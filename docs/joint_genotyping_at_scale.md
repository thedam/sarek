# Joint genotyping at scale — design notes

How to run the DeepVariant → GLnexus joint-genotyping path on ~1000 WGS samples
reliably and cost-efficiently on cloud batch (GCP Batch / AWS Batch) with object
storage (GCS / S3) and spot/preemptible VMs.

This complements the Part 1 implementation
([`bam_joint_genotyping_deepvariant`](../subworkflows/local/bam_joint_genotyping_deepvariant/main.nf))
and the prototype [`conf/gcp_batch.config`](../conf/gcp_batch.config).

## 1. The numbers that drive the design

| Quantity            | Per sample | Cohort (×1000) |
| ------------------- | ---------- | -------------- |
| CRAM                | 20–30 GB   | 20–30 TB       |
| gVCF (DeepVariant)  | ~5 GB      | ~5 TB          |
| DeepVariant compute | independent per sample | embarrassingly parallel |
| GLnexus compute     | — | single merge over **all** gVCFs |

Two qualitatively different problems:

- **DeepVariant** is *embarrassingly parallel* — 1000 independent jobs. It scales
  trivially on spot VMs; the only concern is throughput/quota and per-job cost.
- **GLnexus** is a *gather* — it must see every sample at a given locus. This is the
  hard part: one naive job that reads 5 TB and holds the cohort in memory does not
  fit on a single VM, cannot tolerate a spot preemption (it would restart from
  scratch after hours), and is the dominant cost/risk.

The design therefore focuses almost entirely on **how to shard and harden GLnexus**.

## 2. Pipeline shape

```
CRAM (GCS) ──▶ DeepVariant (×1000, spot)  ──▶ per-sample gVCF (GCS)
                                                     │
                              shard by genomic region (N shards)
                                                     ▼
                              GLnexus per shard (×N, spot)  ──▶ per-shard multi-sample BCF
                                                     │
                                          bcftools concat (gather)
                                                     ▼
                              cohort VCF.gz ──▶ VEP (sharded) ──▶ annotated cohort VCF
```

## 3. Challenge → solution

### 3.1 GLnexus is single-node and memory + local-disk bound

GLnexus builds a RocksDB on local disk and keeps working state in RAM; both grow
with `#samples × region_size`. For 1000 WGS samples a whole-genome single run needs
hundreds of GB of RAM and fast local scratch.

**Solution — scatter by genomic region.** Run GLnexus once per region shard (e.g.
per chromosome, or fixed ~10–50 Mb windows), then `bcftools concat` the per-shard
multi-sample VCFs into the cohort VCF. Variants are independent across regions, so
this is exact (no approximation). Benefits:

- **Bounded memory** per shard → fits a standard high-mem VM, no special hardware.
- **Parallelism** → N shards run concurrently instead of one multi-day job.
- **Cheap failure** → a lost shard re-runs in minutes, not from scratch.

The subworkflow already exposes a `bed` input for exactly this — the scale path
passes one shard BED per GLnexus invocation and gathers with `bcftools/concat`
(already vendored in the repo). Sarek's existing `prepare_intervals` /
`--nucleotides_per_second` machinery produces the shard BEDs.

### 3.2 …but GLnexus reads each gVCF *whole*, even with `--bed`

`--bed` restricts the *output* region, but GLnexus still streams every input gVCF in
full. Naively scattering into N shards would read the 5 TB of gVCFs N times — an I/O
and egress disaster.

**Solution — slice the gVCFs by region first (one extra parallel, spot-friendly
step):**

```
per (sample, shard):  bcftools view -r <shard> sample.g.vcf.gz -Oz -o sample.<shard>.g.vcf.gz   (tabix)
per shard:            GLnexus(all samples' <shard> slices, --bed <shard>)
```

Now each GLnexus shard reads only its slice (≈ 5 TB / N), the slicing itself is
1000×N tiny independent tasks ideal for spot, and total bytes read is ~1× the cohort
instead of N×. This is the standard "shard the gather" pattern and is the single most
important optimization at this scale. (For a first iteration, whole-chromosome shards
with `--bed` and no pre-slicing is acceptable — ~24 shards, ~24× reads of a 5 TB set
is tolerable; pre-slicing becomes worth it for finer windows.)

### 3.3 Spot/preemptible reliability

Spot VMs are reclaimed without notice (~60–91% cheaper). Two safeguards:

- **Sharded + idempotent tasks** — Nextflow re-submits a preempted task; because each
  task writes to a unique work dir and is deterministic, retries are safe. Small
  shards mean a retry wastes minutes, not hours.
- **Spot-aware retries** — `google.batch.maxSpotAttempts` re-runs a preempted task on a
  fresh spot VM (separate from `maxRetries`, which covers other transient failures). The
  expensive GLnexus *gather* is the most painful task to lose, so either keep shards small
  (a re-run then costs minutes) or pin that one process to a **non-spot queue** for
  guaranteed completion — there is no per-task spot toggle in a single config. See
  `conf/gcp_batch.config`.
- **Checkpoint via the object store** — published per-shard BCFs in GCS act as
  natural checkpoints; with `-resume` a re-launch skips completed shards.

### 3.4 Storage & data movement (5 TB gVCF, 20–30 TB CRAM)

- **Co-locate compute and bucket in the same region** — cross-region/internet egress
  is the dominant hidden cost at TB scale. Pin `google.location` (or AWS region) to
  the bucket region.
- **Stream from object storage** — GLnexus and bcftools read gVCFs straight from GCS;
  DeepVariant localizes one CRAM per task. Keep `work-dir` on GCS so no control node
  ever holds TBs.
- **Local SSD for GLnexus scratch** — the RocksDB must live on local SSD, not network
  persistent disk, or the DB becomes the bottleneck (`scratch = true` + a large local
  `disk`).
- **Lifecycle policy on intermediates** — gVCF slices and per-shard BCFs are
  reproducible; put them under a GCS prefix with an auto-delete lifecycle rule, or
  clean the work dir after the run.

### 3.5 Cost optimization summary

1. **Spot everywhere** the work is short/idempotent (DeepVariant, slicing, per-shard
   GLnexus) — biggest single lever.
2. **Right-size** each step (DeepVariant: many CPUs; GLnexus: RAM + local SSD;
   bcftools: cheap) instead of one fat machine type for all.
3. **Scatter** to trade a few cheap VMs running in parallel for one expensive
   long-lived VM.
4. **Read the cohort once** via gVCF pre-slicing instead of N× via naive sharding.
5. **Delete reproducible intermediates** and keep everything in-region.

### 3.6 VEP at scale

VEP is also a per-record transform → shard the cohort VCF by region, annotate in
parallel, `bcftools concat`. Use a VEP cache staged once in the bucket (not
downloaded per task). Sarek's `ANNOTATE` path is reused unchanged; only the input is
the single multi-sample joint VCF.

## 4. Recommended starting topology (~1000 WGS)

- DeepVariant: 1000 × spot `n2-standard-32`, ~200 GB disk, retry-on-spot.
- gVCF slicing: 1000 × ~24 shards = ~24k tiny spot tasks (or skip for v1).
- GLnexus: ~24 (per-chromosome) × spot high-mem (`n2-highmem-32`, 208 GB RAM) +
  ~1.5 TB local SSD, `--config DeepVariantWGS`.
- Gather: `bcftools concat` → cohort VCF; VEP sharded; final annotated cohort VCF.

## 5. AWS Batch equivalent

The same shape maps directly:

- `process.executor = 'awsbatch'`, `aws.batch.cliPath`, `work-dir = 's3://…'`.
- Spot via a Batch **compute environment** of Spot instances; on-demand CE as the
  fallback queue for final retries.
- Local NVMe SSD instances (e.g. `m5d`/`r5d`) for GLnexus scratch; size the Batch job
  `memory`/`vcpus` per step.
- S3 lifecycle rules for intermediates; keep buckets and compute in one region.

## 6. What is implemented vs. designed here

- **Implemented (Part 1):** the module, subworkflow (with a `bed` input ready for
  sharding), pipeline integration, annotation reuse, tests, and a local test profile.
- **Designed (Part 2, this doc + `conf/gcp_batch.config`):** the scatter-gather over
  region shards, gVCF pre-slicing, spot hardening, storage/egress strategy and the
  cloud executor settings. The scatter-gather is a natural extension of the existing
  subworkflow (swap the single whole-genome GLnexus call for a per-shard call over
  the `bed` channel + a `BCFTOOLS_CONCAT` gather) and is left as the next step because
  it is best validated against real cloud infrastructure, which is out of scope here.
```
