# Solution walkthrough — DeepVariant joint genotyping for sarek

Recruitment task: extend nf-core/sarek so that the DeepVariant path can do **cohort
joint genotyping** and emit a single multi-sample, VEP-annotated VCF — and design how
to run it at ~1000-WGS scale on the cloud.

This document is the reviewer's entry point. PR: **https://github.com/thedam/sarek/pull/1**

---

## TL;DR

```bash
export NXF_VER=25.10.2   # required: >=25.10.2 and <26.04 (see "Gotchas")

# Only needed if your shell has AWS credentials (see "Gotchas"): a tiny override so
# nf-schema does not try to validate sarek's default public-S3 cache paths.
cat > /tmp/sarek_no_s3.config <<'EOF'
params { snpeff_cache = null; vep_cache = null; igenomes_base = null; sentieon_dnascope_model = null }
EOF

# Minimal end-to-end test (3 samples, single chromosome, runs in minutes)
nextflow run . -profile test_joint_genotype,docker --outdir results_jg -c /tmp/sarek_no_s3.config
```

Headline command on real input:

```bash
nextflow run . -profile docker --step variant_calling \
    --tools deepvariant,vep --joint_genotype --input <cohort.csv> --genome GATK.GRCh38 \
    --outdir results
```

`--joint_genotype` merges every sample's DeepVariant gVCF with **GLnexus** into one
multi-sample VCF, which is then VEP-annotated through sarek's existing `ANNOTATE` path.

---

## Part 1 — what was implemented

### Tool choice: GLnexus
GLnexus is the joint genotyper recommended for DeepVariant gVCFs (it ships a
`DeepVariant*` config preset). GATK `GenomicsDBImport`/`GenotypeGVCFs` assumes
HaplotypeCaller-style gVCFs and does not handle DeepVariant likelihoods well.

### Data flow
```
CRAM ──▶ DeepVariant (per sample) ──▶ per-sample gVCF
                                            │  collect whole cohort
                                            ▼
                                        GLnexus  ──▶ multi-sample BCF
                                            │  bcftools view -Oz
                                            ▼
                                  multi-sample VCF.gz + tbi
                                            │  (reuses sarek ANNOTATE)
                                            ▼
                                  VEP-annotated cohort VCF
```

### File map
| Piece | Path |
| --- | --- |
| GLnexus module (vendored from nf-core) | [`../modules/nf-core/glnexus/main.nf`](../modules/nf-core/glnexus/main.nf) |
| **Joint genotyping subworkflow** | [`../subworkflows/local/bam_joint_genotyping_deepvariant/main.nf`](../subworkflows/local/bam_joint_genotyping_deepvariant/main.nf) |
| `--joint_genotype` param | [`../nextflow.config`](../nextflow.config) + [`../nextflow_schema.json`](../nextflow_schema.json) |
| Wiring into variant calling (the call site) | [`../subworkflows/local/bam_variant_calling_germline_all/main.nf#L125-L132`](../subworkflows/local/bam_variant_calling_germline_all/main.nf#L125-L132) |
| Param threaded from the pipeline | [`../workflows/sarek/main.nf#L427`](../workflows/sarek/main.nf#L427) |
| Validation (requires `deepvariant`) | [`../subworkflows/local/samplesheet_to_channel/main.nf`](../subworkflows/local/samplesheet_to_channel/main.nf) |
| Module-level config (GLnexus preset, output naming) | [`../conf/modules/joint_genotype.config`](../conf/modules/joint_genotype.config) |
| Annotation reuse (joint VCF → VEP) | `vcf_all` → `POST_VARIANTCALLING` → `VCF_ANNOTATE_ALL` in [`../workflows/sarek/main.nf`](../workflows/sarek/main.nf) |

### Design decisions
- **Separate `--joint_genotype` flag** (not the existing `--joint_germline`, which is
  GATK/Sentieon-specific): keeps the two joint-calling paths independent.
- **Joint VCF replaces the single-sample VCFs** in the output channel (mirrors how the
  GATK joint path replaces `vcf_haplotypecaller`), so the cohort VCF becomes the thing
  that flows into annotation and the variants CSV.
- **Optional `bed` input** on the subworkflow — already in place so the same path can be
  sharded by region at scale (see Part 2).

---

## Part 2 — running at scale (~1000 WGS samples)

- Prototype cloud profile: [`../conf/gcp_batch.config`](../conf/gcp_batch.config) — GCP
  Batch, spot/preemptible VMs, **local-SSD scratch for the GLnexus RocksDB**,
  retry-on-preemption, right-sized resources, GCS work dir.
- Design notes: [`joint_genotyping_at_scale.md`](joint_genotyping_at_scale.md) — the core
  ideas: **scatter GLnexus by genomic region + pre-slice the gVCFs** (so the 5 TB cohort
  is read once, not N times), spot hardening, storage/egress strategy, and the AWS
  equivalent.

---

## How to run

`export NXF_VER=25.10.2` first in every shell.

### 1. Minimal test (CI-style, minutes) — verified ✅
```bash
nextflow run . -profile test_joint_genotype,docker --outdir results_jg -c /tmp/sarek_no_s3.config
```
(`-c /tmp/sarek_no_s3.config` is only needed if your shell has AWS credentials — create it as
shown in the [TL;DR](#tldr) / [Gotchas](#known-environment-gotchas-not-feature-bugs); on a clean
environment you can drop it.)

Produces `results_jg/variant_calling/deepvariant/joint_genotyping/joint_genotyping.deepvariant.vcf.gz`
with 3 sample columns.

### 2. Real data — 1000 Genomes CEU trio (chr21, GRCh37) — verified ✅

**Step 2a — download & build the cohort.** Input data is not committed; regenerate it
with the helper script (needs `samtools` + `curl` on PATH):

```bash
bash tests/download_testdata_3crams_GRCh37.sh            # fast 2 Mb window (default)
# REGION=21 bash tests/download_testdata_3crams_GRCh37.sh   # whole chr21 (slower)
```

It streams the CEU trio (NA12878 / NA12891 / NA12892) high-coverage alignments from the
1000 Genomes **phase3** FTP, subsets each to chr21 (**GRCh37**, contig `21`), writes
**CRAM 3.0**, fetches the GRCh37 chr21 reference, and emits a ready-to-use samplesheet.
Everything lands under `tests/testdata/` (gitignored):

```
tests/testdata/
├── cram/             NA12878.chr21.cram, NA12891.chr21.cram, NA12892.chr21.cram (+ .crai)
├── reference/        Homo_sapiens.GRCh37.chr21.fa (+ .fai / .dict)
└── samplesheet.csv   <-- sarek input  (patient,sample,sex,status,cram,crai)
```

**Step 2b — run joint genotyping + VEP** on the generated samplesheet. This is the full
deliverable: one multi-sample, **VEP-annotated** cohort VCF. The trio is GRCh37 (contig
`21`), so point sarek at the matching chr21 reference (not `--genome GATK.GRCh38`) and at a
GRCh37 VEP cache:

```bash
nextflow run . -profile docker --step variant_calling --tools deepvariant,vep --joint_genotype \
    --input     tests/testdata/samplesheet.csv \
    --fasta     $PWD/tests/testdata/reference/Homo_sapiens.GRCh37.chr21.fa \
    --fasta_fai $PWD/tests/testdata/reference/Homo_sapiens.GRCh37.chr21.fa.fai \
    --dict      $PWD/tests/testdata/reference/Homo_sapiens.GRCh37.chr21.dict \
    --igenomes_ignore --skip_tools vcftools \
    --download_cache --vep_genome GRCh37 --vep_species homo_sapiens --vep_cache_version 115 \
    --outdir results_trio -c /tmp/sarek_no_s3.config
```
Notes:
- `--download_cache` fetches the GRCh37 v115 VEP cache from Ensembl; pass `--vep_cache <dir>`
  instead to reuse an existing cache.
- `-c /tmp/sarek_no_s3.config` only if your shell has AWS credentials (see [Gotchas](#known-environment-gotchas-not-feature-bugs)).
- To run **joint genotyping only** (skip the VEP cache download), drop `vep` and the
  `--vep_*` / `--download_cache` flags.
- On a re-run add `-resume`: DeepVariant + GLnexus are reused and only VEP re-runs.

Result — a single multi-sample, **VEP-annotated** cohort VCF (NA12878/NA12891/NA12892):
```bash
V=$(find results_trio/annotation -name "*joint_genotyping*VEP*vcf.gz" | head -1)
zcat "$V" | grep -m1 '^#CHROM'    # 3 sample columns: NA12878 NA12891 NA12892
zcat "$V" | grep -m1 'ID=CSQ'     # VEP annotation header (CSQ)
```
(The intermediate un-annotated joint VCF — **89,486 variants on chr21** with `REGION=21` —
is at `results_trio/variant_calling/deepvariant/joint_genotyping/joint_genotyping.deepvariant.vcf.gz`.)

---

## Test data (not committed — regenerate with a script)

Per the task, input data is **not** committed; scripts regenerate it.

| Script | What it produces |
| --- | --- |
| [`../tests/download_testdata_3crams_GRCh37.sh`](../tests/download_testdata_3crams_GRCh37.sh) | **(used for the verified real-data run, §2)** streams the 1000G **phase3** CEU trio (NA12878/91/92), subsets to chr21 (**GRCh37**), writes CRAM 3.0 + GRCh37 chr21 reference + samplesheet → `tests/testdata/` |
| [`../tests/1000G/build_testdata.sh`](../tests/1000G/build_testdata.sh) | alternative: streams 30x high-coverage **GRCh38** CRAMs, subsets to one chromosome, → `tests/1000G/data/` (matches sarek's default `--genome GATK.GRCh38`) |

Both write **CRAM 3.0** (recent samtools defaults to 3.1, which the older htslib in some
tool containers cannot decode). Generated data is gitignored.

Samplesheets:
- Minimal test (committed): [`../tests/csv/3.0/joint_genotype_cram.csv`](../tests/csv/3.0/joint_genotype_cram.csv) — 3 samples from nf-core test-datasets.
- Real trio (generated by the script above): `tests/testdata/samplesheet.csv` (gitignored).

---

## Tests

| Test | Location | Note |
| --- | --- | --- |
| Subworkflow nf-test (+ snapshot) | [`../subworkflows/local/bam_joint_genotyping_deepvariant/tests/main.nf.test`](../subworkflows/local/bam_joint_genotyping_deepvariant/tests/main.nf.test) | runs anywhere, no cloud creds — **passing** |
| Pipeline nf-test | [`../tests/joint_genotype_deepvariant.nf.test`](../tests/joint_genotype_deepvariant.nf.test) | CI-oriented (uses `-profile test`, public S3) |
| GLnexus module test | `modules/nf-core/glnexus/tests/` | from nf-core (gitignored per convention); run `nf-core modules test glnexus` |

```bash
nf-test test subworkflows/local/bam_joint_genotyping_deepvariant/tests/main.nf.test --profile docker
```

---

## Known environment gotchas (not feature bugs)

- **Nextflow version**: pipeline requires `>=25.10.2 <26.04`. The 26.04 edge build's strict
  config parser rejects sarek's conditional `withName` blocks. Pin `NXF_VER=25.10.2`.
- **S3 403 on default cache paths** (`--snpeff_cache`, `--vep_cache`, `--igenomes_base`,
  `--sentieon_dnascope_model`): sarek defaults these to **public** S3 / igenomes buckets.
  The 403 appears only when the shell has an **existing AWS profile / credentials**
  (`AWS_PROFILE`, `~/.aws/...`, `AWS_ACCESS_KEY_ID`): the S3 client then *authenticates* with
  those credentials and is denied, instead of reading the public buckets anonymously. Three
  ways out:
  - run with the AWS profile unset:
    `env -u AWS_PROFILE -u AWS_ACCESS_KEY_ID -u AWS_SECRET_ACCESS_KEY nextflow run ...`
  - force anonymous S3: add `aws.client.anonymous = true` via `-c`
  - drop the unused params with a tiny `no_s3.config` and add `-c no_s3.config` (these params
    are only consumed by annotation/Sentieon, not by `--igenomes_ignore` runs):
    ```groovy
    params { snpeff_cache = null; vep_cache = null; igenomes_base = null; sentieon_dnascope_model = null }
    ```
- **vcftools QC** (`--skip_tools vcftools`): vcftools 0.1.16 SIGABRTs on minimal multi-sample
  VCFs; unrelated to joint genotyping.
- **mosdepth coverage = 0 on CRAM 3.1**: recent samtools writes CRAM 3.1, whose codecs the
  older htslib in the mosdepth container can't decode (DeepVariant reads it fine — calling is
  unaffected). The test-data script writes **CRAM 3.0** to avoid this.

---

## Status

| Deliverable | Status |
| --- | --- |
| GLnexus module (container/io/versions) | ✅ |
| Joint genotyping subworkflow across the cohort | ✅ |
| `--joint_genotype` integration into VARIANT_CALLING | ✅ |
| ANNOTATE integration → VEP-annotated multi-sample VCF | ✅ implemented & wired |
| nf-test for module + subworkflow | ✅ |
| Local end-to-end test profile | ✅ verified |
| Fork / PR-ready branch | ✅ PR #1 |
| No input data committed + regeneration script (1000G single-chr) | ✅ |
| Part 2: GCP/AWS config prototype + scale design | ✅ |
