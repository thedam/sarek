#!/usr/bin/env bash
#
# download_testdata.sh
# ---------------------------------------------------------------------------
# Fetch a tiny, self-contained test cohort for the DeepVariant joint-genotyping
# extension of sarek (`--tools deepvariant --joint_genotype`).
#
# It downloads 3 unrelated-as-cohort WGS samples (the CEU trio) from the
# 1000 Genomes Project, subsets each alignment to a single chromosome (chr21),
# rewrites the header to that single contig, and emits CRAM + index plus a
# GRCh37 chr21 reference and a ready-to-use sarek samplesheet.
#
# NOTE: no large input data is committed to the repo — run this script to
# (re)create the bundle locally for review.
#
# Requirements: bash, curl, samtools (>=1.10, built with libcurl), gzip.
#   conda create -n testdata -c bioconda -c conda-forge samtools curl
#
# Usage:
#   bash tests/download_testdata.sh
#   REGION=21 bash tests/download_testdata.sh        # whole chromosome (slower)
# ---------------------------------------------------------------------------
set -euo pipefail

# --- configuration ---------------------------------------------------------
# CEU trio — three samples used here as a small joint-genotyping cohort.
SAMPLES=(NA12878 NA12891 NA12892)
declare -A SEX=( [NA12878]=XX [NA12891]=XY [NA12892]=XX )

CHR="21"                                   # contig name as used in 1000G GRCh37 BAMs
REGION="${REGION:-21:14000000-16000000}"   # subset window; set REGION=21 for full chr
OUTDIR="${OUTDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/testdata}"

G1K_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"
REF_URL="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.21.fa.gz"

ALN_DIR="$OUTDIR/cram"
REF_DIR="$OUTDIR/reference"
REF_FA="$REF_DIR/Homo_sapiens.GRCh37.chr21.fa"
SAMPLESHEET="$OUTDIR/samplesheet.csv"

mkdir -p "$ALN_DIR" "$REF_DIR"

command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH"; exit 1; }
command -v curl     >/dev/null || { echo "ERROR: curl not found in PATH";     exit 1; }

# --- 1. reference (GRCh37 chr21) -------------------------------------------
if [[ ! -f "$REF_FA" ]]; then
  echo ">> Downloading GRCh37 chr21 reference ..."
  curl -fsSL "$REF_URL" | gzip -dc > "$REF_FA"
fi
[[ -f "$REF_FA.fai"  ]] || samtools faidx "$REF_FA"
[[ -f "${REF_FA%.fa}.dict" ]] || samtools dict "$REF_FA" -o "${REF_FA%.fa}.dict"

# --- 2. per-sample alignment subsets ---------------------------------------
for s in "${SAMPLES[@]}"; do
  out_cram="$ALN_DIR/${s}.chr21.cram"
  if [[ -f "$out_cram" ]]; then
    echo ">> $s already present, skipping."
    continue
  fi

  # Resolve the BAM URL from the directory listing (robust to date suffixes).
  # We use the high-coverage PCR-free alignments: they exist for all three
  # samples (low-coverage alignment/ is missing for NA12891/NA12892) and give
  # depth suitable for DeepVariant. Still GRCh37 (hs37d5).
  echo ">> Resolving alignment URL for $s ..."
  listing="$(curl -fsSL "$G1K_BASE/$s/high_coverage_alignment/")"
  bam="$(printf '%s\n' "$listing" \
          | grep -oE "${s}\.mapped[^\"'<> ]*high_coverage[^\"'<> ]*\.bam" \
          | head -n1)"
  [[ -n "$bam" ]] || { echo "ERROR: could not find a high_coverage BAM for $s"; exit 1; }
  bam_url="$G1K_BASE/$s/high_coverage_alignment/$bam"

  echo ">> Fetching $REGION for $s ($bam) ..."
  # Stream only the requested region, then:
  #   - keep @SQ for chr21 only (drop other contigs),
  #   - sanitise mate references that point to other chromosomes (-> '*'),
  # so the result is a clean single-contig alignment DeepVariant accepts.
  samtools view -h "$bam_url" "$REGION" \
    | awk -v c="$CHR" 'BEGIN{OFS="\t"}
        /^@SQ/ { if ($0 ~ "\tSN:"c"\t" || $0 ~ "\tSN:"c"$") print; next }
        /^@/   { print; next }
        { if ($7 != "=" && $7 != c) { $7="*"; $8="0" } print }' \
    | samtools sort -O bam -o "$ALN_DIR/${s}.chr21.bam" -

  # Convert to CRAM against our chr21 reference (sarek --step variant_calling input).
  # Force CRAM 3.0: recent samtools defaults to CRAM 3.1, whose new codecs cannot be
  # decoded by the older htslib in some tool containers (e.g. mosdepth -> empty coverage).
  samtools view -C --output-fmt-option version=3.0 -T "$REF_FA" "$ALN_DIR/${s}.chr21.bam" -o "$out_cram"
  samtools index "$out_cram"
  rm -f "$ALN_DIR/${s}.chr21.bam"
  echo ">> Wrote $out_cram"
done

# --- 3. sarek samplesheet --------------------------------------------------
echo ">> Writing $SAMPLESHEET"
{
  echo "patient,sample,sex,status,cram,crai"
  for s in "${SAMPLES[@]}"; do
    echo "${s},${s},${SEX[$s]},0,${ALN_DIR}/${s}.chr21.cram,${ALN_DIR}/${s}.chr21.cram.crai"
  done
} > "$SAMPLESHEET"

cat <<EOF

Done. Test bundle in: $OUTDIR
  reference : $REF_FA
  alignments: $ALN_DIR/*.chr21.cram
  samplesheet: $SAMPLESHEET

Run sarek with, e.g.:
  nextflow run . -profile docker \\
    --input $SAMPLESHEET \\
    --fasta $REF_FA \\
    --step variant_calling \\
    --tools deepvariant --joint_genotype \\
    --outdir results
EOF
