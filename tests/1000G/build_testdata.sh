#!/usr/bin/env bash
#
# build_testdata.sh — create a small joint-genotyping test cohort from 1000 Genomes
# ---------------------------------------------------------------------------------
# Downloads the high-coverage (30x) GRCh38 CRAMs of a handful of 1000 Genomes
# samples, subsets each to a single chromosome (or region) and writes a sarek
# samplesheet. The full CRAMs are 20-30 GB each, so we never keep them: samtools
# streams the remote CRAM over HTTPS and only the requested region is downloaded
# (a few hundred MB per sample for a whole chromosome, far less for a sub-region).
#
# The input data is intentionally NOT committed to the repository — run this
# script once to materialise it for review.
#
# Requirements: samtools >= 1.13, curl, awk  (e.g. `conda install -c bioconda samtools`)
#
# Usage:
#   tests/1000G/build_testdata.sh -r /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa
#
# Options (env vars or flags):
#   -r REFERENCE   GRCh38 analysis-set FASTA the CRAMs were aligned against (REQUIRED).
#                  This MUST be the same assembly used by the 1000G high-coverage CRAMs.
#                  Download once with:  (see ${REF_URL} below)
#   -s SAMPLES     Space-separated 1000G sample IDs   (default: "HG00096 HG00097 HG00099")
#   -c REGION      Chromosome or region to subset to  (default: "chr21")
#                  Use e.g. "chr21:1-5000000" for an even faster, tiny test.
#   -o OUTDIR      Output directory                   (default: tests/1000G/data)
#   -t THREADS     samtools threads                   (default: 4)
#
# Output:
#   ${OUTDIR}/<sample>.<region>.cram(+.crai)   subset CRAMs
#   ${OUTDIR}/joint_genotype_1000G.csv         sarek samplesheet (step: variant_calling)
#
# Run sarek on the result:
#   nextflow run . -profile docker \
#       --step variant_calling --tools deepvariant,vep --joint_genotype \
#       --input ${OUTDIR}/joint_genotype_1000G.csv \
#       --fasta /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa --igenomes_ignore \
#       --outdir results_1000G
#
set -euo pipefail

# --- defaults ---------------------------------------------------------------
SAMPLES=${SAMPLES:-"HG00096 HG00097 HG00099"}
REGION=${REGION:-"chr21"}
OUTDIR=${OUTDIR:-"$(dirname "$0")/data"}
THREADS=${THREADS:-4}
REFERENCE=${REFERENCE:-""}

# 1000G high-coverage (30x, GRCh38) data collection + its index (run -> sample -> cram path)
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
INDEX_URL="${BASE_URL}/1000G_2504_high_coverage.sequence.index"
REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# --- args -------------------------------------------------------------------
while getopts "r:s:c:o:t:h" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        s) SAMPLES="$OPTARG" ;;
        c) REGION="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Invalid option; run with -h for help" >&2; exit 1 ;;
    esac
done

if [[ -z "${REFERENCE}" ]]; then
    echo "ERROR: a GRCh38 analysis-set FASTA is required (-r). Download it once with:" >&2
    echo "  curl -O ${REF_URL} && curl -O ${REF_URL}.fai && samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa" >&2
    exit 1
fi
command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

mkdir -p "${OUTDIR}"
SAMPLESHEET="${OUTDIR}/joint_genotype_1000G.csv"
REGION_TAG=${REGION//[:]/_}

echo "Resolving CRAM URLs from ${INDEX_URL} ..."
INDEX_FILE="${OUTDIR}/1000G_2504_high_coverage.sequence.index"
[[ -f "${INDEX_FILE}" ]] || curl -fsSL "${INDEX_URL}" -o "${INDEX_FILE}"

echo "patient,status,sample,cram,crai" > "${SAMPLESHEET}"

for SAMPLE in ${SAMPLES}; do
    # Column 1 of the index is the ENA/analysis path; column with SAMPLE_NAME identifies the sample.
    # Pick the analysis CRAM (path ending in .cram) for this sample.
    REL_PATH=$(awk -v s="${SAMPLE}" -F'\t' '$0 ~ s && $1 ~ /\.cram$/ {print $1; exit}' "${INDEX_FILE}")
    if [[ -z "${REL_PATH}" ]]; then
        echo "  WARNING: no CRAM found for ${SAMPLE} in the index, skipping" >&2
        continue
    fi
    CRAM_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/${REL_PATH#ftp.1000genomes.ebi.ac.uk/vol1/ftp/}"
    OUT_CRAM="${OUTDIR}/${SAMPLE}.${REGION_TAG}.cram"

    echo "Subsetting ${SAMPLE} ${REGION} (streaming from ${CRAM_URL}) ..."
    # Force CRAM 3.0: recent samtools defaults to CRAM 3.1, whose new codecs cannot be
    # decoded by the older htslib bundled in some tool containers (e.g. mosdepth -> empty
    # coverage). CRAM 3.0 is readable everywhere.
    samtools view -@ "${THREADS}" -C --output-fmt-option version=3.0 -T "${REFERENCE}" -o "${OUT_CRAM}" "${CRAM_URL}" "${REGION}"
    samtools index -@ "${THREADS}" "${OUT_CRAM}"

    echo "${SAMPLE},0,${SAMPLE},$(readlink -f "${OUT_CRAM}"),$(readlink -f "${OUT_CRAM}").crai" >> "${SAMPLESHEET}"
done

echo
echo "Done. Samplesheet: ${SAMPLESHEET}"
cat "${SAMPLESHEET}"
