#!/usr/bin/env bash
set -euo pipefail

# Phase 2 terminal command reference for manual runs.
# This script intentionally prints commands and executes them in order.
# Edit variables below as needed for your environment.

LOCUS_ID="${LOCUS_ID:-chr7_locus10}"
ROOT_DIR="${ROOT_DIR:-$(pwd)}"
VEP_IMAGE="${VEP_IMAGE:-ensemblorg/ensembl-vep:release_115.2}"
VEP_CACHE_URL="${VEP_CACHE_URL:-https://ftp.ensembl.org/pub/grch37/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh37.tar.gz}"

VEP_RAW_DIR="results/annotations/vep_raw"
BUNDLE_DIR="results/fine_mapping/loci_inputs/${LOCUS_ID}"
SUSIE_RAW_DIR="results/fine_mapping/susie_raw/${LOCUS_ID}"
VEP_INPUT_LOCUS="${VEP_RAW_DIR}/${LOCUS_ID}_vep_input_ensembl.tsv"
VEP_OUTPUT_LOCUS="${VEP_RAW_DIR}/vep_output_${LOCUS_ID}.tsv"

run_cmd() {
  printf '\n>>> %s\n' "$*"
  "$@"
}

run_shell() {
  printf '\n>>> %s\n' "$*"
  bash -lc "$*"
}

# 0) Basic checks
run_shell "command -v docker && docker --version && docker info >/dev/null"
run_shell "command -v tabix || true"
run_shell "command -v bcftools || true"
run_shell "command -v plink || true"
run_shell "command -v Rscript || true"

# 1) Pull VEP image (ARM64 compatible tag)
run_cmd docker pull "${VEP_IMAGE}"

# 2) Build per-locus VEP input from summary stats
run_shell "awk -F '\t' 'BEGIN{OFS=\"\\t\"} NR==1{for(i=1;i<=NF;i++){if(\$i==\"CHR\")c=i; if(\$i==\"BP\")b=i; if(\$i==\"A1\")a1=i; if(\$i==\"A2\")a2=i; if(\$i==\"SNP\")s=i;} next} {id=(\$s!=\"\"?\$s:\$c\":\"\$b); print \$c,\$b,\$b,\$a1\"/\"\$a2,\"+\",id}' \"${BUNDLE_DIR}/summary_stats.tsv\" > \"${VEP_INPUT_LOCUS}\""
run_shell "test -s \"${VEP_INPUT_LOCUS}\" && head -n 3 \"${VEP_INPUT_LOCUS}\""

# 3) Rebuild GRCh37 cache (optional/manual; very large download)
# NOTE: This is intentionally explicit because cache integrity is critical.
run_shell "mkdir -p data/reference/vep/tmp"
run_shell "rm -rf data/reference/vep/homo_sapiens/115_GRCh37"
run_shell "curl -fL -C - \"${VEP_CACHE_URL}\" -o data/reference/vep/tmp/homo_sapiens_vep_115_GRCh37.tar.gz"
run_shell "stat -f '%z' data/reference/vep/tmp/homo_sapiens_vep_115_GRCh37.tar.gz"
run_shell "gzip -t data/reference/vep/tmp/homo_sapiens_vep_115_GRCh37.tar.gz"
run_shell "tar -xzf data/reference/vep/tmp/homo_sapiens_vep_115_GRCh37.tar.gz -C data/reference/vep"
run_shell "test -d data/reference/vep/homo_sapiens/115_GRCh37/7"

# 4) Docker VEP run (offline cache mode, GRCh37)
run_cmd docker run --rm \
  -v "${ROOT_DIR}":/work \
  -v "${ROOT_DIR}/data/reference/vep":/opt/vep/.vep \
  "${VEP_IMAGE}" \
  vep \
  --input_file "/work/${VEP_INPUT_LOCUS}" \
  --format ensembl \
  --output_file "/work/${VEP_OUTPUT_LOCUS}" \
  --tab \
  --everything \
  --cache \
  --offline \
  --dir_cache /opt/vep/.vep \
  --assembly GRCh37 \
  --force_overwrite \
  --fork 4 \
  --stats_file "/work/${VEP_RAW_DIR}/vep_summary_${LOCUS_ID}.html"
run_shell "test -s \"${VEP_OUTPUT_LOCUS}\" && wc -l \"${VEP_OUTPUT_LOCUS}\""

# 5) Parse VEP output
run_cmd python3 scripts/parse_vep_results.py --vep-raw-output "${VEP_OUTPUT_LOCUS}"
run_shell "test -s results/annotations/variant_annotations.tsv && wc -l results/annotations/variant_annotations.tsv"
run_shell "test -s results/annotations/variant_consequence_scores.tsv && wc -l results/annotations/variant_consequence_scores.tsv"

# 6) Pilot pipeline resume (starts from dependency checks and includes LD/SuSiE)
run_cmd bash scripts/run_pilot_chr7_locus10.sh
