# VEP Annotation Guide

## Input contract

Phase 1 input table:
- `results/tables/vep_input.tsv`

Conversion command:
- `python scripts/parse_vep_results.py`

This writes:
- `results/annotations/vep_raw/vep_input_ensembl.tsv`
- `results/annotations/vep_raw/vep_command.sh`

## Local VEP requirements

- `vep` available in PATH
- Local cache under `data/reference/vep/`
- Assembly alignment with project build (GRCh37 by default)

## Raw output location

- `results/annotations/vep_raw/vep_output.tsv`

## Parsed outputs

- `results/annotations/variant_annotations.tsv`
- `results/annotations/variant_consequence_scores.tsv`

Multi-consequence rows are exploded to preserve transcript-aware detail.
