# Fine-Mapping (SuSiE) Guide

## Inputs reused from Phase 1

- `data/interim/ckdgen_egfr_cleaned.tsv`
- `results/tables/lead_variants.tsv`
- `results/tables/lead_loci.tsv`
- `results/loci/*_window.tsv` (if already generated)

## Bundle outputs

`python scripts/prepare_fine_mapping_inputs.py`

Writes bundles under `results/fine_mapping/loci_inputs/<locus_id>/`:
- `summary_stats.tsv`
- `zscores.tsv`
- `metadata.yaml`
- `ld_matrix.tsv` (placeholder contract path)

## LD contract

- Chromosome VCF resources expected in `data/reference/ld/` or `data/reference/Ld/`
- Per-locus LD matrix must be square, SNP-aligned, and match summary SNP ordering.

## SuSiE scaffold

`python scripts/run_susie_scaffold.py`

Creates:
- `scripts/susie_run_template.R`
- `results/fine_mapping/susie_commands.tsv`

Run per locus using generated command entries.

## Parse real SuSiE outputs

Place outputs under `results/fine_mapping/susie_raw/<locus_id>/`:
- `pip.tsv`
- `credible_sets.tsv`

Then run:
- `python scripts/parse_susie_results.py`
