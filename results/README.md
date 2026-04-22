# Results Index

This directory stores generated project outputs.

## Final/primary outputs
- `reports/`: narrative summaries and interpretation reports.
- `fine_mapping/`: aggregate PIP and credible set outputs.
- `annotations/`: parsed variant annotation tables.
- `target_prioritization/`: variant and gene ranking outputs.
- `benchmarking/`: SuSiE vs FINEMAP comparison tables.
- `database/`: DuckDB file and exported integration tables.
- `plots/benchmarking/`: benchmark plots referenced in the main README.

## Extended outputs
- `coloc/`: coloc-ready representative-locus inputs and status table.
- `loci/`: per-locus organized artifacts.

## Debug/intermediate-heavy outputs
- `debug/`, `logs/`, `phase2_final/`, `plots/phase2/`
- `fine_mapping/loci_inputs/`, `fine_mapping/susie_raw/`, `annotations/vep_raw/`

These are useful locally for diagnostics and reproducibility but are often too large/noisy for a minimal public snapshot.
