# Benchmarking Summary (SuSiE vs FINEMAP)

## Scope
- Successful loci benchmark cohort from `multi_locus_status.tsv`: **13**
- GWAS/VEP/SuSiE were not rerun; existing outputs were reused.

## Benchmark Outcome
- Loci attempted: **13**
- Loci benchmarked successfully (FINEMAP completed): **13**
- Loci failed: **0**
- Top-variant agreement (among comparable loci): **1/13**

## Method Stability
- Stable loci: loci with matching SuSiE and FINEMAP top variants.
- Method-sensitive loci: loci where top variants differ, or FINEMAP failed.
- Runtime comparison plotted where runtime values were available.

## Runtime Notes
- SuSiE runtime was not consistently logged in a parseable per-locus format; values may be null.
- FINEMAP runtime is measured from command execution when FINEMAP runs.

## Top Blockers
- None

## Output Files
- `results/benchmarking/susie_vs_finemap_variant_comparison.tsv`
- `results/benchmarking/susie_vs_finemap_locus_summary.tsv`
- `results/plots/benchmarking/susie_vs_finemap_runtime.png`
- `results/plots/benchmarking/susie_vs_finemap_top_variant_overlap.png`
- `results/plots/benchmarking/susie_vs_finemap_credible_set_sizes.png`