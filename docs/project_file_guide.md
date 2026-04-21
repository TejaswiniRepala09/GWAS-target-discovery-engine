# Project File Guide

This is a practical map of the repository for collaborators, reviewers, and interview discussion.

## Top-Level Folders

- `config/`
  - YAML settings for column mappings, fine-mapping, VEP, Open Targets, and structure configuration.
- `data/`
  - `raw/`: original input files.
  - `interim/`: cleaned intermediate files.
  - `reference/`: LD references, gene annotations, VEP cache.
  - `external/`: cached API responses and external resources.
- `src/`
  - Reusable modules grouped by domain (`fine_mapping`, `annotation`, `target_prioritization`, `structure`, `utils`).
- `scripts/`
  - Pipeline and utility entry points.
- `results/`
  - All generated outputs (tables, plots, reports, per-locus artifacts, logs, benchmark artifacts, databases).
- `docs/`
  - Methodology, architecture, and interpretive guides.

## Key Execution Scripts

- `scripts/run_pipeline.py`: Phase 1 pipeline.
- `scripts/run_phase2.py`: single-pass Phase 2 integration.
- `scripts/run_phase2_multi_locus.py`: multi-locus run orchestration.
- `scripts/integrate_system_layer.py`: DuckDB + optional PostgreSQL export.
- `scripts/run_finemap_benchmark.py`: SuSiE vs FINEMAP benchmarking.
- `scripts/setup_finemap_docker.sh`: FINEMAP Docker wrapper setup.
- `scripts/run_benchmark.sh`: benchmark runner with required environment defaults.

## Most Important Output Files

### Fine-mapping
- `results/fine_mapping/pip_summary.tsv`
- `results/fine_mapping/credible_sets.tsv`
- `results/fine_mapping/loci_inputs/<locus_id>/summary_stats_harmonized.tsv`
- `results/fine_mapping/loci_inputs/<locus_id>/ld_matrix.tsv`

### Annotation and prioritization
- `results/annotations/variant_annotations.tsv`
- `results/annotations/variant_consequence_scores.tsv`
- `results/target_prioritization/variant_priority_scores.tsv`
- `results/target_prioritization/gene_prioritization.tsv`

### Per-locus organization
- `results/loci/<locus_id>/inputs/`
- `results/loci/<locus_id>/annotations/`
- `results/loci/<locus_id>/ld/`
- `results/loci/<locus_id>/susie/`
- `results/loci/<locus_id>/prioritization/`
- `results/loci/<locus_id>/plots/`
- `results/loci/<locus_id>/logs/`

### Benchmarking
- `results/benchmarking/susie_vs_finemap_variant_comparison.tsv`
- `results/benchmarking/susie_vs_finemap_locus_summary.tsv`
- `results/plots/benchmarking/susie_vs_finemap_runtime.png`
- `results/plots/benchmarking/susie_vs_finemap_top_variant_overlap.png`
- `results/plots/benchmarking/susie_vs_finemap_credible_set_sizes.png`

### Reports
- `results/reports/final_project_summary.md`
- `results/reports/benchmarking_summary.md`
- `results/reports/representative_loci_summary.md`
- `results/reports/chr7_locus19_interpretation.md`
- `results/reports/chr7_locus14_interpretation.md`
- `results/reports/chr7_locus10_interpretation.md`

### System integration
- `results/database/target_discovery.duckdb`
- `results/database/target_variant_table.tsv`
- `results/database/target_gene_table.tsv`

## Debug vs Final Output Convention

- Final presentation plots:
  - `results/plots/phase2/`
  - `results/plots/benchmarking/`
  - `results/phase2_final/plots/`
- Debug/intermediate artifacts:
  - `results/debug/`
  - `results/fine_mapping/diagnostics/`
  - `results/logs/`

## Suggested Entry Order for New Readers

1. `README.md`
2. `docs/project_architecture.md`
3. `results/reports/final_project_summary.md`
4. `results/reports/representative_loci_summary.md`
5. Benchmarking outputs in `results/benchmarking/`
