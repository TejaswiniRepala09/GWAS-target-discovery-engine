# Scripts Index

This folder contains executable entry points and utility scripts.

## Main entry points
- `run_pipeline.py`: Phase 1 GWAS preprocessing + lead loci extraction.
- `run_phase2.py`: Phase 2 single-pass integration.
- `run_phase2_multi_locus.py`: controlled multi-locus batch execution.
- `integrate_system_layer.py`: DuckDB integration + optional PostgreSQL export.
- `run_finemap_benchmark.py`: SuSiE vs FINEMAP benchmark workflow.
- `run_coloc_representative_loci.py`: coloc-ready layer for representative loci.

## Utilities
- `parse_vep_results.py`, `parse_susie_results.py`: parser helpers.
- `prepare_fine_mapping_inputs.py`, `run_susie_scaffold.py`: fine-mapping scaffolding.
- `fetch_open_targets.py`, `prepare_structure_candidates.py`, `rank_targets.py`: downstream layers.
- `check_output_contracts.py`: lightweight output schema/contract validation.

## Notes
- Scripts are designed to reuse existing project outputs where possible.
- Heavy computations are not required for documentation/validation passes.
