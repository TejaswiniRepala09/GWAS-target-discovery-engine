# GWAS Target Discovery Engine

Genetics-driven target discovery workflow for CKD that extends GWAS signal discovery into fine-mapping, annotation, gene prioritization, systems integration, and benchmark comparison.

## What This Project Does

This repository implements an end-to-end translational pipeline:

1. GWAS summary-statistics QC and locus discovery (Phase 1)
2. Per-locus fine-mapping input preparation (LD-aware contracts)
3. VEP annotation parsing and consequence scoring
4. SuSiE fine-mapping parsing (PIP + credible sets)
5. Variant-to-gene and gene-level prioritization
6. Multi-locus execution with per-locus logs/status
7. DuckDB analytical integration and optional PostgreSQL export
8. SuSiE vs FINEMAP benchmark on successful loci

## Repository Layout

```text
config/                     # YAML settings for all stages
src/                        # core reusable pipeline modules
scripts/                    # executable entry points
docs/                       # methods, architecture, file guide
data/                       # raw/interim/reference/external (local resources)
results/
  tables/                   # lead variants/loci and core tables
  annotations/              # parsed VEP outputs
  fine_mapping/             # pip_summary, credible_sets, diagnostics
  target_prioritization/    # variant/gene prioritization tables
  reports/                  # final summaries, interpretation reports
  plots/                    # phase2 + benchmarking figures
  database/                 # DuckDB + exported integration tables
  loci/<locus_id>/          # per-locus organized artifacts
```

## Where To See Final Results

Start here:

- Final project summary: `results/reports/final_project_summary.md`
- Representative loci interpretation: `results/reports/representative_loci_summary.md`
- Benchmark summary: `results/reports/benchmarking_summary.md`

Core output tables:

- `results/fine_mapping/pip_summary.tsv`
- `results/fine_mapping/credible_sets.tsv`
- `results/annotations/variant_annotations.tsv`
- `results/target_prioritization/variant_priority_scores.tsv`
- `results/target_prioritization/gene_prioritization.tsv`
- `results/benchmarking/susie_vs_finemap_locus_summary.tsv`
- `results/database/target_discovery.duckdb`

Key plots:

- `results/plots/benchmarking/susie_vs_finemap_runtime.png`
- `results/plots/benchmarking/susie_vs_finemap_top_variant_overlap.png`
- `results/plots/benchmarking/susie_vs_finemap_credible_set_sizes.png`

## Reproducibility Notes

- Assembly/build used in this project is **GRCh37**.
- LD references are 1000 Genomes Phase 3 chromosome VCFs.
- VEP is run via Docker in offline cache mode (GRCh37 cache).
- SuSiE and FINEMAP require external binaries/runtime (R + FINEMAP setup).
- Association/fine-mapping outputs are prioritization evidence, not mechanistic proof.

## Environment Setup

```bash
cd /path/to/GWAS-target-discovery-engine
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Exact Pipeline Run Order (Replicate Project)

### 1) Phase 1 GWAS discovery

```bash
source .venv/bin/activate
python scripts/run_pipeline.py
```

Expected key outputs:

- `data/interim/ckdgen_egfr_cleaned.tsv`
- `results/tables/lead_variants.tsv`
- `results/tables/lead_loci.tsv`
- `results/tables/vep_input.tsv`

### 2) Phase 2 single pass integration

```bash
source .venv/bin/activate
python scripts/run_phase2.py
```

### 3) Multi-locus controlled batch run

```bash
source .venv/bin/activate
python scripts/run_phase2_multi_locus.py --mode top_n --max-loci 20 --continue-on-error
```

### 4) System integration layer (DuckDB, optional PostgreSQL)

```bash
source .venv/bin/activate
python scripts/integrate_system_layer.py
```

Optional PostgreSQL export:

```bash
export TARGET_DISCOVERY_PG_DSN='postgresql://user:password@localhost:5432/target_discovery'
python scripts/integrate_system_layer.py --postgres-dsn "$TARGET_DISCOVERY_PG_DSN"
```

### 5) Benchmark (SuSiE vs FINEMAP)

```bash
source .venv/bin/activate
bash scripts/run_benchmark.sh
```

## VEP (Docker) and SuSiE External Execution

This repo includes helpers and contracts for external tools.

- VEP parsing: `scripts/parse_vep_results.py`
- SuSiE parsing: `scripts/parse_susie_results.py`
- Pilot and utility scripts are in `scripts/` for reproducible execution traces.

See methods docs:

- `docs/vep_annotation.md`
- `docs/fine_mapping.md`
- `docs/phase2_methodology.md`

## System / Architecture Docs

- Architecture: `docs/project_architecture.md`
- File map: `docs/project_file_guide.md`
- Final report: `results/reports/final_project_summary.md`

## Limitations

- Fine-mapping quality depends on LD consistency and conditioning.
- Some loci remain numerically unstable and are flagged in status reports.
- Non-coding loci have gene-assignment uncertainty.
- eQTL/regulatory context is supportive, not definitive causality proof.

## License / Use

This repository is intended for research, reproducible methods demonstration, and portfolio/interview presentation.
