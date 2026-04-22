# GWAS Target Discovery Engine

Genetics-driven target discovery workflow for CKD that extends GWAS signal discovery into fine-mapping, annotation, gene prioritization, systems integration, and benchmark comparison.

## Quickstart in 5 Commands

```bash
git clone https://github.com/TejaswiniRepala09/GWAS-target-discovery-engine.git
cd GWAS-target-discovery-engine
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
python scripts/run_phase2_multi_locus.py --mode top_n --max-loci 20 --continue-on-error
python scripts/integrate_system_layer.py && bash scripts/run_benchmark.sh
```

Notes:
- External prerequisites for full reproducibility include PLINK, bcftools/tabix, R + `susieR`, Docker (for VEP), and FINEMAP wrapper/binary.
- If you only want to inspect completed outputs, skip compute and open files under `results/reports/`, `results/fine_mapping/`, `results/target_prioritization/`, and `results/database/`.

## Key Results at a Glance

### 1) Representative locus association pattern (chr7_locus19)
![chr7_locus19 locus plot](results/plots/phase2/chr7_locus19_locus_plot.png)

### 2) Fine-mapping posterior support (chr7_locus19 PIP)
![chr7_locus19 PIP plot](results/plots/phase2/chr7_locus19_pip_plot.png)

### 3) Cross-method benchmark agreement (SuSiE vs FINEMAP)
![SuSiE vs FINEMAP top-variant overlap](results/plots/benchmarking/susie_vs_finemap_top_variant_overlap.png)

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

## Reproducibility Matrix

| Step | Runs locally? | Requires external dependency? | Main tools | Notes |
|---|---|---|---|---|
| GWAS preprocessing / lead loci | Yes | No (if input TSV already present) | Python, Polars | Requires GWAS summary stats file placement under `data/raw/`. |
| VEP annotation | Yes | Yes | Docker VEP (`ensemblorg/ensembl-vep`), GRCh37 cache | Offline cache mode expected; reference cache under `data/reference/vep/`. |
| LD extraction | Yes | Yes | bcftools, tabix, 1000G VCFs | Region-based extraction from chromosome VCF + `.tbi` index. |
| SuSiE fine-mapping | Yes | Yes | Rscript, `susieR`, PLINK outputs | Requires harmonized summary stats + LD matrix consistency. |
| DuckDB integration | Yes | No | Python, DuckDB, Polars | Reads existing TSV outputs; does not rerun heavy analysis. |
| PostgreSQL export | Optional | Yes | `psycopg`, PostgreSQL instance | Triggered via `TARGET_DISCOVERY_PG_DSN` or `--postgres-dsn`. |
| FINEMAP benchmarking | Yes | Yes | FINEMAP wrapper/binary, Python, plotting libs | Uses existing successful-locus cohort and continue-on-error behavior. |

## Data Model / Database Schema

### DuckDB integration tables
- `target_variant_table` (`results/database/target_variant_table.tsv`)
  - Variant-centric integrated view from prioritization + fine-mapping + status.
  - Key fields: `locus_id`, `variant_id`, `CHR`, `BP`, `PIP`, `credible_set_id`, `consequence`, `variant_priority_score`, `candidate_gene`, `gene_score`, `success`, `susie_converged`, `locus_type`, `eqtl_support_flag`.
- `target_gene_table` (`results/database/target_gene_table.tsv`)
  - Gene-centric integrated view for downstream ranking.
  - Key fields: `locus_id`, `gene`, `max_PIP`, `gene_score`, `eqtl_support_flag`, `locus_type`.

### PostgreSQL export tables
Defined by `scripts/integrate_system_layer.py`:
- `loci(locus_id, chromosome, locus_start, locus_end, locus_type, success)`
- `variants(variant_id, locus_id, chr, bp, pip, credible_set_id, consequence, variant_priority_score)`
- `genes(gene_symbol, locus_id, gene_score, eqtl_support_flag, candidate_flag)`

## Example Queries

```sql
-- 1) High-confidence variants across successful loci
SELECT locus_id, variant_id, pip, consequence, variant_priority_score
FROM target_variant_table
WHERE success = TRUE AND pip >= 0.80
ORDER BY pip DESC, variant_priority_score DESC
LIMIT 20;
```

```sql
-- 2) Top candidate genes by score
SELECT locus_id, gene, max_PIP, gene_score, locus_type
FROM target_gene_table
ORDER BY gene_score DESC, max_PIP DESC
LIMIT 20;
```

```sql
-- 3) Loci with converged SuSiE and strong posterior support
SELECT locus_id,
       MAX(pip) AS max_pip,
       COUNT(*) FILTER (WHERE pip >= 0.10) AS n_supporting_variants
FROM target_variant_table
WHERE success = TRUE AND susie_converged = TRUE
GROUP BY locus_id
HAVING MAX(pip) >= 0.50
ORDER BY max_pip DESC;
```

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

- Pipeline architecture: `docs/pipeline_architecture.md`
- Data lineage: `docs/data_lineage.md`
- Architecture narrative: `docs/project_architecture.md`
- File map: `docs/project_file_guide.md`
- Final report: `results/reports/final_project_summary.md`

## Limitations

- Fine-mapping quality depends on LD consistency and conditioning.
- Some loci remain numerically unstable and are flagged in status reports.
- Non-coding loci have gene-assignment uncertainty.
- eQTL/regulatory context is supportive, not definitive causality proof.

## License / Use

This repository is intended for research, reproducible methods demonstration, and portfolio/interview presentation.
