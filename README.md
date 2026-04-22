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

## Nextflow Orchestration (Lightweight, Real Wrapper)

This repository includes a DSL2 Nextflow wrapper that orchestrates existing working scripts without changing scientific logic.

Verified local prerequisites (macOS):
- Java 17+ (`openjdk@17` recommended)
- Nextflow 25.x

Install (verified):

```bash
brew install openjdk@17 nextflow
```

Environment setup (recommended before running Nextflow):

```bash
export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"
export NXF_HOME="$PWD/.nxf_home"
```

Files:
- `nextflow.config`
- `main.nf`
- `modules/pipeline_processes.nf`
- `workflows/phase2_orchestration.nf`

Stage mapping:
1. Locus selection / batch setup -> `scripts/run_phase2_multi_locus.py`
2. VEP annotation -> internal in `run_phase2_multi_locus.py`
3. LD extraction + harmonization -> internal in `run_phase2_multi_locus.py`
4. SuSiE fine-mapping -> internal in `run_phase2_multi_locus.py`
5. Prioritization -> internal in `run_phase2_multi_locus.py`
6. Integration/reporting -> `scripts/integrate_system_layer.py`
7. Benchmarking (optional) -> `scripts/run_benchmark.sh`

Run examples:

```bash
# Lightweight validation (no heavy execution)
nextflow run main.nf -preview --mode top_n --max_loci 1 --locus_ids 'chr7_locus10' --continue_on_error true

# One locus (explicit locus ID)
nextflow run main.nf --mode top_n --max_loci 1 --locus_ids 'chr7_locus10' --continue_on_error true

# Chromosome subset
nextflow run main.nf --mode chromosome --chromosomes '7,8' --max_loci 10 --continue_on_error true

# Top 10 loci genome-wide with benchmark
nextflow run main.nf --mode top_n --max_loci 10 --continue_on_error true --run_benchmark true
```

Validation notes:
- `-preview` validates wiring/config without running heavy pipeline steps.
- In restricted-network environments, Nextflow may show a harmless update-check warning (`Could not resolve host: www.nextflow.io`) while still validating the workflow graph.

## Key Results at a Glance

### 1) Runtime comparison (SuSiE vs FINEMAP)
![SuSiE vs FINEMAP runtime](results/plots/benchmarking/susie_vs_finemap_runtime.png)

### 2) Cross-method benchmark agreement (top variant overlap)
![SuSiE vs FINEMAP top-variant overlap](results/plots/benchmarking/susie_vs_finemap_top_variant_overlap.png)

### 3) Credible-set size comparison
![SuSiE vs FINEMAP credible set sizes](results/plots/benchmarking/susie_vs_finemap_credible_set_sizes.png)

## Pipeline Architecture Diagram

```mermaid
flowchart LR
    A[GWAS Summary Statistics] --> B[QC + Lead Loci]
    B --> C[Per-Locus Input Bundles]
    C --> D[VEP Annotation Parsing]
    C --> E[Region-Based LD Extraction]
    E --> F[SNP Harmonization + LD Validation]
    F --> G[SuSiE Fine-Mapping]
    D --> H[Variant Functional Scoring]
    G --> H
    H --> I[Variant + Gene Prioritization]
    I --> J[DuckDB Integration]
    J --> K[PostgreSQL Export]
    G --> L[FINEMAP Benchmark]
    I --> M[Reports + Plots]
    L --> M
```

## Project Architecture Diagram

```mermaid
flowchart TB
    subgraph Compute["Computation Layer"]
      S1[scripts/run_pipeline.py]
      S2[scripts/run_phase2_multi_locus.py]
      S3[scripts/run_benchmark.sh]
    end

    subgraph Data["Data + Results Layer"]
      D1[results/fine_mapping/*.tsv]
      D2[results/annotations/*.tsv]
      D3[results/target_prioritization/*.tsv]
      D4[results/reports/*.md]
      D5[results/plots/benchmarking/*.png]
    end

    subgraph System["System Integration Layer"]
      DB1[(DuckDB: target_discovery.duckdb)]
      DB2[(PostgreSQL: loci / variants / genes)]
    end

    S1 --> D1
    S2 --> D1
    S2 --> D2
    S2 --> D3
    S2 --> D4
    S3 --> D5
    D1 --> DB1
    D2 --> DB1
    D3 --> DB1
    DB1 --> DB2
```

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

## Colocalization / Target-Gene Confidence Layer

To strengthen non-coding locus interpretation, this repo now includes a coloc-ready representative-loci layer:

```bash
python scripts/run_coloc_representative_loci.py
```

Outputs:
- `results/coloc/<locus_id>/gwas_coloc_input.tsv`
- `results/coloc/<locus_id>/eqtl_input_template.tsv`
- `results/coloc/coloc_status.tsv`
- `results/reports/coloc_summary.md`

Current status:
- Full coloc posterior inference is not run by default because matched eQTL summary statistics are not yet committed.
- The layer is designed to be honest and reproducible: it prepares exact GWAS-side inputs and documents what is missing for production-grade shared-signal testing.

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

## Testing and CI

This repo includes lightweight reliability checks on every push/PR via GitHub Actions:

- Python syntax validation for `scripts/` and `src/`
- Output contract validation for key committed TSV artifacts

Contract checks verify:
- file exists
- TSV is readable
- expected columns are present
- key columns are not entirely null/empty

Checked files:
- `results/fine_mapping/pip_summary.tsv`
- `results/fine_mapping/credible_sets.tsv`
- `results/target_prioritization/variant_priority_scores.tsv`
- `results/target_prioritization/gene_prioritization.tsv`
- `results/reports/multi_locus_status.tsv`
- `results/benchmarking/susie_vs_finemap_locus_summary.tsv`

Run locally:

```bash
python scripts/check_output_contracts.py
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

Folder guides:
- `docs/README.md`
- `scripts/README.md`
- `results/README.md`

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
