# Methodology

## 1) GWAS Summary Statistics Ingestion

The pipeline ingests CKDGen eGFR summary statistics (`.gz` tab-delimited) directly with Polars.

Why summary statistics?
- They provide variant-level association evidence at scale.
- They are portable for reproducible secondary analysis without individual-level genotype data.

The workflow applies explicit column mapping to canonical names (for example, `chr` -> `CHR`, `pos` -> `BP`, `P.value.GC` -> `P_GC`) using config-driven rules.

## 2) Cleaning And Validation

Rows are retained only if they satisfy:
- valid chromosome and base-pair position,
- valid analysis p-value in `(0, 1]`,
- autosomal chromosome filter (`1` to `22`) for default Manhattan/Q-Q analyses.

The pipeline prioritizes corrected p-values (`P_GC`) when available, while preserving both corrected and uncorrected fields.

## 3) QC And Exploratory Outputs

Generated metrics:
- raw and cleaned variant counts,
- per-chromosome counts,
- counts at `p < 5e-8` and `p < 1e-5`,
- top hit table,
- summary statistics for effect-size and uncertainty fields when present.

## 4) GWAS Visualization Logic

### Manhattan Plot
- Uses cumulative genomic coordinates so all chromosomes appear on one continuous x-axis.
- Uses `-log10(p)` on the y-axis to convert tiny p-values into readable magnitudes.
- Includes horizontal reference lines at `5e-8` (genome-wide significance) and `1e-5` (suggestive).

### Q-Q Plot
- Compares observed vs expected `-log10(p)` under null expectation.
- Includes identity line for quick inflation/deviation interpretation.
- Clips extreme p-values to avoid numerical instability in log-transform.

## 5) Distance-Based Lead Locus Extraction

Lead variants are selected from genome-wide significant variants and collapsed by genomic distance threshold (configurable, default 500 kb).

Important caveat:
- This is not LD-based clumping and not fine-mapping.
- Nearby SNPs can be correlated by linkage disequilibrium (LD), so this step is a practical approximation for prioritization only.

## 6) Fine-Mapping Scaffolding (SuSiE Preparation)

Per-locus window files are exported around each lead variant for downstream fine-mapping.

SuSiE is not run here. In a later stage, these locus files should be joined with ancestry-matched LD matrices and passed to a proper fine-mapping model.

## 7) VEP Preparation

A prioritized variant table (`CHR`, `BP`, alleles, IDs where available) is exported for downstream Ensembl VEP workflows.

Exact final VEP input format depends on:
- local vs web VEP mode,
- reference genome build,
- representation requirements (VCF-like, tab format, or HGVS).

## 8) DuckDB Analytics Layer

A lightweight DuckDB database stores key outputs and metadata.

Why this matters:
- simplifies reproducible result tracking,
- enables fast SQL-based comparisons across reruns,
- supports production-style extensions for dashboards, APIs, or cloud execution.
