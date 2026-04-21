# Data Dictionary

## Input Dataset

File:
- `data/raw/metal_eGFR_meta1.TBL.map.annot.gc.gz`

Source notes:
- `readme.txt` (retained from CKDGen distribution metadata)

## Canonical Columns Used In This Project

- `CHR`: Chromosome
- `BP`: Base-pair position (b37 from source notes)
- `SNP`: Variant identifier (usually `RSID`)
- `A1`: Effect allele
- `A2`: Non-effect allele
- `EAF`: Effect allele frequency
- `BETA`: Effect size
- `SE`: Standard error of effect
- `P`: Uncorrected p-value
- `P_GC`: Genomic-control corrected p-value
- `SE_GC`: Genomic-control corrected standard error
- `N`: Sample size
- `DIRECTION`: Direction of effects across cohorts
- `MAC`: Minor allele count
- `P_ANALYSIS`: Selected p-value used for analysis (configured preference order)

## Observed Source-to-Canonical Mapping

- `chr` or `Chr` -> `CHR`
- `pos` or `Pos` -> `BP`
- `RSID` -> `SNP`
- `Allele1` -> `A1`
- `Allele2` -> `A2`
- `Freq1` -> `EAF`
- `Effect` -> `BETA`
- `StdErr` -> `SE`
- `P-value` -> `P`
- `P.value.GC` -> `P_GC`
- `StdErr.GC` -> `SE_GC`
- `n` -> `N`
- `Direction` -> `DIRECTION`
- `mac` -> `MAC`

## Output Tables

- `results/tables/qc_counts.tsv`: core QC counts
- `results/tables/qc_variants_per_chromosome.tsv`: variant counts by chromosome
- `results/tables/qc_numeric_summary.tsv`: numeric summaries for core quantitative columns
- `results/tables/top_hits.tsv`: strongest associations by analysis p-value
- `results/tables/lead_variants.tsv`: lead variants from distance-based collapse
- `results/tables/lead_loci.tsv`: corresponding lead locus intervals
- `results/tables/vep_input.tsv`: prioritized variants for VEP preparation
