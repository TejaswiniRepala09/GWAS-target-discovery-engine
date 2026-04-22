# Data Lineage

This document tracks how key project artifacts are created and connected.

```mermaid
flowchart TB
    R1[data/raw GWAS input] --> I1[data/interim/ckdgen_egfr_cleaned.tsv]
    I1 --> T1[results/tables/lead_variants.tsv]
    I1 --> T2[results/tables/lead_loci.tsv]
    T1 --> T3[results/tables/vep_input.tsv]
    T2 --> FM1[results/fine_mapping/loci_inputs/<locus_id>/summary_stats.tsv]
    FM1 --> FM2[results/fine_mapping/pip_summary.tsv]
    FM1 --> FM3[results/fine_mapping/credible_sets.tsv]
    T3 --> A1[results/annotations/variant_annotations.tsv]
    A1 --> A2[results/annotations/variant_consequence_scores.tsv]
    FM2 --> P1[results/target_prioritization/variant_priority_scores.tsv]
    A2 --> P1
    P1 --> P2[results/target_prioritization/gene_prioritization.tsv]
    FM2 --> DB1[results/database/target_discovery.duckdb]
    FM3 --> DB1
    P1 --> DB1
    P2 --> DB1
    DB1 --> DB2[results/database/target_variant_table.tsv]
    DB1 --> DB3[results/database/target_gene_table.tsv]
    DB1 --> PG[(PostgreSQL: loci / variants / genes)]
```

## Key lineage checkpoints

- `lead_loci.tsv` defines locus windows.
- `pip_summary.tsv` and `credible_sets.tsv` are the primary fine-mapping outputs.
- `variant_priority_scores.tsv` and `gene_prioritization.tsv` integrate statistical + functional evidence.
- `target_discovery.duckdb` is the analytical source of truth for integrated querying.
