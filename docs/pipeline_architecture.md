# Pipeline Architecture

This diagram summarizes the execution architecture from GWAS inputs to integrated outputs and benchmarking.

```mermaid
flowchart LR
    A[GWAS Summary Statistics] --> B[QC + Lead Loci Extraction]
    B --> C[Per-Locus Inputs]
    C --> D[VEP Input + Annotation Parsing]
    C --> E[Region-based LD Extraction]
    E --> F[SNP Harmonization + LD Validation]
    F --> G[SuSiE Fine-mapping]
    D --> H[Variant Functional Scoring]
    G --> H
    H --> I[Variant/Gene Prioritization]
    I --> J[DuckDB Integration]
    J --> K[PostgreSQL Export]
    G --> L[FINEMAP Benchmark]
    I --> M[Reports + Plots]
    L --> M
```

## Practical Notes

- Assembly handling is explicit as GRCh37 in the current project workflow.
- VEP is executed through Docker + offline cache mode.
- LD extraction is region-based from 1000 Genomes chromosome VCFs.
- Multi-locus runs use continue-on-error behavior so one failed locus does not stop the full batch.
