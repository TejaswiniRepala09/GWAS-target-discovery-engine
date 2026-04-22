# Colocalization Summary (Representative Loci)

This layer prepares **coloc-ready inputs** and records missing requirements for full coloc.

## What was run
- GWAS-side coloc input preparation for representative loci.
- Candidate-gene extraction from existing prioritization outputs.
- Missing-data manifest generation for production coloc.

## Why full coloc was not run
- The repository currently does not contain matched eQTL summary statistics per representative locus and tissue.
- Without matched eQTL effect/SE/p-value data on overlapping SNPs, coloc posterior probabilities (PP.H0-H4) cannot be estimated reliably.

## Locus-level status

| locus_id | candidate_genes | coloc_run | shared_signal_support | confidence_interpretation |
|---|---|---|---|---|
| chr7_locus19 | NA | False | weak_to_moderate_indirect_support | Coloc not run. Existing project notes suggest broader cross-tissue regulatory support for PRKAG2, but no definitive kidney-specific shared-signal test yet. |
| chr7_locus14 | NA | False | inconclusive | Coloc not run. Current evidence is predominantly non-coding and gene assignment remains uncertain without matched eQTL summary statistics. |
| chr7_locus10 | NA | False | inconclusive | Coloc not run. Locus remains complex/diffuse; no matched eQTL summary statistics available to test shared causal signal. |

## Production-grade coloc requirements
For each locus and candidate gene/tissue, provide eQTL summary statistics with:
- same genome build (GRCh37)
- SNP-level fields: SNP, CHR, BP, effect allele, other allele, beta, se, p, n, maf
- tissue and gene metadata
- robust harmonization to GWAS alleles and SNP overlap.

Generated files:
- `results/coloc/<locus_id>/gwas_coloc_input.tsv`
- `results/coloc/<locus_id>/eqtl_input_template.tsv`
- `results/coloc/coloc_status.tsv`
