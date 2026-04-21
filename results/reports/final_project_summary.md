# Final Project Summary

## Project Goal
Build a genetics-driven target discovery workflow for CKD that starts from GWAS summary statistics and produces an interpretable shortlist of candidate variants and genes for downstream biological follow-up.

## Data Used
- CKDGen eGFR GWAS summary statistics (Phase 1 and Phase 2 starting point).
- 1000 Genomes Phase 3 LD references (per-chromosome VCF + index).
- Ensembl VEP cache/resources for GRCh37 annotation.
- Optional external evidence layers (Open Targets, UniProt/PDB/AlphaFold metadata checks).

## Pipeline Overview
1. **Phase 1 (signal discovery):** QC, association visualization, top hits, lead variants, lead loci, and per-locus exports.
2. **Phase 2 (target prioritization):**
   - VEP annotation and consequence scoring
   - Per-locus fine-mapping bundle preparation with LD contract checks
   - SuSiE fine-mapping and credible-set summaries
   - Variant-to-gene mapping (annotation-first with fallback logic)
   - Open Targets and structure-prep support layers
   - Transparent variant and gene prioritization
3. **System layer:** DuckDB analytical integration and PostgreSQL relational export.
4. **Benchmarking layer:** SuSiE vs FINEMAP comparison on completed loci.

## Why These Methods Were Chosen
- **SuSiE**: practical Bayesian fine-mapping with interpretable PIP and credible-set outputs.
- **VEP**: standard, transcript-aware consequence annotation for functional interpretation.
- **Polars + DuckDB**: fast, explicit tabular pipelines and local analytical integration.
- **PostgreSQL export**: production-style relational layer for downstream apps/use cases.
- **FINEMAP benchmark**: second fine-mapping method to evaluate robustness and method sensitivity.

## Representative Loci Story
The project intentionally uses a three-locus narrative to show both strength and uncertainty:

- **Strong/clean example:** `chr7_locus19`
  - Highly concentrated posterior signal.
  - Candidate gene interpretation is biologically coherent (PRKAG2 as a candidate, not proven causal).
- **Moderate/complex example:** `chr7_locus14`
  - Signal is interpretable but broader, with more uncertainty in gene assignment.
- **Weak/diffuse example:** `chr7_locus10`
  - Diffuse credible-set structure and higher ambiguity.
  - Useful as an honest limitations case.

This combination makes the project presentation stronger than showing only best-case loci.

## Key Benchmark Findings (SuSiE vs FINEMAP)
Benchmark cohort: **13 successful loci** (from `results/reports/multi_locus_status.tsv` where `success=TRUE`).

- FINEMAP executed for all 13 loci.
- Exact top-variant agreement: **1/13** loci.
- Despite low exact overlap, disagreement is often within the same locus neighborhood:
  - Top-hit distance between methods:
    - median: **24.4 kb**
    - 7/13 loci within **25 kb**
    - 8/13 loci within **50 kb**
    - 4/13 loci > **100 kb**
- FINEMAP runtime (available in benchmark table): median ~**4.34 sec** per locus.

Interpretation:
- In this dataset, method choice materially affects the **single top SNP** call.
- Much of the disagreement looks like **in-locus ambiguity** rather than complete locus-level contradiction.
- Practically, this supports reporting a **credible-variant region** and method-consensus context, not only one winner SNP.

## What Worked Well
- End-to-end reproducibility across multiple loci with per-locus output organization.
- Clear distinction between executed analysis vs scaffolded external-tool contracts.
- Transparent scoring and interpretation tables that are easy to audit.
- System integration (DuckDB + PostgreSQL) that converts analysis outputs into queryable data products.

## What Was Difficult
- Genome-build and cache consistency across GWAS, LD references, and VEP.
- SNP harmonization across GWAS/VCF/LD/fine-mapping formats.
- LD conditioning and numerical stability in some loci (SuSiE failures).
- Method disagreement at top-SNP level during benchmarking.

## Structural Applicability Conclusion
- Most prioritized loci are non-coding dominant.
- Structural follow-up is **conditional**, not universal:
  - most useful when coding variants carry meaningful posterior support and protein-change annotation.
- Structural interpretation is supportive prioritization context, not causal proof.

## Limitations
- Statistical fine-mapping is sensitive to LD quality and model assumptions.
- Several genome-wide loci failed SuSiE due to numerical instability (negative residual variance), highlighting real-world robustness constraints.
- Gene assignment remains uncertain for non-coding signals.
- eQTL/regulatory evidence is helpful but not definitive for causal gene assignment.
- Benchmarking currently focuses on 13 successful loci rather than all lead loci.

## Future Work
1. Harden fine-mapping numerics (LD shrinkage defaults, diagnostics-first rerun policy).
2. Extend FINEMAP/SuSiE benchmarking to larger cohorts with standardized runtime capture for both methods.
3. Add colocalization/eQTL integration as a first-class scoring input.
4. Expand structure branch for high-confidence coding loci with explicit modeling contracts.
5. Add CI checks for key output contracts and schema consistency.

## Final Status
This project is now end-to-end enough for portfolio/recruiter review:
- full association-to-prioritization flow,
- per-locus artifacts and plots,
- system integration (DuckDB + PostgreSQL),
- benchmarking against a second fine-mapping method,
- explicit uncertainty and limitation handling.

## Relevance for Bioinformatics / Pharma Roles
- Demonstrates practical translational genetics workflow design.
- Shows ability to bridge methods, data engineering, and scientific interpretation.
- Emphasizes production habits: modular scripts, robust failure handling, integration layers, and honest evidence language.
