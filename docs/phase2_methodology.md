# Phase 2 Methodology

Phase 2 extends Phase 1 association discovery into target prioritization while preserving explicit uncertainty.

Pipeline flow:
1. Prepare VEP-ready input and parse transcript-aware consequence annotations.
2. Build locus-level fine-mapping bundles and validate LD contracts.
3. Scaffold and parse SuSiE outputs (PIP + credible sets) when available.
4. Build variant-to-gene map (VEP primary, GENCODE nearest fallback).
5. Integrate Open Targets evidence context.
6. Build protein mapping and structure candidate availability checks.
7. Aggregate transparent variant and gene prioritization tables.

## Build and Coordinate Policy

- Required assembly is set in `config/settings.yaml` (`project.assembly`) and must match `phase2.required_assembly`.
- Current repository evidence is GRCh37/b37 (raw source notes, GENCODE file naming, and 1000G Phase3 references).
- The pipeline does not silently liftover coordinates.

## Scientific interpretation guardrails

- Fine-mapping supports statistical prioritization, not mechanism proof.
- Nearest-gene fallback is heuristic.
- Open Targets and structure checks provide translational context, not causal validation.

## Pilot locus note (`chr7_locus10`)

- Pilot locus region: `chr7:65200393-66437208` (GRCh37).
- Observed SuSiE outputs include non-zero PIP variants and two credible sets (`CS1`, `CS2`), with a single-variant CS at PIP 1.0.
- The run produced an IBSS non-convergence warning at 200 iterations. Results are still useful for exploratory prioritization, but should be interpreted with caution.

Recommended production follow-up for this warning:
1. Run `susie_rss` diagnostics (`kriging_rss`, allele/LD consistency checks) on the same locus.
2. Recompute LD from ancestry-matched samples and strict SNP harmonization (including ambiguous allele handling and strand checks).
3. Increase iterations and evaluate residual variance settings only after LD/summary consistency is confirmed.
