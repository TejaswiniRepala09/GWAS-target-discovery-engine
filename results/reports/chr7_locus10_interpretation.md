# chr7_locus10 Interpretation

## Statistical picture
This is the messiest of the representative loci. It has strong peaks, but the posterior mass is spread over many variants and larger credible sets.

- Max PIP: `1.0`
- Variants with nonzero PIP: `81`
- Variants with PIP >= 0.1: `23`
- Credible sets: `10`
- Credible set sizes: `[1, 1, 1, 2, 2, 3, 4, 4, 4, 6]`

This is a realistic example of a locus where fine-mapping helps, but uncertainty remains substantial.

## LD context
The LD matrix passed conditioning checks (`ld_well_conditioned=TRUE`) and SuSiE converged, but RSS consistency remains imperfect (`rss_consistency=FALSE`). That pattern is consistent with a complex regional architecture and/or summary-LD mismatch effects.

## Candidate gene interpretation
Top **candidate genes** from prioritization include:
- `GS1-124K5.12`
- `CRCP`
- `TMEM248`
- `AC006001.1`
- `TPST1`

These are candidate assignments, not confirmed causal genes.

Important caveat:
- Most support is non-coding/regulatory.
- Because the signal is predominantly non-coding, the true target gene may differ due to distal regulatory effects.

## Variant interpretation (coding vs regulatory)
This locus includes coding annotations, but they are not driving the fine-mapped signal under current posterior support.

- Coding variant count (selected classes): `1`
- Max PIP among coding variants: `0`
- Example missense annotations:
  - `rs9530` (GUSB), missense, PIP `0.0`
  - `rs17138089` (CRCP), missense, PIP missing in matched set

So, this is currently better interpreted as a regulatory-style locus than a coding-mechanism locus.

## eQTL / Regulatory Context
External quick checks:

1. **GTEx v8 single-tissue eQTL checks**
- No strong direct eQTL signal was captured in the lightweight lookup for the top candidate genes tested from this locus.

2. **Ensembl regulatory-overlap query (GRCh37 coordinates for this locus)**
- Regulatory features in-region: `164`
- Dominant feature types: `enhancer`, `CTCF_binding_site`, `promoter`, `open_chromatin_region`

Interpretation:
- The heavy regulatory-feature density supports a non-coding regulatory mechanism.
- eQTL-to-gene assignment is still unresolved in this quick pass, so gene mapping remains uncertain.

## Working biological hypothesis
Variants in this region may alter regulatory elements that control one or more nearby genes involved in kidney-relevant biology. The locus is biologically interesting, but target-gene resolution is still incomplete.

## Confidence
- **Statistical confidence:** moderate (with complexity)
- **Gene-level mechanistic confidence:** low
- **Causal-gene certainty:** not established

## Limitations
- Fine-mapping narrows variants but does not prove mechanism.
- Coding annotations here do not currently align with high posterior support.
- Regulatory and eQTL effects may be context-specific and not fully captured in current public resources.
