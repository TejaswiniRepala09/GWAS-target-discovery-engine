# chr7_locus14 Interpretation

## Statistical picture
This locus is interpretable, but it is not as clean as the strongest signals. Posterior support is spread across more variants and set sizes.

- Max PIP: `1.0`
- Variants with nonzero PIP: `35`
- Variants with PIP >= 0.1: `15`
- Credible sets: `10`
- Credible set sizes: `[1, 1, 1, 1, 1, 1, 1, 1, 2, 6]`

This looks like a moderate-to-complex locus rather than a single-variant story.

## LD context
LD diagnostics are usable (`ld_well_conditioned=TRUE`, `susie_converged=TRUE`) but not perfect (`rss_consistency=FALSE`). So the fine-mapping output is informative, but we should treat mechanistic interpretation with caution.

## Candidate gene interpretation
`PPP1R3A` is the top **candidate gene** from prioritization at this locus. It should not be treated as a confirmed causal gene.

Biological context:
- PPP1R3A is involved in glycogen-related metabolic regulation.
- Metabolic signaling is relevant to CKD risk biology, especially in energy-stress settings.

Important caveat:
- This locus is predominantly non-coding in the current annotation set.
- Because the signal is predominantly non-coding, the true target gene may differ due to distal regulatory effects.

## Variant interpretation (coding vs regulatory)
The leading variants are non-coding under current filters, and no high-confidence coding mechanism emerges here.

- Coding variant count (selected classes): `0`
- Max PIP among coding variants: `NA`
- Structural follow-up suitability: currently low

## eQTL / Regulatory Context
External quick checks:

1. **GTEx v8 single-tissue eQTL query (PPP1R3A)**
- No robust signal was captured in this lightweight lookup for PPP1R3A.

2. **Ensembl regulatory-overlap query (GRCh37 coordinates for this locus)**
- No strong regulatory-overlap signal was returned by this quick endpoint query for the tested window.

Interpretation:
- At this point, there is no strong direct eQTL/regulatory support that cleanly links this locus to one gene.
- That adds uncertainty to gene assignment and keeps this locus in the exploratory bucket.

## Working biological hypothesis
A plausible model is that non-coding variants in this region perturb local regulatory control of metabolic genes, with PPP1R3A as one candidate. The data support prioritization, not confirmation.

## Confidence
- **Statistical confidence:** moderate
- **Gene-level mechanistic confidence:** low-to-moderate
- **Causal-gene certainty:** not established

## Limitations
- Fine-mapping prioritizes variants, not definitive target genes.
- Non-coding loci are difficult to map to one gene with confidence.
- eQTL resources may miss disease-context or kidney-state specific effects.
