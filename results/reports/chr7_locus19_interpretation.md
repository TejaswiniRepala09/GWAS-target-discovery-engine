# chr7_locus19 Interpretation

## Statistical picture
This is the cleanest locus in the current set. The posterior mass is concentrated in a small group of variants, and each credible set is a singleton.

- Max PIP: `1.0`
- Variants with nonzero PIP: `9`
- Variants with PIP >= 0.1: `9`
- Credible sets: `9`
- Credible set sizes: `[1, 1, 1, 1, 1, 1, 1, 1, 1]`

In short, the fine-mapping signal is focused and statistically strong.

## LD context
The LD diagnostics indicate a stable matrix and successful SuSiE fit (`ld_well_conditioned=TRUE`, `susie_converged=TRUE`), but RSS consistency is not perfect (`rss_consistency=FALSE`). That means the localization is strong, but we should still avoid over-interpreting mechanism from statistics alone.

## Candidate gene interpretation
`PRKAG2` is the top **candidate gene** from this locus based on prioritization and local mapping. It is not confirmed as the causal gene.

Why PRKAG2 is biologically plausible:
- PRKAG2 encodes the gamma-2 subunit of AMPK, a core energy-sensing pathway.
- AMPK signaling is relevant to renal metabolic stress responses.

Important caveat:
- The locus signal is predominantly non-coding in current annotations.
- Because the signal is predominantly non-coding, the true target gene may differ due to distal regulatory effects.

## Variant interpretation (coding vs regulatory)
Current locus annotations are dominated by:
- `intron_variant`
- `non_coding_transcript_variant`

No strong protein-altering lead was identified in this locus under the current consequence filters. So this is currently a regulatory-style signal, not a coding-mechanism signal.

## eQTL / Regulatory Context
Lightweight external checks were performed:

1. **GTEx v8 single-tissue eQTL query (PRKAG2)**
- PRKAG2 has eQTL signal across multiple tissues (`total_eqtl=360` in this query context).
- In this quick check, there was **no Kidney Cortex hit** (`kidney_eqtl=0`).

2. **Ensembl regulatory-overlap query (GRCh37 coordinates for this locus)**
- Regulatory features in-region: `31`
- Dominant feature types: `enhancer`, `CTCF_binding_site`, `promoter`

Interpretation:
- The regulatory annotation density supports a non-coding regulatory mechanism.
- eQTL support exists for PRKAG2 in non-kidney tissues, which keeps PRKAG2 as a plausible candidate.
- This is supportive evidence, not proof of target-gene causality in CKD.

## Working biological hypothesis
If this locus changes local regulation near PRKAG2, it could alter AMPK-related energy signaling and stress handling. That is biologically relevant to kidney function, but the exact gene-tissue-direction chain is still unresolved.

## Confidence
- **Statistical confidence:** high
- **Gene-level mechanistic confidence:** moderate
- **Causal-gene certainty:** not established

## Limitations
- Fine-mapping identifies likely variants; it does not identify causal genes with certainty.
- Non-coding loci are hard to map to one gene.
- eQTL signals are tissue- and context-dependent; kidney-relevant effects may be missing in the available panel.
