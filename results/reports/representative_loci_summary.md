# Representative Loci Summary

## Why These 3 Loci
- `chr7_locus19`: clean/strong case with compact posterior support.
- `chr7_locus14`: moderate complexity with broader credible-set structure.
- `chr7_locus10`: diffuse/complex case that highlights limitations and caution.

## What Each Locus Teaches
### chr7_locus19
- Max PIP: `1`
- Credible sets: `9` with sizes `[1, 1, 1, 1, 1, 1, 1, 1, 1]`
- LD interpretation: `simple` (`ld_well_conditioned=TRUE`, `rss_consistency=FALSE`)
- Top candidate gene: PRKAG2
- Coding signal: `0` coding variants; structure applicability: `no`
- Confidence statement: `high confidence`

### chr7_locus14
- Max PIP: `1`
- Credible sets: `10` with sizes `[1, 1, 1, 1, 1, 1, 1, 1, 2, 6]`
- LD interpretation: `complex` (`ld_well_conditioned=TRUE`, `rss_consistency=FALSE`)
- Top candidate gene: PPP1R3A
- Coding signal: `0` coding variants; structure applicability: `no`
- Confidence statement: `caution / exploratory`

### chr7_locus10
- Max PIP: `1`
- Credible sets: `10` with sizes `[1, 1, 1, 2, 2, 3, 4, 4, 4, 6]`
- LD interpretation: `complex` (`ld_well_conditioned=TRUE`, `rss_consistency=FALSE`)
- Top candidate genes: GS1-124K5.12, CRCP, TMEM248, AC006001.1, TPST1
- Coding signal: `1` coding variants; structure applicability: `yes`
- Confidence statement: `caution / exploratory`

## Why Together They Tell a Strong Story
- They show a full translational range: clean locus, realistic complexity, and a cautionary diffuse locus.
- This supports honest target discovery communication: prioritization strength, uncertainty, and limitations.

## Gene Assignment Uncertainty
- In these data, fine-mapping identifies likely causal **variants**, not definitive causal genes.
- Gene links come from annotation and prioritization layers, so they should be read as candidate-gene evidence.
- This matters most for non-coding loci, where regulatory effects can act at a distance.
- For this reason, PRKAG2 and PPP1R3A are currently prioritized candidates, not confirmed CKD target genes.

## eQTL / Regulatory Context
- `chr7_locus19` (PRKAG2 candidate): external checks showed broad regulatory feature overlap and cross-tissue eQTL support for PRKAG2, but no clear kidney-cortex hit in the quick GTEx query.
- `chr7_locus14` (PPP1R3A candidate): no strong eQTL/regulatory signal was recovered in the lightweight checks, which increases gene-mapping uncertainty.
- `chr7_locus10` (complex regulatory locus): high density of regulatory features was detected, but no clean single-gene eQTL assignment emerged in this quick pass.
- eQTL evidence is useful for prioritization, but by itself it does not prove causal direction or mechanism.

## Where AlphaFold Fits
- AlphaFold-style structural follow-up is most appropriate when coding variants have protein-change annotations (for example HGVSp/protein_position).
- It is supportive evidence only, not proof of causal mechanism.
- Loci dominated by non-coding/regulatory signatures should prioritize functional genomics before structural claims.

## When Structural Analysis Is Applicable
- Most loci in this batch are non-coding-dominant, so structural modeling is not a universal next step.
- Structural follow-up is best reserved for loci where coding variants carry meaningful posterior support.
- In this representative set, structural interpretation is conditional and exploratory, not central to the strongest signals.
