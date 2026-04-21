# Target Prioritization

## Variant score components

- GWAS significance component (`-log10 p`)
- Consequence severity component (VEP impact + keyword rule)
- SuSiE PIP component (if available)

Output:
- `results/target_prioritization/variant_priority_scores.tsv`

## Gene-level aggregation

Aggregates:
- strongest p-value
- max PIP
- best consequence severity
- supporting locus count
- best variant
- prioritization_reason

Output:
- `results/target_prioritization/gene_prioritization.tsv`

## Open Targets integration

- API: `https://api.platform.opentarget.org/v4/graphql`
- cache: `data/external/opentargets_cache/`
- output: `results/target_prioritization/opentargets_evidence.tsv`
