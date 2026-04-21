# Structural Interpretation Layer

This phase performs lightweight mapping and availability checks, not full modeling.

## Data sources

- UniProt REST (`https://rest.uniprot.org`)
- PDB core entry API (`https://data.rcsb.org/rest/v1/core/entry/`)
- AlphaFold API (`https://alphafold.ebi.ac.uk/api/`)

## Outputs

- `results/structure/protein_mapping.tsv`
- `results/structure/structure_candidates.tsv`
- `results/structure/structural_interpretation_summary.tsv`

## Interpretation

These are **candidate structural follow-up** outputs only.
No stability deltas or causal claims are inferred by this layer.
