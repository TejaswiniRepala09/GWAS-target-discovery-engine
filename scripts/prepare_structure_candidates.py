#!/usr/bin/env python3
"""Build protein mapping and structure-candidate summary tables."""

from pathlib import Path
import sys

import polars as pl

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.io_utils import load_yaml
from src.structure.protein_mapping import build_protein_mapping
from src.structure.structural_summary import build_structure_candidates
from src.utils.logging_utils import setup_phase2_logging


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    struct_cfg = load_yaml("config/structure.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    var_path = Path("results/target_prioritization/variant_priority_scores.tsv")
    gene_path = Path("results/target_prioritization/gene_prioritization.tsv")
    if not var_path.exists() or not gene_path.exists():
        raise FileNotFoundError("Missing variant_priority_scores.tsv or gene_prioritization.tsv")

    variants = pl.read_csv(var_path, separator="\t")
    genes = pl.read_csv(gene_path, separator="\t")

    protein_map = build_protein_mapping(genes, struct_cfg)
    struct = build_structure_candidates(variants, protein_map, struct_cfg)

    Path(struct_cfg["paths"]["protein_mapping_output"]).parent.mkdir(parents=True, exist_ok=True)
    protein_map.write_csv(struct_cfg["paths"]["protein_mapping_output"], separator="\t")
    struct["structure_candidates"].write_csv(struct_cfg["paths"]["structure_candidates_output"], separator="\t")
    struct["structural_summary"].write_csv(struct_cfg["paths"]["structural_summary_output"], separator="\t")
    print(f"Protein mappings: {protein_map.height}")
    print(f"Structure candidates: {struct['structure_candidates'].height}")
