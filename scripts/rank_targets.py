#!/usr/bin/env python3
"""Build phase 2 variant and gene prioritization tables."""

from pathlib import Path
import sys

import polars as pl

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.annotation.gene_mapping import build_variant_gene_map
from src.annotation.variant_prioritization import build_gene_priority_table, build_variant_priority_table
from src.io_utils import load_yaml
from src.target_prioritization.ranking import merge_gene_evidence
from src.utils.logging_utils import setup_phase2_logging


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    vep_cfg = load_yaml("config/vep.yaml")
    ot_cfg = load_yaml("config/opentargets.yaml")
    fine_cfg = load_yaml("config/fine_mapping.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    lead_variants = pl.read_csv(settings["paths"]["lead_variants_output"], separator="\t")

    consequence = pl.DataFrame()
    cons_path = Path(vep_cfg["paths"]["variant_consequence_scores_output"])
    if cons_path.exists():
        consequence = pl.read_csv(cons_path, separator="\t")

    pip_summary = pl.DataFrame()
    pip_path = Path(fine_cfg["paths"]["fine_mapping_root"]) / "pip_summary.tsv"
    if pip_path.exists():
        pip_summary = pl.read_csv(pip_path, separator="\t")

    variant_scores = build_variant_priority_table(lead_variants, consequence, pip_summary)

    variant_gene_map = build_variant_gene_map(
        variant_scores,
        consequence,
        settings["phase2"]["references"]["genes_dir"],
        "results/annotations/variant_gene_map.tsv",
    )

    gene_priority = build_gene_priority_table(variant_scores, variant_gene_map)

    ot_path = Path(ot_cfg["paths"]["output_tsv"])
    if ot_path.exists():
        ot_df = pl.read_csv(ot_path, separator="\t")
        gene_priority = merge_gene_evidence(gene_priority, ot_df)

    out_dir = Path("results/target_prioritization")
    out_dir.mkdir(parents=True, exist_ok=True)
    variant_scores.write_csv(out_dir / "variant_priority_scores.tsv", separator="\t")
    gene_priority.write_csv(out_dir / "gene_prioritization.tsv", separator="\t")

    print(f"Variant priority rows: {variant_scores.height}")
    print(f"Gene priority rows: {gene_priority.height}")
