#!/usr/bin/env python3
"""Fetch Open Targets evidence for gene prioritization table."""

from pathlib import Path
import sys

import polars as pl

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.io_utils import load_yaml
from src.target_prioritization.opentargets_api import fetch_opentargets_evidence
from src.utils.logging_utils import setup_phase2_logging


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    ot_cfg = load_yaml("config/opentargets.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    gene_path = Path("results/target_prioritization/gene_prioritization.tsv")
    if not gene_path.exists():
        raise FileNotFoundError(f"Missing gene prioritization table: {gene_path}")
    genes = pl.read_csv(gene_path, separator="\t")
    evidence = fetch_opentargets_evidence(genes, ot_cfg)
    out = Path(ot_cfg["paths"]["output_tsv"])
    out.parent.mkdir(parents=True, exist_ok=True)
    evidence.write_csv(out, separator="\t")
    print(f"Open Targets evidence rows: {evidence.height}")
