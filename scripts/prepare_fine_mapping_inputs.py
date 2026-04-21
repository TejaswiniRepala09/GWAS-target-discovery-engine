#!/usr/bin/env python3
"""Prepare fine-mapping locus bundles and LD resource checks."""

from pathlib import Path
import sys

import polars as pl

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.fine_mapping.fine_mapping_pipeline import run_fine_mapping_preparation
from src.io_utils import load_yaml
from src.utils.logging_utils import setup_phase2_logging


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    fine_cfg = load_yaml("config/fine_mapping.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    cleaned = pl.read_csv(fine_cfg["paths"]["cleaned_summary_stats"], separator="\t")
    lead = pl.read_csv(fine_cfg["paths"]["lead_variants"], separator="\t")
    outputs = run_fine_mapping_preparation(cleaned, lead, fine_cfg)
    print(f"Prepared bundles: {outputs['locus_metadata'].height}")
