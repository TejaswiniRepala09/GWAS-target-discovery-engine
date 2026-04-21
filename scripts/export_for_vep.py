#!/usr/bin/env python3
"""Export VEP-ready prioritized variant table."""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.export import export_vep_input
from src.io_utils import load_summary_stats, load_yaml


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    lead = load_summary_stats(settings["paths"]["lead_variants_output"])
    export_vep_input(lead, settings)
