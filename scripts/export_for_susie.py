#!/usr/bin/env python3
"""Export per-locus windows for future SuSiE fine-mapping."""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.export import export_locus_windows
from src.io_utils import load_summary_stats, load_yaml


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    cleaned = load_summary_stats(settings["paths"]["cleaned_summary_stats"])
    lead = load_summary_stats(settings["paths"]["lead_variants_output"])
    export_locus_windows(cleaned, lead, settings)
