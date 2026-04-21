#!/usr/bin/env python3
"""Extract distance-based lead loci from cleaned summary statistics."""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.io_utils import load_summary_stats, load_yaml, write_tsv
from src.loci import extract_distance_based_lead_loci


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    df = load_summary_stats(settings["paths"]["cleaned_summary_stats"])
    lead_variants, lead_loci = extract_distance_based_lead_loci(df, settings)
    write_tsv(lead_variants, settings["paths"]["lead_variants_output"])
    write_tsv(lead_loci, settings["paths"]["lead_loci_output"])
