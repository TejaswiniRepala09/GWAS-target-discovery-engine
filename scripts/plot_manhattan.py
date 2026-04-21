#!/usr/bin/env python3
"""Generate Manhattan plot from cleaned summary statistics."""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.io_utils import load_summary_stats, load_yaml
from src.plotting import plot_manhattan


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    df = load_summary_stats(settings["paths"]["cleaned_summary_stats"])
    plot_manhattan(df, settings, output_prefix="manhattan")
