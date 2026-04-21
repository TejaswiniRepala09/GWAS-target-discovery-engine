"""Pipeline orchestration for CKD target discovery workflow."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import polars as pl

from src.database import write_results_database
from src.export import export_locus_windows, export_vep_input
from src.io_utils import load_summary_stats, load_yaml, setup_logging, write_tsv
from src.loci import extract_distance_based_lead_loci
from src.plotting import (
    plot_beta_distribution,
    plot_effect_vs_significance,
    plot_manhattan,
    plot_per_chromosome_counts,
    plot_qq,
    plot_sample_size_distribution,
)
from src.preprocess import clean_summary_stats
from src.qc import compute_qc_tables
from src.schema import apply_column_mappings, validate_schema


def load_configs(settings_path: str | Path, mappings_path: str | Path) -> tuple[dict[str, Any], dict[str, Any]]:
    return load_yaml(settings_path), load_yaml(mappings_path)


def run_pipeline(settings_path: str | Path = "config/settings.yaml", mappings_path: str | Path = "config/column_mappings.yaml") -> None:
    settings, mappings = load_configs(settings_path, mappings_path)
    setup_logging(settings["paths"]["log_file"])
    raw_df = load_summary_stats(settings["paths"]["input_summary_stats"])
    mapped_df = apply_column_mappings(raw_df, mappings)
    schema_info = validate_schema(mapped_df, mappings)
    cleaned_df = clean_summary_stats(mapped_df, settings, schema_info.pvalue_column)
    write_tsv(cleaned_df, settings["paths"]["cleaned_summary_stats"])

    qc = compute_qc_tables(raw_df, cleaned_df, settings)
    write_tsv(qc["counts"], "results/tables/qc_counts.tsv")
    write_tsv(qc["per_chromosome"], "results/tables/qc_variants_per_chromosome.tsv")
    if qc["numeric_summary"].height:
        write_tsv(qc["numeric_summary"], "results/tables/qc_numeric_summary.tsv")
    write_tsv(qc["top_hits"], settings["paths"]["top_hits_output"])

    plot_manhattan(cleaned_df, settings)
    plot_qq(cleaned_df, settings)
    plot_per_chromosome_counts(qc["per_chromosome"], settings)
    plot_effect_vs_significance(cleaned_df, settings)
    plot_beta_distribution(cleaned_df, settings)
    plot_sample_size_distribution(cleaned_df, settings)

    lead_variants, lead_loci = extract_distance_based_lead_loci(cleaned_df, settings)
    if lead_variants.height:
        write_tsv(lead_variants, settings["paths"]["lead_variants_output"])
    if lead_loci.height:
        write_tsv(lead_loci, settings["paths"]["lead_loci_output"])
    exported_loci = export_locus_windows(cleaned_df, lead_variants, settings) if lead_variants.height else pl.DataFrame()
    export_vep_input(lead_variants if lead_variants.height else qc["top_hits"], settings)

    metadata = pl.DataFrame(
        {
            "metric": ["raw_variant_count", "cleaned_variant_count", "selected_pvalue_column"],
            "value": [str(raw_df.height), str(cleaned_df.height), schema_info.pvalue_column],
        }
    )
    write_results_database(
        settings["paths"]["database_path"],
        metadata,
        lead_variants if lead_variants.height else pl.DataFrame(),
        lead_loci if lead_loci.height else pl.DataFrame(),
        qc["top_hits"],
        exported_loci,
    )
    logging.info("Pipeline completed successfully.")
