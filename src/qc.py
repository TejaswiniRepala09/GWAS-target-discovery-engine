"""GWAS QC/exploration utilities."""

from __future__ import annotations

from typing import Any

import polars as pl


def compute_qc_tables(
    raw_df: pl.DataFrame, cleaned_df: pl.DataFrame, settings: dict[str, Any]
) -> dict[str, pl.DataFrame]:
    """Compute core QC tables for summary reporting."""
    analysis_cfg = settings["analysis"]
    p_col = analysis_cfg["selected_pvalue_column"]

    counts = pl.DataFrame(
        {
            "metric": [
                "total_variants_raw",
                "total_variants_cleaned",
                "genome_wide_significant_variants",
                "suggestive_variants",
            ],
            "value": [
                raw_df.height,
                cleaned_df.height,
                cleaned_df.filter(pl.col(p_col) < analysis_cfg["genome_wide_significance_threshold"]).height,
                cleaned_df.filter(pl.col(p_col) < analysis_cfg["suggestive_significance_threshold"]).height,
            ],
        }
    )
    per_chr = cleaned_df.group_by("CHR").len().rename({"len": "variant_count"}).sort("CHR")
    top_hits = cleaned_df.sort(p_col).head(analysis_cfg["top_hits_n"])

    summary_tables: list[pl.DataFrame] = []
    for col in ["BETA", "SE", "EAF", "N"]:
        if col in cleaned_df.columns:
            summary_tables.append(
                cleaned_df.select(
                    pl.lit(col).alias("variable"),
                    pl.col(col).mean().alias("mean"),
                    pl.col(col).std().alias("std"),
                    pl.col(col).min().alias("min"),
                    pl.col(col).max().alias("max"),
                )
            )
    numeric_summary = pl.concat(summary_tables) if summary_tables else pl.DataFrame()

    return {"counts": counts, "per_chromosome": per_chr, "top_hits": top_hits, "numeric_summary": numeric_summary}
