"""Preprocessing for GWAS summary statistics."""

from __future__ import annotations

import logging
from typing import Any

import polars as pl


def clean_summary_stats(df: pl.DataFrame, settings: dict[str, Any], pvalue_column: str) -> pl.DataFrame:
    """
    Clean summary statistics for downstream GWAS visualization and locus extraction.

    GWAS summary statistics can contain malformed rows or mixed-type columns from large meta-analysis exports.
    This function applies consistent coercion and quality filters while preserving non-core columns.
    """
    analysis_cfg = settings["analysis"]
    selected_p_col = analysis_cfg["selected_pvalue_column"]

    cast_exprs = [
        pl.col("CHR").cast(pl.Int64, strict=False),
        pl.col("BP").cast(pl.Int64, strict=False),
        pl.col(pvalue_column).cast(pl.Float64, strict=False).alias(selected_p_col),
    ]
    for col in ["BETA", "SE", "EAF", "N", "P", "P_GC", "SE_GC", "MAC"]:
        if col in df.columns:
            cast_exprs.append(pl.col(col).cast(pl.Float64, strict=False))

    cleaned = df.with_columns(cast_exprs)
    filter_expr = (
        pl.col("CHR").is_not_null()
        & pl.col("BP").is_not_null()
        & pl.col(selected_p_col).is_not_null()
        & (pl.col("BP") > 0)
        & (pl.col(selected_p_col) > 0.0)
        & (pl.col(selected_p_col) <= 1.0)
    )
    if analysis_cfg.get("autosomes_only", True):
        filter_expr = filter_expr & pl.col("CHR").is_between(
            analysis_cfg["chromosome_min"], analysis_cfg["chromosome_max"]
        )
    cleaned = cleaned.filter(filter_expr)
    logging.info("Rows after cleaning: %s", cleaned.height)
    return cleaned
