"""Merge evidence into final gene prioritization."""

from __future__ import annotations

import polars as pl


def merge_gene_evidence(gene_priority: pl.DataFrame, opentargets_evidence: pl.DataFrame) -> pl.DataFrame:
    """Merge Open Targets evidence and update final score transparently."""
    if gene_priority.height == 0:
        return gene_priority

    merged = gene_priority
    if opentargets_evidence.height:
        merged = gene_priority.join(
            opentargets_evidence,
            on=[c for c in ["gene_symbol", "ensembl_gene_id"] if c in gene_priority.columns and c in opentargets_evidence.columns],
            how="left",
        )

    merged = merged.with_columns(
        pl.col("opentargets_support").cast(pl.Float64, strict=False).fill_null(0.0).alias("component_opentargets")
        if "opentargets_support" in merged.columns else pl.lit(0.0).alias("component_opentargets")
    ).with_columns(
        (pl.col("gene_prioritization_score") + pl.col("component_opentargets") * 0.15).alias("gene_prioritization_score_final")
    ).with_columns(
        (
            pl.col("prioritization_reason")
            + pl.lit("; opentargets_support=")
            + pl.col("component_opentargets").cast(pl.Utf8)
        ).alias("prioritization_reason")
    )

    return merged.sort("gene_prioritization_score_final", descending=True)
