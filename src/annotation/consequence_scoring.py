"""Consequence severity scoring with explicit transparent rules."""

from __future__ import annotations

from typing import Any

import polars as pl


def score_consequences(annotation_df: pl.DataFrame, vep_cfg: dict[str, Any]) -> pl.DataFrame:
    """Attach consequence severity component scores."""
    if annotation_df.height == 0:
        return annotation_df

    impact_weights = vep_cfg.get("consequence_weights", {})
    keyword_scores = vep_cfg.get("consequence_keyword_scores", {})

    cons = pl.col("Consequence_term").cast(pl.Utf8).fill_null("")

    score_expr = pl.lit(0.0)
    for keyword, score in keyword_scores.items():
        score_expr = pl.when(cons.str.contains(keyword, literal=True)).then(pl.lit(float(score))).otherwise(score_expr)

    scored = annotation_df.with_columns(
        pl.col("IMPACT").cast(pl.Utf8).replace(impact_weights, default=0.0).cast(pl.Float64).alias("impact_severity_score"),
        score_expr.alias("consequence_keyword_score"),
    ).with_columns(
        (pl.col("impact_severity_score") * 0.6 + pl.col("consequence_keyword_score") * 0.4).alias("consequence_severity_score")
    )

    return scored
