"""Credible set post-processing helpers."""

from __future__ import annotations

import polars as pl


def normalize_credible_sets(credible_sets: pl.DataFrame, pip_summary: pl.DataFrame) -> pl.DataFrame:
    """Attach rank and optional PIP to credible set members."""
    if credible_sets.height == 0:
        return credible_sets

    cs = credible_sets
    if "PIP" in pip_summary.columns and "SNP" in pip_summary.columns:
        cs = cs.join(pip_summary.select([c for c in ["locus_id", "SNP", "PIP"] if c in pip_summary.columns]), on=[c for c in ["locus_id", "SNP"] if c in cs.columns], how="left")
    if "PIP" in cs.columns:
        cs = cs.with_columns(pl.col("PIP").rank(method="ordinal", descending=True).over([c for c in ["locus_id", "credible_set_id"] if c in cs.columns]).alias("posterior_rank"))
    return cs
