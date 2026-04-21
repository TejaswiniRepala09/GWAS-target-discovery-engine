"""AlphaFold follow-up annotations."""

from __future__ import annotations

import polars as pl


def annotate_alphafold_followup(df: pl.DataFrame) -> pl.DataFrame:
    """Add explicit non-causal structural follow-up language."""
    if df.height == 0:
        return df
    return df.with_columns(
        pl.lit("Candidate structural impact assessment; modeling/experimental validation required.").alias("structural_followup_note")
    )
