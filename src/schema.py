"""Schema validation and column normalization utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import polars as pl


@dataclass(frozen=True)
class SchemaInfo:
    """Container describing resolved canonical columns."""

    required_present: list[str]
    optional_present: list[str]
    pvalue_column: str


def apply_column_mappings(df: pl.DataFrame, mappings_cfg: dict[str, Any]) -> pl.DataFrame:
    """Rename source columns to canonical names using configuration mappings."""
    mappings: dict[str, str] = mappings_cfg.get("canonical_mappings", {})
    rename_map = {source: target for source, target in mappings.items() if source in df.columns}
    if rename_map:
        df = df.rename(rename_map)
    return df


def choose_pvalue_column(df: pl.DataFrame, mappings_cfg: dict[str, Any]) -> str:
    """Select p-value column by configured priority."""
    priority: list[str] = mappings_cfg.get("pvalue_priority", ["P_GC", "P"])
    for col in priority:
        if col in df.columns:
            return col
    raise ValueError(f"None of configured p-value columns were found: {priority}")


def validate_schema(df: pl.DataFrame, mappings_cfg: dict[str, Any]) -> SchemaInfo:
    """Validate required columns and report optional ones available."""
    required_cols: list[str] = mappings_cfg.get("required_columns", [])
    optional_cols: list[str] = mappings_cfg.get("optional_columns", [])
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns after mapping: {missing}. Found: {df.columns}")
    pval_col = choose_pvalue_column(df, mappings_cfg)
    return SchemaInfo(
        required_present=[col for col in required_cols if col in df.columns],
        optional_present=[col for col in optional_cols if col in df.columns],
        pvalue_column=pval_col,
    )
