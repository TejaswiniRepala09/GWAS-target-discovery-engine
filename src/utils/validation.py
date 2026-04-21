"""Validation helpers for files, schema, and assembly consistency."""

from __future__ import annotations

from pathlib import Path

import polars as pl


def require_columns(df: pl.DataFrame, columns: list[str], context: str) -> None:
    """Raise readable error if required columns are missing."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {context}: {missing}")


def optional_file(path: str | Path) -> Path | None:
    """Return path if exists, else None."""
    p = Path(path)
    return p if p.exists() else None


def require_file(path: str | Path, context: str) -> Path:
    """Require file exists."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing required file for {context}: {p}")
    return p


def assert_assembly(expected: str, observed: str | None, context: str) -> None:
    """Validate assembly/build where observed metadata exists."""
    if observed is None:
        raise ValueError(
            f"Assembly for {context} could not be automatically verified. "
            f"Set explicit assembly in config and confirm coordinates manually. Expected: {expected}"
        )
    if expected.upper() != observed.upper():
        raise ValueError(f"Assembly mismatch for {context}: expected {expected}, observed {observed}")
