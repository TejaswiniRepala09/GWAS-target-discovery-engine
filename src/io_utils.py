"""I/O utilities for CKD GWAS summary statistics workflow."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import polars as pl
import yaml


def setup_logging(log_file: str | Path | None = None, level: int = logging.INFO) -> None:
    """Configure console and optional file logging."""
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path))
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML file into a dictionary."""
    yaml_path = Path(path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")
    with yaml_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def ensure_path_exists(path: str | Path) -> Path:
    """Raise a clear error if a required path is missing."""
    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Required input file does not exist: {resolved}")
    return resolved


def load_summary_stats(path: str | Path) -> pl.DataFrame:
    """Read a gzipped or plain TSV summary statistics file with Polars."""
    file_path = ensure_path_exists(path)
    logging.info("Loading summary statistics: %s", file_path)
    numeric_overrides = {
        "n": pl.Float64,
        "N": pl.Float64,
        "mac": pl.Float64,
        "Freq1": pl.Float64,
        "Effect": pl.Float64,
        "StdErr": pl.Float64,
        "P-value": pl.Float64,
        "P.value": pl.Float64,
        "P.value.GC": pl.Float64,
        "StdErr.GC": pl.Float64,
        "chr": pl.Int64,
        "pos": pl.Int64,
    }
    df = pl.read_csv(
        file_path,
        separator="\t",
        infer_schema_length=10000,
        null_values=["NA", ""],
        schema_overrides=numeric_overrides,
    )
    logging.info("Loaded %s rows and %s columns", df.height, df.width)
    logging.info("Detected columns: %s", ", ".join(df.columns))
    return df


def write_tsv(df: pl.DataFrame, path: str | Path) -> None:
    """Write dataframe as tab-delimited text."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(output_path, separator="\t")
    logging.info("Wrote table: %s (%s rows)", output_path, df.height)
