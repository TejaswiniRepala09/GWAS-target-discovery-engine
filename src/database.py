"""DuckDB integration for storing analysis artifacts."""

from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl


def _connect(db_path: str | Path) -> duckdb.DuckDBPyConnection:
    path = Path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(str(path))


def upsert_table(conn: duckdb.DuckDBPyConnection, table_name: str, df: pl.DataFrame) -> None:
    """Replace table contents with current analysis output."""
    if df.width == 0:
        conn.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT 1 AS placeholder WHERE 1=0")
        return
    conn.register("tmp_df", df.to_pandas())
    conn.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM tmp_df")


def write_results_database(
    db_path: str | Path,
    metadata_df: pl.DataFrame,
    lead_variants_df: pl.DataFrame,
    lead_loci_df: pl.DataFrame,
    top_hits_df: pl.DataFrame,
    exported_locus_files_df: pl.DataFrame,
) -> None:
    """Write core project result tables to DuckDB."""
    conn = _connect(db_path)
    try:
        upsert_table(conn, "cleaned_summary_stats_metadata", metadata_df)
        upsert_table(conn, "lead_variants", lead_variants_df)
        upsert_table(conn, "lead_loci", lead_loci_df)
        upsert_table(conn, "top_hits", top_hits_df)
        upsert_table(conn, "exported_locus_files", exported_locus_files_df)
    finally:
        conn.close()


def write_phase2_tables(db_path: str | Path, tables: dict[str, pl.DataFrame]) -> None:
    """Write phase 2 analysis tables to DuckDB without affecting phase 1 function signatures."""
    conn = _connect(db_path)
    try:
        for table_name, df in tables.items():
            if df is None:
                continue
            upsert_table(conn, table_name, df)
    finally:
        conn.close()
