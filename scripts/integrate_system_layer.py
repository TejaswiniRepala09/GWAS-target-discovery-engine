#!/usr/bin/env python3
"""Build system integration artifacts from existing Phase 2 outputs.

This script does not rerun analysis. It only integrates existing TSV outputs into:
1) DuckDB analytical layer
2) Optional PostgreSQL export layer
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path
from typing import Iterable

import duckdb
import polars as pl


ROOT = Path(__file__).resolve().parent.parent


def _print_step(cmd: str) -> None:
    print(f"$ {cmd}")


def _ensure_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input missing: {path}")


def _load_tsv_table(conn: duckdb.DuckDBPyConnection, table: str, path: Path) -> None:
    _ensure_file(path)
    conn.execute(
        f"""
        CREATE OR REPLACE TABLE {table} AS
        SELECT * FROM read_csv_auto(?, delim='\\t', header=true, sample_size=-1, ignore_errors=true);
        """,
        [str(path)],
    )


def _read_per_locus_tables(pattern: str, locus_from_parent: bool = False) -> pl.DataFrame:
    paths = sorted(glob.glob(pattern))
    frames: list[pl.DataFrame] = []
    for p in paths:
        path = Path(p)
        locus_id = path.parent.name
        if locus_id.endswith("_regularized") or locus_id.endswith("_corrected"):
            continue
        df = pl.read_csv(path, separator="\t")
        if locus_from_parent or "locus_id" not in df.columns:
            df = df.with_columns(pl.lit(locus_id).alias("locus_id"))
        frames.append(df)
    if not frames:
        return pl.DataFrame()
    return pl.concat(frames, how="diagonal_relaxed")


def _unique_locus_count(conn: duckdb.DuckDBPyConnection, table: str) -> int:
    return int(conn.execute(f"SELECT COUNT(DISTINCT CAST(locus_id AS VARCHAR)) FROM {table};").fetchone()[0] or 0)


def _replace_table_from_polars(conn: duckdb.DuckDBPyConnection, table: str, df: pl.DataFrame) -> None:
    conn.register("tmp_polars_df", df.to_pandas())
    conn.execute(f"CREATE OR REPLACE TABLE {table} AS SELECT * FROM tmp_polars_df;")


def _create_duckdb_tables(db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(db_path))
    try:
        _load_tsv_table(conn, "pip_summary", ROOT / "results/fine_mapping/pip_summary.tsv")
        _load_tsv_table(conn, "credible_sets", ROOT / "results/fine_mapping/credible_sets.tsv")
        _load_tsv_table(conn, "variant_priority", ROOT / "results/target_prioritization/variant_priority_scores.tsv")
        _load_tsv_table(conn, "gene_prioritization", ROOT / "results/target_prioritization/gene_prioritization.tsv")
        _load_tsv_table(conn, "variant_annotations", ROOT / "results/annotations/variant_annotations.tsv")
        _load_tsv_table(conn, "locus_status", ROOT / "results/reports/multi_locus_status.tsv")

        # If aggregate TSVs are stale, rebuild from per-locus outputs.
        per_locus_variant_priority = _read_per_locus_tables(
            str(ROOT / "results/target_prioritization/chr*_locus*/variant_priority_scores.tsv")
        )
        if per_locus_variant_priority.height:
            if _unique_locus_count(conn, "variant_priority") < per_locus_variant_priority.select(pl.col("locus_id").cast(pl.Utf8).n_unique()).item():
                _replace_table_from_polars(conn, "variant_priority", per_locus_variant_priority)

        per_locus_gene_priority = _read_per_locus_tables(
            str(ROOT / "results/target_prioritization/chr*_locus*/gene_prioritization.tsv")
        )
        if per_locus_gene_priority.height:
            if _unique_locus_count(conn, "gene_prioritization") < per_locus_gene_priority.select(pl.col("locus_id").cast(pl.Utf8).n_unique()).item():
                _replace_table_from_polars(conn, "gene_prioritization", per_locus_gene_priority)

        per_locus_pip = _read_per_locus_tables(
            str(ROOT / "results/fine_mapping/susie_raw/chr*_locus*/pip.tsv"),
            locus_from_parent=True,
        )
        if per_locus_pip.height:
            if _unique_locus_count(conn, "pip_summary") < per_locus_pip.select(pl.col("locus_id").cast(pl.Utf8).n_unique()).item():
                _replace_table_from_polars(conn, "pip_summary", per_locus_pip)

        per_locus_cs = _read_per_locus_tables(
            str(ROOT / "results/fine_mapping/susie_raw/chr*_locus*/credible_sets.tsv"),
            locus_from_parent=True,
        )
        if per_locus_cs.height:
            if _unique_locus_count(conn, "credible_sets") < per_locus_cs.select(pl.col("locus_id").cast(pl.Utf8).n_unique()).item():
                _replace_table_from_polars(conn, "credible_sets", per_locus_cs)

        # Optional enrichment sources for derived fields.
        ot_path = ROOT / "results/target_prioritization/opentargets_evidence.tsv"
        if ot_path.exists():
            _load_tsv_table(conn, "opentargets_evidence", ot_path)
        else:
            conn.execute(
                """
                CREATE OR REPLACE TABLE opentargets_evidence AS
                SELECT NULL::VARCHAR AS gene_symbol, NULL::VARCHAR AS ensembl_gene_id,
                       NULL::VARCHAR AS opentargets_status, 0.0::DOUBLE AS opentargets_support,
                       NULL::VARCHAR AS opentargets_note
                WHERE 1=0;
                """
            )

        lead_loci_path = ROOT / "results/tables/lead_loci.tsv"
        if lead_loci_path.exists():
            _load_tsv_table(conn, "lead_loci", lead_loci_path)
        else:
            conn.execute(
                """
                CREATE OR REPLACE TABLE lead_loci AS
                SELECT NULL::VARCHAR AS locus_id, NULL::INTEGER AS CHR,
                       NULL::BIGINT AS locus_start, NULL::BIGINT AS locus_end
                WHERE 1=0;
                """
            )

        conn.execute(
            """
            CREATE OR REPLACE TEMP TABLE locus_signal AS
            SELECT
                locus_id,
                MAX(CAST(PIP AS DOUBLE)) AS max_pip,
                CASE
                    WHEN MAX(CAST(PIP AS DOUBLE)) >= 0.9 THEN 'strong'
                    WHEN MAX(CAST(PIP AS DOUBLE)) >= 0.3 THEN 'moderate'
                    ELSE 'weak'
                END AS locus_type
            FROM pip_summary
            GROUP BY locus_id;
            """
        )

        conn.execute(
            """
            CREATE OR REPLACE TEMP TABLE va_variant_summary AS
            SELECT
                CAST(Uploaded_variation AS VARCHAR) AS Uploaded_variation,
                string_agg(DISTINCT CAST(Consequence AS VARCHAR), '|' ORDER BY CAST(Consequence AS VARCHAR)) AS consequence,
                MIN(NULLIF(CAST(SYMBOL AS VARCHAR), '')) AS symbol_any,
                MIN(NULLIF(CAST(Gene AS VARCHAR), '')) AS gene_any
            FROM variant_annotations
            GROUP BY 1;
            """
        )

        conn.execute(
            """
            CREATE OR REPLACE TEMP TABLE var_base AS
            SELECT DISTINCT
                vp.locus_id,
                COALESCE(NULLIF(CAST(vp.SNP AS VARCHAR), ''), NULLIF(CAST(vp.Uploaded_variation AS VARCHAR), ''), NULLIF(CAST(vp.MARKERNAME AS VARCHAR), '')) AS variant_id,
                TRY_CAST(vp.CHR AS INTEGER) AS CHR,
                TRY_CAST(vp.BP AS BIGINT) AS BP,
                CAST(vp.variant_priority_score AS DOUBLE) AS variant_priority_score,
                COALESCE(NULLIF(CAST(vp.gene_symbol AS VARCHAR), ''), vas.symbol_any, vas.gene_any) AS candidate_gene,
                NULLIF(CAST(vp.ensembl_gene_id AS VARCHAR), '') AS ensembl_gene_id,
                vas.consequence
            FROM variant_priority vp
            LEFT JOIN va_variant_summary vas
                ON CAST(vp.Uploaded_variation AS VARCHAR) = vas.Uploaded_variation;
            """
        )

        conn.execute(
            """
            CREATE OR REPLACE TEMP TABLE ot_eqtl AS
            SELECT
                COALESCE(NULLIF(CAST(gene_symbol AS VARCHAR), ''), NULLIF(CAST(ensembl_gene_id AS VARCHAR), '')) AS gene_key,
                MAX(CASE WHEN TRY_CAST(opentargets_support AS DOUBLE) > 0 THEN TRUE ELSE FALSE END) AS eqtl_support_flag
            FROM opentargets_evidence
            GROUP BY 1;
            """
        )

        conn.execute(
            """
            CREATE OR REPLACE TABLE target_variant_table AS
            WITH pip_norm AS (
                SELECT CAST(locus_id AS VARCHAR) AS locus_id, CAST(SNP AS VARCHAR) AS variant_id, CAST(PIP AS DOUBLE) AS PIP
                FROM pip_summary
            ),
            cs_norm AS (
                SELECT CAST(locus_id AS VARCHAR) AS locus_id, CAST(SNP AS VARCHAR) AS variant_id, CAST(credible_set_id AS VARCHAR) AS credible_set_id
                FROM credible_sets
            )
            SELECT
                vb.locus_id,
                vb.variant_id,
                vb.CHR,
                vb.BP,
                p.PIP,
                cs.credible_set_id,
                vb.consequence,
                vb.variant_priority_score,
                vb.candidate_gene,
                CAST(gp.gene_prioritization_score AS DOUBLE) AS gene_score,
                CAST(ls.success AS BOOLEAN) AS success,
                CAST(ls.susie_converged AS BOOLEAN) AS susie_converged,
                lsg.locus_type,
                COALESCE(ote.eqtl_support_flag, FALSE) AS eqtl_support_flag
            FROM var_base vb
            LEFT JOIN pip_norm p ON vb.locus_id = p.locus_id AND vb.variant_id = p.variant_id
            LEFT JOIN cs_norm cs ON vb.locus_id = cs.locus_id AND vb.variant_id = cs.variant_id
            LEFT JOIN gene_prioritization gp
                   ON vb.locus_id = CAST(gp.locus_id AS VARCHAR)
                  AND vb.candidate_gene = COALESCE(NULLIF(CAST(gp.gene_symbol AS VARCHAR), ''), NULLIF(CAST(gp.ensembl_gene_id AS VARCHAR), ''))
            LEFT JOIN locus_status ls ON vb.locus_id = CAST(ls.locus_id AS VARCHAR)
            LEFT JOIN locus_signal lsg ON vb.locus_id = lsg.locus_id
            LEFT JOIN ot_eqtl ote ON COALESCE(vb.candidate_gene, vb.ensembl_gene_id) = ote.gene_key;
            """
        )

        conn.execute(
            """
            CREATE OR REPLACE TABLE target_gene_table AS
            SELECT *
            FROM (
            SELECT
                CAST(gp.locus_id AS VARCHAR) AS locus_id,
                COALESCE(NULLIF(CAST(gp.gene_symbol AS VARCHAR), ''), NULLIF(CAST(gp.ensembl_gene_id AS VARCHAR), '')) AS gene,
                CAST(gp.max_pip AS DOUBLE) AS max_PIP,
                CAST(gp.gene_prioritization_score AS DOUBLE) AS gene_score,
                COALESCE(ote.eqtl_support_flag, FALSE) AS eqtl_support_flag,
                lsg.locus_type
            FROM gene_prioritization gp
            LEFT JOIN locus_signal lsg ON CAST(gp.locus_id AS VARCHAR) = lsg.locus_id
            LEFT JOIN ot_eqtl ote
                   ON COALESCE(NULLIF(CAST(gp.gene_symbol AS VARCHAR), ''), NULLIF(CAST(gp.ensembl_gene_id AS VARCHAR), '')) = ote.gene_key
            ) t
            WHERE gene IS NOT NULL;
            """
        )

        out_dir = ROOT / "results/database"
        out_dir.mkdir(parents=True, exist_ok=True)
        conn.execute(
            "COPY target_variant_table TO ? (DELIMITER '\\t', HEADER, FORMAT CSV);",
            [str(out_dir / "target_variant_table.tsv")],
        )
        conn.execute(
            "COPY target_gene_table TO ? (DELIMITER '\\t', HEADER, FORMAT CSV);",
            [str(out_dir / "target_gene_table.tsv")],
        )
    finally:
        conn.close()


def _get_counts_and_sample(db_path: Path) -> tuple[int, int, int, pl.DataFrame]:
    conn = duckdb.connect(str(db_path), read_only=True)
    try:
        loci_n = conn.execute("SELECT COUNT(DISTINCT locus_id) FROM target_variant_table;").fetchone()[0]
        variant_n = conn.execute("SELECT COUNT(*) FROM target_variant_table;").fetchone()[0]
        gene_n = conn.execute("SELECT COUNT(*) FROM target_gene_table;").fetchone()[0]
        sample = pl.DataFrame(
            conn.execute(
                "SELECT * FROM target_variant_table WHERE PIP > 0.8 ORDER BY PIP DESC, variant_priority_score DESC LIMIT 10;"
            ).fetch_df()
        )
    finally:
        conn.close()
    return int(loci_n or 0), int(variant_n or 0), int(gene_n or 0), sample


def _rows(records: Iterable[dict]) -> list[tuple]:
    return [tuple(rec.values()) for rec in records]


def _export_postgres(db_path: Path, dsn: str) -> list[str]:
    try:
        import psycopg
    except ImportError as exc:
        raise RuntimeError(
            "PostgreSQL export requires `psycopg`. Install it with: python3 -m pip install psycopg[binary]"
        ) from exc

    dd = duckdb.connect(str(db_path), read_only=True)
    try:
        loci_df = pl.DataFrame(
            dd.execute(
                """
                WITH successful_loci AS (
                    SELECT DISTINCT CAST(locus_id AS VARCHAR) AS locus_id
                    FROM locus_status
                    WHERE CAST(success AS BOOLEAN) IS TRUE
                )
                SELECT
                    tv.locus_id,
                    COALESCE(TRY_CAST(ll.CHR AS INTEGER), TRY_CAST(ls.chromosome AS INTEGER), TRY_CAST(tv.CHR AS INTEGER)) AS chromosome,
                    TRY_CAST(ll.locus_start AS BIGINT) AS locus_start,
                    TRY_CAST(ll.locus_end AS BIGINT) AS locus_end,
                    tv.locus_type,
                    COALESCE(CAST(ls.success AS BOOLEAN), FALSE) AS success
                FROM (SELECT DISTINCT locus_id, CHR, locus_type FROM target_variant_table) tv
                INNER JOIN successful_loci sl ON tv.locus_id = sl.locus_id
                LEFT JOIN lead_loci ll ON tv.locus_id = CAST(ll.locus_id AS VARCHAR)
                LEFT JOIN locus_status ls ON tv.locus_id = CAST(ls.locus_id AS VARCHAR);
                """
            ).fetch_df()
        )
        variants_df = pl.DataFrame(
            dd.execute(
                """
                WITH successful_loci AS (
                    SELECT DISTINCT CAST(locus_id AS VARCHAR) AS locus_id
                    FROM locus_status
                    WHERE CAST(success AS BOOLEAN) IS TRUE
                )
                SELECT DISTINCT
                    variant_id,
                    locus_id,
                    TRY_CAST(CHR AS INTEGER) AS chr,
                    TRY_CAST(BP AS BIGINT) AS bp,
                    CAST(PIP AS DOUBLE) AS pip,
                    CAST(credible_set_id AS VARCHAR) AS credible_set_id,
                    CAST(consequence AS VARCHAR) AS consequence,
                    CAST(variant_priority_score AS DOUBLE) AS variant_priority_score
                FROM target_variant_table
                WHERE variant_id IS NOT NULL
                  AND locus_id IN (SELECT locus_id FROM successful_loci);
                """
            ).fetch_df()
        )
        genes_df = pl.DataFrame(
            dd.execute(
                """
                WITH successful_loci AS (
                    SELECT DISTINCT CAST(locus_id AS VARCHAR) AS locus_id
                    FROM locus_status
                    WHERE CAST(success AS BOOLEAN) IS TRUE
                )
                SELECT DISTINCT
                    gene AS gene_symbol,
                    locus_id,
                    CAST(gene_score AS DOUBLE) AS gene_score,
                    CAST(eqtl_support_flag AS BOOLEAN) AS eqtl_support_flag,
                    TRUE AS candidate_flag
                FROM target_gene_table
                WHERE gene IS NOT NULL
                  AND locus_id IN (SELECT locus_id FROM successful_loci);
                """
            ).fetch_df()
        )
    finally:
        dd.close()

    created = []
    with psycopg.connect(dsn) as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS loci (
                    locus_id TEXT PRIMARY KEY,
                    chromosome INT,
                    locus_start BIGINT,
                    locus_end BIGINT,
                    locus_type TEXT,
                    success BOOLEAN
                );
                """
            )
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS variants (
                    variant_id TEXT PRIMARY KEY,
                    locus_id TEXT REFERENCES loci(locus_id),
                    chr INT,
                    bp BIGINT,
                    pip DOUBLE PRECISION,
                    credible_set_id TEXT,
                    consequence TEXT,
                    variant_priority_score DOUBLE PRECISION
                );
                """
            )
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS genes (
                    gene_symbol TEXT,
                    locus_id TEXT REFERENCES loci(locus_id),
                    gene_score DOUBLE PRECISION,
                    eqtl_support_flag BOOLEAN,
                    candidate_flag BOOLEAN
                );
                """
            )

            # Truncate related tables in a single statement to satisfy FK constraints.
            cur.execute("TRUNCATE TABLE variants, genes, loci;")

            if loci_df.height:
                cur.executemany(
                    "INSERT INTO loci (locus_id, chromosome, locus_start, locus_end, locus_type, success) VALUES (%s, %s, %s, %s, %s, %s)",
                    _rows(loci_df.to_dicts()),
                )
            if variants_df.height:
                cur.executemany(
                    "INSERT INTO variants (variant_id, locus_id, chr, bp, pip, credible_set_id, consequence, variant_priority_score) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)",
                    _rows(variants_df.to_dicts()),
                )
            if genes_df.height:
                cur.executemany(
                    "INSERT INTO genes (gene_symbol, locus_id, gene_score, eqtl_support_flag, candidate_flag) VALUES (%s, %s, %s, %s, %s)",
                    _rows(genes_df.to_dicts()),
                )
        conn.commit()
    created.extend(["loci", "variants", "genes"])
    return created


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Integrate Phase 2 outputs into DuckDB + optional PostgreSQL.")
    parser.add_argument(
        "--duckdb-path",
        default=str(ROOT / "results/database/target_discovery.duckdb"),
        help="DuckDB file path.",
    )
    parser.add_argument(
        "--postgres-dsn",
        default=os.getenv("TARGET_DISCOVERY_PG_DSN", ""),
        help="PostgreSQL DSN. If omitted, PostgreSQL export is skipped.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    db_path = Path(args.duckdb_path)

    _print_step(f"python3 {Path(__file__).name} --duckdb-path {db_path}")
    _create_duckdb_tables(db_path)

    loci_n, var_n, gene_n, sample = _get_counts_and_sample(db_path)
    print(f"DuckDB validation: loci={loci_n}, variants={var_n}, genes={gene_n}")
    print("Sample query: SELECT * FROM target_variant_table WHERE PIP > 0.8 LIMIT 10;")
    if sample.height:
        print(sample)
    else:
        print("No rows returned for PIP > 0.8.")

    if args.postgres_dsn:
        _print_step("postgres export (DSN provided)")
        tables = _export_postgres(db_path, args.postgres_dsn)
        print(f"PostgreSQL tables created/refreshed: {', '.join(tables)}")
    else:
        print("PostgreSQL export skipped: no DSN provided. Set TARGET_DISCOVERY_PG_DSN or pass --postgres-dsn.")


if __name__ == "__main__":
    main()
