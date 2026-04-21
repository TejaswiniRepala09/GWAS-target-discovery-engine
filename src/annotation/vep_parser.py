"""VEP input formatting, command generation, and output parsing."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any

import polars as pl

from src.utils.path_utils import ensure_dir
from src.utils.validation import require_columns, require_file

LOGGER = logging.getLogger(__name__)


def build_vep_input_file(vep_input_df: pl.DataFrame, vep_cfg: dict[str, Any]) -> Path:
    """
    Convert project VEP input table into Ensembl tab format.

    Output columns:
    CHR START END ALLELE STRAND ID
    """
    require_columns(vep_input_df, vep_cfg.get("required_input_columns", ["CHR", "BP", "A1", "A2"]), "vep_input.tsv")

    out_path = Path(vep_cfg["paths"]["vep_formatted_input"])
    ensure_dir(out_path.parent)

    formatted = vep_input_df.with_columns(
        pl.col("CHR").cast(pl.Utf8).alias("CHR"),
        pl.col("BP").cast(pl.Int64).alias("START"),
        pl.col("BP").cast(pl.Int64).alias("END"),
        (pl.col("A1").cast(pl.Utf8) + pl.lit("/") + pl.col("A2").cast(pl.Utf8)).alias("ALLELE"),
        pl.lit("+").alias("STRAND"),
        pl.when(pl.col("SNP").is_not_null())
        .then(pl.col("SNP").cast(pl.Utf8))
        .otherwise(pl.concat_str([pl.col("CHR").cast(pl.Utf8), pl.lit(":"), pl.col("BP").cast(pl.Utf8)]))
        .alias("ID"),
    ).select(["CHR", "START", "END", "ALLELE", "STRAND", "ID"])

    formatted.write_csv(out_path, separator="\t", include_header=False)
    LOGGER.info("Wrote VEP-formatted input: %s (%s rows)", out_path, formatted.height)
    return out_path


def _check_vep_runtime(vep_cfg: dict[str, Any]) -> dict[str, str | bool]:
    exe = vep_cfg["local_vep"].get("executable", "vep")
    cache_dir = Path(vep_cfg["local_vep"].get("cache_dir", "data/reference/vep"))
    return {
        "vep_in_path": shutil.which(exe) is not None,
        "cache_present": cache_dir.exists(),
        "cache_dir": str(cache_dir),
        "executable": exe,
    }


def generate_local_vep_command(vep_cfg: dict[str, Any]) -> tuple[str, dict[str, str | bool]]:
    """Generate exact local VEP command string and validation status."""
    status = _check_vep_runtime(vep_cfg)
    paths = vep_cfg["paths"]
    local = vep_cfg["local_vep"]
    assembly = vep_cfg.get("assembly", "GRCh37")

    cmd = " \\\n  ".join(
        [
            f"{status['executable']}",
            f"--input_file {paths['vep_formatted_input']}",
            "--format ensembl",
            f"--output_file {paths['vep_raw_output']}",
            "--tab",
            "--everything",
            "--cache",
            "--offline" if bool(local.get("offline", True)) else "",
            f"--dir_cache {local.get('cache_dir', 'data/reference/vep')}",
            f"--assembly {assembly}",
            "--force_overwrite",
            f"--fork {int(local.get('fork', 4))}",
            f"--stats_file {local.get('stats_file', 'results/annotations/vep_raw/vep_summary.html')}",
            f"2>&1 | tee {local.get('log_file', 'results/logs/phase2_vep.log')}",
        ]
    ).replace("  \\\n  ", " \\\n  ")
    return cmd, status


def parse_vep_output(vep_cfg: dict[str, Any], raw_output_path: str | Path | None = None) -> pl.DataFrame:
    """Parse tabular VEP output into normalized long-form annotation table."""
    source_path = raw_output_path if raw_output_path is not None else vep_cfg["paths"]["vep_raw_output"]
    raw_path = require_file(source_path, "VEP raw output")

    # VEP tab output includes many "##" metadata lines and a "#Uploaded_variation" header line.
    # We must keep that header (not treat it as a comment), otherwise columns are misparsed.
    header_idx: int | None = None
    with raw_path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            if line.startswith("#Uploaded_variation") or line.startswith("Uploaded_variation\t"):
                header_idx = idx
                break
    if header_idx is None:
        raise ValueError(
            f"Could not find VEP header line (#Uploaded_variation) in {raw_path}. "
            "Ensure the file is tabular VEP output produced with '--tab'."
        )

    df = pl.read_csv(raw_path, separator="\t", skip_rows=header_idx, null_values=["", "-", "."], infer_schema_length=10000)
    rename_map = {c: c.lstrip("#") for c in df.columns if c.startswith("#")}
    if rename_map:
        df = df.rename(rename_map)
    if "Uploaded_variation" not in df.columns:
        raise ValueError(
            f"Parsed VEP table from {raw_path} is missing Uploaded_variation column after header normalization."
        )
    df = df.with_columns(pl.col("Uploaded_variation").cast(pl.Utf8, strict=False))
    if df.select(pl.col("Uploaded_variation").drop_nulls().len()).item() == 0:
        raise ValueError(
            f"Uploaded_variation is entirely null/empty in {raw_path}. "
            "Check VEP input IDs and assembly/cache compatibility."
        )

    expected = vep_cfg.get("expected_output_columns", [])
    for col in expected:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))

    if "Consequence" in df.columns:
        df = df.with_columns(pl.col("Consequence").cast(pl.Utf8).str.split(",").alias("_cons"))
        df = df.explode("_cons").with_columns(
            pl.col("_cons").cast(pl.Utf8).str.replace_all("&", ",").str.split(",").alias("_cons2")
        ).explode("_cons2").with_columns(
            pl.col("_cons2").cast(pl.Utf8).str.strip_chars().alias("Consequence_term")
        ).drop(["_cons", "_cons2"])
    else:
        df = df.with_columns(pl.lit(None).alias("Consequence_term"))

    if "Location" in df.columns:
        df = df.with_columns(
            pl.col("Location").cast(pl.Utf8).str.split(":").list.get(0).alias("CHR"),
            pl.col("Location").cast(pl.Utf8).str.split(":").list.get(1).alias("_pos"),
        ).with_columns(
            pl.col("_pos").cast(pl.Utf8).str.split("-").list.get(0).cast(pl.Int64, strict=False).alias("BP")
        ).drop("_pos")

    out_cols = [
        "Uploaded_variation",
        "CHR",
        "BP",
        "Location",
        "Allele",
        "Gene",
        "SYMBOL",
        "Feature",
        "Consequence",
        "Consequence_term",
        "IMPACT",
        "HGVSc",
        "HGVSp",
        "Amino_acids",
        "Protein_position",
        "Existing_variation",
    ]
    out = df.select([c for c in out_cols if c in df.columns])
    out_path = Path(vep_cfg["paths"]["variant_annotations_output"])
    ensure_dir(out_path.parent)
    out.write_csv(out_path, separator="\t")
    LOGGER.info("Parsed VEP annotations: %s (%s rows)", out_path, out.height)
    return out
