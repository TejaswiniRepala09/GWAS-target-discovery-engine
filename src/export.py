"""Export helpers for downstream fine-mapping and VEP annotation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl


def export_locus_windows(df: pl.DataFrame, lead_df: pl.DataFrame, settings: dict[str, Any]) -> pl.DataFrame:
    """Export one per-locus variant table around each lead variant for future SuSiE input."""
    out_dir = Path(settings["paths"]["locus_output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    window = int(settings["analysis"]["susie_window_bp"])
    exported_rows: list[dict[str, Any]] = []

    for row in lead_df.iter_rows(named=True):
        chrom = int(row["CHR"])
        bp = int(row["BP"])
        locus_id = row["locus_id"]
        start, end = max(1, bp - window), bp + window
        locus_df = df.filter((pl.col("CHR") == chrom) & (pl.col("BP") >= start) & (pl.col("BP") <= end))
        out_path = out_dir / f"{locus_id}_window.tsv"
        locus_df.write_csv(out_path, separator="\t")
        exported_rows.append(
            {"locus_id": locus_id, "chr": chrom, "lead_bp": bp, "window_start": start, "window_end": end, "path": str(out_path)}
        )
    return pl.from_dicts(exported_rows) if exported_rows else pl.DataFrame()


def export_vep_input(lead_df: pl.DataFrame, settings: dict[str, Any]) -> pl.DataFrame:
    """
    Build a VEP-preparation table.

    VEP exact input format can depend on local installation and genome build settings.
    This export keeps CHR/BP/alleles with optional SNP id to simplify conversion to the final VEP input form.
    """
    cols = [c for c in ["CHR", "BP", "SNP", "A1", "A2", "locus_id"] if c in lead_df.columns]
    vep_df = lead_df.select(cols).sort(["CHR", "BP"])
    output = Path(settings["paths"]["vep_output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    vep_df.write_csv(output, separator="\t")
    return vep_df
