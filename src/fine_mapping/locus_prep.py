"""Build per-locus fine-mapping input bundles."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from src.utils.path_utils import ensure_dir
from src.utils.validation import require_columns

LOGGER = logging.getLogger(__name__)


def _bundle_paths(bundle_dir: Path, cfg: dict[str, Any]) -> dict[str, Path]:
    contract = cfg["output_contract"]
    return {
        "summary": bundle_dir / contract["summary_stats_filename"],
        "zscores": bundle_dir / contract["zscores_filename"],
        "metadata": bundle_dir / contract["metadata_filename"],
        "ld_placeholder": bundle_dir / contract["ld_matrix_filename"],
    }


def create_locus_bundles(
    cleaned_df: pl.DataFrame,
    lead_df: pl.DataFrame,
    fine_cfg: dict[str, Any],
) -> pl.DataFrame:
    """Create one input bundle per locus for downstream susieR execution."""
    require_columns(cleaned_df, fine_cfg.get("required_columns", ["CHR", "BP", "SNP", "BETA", "SE"]), "cleaned GWAS")
    require_columns(lead_df, ["locus_id", "CHR", "BP"], "lead variants")

    bundles_root = ensure_dir(fine_cfg["paths"]["bundles_dir"])
    window_bp = int(fine_cfg.get("window_bp", 500000))

    rows: list[dict[str, Any]] = []
    for lead in lead_df.iter_rows(named=True):
        locus_id = str(lead["locus_id"])
        chrom = int(lead["CHR"])
        bp = int(lead["BP"])
        start = max(1, bp - window_bp)
        end = bp + window_bp

        locus = cleaned_df.filter(
            (pl.col("CHR") == chrom) & (pl.col("BP") >= start) & (pl.col("BP") <= end)
        ).sort(["CHR", "BP"])
        if locus.height == 0:
            continue

        bundle_dir = ensure_dir(bundles_root / locus_id)
        paths = _bundle_paths(bundle_dir, fine_cfg)

        summary_cols = [c for c in ["SNP", "CHR", "BP", "A1", "A2", "EAF", "BETA", "SE", "P", "P_GC", "P_ANALYSIS", "N"] if c in locus.columns]
        summary = locus.select(summary_cols)
        summary.write_csv(paths["summary"], separator="\t")

        if "BETA" in summary.columns and "SE" in summary.columns:
            zscores = summary.with_columns((pl.col("BETA") / pl.col("SE")).alias("Z")).select([c for c in ["SNP", "Z"] if c in summary.columns or c == "Z"])
            zscores.write_csv(paths["zscores"], separator="\t")

        if not paths["ld_placeholder"].exists():
            paths["ld_placeholder"].write_text(
                "# Provide locus LD matrix here. Expected first column SNP and square numeric matrix.\n",
                encoding="utf-8",
            )

        meta = {
            "locus_id": locus_id,
            "chr": chrom,
            "lead_bp": bp,
            "window_start": start,
            "window_end": end,
            "assembly": fine_cfg.get("assembly", "GRCh37"),
            "variant_count": int(summary.height),
            "summary_stats": str(paths["summary"]),
            "zscores": str(paths["zscores"]),
            "ld_matrix_expected": str(paths["ld_placeholder"]),
        }
        with paths["metadata"].open("w", encoding="utf-8") as handle:
            yaml.safe_dump(meta, handle, sort_keys=False)

        rows.append(meta)

    manifest = pl.from_dicts(rows) if rows else pl.DataFrame()
    LOGGER.info("Created %s fine-mapping bundles", manifest.height)
    return manifest
