"""LD resource and matrix contract validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl

from src.utils.path_utils import first_existing_dir


def detect_ld_resource_dir(fine_cfg: dict[str, Any]) -> Path | None:
    """Find first existing LD reference dir from config candidates."""
    candidates = fine_cfg["paths"].get("ld_matrix_dir_candidates", [])
    return first_existing_dir(candidates)


def validate_ld_resources_for_chromosomes(chromosomes: list[int], fine_cfg: dict[str, Any]) -> pl.DataFrame:
    """Validate chromosome-level 1000G VCF/TBI availability."""
    ld_dir = detect_ld_resource_dir(fine_cfg)
    rows: list[dict[str, Any]] = []
    for chrom in chromosomes:
        prefix = f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        vcf = ld_dir / prefix if ld_dir else Path(prefix)
        tbi = ld_dir / f"{prefix}.tbi" if ld_dir else Path(f"{prefix}.tbi")
        rows.append(
            {
                "CHR": int(chrom),
                "ld_dir": str(ld_dir) if ld_dir else None,
                "vcf_expected": str(vcf),
                "tbi_expected": str(tbi),
                "vcf_exists": bool(vcf.exists()) if ld_dir else False,
                "tbi_exists": bool(tbi.exists()) if ld_dir else False,
            }
        )
    return pl.from_dicts(rows) if rows else pl.DataFrame()


def validate_locus_ld_matrix_contract(bundle_manifest: pl.DataFrame, fine_cfg: dict[str, Any]) -> pl.DataFrame:
    """Validate whether per-locus LD matrices are present in bundle directories."""
    rows: list[dict[str, Any]] = []
    filename = fine_cfg["output_contract"].get("ld_matrix_filename", "ld_matrix.tsv")
    for row in bundle_manifest.iter_rows(named=True):
        bundle_dir = Path(fine_cfg["paths"]["bundles_dir"]) / str(row["locus_id"])
        ld_path = bundle_dir / filename
        rows.append(
            {
                "locus_id": row["locus_id"],
                "ld_matrix_path": str(ld_path),
                "ld_matrix_present": ld_path.exists() and ld_path.stat().st_size > 0,
            }
        )
    return pl.from_dicts(rows) if rows else pl.DataFrame()
