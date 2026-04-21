"""Distance-based lead locus extraction for GWAS summary statistics."""

from __future__ import annotations

from typing import Any

import polars as pl


def extract_distance_based_lead_loci(df: pl.DataFrame, settings: dict[str, Any]) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Extract lead variants using physical distance, not LD clumping.

    This is a practical approximation for prioritization. Nearby significant SNPs can still represent
    the same signal due to linkage disequilibrium, which true fine-mapping addresses later.
    """
    cfg = settings["analysis"]
    p_col = cfg["selected_pvalue_column"]
    sig = df.filter(pl.col(p_col) < cfg["genome_wide_significance_threshold"]).sort(["CHR", p_col, "BP"])
    if sig.height == 0:
        return pl.DataFrame(), pl.DataFrame()

    lead_rows: list[dict[str, Any]] = []
    loci_rows: list[dict[str, Any]] = []
    distance_bp = int(cfg["lead_locus_distance_bp"])

    for chrom in sig["CHR"].unique().sort().to_list():
        chr_df = sig.filter(pl.col("CHR") == chrom).sort("BP")
        current_end = -1
        locus_idx = 0
        locus_variants: list[dict[str, Any]] = []
        for row in chr_df.iter_rows(named=True):
            bp = int(row["BP"])
            if current_end >= 0 and bp > current_end:
                locus_id = f"chr{chrom}_locus{locus_idx}"
                best = min(locus_variants, key=lambda x: x[p_col])
                lead_rows.append({**best, "locus_id": locus_id})
                positions = [int(x["BP"]) for x in locus_variants]
                loci_rows.append(
                    {"locus_id": locus_id, "CHR": chrom, "locus_start": min(positions), "locus_end": max(positions)}
                )
                locus_idx += 1
                locus_variants = []
            locus_variants.append(row)
            current_end = bp + distance_bp
        if locus_variants:
            locus_id = f"chr{chrom}_locus{locus_idx}"
            best = min(locus_variants, key=lambda x: x[p_col])
            lead_rows.append({**best, "locus_id": locus_id})
            positions = [int(x["BP"]) for x in locus_variants]
            loci_rows.append({"locus_id": locus_id, "CHR": chrom, "locus_start": min(positions), "locus_end": max(positions)})

    lead_df = pl.from_dicts(lead_rows).sort(["CHR", "BP"])
    loci_df = pl.from_dicts(loci_rows).sort(["CHR", "locus_start"])
    return lead_df, loci_df
