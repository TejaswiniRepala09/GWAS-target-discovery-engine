#!/usr/bin/env python3
"""Select top-ranked loci chromosomes for targeted 1000 Genomes downloads."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import polars as pl


def load_tsv(path: str | Path) -> pl.DataFrame:
    """Load a tab-delimited file with basic existence checks."""
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    return pl.read_csv(file_path, separator="\t")


def detect_pvalue_column(columns: Iterable[str]) -> str:
    """
    Detect the p-value column with preference order:
    1) P_ANALYSIS
    2) P
    3) common alternatives (including P_GC)
    """
    columns_list = list(columns)
    preferred = ["P_ANALYSIS", "P", "P_GC", "P.value.GC", "P.value", "P-value"]
    for col in preferred:
        if col in columns_list:
            return col

    # Fallback to "p-value like" columns if standard names are absent.
    fallback = [
        col
        for col in columns_list
        if ("p" in col.lower()) and ("se" not in col.lower()) and ("position" not in col.lower())
    ]
    if fallback:
        return fallback[0]

    raise ValueError(
        "Could not find a p-value column. Expected one of: "
        "P_ANALYSIS, P, P_GC, P.value.GC, P.value, P-value."
    )


def get_top_ranked_variants(lead_variants: pl.DataFrame, pvalue_col: str, top_n: int) -> pl.DataFrame:
    """Cast p-values to numeric, sort ascending, and keep top N strongest signals."""
    required = ["CHR", "BP"]
    missing = [c for c in required if c not in lead_variants.columns]
    if missing:
        raise ValueError(f"Missing required columns in lead_variants.tsv: {missing}")

    ranked = (
        lead_variants.with_columns(pl.col(pvalue_col).cast(pl.Float64, strict=False).alias("_pvalue_numeric"))
        .filter(pl.col("_pvalue_numeric").is_not_null())
        .sort("_pvalue_numeric")
        .head(top_n)
    )
    return ranked


def maybe_join_locus_bounds(top_variants: pl.DataFrame, lead_loci: pl.DataFrame | None) -> pl.DataFrame:
    """Optionally join locus_start/locus_end from lead_loci by locus_id when available."""
    if lead_loci is None:
        return top_variants
    if "locus_id" not in top_variants.columns or "locus_id" not in lead_loci.columns:
        return top_variants

    loci_cols = [c for c in ["locus_id", "locus_start", "locus_end"] if c in lead_loci.columns]
    return top_variants.join(lead_loci.select(loci_cols), on="locus_id", how="left")


def build_1000g_download_table(chromosomes: pl.Series) -> pl.DataFrame:
    """Build download-ready 1000 Genomes VCF and index filenames for each chromosome."""
    chrom_list = [str(x) for x in chromosomes.to_list()]
    vcf = [
        f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        for chrom in chrom_list
    ]
    tbi = [f"{name}.tbi" for name in vcf]
    return pl.DataFrame({"CHR": chrom_list, "vcf_file": vcf, "tbi_file": tbi})


def main() -> None:
    """CLI entrypoint."""
    parser = argparse.ArgumentParser(
        description="Select top chromosomes from lead variants for 1000G download planning."
    )
    parser.add_argument(
        "--lead-variants",
        default="results/tables/lead_variants.tsv",
        help="Path to lead_variants.tsv",
    )
    parser.add_argument(
        "--lead-loci",
        default="results/tables/lead_loci.tsv",
        help="Path to lead_loci.tsv (optional join source)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top variants/loci to retain (default: 5)",
    )
    parser.add_argument(
        "--skip-loci-join",
        action="store_true",
        help="Skip optional locus_id join with lead_loci.tsv",
    )
    parser.add_argument(
        "--output-top-variants",
        default="results/tables/top_ranked_variants.tsv",
        help="Output path for ranked top variants",
    )
    parser.add_argument(
        "--output-chromosomes",
        default="results/tables/top_loci_chromosomes.tsv",
        help="Output path for unique chromosomes + 1000G files",
    )
    args = parser.parse_args()

    if args.top_n <= 0:
        raise ValueError("--top-n must be >= 1")

    # 1) Load lead variants and detect p-value column with preference order.
    lead_variants = load_tsv(args.lead_variants)
    pvalue_col = detect_pvalue_column(lead_variants.columns)

    # 2) Rank by p-value and keep top N rows.
    top_variants = get_top_ranked_variants(lead_variants, pvalue_col, args.top_n)

    # 3) Optional locus_id join to add locus_start/locus_end context.
    lead_loci_df = None
    if not args.skip_loci_join:
        try:
            lead_loci_df = load_tsv(args.lead_loci)
        except FileNotFoundError:
            lead_loci_df = None
    top_variants = maybe_join_locus_bounds(top_variants, lead_loci_df)

    # 4) Build unique chromosome list and file names for 1000G downloads.
    unique_chr = top_variants.select("CHR").unique().sort("CHR")["CHR"]
    download_table = build_1000g_download_table(unique_chr)

    # 5) Write outputs.
    out_top = Path(args.output_top_variants)
    out_chr = Path(args.output_chromosomes)
    out_top.parent.mkdir(parents=True, exist_ok=True)
    out_chr.parent.mkdir(parents=True, exist_ok=True)

    # Keep a clean, reviewer-friendly top table.
    preferred_cols = ["locus_id", "CHR", "BP", pvalue_col, "locus_start", "locus_end", "SNP"]
    cols_to_write = [c for c in preferred_cols if c in top_variants.columns]
    top_variants.select(cols_to_write).write_csv(out_top, separator="\t")
    download_table.write_csv(out_chr, separator="\t")

    # 6) Print requested summaries.
    print(f"Using p-value column: {pvalue_col}")
    print("\nTop loci / variants:")
    print(top_variants.select(cols_to_write))
    print("\nUnique chromosomes:")
    print(unique_chr.to_list())
    print("\n1000 Genomes download files:")
    for row in download_table.iter_rows(named=True):
        print(row["vcf_file"])
        print(row["tbi_file"])

    print(f"\nWrote: {out_top}")
    print(f"Wrote: {out_chr}")


if __name__ == "__main__":
    main()

