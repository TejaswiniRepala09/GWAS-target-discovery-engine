"""Variant-level and gene-level prioritization tables."""

from __future__ import annotations

from typing import Any

import polars as pl


def _p_col(df: pl.DataFrame) -> str:
    for col in ["P_ANALYSIS", "P_GC", "P"]:
        if col in df.columns:
            return col
    raise ValueError("Missing p-value column; expected one of P_ANALYSIS, P_GC, P")


def build_variant_priority_table(
    variants: pl.DataFrame,
    consequence_scores: pl.DataFrame,
    pip_summary: pl.DataFrame,
) -> pl.DataFrame:
    """Merge GWAS + consequence + PIP with transparent score components."""
    pcol = _p_col(variants)

    out = variants.with_columns((-pl.col(pcol).clip(1e-300, 1.0).log10()).alias("component_gwas_log10p"))

    if consequence_scores.height:
        join_key = "SNP" if "SNP" in out.columns else "MARKERNAME"
        cons_key = "Uploaded_variation" if "Uploaded_variation" in consequence_scores.columns else join_key
        out = out.join(
            consequence_scores.select([c for c in [cons_key, "consequence_severity_score", "Consequence_term", "gene_symbol", "ensembl_gene_id"] if c in consequence_scores.columns]),
            left_on=join_key,
            right_on=cons_key,
            how="left",
        )

    if pip_summary.height and "SNP" in out.columns and "SNP" in pip_summary.columns:
        pip_cols = ["SNP", "PIP"]
        if "locus_id" not in out.columns and "locus_id" in pip_summary.columns:
            pip_cols.append("locus_id")
        out = out.join(pip_summary.select([c for c in pip_cols if c in pip_summary.columns]), on="SNP", how="left")

    if "consequence_severity_score" in out.columns:
        out = out.with_columns(pl.col("consequence_severity_score").cast(pl.Float64, strict=False).fill_null(0.0).alias("component_consequence"))
    else:
        out = out.with_columns(pl.lit(0.0).alias("component_consequence"))

    if "PIP" in out.columns:
        out = out.with_columns(pl.col("PIP").cast(pl.Float64, strict=False).fill_null(0.0).alias("component_pip"))
    else:
        out = out.with_columns(pl.lit(0.0).alias("component_pip"))

    out = out.with_columns(
        (pl.col("component_gwas_log10p") * 0.4 + pl.col("component_consequence") * 0.3 + pl.col("component_pip") * 0.3).alias("variant_priority_score")
    )

    return out.sort("variant_priority_score", descending=True)


def build_gene_priority_table(variant_scores: pl.DataFrame, variant_gene_map: pl.DataFrame) -> pl.DataFrame:
    """Aggregate variant evidence to gene-level priorities."""
    if variant_scores.height == 0:
        return pl.DataFrame()

    merged = variant_scores
    if variant_gene_map.height and "SNP" in variant_scores.columns and "SNP" in variant_gene_map.columns:
        merged = variant_scores.join(
            variant_gene_map.select([c for c in ["SNP", "gene_symbol", "ensembl_gene_id"] if c in variant_gene_map.columns]),
            on="SNP",
            how="left",
        )

    pcol = "P_ANALYSIS" if "P_ANALYSIS" in merged.columns else ("P_GC" if "P_GC" in merged.columns else "P")

    grouped = merged.group_by([c for c in ["gene_symbol", "ensembl_gene_id"] if c in merged.columns]).agg(
        pl.min(pcol).alias("strongest_pvalue"),
        pl.max("PIP").alias("max_pip") if "PIP" in merged.columns else pl.lit(None).alias("max_pip"),
        pl.max("consequence_severity_score").alias("best_consequence_severity") if "consequence_severity_score" in merged.columns else pl.lit(None).alias("best_consequence_severity"),
        pl.n_unique("locus_id").alias("supporting_locus_count") if "locus_id" in merged.columns else pl.lit(None).alias("supporting_locus_count"),
        pl.col("SNP").drop_nulls().first().alias("best_variant") if "SNP" in merged.columns else pl.lit(None).alias("best_variant"),
        pl.max("variant_priority_score").alias("gene_prioritization_score"),
    ).with_columns(
        (
            pl.lit("strongest_pvalue=")
            + pl.col("strongest_pvalue").cast(pl.Utf8)
            + pl.lit("; max_pip=")
            + pl.col("max_pip").cast(pl.Utf8)
            + pl.lit("; best_variant=")
            + pl.col("best_variant").cast(pl.Utf8)
        ).alias("prioritization_reason")
    )

    return grouped.sort("gene_prioritization_score", descending=True)
