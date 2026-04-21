"""Variant-to-gene mapping using VEP primary mapping and GENCODE fallback."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Any

import polars as pl

LOGGER = logging.getLogger(__name__)


def _read_gtf(gtf_path: Path) -> pl.DataFrame:
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    rows: list[dict[str, Any]] = []
    with opener(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = parts
            if feature != "gene":
                continue
            gene_id = None
            gene_name = None
            for item in attributes.split(";"):
                item = item.strip()
                if item.startswith("gene_id"):
                    gene_id = item.split(" ", 1)[1].replace('"', "").strip()
                elif item.startswith("gene_name"):
                    gene_name = item.split(" ", 1)[1].replace('"', "").strip()
            rows.append(
                {
                    "CHR": chrom.replace("chr", ""),
                    "gene_start": int(start),
                    "gene_end": int(end),
                    "ensembl_gene_id": gene_id,
                    "gene_symbol": gene_name,
                }
            )
    return pl.from_dicts(rows) if rows else pl.DataFrame()


def _nearest_gene_lookup(variants: pl.DataFrame, genes: pl.DataFrame) -> pl.DataFrame:
    if variants.height == 0 or genes.height == 0:
        return pl.DataFrame()

    genes_by_chr = {str(ch): genes.filter(pl.col("CHR") == str(ch)) for ch in genes["CHR"].unique().to_list()}
    out_rows: list[dict[str, Any]] = []

    for row in variants.select([c for c in ["SNP", "CHR", "BP", "locus_id"] if c in variants.columns]).iter_rows(named=True):
        chrom = str(row["CHR"])
        if chrom not in genes_by_chr or genes_by_chr[chrom].height == 0:
            continue
        bp = int(row["BP"])
        chr_genes = genes_by_chr[chrom].with_columns(
            pl.when((pl.col("gene_start") <= bp) & (pl.col("gene_end") >= bp))
            .then(pl.lit(0))
            .otherwise(pl.min_horizontal((pl.col("gene_start") - bp).abs(), (pl.col("gene_end") - bp).abs()))
            .alias("distance_bp")
        ).sort("distance_bp")
        best = chr_genes.row(0, named=True)
        out_rows.append(
            {
                "SNP": row.get("SNP"),
                "locus_id": row.get("locus_id"),
                "CHR": chrom,
                "BP": bp,
                "gene_symbol": best.get("gene_symbol"),
                "ensembl_gene_id": best.get("ensembl_gene_id"),
                "mapping_source": "nearest_gene_gencode",
                "distance_bp": int(best.get("distance_bp", 0)),
                "mapping_warning": "Nearest-gene fallback does not imply causal assignment.",
            }
        )

    return pl.from_dicts(out_rows) if out_rows else pl.DataFrame()


def build_variant_gene_map(
    variant_table: pl.DataFrame,
    annotation_table: pl.DataFrame,
    genes_dir: str | Path,
    output_path: str | Path,
) -> pl.DataFrame:
    """Build variant-gene mapping table with VEP-first + GENCODE fallback behavior."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    vep_map = pl.DataFrame()
    if annotation_table.height:
        rename_map = {}
        if "SYMBOL" in annotation_table.columns:
            rename_map["SYMBOL"] = "gene_symbol"
        if "Gene" in annotation_table.columns:
            rename_map["Gene"] = "ensembl_gene_id"
        ann = annotation_table.rename(rename_map)

        right_key = "Uploaded_variation" if "Uploaded_variation" in ann.columns else "SNP"
        left_key = "SNP" if "SNP" in variant_table.columns else right_key
        cols = [c for c in [right_key, "gene_symbol", "ensembl_gene_id", "Consequence_term"] if c in ann.columns]

        if cols:
            vep_map = variant_table.join(
                ann.select(cols),
                left_on=left_key,
                right_on=right_key,
                how="left",
            ).with_columns(
                pl.when(pl.col("gene_symbol").is_not_null() | pl.col("ensembl_gene_id").is_not_null())
                .then(pl.lit("vep"))
                .otherwise(pl.lit(None))
                .alias("mapping_source"),
                pl.lit(None).alias("distance_bp"),
                pl.lit(None).alias("mapping_warning"),
            )

    if vep_map.height == 0:
        vep_map = variant_table.with_columns(
            pl.lit(None).alias("gene_symbol"),
            pl.lit(None).alias("ensembl_gene_id"),
            pl.lit(None).alias("mapping_source"),
            pl.lit(None).alias("distance_bp"),
            pl.lit(None).alias("mapping_warning"),
        )

    unresolved = vep_map.filter(pl.col("gene_symbol").is_null() & pl.col("ensembl_gene_id").is_null())

    gtf_candidates = sorted(Path(genes_dir).glob("*.gtf*")) if Path(genes_dir).exists() else []
    nearest = pl.DataFrame()
    if unresolved.height and gtf_candidates:
        genes = _read_gtf(gtf_candidates[0])
        nearest = _nearest_gene_lookup(unresolved, genes)
    elif unresolved.height:
        LOGGER.warning("GENCODE GTF missing under %s; nearest-gene fallback skipped", genes_dir)

    out = vep_map.filter(pl.col("mapping_source").is_not_null()) if vep_map.height else pl.DataFrame()
    if nearest.height:
        out = pl.concat([out, nearest], how="diagonal_relaxed") if out.height else nearest

    out.write_csv(output_path, separator="\t")
    return out
