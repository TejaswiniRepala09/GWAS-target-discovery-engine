#!/usr/bin/env python3
"""Prepare coloc-ready artifacts for representative loci.

This script is intentionally lightweight and honest:
- It DOES NOT fabricate coloc posterior probabilities.
- It prepares GWAS-side coloc inputs for representative loci.
- It records exactly what eQTL-side data are missing for full coloc runs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import polars as pl


ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results" / "coloc"
REP_LOCI = ["chr7_locus19", "chr7_locus14", "chr7_locus10"]


@dataclass
class ColocStatus:
    locus_id: str
    candidate_genes: str
    coloc_run: bool
    data_used: str
    missing_data: str
    shared_signal_support: str
    confidence_interpretation: str


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def pick_summary_stats_path(locus_id: str) -> Path:
    base = ROOT / "results" / "fine_mapping" / "loci_inputs" / locus_id
    harmonized = base / "summary_stats_harmonized.tsv"
    standard = base / "summary_stats.tsv"
    raw = base / "summary_stats_raw.tsv"
    for p in [harmonized, standard, raw]:
        if p.exists():
            return p
    raise FileNotFoundError(f"No summary stats file found for {locus_id} under {base}")


def extract_candidate_genes(gene_prioritization: pl.DataFrame, locus_id: str, top_n: int = 5) -> list[str]:
    if "locus_id" not in gene_prioritization.columns:
        return []

    df = gene_prioritization.filter(pl.col("locus_id") == locus_id)
    if df.is_empty():
        return []

    score_col = "gene_prioritization_score_final" if "gene_prioritization_score_final" in df.columns else (
        "gene_prioritization_score" if "gene_prioritization_score" in df.columns else None
    )
    if score_col:
        df = df.sort(score_col, descending=True)

    genes = (
        df.select(pl.col("gene_symbol"))
        .drop_nulls()
        .head(top_n)
        .to_series()
        .to_list()
    )
    return [str(g) for g in genes if str(g).strip()]


def build_gwas_coloc_table(summary_df: pl.DataFrame) -> pl.DataFrame:
    required = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]
    missing = [c for c in required if c not in summary_df.columns]
    if missing:
        raise ValueError(f"Missing required GWAS columns for coloc input: {missing}")

    out = (
        summary_df.select(required + [c for c in ["N", "EAF"] if c in summary_df.columns])
        .with_columns(
            [
                pl.col("SNP").cast(pl.Utf8),
                pl.col("CHR").cast(pl.Utf8).str.replace("chr", ""),
                pl.col("BP").cast(pl.Int64),
                pl.col("A1").cast(pl.Utf8),
                pl.col("A2").cast(pl.Utf8),
                pl.col("BETA").cast(pl.Float64),
                pl.col("SE").cast(pl.Float64),
                pl.col("P").cast(pl.Float64),
            ]
        )
        .with_columns(
            [
                (pl.col("SE") ** 2).alias("varbeta"),
                (pl.col("BETA") / pl.col("SE")).alias("z"),
            ]
        )
    )
    return out


def write_eqtl_template(path: Path, locus_id: str, genes: Iterable[str]) -> None:
    template = pl.DataFrame(
        {
            "locus_id": [locus_id],
            "gene_symbol": [",".join(genes) if genes else ""],
            "required_columns": [
                "SNP,CHR,BP,effect_allele,other_allele,beta,se,p,n,maf,tissue,gene_symbol,build"
            ],
            "notes": [
                "Provide eQTL summary statistics for the same build (GRCh37), harmonized alleles, and overlapping SNP set."
            ],
        }
    )
    template.write_csv(path, separator="\t")


def classify_support(locus_id: str, genes: list[str]) -> tuple[str, str]:
    # Keep this conservative and report-style, based on existing project interpretation language.
    if locus_id == "chr7_locus19":
        return (
            "weak_to_moderate_indirect_support",
            "Coloc not run. Existing project notes suggest broader cross-tissue regulatory support for PRKAG2, but no definitive kidney-specific shared-signal test yet.",
        )
    if locus_id == "chr7_locus14":
        return (
            "inconclusive",
            "Coloc not run. Current evidence is predominantly non-coding and gene assignment remains uncertain without matched eQTL summary statistics.",
        )
    return (
        "inconclusive",
        "Coloc not run. Locus remains complex/diffuse; no matched eQTL summary statistics available to test shared causal signal.",
    )


def main() -> None:
    ensure_dir(RESULTS_DIR)

    gene_path = ROOT / "results" / "target_prioritization" / "gene_prioritization.tsv"
    gene_prioritization = pl.read_csv(gene_path, separator="\t") if gene_path.exists() else pl.DataFrame()

    statuses: list[ColocStatus] = []

    for locus_id in REP_LOCI:
        locus_dir = RESULTS_DIR / locus_id
        ensure_dir(locus_dir)

        summary_path = pick_summary_stats_path(locus_id)
        summary_df = pl.read_csv(summary_path, separator="\t")
        gwas_coloc = build_gwas_coloc_table(summary_df)
        gwas_coloc_path = locus_dir / "gwas_coloc_input.tsv"
        gwas_coloc.write_csv(gwas_coloc_path, separator="\t")

        genes = extract_candidate_genes(gene_prioritization, locus_id)
        write_eqtl_template(locus_dir / "eqtl_input_template.tsv", locus_id, genes)

        support, interpretation = classify_support(locus_id, genes)
        statuses.append(
            ColocStatus(
                locus_id=locus_id,
                candidate_genes=",".join(genes) if genes else "",
                coloc_run=False,
                data_used=f"GWAS summary stats from {summary_path.relative_to(ROOT)}",
                missing_data=(
                    "Matched eQTL summary statistics (same build, tissue-resolved, harmonized alleles, overlapping SNPs)"
                ),
                shared_signal_support=support,
                confidence_interpretation=interpretation,
            )
        )

    status_df = pl.DataFrame([s.__dict__ for s in statuses])
    status_df.write_csv(RESULTS_DIR / "coloc_status.tsv", separator="\t")

    md_lines = [
        "# Colocalization Summary (Representative Loci)",
        "",
        "This layer prepares **coloc-ready inputs** and records missing requirements for full coloc.",
        "",
        "## What was run",
        "- GWAS-side coloc input preparation for representative loci.",
        "- Candidate-gene extraction from existing prioritization outputs.",
        "- Missing-data manifest generation for production coloc.",
        "",
        "## Why full coloc was not run",
        "- The repository currently does not contain matched eQTL summary statistics per representative locus and tissue.",
        "- Without matched eQTL effect/SE/p-value data on overlapping SNPs, coloc posterior probabilities (PP.H0-H4) cannot be estimated reliably.",
        "",
        "## Locus-level status",
        "",
        "| locus_id | candidate_genes | coloc_run | shared_signal_support | confidence_interpretation |",
        "|---|---|---|---|---|",
    ]
    for row in status_df.to_dicts():
        md_lines.append(
            f"| {row['locus_id']} | {row['candidate_genes'] or 'NA'} | {row['coloc_run']} | {row['shared_signal_support']} | {row['confidence_interpretation']} |"
        )

    md_lines.extend(
        [
            "",
            "## Production-grade coloc requirements",
            "For each locus and candidate gene/tissue, provide eQTL summary statistics with:",
            "- same genome build (GRCh37)",
            "- SNP-level fields: SNP, CHR, BP, effect allele, other allele, beta, se, p, n, maf",
            "- tissue and gene metadata",
            "- robust harmonization to GWAS alleles and SNP overlap.",
            "",
            "Generated files:",
            "- `results/coloc/<locus_id>/gwas_coloc_input.tsv`",
            "- `results/coloc/<locus_id>/eqtl_input_template.tsv`",
            "- `results/coloc/coloc_status.tsv`",
        ]
    )

    (ROOT / "results" / "reports").mkdir(parents=True, exist_ok=True)
    (ROOT / "results" / "reports" / "coloc_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print("Coloc-ready layer completed.")
    print(f"Wrote: {RESULTS_DIR / 'coloc_status.tsv'}")
    print(f"Wrote: {ROOT / 'results' / 'reports' / 'coloc_summary.md'}")


if __name__ == "__main__":
    main()

