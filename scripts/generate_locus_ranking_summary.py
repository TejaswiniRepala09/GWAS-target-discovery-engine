#!/usr/bin/env python3
"""Generate locus-level ranking summary from completed fine-mapping outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _require_columns(df: pd.DataFrame, required: set[str], name: str) -> None:
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{name} missing required columns: {sorted(missing)}")


def _to_markdown_table(df: pd.DataFrame) -> str:
    """Render a markdown table without optional dependencies like tabulate."""
    if df.empty:
        return "_No rows._"
    cols = [str(c) for c in df.columns]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    body_lines = []
    for _, row in df.iterrows():
        vals = []
        for c in df.columns:
            v = row[c]
            if pd.isna(v):
                vals.append("")
            else:
                vals.append(str(v))
        body_lines.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, sep] + body_lines)


def _classify_locus(
    max_pip: float,
    max_cs_size: float,
    n_cs: int,
    n_nonzero: int,
) -> str:
    max_cs_size = float(max_cs_size) if pd.notna(max_cs_size) else np.nan
    diffuse = (pd.notna(max_cs_size) and max_cs_size > 10) or (n_nonzero > 50)

    if (max_pip >= 0.9) and (pd.notna(max_cs_size) and max_cs_size <= 3):
        return "strong"
    if (max_pip < 0.3) or diffuse:
        return "weak"
    if (0.3 <= max_pip < 0.9) or (n_cs > 1):
        return "moderate"
    return "moderate"


def main() -> None:
    root = Path.cwd()
    pip_path = root / "results/fine_mapping/pip_summary.tsv"
    cs_path = root / "results/fine_mapping/credible_sets.tsv"
    vp_path = root / "results/target_prioritization/variant_priority_scores.tsv"
    gp_path = root / "results/target_prioritization/gene_prioritization.tsv"
    status_path = root / "results/reports/multi_locus_status.tsv"

    for p in [pip_path, cs_path, vp_path, gp_path, status_path]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input file: {p}")

    pip = pd.read_csv(pip_path, sep="\t")
    cs = pd.read_csv(cs_path, sep="\t")
    vp = pd.read_csv(vp_path, sep="\t")
    gp = pd.read_csv(gp_path, sep="\t")
    status = pd.read_csv(status_path, sep="\t")

    _require_columns(pip, {"locus_id", "PIP", "SNP"}, "pip_summary.tsv")
    _require_columns(cs, {"locus_id", "credible_set_id", "SNP"}, "credible_sets.tsv")
    _require_columns(gp, {"locus_id", "gene_symbol"}, "gene_prioritization.tsv")
    _require_columns(status, {"locus_id", "success", "susie_converged"}, "multi_locus_status.tsv")

    pip["PIP"] = pd.to_numeric(pip["PIP"], errors="coerce").fillna(0.0)

    pip_metrics = (
        pip.groupby("locus_id", as_index=False)
        .agg(
            max_pip=("PIP", "max"),
            mean_pip=("PIP", "mean"),
            number_of_nonzero_pip_variants=("PIP", lambda s: int((s > 0).sum())),
            number_of_high_pip_variants_gt_01=("PIP", lambda s: int((s > 0.1).sum())),
            number_of_high_pip_variants_gt_05=("PIP", lambda s: int((s > 0.5).sum())),
        )
    )

    if len(cs):
        cs_sizes = cs.groupby(["locus_id", "credible_set_id"], as_index=False).size().rename(columns={"size": "credible_set_size"})
        cs_metrics = (
            cs_sizes.groupby("locus_id", as_index=False)
            .agg(
                number_of_credible_sets=("credible_set_id", "nunique"),
                min_credible_set_size=("credible_set_size", "min"),
                max_credible_set_size=("credible_set_size", "max"),
                total_variants_in_credible_sets=("credible_set_size", "sum"),
            )
        )
    else:
        cs_metrics = pd.DataFrame(
            columns=[
                "locus_id",
                "number_of_credible_sets",
                "min_credible_set_size",
                "max_credible_set_size",
                "total_variants_in_credible_sets",
            ]
        )

    score_col = "gene_prioritization_score_final" if "gene_prioritization_score_final" in gp.columns else "gene_prioritization_score"
    if score_col not in gp.columns:
        gp[score_col] = np.nan
    gp[score_col] = pd.to_numeric(gp[score_col], errors="coerce")

    n_genes = gp.groupby("locus_id", as_index=False).agg(number_of_genes_prioritized=("gene_symbol", "count"))
    gp_sorted = gp.sort_values(["locus_id", score_col], ascending=[True, False])
    top_gene = (
        gp_sorted.groupby("locus_id", as_index=False)
        .first()[["locus_id", "gene_symbol", score_col]]
        .rename(columns={"gene_symbol": "top_gene", score_col: "top_gene_score"})
    )
    gene_metrics = n_genes.merge(top_gene, on="locus_id", how="left")

    status_metrics = status[["locus_id", "success", "susie_converged"]].copy()
    status_metrics["diagnostics_flag"] = status["diagnostics_flag"] if "diagnostics_flag" in status.columns else np.nan

    summary = (
        pip_metrics.merge(cs_metrics, on="locus_id", how="left")
        .merge(gene_metrics, on="locus_id", how="left")
        .merge(status_metrics, on="locus_id", how="left")
    )

    for col in [
        "number_of_credible_sets",
        "min_credible_set_size",
        "max_credible_set_size",
        "total_variants_in_credible_sets",
        "number_of_genes_prioritized",
    ]:
        if col in summary.columns:
            summary[col] = pd.to_numeric(summary[col], errors="coerce")

    summary["locus_type"] = summary.apply(
        lambda r: _classify_locus(
            max_pip=float(r["max_pip"]),
            max_cs_size=float(r["max_credible_set_size"]) if pd.notna(r["max_credible_set_size"]) else np.nan,
            n_cs=int(r["number_of_credible_sets"]) if pd.notna(r["number_of_credible_sets"]) else 0,
            n_nonzero=int(r["number_of_nonzero_pip_variants"]),
        ),
        axis=1,
    )

    summary["interpretability_score"] = (
        (summary["max_pip"] * 2.0)
        - np.log(pd.to_numeric(summary["max_credible_set_size"], errors="coerce").fillna(0.0) + 1.0)
        - (pd.to_numeric(summary["number_of_credible_sets"], errors="coerce").fillna(0.0) * 0.5)
    )

    summary = summary.sort_values(["interpretability_score", "max_pip"], ascending=[False, False]).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)

    output = summary[
        [
            "locus_id",
            "max_pip",
            "mean_pip",
            "number_of_nonzero_pip_variants",
            "number_of_high_pip_variants_gt_01",
            "number_of_high_pip_variants_gt_05",
            "number_of_credible_sets",
            "min_credible_set_size",
            "max_credible_set_size",
            "total_variants_in_credible_sets",
            "top_gene",
            "top_gene_score",
            "number_of_genes_prioritized",
            "success",
            "susie_converged",
            "diagnostics_flag",
            "locus_type",
            "interpretability_score",
            "rank",
        ]
    ].copy()

    # Required concise output schema requested by user.
    output["number_of_high_pip_variants"] = output["number_of_high_pip_variants_gt_01"]

    out_tsv = root / "results/reports/locus_ranking_summary.tsv"
    out_md = root / "results/reports/locus_ranking_summary.md"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(out_tsv, sep="\t", index=False)

    top_strong = output[output["locus_type"] == "strong"].head(5)
    top_moderate = output[output["locus_type"] == "moderate"].head(5)
    top_weak = output[output["locus_type"] == "weak"].head(5)

    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Locus Ranking Summary\n\n")
        f.write("This table ranks loci by a simple interpretability heuristic combining high PIP and compact credible sets.\n")
        f.write("Use this as a presentation aid, not as a causal claim.\n\n")
        f.write("## Top 5 Strong Loci\n\n")
        if len(top_strong):
            f.write(_to_markdown_table(top_strong))
            f.write("\n\n")
        else:
            f.write("No strong loci in current outputs.\n\n")

        f.write("## Top 5 Moderate Loci\n\n")
        if len(top_moderate):
            f.write(_to_markdown_table(top_moderate))
            f.write("\n\n")
        else:
            f.write("No moderate loci in current outputs.\n\n")

        f.write("## Top 5 Weak Loci\n\n")
        if len(top_weak):
            f.write(_to_markdown_table(top_weak))
            f.write("\n\n")
        else:
            f.write("No weak loci in current outputs.\n\n")

    class_counts = output["locus_type"].value_counts().to_dict()
    top3 = output.head(3)[["rank", "locus_id", "interpretability_score", "max_pip", "locus_type"]]

    print(f"total_loci_analyzed={len(output)}")
    print(f"class_counts={class_counts}")
    print("top_3_loci_overall:")
    print(top3.to_string(index=False))
    print(f"wrote_tsv={out_tsv}")
    print(f"wrote_md={out_md}")


if __name__ == "__main__":
    main()
