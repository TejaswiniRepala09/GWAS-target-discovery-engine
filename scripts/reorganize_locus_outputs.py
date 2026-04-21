#!/usr/bin/env python3
"""Backfill per-locus output folders from existing multi-locus outputs.

This script does not rerun scientific analysis. It copies existing artifacts into:
results/loci/<locus_id>/{inputs,annotations,ld,susie,prioritization,plots,logs}
and regenerates only missing core plots when possible.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


PLOT_SUFFIXES = [
    "consequence_barplot.png",
    "top_genes_barplot.png",
    "locus_plot.png",
    "ld_heatmap.png",
    "pip_plot.png",
    "credible_sets_plot.png",
    "gene_prioritization_plot.png",
]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def locus_dir_map(repo_root: Path, locus_id: str) -> dict[str, Path]:
    root = repo_root / "results" / "loci" / locus_id
    return {
        "root": root,
        "inputs": root / "inputs",
        "annotations": root / "annotations",
        "ld": root / "ld",
        "susie": root / "susie",
        "prioritization": root / "prioritization",
        "plots": root / "plots",
        "logs": root / "logs",
    }


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return True


def plot_from_existing_data(locus_id: str, repo_root: Path, out_dir: Path) -> list[Path]:
    """Regenerate only core missing per-locus plots from already-computed outputs."""
    ensure_dir(out_dir)
    out_paths: list[Path] = []

    ann_path = repo_root / f"results/annotations/{locus_id}/variant_annotations.tsv"
    sum_path = repo_root / f"results/fine_mapping/loci_inputs/{locus_id}/summary_stats_harmonized.tsv"
    ld_path = repo_root / f"results/fine_mapping/loci_inputs/{locus_id}/ld_matrix.tsv"
    pip_path = repo_root / f"results/fine_mapping/susie_raw/{locus_id}_regularized/pip.tsv"
    cs_path = repo_root / f"results/fine_mapping/susie_raw/{locus_id}_regularized/credible_sets.tsv"
    if not pip_path.exists():
        pip_path = repo_root / f"results/fine_mapping/susie_raw/{locus_id}/pip.tsv"
    if not cs_path.exists():
        cs_path = repo_root / f"results/fine_mapping/susie_raw/{locus_id}/credible_sets.tsv"
    gene_path = repo_root / f"results/target_prioritization/{locus_id}/gene_prioritization.tsv"

    ann = pd.read_csv(ann_path, sep="\t") if ann_path.exists() else pd.DataFrame()
    summary_df = pd.read_csv(sum_path, sep="\t") if sum_path.exists() else pd.DataFrame()
    pip_df = pd.read_csv(pip_path, sep="\t") if pip_path.exists() else pd.DataFrame()
    cs_df = pd.read_csv(cs_path, sep="\t") if cs_path.exists() else pd.DataFrame()
    gene_df = pd.read_csv(gene_path, sep="\t") if gene_path.exists() else pd.DataFrame()

    if "Consequence_term" in ann.columns and ann["Consequence_term"].dropna().shape[0] > 0:
        p = out_dir / f"{locus_id}_consequence_barplot.png"
        if not p.exists():
            counts = ann["Consequence_term"].dropna().astype(str).value_counts().head(20)
            fig = plt.figure(figsize=(12, 6), dpi=200)
            counts.sort_values().plot(kind="barh", color="#2a9d8f")
            plt.title(f"{locus_id}: Top Variant Consequences")
            plt.tight_layout()
            fig.savefig(p)
            plt.close(fig)
            out_paths.append(p)

    if "SYMBOL" in ann.columns:
        genes = ann["SYMBOL"].dropna().astype(str)
        genes = genes[(genes != "") & (genes != "-")]
        if len(genes):
            p = out_dir / f"{locus_id}_top_genes_barplot.png"
            if not p.exists():
                gc = genes.value_counts().head(15)
                fig = plt.figure(figsize=(12, 6), dpi=200)
                gc.sort_values().plot(kind="barh", color="#264653")
                plt.title(f"{locus_id}: Top Genes by Annotated Variant Count")
                plt.tight_layout()
                fig.savefig(p)
                plt.close(fig)
                out_paths.append(p)

    if not summary_df.empty and "BP" in summary_df.columns:
        pcol = "P_ANALYSIS" if "P_ANALYSIS" in summary_df.columns else ("P_GC" if "P_GC" in summary_df.columns else "P")
        tmp = summary_df.copy()
        tmp["BP"] = pd.to_numeric(tmp["BP"], errors="coerce")
        tmp[pcol] = pd.to_numeric(tmp[pcol], errors="coerce")
        tmp = tmp.dropna(subset=["BP", pcol])
        if len(tmp):
            p = out_dir / f"{locus_id}_locus_plot.png"
            if not p.exists():
                tmp["minus_log10p"] = -np.log10(np.clip(tmp[pcol].to_numpy(dtype=float), 1e-300, 1.0))
                lead_idx = tmp[pcol].idxmin()
                fig = plt.figure(figsize=(12, 5), dpi=200)
                plt.scatter(tmp["BP"], tmp["minus_log10p"], s=10, alpha=0.7, color="#3a86ff")
                plt.scatter([tmp.loc[lead_idx, "BP"]], [tmp.loc[lead_idx, "minus_log10p"]], color="#d62828", s=45)
                plt.title(f"{locus_id}: Locus Association Plot")
                plt.xlabel("Genomic position (bp)")
                plt.ylabel("-log10(P)")
                plt.tight_layout()
                fig.savefig(p)
                plt.close(fig)
                out_paths.append(p)

    if ld_path.exists():
        p = out_dir / f"{locus_id}_ld_heatmap.png"
        if not p.exists():
            lddf = pd.read_csv(ld_path, sep="\t")
            if "SNP" in lddf.columns:
                M = lddf.drop(columns=["SNP"]).to_numpy(dtype=float)
            else:
                # Some earlier files are saved as pure numeric square matrices.
                M = lddf.to_numpy(dtype=float)
            max_n = 600
            if M.shape[0] > max_n:
                idx = np.linspace(0, M.shape[0] - 1, max_n).astype(int)
                Mplot = M[np.ix_(idx, idx)]
            else:
                Mplot = M
            fig = plt.figure(figsize=(8, 7), dpi=220)
            im = plt.imshow(Mplot, cmap="viridis", vmin=-1, vmax=1, interpolation="nearest", aspect="auto")
            plt.colorbar(im, fraction=0.046, pad=0.04, label="LD (r)")
            plt.title(f"{locus_id}: LD Heatmap")
            plt.tight_layout()
            fig.savefig(p)
            plt.close(fig)
            out_paths.append(p)

    if not pip_df.empty and "SNP" in pip_df.columns and "PIP" in pip_df.columns:
        p = out_dir / f"{locus_id}_pip_plot.png"
        if not p.exists():
            pip = pip_df.copy()
            pip["BP"] = pd.to_numeric(pip["SNP"].astype(str).str.split(":").str[1], errors="coerce")
            if pip["BP"].isna().all():
                pip["BP"] = np.arange(len(pip))
            pip["PIP"] = pd.to_numeric(pip["PIP"], errors="coerce").fillna(0.0)
            pip = pip.sort_values("BP")
            fig = plt.figure(figsize=(12, 5), dpi=200)
            plt.scatter(pip["BP"], pip["PIP"], s=10, alpha=0.7, color="#8338ec")
            top = pip.nlargest(5, "PIP")
            plt.scatter(top["BP"], top["PIP"], color="#ff006e", s=30)
            plt.title(f"{locus_id}: SuSiE PIP by Position")
            plt.xlabel("Genomic position")
            plt.ylabel("PIP")
            plt.ylim(-0.02, 1.02)
            plt.tight_layout()
            fig.savefig(p)
            plt.close(fig)
            out_paths.append(p)

    if not cs_df.empty and "credible_set_id" in cs_df.columns:
        p = out_dir / f"{locus_id}_credible_sets_plot.png"
        if not p.exists():
            sizes = cs_df["credible_set_id"].astype(str).value_counts().sort_index()
            fig = plt.figure(figsize=(8, 5), dpi=200)
            sizes.plot(kind="bar", color="#fb8500")
            plt.title(f"{locus_id}: Credible Set Sizes")
            plt.tight_layout()
            fig.savefig(p)
            plt.close(fig)
            out_paths.append(p)

    if not gene_df.empty and "gene_symbol" in gene_df.columns:
        score_col = "gene_prioritization_score_final" if "gene_prioritization_score_final" in gene_df.columns else "gene_prioritization_score"
        if score_col in gene_df.columns:
            p = out_dir / f"{locus_id}_gene_prioritization_plot.png"
            if not p.exists():
                g = gene_df.copy()
                g["gene_symbol"] = g["gene_symbol"].fillna("").astype(str)
                g = g[(g["gene_symbol"] != "") & (g["gene_symbol"] != "-")].sort_values(score_col, ascending=False).head(15)
                if len(g):
                    fig = plt.figure(figsize=(12, 6), dpi=200)
                    plt.barh(g["gene_symbol"][::-1], g[score_col][::-1], color="#e76f51")
                    plt.title(f"{locus_id}: Top Prioritized Genes")
                    plt.xlabel(score_col)
                    plt.tight_layout()
                    fig.savefig(p)
                    plt.close(fig)
                    out_paths.append(p)

    return out_paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Reorganize existing outputs into per-locus folder structure.")
    parser.add_argument("--locus-prefix", default="", help="Optional prefix filter (e.g. chr7_).")
    parser.add_argument("--successful-only", action="store_true", default=True, help="Use only success=True loci from status table.")
    parser.add_argument(
        "--include-existing-folders",
        action="store_true",
        default=True,
        help="Also discover loci from existing per-locus output folders (annotations/prioritization/fine_mapping).",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    status_path = repo_root / "results/reports/multi_locus_status.tsv"
    if not status_path.exists():
        raise FileNotFoundError(f"Missing status file: {status_path}")

    status = pd.read_csv(status_path, sep="\t")
    if args.successful_only and "success" in status.columns:
        status = status[status["success"] == True].copy()  # noqa: E712
    loci_from_status = set(status["locus_id"].astype(str).unique().tolist())

    loci_existing: set[str] = set()
    if args.include_existing_folders:
        for base in [
            repo_root / "results/annotations",
            repo_root / "results/target_prioritization",
            repo_root / "results/fine_mapping/loci_inputs",
            repo_root / "results/fine_mapping/susie_raw",
        ]:
            if base.exists():
                for p in base.iterdir():
                    if p.is_dir() and "locus" in p.name and not p.name.endswith("_regularized"):
                        loci_existing.add(p.name)

    loci = sorted(loci_from_status.union(loci_existing))
    if args.locus_prefix:
        loci = [x for x in loci if x.startswith(args.locus_prefix)]

    copied_files: list[str] = []
    found_plots: list[str] = []
    regenerated_plots: list[str] = []

    for locus_id in loci:
        d = locus_dir_map(repo_root, locus_id)
        for p in d.values():
            ensure_dir(p)

        src_inputs = repo_root / "results/fine_mapping/loci_inputs" / locus_id
        src_ann = repo_root / "results/annotations" / locus_id
        src_susie = repo_root / "results/fine_mapping/susie_raw" / locus_id
        src_susie_reg = repo_root / "results/fine_mapping/susie_raw" / f"{locus_id}_regularized"
        src_diag = repo_root / "results/fine_mapping/diagnostics"
        src_prior = repo_root / "results/target_prioritization" / locus_id
        src_vep_raw = repo_root / "results/annotations/vep_raw"
        src_plot = repo_root / "results/plots/phase2"
        src_log = repo_root / "results/logs" / f"{locus_id}.log"

        if src_inputs.exists():
            for p in src_inputs.glob("*"):
                if p.is_file() and copy_if_exists(p, d["inputs"] / p.name):
                    copied_files.append(str(d["inputs"] / p.name))

        if src_ann.exists():
            for p in src_ann.glob("*"):
                if p.is_file() and copy_if_exists(p, d["annotations"] / p.name):
                    copied_files.append(str(d["annotations"] / p.name))

        for raw_name in [f"{locus_id}_vep_input_ensembl.tsv", f"vep_output_{locus_id}.tsv", f"vep_summary_{locus_id}.html"]:
            p = src_vep_raw / raw_name
            if copy_if_exists(p, d["annotations"] / raw_name):
                copied_files.append(str(d["annotations"] / raw_name))

        if src_susie.exists():
            for p in src_susie.glob("*"):
                if p.is_file():
                    target = d["ld"] / p.name if any(
                        p.name.endswith(x)
                        for x in [".bed", ".bim", ".fam", ".ld", ".vcf.gz", ".vcf.gz.tbi", ".nosex", ".log", "region.txt", "ld_snp_mismatch.tsv", "ld_nonfinite_dropped.tsv"]
                    ) else d["susie"] / p.name
                    if copy_if_exists(p, target):
                        copied_files.append(str(target))

        if src_susie_reg.exists():
            for p in src_susie_reg.glob("*"):
                if p.is_file() and copy_if_exists(p, d["susie"] / p.name):
                    copied_files.append(str(d["susie"] / p.name))

        for diag_name in [f"{locus_id}_diagnostics.txt", f"{locus_id}_regularized_diagnostics.txt"]:
            p = src_diag / diag_name
            if copy_if_exists(p, d["susie"] / diag_name):
                copied_files.append(str(d["susie"] / diag_name))

        if src_prior.exists():
            for p in src_prior.glob("*"):
                if p.is_file() and copy_if_exists(p, d["prioritization"] / p.name):
                    copied_files.append(str(d["prioritization"] / p.name))

        if copy_if_exists(src_log, d["logs"] / src_log.name):
            copied_files.append(str(d["logs"] / src_log.name))

        for suffix in PLOT_SUFFIXES + ["ld_eigenvalues.png", "residuals_plot.png", "regularized_ld_eigenvalues.png", "regularized_residuals_plot.png"]:
            plot_name = f"{locus_id}_{suffix}"
            if copy_if_exists(src_plot / plot_name, d["plots"] / plot_name):
                found_plots.append(str(d["plots"] / plot_name))

        missing_core = [d["plots"] / f"{locus_id}_{suffix}" for suffix in PLOT_SUFFIXES if not (d["plots"] / f"{locus_id}_{suffix}").exists()]
        if missing_core:
            regen = plot_from_existing_data(locus_id, repo_root, d["plots"])
            regenerated_plots.extend(str(p) for p in regen)

    index_rows = []
    for locus_id in loci:
        d = locus_dir_map(repo_root, locus_id)
        index_rows.append(
            {
                "locus_id": locus_id,
                "inputs_dir": str(d["inputs"]),
                "annotations_dir": str(d["annotations"]),
                "ld_dir": str(d["ld"]),
                "susie_dir": str(d["susie"]),
                "prioritization_dir": str(d["prioritization"]),
                "plots_dir": str(d["plots"]),
                "logs_dir": str(d["logs"]),
            }
        )
    index_path = repo_root / "results/reports/locus_output_index.tsv"
    pd.DataFrame(index_rows).to_csv(index_path, sep="\t", index=False)

    print(
        json.dumps(
            {
                "loci_processed": len(loci),
                "copied_file_count": len(set(copied_files)),
                "found_plot_count": len(set(found_plots)),
                "regenerated_plot_count": len(set(regenerated_plots)),
                "locus_output_index": str(index_path),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
