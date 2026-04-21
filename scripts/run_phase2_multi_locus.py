#!/usr/bin/env python3
"""Controlled multi-locus Phase 2 batch runner (pilot logic generalized)."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import polars as pl
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(REPO_ROOT / "results" / "logs" / "mplconfig"))

import matplotlib

from src.annotation.consequence_scoring import score_consequences
from src.annotation.vep_parser import parse_vep_output

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CHR7_VCF_NAME = "ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
VEP_DOCKER_IMAGE = "ensemblorg/ensembl-vep:release_115.2"
PLOT_SUFFIXES = [
    "consequence_barplot.png",
    "top_genes_barplot.png",
    "locus_plot.png",
    "ld_heatmap.png",
    "pip_plot.png",
    "credible_sets_plot.png",
    "gene_prioritization_plot.png",
]


@dataclass
class LocusRunResult:
    locus_id: str
    chromosome: str
    step_reached: str
    success: bool
    error_message: str
    n_variants_input: int
    n_variants_after_harmonization: int
    susie_converged: bool
    diagnostics_flag: bool


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _format_cmd(cmd: list[str]) -> str:
    return " ".join(shlex.quote(c) for c in cmd)


def run_cmd(cmd: list[str], log_file: Path, cwd: Path) -> None:
    print(f">>> {_format_cmd(cmd)}")
    with log_file.open("a", encoding="utf-8") as log:
        log.write(f"\n>>> {_format_cmd(cmd)}\n")
        proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
        if proc.stdout:
            print(proc.stdout, end="")
            log.write(proc.stdout)
        if proc.stderr:
            print(proc.stderr, end="")
            log.write(proc.stderr)
        if proc.returncode != 0:
            raise RuntimeError(f"Command failed ({proc.returncode}): {_format_cmd(cmd)}")


def _docker_image_available(image: str, cwd: Path) -> bool:
    proc = subprocess.run(["docker", "image", "inspect", image], cwd=str(cwd), capture_output=True, text=True)
    return proc.returncode == 0


def _find_chr_vcf(chr_value: str, repo_root: Path) -> Path:
    chr_clean = str(chr_value).replace("chr", "")
    candidates = [
        repo_root / "data/reference/Ld" / f"ALL.chr{chr_clean}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        repo_root / "data/reference/ld" / f"ALL.chr{chr_clean}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"Missing chr{chr_clean} 1000G VCF in expected paths: {candidates}")


def _make_vep_input(summary_df: pd.DataFrame, out_path: Path) -> None:
    required = {"CHR", "BP", "A1", "A2", "SNP"}
    missing = required - set(summary_df.columns)
    if missing:
        raise ValueError(f"summary stats missing required cols for VEP input: {sorted(missing)}")
    tmp = summary_df.copy()
    tmp["CHR"] = tmp["CHR"].astype(str).str.replace("chr", "", regex=False)
    tmp["BP"] = pd.to_numeric(tmp["BP"], errors="coerce").astype("Int64")
    tmp["A1"] = tmp["A1"].astype(str)
    tmp["A2"] = tmp["A2"].astype(str)
    tmp["SNP"] = tmp["SNP"].astype(str)
    out = pd.DataFrame(
        {
            "CHR": tmp["CHR"],
            "START": tmp["BP"],
            "END": tmp["BP"],
            "ALLELE": tmp["A1"] + "/" + tmp["A2"],
            "STRAND": "+",
            "ID": tmp["SNP"],
        }
    )
    out.to_csv(out_path, sep="\t", header=False, index=False)


def _harmonize_summary_with_bim(
    summary_path: Path,
    bim_path: Path,
    ld_path: Path,
    out_summary_path: Path,
    out_ld_path: Path,
    mismatch_path: Path,
) -> tuple[int, int]:
    summary = pd.read_csv(summary_path, sep="\t")
    req = {"SNP", "CHR", "BP", "A1", "A2"}
    missing = req - set(summary.columns)
    if missing:
        raise ValueError(f"summary stats missing required columns: {sorted(missing)}")

    bim = pd.read_csv(bim_path, sep=r"\s+", header=None, names=["CHR", "ID", "CM", "BP", "A1", "A2"])
    M = np.loadtxt(ld_path)
    if M.shape[0] != M.shape[1]:
        raise ValueError(f"LD not square: {M.shape}")
    if M.shape[0] != len(bim):
        raise ValueError(f"LD dim and BIM mismatch: {M.shape[0]} vs {len(bim)}")

    summary["CHR"] = summary["CHR"].astype(str).str.replace("chr", "", regex=False)
    bim["CHR"] = bim["CHR"].astype(str).str.replace("chr", "", regex=False)
    summary["BP"] = pd.to_numeric(summary["BP"], errors="coerce").astype("Int64")
    bim["BP"] = pd.to_numeric(bim["BP"], errors="coerce").astype("Int64")
    summary["A1u"] = summary["A1"].astype(str).str.upper()
    summary["A2u"] = summary["A2"].astype(str).str.upper()
    bim["A1u"] = bim["A1"].astype(str).str.upper()
    bim["A2u"] = bim["A2"].astype(str).str.upper()
    summary["PAIR"] = summary.apply(lambda r: "|".join(sorted([r["A1u"], r["A2u"]])), axis=1)
    bim["PAIR"] = bim.apply(lambda r: "|".join(sorted([r["A1u"], r["A2u"]])), axis=1)
    bim["LD_INDEX"] = np.arange(len(bim), dtype=int)

    pcol = "P_ANALYSIS" if "P_ANALYSIS" in summary.columns else ("P_GC" if "P_GC" in summary.columns else "P")
    summary["P_NUM"] = pd.to_numeric(summary[pcol], errors="coerce")
    summary_keyed = summary.sort_values(["P_NUM", "SNP"], na_position="last").drop_duplicates(
        subset=["CHR", "BP", "PAIR"], keep="first"
    )

    joined = summary_keyed.merge(
        bim[["CHR", "BP", "PAIR", "A1u", "A2u", "LD_INDEX"]],
        on=["CHR", "BP", "PAIR"],
        how="inner",
        suffixes=("_SUM", "_LD"),
    )
    if joined.empty:
        raise ValueError("No overlapping variants after CHR/BP/allele harmonization")

    joined["ALLELE_ORIENTATION"] = np.where(
        (joined["A1u_SUM"] == joined["A1u_LD"]) & (joined["A2u_SUM"] == joined["A2u_LD"]),
        "same",
        np.where(
            (joined["A1u_SUM"] == joined["A2u_LD"]) & (joined["A2u_SUM"] == joined["A1u_LD"]),
            "swapped",
            "other",
        ),
    )

    summary_keys = set(summary_keyed[["CHR", "BP", "PAIR"]].astype(str).agg(":".join, axis=1).tolist())
    bim_keys = set(bim[["CHR", "BP", "PAIR"]].astype(str).agg(":".join, axis=1).tolist())
    rows: list[dict[str, str]] = []
    rows += [{"type": "summary_not_in_ld", "key": x} for x in sorted(summary_keys - bim_keys)]
    rows += [{"type": "ld_not_in_summary", "key": x} for x in sorted(bim_keys - summary_keys)]
    rows += [
        {"type": "allele_orientation", "key": k, "status": s}
        for k, s in joined[["SNP", "ALLELE_ORIENTATION"]].astype(str).itertuples(index=False, name=None)
        if s != "same"
    ]
    pd.DataFrame(rows).to_csv(mismatch_path, sep="\t", index=False)

    # SNP quality filters
    def _is_ambiguous(a: str, b: str) -> bool:
        return {a, b} in ({"A", "T"}, {"C", "G"})

    joined["ambiguous"] = joined.apply(lambda r: _is_ambiguous(r["A1u_SUM"], r["A2u_SUM"]), axis=1)
    joined = joined.loc[~joined["ambiguous"]].copy()

    joined = joined.sort_values(["BP", "SNP", "LD_INDEX"]).drop_duplicates(subset=["LD_INDEX"], keep="first")
    common = [
        f"{c}:{bp}:{a1}:{a2}"
        for c, bp, a1, a2 in joined[["CHR", "BP", "A1u_SUM", "A2u_SUM"]].itertuples(index=False, name=None)
    ]
    if len(common) < 10:
        raise ValueError(f"Too few overlapping SNPs after harmonization/filtering: {len(common)}")

    ord_idx = joined["LD_INDEX"].astype(int).tolist()
    A = M[np.ix_(ord_idx, ord_idx)]
    finite_mask = np.isfinite(A).all(axis=1) & np.isfinite(A).all(axis=0)
    if not finite_mask.all():
        dropped = joined.loc[~finite_mask, ["SNP", "CHR", "BP", "A1u_SUM", "A2u_SUM"]].copy()
        dropped["reason"] = "non_finite_ld_row_or_col"
        drop_log = mismatch_path.with_name("ld_nonfinite_dropped.tsv")
        dropped.to_csv(drop_log, sep="\t", index=False)
        joined = joined.loc[finite_mask].copy()
        A = A[np.ix_(finite_mask, finite_mask)]
        common = [
            f"{c}:{bp}:{a1}:{a2}"
            for c, bp, a1, a2 in joined[["CHR", "BP", "A1u_SUM", "A2u_SUM"]].itertuples(index=False, name=None)
        ]
        if len(common) < 10:
            raise ValueError(
                f"Too few overlapping SNPs after dropping non-finite LD rows/cols: {len(common)}"
            )

    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 1.0)
    bad = ~np.isfinite(A)
    if bad.any():
        # Final guard for numerical safety before susie_rss; should be rare after row/col drop.
        A[bad] = 0.0
        np.fill_diagonal(A, 1.0)
    if np.max(np.abs(A - A.T)) > 1e-8:
        raise ValueError("LD matrix asymmetric after harmonization")

    out_ld = pd.DataFrame(A, columns=common, index=common)
    out_ld.insert(0, "SNP", common)
    out_ld.to_csv(out_ld_path, sep="\t", index=False)

    summary_cols = [c for c in summary_keyed.columns if c not in {"A1u", "A2u", "PAIR", "P_NUM"}]
    out_sum = joined[summary_cols].copy()
    out_sum["SNP"] = common
    out_sum.to_csv(out_summary_path, sep="\t", index=False)

    return int(len(summary)), int(len(out_sum))


def _regularize_ld(ld_path: Path, out_path: Path, shrink_lambda: float = 0.1) -> dict[str, float]:
    ld = pd.read_csv(ld_path, sep="\t")
    snps = ld["SNP"].astype(str).tolist()
    R = ld.drop(columns=["SNP"]).to_numpy(dtype=float)
    R = (R + R.T) / 2.0
    np.fill_diagonal(R, 1.0)

    eig0 = np.linalg.eigvalsh(R)
    min_before = float(np.min(eig0))

    w, v = np.linalg.eigh(R)
    w = np.clip(w, 1e-6, None)
    R_psd = v @ np.diag(w) @ v.T
    R_psd = (R_psd + R_psd.T) / 2.0
    d = np.sqrt(np.clip(np.diag(R_psd), 1e-12, None))
    R_psd = R_psd / np.outer(d, d)
    R_psd = (R_psd + R_psd.T) / 2.0
    np.fill_diagonal(R_psd, 1.0)

    R_reg = (1.0 - shrink_lambda) * R_psd + shrink_lambda * np.eye(R_psd.shape[0])
    R_reg = (R_reg + R_reg.T) / 2.0
    np.fill_diagonal(R_reg, 1.0)

    eig1 = np.linalg.eigvalsh(R_reg)
    out = pd.DataFrame(R_reg, columns=snps)
    out.insert(0, "SNP", snps)
    out.to_csv(out_path, sep="\t", index=False)
    return {"min_eig_before": min_before, "min_eig_after": float(np.min(eig1)), "shrink_lambda": shrink_lambda}


def _parse_diag_txt(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ": " in line and not line.startswith("["):
            k, v = line.split(": ", 1)
            out[k.strip()] = v.strip()
    return out


def _build_variant_and_gene_priority(
    summary_harmonized: pd.DataFrame,
    pip_df: pd.DataFrame,
    vep_scores: pl.DataFrame,
    locus_id: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    variants = summary_harmonized.copy()
    if "locus_id" not in variants.columns:
        variants["locus_id"] = locus_id

    pcol = "P_ANALYSIS" if "P_ANALYSIS" in variants.columns else ("P_GC" if "P_GC" in variants.columns else "P")
    variants["component_gwas_log10p"] = -np.log10(np.clip(pd.to_numeric(variants[pcol], errors="coerce"), 1e-300, 1.0))

    pip_m = pip_df[["SNP", "PIP"]].copy()
    pip_m["PIP"] = pd.to_numeric(pip_m["PIP"], errors="coerce")
    variants = variants.merge(pip_m, on="SNP", how="left")
    variants["component_pip"] = variants["PIP"].fillna(0.0)

    vep_pd = vep_scores.to_pandas() if vep_scores.height else pd.DataFrame()
    if not vep_pd.empty:
        vep_pd["CHR"] = vep_pd["Location"].astype(str).str.split(":").str[0].str.replace("chr", "", regex=False)
        vep_pd["BP"] = pd.to_numeric(vep_pd["Location"].astype(str).str.split(":").str[1].str.split("-").str[0], errors="coerce")
        cons = (
            vep_pd.groupby(["CHR", "BP"], as_index=False)
            .agg(
                consequence_severity_score=("consequence_severity_score", "max"),
                Consequence_term=("Consequence_term", "first"),
                gene_symbol=("SYMBOL", "first"),
                ensembl_gene_id=("Gene", "first"),
                Uploaded_variation=("Uploaded_variation", "first"),
            )
        )
        variants["CHR"] = variants["CHR"].astype(str).str.replace("chr", "", regex=False)
        variants["BP"] = pd.to_numeric(variants["BP"], errors="coerce")
        variants = variants.merge(cons, on=["CHR", "BP"], how="left")
    else:
        variants["consequence_severity_score"] = 0.0
        variants["Consequence_term"] = np.nan
        variants["gene_symbol"] = np.nan
        variants["ensembl_gene_id"] = np.nan
        variants["Uploaded_variation"] = np.nan

    variants["component_consequence"] = pd.to_numeric(variants["consequence_severity_score"], errors="coerce").fillna(0.0)
    variants["variant_priority_score"] = (
        variants["component_gwas_log10p"] * 0.4 + variants["component_consequence"] * 0.3 + variants["component_pip"] * 0.3
    )
    variants = variants.sort_values("variant_priority_score", ascending=False)

    g = variants.copy()
    g["gene_symbol"] = g["gene_symbol"].fillna("").astype(str)
    g["ensembl_gene_id"] = g["ensembl_gene_id"].fillna("").astype(str)
    gene = (
        g.groupby(["gene_symbol", "ensembl_gene_id"], as_index=False)
        .agg(
            strongest_pvalue=(pcol, "min"),
            max_pip=("PIP", "max"),
            best_consequence_severity=("consequence_severity_score", "max"),
            supporting_variant_count=("SNP", "nunique"),
            supporting_locus_count=("locus_id", "nunique"),
            best_variant=("SNP", "first"),
            gene_prioritization_score=("variant_priority_score", "max"),
        )
        .sort_values("gene_prioritization_score", ascending=False)
    )
    gene["prioritization_reason"] = (
        "strongest_pvalue="
        + gene["strongest_pvalue"].astype(str)
        + "; max_pip="
        + gene["max_pip"].astype(str)
        + "; best_variant="
        + gene["best_variant"].astype(str)
    )
    gene["locus_id"] = locus_id

    variants["locus_id"] = locus_id
    return variants, gene


def _plot_locus_outputs(
    locus_id: str,
    summary_df: pd.DataFrame,
    ld_path: Path,
    pip_df: pd.DataFrame,
    cs_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    vep_ann: pl.DataFrame,
    out_dir: Path,
) -> list[Path]:
    out_paths: list[Path] = []
    _ensure_dir(out_dir)
    ann = vep_ann.to_pandas() if vep_ann.height else pd.DataFrame()

    # consequence barplot
    if "Consequence_term" in ann.columns and ann["Consequence_term"].dropna().shape[0] > 0:
        cons_counts = ann["Consequence_term"].dropna().astype(str).value_counts().head(20)
        fig = plt.figure(figsize=(12, 6), dpi=200)
        cons_counts.sort_values().plot(kind="barh", color="#2a9d8f")
        plt.title(f"{locus_id}: Top Variant Consequences")
        plt.xlabel("Count")
        plt.ylabel("Consequence term")
        plt.tight_layout()
        p = out_dir / f"{locus_id}_consequence_barplot.png"
        fig.savefig(p)
        plt.close(fig)
        out_paths.append(p)

    # top genes barplot from annotations
    if "SYMBOL" in ann.columns:
        genes = ann["SYMBOL"].dropna().astype(str)
        genes = genes[(genes != "") & (genes != "-")]
        if len(genes):
            gc = genes.value_counts().head(15)
            fig = plt.figure(figsize=(12, 6), dpi=200)
            gc.sort_values().plot(kind="barh", color="#264653")
            plt.title(f"{locus_id}: Top Genes by Annotated Variant Count")
            plt.tight_layout()
            p = out_dir / f"{locus_id}_top_genes_barplot.png"
            fig.savefig(p)
            plt.close(fig)
            out_paths.append(p)

    # locus association plot
    pcol = "P_ANALYSIS" if "P_ANALYSIS" in summary_df.columns else ("P_GC" if "P_GC" in summary_df.columns else "P")
    tmp = summary_df.copy()
    tmp["BP"] = pd.to_numeric(tmp["BP"], errors="coerce")
    tmp[pcol] = pd.to_numeric(tmp[pcol], errors="coerce")
    tmp = tmp.dropna(subset=["BP", pcol])
    if len(tmp):
        tmp["minus_log10p"] = -np.log10(np.clip(tmp[pcol].to_numpy(dtype=float), 1e-300, 1.0))
        lead_idx = tmp[pcol].idxmin()
        fig = plt.figure(figsize=(12, 5), dpi=200)
        plt.scatter(tmp["BP"], tmp["minus_log10p"], s=10, alpha=0.7, color="#3a86ff")
        plt.scatter([tmp.loc[lead_idx, "BP"]], [tmp.loc[lead_idx, "minus_log10p"]], color="#d62828", s=45)
        plt.title(f"{locus_id}: Locus Association Plot")
        plt.xlabel("Genomic position (bp)")
        plt.ylabel("-log10(P)")
        plt.tight_layout()
        p = out_dir / f"{locus_id}_locus_plot.png"
        fig.savefig(p)
        plt.close(fig)
        out_paths.append(p)

    # LD heatmap
    lddf = pd.read_csv(ld_path, sep="\t")
    M = lddf.drop(columns=["SNP"]).to_numpy(dtype=float)
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
    p = out_dir / f"{locus_id}_ld_heatmap.png"
    fig.savefig(p)
    plt.close(fig)
    out_paths.append(p)

    # PIP plot
    pip = pip_df.copy()
    if len(pip):
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
        p = out_dir / f"{locus_id}_pip_plot.png"
        fig.savefig(p)
        plt.close(fig)
        out_paths.append(p)

    # credible sets plot
    if len(cs_df):
        cs_sizes = cs_df["credible_set_id"].astype(str).value_counts().sort_index()
        fig = plt.figure(figsize=(8, 5), dpi=200)
        cs_sizes.plot(kind="bar", color="#fb8500")
        plt.title(f"{locus_id}: Credible Set Sizes")
        plt.xlabel("Credible set")
        plt.ylabel("Variants")
        plt.tight_layout()
        p = out_dir / f"{locus_id}_credible_sets_plot.png"
        fig.savefig(p)
        plt.close(fig)
        out_paths.append(p)

    # gene prioritization plot
    if len(gene_df):
        score_col = "gene_prioritization_score_final" if "gene_prioritization_score_final" in gene_df.columns else "gene_prioritization_score"
        g = gene_df.copy()
        g["gene_symbol"] = g["gene_symbol"].fillna("").astype(str)
        g = g[(g["gene_symbol"] != "") & (g["gene_symbol"] != "-")]
        g = g.sort_values(score_col, ascending=False).head(15)
        if len(g):
            fig = plt.figure(figsize=(12, 6), dpi=200)
            plt.barh(g["gene_symbol"][::-1], g[score_col][::-1], color="#e76f51")
            plt.title(f"{locus_id}: Top Prioritized Genes")
            plt.xlabel(score_col)
            plt.tight_layout()
            p = out_dir / f"{locus_id}_gene_prioritization_plot.png"
            fig.savefig(p)
            plt.close(fig)
            out_paths.append(p)

    return out_paths


def _locus_dir_map(repo_root: Path, locus_id: str) -> dict[str, Path]:
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


def _copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    _ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return True


def _materialize_locus_outputs(repo_root: Path, locus_id: str, locus_log: Path) -> list[Path]:
    d = _locus_dir_map(repo_root, locus_id)
    for p in d.values():
        _ensure_dir(p)

    copied: list[Path] = []
    src_inputs = repo_root / "results/fine_mapping/loci_inputs" / locus_id
    src_ann = repo_root / "results/annotations" / locus_id
    src_susie = repo_root / "results/fine_mapping/susie_raw" / locus_id
    src_susie_reg = repo_root / "results/fine_mapping/susie_raw" / f"{locus_id}_regularized"
    src_diag = repo_root / "results/fine_mapping/diagnostics"
    src_prior = repo_root / "results/target_prioritization" / locus_id
    src_vep_raw = repo_root / "results/annotations/vep_raw"
    src_plot_global = repo_root / "results/plots/phase2"

    if src_inputs.exists():
        for p in src_inputs.glob("*"):
            if p.is_file() and _copy_if_exists(p, d["inputs"] / p.name):
                copied.append(d["inputs"] / p.name)

    if src_ann.exists():
        for p in src_ann.glob("*"):
            if p.is_file() and _copy_if_exists(p, d["annotations"] / p.name):
                copied.append(d["annotations"] / p.name)

    # VEP raw artifacts for this locus.
    for raw_name in [
        f"{locus_id}_vep_input_ensembl.tsv",
        f"vep_output_{locus_id}.tsv",
        f"vep_summary_{locus_id}.html",
    ]:
        p = src_vep_raw / raw_name
        if _copy_if_exists(p, d["annotations"] / raw_name):
            copied.append(d["annotations"] / raw_name)

    if src_susie.exists():
        for p in src_susie.glob("*"):
            if p.is_file():
                target = d["ld"] / p.name if any(
                    p.name.endswith(x)
                    for x in [".bed", ".bim", ".fam", ".ld", ".vcf.gz", ".vcf.gz.tbi", ".nosex", ".log", "region.txt", "ld_snp_mismatch.tsv", "ld_nonfinite_dropped.tsv"]
                ) else d["susie"] / p.name
                if _copy_if_exists(p, target):
                    copied.append(target)

    if src_susie_reg.exists():
        for p in src_susie_reg.glob("*"):
            if p.is_file() and _copy_if_exists(p, d["susie"] / p.name):
                copied.append(d["susie"] / p.name)

    for diag_name in [f"{locus_id}_diagnostics.txt", f"{locus_id}_regularized_diagnostics.txt"]:
        p = src_diag / diag_name
        if _copy_if_exists(p, d["susie"] / diag_name):
            copied.append(d["susie"] / diag_name)

    if src_prior.exists():
        for p in src_prior.glob("*"):
            if p.is_file() and _copy_if_exists(p, d["prioritization"] / p.name):
                copied.append(d["prioritization"] / p.name)

    # Preserve centralized plot compatibility while ensuring per-locus plot location.
    for suffix in PLOT_SUFFIXES + ["ld_eigenvalues.png", "residuals_plot.png", "regularized_ld_eigenvalues.png", "regularized_residuals_plot.png"]:
        name = f"{locus_id}_{suffix}"
        p = src_plot_global / name
        if _copy_if_exists(p, d["plots"] / name):
            copied.append(d["plots"] / name)

    if _copy_if_exists(locus_log, d["logs"] / locus_log.name):
        copied.append(d["logs"] / locus_log.name)

    return copied


def _select_loci(
    lead_loci_path: Path,
    lead_variants_path: Path,
    chromosome_filter: str,
    max_loci: int,
    mode: str,
) -> pd.DataFrame:
    loci = pd.read_csv(lead_loci_path, sep="\t")
    lead = pd.read_csv(lead_variants_path, sep="\t")

    chr_col = "CHR" if "CHR" in loci.columns else "chr"
    loci["_CHR"] = loci[chr_col].astype(str).str.replace("chr", "", regex=False)
    if chromosome_filter:
        wanted = {c.strip().replace("chr", "") for c in chromosome_filter.split(",")}
        loci = loci[loci["_CHR"].isin(wanted)].copy()

    if mode == "top_n":
        pcol = "P_ANALYSIS" if "P_ANALYSIS" in lead.columns else ("P_GC" if "P_GC" in lead.columns else "P")
        top_ids = (
            lead.assign(_p=pd.to_numeric(lead[pcol], errors="coerce"))
            .sort_values("_p")
            .drop_duplicates("locus_id")
            .head(max_loci)["locus_id"]
            .tolist()
        )
        loci = loci[loci["locus_id"].isin(top_ids)].copy()
    else:
        loci = loci.sort_values(["_CHR", "locus_id"]).head(max_loci).copy()
    return loci


def main() -> None:
    parser = argparse.ArgumentParser(description="Controlled multi-locus Phase 2 batch runner.")
    parser.add_argument("--chromosomes", default="7", help="Comma-separated chromosomes to run. Default: 7")
    parser.add_argument("--max-loci", type=int, default=10, help="Maximum loci to run. Default: 10")
    parser.add_argument("--mode", choices=["chromosome", "top_n"], default="chromosome", help="Locus selection mode.")
    parser.add_argument("--continue-on-error", action="store_true", default=True, help="Continue after locus failure.")
    parser.add_argument("--docker-image", default=VEP_DOCKER_IMAGE, help="Docker VEP image tag.")
    parser.add_argument(
        "--locus-ids",
        default="",
        help="Optional comma-separated locus IDs to run (e.g. chr7_locus10).",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    settings = _read_yaml(repo_root / "config/settings.yaml")
    vep_cfg = _read_yaml(repo_root / "config/vep.yaml")

    if settings["project"].get("assembly") != "GRCh37":
        raise ValueError(f"Batch runner requires GRCh37. Found: {settings['project'].get('assembly')}")

    required_bins = ["docker", "bcftools", "tabix", "plink", "Rscript", "python3"]
    for b in required_bins:
        if not shutil.which(b):
            raise RuntimeError(f"Missing required binary on PATH: {b}")

    lead_loci_path = repo_root / "results/tables/lead_loci.tsv"
    lead_variants_path = repo_root / "results/tables/lead_variants.tsv"
    cleaned_path = repo_root / "data/interim/ckdgen_egfr_cleaned.tsv"
    if not lead_loci_path.exists() or not lead_variants_path.exists() or not cleaned_path.exists():
        raise FileNotFoundError("Missing lead_loci.tsv, lead_variants.tsv, or cleaned GWAS file")

    loci_df = _select_loci(lead_loci_path, lead_variants_path, args.chromosomes, args.max_loci, args.mode)
    if args.locus_ids.strip():
        wanted = {x.strip() for x in args.locus_ids.split(",") if x.strip()}
        loci_df = loci_df[loci_df["locus_id"].astype(str).isin(wanted)].copy()
        if loci_df.empty:
            raise ValueError(f"No loci matched --locus-ids={sorted(wanted)}")
    cleaned = pd.read_csv(cleaned_path, sep="\t")
    cleaned["CHR"] = cleaned["CHR"].astype(str).str.replace("chr", "", regex=False)
    cleaned["BP"] = pd.to_numeric(cleaned["BP"], errors="coerce").astype("Int64")

    status_rows: list[LocusRunResult] = []
    pip_all: list[pd.DataFrame] = []
    cs_all: list[pd.DataFrame] = []
    var_all: list[pd.DataFrame] = []
    gene_all: list[pd.DataFrame] = []
    all_plot_paths: list[Path] = []

    _ensure_dir(repo_root / "results/logs")
    _ensure_dir(repo_root / "results/reports")

    # Ensure image is available up front; tolerate pull-helper issues if image already local.
    bootstrap_log = repo_root / "results/logs/multi_locus_batch.log"
    if not _docker_image_available(args.docker_image, repo_root):
        try:
            run_cmd(["docker", "pull", args.docker_image], bootstrap_log, repo_root)
        except Exception:  # noqa: BLE001
            if not _docker_image_available(args.docker_image, repo_root):
                raise RuntimeError(
                    f"Docker image unavailable and pull failed: {args.docker_image}. "
                    "Fix Docker credential helper or pre-pull the image manually."
                )

    for _, locus in loci_df.iterrows():
        locus_id = str(locus["locus_id"])
        chr_raw = str(locus["CHR"]).replace("chr", "")
        start = int(locus["locus_start"])
        end = int(locus["locus_end"])

        locus_log = repo_root / f"results/logs/{locus_id}.log"
        step = "init"
        n_in = 0
        n_h = 0
        converged = False
        diag_flag = False
        err = ""

        try:
            print(f"\n===== LOCUS {locus_id} (chr{chr_raw}:{start}-{end}) =====")
            _ensure_dir(repo_root / f"results/fine_mapping/loci_inputs/{locus_id}")
            _ensure_dir(repo_root / f"results/fine_mapping/susie_raw/{locus_id}")
            _ensure_dir(repo_root / "results/annotations/vep_raw")
            _ensure_dir(repo_root / f"results/annotations/{locus_id}")
            _ensure_dir(repo_root / f"results/target_prioritization/{locus_id}")
            locus_dirs = _locus_dir_map(repo_root, locus_id)
            for p in locus_dirs.values():
                _ensure_dir(p)

            bundle_dir = repo_root / f"results/fine_mapping/loci_inputs/{locus_id}"
            susie_dir = repo_root / f"results/fine_mapping/susie_raw/{locus_id}"
            vep_raw_out = repo_root / f"results/annotations/vep_raw/vep_output_{locus_id}.tsv"
            vep_input_path = repo_root / f"results/annotations/vep_raw/{locus_id}_vep_input_ensembl.tsv"
            region_vcf = susie_dir / f"{locus_id}.region.vcf.gz"
            bim_path = susie_dir / f"{locus_id}.bim"
            ld_raw_path = susie_dir / f"{locus_id}.ld"
            summary_raw_path = bundle_dir / "summary_stats_raw.tsv"
            summary_harm_path = bundle_dir / "summary_stats_harmonized.tsv"
            ld_harm_path = bundle_dir / "ld_matrix.tsv"
            mismatch_path = susie_dir / "ld_snp_mismatch.tsv"
            meta_path = bundle_dir / "metadata.yaml"

            # 1) prepare locus input
            step = "prepare_inputs"
            sub = cleaned[(cleaned["CHR"] == chr_raw) & (cleaned["BP"] >= start) & (cleaned["BP"] <= end)].copy()
            if sub.empty:
                raise ValueError("No cleaned GWAS variants in locus window")
            n_in = int(len(sub))
            sub.to_csv(summary_raw_path, sep="\t", index=False)
            meta = {
                "locus_id": locus_id,
                "chr": chr_raw,
                "window_start": start,
                "window_end": end,
                "assembly": "GRCh37",
                "variant_count_input": n_in,
            }
            meta_path.write_text(yaml.safe_dump(meta, sort_keys=False), encoding="utf-8")

            # 2) VEP input + Docker VEP
            step = "vep_run"
            _make_vep_input(sub, vep_input_path)
            run_cmd(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{repo_root}:/work",
                    "-v",
                    f"{repo_root / 'data/reference/vep'}:/opt/vep/.vep",
                    args.docker_image,
                    "vep",
                    "--input_file",
                    f"/work/{vep_input_path.relative_to(repo_root)}",
                    "--format",
                    "ensembl",
                    "--output_file",
                    f"/work/{vep_raw_out.relative_to(repo_root)}",
                    "--tab",
                    "--everything",
                    "--cache",
                    "--offline",
                    "--dir_cache",
                    "/opt/vep/.vep",
                    "--assembly",
                    "GRCh37",
                    "--force_overwrite",
                    "--fork",
                    "4",
                    "--stats_file",
                    f"/work/results/annotations/vep_raw/vep_summary_{locus_id}.html",
                ],
                locus_log,
                repo_root,
            )
            if not vep_raw_out.exists() or vep_raw_out.stat().st_size == 0:
                raise RuntimeError("VEP output missing/empty")

            # 3) VEP parse for locus
            step = "vep_parse"
            ann = parse_vep_output(vep_cfg, raw_output_path=vep_raw_out)
            vep_scored = score_consequences(ann, vep_cfg)
            ann_path_locus = repo_root / f"results/annotations/{locus_id}/variant_annotations.tsv"
            score_path_locus = repo_root / f"results/annotations/{locus_id}/variant_consequence_scores.tsv"
            ann.write_csv(ann_path_locus, separator="\t")
            vep_scored.write_csv(score_path_locus, separator="\t")

            # 4) region-based extraction + LD raw
            step = "ld_extract"
            chr_vcf = _find_chr_vcf(chr_raw, repo_root)
            run_cmd(
                [
                    "bcftools",
                    "view",
                    "-r",
                    f"{chr_raw}:{start}-{end}",
                    str(chr_vcf),
                    "-Oz",
                    "-o",
                    str(region_vcf),
                ],
                locus_log,
                repo_root,
            )
            run_cmd(["tabix", "-f", "-p", "vcf", str(region_vcf)], locus_log, repo_root)
            run_cmd(
                [
                    "plink",
                    "--vcf",
                    str(region_vcf),
                    "--double-id",
                    "--allow-extra-chr",
                    "--snps-only",
                    "just-acgt",
                    "--make-bed",
                    "--out",
                    str(susie_dir / locus_id),
                ],
                locus_log,
                repo_root,
            )
            run_cmd(
                [
                    "plink",
                    "--bfile",
                    str(susie_dir / locus_id),
                    "--r",
                    "square",
                    "--out",
                    str(susie_dir / locus_id),
                ],
                locus_log,
                repo_root,
            )

            # 5) harmonization
            step = "harmonize"
            n_in, n_h = _harmonize_summary_with_bim(
                summary_raw_path,
                bim_path,
                ld_raw_path,
                summary_harm_path,
                ld_harm_path,
                mismatch_path,
            )
            if n_h < 10:
                raise RuntimeError(f"Too few harmonized variants: {n_h}")

            # 6) SuSiE run
            step = "susie_run"
            run_cmd(
                [
                    "Rscript",
                    "scripts/susie_run_template.R",
                    str(summary_harm_path),
                    str(ld_harm_path),
                    str(susie_dir),
                ],
                locus_log,
                repo_root,
            )
            pip_path = susie_dir / "pip.tsv"
            cs_path = susie_dir / "credible_sets.tsv"
            fit_path = susie_dir / "susie_fit.rds"
            z_path = susie_dir / "zscores.tsv"
            for p in [pip_path, cs_path, fit_path, z_path]:
                if not p.exists() or p.stat().st_size == 0:
                    raise RuntimeError(f"Missing SuSiE output: {p}")

            # 7) diagnostics
            step = "diagnostics"
            run_cmd(
                [
                    "Rscript",
                    "scripts/run_susie_diagnostics_pilot.R",
                    locus_id,
                    str(summary_harm_path),
                    str(ld_harm_path),
                    str(fit_path),
                    str(z_path),
                    locus_id,
                ],
                locus_log,
                repo_root,
            )
            diag_path = repo_root / f"results/fine_mapping/diagnostics/{locus_id}_diagnostics.txt"
            diag = _parse_diag_txt(diag_path)
            converged = diag.get("converged", "FALSE").upper() == "TRUE"
            ld_ok = diag.get("ld_well_conditioned_flag", "FALSE").upper() == "TRUE"
            rss_ok = diag.get("rss_consistency_flag", "FALSE").upper() == "TRUE"
            diag_flag = ld_ok and rss_ok and converged

            # optional LD regularization rerun when diagnostics fail
            final_susie_dir = susie_dir
            final_summary = summary_harm_path
            final_ld = ld_harm_path
            if not diag_flag:
                step = "ld_regularization"
                ld_reg_path = bundle_dir / "ld_matrix_regularized.tsv"
                _regularize_ld(ld_harm_path, ld_reg_path, shrink_lambda=0.1)
                reg_dir = repo_root / f"results/fine_mapping/susie_raw/{locus_id}_regularized"
                _ensure_dir(reg_dir)
                run_cmd(
                    [
                        "Rscript",
                        "scripts/susie_run_template.R",
                        str(summary_harm_path),
                        str(ld_reg_path),
                        str(reg_dir),
                    ],
                    locus_log,
                    repo_root,
                )
                run_cmd(
                    [
                        "Rscript",
                        "scripts/run_susie_diagnostics_pilot.R",
                        locus_id,
                        str(summary_harm_path),
                        str(ld_reg_path),
                        str(reg_dir / "susie_fit.rds"),
                        str(reg_dir / "zscores.tsv"),
                        f"{locus_id}_regularized",
                    ],
                    locus_log,
                    repo_root,
                )
                diag2 = _parse_diag_txt(repo_root / f"results/fine_mapping/diagnostics/{locus_id}_regularized_diagnostics.txt")
                converged = diag2.get("converged", "FALSE").upper() == "TRUE"
                ld_ok = diag2.get("ld_well_conditioned_flag", "FALSE").upper() == "TRUE"
                rss_ok = diag2.get("rss_consistency_flag", "FALSE").upper() == "TRUE"
                diag_flag = ld_ok and rss_ok and converged
                final_susie_dir = reg_dir
                final_ld = ld_reg_path

            # 8) per-locus parse + prioritization + plots
            step = "prioritization"
            pip = pd.read_csv(final_susie_dir / "pip.tsv", sep="\t")
            cs = pd.read_csv(final_susie_dir / "credible_sets.tsv", sep="\t")
            if "locus_id" not in pip.columns:
                pip["locus_id"] = locus_id
            if "locus_id" not in cs.columns:
                cs["locus_id"] = locus_id
            if "posterior_rank" not in cs.columns:
                cs["posterior_rank"] = cs.groupby("credible_set_id").cumcount() + 1

            sum_h = pd.read_csv(final_summary, sep="\t")
            var_scores, gene_scores = _build_variant_and_gene_priority(sum_h, pip, vep_scored, locus_id)
            var_scores.to_csv(repo_root / f"results/target_prioritization/{locus_id}/variant_priority_scores.tsv", sep="\t", index=False)
            gene_scores.to_csv(repo_root / f"results/target_prioritization/{locus_id}/gene_prioritization.tsv", sep="\t", index=False)

            step = "plots"
            plots = _plot_locus_outputs(
                locus_id=locus_id,
                summary_df=sum_h,
                ld_path=final_ld,
                pip_df=pip,
                cs_df=cs,
                gene_df=gene_scores,
                vep_ann=ann,
                out_dir=locus_dirs["plots"],
            )
            all_plot_paths.extend(plots)

            # Keep a centralized phase2 plots folder for backward compatibility.
            _ensure_dir(repo_root / "results/plots/phase2")
            for p in plots:
                _copy_if_exists(p, repo_root / "results/plots/phase2" / p.name)

            _materialize_locus_outputs(repo_root, locus_id, locus_log)

            pip_all.append(pip)
            cs_all.append(cs)
            var_all.append(var_scores)
            gene_all.append(gene_scores)
            status_rows.append(
                LocusRunResult(
                    locus_id=locus_id,
                    chromosome=chr_raw,
                    step_reached=step,
                    success=True,
                    error_message="",
                    n_variants_input=n_in,
                    n_variants_after_harmonization=n_h,
                    susie_converged=converged,
                    diagnostics_flag=diag_flag,
                )
            )
        except Exception as e:  # noqa: BLE001
            err = f"{type(e).__name__}: {e}"
            tb = traceback.format_exc()
            with locus_log.open("a", encoding="utf-8") as log:
                log.write("\nERROR\n")
                log.write(err + "\n")
                log.write(tb + "\n")
            status_rows.append(
                LocusRunResult(
                    locus_id=locus_id,
                    chromosome=chr_raw,
                    step_reached=step,
                    success=False,
                    error_message=err,
                    n_variants_input=n_in,
                    n_variants_after_harmonization=n_h,
                    susie_converged=converged,
                    diagnostics_flag=diag_flag,
                )
            )
            print(f"[{locus_id}] FAILED at step={step}: {err}")
            _materialize_locus_outputs(repo_root, locus_id, locus_log)
            if not args.continue_on_error:
                raise

    # Aggregate successful loci
    pip_out = repo_root / "results/fine_mapping/pip_summary.tsv"
    cs_out = repo_root / "results/fine_mapping/credible_sets.tsv"
    vp_out = repo_root / "results/target_prioritization/variant_priority_scores.tsv"
    gp_out = repo_root / "results/target_prioritization/gene_prioritization.tsv"
    _ensure_dir(pip_out.parent)
    _ensure_dir(vp_out.parent)

    if pip_all:
        pd.concat(pip_all, ignore_index=True).to_csv(pip_out, sep="\t", index=False)
    if cs_all:
        pd.concat(cs_all, ignore_index=True).to_csv(cs_out, sep="\t", index=False)
    if var_all:
        pd.concat(var_all, ignore_index=True).to_csv(vp_out, sep="\t", index=False)
    if gene_all:
        pd.concat(gene_all, ignore_index=True).to_csv(gp_out, sep="\t", index=False)

    status_df = pd.DataFrame([r.__dict__ for r in status_rows])
    status_path = repo_root / "results/reports/multi_locus_status.tsv"
    status_df.to_csv(status_path, sep="\t", index=False)

    # Per-locus output directory index for quick navigation.
    index_rows = []
    for locus_id in status_df["locus_id"].astype(str).tolist():
        d = _locus_dir_map(repo_root, locus_id)
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
    pd.DataFrame(index_rows).to_csv(repo_root / "results/reports/locus_output_index.tsv", sep="\t", index=False)

    attempted = int(len(status_df))
    succeeded = int(status_df["success"].sum()) if attempted else 0
    failed = attempted - succeeded
    top_blockers = (
        status_df.loc[~status_df["success"], "error_message"].value_counts().head(5).to_dict()
        if failed
        else {}
    )

    summary_md = repo_root / "results/reports/multi_locus_summary.md"
    with summary_md.open("w", encoding="utf-8") as f:
        f.write("# Multi-Locus Batch Summary\n\n")
        f.write(f"- Attempted loci: {attempted}\n")
        f.write(f"- Succeeded loci: {succeeded}\n")
        f.write(f"- Failed loci: {failed}\n")
        f.write(f"- Chromosome filter: {args.chromosomes}\n")
        f.write(f"- Max loci: {args.max_loci}\n")
        f.write(f"- Mode: {args.mode}\n\n")
        f.write("## Aggregate Outputs\n")
        f.write(f"- `{pip_out}`\n")
        f.write(f"- `{cs_out}`\n")
        f.write(f"- `{vp_out}`\n")
        f.write(f"- `{gp_out}`\n\n")
        f.write("## Status Table\n")
        f.write(f"- `{status_path}`\n\n")
        f.write("## Per-Locus Output Index\n")
        f.write("- `results/reports/locus_output_index.tsv`\n\n")
        f.write("## Top Blockers\n")
        if top_blockers:
            for k, v in top_blockers.items():
                f.write(f"- {k} ({v})\n")
        else:
            f.write("- None\n")
        f.write("\n## Generated Plots\n")
        for p in sorted(set(all_plot_paths)):
            f.write(f"- `{p}`\n")

    print(json.dumps(
        {
            "attempted": attempted,
            "succeeded": succeeded,
            "failed": failed,
            "top_blockers": top_blockers,
            "aggregate_outputs": [str(pip_out), str(cs_out), str(vp_out), str(gp_out)],
            "plots_dir": str(repo_root / "results/plots/phase2"),
            "status_table": str(status_path),
            "summary_report": str(summary_md),
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
