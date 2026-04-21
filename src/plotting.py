"""Matplotlib plotting utilities for GWAS summary statistics."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import polars as pl


def _save_plot(fig: plt.Figure, base_path: Path, dpi: int) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(base_path.with_suffix(".png"), dpi=dpi)
    fig.savefig(base_path.with_suffix(".pdf"), dpi=dpi)
    plt.close(fig)


def _minus_log10(p: np.ndarray) -> np.ndarray:
    clipped = np.clip(p, 1e-300, 1.0)
    return -np.log10(clipped)


def _manhattan_frame(df: pl.DataFrame, p_col: str) -> pl.DataFrame:
    chr_max = df.group_by("CHR").agg(pl.max("BP").alias("chr_max")).sort("CHR")
    offsets: dict[int, int] = {}
    running = 0
    for row in chr_max.iter_rows(named=True):
        offsets[int(row["CHR"])] = running
        running += int(row["chr_max"])
    return df.with_columns(
        (pl.col("BP") + pl.col("CHR").replace(offsets, default=0)).alias("BP_CUM"),
        (-pl.col(p_col).clip(1e-300, 1.0).log10()).alias("NEG_LOG10_P"),
    )


def plot_manhattan(df: pl.DataFrame, settings: dict[str, Any], output_prefix: str = "manhattan") -> None:
    """Create Manhattan plot using cumulative coordinates across chromosomes."""
    cfg = settings["analysis"]
    p_col = cfg["selected_pvalue_column"]
    plot_df = _manhattan_frame(df, p_col)
    if cfg["downsample_for_manhattan"]:
        sig = plot_df.filter(pl.col(p_col) < cfg["suggestive_significance_threshold"])
        nonsig = plot_df.filter(pl.col(p_col) >= cfg["suggestive_significance_threshold"])
        if nonsig.height > cfg["manhattan_non_sig_max_points"]:
            nonsig = nonsig.sample(n=cfg["manhattan_non_sig_max_points"], seed=42)
        plot_df = pl.concat([sig, nonsig])
    plot_df = plot_df.sort(["CHR", "BP"])
    x = plot_df["BP_CUM"].to_numpy()
    y = plot_df["NEG_LOG10_P"].to_numpy()
    chrs = plot_df["CHR"].to_numpy()
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.scatter(x, y, c=(chrs % 2), s=3, cmap="coolwarm", alpha=0.65, linewidths=0)
    ax.axhline(-np.log10(cfg["genome_wide_significance_threshold"]), color="red", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(cfg["suggestive_significance_threshold"]), color="gray", linestyle=":", linewidth=1)
    ticks = plot_df.group_by("CHR").agg(pl.mean("BP_CUM").alias("center")).sort("CHR")
    ax.set_xticks(ticks["center"].to_numpy())
    ax.set_xticklabels([str(v) for v in ticks["CHR"].to_list()])
    ax.set_xlabel("Chromosome (cumulative coordinate)")
    ax.set_ylabel("-log10(p)")
    ax.set_title("Manhattan plot: CKDGen eGFR GWAS")
    _save_plot(fig, Path("results/plots") / output_prefix, settings["plotting"]["dpi"])


def plot_qq(df: pl.DataFrame, settings: dict[str, Any], output_prefix: str = "qq_plot") -> None:
    """Create QQ plot comparing observed and expected GWAS p-value distributions."""
    p_col = settings["analysis"]["selected_pvalue_column"]
    pvals = np.sort(df[p_col].to_numpy())
    n = pvals.size
    expected = _minus_log10(np.arange(1, n + 1) / (n + 1))
    observed = _minus_log10(pvals)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=5, alpha=0.6, color="navy", linewidths=0)
    max_val = max(expected.max(), observed.max())
    ax.plot([0, max_val], [0, max_val], linestyle="--", color="red")
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("Q-Q plot: CKDGen eGFR GWAS")
    _save_plot(fig, Path("results/plots") / output_prefix, settings["plotting"]["dpi"])


def plot_per_chromosome_counts(per_chr: pl.DataFrame, settings: dict[str, Any]) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(per_chr["CHR"].to_numpy().astype(int), per_chr["variant_count"].to_numpy(), color="slategray")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Variant count")
    ax.set_title("Variants per chromosome")
    _save_plot(fig, Path("results/plots/variants_per_chromosome"), settings["plotting"]["dpi"])


def plot_effect_vs_significance(df: pl.DataFrame, settings: dict[str, Any]) -> None:
    if "BETA" not in df.columns:
        return
    p_col = settings["analysis"]["selected_pvalue_column"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df["BETA"].to_numpy(), _minus_log10(df[p_col].to_numpy()), s=4, alpha=0.5, color="teal")
    ax.set_xlabel("Effect size (BETA)")
    ax.set_ylabel("-log10(p)")
    ax.set_title("Effect size vs statistical significance")
    _save_plot(fig, Path("results/plots/effect_vs_significance"), settings["plotting"]["dpi"])


def plot_beta_distribution(df: pl.DataFrame, settings: dict[str, Any]) -> None:
    if "BETA" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["BETA"].drop_nulls().to_numpy(), bins=80, color="steelblue", alpha=0.85)
    ax.set_xlabel("BETA")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of GWAS effect sizes")
    _save_plot(fig, Path("results/plots/beta_distribution"), settings["plotting"]["dpi"])


def plot_sample_size_distribution(df: pl.DataFrame, settings: dict[str, Any]) -> None:
    if "N" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["N"].drop_nulls().to_numpy(), bins=60, color="darkolivegreen", alpha=0.8)
    ax.set_xlabel("Sample size (N)")
    ax.set_ylabel("Count")
    ax.set_title("Sample size distribution")
    _save_plot(fig, Path("results/plots/sample_size_distribution"), settings["plotting"]["dpi"])
