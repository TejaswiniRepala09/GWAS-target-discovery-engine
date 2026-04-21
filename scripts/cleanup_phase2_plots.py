#!/usr/bin/env python3
"""Clean up multi-locus plot outputs into final vs debug buckets.

Rules:
- Keep final plots:
  pip, credible_sets, gene_prioritization, top_genes, consequence, locus, ld_heatmap
- Keep corrected diagnostics:
  corrected_ld_eigenvalues, corrected_residuals_plot
- Move all other per-locus/intermediate plots to results/debug/<locus_id>/
- Copy kept final/corrected plots to results/phase2_final/plots/<locus_id>/
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

import pandas as pd


FINAL_SUFFIXES = {
    "pip_plot.png",
    "credible_sets_plot.png",
    "gene_prioritization_plot.png",
    "top_genes_barplot.png",
    "consequence_barplot.png",
    "locus_plot.png",
    "ld_heatmap.png",
}
CORRECTED_DIAG_SUFFIXES = {
    "corrected_ld_eigenvalues.png",
    "corrected_residuals_plot.png",
}
KEEP_SUFFIXES = FINAL_SUFFIXES | CORRECTED_DIAG_SUFFIXES


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def parse_locus_and_suffix(filename: str) -> tuple[str | None, str | None]:
    """Infer locus_id and plot suffix from plot filename."""
    m = re.match(r"^(chr[^_]+_locus\d+(?:_corrected)?)_(.+\.png)$", filename)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def collect_plot_files(repo_root: Path) -> list[Path]:
    paths: list[Path] = []
    paths.extend(sorted((repo_root / "results" / "plots" / "phase2").glob("*.png")))
    paths.extend(sorted((repo_root / "results" / "loci").glob("*/plots/*.png")))
    # de-duplicate by absolute path
    seen = set()
    unique: list[Path] = []
    for p in paths:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique.append(p)
    return unique


def main() -> None:
    repo_root = Path.cwd()
    debug_root = repo_root / "results" / "debug"
    final_root = repo_root / "results" / "phase2_final" / "plots"
    report_dir = repo_root / "results" / "reports"
    ensure_dir(debug_root)
    ensure_dir(final_root)
    ensure_dir(report_dir)

    plot_files = collect_plot_files(repo_root)
    moves: list[dict[str, str]] = []
    keeps: list[dict[str, str]] = []
    copies: list[dict[str, str]] = []

    print(">>> Classifying and organizing plots")
    for src in plot_files:
        locus_id, suffix = parse_locus_and_suffix(src.name)
        if not locus_id or not suffix:
            continue

        # keep/copy final + corrected diagnostics
        if suffix in KEEP_SUFFIXES:
            final_dst = final_root / locus_id / src.name
            ensure_dir(final_dst.parent)
            shutil.copy2(src, final_dst)
            keeps.append(
                {
                    "locus_id": locus_id,
                    "plot_name": src.name,
                    "source_path": str(src),
                    "classification": "keep",
                }
            )
            copies.append(
                {
                    "locus_id": locus_id,
                    "plot_name": src.name,
                    "target_path": str(final_dst),
                    "copy_type": "phase2_final",
                }
            )
            continue

        # move intermediate plots only from per-locus folder or global phase2 folder
        if "/results/loci/" in str(src) or "/results/plots/phase2/" in str(src):
            dbg_dst = debug_root / locus_id / src.name
            ensure_dir(dbg_dst.parent)
            print(f">>> moving {src} -> {dbg_dst}")
            shutil.move(str(src), str(dbg_dst))
            moves.append(
                {
                    "locus_id": locus_id,
                    "plot_name": src.name,
                    "source_path": str(src),
                    "target_path": str(dbg_dst),
                    "classification": "moved_to_debug",
                }
            )

    moves_df = pd.DataFrame(moves)
    keeps_df = pd.DataFrame(keeps)
    copies_df = pd.DataFrame(copies)

    if len(moves_df):
        moves_df.to_csv(report_dir / "phase2_plot_moves.tsv", sep="\t", index=False)
    if len(keeps_df):
        keeps_df.to_csv(report_dir / "phase2_plot_kept.tsv", sep="\t", index=False)
    if len(copies_df):
        copies_df.to_csv(report_dir / "phase2_plot_final_index.tsv", sep="\t", index=False)

    summary = {
        "plots_scanned": len(plot_files),
        "plots_moved_to_debug": int(len(moves_df)),
        "plots_kept": int(len(keeps_df)),
        "plots_copied_to_phase2_final": int(len(copies_df)),
        "debug_root": str(debug_root),
        "phase2_final_plots_root": str(final_root),
        "reports": [
            str(report_dir / "phase2_plot_moves.tsv"),
            str(report_dir / "phase2_plot_kept.tsv"),
            str(report_dir / "phase2_plot_final_index.tsv"),
        ],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

