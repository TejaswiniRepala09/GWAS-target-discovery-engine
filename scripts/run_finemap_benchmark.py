#!/usr/bin/env python3
"""Benchmark SuSiE outputs against FINEMAP on successful loci only.

This script does not rerun GWAS/VEP/SuSiE. It uses existing SuSiE outputs and
harmonized per-locus inputs to prepare FINEMAP contracts, run FINEMAP when
available, and aggregate comparison outputs.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
STATUS_TSV = ROOT / "results/reports/multi_locus_status.tsv"

BENCHMARK_DIR = ROOT / "results/benchmarking"
FINEMAP_INPUT_DIR = BENCHMARK_DIR / "finemap_inputs"
FINEMAP_RAW_DIR = BENCHMARK_DIR / "finemap_raw"
PLOTS_DIR = ROOT / "results/plots/benchmarking"
REPORT_PATH = ROOT / "results/reports/benchmarking_summary.md"

VARIANT_COMPARISON = BENCHMARK_DIR / "susie_vs_finemap_variant_comparison.tsv"
LOCUS_SUMMARY = BENCHMARK_DIR / "susie_vs_finemap_locus_summary.tsv"
LOCUS_STATUS = BENCHMARK_DIR / "finemap_benchmark_status.tsv"


@dataclass
class LocusBenchmark:
    locus_id: str
    benchmark_status: str
    error_message: str
    susie_top_variant: str | None
    susie_max_pip: float | None
    susie_credible_set_count: int | None
    susie_runtime_sec: float | None
    finemap_top_variant: str | None
    finemap_max_posterior: float | None
    finemap_credible_set_count: int | None
    finemap_runtime_sec: float | None
    overlap_flag: bool | None


def _print_cmd(cmd: list[str]) -> None:
    print(">>>", " ".join(cmd))


def _run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    _print_cmd(cmd)
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)


def _safe_float(v: Any) -> float | None:
    if v is None:
        return None
    try:
        out = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def _detect_col(df: pl.DataFrame, names: list[str]) -> str | None:
    cols_lower = {c.lower(): c for c in df.columns}
    for name in names:
        if name.lower() in cols_lower:
            return cols_lower[name.lower()]
    return None


def _load_successful_loci() -> list[str]:
    status = pl.read_csv(STATUS_TSV, separator="\t")
    if "success" not in status.columns:
        raise RuntimeError(f"'success' column missing in {STATUS_TSV}")
    loci = (
        status.filter(pl.col("success") == True)
        .select("locus_id")
        .to_series()
        .to_list()
    )
    return sorted(set(loci))


def _read_susie(locus_id: str) -> tuple[pl.DataFrame, int, str | None, float | None]:
    pip_path = ROOT / f"results/fine_mapping/susie_raw/{locus_id}/pip.tsv"
    cs_path = ROOT / f"results/fine_mapping/susie_raw/{locus_id}/credible_sets.tsv"
    if not pip_path.exists():
        raise FileNotFoundError(f"Missing SuSiE PIP file: {pip_path}")
    pip_df = pl.read_csv(pip_path, separator="\t")
    snp_col = _detect_col(pip_df, ["SNP", "rsid", "variant_id"])
    pip_col = _detect_col(pip_df, ["PIP", "pip", "posterior"])
    if not snp_col or not pip_col:
        raise RuntimeError(f"Missing required SuSiE columns in {pip_path}")
    pip_df = (
        pip_df.rename({snp_col: "SNP", pip_col: "susie_pip"})
        .with_columns(
            pl.col("SNP").cast(pl.Utf8),
            pl.col("susie_pip").cast(pl.Float64),
            pl.lit(locus_id).alias("locus_id"),
        )
        .select(["locus_id", "SNP", "susie_pip"])
    )
    top = (
        pip_df.sort("susie_pip", descending=True)
        .head(1)
        .to_dicts()
    )
    top_variant = top[0]["SNP"] if top else None
    max_pip = _safe_float(top[0]["susie_pip"]) if top else None
    cs_count = 0
    if cs_path.exists():
        cs_df = pl.read_csv(cs_path, separator="\t")
        cs_col = _detect_col(cs_df, ["credible_set_id", "cs_id", "cs"])
        if cs_col:
            cs_count = cs_df.select(pl.col(cs_col).cast(pl.Utf8).n_unique()).item() or 0
    return pip_df, cs_count, top_variant, max_pip


def _parse_stable_variant(s: str) -> tuple[str, int, str, str] | None:
    # expected stable format CHR:BP:REF:ALT
    p = s.split(":")
    if len(p) < 4:
        return None
    try:
        bp = int(p[1])
    except ValueError:
        return None
    return p[0], bp, p[2], p[3]


def _prepare_finemap_inputs(locus_id: str, out_dir: Path) -> tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    ss_h = ROOT / f"results/fine_mapping/loci_inputs/{locus_id}/summary_stats_harmonized.tsv"
    ss = ROOT / f"results/fine_mapping/loci_inputs/{locus_id}/summary_stats.tsv"
    ld_path = ROOT / f"results/fine_mapping/loci_inputs/{locus_id}/ld_matrix.tsv"
    summary_path = ss_h if ss_h.exists() else ss
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing summary stats for {locus_id}: {summary_path}")
    if not ld_path.exists():
        raise FileNotFoundError(f"Missing LD matrix for {locus_id}: {ld_path}")

    sdf = pl.read_csv(summary_path, separator="\t")
    snp_col = _detect_col(sdf, ["SNP", "variant_id", "MARKERNAME"])
    chr_col = _detect_col(sdf, ["CHR", "chromosome"])
    bp_col = _detect_col(sdf, ["BP", "POS", "position"])
    a1_col = _detect_col(sdf, ["A1", "ALT", "effect_allele"])
    a2_col = _detect_col(sdf, ["A2", "REF", "other_allele"])
    eaf_col = _detect_col(sdf, ["EAF", "maf", "MAF"])
    beta_col = _detect_col(sdf, ["BETA", "beta"])
    se_col = _detect_col(sdf, ["SE", "se"])
    required = [snp_col, chr_col, bp_col, a1_col, a2_col, beta_col, se_col]
    if any(c is None for c in required):
        raise RuntimeError(f"Missing FINEMAP-required columns in {summary_path}")

    z_rows: list[dict[str, Any]] = []
    for rec in sdf.select([snp_col, chr_col, bp_col, a1_col, a2_col, beta_col, se_col] + ([eaf_col] if eaf_col else [])).to_dicts():
        raw_id = str(rec[snp_col])
        parsed = _parse_stable_variant(raw_id)
        chr_v = rec[chr_col]
        bp_v = rec[bp_col]
        a1 = str(rec[a1_col]).upper()
        a2 = str(rec[a2_col]).upper()
        if parsed is not None:
            chr_text, bp_int, ref, alt = parsed
            chr_out = int(chr_text.replace("chr", "")) if chr_text.lower().startswith("chr") else int(chr_text)
            bp_out = bp_int
            # keep alleles from summary stats to preserve model consistency
            a1_out, a2_out = a1, a2
        else:
            chr_out = int(chr_v)
            bp_out = int(bp_v)
            a1_out, a2_out = a1, a2

        maf = None
        if eaf_col:
            eaf = _safe_float(rec[eaf_col])
            if eaf is not None:
                maf = min(eaf, 1 - eaf)
        z_rows.append(
            {
                "rsid": raw_id,
                "chromosome": chr_out,
                "position": bp_out,
                "allele1": a1_out,
                "allele2": a2_out,
                "maf": maf if maf is not None else 0.05,
                "beta": float(rec[beta_col]),
                "se": float(rec[se_col]),
            }
        )

    z_df = pl.DataFrame(z_rows)
    z_file = out_dir / f"{locus_id}.z"
    z_df.write_csv(z_file, separator=" ")

    ld_df = pl.read_csv(ld_path, separator="\t")
    if "SNP" in ld_df.columns:
        ld_df = ld_df.drop("SNP")
    ld_file = out_dir / f"{locus_id}.ld"
    np.savetxt(ld_file, ld_df.to_numpy(), fmt="%.10g")

    n_samples = int(_safe_float(sdf["N"][0]) or 10000) if "N" in sdf.columns and sdf.height else 10000
    master_file = out_dir / f"{locus_id}.master"
    snp_file = out_dir / f"{locus_id}.snp"
    config_file = out_dir / f"{locus_id}.config"
    cred_file = out_dir / f"{locus_id}.cred"
    log_file = out_dir / f"{locus_id}.log"
    master_file.write_text(
        "z;ld;snp;config;cred;log;n_samples\n"
        f"{z_file};{ld_file};{snp_file};{config_file};{cred_file};{log_file};{n_samples}\n"
    )
    return master_file, snp_file, cred_file


def _read_whitespace_table(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        return rows
    header = lines[0].split()
    for ln in lines[1:]:
        parts = ln.split()
        if len(parts) < len(header):
            continue
        rows.append({header[i]: parts[i] for i in range(len(header))})
    return rows


def _parse_finemap_outputs(
    locus_id: str,
    snp_file: Path,
    cred_file: Path,
) -> tuple[pl.DataFrame, str | None, float | None, int | None]:
    if not snp_file.exists() or snp_file.stat().st_size == 0:
        raise FileNotFoundError(f"Missing FINEMAP SNP output: {snp_file}")
    rows = _read_whitespace_table(snp_file)
    if not rows:
        raise RuntimeError(f"Empty FINEMAP SNP output: {snp_file}")
    # infer column names
    sample_keys = {k.lower(): k for k in rows[0].keys()}
    rsid_col = sample_keys.get("rsid") or sample_keys.get("snp") or list(rows[0].keys())[0]
    prob_col = None
    for cand in ["prob", "posterior", "post", "pp", "probability"]:
        if cand in sample_keys:
            prob_col = sample_keys[cand]
            break
    if prob_col is None:
        # fallback to any column containing "prob"
        for k in rows[0].keys():
            if "prob" in k.lower() or "post" in k.lower():
                prob_col = k
                break
    if prob_col is None:
        raise RuntimeError(f"Could not detect FINEMAP posterior column in {snp_file}")

    data = []
    for r in rows:
        pp = _safe_float(r.get(prob_col))
        if pp is None:
            continue
        data.append({"locus_id": locus_id, "SNP": str(r.get(rsid_col)), "finemap_posterior": pp})
    if not data:
        raise RuntimeError(f"No numeric FINEMAP posterior values in {snp_file}")
    finemap_df = pl.DataFrame(data)
    top = finemap_df.sort("finemap_posterior", descending=True).head(1).to_dicts()[0]
    top_variant = str(top["SNP"])
    max_pp = float(top["finemap_posterior"])

    cred_count: int | None = None
    if cred_file.exists() and cred_file.stat().st_size > 0:
        cred_rows = _read_whitespace_table(cred_file)
        cred_count = len(cred_rows) if cred_rows else 0

    return finemap_df, top_variant, max_pp, cred_count


def _plot_runtime(df: pl.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_df = df.select(["locus_id", "susie_runtime_sec", "finemap_runtime_sec"]).sort("locus_id")
    s = [v if v is not None else np.nan for v in plot_df["susie_runtime_sec"].to_list()]
    f = [v if v is not None else np.nan for v in plot_df["finemap_runtime_sec"].to_list()]
    if np.isfinite(np.array(f, dtype=float)).sum() == 0 and np.isfinite(np.array(s, dtype=float)).sum() == 0:
        ax.text(0.5, 0.5, "No runtime values available", ha="center", va="center", fontsize=12)
        ax.set_axis_off()
    else:
        x = np.arange(len(plot_df))
        ax.bar(x - 0.2, np.nan_to_num(s, nan=0.0), width=0.4, label="SuSiE")
        ax.bar(x + 0.2, np.nan_to_num(f, nan=0.0), width=0.4, label="FINEMAP")
        ax.set_xticks(x)
        ax.set_xticklabels(plot_df["locus_id"].to_list(), rotation=70, ha="right")
        ax.set_ylabel("Runtime (sec)")
        ax.set_title("SuSiE vs FINEMAP Runtime")
        ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_overlap(df: pl.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    overlap = int(df.filter(pl.col("overlap_flag") == True).height)
    no_overlap = int(df.filter(pl.col("overlap_flag") == False).height)
    unresolved = int(df.filter(pl.col("overlap_flag").is_null()).height)
    vals = [overlap, no_overlap, unresolved]
    labels = ["Top-variant overlap", "No overlap", "Unavailable"]
    ax.bar(labels, vals)
    ax.set_ylabel("Locus count")
    ax.set_title("SuSiE vs FINEMAP Top Variant Agreement")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.1, str(v), ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_credible_counts(df: pl.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_df = df.select(["locus_id", "susie_credible_set_count", "finemap_credible_set_count"]).sort("locus_id")
    x = np.arange(plot_df.height)
    s = np.array([v if v is not None else np.nan for v in plot_df["susie_credible_set_count"].to_list()], dtype=float)
    f = np.array([v if v is not None else np.nan for v in plot_df["finemap_credible_set_count"].to_list()], dtype=float)
    if np.isfinite(f).sum() == 0 and np.isfinite(s).sum() == 0:
        ax.text(0.5, 0.5, "No credible set counts available", ha="center", va="center", fontsize=12)
        ax.set_axis_off()
    else:
        ax.bar(x - 0.2, np.nan_to_num(s, nan=0.0), width=0.4, label="SuSiE")
        ax.bar(x + 0.2, np.nan_to_num(f, nan=0.0), width=0.4, label="FINEMAP")
        ax.set_xticks(x)
        ax.set_xticklabels(plot_df["locus_id"].to_list(), rotation=70, ha="right")
        ax.set_ylabel("Credible set count")
        ax.set_title("Credible Set Count Comparison")
        ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    FINEMAP_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    FINEMAP_RAW_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    loci = _load_successful_loci()
    # Allow explicit FINEMAP path override for reproducible environments.
    finemap_bin = os.getenv("FINEMAP_BIN", "").strip() or shutil.which("finemap")
    records: list[LocusBenchmark] = []
    variant_tables: list[pl.DataFrame] = []
    top_blockers: dict[str, int] = {}

    for locus_id in loci:
        print(f"\n===== BENCHMARK LOCUS {locus_id} =====")
        locus_input_dir = FINEMAP_INPUT_DIR / locus_id
        locus_raw_dir = FINEMAP_RAW_DIR / locus_id
        locus_raw_dir.mkdir(parents=True, exist_ok=True)

        try:
            susie_df, susie_cs_count, susie_top, susie_max = _read_susie(locus_id)
            master, snp_file, cred_file = _prepare_finemap_inputs(locus_id, locus_input_dir)

            finemap_df: pl.DataFrame | None = None
            finemap_top = None
            finemap_max = None
            finemap_cs_count = None
            finemap_runtime = None
            status = "success"
            err = ""

            if finemap_bin is None:
                status = "failed_missing_finemap_binary"
                err = "FINEMAP binary not found on PATH"
            else:
                cmd = [finemap_bin, "--sss", "--in-files", str(master)]
                t0 = time.perf_counter()
                proc = _run(cmd, cwd=locus_raw_dir)
                finemap_runtime = time.perf_counter() - t0
                if proc.returncode != 0:
                    status = "failed_finemap_run"
                    err = (proc.stderr or proc.stdout or "FINEMAP run failed").strip()[:500]
                else:
                    # copy outputs to raw locus dir if generated next to master
                    for p in [snp_file, cred_file, locus_input_dir / f"{locus_id}.log"]:
                        if p.exists():
                            dst = locus_raw_dir / p.name
                            if dst != p:
                                dst.write_bytes(p.read_bytes())
                    finemap_df, finemap_top, finemap_max, finemap_cs_count = _parse_finemap_outputs(
                        locus_id=locus_id,
                        snp_file=locus_raw_dir / snp_file.name if (locus_raw_dir / snp_file.name).exists() else snp_file,
                        cred_file=locus_raw_dir / cred_file.name if (locus_raw_dir / cred_file.name).exists() else cred_file,
                    )

            if finemap_df is not None:
                joined = susie_df.join(finemap_df, on=["locus_id", "SNP"], how="full")
            else:
                joined = susie_df.with_columns(pl.lit(None, dtype=pl.Float64).alias("finemap_posterior"))
            variant_tables.append(joined)

            overlap = None
            if susie_top and finemap_top:
                overlap = susie_top == finemap_top

            records.append(
                LocusBenchmark(
                    locus_id=locus_id,
                    benchmark_status=status,
                    error_message=err,
                    susie_top_variant=susie_top,
                    susie_max_pip=susie_max,
                    susie_credible_set_count=susie_cs_count,
                    susie_runtime_sec=None,
                    finemap_top_variant=finemap_top,
                    finemap_max_posterior=finemap_max,
                    finemap_credible_set_count=finemap_cs_count,
                    finemap_runtime_sec=finemap_runtime,
                    overlap_flag=overlap,
                )
            )
            if status != "success":
                top_blockers[status] = top_blockers.get(status, 0) + 1

        except Exception as exc:  # continue-on-error per locus
            msg = str(exc)[:500]
            records.append(
                LocusBenchmark(
                    locus_id=locus_id,
                    benchmark_status="failed_exception",
                    error_message=msg,
                    susie_top_variant=None,
                    susie_max_pip=None,
                    susie_credible_set_count=None,
                    susie_runtime_sec=None,
                    finemap_top_variant=None,
                    finemap_max_posterior=None,
                    finemap_credible_set_count=None,
                    finemap_runtime_sec=None,
                    overlap_flag=None,
                )
            )
            top_blockers["failed_exception"] = top_blockers.get("failed_exception", 0) + 1

    summary_df = pl.DataFrame([r.__dict__ for r in records]).with_columns(
        pl.col("locus_id").cast(pl.Utf8)
    )
    summary_df.write_csv(LOCUS_SUMMARY, separator="\t")
    summary_df.write_csv(LOCUS_STATUS, separator="\t")

    if variant_tables:
        var_df = pl.concat(variant_tables, how="diagonal_relaxed").sort(["locus_id", "SNP"])
    else:
        var_df = pl.DataFrame(
            schema={
                "locus_id": pl.Utf8,
                "SNP": pl.Utf8,
                "susie_pip": pl.Float64,
                "finemap_posterior": pl.Float64,
            }
        )
    var_df.write_csv(VARIANT_COMPARISON, separator="\t")

    _plot_runtime(summary_df, PLOTS_DIR / "susie_vs_finemap_runtime.png")
    _plot_overlap(summary_df, PLOTS_DIR / "susie_vs_finemap_top_variant_overlap.png")
    _plot_credible_counts(summary_df, PLOTS_DIR / "susie_vs_finemap_credible_set_sizes.png")

    attempted = len(loci)
    success_n = summary_df.filter(pl.col("benchmark_status") == "success").height
    failed_n = attempted - success_n
    overlap_n = summary_df.filter(pl.col("overlap_flag") == True).height
    comparable_n = summary_df.filter(pl.col("overlap_flag").is_not_null()).height

    report = []
    report.append("# Benchmarking Summary (SuSiE vs FINEMAP)")
    report.append("")
    report.append("## Scope")
    report.append(f"- Successful loci benchmark cohort from `multi_locus_status.tsv`: **{attempted}**")
    report.append("- GWAS/VEP/SuSiE were not rerun; existing outputs were reused.")
    report.append("")
    report.append("## Benchmark Outcome")
    report.append(f"- Loci attempted: **{attempted}**")
    report.append(f"- Loci benchmarked successfully (FINEMAP completed): **{success_n}**")
    report.append(f"- Loci failed: **{failed_n}**")
    report.append(f"- Top-variant agreement (among comparable loci): **{overlap_n}/{comparable_n}**")
    report.append("")
    report.append("## Method Stability")
    report.append("- Stable loci: loci with matching SuSiE and FINEMAP top variants.")
    report.append("- Method-sensitive loci: loci where top variants differ, or FINEMAP failed.")
    report.append("- Runtime comparison plotted where runtime values were available.")
    report.append("")
    report.append("## Runtime Notes")
    report.append("- SuSiE runtime was not consistently logged in a parseable per-locus format; values may be null.")
    report.append("- FINEMAP runtime is measured from command execution when FINEMAP runs.")
    report.append("")
    report.append("## Top Blockers")
    if top_blockers:
        for k, v in sorted(top_blockers.items(), key=lambda x: (-x[1], x[0])):
            report.append(f"- {k}: {v}")
    else:
        report.append("- None")
    report.append("")
    report.append("## Output Files")
    report.append(f"- `{VARIANT_COMPARISON.relative_to(ROOT)}`")
    report.append(f"- `{LOCUS_SUMMARY.relative_to(ROOT)}`")
    report.append(f"- `{PLOTS_DIR.relative_to(ROOT) / 'susie_vs_finemap_runtime.png'}`")
    report.append(f"- `{PLOTS_DIR.relative_to(ROOT) / 'susie_vs_finemap_top_variant_overlap.png'}`")
    report.append(f"- `{PLOTS_DIR.relative_to(ROOT) / 'susie_vs_finemap_credible_set_sizes.png'}`")
    REPORT_PATH.write_text("\n".join(report))

    print("\n===== BENCHMARK COMPLETE =====")
    print(f"attempted={attempted}")
    print(f"succeeded={success_n}")
    print(f"failed={failed_n}")
    print(f"top_blockers={top_blockers}")
    print(f"variant_comparison={VARIANT_COMPARISON}")
    print(f"locus_summary={LOCUS_SUMMARY}")
    print(f"report={REPORT_PATH}")


if __name__ == "__main__":
    main()
