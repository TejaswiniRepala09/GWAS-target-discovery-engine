#!/usr/bin/env python3
"""Find high-confidence coding loci for AlphaFold-style structural follow-up."""

from __future__ import annotations

from pathlib import Path
import re

import pandas as pd

TARGET_TERMS = {"missense_variant", "protein_altering_variant", "stop_gained"}


def is_target_consequence(term: str) -> bool:
    if pd.isna(term):
        return False
    t = str(term).strip()
    if not t:
        return False
    parts = re.split(r"[,&;|]", t)
    return any(p.strip() in TARGET_TERMS for p in parts)


def pick_pip_path(root: Path, locus_id: str) -> Path | None:
    reg = root / f"results/fine_mapping/susie_raw/{locus_id}_regularized/pip.tsv"
    raw = root / f"results/fine_mapping/susie_raw/{locus_id}/pip.tsv"
    if reg.exists():
        return reg
    if raw.exists():
        return raw
    return None


def main() -> None:
    root = Path.cwd()
    ann_root = root / "results/annotations"
    if not ann_root.exists():
        raise FileNotFoundError(f"Missing annotation root: {ann_root}")

    rows: list[pd.DataFrame] = []

    loci = sorted([p.name for p in ann_root.glob("chr*_locus*") if p.is_dir()])
    if not loci:
        raise RuntimeError("No per-locus annotation directories found under results/annotations/")

    for locus_id in loci:
        ann_path = ann_root / locus_id / "variant_annotations.tsv"
        pip_path = pick_pip_path(root, locus_id)
        if not ann_path.exists() or pip_path is None:
            continue

        ann = pd.read_csv(ann_path, sep="\t")
        pip = pd.read_csv(pip_path, sep="\t")

        if not {"Consequence_term", "Uploaded_variation", "SYMBOL", "CHR", "BP"}.issubset(ann.columns):
            continue
        if not {"SNP", "PIP"}.issubset(pip.columns):
            continue

        ann = ann.copy()
        ann["is_target_coding"] = ann["Consequence_term"].apply(is_target_consequence)
        ann = ann[ann["is_target_coding"]].copy()
        if ann.empty:
            continue

        ann["CHR"] = ann["CHR"].astype(str).str.replace("chr", "", regex=False)
        ann["BP"] = pd.to_numeric(ann["BP"], errors="coerce")

        pip = pip.copy()
        pip["SNP"] = pip["SNP"].astype(str)
        pip["PIP"] = pd.to_numeric(pip["PIP"], errors="coerce")

        # Parse CHR/BP from SNP when in CHR:BP:REF:ALT format.
        pip["CHR"] = pip["SNP"].str.split(":").str[0].str.replace("chr", "", regex=False)
        pip["BP"] = pd.to_numeric(pip["SNP"].str.split(":").str[1], errors="coerce")

        # 1) direct ID join
        m1 = ann.merge(
            pip[["SNP", "PIP", "CHR", "BP"]],
            left_on="Uploaded_variation",
            right_on="SNP",
            how="left",
            suffixes=("", "_pip"),
        )

        # 2) coordinate fallback for missing PIP
        need = m1["PIP"].isna()
        if need.any():
            m2 = m1.loc[need, ["Uploaded_variation", "SYMBOL", "Consequence_term", "CHR", "BP"]].merge(
                pip[["SNP", "PIP", "CHR", "BP"]],
                on=["CHR", "BP"],
                how="left",
                suffixes=("", "_pip"),
            )
            m1.loc[need, "SNP"] = m2["SNP"].values
            m1.loc[need, "PIP"] = m2["PIP"].values

        m1["PIP"] = pd.to_numeric(m1["PIP"], errors="coerce")
        m1 = m1.dropna(subset=["PIP"])
        if m1.empty:
            continue

        m1 = m1.rename(
            columns={
                "SYMBOL": "gene",
                "Consequence_term": "consequence",
                "SNP": "variant_id",
            }
        )
        m1["locus_id"] = locus_id

        # keep strongest row per (locus, variant, gene, consequence)
        m1 = (
            m1.sort_values("PIP", ascending=False)
            .drop_duplicates(subset=["locus_id", "variant_id", "gene", "consequence"], keep="first")
        )

        rows.append(m1[["locus_id", "variant_id", "gene", "consequence", "PIP"]])

    if not rows:
        out = root / "results/reports/alphafold_high_confidence_candidates.tsv"
        pd.DataFrame(columns=["locus_id", "variant_id", "gene", "consequence", "PIP", "rank"]).to_csv(out, sep="\t", index=False)
        print("No coding variants matched to PIP across available loci.")
        print(f"wrote={out}")
        return

    df = pd.concat(rows, ignore_index=True)
    df = df[df["PIP"] >= 0.5].copy()

    if df.empty:
        out = root / "results/reports/alphafold_high_confidence_candidates.tsv"
        pd.DataFrame(columns=["locus_id", "variant_id", "gene", "consequence", "PIP", "rank"]).to_csv(out, sep="\t", index=False)
        print("No strong coding locus exists under current strict filter (PIP >= 0.5).")
        print(f"wrote={out}")
        return

    locus_stats = (
        df.groupby("locus_id", as_index=False)
        .agg(
            max_coding_pip=("PIP", "max"),
            n_coding_pip_ge_05=("PIP", "count"),
            n_coding_pip_ge_08=("PIP", lambda s: int((s >= 0.8).sum())),
        )
        .sort_values(["max_coding_pip", "n_coding_pip_ge_05"], ascending=[False, False])
        .reset_index(drop=True)
    )
    locus_stats["locus_rank"] = locus_stats.index + 1

    out_df = df.merge(locus_stats, on="locus_id", how="left")
    out_df = out_df.sort_values(["locus_rank", "PIP"], ascending=[True, False]).reset_index(drop=True)
    out_df["rank"] = out_df.index + 1

    out = root / "results/reports/alphafold_high_confidence_candidates.tsv"
    out_df[[
        "locus_id",
        "variant_id",
        "gene",
        "consequence",
        "PIP",
        "locus_rank",
        "max_coding_pip",
        "n_coding_pip_ge_05",
        "n_coding_pip_ge_08",
        "rank",
    ]].to_csv(out, sep="\t", index=False)

    top3 = locus_stats.head(3)
    print("top_3_loci_suitable_for_alphafold:")
    print(top3.to_string(index=False))
    print(f"wrote={out}")


if __name__ == "__main__":
    main()
