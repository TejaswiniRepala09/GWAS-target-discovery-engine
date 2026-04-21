#!/usr/bin/env python3
"""Generate representative-locus interpretation and AlphaFold candidate reports."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import pandas as pd
import numpy as np

REP_LOCI = ["chr7_locus19", "chr7_locus14", "chr7_locus10"]
CODING_TERMS = {
    "missense_variant",
    "protein_altering_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
}


@dataclass
class LocusSummary:
    locus_id: str
    max_pip: float
    n_nonzero: int
    n_pip_ge_01: int
    cs_count: int
    cs_sizes: list[int]
    cs_max: int
    cs_min: int
    top_variants: list[dict]
    top_genes: list[dict]
    ld_quality: str
    ld_well_conditioned: str
    rss_consistency: str
    converged: str
    coding_variant_count: int
    coding_top_pip: float | None
    coding_examples: list[dict]
    has_protein_change: bool


def read_diag_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if ": " in line and not line.startswith("["):
            k, v = line.split(": ", 1)
            out[k.strip()] = v.strip()
    return out


def infer_ld_complexity(cs_max: int, n_pip_ge_01: int, n_nonzero: int, diag: dict[str, str]) -> str:
    ld_flag = diag.get("ld_well_conditioned_flag", "UNKNOWN")
    if cs_max <= 2 and n_pip_ge_01 <= 10 and n_nonzero <= 20 and ld_flag == "TRUE":
        return "simple"
    if cs_max >= 5 or n_pip_ge_01 >= 20 or n_nonzero >= 60:
        return "complex"
    return "moderate"


def load_locus_summary(root: Path, locus_id: str) -> LocusSummary:
    pip_path = root / f"results/fine_mapping/susie_raw/{locus_id}_regularized/pip.tsv"
    cs_path = root / f"results/fine_mapping/susie_raw/{locus_id}_regularized/credible_sets.tsv"
    if not pip_path.exists():
        pip_path = root / f"results/fine_mapping/susie_raw/{locus_id}/pip.tsv"
    if not cs_path.exists():
        cs_path = root / f"results/fine_mapping/susie_raw/{locus_id}/credible_sets.tsv"
    vp_path = root / f"results/target_prioritization/{locus_id}/variant_priority_scores.tsv"
    gp_path = root / f"results/target_prioritization/{locus_id}/gene_prioritization.tsv"
    ann_path = root / f"results/annotations/{locus_id}/variant_annotations.tsv"
    diag_path = root / f"results/fine_mapping/diagnostics/{locus_id}_regularized_diagnostics.txt"
    if not diag_path.exists():
        diag_path = root / f"results/fine_mapping/diagnostics/{locus_id}_diagnostics.txt"

    for p in [pip_path, cs_path, vp_path, gp_path, ann_path]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required per-locus file: {p}")

    pip = pd.read_csv(pip_path, sep="\t")
    cs = pd.read_csv(cs_path, sep="\t")
    vp = pd.read_csv(vp_path, sep="\t")
    gp = pd.read_csv(gp_path, sep="\t")
    ann = pd.read_csv(ann_path, sep="\t")
    diag = read_diag_kv(diag_path)

    pip["PIP"] = pd.to_numeric(pip["PIP"], errors="coerce").fillna(0.0)
    max_pip = float(pip["PIP"].max())
    n_nonzero = int((pip["PIP"] > 0).sum())
    n_pip_ge_01 = int((pip["PIP"] >= 0.1).sum())
    top_variants = pip.sort_values("PIP", ascending=False)[["SNP", "PIP"]].head(5).to_dict("records")

    cs_sizes = cs.groupby("credible_set_id").size().sort_values().tolist() if len(cs) else []
    cs_count = len(cs_sizes)
    cs_max = int(max(cs_sizes)) if cs_sizes else 0
    cs_min = int(min(cs_sizes)) if cs_sizes else 0

    gscore_col = "gene_prioritization_score_final" if "gene_prioritization_score_final" in gp.columns else "gene_prioritization_score"
    gp[gscore_col] = pd.to_numeric(gp[gscore_col], errors="coerce")
    top_genes = (
        gp.sort_values(gscore_col, ascending=False)[["gene_symbol", gscore_col]]
        .dropna(subset=["gene_symbol"])
        .query("gene_symbol != '-'")
        .head(5)
        .rename(columns={gscore_col: "score"})
        .to_dict("records")
    )

    # coding evidence:
    # - use prioritized variant table for coding-PIP alignment (same locus/SNP namespace)
    # - use annotation table for protein-change descriptors (HGVSp/Protein_position)
    vp["Consequence_term"] = vp["Consequence_term"].astype(str)
    coding_vp = vp[vp["Consequence_term"].isin(CODING_TERMS)].copy()
    ann["Consequence_term"] = ann["Consequence_term"].astype(str)
    coding_ann = ann[ann["Consequence_term"].isin(CODING_TERMS)].copy()

    coding_top_pip = None
    coding_examples: list[dict] = []
    has_protein_change = False
    coding_count = 0
    if len(coding_vp):
        coding_vp["PIP"] = pd.to_numeric(coding_vp["PIP"], errors="coerce")
        coding_count = int(coding_vp[["SNP"]].drop_duplicates().shape[0])
        coding_top_pip = float(coding_vp["PIP"].max()) if coding_vp["PIP"].notna().any() else None
    if len(coding_ann):
        m = coding_ann.copy()
        has_protein_change = bool(
            m["HGVSp"].astype(str).str.strip().replace({"nan": ""}).ne("").any()
            or m["Protein_position"].astype(str).str.strip().replace({"nan": ""}).ne("").any()
        )
        # enrich with PIP via CHR/BP lookup from prioritized table when available
        if len(coding_vp) and {"CHR", "BP"}.issubset(coding_vp.columns):
            map_pip = (
                coding_vp.assign(
                    CHR=coding_vp["CHR"].astype(str).str.replace("chr", "", regex=False),
                    BP=pd.to_numeric(coding_vp["BP"], errors="coerce"),
                )[["CHR", "BP", "PIP"]]
                .dropna(subset=["CHR", "BP"])
            )
            m = m.assign(
                CHR=m["CHR"].astype(str).str.replace("chr", "", regex=False),
                BP=pd.to_numeric(m["BP"], errors="coerce"),
            ).merge(map_pip, on=["CHR", "BP"], how="left")
        else:
            m["PIP"] = np.nan
        coding_examples = (
            m.sort_values("PIP", ascending=False)[
                ["Uploaded_variation", "Consequence_term", "SYMBOL", "PIP", "HGVSp", "Protein_position"]
            ]
            .drop_duplicates(subset=["Uploaded_variation", "Consequence_term", "SYMBOL"])
            .head(5)
            .to_dict("records")
        )

    ld_quality = infer_ld_complexity(cs_max, n_pip_ge_01, n_nonzero, diag)

    return LocusSummary(
        locus_id=locus_id,
        max_pip=max_pip,
        n_nonzero=n_nonzero,
        n_pip_ge_01=n_pip_ge_01,
        cs_count=cs_count,
        cs_sizes=cs_sizes,
        cs_max=cs_max,
        cs_min=cs_min,
        top_variants=top_variants,
        top_genes=top_genes,
        ld_quality=ld_quality,
        ld_well_conditioned=diag.get("ld_well_conditioned_flag", "UNKNOWN"),
        rss_consistency=diag.get("rss_consistency_flag", "UNKNOWN"),
        converged=diag.get("converged", "UNKNOWN"),
        coding_variant_count=coding_count,
        coding_top_pip=coding_top_pip,
        coding_examples=coding_examples,
        has_protein_change=has_protein_change,
    )


def confidence_label(s: LocusSummary) -> str:
    if s.max_pip >= 0.99 and s.cs_max <= 2 and s.ld_quality == "simple":
        return "high confidence"
    if s.ld_quality == "complex" or s.n_nonzero >= 60 or s.cs_max >= 5:
        return "caution / exploratory"
    return "moderate / complex"


def write_locus_report(root: Path, s: LocusSummary) -> Path:
    out = root / f"results/reports/{s.locus_id}_interpretation.md"
    lines = []
    lines.append(f"# {s.locus_id} Interpretation")
    lines.append("")
    lines.append("## 1) Statistical Pattern")
    lines.append(f"- Max PIP: `{s.max_pip:.6g}`")
    lines.append(f"- Credible sets: `{s.cs_count}`; sizes: `{s.cs_sizes}`")
    lines.append(f"- Variants with nonzero PIP: `{s.n_nonzero}`")
    lines.append(f"- Variants with PIP >= 0.1: `{s.n_pip_ge_01}`")
    focus = "focused" if (s.cs_max <= 2 and s.n_pip_ge_01 <= 10) else "diffuse/complex"
    lines.append(f"- Signal pattern: `{focus}`")
    lines.append("")
    lines.append("## 2) LD Interpretation")
    lines.append(f"- LD complexity: `{s.ld_quality}`")
    lines.append(f"- Diagnostic flags: `ld_well_conditioned={s.ld_well_conditioned}`, `rss_consistency={s.rss_consistency}`, `susie_converged={s.converged}`")
    lines.append("")
    lines.append("## 3) Biological Interpretation")
    tg = ", ".join([f"{x['gene_symbol']} ({x['score']:.3g})" for x in s.top_genes]) if s.top_genes else "No prioritized gene available"
    lines.append(f"- Top candidate genes: {tg}")
    tv = ", ".join([f"{x['SNP']} (PIP={x['PIP']:.3g})" for x in s.top_variants]) if s.top_variants else "No top variants available"
    lines.append(f"- Top candidate variants: {tv}")
    lines.append(f"- Coding variant count (selected consequence classes): `{s.coding_variant_count}`")
    if s.coding_top_pip is not None:
        lines.append(f"- Max PIP among coding variants: `{s.coding_top_pip:.6g}`")
    else:
        lines.append("- Max PIP among coding variants: `NA`")
    model_app = "applicable" if (s.coding_variant_count > 0 and s.has_protein_change) else "not currently applicable"
    lines.append(f"- Structural modeling applicability: `{model_app}`")
    if s.coding_examples:
        lines.append("- Coding examples:")
        for x in s.coding_examples:
            lines.append(
                f"  - {x.get('Uploaded_variation','')}: {x.get('Consequence_term','')} | gene={x.get('SYMBOL','')} | PIP={x.get('PIP','')} | HGVSp={x.get('HGVSp','')}"
            )
    lines.append("")
    lines.append("## 4) Confidence Statement")
    lines.append(f"- Confidence: `{confidence_label(s)}`")
    lines.append("- Interpretation caveat: fine-mapping + annotation prioritize candidates but do not prove mechanism/causality.")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> None:
    root = Path.cwd()

    # A) coding-locus search across available loci
    loci = []
    for p in (root / "results/annotations").glob("chr*_locus*"):
        if p.is_dir():
            loci.append(p.name)
    loci = sorted(set(loci))

    records = []
    for locus in loci:
        try:
            s = load_locus_summary(root, locus)
        except FileNotFoundError:
            continue
        # pull coding-support stats via per-locus variant_priority_scores
        vp_path = root / f"results/target_prioritization/{locus}/variant_priority_scores.tsv"
        vp = pd.read_csv(vp_path, sep="\t")
        vp["Consequence_term"] = vp["Consequence_term"].astype(str)
        vpc = vp[vp["Consequence_term"].isin(CODING_TERMS)].copy()
        vpc["variant_priority_score"] = pd.to_numeric(vpc["variant_priority_score"], errors="coerce")
        top_vscore = float(vpc["variant_priority_score"].max()) if len(vpc) else np.nan

        records.append(
            {
                "locus_id": locus,
                "success": True,
                "max_pip_all_variants": s.max_pip,
                "max_pip_coding_variants": float(pd.to_numeric(vpc["PIP"], errors="coerce").max()) if len(vpc) else np.nan,
                "coding_variant_count": s.coding_variant_count,
                "top_coding_variant_priority_score": top_vscore,
                "susie_converged": s.converged,
                "ld_well_conditioned": s.ld_well_conditioned,
                "rss_consistency_flag": s.rss_consistency,
                "has_protein_change_annotation": s.has_protein_change,
            }
        )

    cand = pd.DataFrame(records)
    if cand.empty:
        raise RuntimeError("No loci found with complete per-locus outputs for coding-locus scan.")

    cand["max_pip_coding_variants"] = pd.to_numeric(cand["max_pip_coding_variants"], errors="coerce")
    cand["top_coding_variant_priority_score"] = pd.to_numeric(cand["top_coding_variant_priority_score"], errors="coerce")
    cand["coding_variant_count"] = pd.to_numeric(cand["coding_variant_count"], errors="coerce").fillna(0).astype(int)

    cand = cand.sort_values(
        [
            "max_pip_coding_variants",
            "top_coding_variant_priority_score",
            "coding_variant_count",
        ],
        ascending=[False, False, False],
        na_position="last",
    )
    out_candidates = root / "results/reports/alphafold_candidate_loci.tsv"
    cand.to_csv(out_candidates, sep="\t", index=False)

    # B) representative loci deep interpretation
    reps = []
    rep_paths = []
    for locus in REP_LOCI:
        s = load_locus_summary(root, locus)
        reps.append(s)
        rep_paths.append(str(write_locus_report(root, s)))

    # C) combined summary
    combined = root / "results/reports/representative_loci_summary.md"
    lines = []
    lines.append("# Representative Loci Summary")
    lines.append("")
    lines.append("## Why These 3 Loci")
    lines.append("- `chr7_locus19`: clean/strong case with compact posterior support.")
    lines.append("- `chr7_locus14`: moderate complexity with broader credible-set structure.")
    lines.append("- `chr7_locus10`: diffuse/complex case that highlights limitations and caution.")
    lines.append("")
    lines.append("## What Each Locus Teaches")
    for s in reps:
        lines.append(f"### {s.locus_id}")
        lines.append(f"- Max PIP: `{s.max_pip:.6g}`")
        lines.append(f"- Credible sets: `{s.cs_count}` with sizes `{s.cs_sizes}`")
        lines.append(f"- LD interpretation: `{s.ld_quality}` (`ld_well_conditioned={s.ld_well_conditioned}`, `rss_consistency={s.rss_consistency}`)")
        lines.append(f"- Top genes: {', '.join([x['gene_symbol'] for x in s.top_genes]) if s.top_genes else 'NA'}")
        lines.append(f"- Coding signal: `{s.coding_variant_count}` coding variants; structure applicability: `{'yes' if (s.coding_variant_count>0 and s.has_protein_change) else 'no'}`")
        lines.append(f"- Confidence statement: `{confidence_label(s)}`")
        lines.append("")
    lines.append("## Why Together They Tell a Strong Story")
    lines.append("- They show a full translational range: clean locus, realistic complexity, and a cautionary diffuse locus.")
    lines.append("- This supports honest target discovery communication: prioritization strength, uncertainty, and limitations.")
    lines.append("")
    lines.append("## Where AlphaFold Fits")
    lines.append("- AlphaFold-style structural follow-up is most appropriate when coding variants have protein-change annotations (for example HGVSp/protein_position).")
    lines.append("- It is supportive evidence only, not proof of causal mechanism.")
    lines.append("- Loci dominated by non-coding/regulatory signatures should prioritize functional genomics before structural claims.")
    combined.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "alphafold_candidate_table": str(out_candidates),
        "best_alphafold_locus": cand.iloc[0]["locus_id"],
        "representative_reports": rep_paths,
        "combined_summary": str(combined),
    }, indent=2))


if __name__ == "__main__":
    main()
