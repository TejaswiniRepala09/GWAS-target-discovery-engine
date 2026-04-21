"""Phase 2 markdown reporting."""

from __future__ import annotations

from pathlib import Path

import polars as pl


def _preview(df: pl.DataFrame, n: int = 10) -> str:
    if df.height == 0:
        return "No rows available."
    return str(df.head(n))


def write_phase2_report(
    output_path: str | Path,
    variant_scores: pl.DataFrame,
    gene_scores: pl.DataFrame,
    credible_sets: pl.DataFrame,
    structure_candidates: pl.DataFrame,
) -> Path:
    """Write concise phase 2 summary report."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(
        [
            "# Phase 2 Summary Report",
            "",
            "This report summarizes Phase 2 prioritization outputs.",
            "",
            "## Caveats",
            "- GWAS association is not causality.",
            "- Fine-mapping supports statistical prioritization, not mechanism proof.",
            "- Nearest-gene fallback is heuristic.",
            "- Structure candidate checks are follow-up pointers, not biophysical estimates.",
            "",
            "## Top Variants",
            _preview(variant_scores),
            "",
            "## Top Genes",
            _preview(gene_scores),
            "",
            "## Credible Sets",
            _preview(credible_sets),
            "",
            "## Structure Candidates",
            _preview(structure_candidates),
        ]
    )
    path.write_text(text, encoding="utf-8")
    return path
