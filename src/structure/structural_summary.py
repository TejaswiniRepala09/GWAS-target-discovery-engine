"""Build structure candidates from variant/gene tables and API checks."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from src.structure.alphafold_scaffold import annotate_alphafold_followup
from src.structure.mutation_formatting import short_protein_change
from src.structure.structure_fetch import query_alphafold, query_pdb


def build_structure_candidates(
    variant_scores: pl.DataFrame,
    protein_mapping: pl.DataFrame,
    structure_cfg: dict,
) -> dict[str, pl.DataFrame]:
    """Build structure candidate and summary tables with availability checks only."""
    if variant_scores.height == 0:
        empty = pl.DataFrame()
        return {"structure_candidates": empty, "structural_summary": empty}

    keywords = structure_cfg.get("consequence_keywords", [])
    cons_col = "Consequence_term" if "Consequence_term" in variant_scores.columns else None
    candidates = variant_scores
    if cons_col:
        pattern = "|".join(keywords)
        candidates = candidates.filter(pl.col(cons_col).cast(pl.Utf8).str.contains(pattern, literal=False)) if pattern else candidates

    if protein_mapping.height and "gene_symbol" in candidates.columns and "gene_symbol" in protein_mapping.columns:
        candidates = candidates.join(protein_mapping, on=[c for c in ["gene_symbol", "ensembl_gene_id"] if c in candidates.columns and c in protein_mapping.columns], how="left")

    candidates = candidates.with_columns(
        pl.col("HGVSp").map_elements(short_protein_change, return_dtype=pl.Utf8).alias("protein_change_short") if "HGVSp" in candidates.columns else pl.lit(None).alias("protein_change_short")
    )

    # Lightweight availability checks: if accession-like field exists, query AlphaFold endpoint.
    timeout = int(structure_cfg["apis"].get("timeout_seconds", 20))
    pdb_cache = Path(structure_cfg["paths"]["pdb_cache_dir"])
    af_cache = Path(structure_cfg["paths"]["alphafold_cache_dir"])
    pdb_cache.mkdir(parents=True, exist_ok=True)
    af_cache.mkdir(parents=True, exist_ok=True)

    if "uniprot_accession" in candidates.columns:
        af_exists: list[bool] = []
        for acc in candidates["uniprot_accession"].to_list():
            if acc is None:
                af_exists.append(False)
                continue
            payload = query_alphafold(structure_cfg["apis"]["alphafold_base_url"], str(acc), timeout, af_cache)
            af_exists.append(bool(payload.get("exists", False)))
        candidates = candidates.with_columns(pl.Series("alphafold_model_available", af_exists))

    # PDB id may be unavailable in this stage; keep placeholder logic.
    if "pdb_id" in candidates.columns:
        pdb_exists: list[bool] = []
        for pid in candidates["pdb_id"].to_list():
            if pid is None:
                pdb_exists.append(False)
                continue
            payload = query_pdb(structure_cfg["apis"]["pdb_base_url"], str(pid), timeout, pdb_cache)
            pdb_exists.append(bool(payload.get("exists", False)))
        candidates = candidates.with_columns(pl.Series("pdb_entry_available", pdb_exists))

    candidates = annotate_alphafold_followup(candidates)

    summary_cols = [
        "gene_symbol",
        "ensembl_gene_id",
        "SNP",
        "variant_priority_score",
        "Consequence_term",
        "HGVSp",
        "protein_change_short",
        "uniprot_accession",
        "alphafold_model_available",
        "pdb_entry_available",
        "structural_followup_note",
    ]
    summary = candidates.select([c for c in summary_cols if c in candidates.columns]).sort("variant_priority_score", descending=True) if candidates.height else pl.DataFrame()
    return {"structure_candidates": candidates, "structural_summary": summary}
