"""UniProt mapping for prioritized genes with local cache."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import polars as pl

try:
    import httpx
except Exception:  # noqa: BLE001
    httpx = None


def _cache_path(cache_dir: Path, key: str) -> Path:
    return cache_dir / f"{hashlib.sha256(key.encode('utf-8')).hexdigest()}.json"


def _query_uniprot(base_url: str, timeout: int, gene_symbol: str) -> dict[str, Any]:
    if httpx is None:
        raise RuntimeError("httpx is not installed. Install dependencies from requirements.txt")
    # Search endpoint for reviewed human proteins by gene name
    url = f"{base_url.rstrip('/')}/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606",
        "fields": "accession,gene_names,protein_name",
        "format": "json",
        "size": 1,
    }
    with httpx.Client(timeout=timeout) as client:
        resp = client.get(url, params=params)
        resp.raise_for_status()
        return resp.json()


def build_protein_mapping(gene_scores: pl.DataFrame, structure_cfg: dict[str, Any]) -> pl.DataFrame:
    """Map gene symbols to UniProt accessions and write cache-backed table."""
    if gene_scores.height == 0:
        return pl.DataFrame()

    base_url = structure_cfg["apis"]["uniprot_base_url"]
    timeout = int(structure_cfg["apis"].get("timeout_seconds", 20))
    cache_dir = Path(structure_cfg["paths"]["uniprot_cache_dir"])
    cache_dir.mkdir(parents=True, exist_ok=True)

    max_genes = int(structure_cfg.get("paths", {}).get("max_genes_query", 50))
    rows: list[dict[str, Any]] = []
    for row in gene_scores.select([c for c in ["gene_symbol", "ensembl_gene_id"] if c in gene_scores.columns]).unique().head(max_genes).iter_rows(named=True):
        symbol = row.get("gene_symbol")
        if not symbol:
            continue
        cpath = _cache_path(cache_dir, f"uniprot::{symbol}")
        if cpath.exists():
            payload = json.loads(cpath.read_text(encoding="utf-8"))
        else:
            try:
                payload = _query_uniprot(base_url, timeout, str(symbol))
            except Exception as exc:  # noqa: BLE001
                payload = {"status": "error", "message": str(exc)}
            cpath.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        if payload.get("status") == "error":
            rows.append(
                {
                    "gene_symbol": symbol,
                    "ensembl_gene_id": row.get("ensembl_gene_id"),
                    "uniprot_accession": None,
                    "uniprot_status": "error",
                    "uniprot_note": payload.get("message"),
                }
            )
            continue

        result = (payload.get("results") or [None])[0]
        accession = result.get("primaryAccession") if result else None
        rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_gene_id": row.get("ensembl_gene_id"),
                "uniprot_accession": accession,
                "uniprot_status": "ok" if accession else "not_found",
                "uniprot_note": None if accession else "No UniProt match in query result",
            }
        )

    return pl.from_dicts(rows) if rows else pl.DataFrame()
