"""Open Targets GraphQL integration with local JSON caching."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import polars as pl

LOGGER = logging.getLogger(__name__)

try:
    import httpx
except Exception:  # noqa: BLE001
    httpx = None


QUERY = """
query targetSummary($ensg: String!) {
  target(ensemblId: $ensg) {
    id
    approvedSymbol
    approvedName
    tractability {
      label
      modality
      value
    }
  }
}
"""


def _cache_file(cache_dir: Path, key: str) -> Path:
    return cache_dir / f"{key}.json"


def _hash_key(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _query_target(base_url: str, timeout: int, ensembl_gene_id: str) -> dict[str, Any]:
    if httpx is None:
        raise RuntimeError("httpx is not installed. Install dependencies from requirements.txt")
    with httpx.Client(timeout=timeout) as client:
        response = client.post(base_url, json={"query": QUERY, "variables": {"ensg": ensembl_gene_id}})
        response.raise_for_status()
        return response.json()


def fetch_opentargets_evidence(gene_df: pl.DataFrame, ot_cfg: dict[str, Any]) -> pl.DataFrame:
    """Fetch Open Targets records by Ensembl ID, fallback to ambiguity warning row."""
    if gene_df.height == 0:
        return pl.DataFrame()

    base_url = ot_cfg["opentargets"]["base_url"]
    timeout = int(ot_cfg["opentargets"].get("timeout_seconds", 25))
    enabled = bool(ot_cfg["opentargets"].get("enabled", True))
    cache_dir = Path(ot_cfg["paths"]["cache_dir"])
    cache_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    max_genes = int(ot_cfg.get("query", {}).get("page_size", 25))
    unique = gene_df.select([c for c in ["gene_symbol", "ensembl_gene_id"] if c in gene_df.columns]).unique().head(max_genes)

    for row in unique.iter_rows(named=True):
        gene_symbol = row.get("gene_symbol")
        ensg = row.get("ensembl_gene_id")

        if not ensg:
            rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": None,
                    "opentargets_status": "skipped_missing_ensembl",
                    "opentargets_support": 0.0,
                    "opentargets_note": "No Ensembl ID; symbol-only query is ambiguous and skipped by design.",
                }
            )
            continue

        key = _hash_key(str(ensg))
        cfile = _cache_file(cache_dir, key)

        payload: dict[str, Any]
        if cfile.exists():
            payload = json.loads(cfile.read_text(encoding="utf-8"))
        elif not enabled:
            payload = {"status": "disabled"}
            cfile.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        else:
            try:
                payload = _query_target(base_url, timeout, str(ensg))
            except Exception as exc:  # noqa: BLE001
                payload = {"status": "error", "message": str(exc)}
                LOGGER.warning("Open Targets query failed for %s: %s", ensg, exc)
            cfile.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        if payload.get("status") == "disabled":
            rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": ensg,
                    "opentargets_status": "disabled",
                    "opentargets_support": 0.0,
                    "opentargets_note": "Open Targets client disabled in config.",
                }
            )
            continue

        if payload.get("status") == "error":
            rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": ensg,
                    "opentargets_status": "error",
                    "opentargets_support": 0.0,
                    "opentargets_note": payload.get("message"),
                }
            )
            continue

        target = payload.get("data", {}).get("target") or {}
        tractability = target.get("tractability") or []
        tract_label = tractability[0].get("label") if tractability else None
        rows.append(
            {
                "gene_symbol": gene_symbol,
                "ensembl_gene_id": ensg,
                "opentargets_status": "ok",
                "opentargets_support": 1.0 if target.get("id") else 0.0,
                "opentargets_target_id": target.get("id"),
                "opentargets_approved_symbol": target.get("approvedSymbol"),
                "opentargets_tractability_label": tract_label,
                "opentargets_note": None,
            }
        )

    return pl.from_dicts(rows) if rows else pl.DataFrame()
