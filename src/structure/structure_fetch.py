"""PDB and AlphaFold availability checks with JSON cache."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

try:
    import httpx
except Exception:  # noqa: BLE001
    httpx = None


def _cache(cache_dir: Path, key: str) -> Path:
    return cache_dir / f"{hashlib.sha256(key.encode('utf-8')).hexdigest()}.json"


def query_pdb(pdb_base_url: str, pdb_id: str, timeout: int, cache_dir: Path) -> dict[str, Any]:
    key = f"pdb::{pdb_id}"
    cpath = _cache(cache_dir, key)
    if cpath.exists():
        return json.loads(cpath.read_text(encoding="utf-8"))
    url = f"{pdb_base_url.rstrip('/')}/{pdb_id}"
    try:
        if httpx is None:
            raise RuntimeError("httpx is not installed. Install dependencies from requirements.txt")
        with httpx.Client(timeout=timeout) as client:
            resp = client.get(url)
        payload = {"status_code": resp.status_code, "exists": resp.status_code == 200}
    except Exception as exc:  # noqa: BLE001
        payload = {"status": "error", "message": str(exc), "exists": False}
    cpath.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def query_alphafold(alphafold_base_url: str, uniprot_acc: str, timeout: int, cache_dir: Path) -> dict[str, Any]:
    key = f"alphafold::{uniprot_acc}"
    cpath = _cache(cache_dir, key)
    if cpath.exists():
        return json.loads(cpath.read_text(encoding="utf-8"))
    url = f"{alphafold_base_url.rstrip('/')}/prediction/{uniprot_acc}"
    try:
        if httpx is None:
            raise RuntimeError("httpx is not installed. Install dependencies from requirements.txt")
        with httpx.Client(timeout=timeout) as client:
            resp = client.get(url)
        payload = {"status_code": resp.status_code, "exists": resp.status_code == 200}
    except Exception as exc:  # noqa: BLE001
        payload = {"status": "error", "message": str(exc), "exists": False}
    cpath.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload
