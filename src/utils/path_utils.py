"""Path and directory utilities."""

from __future__ import annotations

from pathlib import Path


def ensure_dir(path: str | Path) -> Path:
    """Create directory if missing."""
    resolved = Path(path)
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def first_existing_dir(candidates: list[str] | tuple[str, ...]) -> Path | None:
    """Return first existing directory from candidate list."""
    for item in candidates:
        p = Path(item)
        if p.exists() and p.is_dir():
            return p
    return None
