"""Logging helpers for Phase 2."""

from __future__ import annotations

import logging
from pathlib import Path


def setup_phase2_logging(log_file: str | Path, level: int = logging.INFO) -> None:
    """Configure phase 2 logger with stream + file handlers."""
    path = Path(log_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(path)],
        force=True,
    )
