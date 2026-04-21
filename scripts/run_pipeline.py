#!/usr/bin/env python3
"""Run full CKD target discovery pipeline."""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.pipeline import run_pipeline


if __name__ == "__main__":
    run_pipeline()
