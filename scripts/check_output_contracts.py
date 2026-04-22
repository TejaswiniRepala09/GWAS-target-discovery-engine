#!/usr/bin/env python3
"""Lightweight contract checks for key committed TSV outputs.

Checks:
1) file exists
2) file is readable as TSV
3) required columns are present
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

CONTRACTS = {
    "results/fine_mapping/pip_summary.tsv": {
        "required_columns": {"SNP", "PIP", "locus_id"},
    },
    "results/fine_mapping/credible_sets.tsv": {
        "required_columns": {"credible_set_id", "SNP", "locus_id", "posterior_rank"},
    },
    "results/target_prioritization/variant_priority_scores.tsv": {
        "required_columns": {"SNP", "locus_id", "PIP", "variant_priority_score", "Uploaded_variation"},
    },
    "results/target_prioritization/gene_prioritization.tsv": {
        "required_columns": {"gene_symbol", "locus_id", "max_pip", "gene_prioritization_score"},
    },
}


def read_header(path: Path) -> list[str]:
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty TSV (no header): {path}") from exc


def main() -> int:
    failures: list[str] = []
    print("Checking output contracts...")
    for rel_path, spec in CONTRACTS.items():
        path = ROOT / rel_path
        if not path.exists():
            failures.append(f"[missing] {rel_path}")
            continue
        try:
            header = read_header(path)
        except Exception as exc:
            failures.append(f"[unreadable] {rel_path}: {exc}")
            continue

        missing_cols = sorted(spec["required_columns"] - set(header))
        if missing_cols:
            failures.append(f"[schema] {rel_path}: missing columns {missing_cols}")
            continue

        print(f"[ok] {rel_path} ({len(header)} columns)")

    if failures:
        print("\nContract check failures:")
        for item in failures:
            print(f"- {item}")
        return 1

    print("\nAll output contracts passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

