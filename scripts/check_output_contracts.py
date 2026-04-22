#!/usr/bin/env python3
"""Lightweight contract checks for key committed TSV outputs.

Checks:
1) file exists
2) file is readable as TSV
3) required columns are present
4) key columns are not entirely null/empty
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

CONTRACTS = {
    "results/fine_mapping/pip_summary.tsv": {
        "required_columns": {"SNP", "PIP", "locus_id"},
        "non_empty_columns": {"SNP", "PIP", "locus_id"},
    },
    "results/fine_mapping/credible_sets.tsv": {
        "required_columns": {"credible_set_id", "SNP", "locus_id", "posterior_rank"},
        "non_empty_columns": {"credible_set_id", "SNP", "locus_id"},
    },
    "results/target_prioritization/variant_priority_scores.tsv": {
        "required_columns": {"SNP", "locus_id", "PIP", "variant_priority_score"},
        "non_empty_columns": {"SNP", "locus_id", "variant_priority_score"},
    },
    "results/target_prioritization/gene_prioritization.tsv": {
        "required_columns": {"gene_symbol"},
        "required_any_columns": [{"gene_prioritization_score", "gene_prioritization_score_final"}],
        "non_empty_columns": {"gene_symbol", "strongest_pvalue"},
    },
    "results/reports/multi_locus_status.tsv": {
        "required_columns": {
            "locus_id",
            "chromosome",
            "step_reached",
            "success",
            "error_message",
            "n_variants_input",
            "n_variants_after_harmonization",
            "susie_converged",
            "diagnostics_flag",
        },
        "non_empty_columns": {"locus_id", "chromosome", "success", "step_reached"},
    },
    "results/benchmarking/susie_vs_finemap_locus_summary.tsv": {
        "required_columns": {"locus_id", "susie_top_variant", "finemap_top_variant", "benchmark_status"},
        "non_empty_columns": {"locus_id", "benchmark_status"},
    },
}


def read_header(path: Path) -> list[str]:
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty TSV (no header): {path}") from exc


def is_populated(value: str) -> bool:
    text = str(value).strip()
    return text not in {"", "NA", "N/A", "NaN", "nan", "null", "NULL", "None"}


def check_non_empty_columns(path: Path, header: list[str], cols: set[str]) -> list[str]:
    idx_map = {name: i for i, name in enumerate(header)}
    has_value = {c: False for c in cols}

    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # skip header
        for row in reader:
            for col in cols:
                i = idx_map[col]
                if i < len(row) and is_populated(row[i]):
                    has_value[col] = True
            if all(has_value.values()):
                break

    return [c for c, ok in has_value.items() if not ok]


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

        any_group_failed = False
        for group in spec.get("required_any_columns", []):
            if not (set(group) & set(header)):
                failures.append(f"[schema] {rel_path}: expected one of {sorted(group)}")
                any_group_failed = True
        if any_group_failed:
            continue

        missing_values = check_non_empty_columns(path, header, set(spec.get("non_empty_columns", set())))
        if missing_values:
            failures.append(f"[empty] {rel_path}: columns entirely null/empty {missing_values}")
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
