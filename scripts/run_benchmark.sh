#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

export PATH="$ROOT_DIR/tools/bin:$PATH"
export MPLBACKEND=Agg
export MPLCONFIGDIR="$ROOT_DIR/results/logs/.mpl"

echo ">>> .venv/bin/python -u scripts/run_finemap_benchmark.py"
.venv/bin/python -u scripts/run_finemap_benchmark.py
