#!/usr/bin/env python3
"""Validate VEP runtime contract, write command, and parse VEP output if present."""

from pathlib import Path
import sys
import argparse

sys.path.append(str(Path(__file__).resolve().parent.parent))

import polars as pl

from src.annotation.consequence_scoring import score_consequences
from src.annotation.vep_parser import build_vep_input_file, generate_local_vep_command, parse_vep_output
from src.io_utils import load_yaml
from src.utils.logging_utils import setup_phase2_logging
from src.utils.validation import require_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse VEP output and score consequences.")
    parser.add_argument(
        "--vep-raw-output",
        default=None,
        help="Optional explicit VEP raw output path (e.g., per-locus file).",
    )
    args = parser.parse_args()

    settings = load_yaml("config/settings.yaml")
    vep_cfg = load_yaml("config/vep.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    vep_input_path = require_file(vep_cfg["paths"]["vep_input"], "vep input")
    vep_input = pl.read_csv(vep_input_path, separator="\t")
    formatted = build_vep_input_file(vep_input, vep_cfg)

    cmd, status = generate_local_vep_command(vep_cfg)
    cmd_path = Path(vep_cfg["paths"]["vep_raw_dir"]) / "vep_command.sh"
    cmd_path.parent.mkdir(parents=True, exist_ok=True)
    cmd_path.write_text(cmd + "\n", encoding="utf-8")

    print(f"VEP formatted input: {formatted}")
    print(f"VEP command script: {cmd_path}")
    print(f"VEP in PATH: {status['vep_in_path']} | cache present: {status['cache_present']} ({status['cache_dir']})")

    raw = Path(args.vep_raw_output) if args.vep_raw_output else Path(vep_cfg["paths"]["vep_raw_output"])
    if raw.exists():
        ann = parse_vep_output(vep_cfg, raw_output_path=raw)
        scored = score_consequences(ann, vep_cfg)
        scored.write_csv(vep_cfg["paths"]["variant_consequence_scores_output"], separator="\t")
        print(f"Parsed VEP rows: {ann.height}")
        print(f"Scored consequence rows: {scored.height}")
    else:
        print(f"VEP raw output not found yet: {raw}")
