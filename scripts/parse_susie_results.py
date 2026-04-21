#!/usr/bin/env python3
"""Parse real SuSiE outputs from results/fine_mapping/susie_raw/."""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.fine_mapping.susie_io import parse_susie_outputs
from src.fine_mapping.credible_sets import normalize_credible_sets
from src.io_utils import load_yaml
from src.utils.logging_utils import setup_phase2_logging


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    fine_cfg = load_yaml("config/fine_mapping.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    parsed = parse_susie_outputs(fine_cfg)
    pip = parsed["pip_summary"]
    cs = normalize_credible_sets(parsed["credible_sets"], pip)

    out_root = Path(fine_cfg["paths"]["fine_mapping_root"])
    out_root.mkdir(parents=True, exist_ok=True)
    if pip.height:
        pip.write_csv(out_root / "pip_summary.tsv", separator="\t")
    if cs.height:
        cs.write_csv(out_root / "credible_sets.tsv", separator="\t")
    if parsed["susie_missing"].height:
        parsed["susie_missing"].write_csv(out_root / "susie_missing_inputs.tsv", separator="\t")

    print(f"Parsed PIP rows: {pip.height}")
    print(f"Parsed credible set rows: {cs.height}")
