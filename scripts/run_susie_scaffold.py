#!/usr/bin/env python3
"""Generate susieR script template and per-locus command manifest."""

from pathlib import Path
import sys

import polars as pl

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.fine_mapping.susie_io import generate_susie_commands, write_susie_r_template
from src.io_utils import load_yaml
from src.utils.logging_utils import setup_phase2_logging


if __name__ == "__main__":
    settings = load_yaml("config/settings.yaml")
    fine_cfg = load_yaml("config/fine_mapping.yaml")
    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    meta_path = Path(fine_cfg["paths"]["fine_mapping_root"]) / "locus_metadata.tsv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing locus metadata. Run prepare_fine_mapping_inputs.py first: {meta_path}")

    manifest = pl.read_csv(meta_path, separator="\t")
    template = write_susie_r_template(fine_cfg["susie"]["rscript_template"])
    cmds = generate_susie_commands(manifest, fine_cfg)
    out_path = Path(fine_cfg["paths"]["fine_mapping_root"]) / "susie_commands.tsv"
    cmds.write_csv(out_path, separator="\t")

    print(f"Wrote template: {template}")
    print(f"Wrote commands: {out_path}")
