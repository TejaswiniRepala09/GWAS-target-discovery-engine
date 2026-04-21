"""Fine-mapping preparation pipeline (bundle + contracts + parse)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import polars as pl

from src.fine_mapping.credible_sets import normalize_credible_sets
from src.fine_mapping.ld import (
    validate_ld_resources_for_chromosomes,
    validate_locus_ld_matrix_contract,
)
from src.fine_mapping.locus_prep import create_locus_bundles
from src.fine_mapping.susie_io import generate_susie_commands, parse_susie_outputs

LOGGER = logging.getLogger(__name__)


def run_fine_mapping_preparation(
    cleaned_df: pl.DataFrame,
    lead_variants_df: pl.DataFrame,
    fine_cfg: dict[str, Any],
) -> dict[str, pl.DataFrame]:
    """Run locus bundle preparation + LD checks + susie parser/scaffold outputs."""
    root = Path(fine_cfg["paths"]["fine_mapping_root"])
    root.mkdir(parents=True, exist_ok=True)

    locus_meta = create_locus_bundles(cleaned_df, lead_variants_df, fine_cfg)
    locus_meta.write_csv(root / "locus_metadata.tsv", separator="\t") if locus_meta.height else None

    chroms = sorted(set(int(x) for x in lead_variants_df["CHR"].to_list())) if lead_variants_df.height and "CHR" in lead_variants_df.columns else []
    ld_chr_status = validate_ld_resources_for_chromosomes(chroms, fine_cfg)
    ld_chr_status.write_csv(root / "ld_chromosome_resource_status.tsv", separator="\t") if ld_chr_status.height else None

    ld_locus_status = validate_locus_ld_matrix_contract(locus_meta, fine_cfg) if locus_meta.height else pl.DataFrame()
    ld_locus_status.write_csv(root / "ld_locus_contract_status.tsv", separator="\t") if ld_locus_status.height else None

    susie_commands = generate_susie_commands(locus_meta, fine_cfg) if locus_meta.height else pl.DataFrame()
    susie_commands.write_csv(root / "susie_commands.tsv", separator="\t") if susie_commands.height else None

    susie = parse_susie_outputs(fine_cfg)
    pip = susie["pip_summary"]
    cs = normalize_credible_sets(susie["credible_sets"], pip)

    pip.write_csv(root / "pip_summary.tsv", separator="\t") if pip.height else None
    cs.write_csv(root / "credible_sets.tsv", separator="\t") if cs.height else None
    susie["susie_missing"].write_csv(root / "susie_missing_inputs.tsv", separator="\t") if susie["susie_missing"].height else None

    return {
        "locus_metadata": locus_meta,
        "ld_chr_status": ld_chr_status,
        "ld_locus_status": ld_locus_status,
        "susie_commands": susie_commands,
        "pip_summary": pip,
        "credible_sets": cs,
        "susie_missing": susie["susie_missing"],
    }
