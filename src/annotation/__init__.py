"""Annotation and gene mapping layer for Phase 2."""

from src.annotation.consequence_scoring import score_consequences
from src.annotation.gene_mapping import build_variant_gene_map
from src.annotation.variant_prioritization import build_variant_priority_table
from src.annotation.vep_parser import (
    build_vep_input_file,
    generate_local_vep_command,
    parse_vep_output,
)

__all__ = [
    "build_vep_input_file",
    "generate_local_vep_command",
    "parse_vep_output",
    "score_consequences",
    "build_variant_gene_map",
    "build_variant_priority_table",
]
