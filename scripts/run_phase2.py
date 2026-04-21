#!/usr/bin/env python3
"""Run Phase 2 in scaffold-safe mode, parsing real outputs when present."""

from pathlib import Path
import sys

import polars as pl

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.annotation.consequence_scoring import score_consequences
from src.annotation.gene_mapping import build_variant_gene_map
from src.annotation.variant_prioritization import build_gene_priority_table, build_variant_priority_table
from src.annotation.vep_parser import build_vep_input_file, generate_local_vep_command, parse_vep_output
from src.database import write_phase2_tables
from src.fine_mapping.fine_mapping_pipeline import run_fine_mapping_preparation
from src.io_utils import load_yaml
from src.structure.protein_mapping import build_protein_mapping
from src.structure.structural_summary import build_structure_candidates
from src.target_prioritization.opentargets_api import fetch_opentargets_evidence
from src.target_prioritization.ranking import merge_gene_evidence
from src.target_prioritization.reporting import write_phase2_report
from src.utils.logging_utils import setup_phase2_logging


def main() -> None:
    settings = load_yaml("config/settings.yaml")
    fine_cfg = load_yaml("config/fine_mapping.yaml")
    vep_cfg = load_yaml("config/vep.yaml")
    ot_cfg = load_yaml("config/opentargets.yaml")
    struct_cfg = load_yaml("config/structure.yaml")

    setup_phase2_logging(settings["paths"]["phase2_log_file"])

    # Assembly consistency guard
    project_assembly = settings["project"].get("assembly")
    required = settings["phase2"].get("required_assembly")
    if not project_assembly:
        raise ValueError("config/settings.yaml requires project.assembly for Phase 2")
    if project_assembly != required:
        raise ValueError(f"Assembly mismatch in settings: project.assembly={project_assembly}, required_assembly={required}")

    cleaned = pl.read_csv(fine_cfg["paths"]["cleaned_summary_stats"], separator="\t")
    lead = pl.read_csv(fine_cfg["paths"]["lead_variants"], separator="\t")

    # Fine-mapping scaffold + parser
    fm = run_fine_mapping_preparation(cleaned, lead, fine_cfg)

    # VEP scaffold + parse when available
    vep_input = pl.read_csv(vep_cfg["paths"]["vep_input"], separator="\t")
    build_vep_input_file(vep_input, vep_cfg)
    vep_cmd, vep_status = generate_local_vep_command(vep_cfg)
    cmd_path = Path(vep_cfg["paths"]["vep_raw_dir"]) / "vep_command.sh"
    cmd_path.parent.mkdir(parents=True, exist_ok=True)
    cmd_path.write_text(vep_cmd + "\n", encoding="utf-8")

    annotations = pl.DataFrame()
    consequence_scores = pl.DataFrame()
    raw_vep = Path(vep_cfg["paths"]["vep_raw_output"])
    if raw_vep.exists():
        annotations = parse_vep_output(vep_cfg)
        consequence_scores = score_consequences(annotations, vep_cfg)
        consequence_scores.write_csv(vep_cfg["paths"]["variant_consequence_scores_output"], separator="\t")

    # Variant + gene prioritization
    variant_scores = build_variant_priority_table(lead, consequence_scores, fm["pip_summary"])
    variant_gene_map = build_variant_gene_map(
        variant_scores,
        consequence_scores,
        settings["phase2"]["references"]["genes_dir"],
        "results/annotations/variant_gene_map.tsv",
    )
    gene_priority = build_gene_priority_table(variant_scores, variant_gene_map)

    # Open Targets
    ot = fetch_opentargets_evidence(gene_priority, ot_cfg) if gene_priority.height else pl.DataFrame()
    if ot.height:
        Path(ot_cfg["paths"]["output_tsv"]).parent.mkdir(parents=True, exist_ok=True)
        ot.write_csv(ot_cfg["paths"]["output_tsv"], separator="\t")
        gene_priority = merge_gene_evidence(gene_priority, ot)

    # Structure candidates
    protein_map = build_protein_mapping(gene_priority, struct_cfg) if gene_priority.height else pl.DataFrame()
    struct = build_structure_candidates(variant_scores, protein_map, struct_cfg) if variant_scores.height else {"structure_candidates": pl.DataFrame(), "structural_summary": pl.DataFrame()}

    # Write outputs
    ann_dir = Path(settings["phase2"]["outputs"]["annotations_dir"])
    tp_dir = Path(settings["phase2"]["outputs"]["target_prioritization_dir"])
    st_dir = Path(settings["phase2"]["outputs"]["structure_dir"])
    for d in [ann_dir, tp_dir, st_dir]:
        d.mkdir(parents=True, exist_ok=True)

    if annotations.height:
        annotations.write_csv(ann_dir / "variant_annotations.tsv", separator="\t")
    if consequence_scores.height:
        consequence_scores.write_csv(ann_dir / "variant_consequence_scores.tsv", separator="\t")
    if variant_gene_map.height:
        variant_gene_map.write_csv(ann_dir / "variant_gene_map.tsv", separator="\t")

    variant_scores.write_csv(tp_dir / "variant_priority_scores.tsv", separator="\t")
    gene_priority.write_csv(tp_dir / "gene_prioritization.tsv", separator="\t")

    if protein_map.height:
        protein_map.write_csv(struct_cfg["paths"]["protein_mapping_output"], separator="\t")
    if struct["structure_candidates"].height:
        struct["structure_candidates"].write_csv(struct_cfg["paths"]["structure_candidates_output"], separator="\t")
    if struct["structural_summary"].height:
        struct["structural_summary"].write_csv(struct_cfg["paths"]["structural_summary_output"], separator="\t")

    report_path = write_phase2_report(
        settings["paths"]["phase2_report"],
        variant_scores,
        gene_priority,
        fm["credible_sets"],
        struct["structural_summary"],
    )

    # Missing-input log table for scaffold mode
    missing_rows = [
        {
            "component": "vep",
            "raw_output_present": raw_vep.exists(),
            "vep_in_path": bool(vep_status["vep_in_path"]),
            "vep_cache_present": bool(vep_status["cache_present"]),
        },
        {
            "component": "susie",
            "raw_outputs_present": bool(Path(fine_cfg["paths"]["susie_raw_dir"]).exists()),
            "parsed_pip_rows": int(fm["pip_summary"].height),
            "parsed_cs_rows": int(fm["credible_sets"].height),
        },
    ]
    pl.from_dicts(missing_rows).write_csv(Path(settings["phase2"]["outputs"]["reports_dir"]) / "phase2_missing_inputs.tsv", separator="\t")

    # DuckDB extension
    write_phase2_tables(
        settings["paths"]["database_path"],
        {
            "variant_annotations": annotations,
            "variant_consequence_scores": consequence_scores,
            "variant_gene_map": variant_gene_map,
            "credible_sets": fm["credible_sets"],
            "pip_summary": fm["pip_summary"],
            "variant_priority_scores": variant_scores,
            "gene_prioritization": gene_priority,
            "opentargets_evidence": ot,
            "protein_mapping": protein_map,
            "structure_candidates": struct["structure_candidates"],
        },
    )

    print(f"Phase 2 completed. Report: {report_path}")


if __name__ == "__main__":
    main()
