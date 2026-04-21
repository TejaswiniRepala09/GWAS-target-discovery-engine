"""Target-prioritization integration modules."""

from src.target_prioritization.opentargets_api import fetch_opentargets_evidence
from src.target_prioritization.ranking import merge_gene_evidence
from src.target_prioritization.reporting import write_phase2_report

__all__ = ["fetch_opentargets_evidence", "merge_gene_evidence", "write_phase2_report"]
