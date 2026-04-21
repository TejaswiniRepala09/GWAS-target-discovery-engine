"""Mutation string helpers from VEP HGVSp-like fields."""

from __future__ import annotations

import re


def short_protein_change(hgvsp: str | None) -> str | None:
    """Convert p.Gly12Asp -> G12D when parseable."""
    if not hgvsp:
        return None
    m = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)", hgvsp)
    if not m:
        return hgvsp
    aa = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
        "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
        "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
    }
    ref, pos, alt = m.groups()
    return f"{aa.get(ref, ref)}{pos}{aa.get(alt, alt)}"
