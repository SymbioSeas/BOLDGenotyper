"""Classify species labels and count species under explicit definitions.

label_level():

    species  two or more words naming a species ("Genus alpha")
    open     open nomenclature ("Genus sp. 1", "Genus cf. alpha", "Genus aff. alpha")
    genus    a single word ("Genus")
    missing  empty or NaN

species_count_summary() returns the species counts defined in the run README:

    n_species_labels_input          unique non-empty species labels among unique input records
    n_species_labels_species_level  of those, labels at level "species"
    n_species_labels_qc_passed      unique non-empty labels among QC-passed records
    n_species_haplotype_assigned    unique species-level assigned_sp across haplotypes
    n_species_confident             unique species among haplotypes whose majority_fraction
                                    is at least the majority threshold
"""

import math
from typing import Dict, Iterable, Optional

import pandas as pd

_OPEN_TOKENS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr.", "gr", "gr."}


def _clean(label) -> str:
    if label is None or (isinstance(label, float) and math.isnan(label)):
        return ""
    text = str(label).strip()
    return "" if text.lower() in ("nan", "none") else text


def label_level(label) -> str:
    """Return 'species', 'open', 'genus' or 'missing' for a species label."""
    parts = _clean(label).split()
    if not parts:
        return "missing"
    if len(parts) == 1:
        return "genus"
    if any(p.lower() in _OPEN_TOKENS for p in parts[1:]):
        return "open"
    return "species"


def genus_of(label) -> Optional[str]:
    """Return the first word of a label, or None for a missing label."""
    parts = _clean(label).split()
    return parts[0] if parts else None


def species_count_summary(input_labels: Iterable, qc_passed_labels: Iterable,
                          haplotype_taxonomy: pd.DataFrame,
                          majority_threshold: float) -> Dict[str, int]:
    """Species counts under the definitions in the module docstring."""
    inputs = {_clean(x) for x in input_labels} - {""}
    qc = {_clean(x) for x in qc_passed_labels} - {""}
    tax = haplotype_taxonomy.copy()
    tax["assigned_sp"] = tax["assigned_sp"].map(_clean)
    species_level = tax[tax["assigned_sp"].map(label_level) == "species"]
    confident = species_level[species_level["majority_fraction"] >= majority_threshold]
    return {
        "n_species_labels_input": len(inputs),
        "n_species_labels_species_level": sum(1 for x in inputs if label_level(x) == "species"),
        "n_species_labels_qc_passed": len(qc),
        "n_species_haplotype_assigned": int(species_level["assigned_sp"].nunique()),
        "n_species_confident": int(confident["assigned_sp"].nunique()),
    }
