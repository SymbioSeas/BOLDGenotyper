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

pool_labels() pools the labels of a set of records (for boldgenotyper-query's
species support): species-level labels count toward support; genus-only and
open-nomenclature labels of the leading species' genus are compatible and left
out of the denominator; genus-only labels of other genera count against;
missing labels are listed separately.
"""

import math
from collections import Counter
from dataclasses import dataclass
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


def haplotype_species_map(haplotype_taxonomy: pd.DataFrame) -> Dict[str, str]:
    """Map each haplotype_id and display label (haplotype_sp) to its assigned_sp."""
    mapping = dict(zip(haplotype_taxonomy["haplotype_id"], haplotype_taxonomy["assigned_sp"]))
    if "haplotype_sp" in haplotype_taxonomy.columns:
        mapping.update(zip(haplotype_taxonomy["haplotype_sp"], haplotype_taxonomy["assigned_sp"]))
    return mapping


MISSING_LABEL_TEXT = "(no species label)"


@dataclass(frozen=True)
class PooledLabels:
    """Pooled labels of a set of records."""
    text: str                         # "Genus alpha 42; Genus beta 4; Genus 3; (no species label) 2"
    leading_species: Optional[str]    # most common species-level label; None if tied or none
    support: Optional[float]          # n_leading / n_species_level; None without species labels
    n_leading: int
    n_species_level: int              # species labels + genus-only labels of other genera
    n_compatible: int                 # genus-only / open labels of the leading genus
    n_missing: int
    n_distinct_species: int


def pool_labels(labels: Iterable) -> PooledLabels:
    """Pool record labels and compute species support (see module docstring)."""
    cleaned = [_clean(x) for x in labels]
    n_missing = sum(1 for x in cleaned if not x)
    counts = Counter(x for x in cleaned if x)
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    parts = [f"{label} {n}" for label, n in ordered]
    if n_missing:
        parts.append(f"{MISSING_LABEL_TEXT} {n_missing}")
    text = "; ".join(parts)

    species = [(label, n) for label, n in ordered if label_level(label) == "species"]
    if not species:
        return PooledLabels(text, None, None, 0, 0, 0, n_missing, 0)

    top_label, n_top = species[0]
    tied = len(species) > 1 and species[1][1] == n_top
    leading_genus = genus_of(top_label)
    n_species_labels = sum(n for _, n in species)
    other = [(label, n) for label, n in ordered if label_level(label) in ("genus", "open")]
    n_compatible = sum(n for label, n in other if genus_of(label) == leading_genus)
    n_conflicting = sum(n for label, n in other if genus_of(label) != leading_genus)
    denominator = n_species_labels + n_conflicting
    return PooledLabels(
        text=text,
        leading_species=None if tied else top_label,
        support=n_top / denominator,
        n_leading=n_top,
        n_species_level=denominator,
        n_compatible=n_compatible,
        n_missing=n_missing,
        n_distinct_species=len(species),
    )
