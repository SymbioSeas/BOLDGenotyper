"""The haplotype taxon rule uses the same labeled-member denominator as majority_fraction.

Members without a species label (genus-only labels of the leading genus, open
nomenclature of that genus, and records with no label) are left out of the
denominator; genus-only labels of another genus count against.
"""

import numpy as np
import pandas as pd
import pytest

from boldgenotyper import quality_control as qc
from boldgenotyper import utils


def _assign(groups, threshold=0.7):
    rows = []
    for hid, labels in groups.items():
        for species, genus in labels:
            rows.append({"haplotype_id": hid, "species": species, "genus": genus})
    df = pd.DataFrame(rows)
    assign, _ = utils.assign_consensus_taxonomy(df, group_col="haplotype_id", majority_threshold=threshold)
    return assign.set_index("haplotype_id"), df


A, G = ("Genus alpha", "Genus"), ("", "Genus")


@pytest.mark.parametrize("labels,sp,level,frac", [
    ([A] * 5 + [G] * 5, "Genus alpha", "species", 1.0),               # v1: 0.5 -> genus
    ([A] * 5 + [("", "")] * 5, "Genus alpha", "species", 1.0),        # no label at all
    ([A] * 8 + [("Genus beta", "Genus")] * 2, "Genus alpha", "species", 0.8),
    ([A] * 7 + [("Genus beta", "Genus")] * 3, "Genus", "genus", 0.7),  # not more than 0.70
    ([A] * 4 + [("Genus beta", "Genus")] * 4, "Genus", "genus", 0.5),  # tie
    ([A] * 6 + [("", "Othergenus")] * 2, "Genus alpha", "species", 0.75),
    ([("Genus sp. 1", "Genus")] * 6, "Genus", "genus", np.nan),       # open nomenclature only
])
def test_species_rule(labels, sp, level, frac):
    assign, _ = _assign({"H0001": labels})
    row = assign.loc["H0001"]
    assert (row.assigned_sp, row.assignment_level) == (sp, level)
    if np.isnan(frac):
        assert pd.isna(row.majority_fraction)
    else:
        assert row.majority_fraction == pytest.approx(frac)


def test_majority_fraction_matches_label_columns():
    groups = {"H0001": [A] * 5 + [G] * 3 + [("", "")] * 2,
              "H0002": [A] * 2 + [("Genus beta", "Genus")] * 6 + [("", "Othergenus")]}
    assign, df = _assign(groups)
    labeled = qc.add_label_columns(df, group_col="haplotype_id", notes_col="notes")
    by_label = labeled.groupby("haplotype_id")["majority_fraction"].first()
    for hid in groups:
        assert assign.loc[hid, "majority_fraction"] == pytest.approx(by_label[hid])
    # H0002: beta 6 of 9 labeled (the Othergenus genus-only label counts against): 0.67
    assert (assign.loc["H0002", "assigned_sp"], assign.loc["H0002", "assignment_level"]) == ("Genus", "genus")


def test_notes_say_what_happened():
    assign, _ = _assign({"H0001": [A] * 5 + [G] * 5, "H0002": [A] * 4 + [("Genus beta", "Genus")] * 4,
                         "H0003": [G] * 3})
    assert assign.loc["H0001", "assignment_notes"] == "majority 1.00 (5 of 5 counted labels)"
    assert "tie" in assign.loc["H0002", "assignment_notes"]
    assert "no species-level label" in assign.loc["H0003", "assignment_notes"]


def test_label_rule_reports_its_leading_species():
    assign, _ = _assign({"H1": [("Genus alpha", "Genus")] * 2 + [("Genus beta", "Genus")],
                         "H2": [("Genus alpha", "Genus"), ("Genus beta", "Genus")],
                         "H3": [("", "Genus")] * 2})
    assert assign.loc["H1", "assignment_level"] == "genus"
    assert assign.loc["H1", "leading_species"] == "Genus alpha" and assign.loc["H1", "n_species_labels"] == 2
    assert assign.loc["H2", "leading_species"] == "" and assign.loc["H2", "n_species_labels"] == 2
    assert assign.loc["H3", "leading_species"] == "" and assign.loc["H3", "n_species_labels"] == 0
