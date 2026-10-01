"""P1-2: barcoding gap only with --extended, computed from record labels."""

import numpy as np
import pandas as pd
import pytest

from boldgenotyper import divergence_analysis as da

IDS = ["H0001", "H0002", "H0003", "H0004"]
D = pd.DataFrame([[0.00, 0.01, 0.05, 0.06],
                  [0.01, 0.00, 0.04, 0.05],
                  [0.05, 0.04, 0.00, 0.02],
                  [0.06, 0.05, 0.02, 0.00]], index=IDS, columns=IDS)


def _records(rows):
    return pd.DataFrame(rows, columns=["haplotype_id", "label", "label_conflict"])


def test_separate_labels_have_a_gap():
    rec = _records([("H0001", "Genus alpha", False)] * 3 + [("H0002", "Genus alpha", False)]
                   + [("H0003", "Genus beta", False)] * 2)
    gap = da.compute_label_based_gap(D, rec).set_index("label")
    assert list(gap.reset_index().columns) == da.LABEL_GAP_COLUMNS
    a = gap.loc["Genus alpha"]
    assert (a.n_records, a.n_haplotypes, a.max_within, a.min_between) == (4, 2, 0.01, 0.04)
    assert a.nearest_other_label == "Genus beta" and bool(a.label_based_gap) is True
    assert gap.loc["Genus beta", "max_within"] == 0.0


def test_label_conflict_records_are_excluded():
    rec = _records([("H0001", "Genus alpha", False)] * 2 + [("H0003", "Genus beta", False)] * 9
                   + [("H0003", "Genus alpha", True)])          # minority label in H0003
    gap = da.compute_label_based_gap(D, rec).set_index("label")
    assert gap.loc["Genus alpha", "n_haplotypes"] == 1
    assert gap.loc["Genus alpha", "min_between"] == 0.05


def test_shared_haplotype_without_conflict_flags_has_no_gap():
    rec = pd.DataFrame({"haplotype_id": ["H0001", "H0001", "H0002"],
                        "label": ["Genus alpha", "Genus beta", "Genus alpha"]})
    gap = da.compute_label_based_gap(D, rec).set_index("label")
    assert gap.loc["Genus alpha", "min_between"] == 0.0
    assert bool(gap.loc["Genus alpha", "label_based_gap"]) is False


def test_genus_only_and_missing_labels_are_not_rows():
    rec = _records([("H0001", "Genus", False), ("H0002", "", False), ("H0003", "Genus beta", False)])
    gap = da.compute_label_based_gap(D, rec)
    assert list(gap["label"]) == ["Genus beta"]
    assert np.isnan(gap.loc[0, "min_between"])
