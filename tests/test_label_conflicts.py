"""P1-6, P1-7, P1-8: label conflicts with evidence, not 'contamination' or 'confidence'."""

import numpy as np
import pandas as pd
import pytest

from boldgenotyper import quality_control as qc
from tests.synthetic_data import A, mutate


def _records():
    rows = []

    def add(hid, species, genus, n, note=""):
        for _ in range(n):
            rows.append({"processid": f"P{len(rows):03d}", "haplotype_id": hid,
                         "species": species, "genus": genus, "notes": note})

    add("H0001", "Genus alpha", "Genus", 13)
    add("H0001", "", "", 1)                          # unlabeled member (v1 hid it in the denominator)
    add("H0002", "Genus alpha", "Genus", 5)
    add("H0002", "", "Genus", 2)                     # genus-only congeners: compatible
    add("H0003", "Genus alpha", "Genus", 4)
    add("H0003", "Genus beta", "Genus", 4)           # tie: no majority
    add("H0004", "Genus beta", "Genus", 9)
    add("H0004", "Genus alpha", "Genus", 1, note="pup of specimen X-12")
    add(np.nan, "Genus gamma", "Genus", 1)           # unassigned record
    return pd.DataFrame(rows)


SEQS = {"H0001": A, "H0002": mutate(A, [1]), "H0003": mutate(A, [1, 2]),
        "H0004": mutate(A, [1, 2, 3, 4])}


@pytest.fixture
def labeled():
    return qc.add_label_columns(_records(), group_col="haplotype_id")


def test_missing_labels_are_reported_not_hidden(labeled):
    h1 = labeled[labeled.haplotype_id == "H0001"].iloc[0]
    assert (h1.n_labeled, h1.n_missing_label, h1.majority_fraction) == (13, 1, 1.0)
    assert not labeled[labeled.haplotype_id == "H0001"].label_conflict.any()


def test_genus_only_congeners_are_compatible(labeled):
    h2 = labeled[labeled.haplotype_id == "H0002"]
    assert not h2.label_conflict.any()
    assert h2.iloc[0].majority_label == "Genus alpha" and h2.iloc[0].majority_fraction == 1.0
    assert set(h2.label) == {"Genus alpha", "Genus"}


def test_tied_haplotype_has_no_majority_and_all_are_conflicts(labeled):
    h3 = labeled[labeled.haplotype_id == "H0003"]
    assert h3.label_conflict.all() and h3.majority_label.isna().all()
    assert h3.iloc[0].majority_fraction == 0.5


def test_minority_label_is_a_conflict(labeled):
    h4 = labeled[labeled.haplotype_id == "H0004"]
    assert h4.label_conflict.sum() == 1
    assert h4[h4.label_conflict].iloc[0].label == "Genus alpha"
    assert h4.iloc[0].majority_fraction == 0.9


def test_unassigned_records_are_not_conflicts(labeled):
    assert not labeled[labeled.haplotype_id.isna()].label_conflict.any()


def test_no_interpretive_columns(labeled):
    for col in ("potential_misidentification", "misidentification_confidence",
                "uncertainty_category", "group_majority_pct"):
        assert col not in labeled.columns


@pytest.mark.parametrize("note,hits", [
    ("pup of specimen ELAME 463-09", ["pup"]),
    ("pupil diameter measured", []),
    ("unlikely to be relevant", []),
    ("likely juvenile, possibly damaged", ["likely", "juvenile", "possibly", "damaged"]),
    ("identified as Genus cf. alpha", ["cf."]),
    ("", []),
])
def test_keyword_hits_are_whole_words(note, hits):
    assert qc.find_keyword_hits(note) == hits


def test_report_files_and_evidence(tmp_path, labeled):
    res = qc.generate_quality_control_report(labeled, tmp_path, "Synth", group_col="haplotype_id",
                                             haplotype_sequences=SEQS)
    conflicts = pd.read_csv(tmp_path / "label_conflicts.csv")
    assert list(conflicts.columns) == qc.LABEL_CONFLICT_COLUMNS
    assert len(conflicts) == 9                       # 8 in the tie + 1 minority
    alpha_in_h4 = conflicts[(conflicts.haplotype_id == "H0004")].iloc[0]
    assert alpha_in_h4.nearest_haplotype_with_label == "H0003"
    assert alpha_in_h4.distance_to_nearest_with_label == pytest.approx(2 / 600)
    assert alpha_in_h4.depositor_keyword_hits == "pup"

    summary = pd.read_csv(tmp_path / "haplotype_label_summary.csv").set_index("haplotype_id")
    assert list(summary.reset_index().columns) == qc.HAPLOTYPE_LABEL_COLUMNS
    assert summary.loc["H0001", "n_missing_label"] == 1
    assert summary.loc["H0003", "n_label_conflicts"] == 8
    assert summary.loc["H0002", "labels"] == "Genus alpha 5; Genus 2"

    hits = pd.read_csv(tmp_path / "depositor_keyword_hits.csv")
    assert list(hits.columns) == ["processid", "haplotype_id", "label", "depositor_note", "keyword_hits"]
    assert len(hits) == 1

    for old in ("potential_misidentifications.csv", "mixed_species_summary.csv",
                "depositor_flags_summary.csv", "contamination_heatmap.pdf"):
        assert not (tmp_path / old).exists()
    readme = (tmp_path / "README.md").read_text().lower()
    assert "contamination" not in readme and "confidence" not in readme
    assert res["summary"]["n_label_conflicts"] == 9
