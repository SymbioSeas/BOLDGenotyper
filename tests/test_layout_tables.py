"""P2: records.csv, haplotypes.csv and haplotype_by_region.csv built from in-memory tables."""

import numpy as np
import pandas as pd
import pytest

from boldgenotyper import layout


def _records():
    return pd.DataFrame({
        "processid": ["P1", "P2", "P3", "P4", "P5", "P6"],
        "sampleid": ["S1", "S2", "S3", "S4", "S5", "S6"],
        "species": ["Genus alpha", "Genus alpha", "", "Genus beta", "Genus alpha", "Genus beta"],
        "genus": ["Genus", "Genus", "Genus", "Genus", "Genus", "Genus"],
        "notes": ["", "pup of X", "", "", "", ""],
        "country/ocean": ["Atlantis", "Atlantis", "", "Lemuria", "Atlantis", ""],
        "nuc": ["ACGT"] * 5 + [""],
        "museumid": ["M1"] * 6,                      # not a key field
        "lat": [10.0, 11.0, np.nan, 5.0, 10.0, np.nan],
        "lon": [-50.0, -51.0, np.nan, -40.0, -50.0, np.nan],
        "coord_status": ["ok", "centroid", "missing", "ok", "ok", "missing"],
        "ocean_basin": ["North Atlantic", np.nan, np.nan, "Indian", "North Atlantic", np.nan],
    })


STATUS = pd.DataFrame({
    "processid": ["P1", "P2", "P3", "P4", "P5", "P6"],
    "status": ["assigned", "assigned", "assigned_corrected", "assigned", "failed_qc", "no_sequence"],
    "haplotype_id": ["H0001", "H0001", "H0001", "H0002", np.nan, np.nan],
    "nearest_haplotype": [np.nan, np.nan, "H0001", np.nan, np.nan, np.nan],
    "distance_to_nearest": [np.nan, np.nan, 0.002, np.nan, np.nan, np.nan],
    "status_note": ["", "", "corrected", "", "too short", "no sequence"],
})
QC = pd.DataFrame({"processid": ["P1", "P2", "P3", "P4", "P5"], "qc_pass": [True, True, True, True, False],
                   "qc_fail_reason": ["", "", "", "", "length"], "raw_length": [650, 650, 640, 652, 90]})
TAXONOMY = pd.DataFrame({"haplotype_id": ["H0001", "H0002"], "assigned_sp": ["Genus alpha", "Genus beta"],
                         "assignment_level": ["species", "species"],
                         "assignment_notes": ["majority 0.67", "majority 1.00"],
                         "haplotype_sp": ["Genus alpha H0001", "Genus beta H0002"]})


@pytest.fixture
def records():
    return layout.build_records_table(_records(), STATUS, QC, TAXONOMY, region_col="ocean_basin")


def test_records_one_row_per_unique_record_in_input_order(records):
    assert list(records["processid"]) == ["P1", "P2", "P3", "P4", "P5", "P6"]
    assert list(records["status"]) == list(STATUS["status"])


def test_records_columns(records):
    expected = layout.record_columns("ocean_basin", ["sampleid", "species", "genus", "notes", "country/ocean"])
    assert list(records.columns) == expected
    assert expected[:3] == ["processid", "status", "status_note"]
    assert "nuc" not in records.columns and "museumid" not in records.columns


def test_records_values(records):
    r = records.set_index("processid")
    assert r.loc["P3", "label"] == "Genus" and r.loc["P3", "haplotype_label"] == "Genus alpha H0001"
    assert r.loc["P3", "assigned_sp"] == "Genus alpha"
    assert r.loc["P5", "qc_fail_reason"] == "length" and pd.isna(r.loc["P6", "qc_pass"])
    assert r.loc["P1", "sequence_length"] == 650
    assert r.loc["P2", "depositor_keyword_hits"] == "pup"
    assert not r["label_conflict"].any()
    assert pd.isna(r.loc["P5", "haplotype_id"])


def test_label_conflict_matches_label_columns():
    recs = _records()
    recs.loc[recs.processid == "P2", "species"] = "Genus beta"      # alpha 2, beta 1
    recs.loc[recs.processid == "P3", "species"] = "Genus alpha"
    out = layout.build_records_table(recs, STATUS, QC, TAXONOMY, region_col="ocean_basin")
    assert out.set_index("processid")["label_conflict"].to_dict() == {
        "P1": False, "P2": True, "P3": False, "P4": False, "P5": False, "P6": False}


@pytest.fixture
def haplotypes(records):
    stats = pd.DataFrame({
        "haplotype_id": ["H0001", "H0002"], "seq_md5": ["m1", "m2"], "core_length": [600, 600],
        "n_exact": [2, 1], "n_corrected": [1, 0], "n_members": [3, 1], "is_singleton": [False, True],
        "nearest_haplotype": ["H0002", "H0001"], "nearest_neighbor_distance": [0.01, 0.01],
        "is_suspect": [False, False], "suspect_reasons": ["", ""],
        "min_orf_coverage": [0.99, 0.98], "max_internal_stops": [0, 0]})
    return layout.build_haplotypes_table(stats, TAXONOMY, records, region_col="ocean_basin",
                                         tree_tips={"H0001"})


def test_haplotype_columns_and_merges(haplotypes):
    assert list(haplotypes.columns) == layout.HAPLOTYPE_COLUMNS
    h = haplotypes.set_index("haplotype_id")
    assert h.loc["H0001", "haplotype_label"] == "Genus alpha H0001"
    assert (h.loc["H0001", "n_members"], h.loc["H0001", "n_corrected"]) == (3, 1)
    assert h.loc["H0001", "labels"] == "Genus alpha 2; Genus 1"
    assert h.loc["H0001", "majority_fraction"] == 1.0 and h.loc["H0001", "n_label_conflicts"] == 0
    assert h.loc["H0001", "nearest_haplotype"] == "H0002" and h.loc["H0001", "nearest_distance"] == 0.01
    assert h.loc["H0001", "length"] == 600


def test_haplotype_geographic_counts(haplotypes):
    h = haplotypes.set_index("haplotype_id")
    # P1 ok, P2 centroid (not geographic quality), P3 missing
    assert (h.loc["H0001", "n_with_coords"], h.loc["H0001", "n_with_region"],
            h.loc["H0001", "n_with_country"]) == (1, 1, 2)


def test_in_tree(haplotypes):
    assert haplotypes.set_index("haplotype_id")["in_tree"].to_dict() == {"H0001": True, "H0002": False}
    stats_only = layout.build_haplotypes_table(
        pd.DataFrame({"haplotype_id": ["H0001"], "seq_md5": ["m1"], "core_length": [600], "n_exact": [1],
                      "n_corrected": [0], "n_members": [1], "is_singleton": [True]}),
        TAXONOMY.iloc[:1], layout.build_records_table(_records(), STATUS, QC, TAXONOMY, "ocean_basin"),
        region_col="ocean_basin", tree_tips=None)
    assert stats_only["in_tree"].isna().all()


def test_haplotype_by_region(records):
    by_region = layout.build_haplotype_by_region(records, region_col="ocean_basin")
    assert list(by_region.columns) == ["haplotype_id", "haplotype_label", "assigned_sp", "ocean_basin",
                                       "n_records", "pct_of_region"]
    rows = {(r.haplotype_id, r.ocean_basin): (r.n_records, r.pct_of_region) for r in by_region.itertuples()}
    assert rows == {("H0001", "North Atlantic"): (1, 100.0), ("H0002", "Indian"): (1, 100.0)}
