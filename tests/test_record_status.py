"""P0-2: one specific, true status per unique input record."""

import pandas as pd
import pytest

from boldgenotyper import record_status as rs


def _inputs():
    unique = ["p1", "p2", "p3", "p4", "p5", "p6"]
    seqs = {"p1", "p2", "p3", "p4", "p5"}
    qc = pd.DataFrame({"processid": ["p1", "p2", "p3", "p4", "p5"],
                       "qc_pass": [True, True, True, True, False],
                       "qc_fail_reason": ["", "", "", "", "N content 10.0% > 5.0%"]})
    exact = pd.DataFrame({"processid": ["p1"], "haplotype_id": ["H0001"]})
    err = pd.DataFrame({"processid": ["p2", "p3"],
                        "status": ["assigned_corrected", "excluded_error_singleton"],
                        "haplotype_id": ["H0001", None],
                        "nearest_haplotype": ["H0001", None],
                        "distance_to_nearest": [0.0017, 0.0017],
                        "tied_haplotypes": ["", "H0001;H0002"]})
    return unique, seqs, qc, exact, ["p4"], err


def test_every_record_gets_exactly_one_status():
    out = rs.build_record_status(*_inputs(), singleton_distance=0.005)
    assert list(out["processid"]) == ["p1", "p2", "p3", "p4", "p5", "p6"]
    assert list(out["status"]) == ["assigned", "assigned_corrected", "excluded_error_singleton",
                                   "excluded_core_coverage", "failed_qc", "no_sequence"]
    assert set(out["status"]) <= set(rs.STATUSES)


def test_notes_are_specific_and_true():
    out = rs.build_record_status(*_inputs(), singleton_distance=0.005).set_index("processid")
    assert "No valid core region" not in " ".join(out["status_note"])
    assert "H0001" in out.loc["p2", "status_note"]
    assert "H0001, H0002" in out.loc["p3", "status_note"]
    assert "N content" in out.loc["p5", "status_note"]
    assert out.loc["p2", "haplotype_id"] == "H0001"
    assert pd.isna(out.loc["p3", "haplotype_id"])
    assert out.loc["p2", "distance_to_nearest"] == 0.0017


def test_correction_disabled_note_names_nearest():
    unique, seqs, qc, exact, core, err = _inputs()
    err.loc[0, ["status", "haplotype_id"]] = ["excluded_error_singleton", None]
    out = rs.build_record_status(unique, seqs, qc, exact, core, err, singleton_distance=0.005)
    note = out.set_index("processid").loc["p2", "status_note"]
    assert "correction disabled" in note and "H0001" in note


def test_unexplained_record_raises():
    unique, seqs, qc, exact, core, err = _inputs()
    with pytest.raises(ValueError, match="p7"):
        rs.build_record_status(unique + ["p7"], seqs | {"p7"}, qc, exact, core, err, 0.005)
