"""P0-3: the assignment summary reports explicit denominators."""

import pandas as pd

from boldgenotyper import reports


def test_summary_denominators_and_percentages():
    status = pd.DataFrame({"status": ["assigned"] * 6 + ["assigned_corrected"] * 2 +
                           ["excluded_error_singleton", "excluded_core_coverage",
                            "failed_qc", "no_sequence"]})
    s = reports.compute_assignment_summary(status, n_input=13, n_haplotypes=3,
                                           species_counts={}, n_excluded_error_tied=0,
                                           n_excluded_error_far=1)
    assert s["n_input"] == 13 and s["n_unique"] == 12
    assert s["n_qc_pass"] == 10          # everything except failed_qc and no_sequence
    assert s["n_core_ok"] == 9           # minus excluded_core_coverage
    assert s["pct_assigned"] == 80.0     # (6 + 2) / 10
    assert s["pct_assigned_exact"] == 60.0
    assert s["n_excluded_error"] == 1
    assert (s["n_excluded_error_tied"], s["n_excluded_error_far"]) == (0, 1)
    assert list(s) == reports.SUMMARY_COLUMNS


def test_summary_csv_has_no_empty_columns(tmp_path):
    status = pd.DataFrame({"processid": ["a", "b"], "status": ["assigned", "failed_qc"],
                           "assigned_haplotype": ["H0001", None]})
    diag = tmp_path / "diag.csv"
    status.to_csv(diag, index=False)
    ann = tmp_path / "ann.csv"
    pd.DataFrame({"processid": ["a"], "haplotype_id": ["H0001"]}).to_csv(ann, index=False)
    out = tmp_path / "summary.csv"
    reports.generate_assignment_summary(str(ann), str(diag), str(out), n_input=2, n_haplotypes=1)
    df = pd.read_csv(out)
    assert list(df.columns) == reports.SUMMARY_COLUMNS
    assert df.notna().all(axis=None)
    assert df.loc[0, "pct_assigned"] == 100.0


def test_update_summary_keeps_integer_counts(tmp_path):
    """Pre-tag fix: values filled in Phase 13 leave the counts written in Phase 12 as integers."""
    out = tmp_path / "summary.csv"
    row = {c: 0 for c in reports.SUMMARY_COLUMNS}
    row.update(n_input=9785, n_unique=9221, pct_assigned=83.26)
    pd.DataFrame([row], columns=reports.SUMMARY_COLUMNS).to_csv(out, index=False)
    updated = reports.update_summary(out, {"n_singleton_label_absent": 39})
    written = pd.read_csv(out, dtype=str).iloc[0]
    assert (written["n_input"], written["n_unique"]) == ("9785", "9221")
    assert written["pct_assigned"] == "83.26"
    assert written["n_singleton_label_absent"] == "39"
    assert updated["n_input"] == 9785 and isinstance(updated["n_input"], int)
