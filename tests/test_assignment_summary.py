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
