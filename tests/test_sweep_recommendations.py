"""P3: the sweep reports the elbow with the share of records whose haplotype changes."""

import pandas as pd
import pytest

from boldgenotyper import parameter_sweep as ps

THRESHOLDS = [0.005, 0.01, 0.02]


def _tracking():
    # md5 of each record's haplotype per threshold ('' unassigned, 'NA' not in the run)
    return pd.DataFrame({
        "processid": ["r1", "r2", "r3", "r4"],
        "md5_0.005": ["a", "a", "s", "b"],
        "md5_0.01": ["a", "a", "a", "b"],     # r3's singleton corrected into a
        "md5_0.02": ["a", "a", "a", ""],      # r4 left unassigned
    })


def test_adjacent_change_rates():
    rates = ps.adjacent_change_rates(_tracking(), THRESHOLDS)
    assert list(rates.columns) == ["threshold_from", "threshold_to", "n_records", "n_changed", "pct_changed"]
    assert rates[["n_records", "n_changed"]].values.tolist() == [[4, 1], [4, 1]]
    assert rates["pct_changed"].tolist() == [25.0, 25.0]


def test_records_not_in_both_runs_are_left_out():
    t = _tracking()
    t.loc[0, "md5_0.01"] = "NA"
    rates = ps.adjacent_change_rates(t, THRESHOLDS)
    assert rates["n_records"].tolist() == [3, 3]


def test_recommendations_report_elbow_and_stability(tmp_path):
    summary = pd.DataFrame({"threshold": THRESHOLDS, "n_haplotypes": [30, 12, 11],
                            "n_singletons": [20, 3, 2], "n_samples": [100] * 3,
                            "pct_assigned": [90.0, 95.0, 95.0]})
    out = tmp_path / "recommendations.txt"
    ps.generate_recommendations(summary, "Synth", out, tracking_df=_tracking())
    text = out.read_text()
    assert "ACCEPTABLE RANGE" not in text and "AVOID" not in text
    assert "elbow" in text.lower()
    assert "0.005 -> 0.010" in text and "25.0%" in text
    assert "boldgenotyper your_data.tsv" in text and "--input" not in text
