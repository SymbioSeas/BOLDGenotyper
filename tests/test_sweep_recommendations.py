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


def test_sweep_runs_use_the_given_genetic_code(tmp_path, monkeypatch):
    from boldgenotyper import cli
    seen = {}

    def fake_pipeline(**kwargs):
        seen["code"] = kwargs["cfg"].coi_validation.mitochondrial_code
        return False

    monkeypatch.setattr(cli, "run_pipeline", fake_pipeline)
    ps.run_single_threshold(tmp_path / "in.tsv", "Synth", 0.01, tmp_path / "run", genetic_code=5)
    assert seen["code"] == 5


def test_sweep_cli_passes_genetic_code(tmp_path, monkeypatch):
    from boldgenotyper import cli
    tsv = tmp_path / "in.tsv"
    tsv.write_text("processid\n")
    seen = {}

    def fake_sweep(**kwargs):
        seen.update(kwargs)
        return {}

    monkeypatch.setattr(ps, "run_parameter_sweep", fake_sweep)
    assert cli.main_sweep([str(tsv), "--genetic-code", "5", "--output", str(tmp_path / "s")]) == 0
    assert seen["genetic_code"] == 5
    assert cli.main_sweep([str(tsv), "--output", str(tmp_path / "s")]) == 0
    assert seen["genetic_code"] == 2
