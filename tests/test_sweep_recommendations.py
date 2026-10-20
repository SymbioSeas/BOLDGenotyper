"""P3: the sweep reports the share of records whose haplotype changes between thresholds.

2.2.0: it recommends the settled threshold (every later step changes fewer than 1% of
records), which does not move with the start of the range; the elbow is kept for reference."""

import re

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


def _summary(thresholds, n_haplotypes):
    n = len(thresholds)
    return pd.DataFrame({"threshold": thresholds, "n_haplotypes": n_haplotypes,
                         "n_singletons": [2] * n, "n_samples": [100] * n,
                         "pct_assigned": [95.0] * n})


def test_recommendations_report_elbow_and_stability_without_settling(tmp_path):
    # every step changes 25% of records: nothing settles, so nothing is recommended
    summary = _summary(THRESHOLDS + [0.03], [30, 12, 11, 10])
    out = tmp_path / "recommendations.txt"
    assert ps.generate_recommendations(summary, "Synth", out, tracking_df=_tracking()) is None
    text = out.read_text()
    assert "ACCEPTABLE RANGE" not in text and "AVOID" not in text
    assert re.search(r"^\s+0\.010\s.*<- elbow$", text, re.M)
    assert "<- settled" not in text
    assert "0.005 -> 0.010" in text and "25.0%" in text
    assert "No threshold settled" in text and "extend the range upward" in text
    assert "boldgenotyper your_data.tsv" in text and "--input" not in text
    assert "--singleton-distance <threshold>" in text


@pytest.mark.parametrize("thresholds,counts", [([0.005], [4]), ([0.005, 0.01], [9, 4]),
                                               ([0.005, 0.01, 0.02], [30, 12, 11])])
def test_no_elbow_below_four_thresholds(tmp_path, thresholds, counts):
    # With 3 thresholds the curvature has one value, so the middle one would always win
    summary = _summary(thresholds, counts)
    assert ps.find_elbow_point(summary) is None
    out = tmp_path / "recommendations.txt"
    assert ps.generate_recommendations(summary, "Synth", out) is None
    text = out.read_text()
    assert "<- elbow" not in text and "<- settled" not in text
    assert "No elbow" in text and "at least 4 thresholds" in text
    assert "No settled threshold" in text and "group_membership_tracking.csv" in text
    assert "--singleton-distance <threshold>" in text


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


def test_sweep_passes_input_flags(tmp_path, monkeypatch):
    from boldgenotyper import cli
    seen = {}

    def fake_pipeline(**kwargs):
        seen.update(kwargs)
        return False

    monkeypatch.setattr(cli, "run_pipeline", fake_pipeline)
    ps.run_single_threshold(tmp_path / "in.tsv", "Synth", 0.01, tmp_path / "run",
                            marker="COI-3P", allow_malformed_rows=True)
    assert seen["marker"] == "COI-3P" and seen["allow_malformed_rows"] is True

    tsv = tmp_path / "in.tsv"
    tsv.write_text("processid\n")
    seen.clear()
    monkeypatch.setattr(ps, "run_parameter_sweep", lambda **kw: seen.update(kw) or {})
    assert cli.main_sweep([str(tsv), "--marker", "COI-3P", "--allow-malformed-rows",
                           "--output", str(tmp_path / "s")]) == 0
    assert seen["marker"] == "COI-3P" and seen["allow_malformed_rows"] is True


def test_sweep_passes_label_map(tmp_path, monkeypatch):
    from boldgenotyper import cli
    seen = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda **kw: seen.update(kw) or False)
    ps.run_single_threshold(tmp_path / "in.tsv", "Synth", 0.01, tmp_path / "run",
                            label_map=tmp_path / "map.csv")
    assert seen["label_map"] == tmp_path / "map.csv"
    tsv = tmp_path / "in.tsv"
    tsv.write_text("processid\n")
    seen.clear()
    monkeypatch.setattr(ps, "run_parameter_sweep", lambda **kw: seen.update(kw) or {})
    assert cli.main_sweep([str(tsv), "--label-map", str(tmp_path / "map.csv"),
                           "--output", str(tmp_path / "s")]) == 0
    assert seen["label_map"] == tmp_path / "map.csv"


# 2.2.0: the settled threshold

SWEEP = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025]
CHANGES = [80, 15, 5, 0, 0]                  # of 1,000 records: 8%, 1.5%, 0.5%, 0%, 0%
COUNTS = [500, 300, 280, 275, 274, 274]      # convex: the elbow is the second threshold


def _tracking_with_changes(thresholds, changes, n=1000):
    """md5 columns where step k changes the haplotype of changes[k] records."""
    data, version, start = {"processid": [f"r{i}" for i in range(n)]}, [0] * n, 0
    for k, t in enumerate(thresholds):
        if k:
            for i in range(start, start + changes[k - 1]):
                version[i] += 1
            start += changes[k - 1]
        data[f"md5_{t}"] = [f"h{i}v{version[i]}" for i in range(n)]
    return pd.DataFrame(data)


def _rates(thresholds, changed, n=1000):
    return pd.DataFrame({"threshold_from": thresholds[:-1], "threshold_to": thresholds[1:],
                         "n_records": [n] * len(changed), "n_changed": changed,
                         "pct_changed": [round(c / n * 100, 2) for c in changed]})


def test_settled_threshold_is_where_every_later_step_stays_under_the_cutoff():
    t = [0.0, 0.005, 0.01, 0.015, 0.02]
    rates = _rates(t, [63, 5, 20, 2])          # 0.005 -> 0.010 is under 1%, but 0.010 -> 0.015 is not
    assert ps.find_settled_threshold(rates) == 0.015
    assert ps.find_settled_threshold(rates, max_pct_changed=2.5) == 0.005


def test_settled_threshold_does_not_move_when_the_first_threshold_is_dropped():
    full = ps.adjacent_change_rates(_tracking_with_changes(SWEEP, CHANGES), SWEEP)
    later = ps.adjacent_change_rates(_tracking_with_changes(SWEEP[1:], CHANGES[1:]), SWEEP[1:])
    assert ps.find_settled_threshold(full) == ps.find_settled_threshold(later) == 0.01
    # the elbow moves with the start of the range
    assert ps.find_elbow_point(_summary(SWEEP, COUNTS)) == 0.005
    assert ps.find_elbow_point(_summary(SWEEP[1:], COUNTS[1:])) == 0.01


def test_settled_share_is_not_rounded():
    under = _rates([0.01, 0.015], [10], n=1004)          # 0.996%, shown as 1.0
    assert under["pct_changed"].tolist() == [1.0]
    assert ps.find_settled_threshold(under) == 0.01
    assert ps.find_settled_threshold(_rates([0.01, 0.015], [10])) is None   # exactly 1%


def test_no_settled_threshold_when_the_last_step_changes_too_many_or_no_steps():
    assert ps.find_settled_threshold(_rates([0.0, 0.005, 0.01], [0, 30])) is None
    assert ps.find_settled_threshold(_rates([0.01], [])) is None


def test_recommendations_recommend_the_settled_threshold_and_caution_on_the_elbow(tmp_path):
    out = tmp_path / "recommendations.txt"
    tracking = _tracking_with_changes(SWEEP, CHANGES)
    assert ps.generate_recommendations(_summary(SWEEP, COUNTS), "Synth", out, tracking_df=tracking) == 0.01
    text = out.read_text()
    assert re.search(r"^\s+0\.010\s.*<- settled$", text, re.M)
    assert re.search(r"^\s+0\.005\s.*<- elbow$", text, re.M)
    assert "fewer than 1% of records" in text
    assert "depends on the range tested" in text and "second threshold tested" in text
    assert "--singleton-distance 0.010" in text
    assert "unevenly spaced" not in text and "sweep lower values" not in text


def test_one_threshold_can_be_settled_and_the_elbow(tmp_path):
    t, out = [0.005, 0.01, 0.015, 0.02, 0.025], tmp_path / "recommendations.txt"
    tracking = _tracking_with_changes(t, [40, 4, 2, 0])
    summary = _summary(t, [300, 280, 275, 274, 274])
    assert ps.generate_recommendations(summary, "Synth", out, tracking_df=tracking) == 0.01
    assert re.search(r"^\s+0\.010\s.*<- settled, elbow$", out.read_text(), re.M)


def test_first_threshold_settled_asks_for_lower_values(tmp_path):
    t, out = [0.01, 0.015, 0.02, 0.025], tmp_path / "recommendations.txt"
    ps.generate_recommendations(_summary(t, [280, 275, 274, 274]), "Synth", out,
                                tracking_df=_tracking_with_changes(t, [2, 0, 0]))
    text = out.read_text()
    assert re.search(r"^\s+0\.010\s.*<- settled", text, re.M)
    assert "sweep lower values" in text


def test_uneven_spacing_is_noted(tmp_path):
    t, out = [0.01, 0.015, 0.02, 0.03, 0.05], tmp_path / "recommendations.txt"   # the default
    ps.generate_recommendations(_summary(t, [280, 275, 274, 274, 273]), "Synth", out,
                                tracking_df=_tracking_with_changes(t, [5, 2, 1, 1]))
    assert "unevenly spaced" in out.read_text()


def test_sweep_recommends_the_settled_threshold_end_to_end(tmp_path, monkeypatch):
    counts = dict(zip(SWEEP, COUNTS))
    monkeypatch.setattr(ps, "run_single_threshold", lambda **kw: {
        "threshold": kw["threshold"], "n_haplotypes": counts[kw["threshold"]], "n_singletons": 2,
        "n_samples": 1000, "pct_assigned": 95.0})
    monkeypatch.setattr(ps, "track_group_membership",
                        lambda runs, organism, thresholds: _tracking_with_changes(SWEEP, CHANGES))
    results = ps.run_parameter_sweep(tmp_path / "Synth.tsv", SWEEP, tmp_path / "sweep", organism="Synth")
    assert results["recommended_threshold"] == 0.01 and results["elbow_threshold"] == 0.005
    readme = (tmp_path / "sweep" / "README.md").read_text()
    assert "**0.010**" in readme and "settled" in readme
    assert (tmp_path / "sweep" / "elbow_plot.pdf").exists()


@pytest.mark.parametrize("result,expected", [
    ({"recommended_threshold": 0.01, "elbow_threshold": 0.005}, "Recommended threshold: 0.010 (settled"),
    ({"recommended_threshold": None, "elbow_threshold": 0.005}, "No threshold settled"),
])
def test_sweep_cli_reports_the_settled_threshold(tmp_path, monkeypatch, capsys, result, expected):
    from boldgenotyper import cli
    tsv = tmp_path / "in.tsv"
    tsv.write_text("processid\n")
    monkeypatch.setattr(ps, "run_parameter_sweep", lambda **kw: result)
    assert cli.main_sweep([str(tsv), "--output", str(tmp_path / "s")]) == 0
    assert expected in capsys.readouterr().out
