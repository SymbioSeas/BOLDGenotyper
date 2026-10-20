"""boldgenotyper-labels: a suggested label map from GBIF or WoRMS (pre-tag A3). No network."""

import json
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd
import pytest

from boldgenotyper import cli, label_backbone as lb
from boldgenotyper.label_variants import load_label_map


class _Response:
    def __init__(self, payload):
        self._body = b"" if payload is None else json.dumps(payload).encode()

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _fake_urlopen(answers):
    """answers: label -> JSON payload (None for an empty body)."""
    def fake(url, timeout=None):
        query = urllib.parse.urlparse(url if isinstance(url, str) else url.full_url).query
        params = urllib.parse.parse_qs(query)
        name = (params.get("name") or params.get("scientificnames[]"))[0]
        return _Response(answers[name])
    return fake


def _run_dir(tmp_path, labels):
    run = tmp_path / "run"
    run.mkdir()
    pd.DataFrame({"processid": [f"P{i}" for i in range(len(labels))], "label": labels}) \
        .to_csv(run / "records.csv", index=False)
    return run


def _gbif(canonical, match_type="EXACT", status="ACCEPTED", species=None):
    return {"canonicalName": canonical, "matchType": match_type, "status": status,
            "species": species or canonical, "rank": "SPECIES"}


@pytest.fixture(autouse=True)
def _no_pause(monkeypatch):
    monkeypatch.setattr(lb, "PAUSE_SECONDS", 0)
    monkeypatch.setattr(lb, "RETRY_PAUSE_SECONDS", 0)


def test_gbif_fuzzy_suggests_spelling(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen({
        "Leptidea sinapsis": _gbif("Leptidea sinapis", "FUZZY"),
        "Leptidea sinapis": _gbif("Leptidea sinapis")}))
    out = lb.suggest_label_map(_run_dir(tmp_path, ["Leptidea sinapsis", "Leptidea sinapis"]))
    assert len(out) == 1
    row = out.iloc[0]
    assert (row["from"], row["to"], row["kind"]) == ("Leptidea sinapsis", "Leptidea sinapis", "spelling")
    assert row["reason"].startswith("GBIF fuzzy match") and row["source"].startswith("GBIF Backbone Taxonomy")


def test_gbif_synonym_marked_synonym(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen({
        "Pontia callidice": _gbif("Pontia callidice", status="SYNONYM", species="Pieris callidice")}))
    out = lb.suggest_label_map(_run_dir(tmp_path, ["Pontia callidice"]))
    assert out.iloc[0]["kind"] == "synonym" and out.iloc[0]["to"] == "Pieris callidice"


def test_accepted_exact_names_not_listed(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen({"Pieris rapae": _gbif("Pieris rapae")}))
    assert lb.suggest_label_map(_run_dir(tmp_path, ["Pieris rapae", "Pieris"])).empty


def test_worms_parsing(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen({
        "Holothuria poli": [[{"match_type": "near_1", "status": "accepted",
                              "scientificname": "Holothuria polii", "valid_name": "Holothuria polii"}]]}))
    out = lb.suggest_label_map(_run_dir(tmp_path, ["Holothuria poli"]), backbone="worms")
    row = out.iloc[0]
    assert (row["from"], row["to"], row["kind"]) == ("Holothuria poli", "Holothuria polii", "spelling")
    assert row["source"].startswith("WoRMS")


def test_network_failure_writes_nothing(tmp_path, monkeypatch, capsys):
    def offline(url, timeout=None):
        raise urllib.error.URLError("no network")
    monkeypatch.setattr(urllib.request, "urlopen", offline)
    run = _run_dir(tmp_path, ["Pieris rapae"])
    assert cli.main_labels(["suggest", str(run)]) == 1
    assert not (run / "label_map_suggested.csv").exists()
    assert "no network" in capsys.readouterr().err


def test_suggestions_load_as_a_label_map(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen({
        "Leptidea sinapsis": _gbif("Leptidea sinapis", "FUZZY")}))
    run = _run_dir(tmp_path, ["Leptidea sinapsis"])
    assert cli.main_labels(["suggest", str(run)]) == 0
    loaded = load_label_map(run / "label_map_suggested.csv")
    assert loaded.loc[0, "from"] == "Leptidea sinapsis" and loaded.loc[0, "to"] == "Leptidea sinapis"


def _worms(scientificname, valid_name, status, match_type="exact"):
    return [[{"match_type": match_type, "status": status, "scientificname": scientificname,
              "valid_name": valid_name}]]


def test_worms_subgenus_form_is_not_a_synonym(tmp_path, monkeypatch):
    """WoRMS gives accepted names with the subgenus and marks the plain binomial an
    'alternative representation': the same name, not a synonym (as WoRMS answered, 20 Oct 2026)."""
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen({
        "Holothuria albiventer": _worms("Holothuria albiventer", "Holothuria (Metriatyla) albiventer",
                                        "alternative representation"),
        "Holothuria polii": _worms("Holothuria polii", "Holothuria (Roweothuria) poli",
                                   "alternative representation", "phonetic"),
        "Holothuria collaris": _worms("Holothuria collaris", "Holothuria (Lessonothuria) verrucosa",
                                      "unaccepted"),
        "Holothuria ocellata": _worms("Holothuria ocellata", "Bohadschia ocellata", "unaccepted",
                                      "exact_subgenus"),
    }))
    labels = ["Holothuria albiventer", "Holothuria polii", "Holothuria collaris", "Holothuria ocellata"]
    out = lb.suggest_label_map(_run_dir(tmp_path, labels), backbone="worms").set_index("from")
    assert "Holothuria albiventer" not in out.index
    assert out.loc["Holothuria polii", ["to", "kind"]].tolist() == ["Holothuria poli", "spelling"]
    assert out.loc["Holothuria collaris", ["to", "kind"]].tolist() == ["Holothuria verrucosa", "synonym"]
    assert out.loc["Holothuria ocellata", ["to", "kind"]].tolist() == ["Bohadschia ocellata", "synonym"]


def test_strip_subgenus():
    assert lb.strip_subgenus("Holothuria (Metriatyla) albiventer") == "Holothuria albiventer"
    assert lb.strip_subgenus("Holothuria albiventer") == "Holothuria albiventer"


def test_a_timeout_is_retried(tmp_path, monkeypatch):
    """One slow answer among hundreds of names should not end the run."""
    calls = []
    answer = _fake_urlopen({"Pieris rapae": _gbif("Pieris rapae")})

    def flaky(url, timeout=None):
        calls.append(url)
        if len(calls) <= 2:
            raise urllib.error.URLError("timed out")
        return answer(url, timeout)
    monkeypatch.setattr(urllib.request, "urlopen", flaky)
    lb.suggest_label_map(_run_dir(tmp_path, ["Pieris rapae"]))
    assert len(calls) == 3


def test_persistent_failure_still_fails(tmp_path, monkeypatch):
    calls = []

    def offline(url, timeout=None):
        calls.append(url)
        raise urllib.error.URLError("no network")
    monkeypatch.setattr(urllib.request, "urlopen", offline)
    with pytest.raises(lb.BackboneError):
        lb.suggest_label_map(_run_dir(tmp_path, ["Pieris rapae"]))
    assert len(calls) == lb.ATTEMPTS
