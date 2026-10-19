"""Label cleanup and spelling-variant detection (pre-tag A3). Nothing is merged automatically."""

import pandas as pd
import pytest

from boldgenotyper import label_variants as lv


def _records(rows):
    """rows: (label, haplotype_id, n) -> one record per count."""
    out = []
    for label, hap, n in rows:
        for _ in range(n):
            out.append({"processid": f"P{len(out):04d}", "haplotype_id": hap, "label": label})
    return pd.DataFrame(out)


def test_clean_label():
    assert lv.clean_label("  Leptidea   sinapis ") == "Leptidea sinapis"
    assert lv.clean_label("Leptidea sinapis_complex") == "Leptidea sinapis complex"
    assert lv.clean_label("leptidea SINAPIS") == "Leptidea sinapis"
    assert lv.clean_label("Abaies albulaDHJ02") == "Abaies albulaDHJ02"
    assert lv.clean_label(None) == ""


def test_detects_misspelling_and_gender_ending():
    rec = _records([("Leptidea sinapis", "H1", 5), ("Leptidea sinapsis", "H2", 1),
                    ("Colias crocea", "H3", 4), ("Colias croceus", "H3", 2)])
    out = lv.detect_label_variants(rec)
    assert list(out.columns) == lv.LABEL_VARIANT_COLUMNS
    lep = out[out.label_a == "Leptidea sinapis"]
    assert len(lep) == 1 and lep.iloc[0]["label_b"] == "Leptidea sinapsis"
    assert lep.iloc[0]["edit_distance"] == 1 and (lep.iloc[0]["n_a"], lep.iloc[0]["n_b"]) == (5, 1)
    col = out[out.label_a == "Colias crocea"].iloc[0]
    assert col["label_b"] == "Colias croceus" and bool(col["ending_change"]) is True


def test_valid_names_two_letters_apart_are_listed_not_merged():
    rec = _records([("Coregonus hoyi", "H1", 4), ("Coregonus kiyi", "H1", 2)])
    before = rec.copy()
    out = lv.detect_label_variants(rec)
    assert len(out) == 1 and out.iloc[0]["n_shared_haplotypes"] == 1
    pd.testing.assert_frame_equal(rec, before)


def test_qualifier_form():
    rec = _records([("Leptidea sinapis", "H1", 3), ("Leptidea sinapis complex", "H1", 1)])
    out = lv.detect_label_variants(rec)
    assert len(out) == 1 and bool(out.iloc[0]["qualifier_form"]) is True


def test_other_genera_and_distant_epithets_ignored():
    rec = _records([("Pieris napi", "H1", 3), ("Pontia napi", "H2", 3),
                    ("Colias crocea", "H3", 3), ("Colias hyale", "H4", 3)])
    assert lv.detect_label_variants(rec).empty


def test_in_label_map_flag():
    rec = _records([("Leptidea sinapis", "H1", 5), ("Leptidea sinapsis", "H1", 1)])
    label_map = pd.DataFrame({"from": ["Leptidea sinapsis"], "to": ["Leptidea sinapis"]})
    out = lv.detect_label_variants(rec, label_map=label_map)
    assert bool(out.iloc[0]["in_label_map"]) is True
    assert bool(lv.detect_label_variants(rec).iloc[0]["in_label_map"]) is False


def test_deposited_label_used_when_mapped():
    rec = _records([("Leptidea sinapis", "H1", 6)])
    rec["label_mapped_from"] = ["", "", "", "", "", "Leptidea sinapsis"]
    out = lv.detect_label_variants(rec)
    assert out.iloc[0]["label_b"] == "Leptidea sinapsis" and out.iloc[0]["n_b"] == 1


# --- Label map (pre-tag A3) --------------------------------------------------------

def _map(tmp_path, text, name="map.csv"):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_load_label_map_requires_columns(tmp_path):
    with pytest.raises(lv.LabelMapError):
        lv.load_label_map(_map(tmp_path, "from,to\nLeptidea sinapsis,Leptidea sinapis\n"))


def test_label_map_rejects_chains_self_and_duplicates(tmp_path):
    chain = _map(tmp_path, "from,to,reason\nGenus a,Genus b,x\nGenus b,Genus c,y\n", "chain.csv")
    with pytest.raises(lv.LabelMapError, match="line"):
        lv.load_label_map(chain)
    self_map = _map(tmp_path, "from,to,reason\nGenus a,Genus a,x\n", "self.csv")
    with pytest.raises(lv.LabelMapError, match="2"):
        lv.load_label_map(self_map)
    dup = _map(tmp_path, "from,to,reason\nGenus a,Genus b,x\nGenus a,Genus c,y\n", "dup.csv")
    with pytest.raises(lv.LabelMapError, match="3"):
        lv.load_label_map(dup)


def test_label_map_matches_after_cleanup(tmp_path):
    label_map = lv.load_label_map(_map(tmp_path, "from,to,reason\nLeptidea  sinapsis,Leptidea sinapis,typo\n"))
    assert label_map.loc[0, "from"] == "Leptidea sinapsis" and label_map.loc[0, "row"] == 2
    df = pd.DataFrame({"processid": ["P1"], "species": ["leptidea_sinapsis"], "genus": ["Leptidea"]})
    out, usage = lv.apply_label_map(df, label_map)
    assert out.loc[0, "species"] == "Leptidea sinapis"
    assert usage.loc[0, "n_records"] == 1


def test_apply_sets_species_genus_and_keeps_deposited(tmp_path):
    label_map = lv.load_label_map(_map(tmp_path, "from,to,reason,source\nPontia callidice,Pieris callidice,genus moved,backbone\n"))
    df = pd.DataFrame({"processid": ["P1", "P2", "P3"],
                       "species": ["Pontia callidice", "Pontia callidice", "Pontia edusa"],
                       "genus": ["Pontia", "Pontia", "Pontia"]})
    out, usage = lv.apply_label_map(df, label_map)
    assert out["species"].tolist() == ["Pieris callidice", "Pieris callidice", "Pontia edusa"]
    assert out["genus"].tolist() == ["Pieris", "Pieris", "Pontia"]
    assert out["species_as_deposited"].tolist() == ["Pontia callidice", "Pontia callidice", "Pontia edusa"]
    assert usage.set_index("from").loc["Pontia callidice", "n_records"] == 2
    assert df["species"].tolist()[0] == "Pontia callidice"          # input unchanged


def test_unused_map_rows_reported(tmp_path):
    label_map = lv.load_label_map(_map(tmp_path, "from,to,reason\nColias croceus,Colias crocea,gender\n"))
    df = pd.DataFrame({"processid": ["P1"], "species": ["Pieris rapae"], "genus": ["Pieris"]})
    _, usage = lv.apply_label_map(df, label_map)
    assert usage.loc[0, "n_records"] == 0
