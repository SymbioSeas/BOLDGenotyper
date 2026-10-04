"""Query output tables and provenance (Phase 10)."""

import json

import pandas as pd

from boldgenotyper import query_assignment as qa, query_outputs as qo
from tests.test_query_assignment import _call, _m

HAPS = pd.DataFrame({"haplotype_id": ["H0001", "H0002", "H0004"],
                     "haplotype_label": ["Genus alpha H0001", "Genus alpha H0002", "Genus beta H0004"],
                     "assigned_sp": ["Genus alpha", "Genus alpha", "Genus beta"],
                     "assignment_level": ["species"] * 3, "assignment_source": ["labels"] * 3,
                     "n_members": [30, 12, 4]})


def test_calls_table_columns_and_values():
    v = _call([_m("H0002"), _m("H0001"), _m("H0004", d=5, region="X")])
    calls = qo.calls_table([v], reference_run={"provenance_md5": "abc"}, record_matches={"q1": 3})
    assert list(calls.columns) == qo.CALL_COLUMNS
    r = calls.iloc[0]
    assert r["equivalence_set"] == "H0001;H0002" and r["species_call"] == "Genus alpha"
    assert r["call_basis"] == "pooled labels" and r["reference_run"] == "abc"
    assert r["n_exact_record_matches"] == 3 and r["identity_pct"] == 100.0


def test_written_tables_have_no_comment_lines(tmp_path):
    v = _call([_m("H0001")])
    path = qo.write_calls(qo.calls_table([v]), tmp_path)
    assert path.name == "query_calls.csv"
    assert path.read_text().splitlines()[0].startswith("query_id,")
    assert pd.isna(pd.read_csv(path).loc[0, "reference_run"])        # no reference run: empty


def test_candidates_mark_set_members_and_use_query_time_taxa():
    matches = [_m("H0001"), _m("H0002"), _m("H0004", d=9, region="X")]
    v = _call(matches)
    taxa = {"H0002": qa.HaplotypeTaxon("Genus beta", "species", "curated")}
    c = qo.candidates_table(matches, [v], haplotypes=HAPS, haplotype_taxa=taxa)
    assert list(c.columns) == qo.CANDIDATE_COLUMNS and c["rank"].tolist() == [1, 2, 3]
    ci = c.set_index("haplotype_id")
    assert ci.loc["H0001", "in_equivalence_set"] and not ci.loc["H0004", "in_equivalence_set"]
    assert ci.loc["H0001", "eligible"] and not ci.loc["H0004", "eligible"]
    assert ci.loc["H0002", "assigned_sp"] == "Genus beta" and ci.loc["H0002", "assignment_source"] == "curated"
    assert ci.loc["H0004", "n_members"] == 4


def test_query_provenance(tmp_path):
    q = tmp_path / "q.fasta"
    q.write_text(">q1\nACGT\n")
    ref = tmp_path / "haplotypes.fasta"
    ref.write_text(">H0001\nACGT\n")
    prov = qo.build_query_provenance(["boldgenotyper-query", "--query", str(q)], q, 1, ref, 1,
                                     {"provenance_md5": "abc"}, {"max_mismatch_delta": 0})
    (tmp_path / "query_calls.csv").write_text("query_id\n")
    path = qo.write_query_provenance(prov, tmp_path)
    d = json.loads(path.read_text())
    assert d["query"]["n_sequences"] == 1 and len(d["query"]["md5"]) == 32
    assert d["reference_run"]["provenance_md5"] == "abc" and d["reference_fasta"]["n_haplotypes"] == 1
    assert d["parameters"]["max_mismatch_delta"] == 0 and d["alignment"]["gap_open"] == -6
    assert "query_calls.csv" in d["outputs"] and d["warnings"] == []
