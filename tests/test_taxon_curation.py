"""Curation file: load, match to haplotypes, apply."""

import pandas as pd
import pytest

from boldgenotyper import taxon_curation as tc

SEQS = {"H0001": "ACGTACGTAA", "H0002": "ACGTACGTCC", "H0003": "TTTTACGTAA"}
HAPS = pd.DataFrame({"haplotype_id": list(SEQS), "seq_md5": ["m1", "m2", "m3"]})
TAX = pd.DataFrame({"haplotype_id": list(SEQS), "assigned_sp": ["Genus", "Genus beta", "Genus beta"],
                    "assignment_level": ["genus", "species", "species"],
                    "assignment_notes": ["tie for most common species; fell back to genus", "", ""]})


def _write(tmp_path, text, name="cur.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_load_requires_columns(tmp_path):
    with pytest.raises(tc.CurationError):
        tc.load_curation(_write(tmp_path, "seq_md5,taxon\nm1,Genus alpha\n"))


def test_load_requires_values(tmp_path):
    with pytest.raises(tc.CurationError):
        tc.load_curation(_write(tmp_path, "seq_md5,taxon,reason\nm1,Genus alpha,\n"))


def test_load_tsv_and_clean_sequences(tmp_path):
    cur = tc.load_curation(_write(tmp_path, "sequence\ttaxon\treason\nacgt-acgt aa\tGenus alpha\tdistance\n", "c.tsv"))
    assert cur.loc[0, "sequence"] == "ACGTACGTAA" and cur.loc[0, "row"] == 2


def test_match_by_md5_sequence_and_containment(tmp_path):
    cur = tc.load_curation(_write(tmp_path,
        "seq_md5,sequence,taxon,reason\n"
        "m1,,Genus alpha,tree\n"            # md5
        ",ACGTACGTCC,Genus beta,check\n"    # exact sequence
        ",GTACGTCC,Genus beta,check\n"      # contained in H0002 only: same haplotype, same taxon
        ",ACGT,Genus gamma,x\n"             # contained in several haplotypes: not applied
        "zz,,Genus delta,x\n"))             # no match
    applied, not_applied = tc.match_curation(cur, HAPS, SEQS)
    assert sorted(applied["haplotype_id"]) == ["H0001", "H0002"]
    reasons = dict(zip(not_applied["row"], not_applied["reason_not_applied"]))
    assert reasons[5].startswith("matches 3 haplotypes") and reasons[6] == "matches no haplotype"


def test_conflicting_rows_for_one_haplotype_are_not_applied(tmp_path):
    cur = tc.load_curation(_write(tmp_path,
        "seq_md5,taxon,reason\nm1,Genus alpha,a\nm1,Genus beta,b\n"))
    applied, not_applied = tc.match_curation(cur, HAPS, SEQS)
    assert applied.empty and len(not_applied) == 2
    assert not_applied["reason_not_applied"].str.contains("disagree").all()


def test_apply_sets_taxon_level_source_and_note(tmp_path):
    cur = tc.load_curation(_write(tmp_path,
        "seq_md5,taxon,reason,curator,date\nm1,Genus alpha,nested in Genus alpha clade,SS,2026-10-04\n"))
    applied, _ = tc.match_curation(cur, HAPS, SEQS)
    out = tc.apply_curation(TAX, applied).set_index("haplotype_id")
    assert out.loc["H0001", "assigned_sp"] == "Genus alpha"
    assert out.loc["H0001", "assignment_level"] == "species"
    assert out.loc["H0001", "assignment_source"] == "curated"
    assert out.loc["H0001", "assignment_notes"].startswith(
        "curated: nested in Genus alpha clade (SS, 2026-10-04); labels: tie")
    assert out.loc["H0002", "assignment_source"] == "labels"
    assert out.loc["H0002", "assigned_sp"] == "Genus beta"


def test_match_cleans_row_sequences():
    cur = pd.DataFrame([{"row": 2, "seq_md5": "", "sequence": "acgt-acgtcc", "taxon": "Genus beta",
                         "reason": "x", "curator": "", "date": ""}])
    applied, _ = tc.match_curation(cur, HAPS, SEQS)
    assert applied["haplotype_id"].tolist() == ["H0002"]


def test_taxa_not_from_labels():
    tax = pd.DataFrame({"haplotype_id": ["H1", "H2", "H3"], "assigned_sp": ["A a", "B b", "C c"],
                        "assignment_source": ["labels", "curated", "evidence"]})
    assert tc.taxa_not_from_labels(tax) == {"H2": "B b", "H3": "C c"}
    assert tc.taxa_not_from_labels(tax.drop(columns="assignment_source")) == {}
