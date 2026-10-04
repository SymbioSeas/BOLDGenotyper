"""Query verdict (Q3-Q6), reference resolution and the summary writer."""

from pathlib import Path
import pytest
from boldgenotyper import query_assignment as qa


from boldgenotyper.haplotype_query import AlignmentResult

QLEN = 220


def _m(hid, d=0, cov=1.0, region="R", identity=None, aligned=QLEN, indel=False, n_gaps=0):
    """A match as align_query_to_haplotype would return it."""
    ident = identity if identity is not None else (aligned - (d - n_gaps)) / aligned * 100
    return AlignmentResult(
        query_id="q1", haplotype_id=hid, score=0.0, identity_pct=ident,
        matches=aligned - (d - n_gaps), aligned_length=aligned, query_length=QLEN,
        haplotype_length=595, divergence=1 - ident / 100, aligned_query="",
        aligned_haplotype="", n_mismatches=d - n_gaps, n_gaps=n_gaps,
        indel_in_query=indel, query_coverage=cov, overlap_coverage=cov,
        differences=d, haplotype_region=region)


LABELS = {"H0001": ["Genus alpha"] * 30, "H0002": ["Genus alpha"] * 10 + ["Genus"] * 2,
          "H0003": ["Genus alpha"] * 2, "H0004": ["Genus beta"] * 4,
          "H0005": ["Genus alpha"] * 4 + ["Genus beta"] * 4}


def _call(matches, labels=LABELS, **kw):
    return qa.assign_query("q1", matches, labels, **kw)


def test_exact_match_to_identical_haplotypes_is_a_haplotype_set():
    v = _call([_m("H0002"), _m("H0001"), _m("H0003"), _m("H0004", d=5, region="X")])
    assert v.call_level == "haplotype_set"
    assert v.equivalence_set == ["H0001", "H0002", "H0003"]
    assert v.n_equivalent == 3 and v.set_identical_over_region is True
    assert v.species_call == "Genus alpha" and v.species_support == 1.0
    assert v.pooled_member_labels == "Genus alpha 42; Genus 2"
    assert v.flags == []


def test_single_best_is_a_haplotype_call():
    v = _call([_m("H0004", d=2, region="X"), _m("H0001")])
    assert v.call_level == "haplotype" and v.equivalence_set == ["H0001"]
    assert v.best_haplotype == "H0001" and v.differences == 0


def test_terminal_difference_and_delta():
    matches = [_m("H0001"), _m("H0002", d=1, region="Y")]
    assert _call(matches).equivalence_set == ["H0001"]
    v = _call(matches, max_mismatch_delta=1)
    assert v.equivalence_set == ["H0001", "H0002"]
    assert v.set_identical_over_region is False and "set_not_identical" in v.flags


def test_low_coverage_and_low_identity_give_no_call():
    v = _call([_m("H0001", cov=0.6)])
    assert v.call_level == "none" and v.flags == ["low_coverage"] and v.equivalence_set == []
    v = _call([_m("H0001", d=20, region="X")])     # 200/220 = 90.9% identity
    assert v.call_level == "none" and v.flags == ["below_identity"]
    assert v.species_call == "" and v.species_support is None


def test_no_matches():
    v = _call([])
    assert v.call_level == "none" and v.best_haplotype is None and v.best_identity == 0.0


def test_ranking_prefers_eligible_then_fewer_differences():
    ranked = qa.rank_matches([_m("H0001", cov=0.5), _m("H0002", d=1, region="Y"), _m("H0003")],
                             min_identity=0.97, min_query_coverage=0.95)
    assert [m.haplotype_id for m in ranked] == ["H0003", "H0002", "H0001"]


def test_mixed_pool_gives_no_species_call():
    v = _call([_m("H0005")])
    assert v.call_level == "haplotype" and v.species_call == ""
    assert v.species_support == 0.5
    assert set(v.flags) == {"species_support_below_threshold", "label_conflict_in_set"}


def test_conflicting_set_still_calls_a_dominant_species():
    v = _call([_m("H0001"), _m("H0004")])          # 30 alpha + 4 beta
    assert v.species_call == "Genus alpha" and round(v.species_support, 2) == 0.88
    assert "label_conflict_in_set" in v.flags


def test_no_reference_labels():
    v = _call([_m("H0001")], labels=None)
    assert v.call_level == "haplotype" and v.species_call == ""
    assert v.flags == ["no_reference_labels"]


def test_indel_flag_comes_from_the_best_match():
    v = _call([_m("H0001", d=1, n_gaps=1, indel=True)])
    assert "indel_in_query" in v.flags and v.n_gaps == 1


def test_resolve_reference_file_unchanged(tmp_path):
    f = tmp_path / "X_haplotypes.fasta"
    f.write_text(">h1\nACGT\n")
    fasta, adir = qa.resolve_reference(f, None)
    assert fasta == f and adir is None


def test_resolve_reference_directory_finds_fasta_and_sets_analysis_dir(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    f = tmp_path / "haplotypes" / "Carcharhiniformes_haplotypes.fasta"
    f.write_text(">h1\nACGT\n")
    fasta, adir = qa.resolve_reference(tmp_path, None)
    assert fasta == f
    assert adir == tmp_path


def test_resolve_reference_directory_keeps_explicit_analysis_dir(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    (tmp_path / "haplotypes" / "X_haplotypes.fasta").write_text(">h\nA\n")
    other = tmp_path / "other"
    fasta, adir = qa.resolve_reference(tmp_path, other)
    assert adir == other


def test_resolve_reference_directory_zero_fastas_errors(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    with pytest.raises(ValueError):
        qa.resolve_reference(tmp_path, None)


def test_resolve_reference_directory_multiple_fastas_errors(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    (tmp_path / "haplotypes" / "A_haplotypes.fasta").write_text(">h\nA\n")
    (tmp_path / "haplotypes" / "B_haplotypes.fasta").write_text(">h\nA\n")
    with pytest.raises(ValueError):
        qa.resolve_reference(tmp_path, None)


# ---------------------------------------------------------------------------
# Task 3: query_summary.csv output
# ---------------------------------------------------------------------------

def test_write_query_summary_roundtrip(tmp_path):
    import pandas as pd
    from boldgenotyper import haplotype_query
    v = _call([_m("H0002"), _m("H0001")])
    out = haplotype_query.write_query_summary([v], tmp_path, parameters={"max_mismatch_delta": 0})
    df = pd.read_csv(out, comment="#")
    assert list(df["query_id"]) == ["q1"]
    assert df.loc[0, "equivalence_set"] == "H0001;H0002"
    assert df.loc[0, "species_call"] == "Genus alpha"
    assert list(df.columns) == haplotype_query.QUERY_SUMMARY_COLUMNS


# --- Curated haplotype taxa (Phase 9) ----------------------------------------

def test_pooled_call_has_pooled_basis():
    v = _call([_m("H0001")])
    assert v.species_call == "Genus alpha" and v.call_basis == "pooled labels"


def test_curated_set_is_called_by_haplotype_taxa():
    taxa = {"H0005": qa.HaplotypeTaxon("Genus alpha", "species", "curated")}
    v = _call([_m("H0005")], haplotype_taxa=taxa)          # labels 4/4: no pooled call
    assert v.species_call == "Genus alpha" and v.call_basis == "haplotype taxa"
    assert v.species_support == 0.5 and "species_support_below_threshold" not in v.flags
    assert "label_conflict_in_set" in v.flags


def test_curated_taxa_that_disagree_give_no_call():
    taxa = {"H0001": qa.HaplotypeTaxon("Genus alpha", "species", "curated"),
            "H0004": qa.HaplotypeTaxon("Genus beta", "species", "labels")}
    v = _call([_m("H0001"), _m("H0004", region="Y")], haplotype_taxa=taxa)
    assert v.species_call == "" and v.call_basis == ""
    assert "curated_taxa_disagree" in v.flags
    assert "H0001 Genus alpha (curated)" in v.notes and "H0004 Genus beta (labels)" in v.notes


def test_curated_with_a_genus_level_member_gives_no_call():
    taxa = {"H0001": qa.HaplotypeTaxon("Genus alpha", "species", "curated"),
            "H0005": qa.HaplotypeTaxon("Genus", "genus", "labels")}
    v = _call([_m("H0001"), _m("H0005")], haplotype_taxa=taxa)
    assert v.species_call == "" and "curated_taxa_disagree" in v.flags


def test_uncurated_sets_still_use_pooled_labels():
    taxa = {"H0005": qa.HaplotypeTaxon("Genus", "genus", "labels")}
    v = _call([_m("H0005")], haplotype_taxa=taxa)
    assert v.species_call == "" and v.call_basis == ""
    assert "species_support_below_threshold" in v.flags


def test_load_haplotype_taxa(tmp_path):
    assert qa.load_haplotype_taxa(None) is None
    assert qa.load_haplotype_taxa(tmp_path) is None
    (tmp_path / "run_provenance.json").write_text("{}")
    (tmp_path / "haplotypes.csv").write_text(
        "haplotype_id,assigned_sp,assignment_level,assignment_source\n"
        "H0001,Genus alpha,species,labels\nH0002,Genus beta,species,curated\n")
    taxa = qa.load_haplotype_taxa(tmp_path)
    assert taxa["H0002"] == qa.HaplotypeTaxon("Genus beta", "species", "curated")
    (tmp_path / "haplotypes.csv").write_text(
        "haplotype_id,assigned_sp,assignment_level\nH0001,Genus alpha,species\n")
    assert qa.load_haplotype_taxa(tmp_path)["H0001"].source == "labels"
