from collections import namedtuple
from pathlib import Path
import pytest
from boldgenotyper import query_assignment as qa

M = namedtuple("M", "haplotype_id identity_pct")


def _matches(*pairs):
    return [M(h, p) for h, p in pairs]


# ---------------------------------------------------------------------------
# Task 1: verdict logic
# ---------------------------------------------------------------------------

def test_confident_single_haplotype_and_species():
    matches = _matches(("h1", 99.6), ("h2", 92.0))
    v = qa.assign_query("q1", matches, {"h1": "Sphyrna lewini"}, 0.97, 0.005)
    assert v.haplotype_call == "confident"
    assert v.assigned_haplotypes == ["h1"]
    assert v.n_tied == 1
    assert v.species_call == "confident"
    assert v.assigned_species == "Sphyrna lewini"
    assert v.candidate_species == []


def test_tied_same_species_is_species_confident():
    matches = _matches(("h1", 99.6), ("h2", 99.6), ("h3", 99.6), ("h4", 90.0))
    sp = {"h1": "Sphyrna lewini", "h2": "Sphyrna lewini", "h3": "Sphyrna lewini"}
    v = qa.assign_query("q1", matches, sp, 0.97, 0.005)
    assert v.haplotype_call == "tied"
    assert v.assigned_haplotypes == ["h1", "h2", "h3"]
    assert v.n_tied == 3
    assert v.species_call == "confident"
    assert v.assigned_species == "Sphyrna lewini"


def test_tied_across_species_is_ambiguous():
    matches = _matches(("h1", 99.6), ("h2", 99.5))
    sp = {"h1": "Sphyrna lewini", "h2": "Sphyrna mokarran"}
    v = qa.assign_query("q1", matches, sp, 0.97, 0.005)
    assert v.haplotype_call == "tied"
    assert v.species_call == "ambiguous"
    assert v.assigned_species is None
    assert v.candidate_species == ["Sphyrna lewini", "Sphyrna mokarran"]


def test_no_confident_match_below_threshold():
    matches = _matches(("h1", 94.2), ("h2", 90.0))
    v = qa.assign_query("q1", matches, {"h1": "Sphyrna lewini"}, 0.97, 0.005)
    assert v.haplotype_call == "no_confident_match"
    assert v.assigned_haplotypes == []
    assert v.species_call == "unassigned"
    assert "94.2" in v.reason


def test_missing_species_metadata_is_unknown():
    matches = _matches(("h1", 99.6),)
    v = qa.assign_query("q1", matches, {}, 0.97, 0.005)
    assert v.haplotype_call == "confident"
    assert v.species_call == "unknown"
    assert v.assigned_species is None


def test_tie_floor_never_dips_below_min_identity():
    # best 97.1, tie_margin 0.5% -> naive floor 96.6, but min_identity 97 wins.
    matches = _matches(("h1", 97.1), ("h2", 96.8))
    v = qa.assign_query("q1", matches, {"h1": "A", "h2": "B"}, 0.97, 0.005)
    assert v.assigned_haplotypes == ["h1"]
    assert v.haplotype_call == "confident"


def test_no_matches_is_no_confident_match():
    v = qa.assign_query("q1", [], {}, 0.97, 0.005)
    assert v.haplotype_call == "no_confident_match"
    assert v.best_identity == 0.0


# ---------------------------------------------------------------------------
# Task 2: reference-by-run-directory resolution
# ---------------------------------------------------------------------------

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
