"""Dereplication keyed by seq_md5, error-singleton correction and the
per-sequence core cutoff (P0-8, P0-2, P0-1)."""

import shutil

import pytest

from boldgenotyper import dereplication, haplotype_ids, utils
from tests.synthetic_data import A, B, C, D, mutate

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")


def test_unique_haplotypes_keyed_by_md5_in_v1_order():
    seqs = [B, A, A, B, A, C]
    heads = [f"p{i}" for i in range(len(seqs))]
    records, mapping = dereplication.identify_unique_haplotypes(seqs, heads)
    assert [r.id for r in records] == [haplotype_ids.seq_md5(s) for s in (A, B, C)]
    assert set(mapping.columns) >= {"processid", "haplotype_id", "n_members", "is_singleton"}
    assert mapping.set_index("processid").loc["p1", "haplotype_id"] == haplotype_ids.seq_md5(A)


def test_unique_haplotypes_ties_keep_first_appearance_order():
    seqs = [C, B, C, B]
    records, _ = dereplication.identify_unique_haplotypes(seqs, [f"p{i}" for i in range(4)])
    assert [r.id for r in records] == [haplotype_ids.seq_md5(C), haplotype_ids.seq_md5(B)]


def test_extract_core_region_reports_cutoff_and_exclusions():
    seqs = [A] * 20 + [A[:480] + "-" * 120]
    heads = [f"p{i}" for i in range(21)]
    core, kept, details = utils.extract_core_region(seqs, heads, return_details=True)
    assert "p20" not in kept
    assert details["excluded_headers"] == ["p20"]
    assert details["min_sequence_length_bp"] == pytest.approx(0.95 * details["core_reference_length_bp"])


def test_extract_core_region_cutoff_is_a_parameter():
    seqs = [A] * 20 + [A[:480] + "-" * 120]
    heads = [f"p{i}" for i in range(21)]
    _, kept, details = utils.extract_core_region(seqs, heads, min_sequence_coverage=0.75,
                                                 return_details=True)
    assert "p20" in kept and details["excluded_headers"] == []


def test_extract_core_region_default_return_unchanged():
    out = utils.extract_core_region([A] * 5, [f"p{i}" for i in range(5)])
    assert isinstance(out, tuple) and len(out) == 2


def _write_fasta(path, seqs):
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return path


@needs_mafft
def test_identify_haplotypes_ids_counts_and_corrections(tmp_path):
    seqs = {f"a{i}": A for i in range(10)}
    seqs.update({"a_short1": A[3:], "a_short2": A[3:]})          # merged into A
    seqs.update({f"b{i}": B for i in range(4)})
    seqs.update({f"c{i}": C for i in range(3)})
    seqs.update({"d0": D, "e_corr": mutate(A, [30]), "e_tie": mutate(A, [10])})
    fasta = _write_fasta(tmp_path / "in.fasta", seqs)
    records, mapping, stats, details = dereplication.identify_haplotypes(
        tsv_path=tmp_path / "Synth.tsv", fasta_path=fasta, output_dir=tmp_path / "out",
        min_singleton_distance=0.005, return_details=True)

    by_md5 = stats.set_index("seq_md5")
    a = by_md5.loc[haplotype_ids.seq_md5(A)]
    assert a["haplotype_id"] == "H0001"
    assert (a["n_exact"], a["n_corrected"], a["n_members"]) == (12, 1, 13)
    assert {haplotype_ids.seq_md5(s) for s in (A, B, C, D)} == set(stats["seq_md5"])
    assert all(r.id.startswith("H") and "_n" not in r.id for r in records)
    assert records[0].description == f"seq_md5={haplotype_ids.seq_md5(A)} n_members=13"

    err = details.error_singletons.set_index("processid")
    assert err.loc["e_corr", "status"] == "assigned_corrected"
    assert err.loc["e_corr", "haplotype_id"] == "H0001"
    assert err.loc["e_tie", "status"] == "excluded_error_singleton"
    assert set(err.loc["e_tie", "tied_haplotypes"].split(";")) == {
        by_md5.loc[haplotype_ids.seq_md5(A), "haplotype_id"],
        by_md5.loc[haplotype_ids.seq_md5(B), "haplotype_id"]}
    m = mapping.set_index("processid")
    assert m.loc["e_corr", "assignment"] == "corrected"
    assert m.loc["a_short1", "n_members"] == 13
    assert "e_tie" not in m.index

    fasta_out = (tmp_path / "out" / "Synth_haplotypes.fasta").read_text()
    assert fasta_out.startswith(f">H0001 seq_md5={haplotype_ids.seq_md5(A)} n_members=13\n")


@needs_mafft
def test_no_assign_corrected_leaves_haplotype_set_unchanged(tmp_path):
    seqs = {f"a{i}": A for i in range(5)}
    seqs.update({f"b{i}": B for i in range(3)})
    seqs["e_corr"] = mutate(A, [30])
    fasta = _write_fasta(tmp_path / "in.fasta", seqs)
    on = dereplication.identify_haplotypes(tmp_path / "S.tsv", fasta, tmp_path / "on", return_details=True)
    off = dereplication.identify_haplotypes(tmp_path / "S.tsv", fasta, tmp_path / "off",
                                            assign_corrected=False, return_details=True)
    assert set(on[2]["seq_md5"]) == set(off[2]["seq_md5"])
    assert off[3].error_singletons.set_index("processid").loc["e_corr", "status"] == "excluded_error_singleton"
    assert off[2].set_index("seq_md5").loc[haplotype_ids.seq_md5(A), "n_corrected"] == 0


def test_vectorized_distances_match_compute_distance():
    """assign_error_singletons uses a numpy p-distance; it must equal _compute_distance."""
    import random
    rng = random.Random(3)
    alphabet = "ACGTacgt-N"
    refs = ["".join(rng.choice(alphabet) for _ in range(90)) for _ in range(25)]
    query = "".join(rng.choice(alphabet) for _ in range(90))
    got = dereplication.distances_to_references(query, refs)
    want = [dereplication._compute_distance(query.upper(), r.upper()) for r in refs]
    assert list(got) == want
    assert dereplication.distances_to_references("----", ["ACGT"])[0] == 1.0
