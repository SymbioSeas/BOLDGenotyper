"""Query input cleaning (Q1) and alignment metrics (Q2) with the overlap aligner."""

import logging

import pytest

from boldgenotyper import haplotype_query as hq
from tests.synthetic_data import A, base_sequence

_NEXT = {"A": "C", "C": "G", "G": "T", "T": "A"}
Q = A[100:320]                      # 220 bp query lying inside the 600 bp "core"


def _mut(seq, positions):
    s = list(seq)
    for p in positions:
        s[p] = _NEXT[s[p]]
    return "".join(s)


def _align(query, hap, hid="H0001"):
    return hq.align_query_to_haplotype(query, "q", hap, hid)


def test_clean_query_sequence_strips_gaps_and_case():
    assert hq.clean_query_sequence("ac-g.t\n a") == ("ACGTA", True)
    assert hq.clean_query_sequence("ACGT") == ("ACGT", False)


def test_load_query_sequences_cleans_and_warns(tmp_path, caplog):
    fasta = tmp_path / "q.fasta"
    fasta.write_text(">aligned\n---" + Q + "..\n>plain\n" + Q + "\n>short\nACGT\n")
    with caplog.at_level(logging.WARNING, logger="boldgenotyper.haplotype_query"):
        records = hq.load_query_sequences(fasta, min_length=100, max_length=2000)
    assert [r.id for r in records] == ["aligned", "plain"]
    assert all(len(r.seq) == 220 for r in records)
    assert any("looks aligned" in m for m in caplog.messages)


def test_identical_query_inside_core():
    r = _align(Q, A)
    assert (r.identity_pct, r.n_mismatches, r.n_gaps, r.differences) == (100.0, 0, 0, 0)
    assert (r.query_start, r.query_end) == (0, 220)
    assert r.query_coverage == 1.0 and r.overlap_coverage == 1.0
    assert r.haplotype_region == Q
    assert not hasattr(r, "match_quality")


def test_terminal_mismatch_is_counted():
    """v1 local alignment dropped the last base and reported 100%."""
    r = _align(Q, _mut(A, [319]))
    assert r.n_mismatches == 1 and r.n_gaps == 0 and r.differences == 1
    assert r.aligned_length == 220
    assert r.identity_pct == pytest.approx(219 / 220 * 100)


@pytest.mark.parametrize("query,n_gaps,indel", [
    (Q[:110] + Q[111:], 1, True),          # 1-bp deletion
    (Q[:110] + Q[113:], 3, False),         # codon deletion
    (Q[:110] + "A" + Q[110:], 1, True),    # 1-bp insertion
])
def test_internal_indels_are_gaps(query, n_gaps, indel):
    r = _align(query, A)
    assert r.n_gaps == n_gaps and r.n_mismatches == 0
    assert r.indel_in_query is indel
    assert r.differences == n_gaps
    assert r.identity_pct == 100.0          # identity ignores gaps; differences do not


def test_mismatches_and_gaps_counted_separately():
    query = _mut(Q, [20, 150])
    query = query[:80] + query[81:]
    r = _align(query, A)
    assert (r.n_mismatches, r.n_gaps, r.differences) == (2, 1, 3)


def test_query_longer_than_core():
    query = base_sequence(30, seed=9) + A + base_sequence(20, seed=8)
    r = _align(query, A)
    assert r.overlap_coverage == 1.0
    assert r.query_coverage == pytest.approx(600 / 750)
    assert r.differences == 0


def test_unrelated_sequence_has_low_coverage():
    r = _align(base_sequence(73, seed=99), A)
    assert r.overlap_coverage < 0.5




# 2.1.0: a query that runs past the haplotype core is compared over the part the reference
# covers: compared_coverage leaves out query bases beyond the haplotype's own ends, and
# eligibility uses it. overlap_coverage keeps its 2.0.0 meaning.

def _bases(n, seed):
    """n random bases (base_sequence counts codons)."""
    return base_sequence(n // 3 + 1, seed=seed)[:n]


def test_query_overhanging_core_start_is_compared_over_overlap():
    query = _bases(53, 5) + A[:73]          # a 5' mini-barcode: 53 bases before the core
    r = _align(query, A)
    assert r.query_overhang == 53
    assert r.aligned_length == 73 and r.differences == 0
    assert r.compared_coverage == 1.0
    assert r.overlap_coverage == pytest.approx(73 / 126)        # unchanged from 2.0.0
    assert r.query_coverage == pytest.approx(73 / 126)


def test_overhang_past_core_end_counts_too():
    query = A[-80:] + _bases(40, 6)
    r = _align(query, A)
    assert r.query_overhang == 40 and r.compared_coverage == 1.0


def test_full_length_query_compared_coverage_unchanged():
    query = base_sequence(30, seed=9) + A + base_sequence(20, seed=8)
    r = _align(query, A)
    assert r.compared_coverage == r.overlap_coverage == 1.0


def test_bases_opposite_the_haplotype_are_compared_not_overhang():
    """Query bases opposite haplotype bases are aligned (as mismatches if they differ), so
    only bases past the haplotype's own ends can be overhang."""
    r = _align("N" * 30 + A[100:320], A)
    assert (r.query_start, r.haplotype_start, r.query_overhang) == (0, 70, 0)
    assert r.n_mismatches == 30 and r.identity_pct < 97
