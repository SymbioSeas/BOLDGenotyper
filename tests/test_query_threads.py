"""2.3.0: boldgenotyper-query runs queries in parallel (--threads); the results are identical."""

import dataclasses
import random

import pytest

from boldgenotyper import haplotype_query

BASE = "".join(random.Random(0).choice("ACGT") for _ in range(300))


def _variants(n, seed):
    """n copies of BASE with 0-6 substitutions each."""
    rng, out = random.Random(seed), []
    for i in range(n):
        s = list(BASE)
        for pos in rng.sample(range(len(BASE)), i % 7):
            s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
        out.append("".join(s))
    return out


def _fastas(tmp_path):
    haps, queries = tmp_path / "haplotypes.fasta", tmp_path / "queries.fasta"
    haps.write_text("".join(f">H{i:04d}\n{s}\n" for i, s in enumerate(_variants(12, 1), 1)))
    queries.write_text("".join(f">q{i}\n{s[20:280]}\n" for i, s in enumerate(_variants(9, 2))))
    return queries, haps


def _plain(items):
    return [dataclasses.asdict(x) for x in items]


def test_threads_give_identical_results(tmp_path):
    queries, haps = _fastas(tmp_path)
    one = haplotype_query.query_against_haplotypes(queries, haps, top_n=5)
    two = haplotype_query.query_against_haplotypes(queries, haps, top_n=5, threads=2)
    assert len(one[2]) == 9 and len(one[0]) == 45
    assert _plain(two[0]) == _plain(one[0])
    assert _plain(two[2]) == _plain(one[2])


def test_more_threads_than_queries(tmp_path):
    queries, haps = _fastas(tmp_path)
    one = haplotype_query.query_against_haplotypes(queries, haps, top_n=3)
    many = haplotype_query.query_against_haplotypes(queries, haps, top_n=3, threads=32)
    assert _plain(many[0]) == _plain(one[0])


def test_threads_must_be_positive(tmp_path):
    queries, haps = _fastas(tmp_path)
    with pytest.raises(ValueError, match="threads"):
        haplotype_query.query_against_haplotypes(queries, haps, threads=0)
