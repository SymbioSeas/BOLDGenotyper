"""P0-7: explicit tree membership, with every exclusion listed and explained."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from boldgenotyper import phylogenetics


def _r(hid, n, length=300):
    return SeqRecord(Seq("A" * length), id=hid, description=f"seq_md5={'0' * 32} n_members={n}")


def test_extract_cluster_size_reads_n_members_header():
    assert phylogenetics._extract_cluster_size(_r("H0001", 42)) == (42, "n_members")


def test_extract_cluster_size_reads_parsed_fasta_description():
    rec = SeqRecord(Seq("ACGT"), id="H0002",
                    description=f"H0002 seq_md5={'a' * 32} n_members=7")
    assert phylogenetics._extract_cluster_size(rec) == (7, "n_members")


def test_select_keeps_everything_at_defaults_except_suspects():
    recs = [_r("H0001", 10), _r("H0002", 1), _r("H0003", 2)]
    kept, excluded = phylogenetics.select_tree_haplotypes(
        recs, min_members=1, min_length=150, suspect={"H0003": "Distance outlier (0.120)"})
    assert [r.id for r in kept] == ["H0001", "H0002"]
    assert excluded == [("H0003", "suspect: Distance outlier (0.120)")]


def test_select_reports_every_reason():
    recs = [_r("H0001", 10), _r("H0002", 2, length=100)]
    kept, excluded = phylogenetics.select_tree_haplotypes(recs, min_members=3, min_length=150, suspect={})
    assert [r.id for r in kept] == ["H0001"]
    assert excluded == [("H0002", "n_members 2 < 3; length 100 < 150 bp")]
