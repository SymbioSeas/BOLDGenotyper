"""P0-8: count-free display IDs and the seq_md5 cross-run key."""

import hashlib

from boldgenotyper import haplotype_ids as hid


def test_seq_md5_ignores_case_and_gaps():
    expected = hashlib.md5(b"ACGT").hexdigest()
    assert hid.seq_md5("acg-t") == expected
    assert hid.seq_md5("AC.GT") == expected


def test_assign_display_ids_orders_by_count_then_md5():
    ids = hid.assign_display_ids({"bbb": 5, "aaa": 5, "ccc": 9, "ddd": 1})
    assert ids == {"ccc": "H0001", "aaa": "H0002", "bbb": "H0003", "ddd": "H0004"}


def test_assign_display_ids_is_deterministic_for_identical_input():
    counts = {f"m{i:02d}": i % 3 for i in range(30)}
    assert hid.assign_display_ids(counts) == hid.assign_display_ids(dict(reversed(list(counts.items()))))


def test_assign_display_ids_widens_past_9999():
    ids = hid.assign_display_ids({f"{i:05d}": 1 for i in range(10000)})
    assert ids["00000"] == "H00001"
    assert ids["09999"] == "H10000"


def test_header_round_trip():
    h = hid.format_fasta_header("H0007", "a" * 32, 42)
    assert h == f"H0007 seq_md5={'a' * 32} n_members=42"
    assert hid.parse_fasta_header(h) == {"haplotype_id": "H0007", "seq_md5": "a" * 32, "n_members": 42}


def test_parse_header_tolerates_v1_headers():
    assert hid.parse_fasta_header("h1_n354 Haplotype 1 with 422 sequences") == {
        "haplotype_id": "h1_n354", "seq_md5": None, "n_members": None}
