"""Synthetic BOLD-portal-style TSV for end-to-end tests.

Sequences are built from codons whose first two bases are not "TA" or "AG".
Any base in the third position of such a codon gives a sense codon under the
vertebrate mitochondrial code (stops: TAA, TAG, AGA, AGG), so every planted
variant passes ORF validation. Variants differ only at third codon positions.

Planted cases (A is a 600 bp stop-free base sequence):

    SYN001-SYN010  A                      exact haplotype
    SYN011-SYN012  A minus first 3 bp     zero-divergence variant, merged into A
    SYN013-SYN016  B = A + codons 10, 20  second haplotype
    SYN017-SYN019  C = A + codons 40-47   1.3% haplotype; one genus-only label
    SYN020         D = A + codons 60-74   distant singleton, retained
    SYN021         A + codon 30           error singleton, nearest A
    SYN022         A + codon 10           error singleton, equidistant from A and B
    SYN023         first 480 bp of A      insufficient core coverage
    SYN024         A with 60 Ns           fails QC (N fraction)
    SYN025         (no sequence)          no sequence
    SYN001 again   duplicate processid row
"""

import csv
import random
from pathlib import Path

_SAFE_CODONS = [
    a + b + c
    for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if (a + b) not in ("TA", "AG")
]
_NEXT_BASE = {"A": "C", "C": "G", "G": "T", "T": "A"}

COLUMNS = [
    "processid", "sampleid", "bin_uri", "species", "genus", "family",
    "country/ocean", "coord", "coord_source", "marker_code", "sex", "life_stage",
    "collection_date_start", "nuc",
]

# Cycled over records: raw BOLD-portal vocabulary variants (P1-5)
_SEX = ["M", "F", " ", "", "U"]
_LIFE_STAGE = ["A", "adult", "I", "", "Immature"]


def base_sequence(n_codons=200, seed=7):
    """Return a reproducible stop-free sequence of n_codons codons."""
    rng = random.Random(seed)
    return "".join(rng.choice(_SAFE_CODONS) for _ in range(n_codons))


def mutate(seq, codon_indices):
    """Change the third base of each listed codon to the next base."""
    s = list(seq)
    for i in codon_indices:
        pos = 3 * i + 2
        s[pos] = _NEXT_BASE[s[pos]]
    return "".join(s)


A = base_sequence()
B = mutate(A, [10, 20])
C = mutate(A, range(40, 48))
D = mutate(A, range(60, 75))

# (processid, sequence, species label)
_RECORDS = [(f"SYN{i:03d}", A, "Genus alpha") for i in range(1, 11)]
_RECORDS += [(f"SYN{i:03d}", A[3:], "Genus alpha") for i in (11, 12)]
_RECORDS += [(f"SYN{i:03d}", B, "Genus beta") for i in range(13, 17)]
_RECORDS += [("SYN017", C, "Genus alpha"), ("SYN018", C, "Genus alpha"), ("SYN019", C, "Genus")]
_RECORDS += [
    ("SYN020", D, "Genus gamma"),
    ("SYN021", mutate(A, [30]), "Genus alpha"),
    ("SYN022", mutate(A, [10]), "Genus beta"),
    ("SYN023", A[:480], "Genus alpha"),
    ("SYN024", A[:270] + "N" * 60 + A[330:], "Genus alpha"),
    ("SYN025", "", "Genus alpha"),
]

EXPECTED = {
    "n_input": len(_RECORDS) + 1,   # one duplicated row
    "n_unique": len(_RECORDS),
    "status": {
        **{f"SYN{i:03d}": "assigned" for i in range(1, 21)},
        "SYN021": "assigned_corrected",
        "SYN022": "excluded_error_singleton",
        "SYN023": "excluded_core_coverage",
        "SYN024": "failed_qc",
        "SYN025": "no_sequence",
    },
    "haplotype_sequences": {A, B, C, D},   # uppercase, ungapped cores
    "members": {"A": 13, "B": 4, "C": 3, "D": 1},   # A: 12 exact + 1 corrected
}

# Cycled over records: ok, ok (GPS), missing, country centroid, zero
_COORDS = [
    ("[12.5, -45.25]", ""),
    ("[-3, 150]", "GPS"),
    ("", ""),
    ("[-20.1, 134.2]", "Coordinates from country centroid"),
    ("[0, 0]", ""),
]


def make_synthetic_bold_tsv(path, duplicate_at=None):
    """
    Write the synthetic TSV to path and return it. The first record is repeated
    at the end, or at row duplicate_at (which leaves a gap in the row index
    after duplicates are removed).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, (pid, seq, label) in enumerate(_RECORDS):
        coord, source = _COORDS[i % len(_COORDS)]
        rows.append({
            "processid": pid,
            "sampleid": f"S-{pid}",
            "bin_uri": "BOLD:SYN0001",
            "species": label if " " in label else "",
            "genus": label.split()[0],
            "family": "Syntheticidae",
            "country/ocean": "Atlantis" if coord else "",
            "coord": coord,
            "coord_source": source,
            "marker_code": "COI-5P",
            "sex": _SEX[i % len(_SEX)],
            "life_stage": _LIFE_STAGE[i % len(_LIFE_STAGE)],
            "collection_date_start": f"20{10 + i % 8:02d}-06-15",
            "nuc": seq,
        })
    # duplicate processid, as in real portal exports
    rows.insert(len(rows) if duplicate_at is None else duplicate_at, dict(rows[0]))
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path
