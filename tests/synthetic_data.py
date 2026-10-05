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
E = mutate(A, range(40, 52))   # 4 sites from C, 12 from A: Genus alpha's spread holds C

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

# split_labels (Phase 14): C's labels become Genus alpha 2 / Genus beta 1 (0.67, genus
# level) and a second Genus alpha haplotype E lies beside it, so distance and BIN name
# Genus alpha, C's leading label
_SPLIT_RECORDS = [("SYN026", E, "Genus alpha"), ("SYN027", E, "Genus alpha")]

# Cycled over records: ok, ok (GPS), missing, country centroid, zero
_COORDS = [
    ("[12.5, -45.25]", ""),
    ("[-3, 150]", "GPS"),
    ("", ""),
    ("[-20.1, 134.2]", "Coordinates from country centroid"),
    ("[0, 0]", ""),
]


def make_synthetic_bold_tsv(path, duplicate_at=None, split_labels=False):
    """
    Write the synthetic TSV to path and return it. The first record is repeated
    at the end, or at row duplicate_at (which leaves a gap in the row index
    after duplicates are removed). split_labels: the Phase 14 variant (above).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = list(_RECORDS)
    if split_labels:
        records = [(p, s, "Genus beta" if p == "SYN019" else lab) for p, s, lab in records]
        records += _SPLIT_RECORDS
    rows = []
    for i, (pid, seq, label) in enumerate(records):
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


# ---------------------------------------------------------------------------
# Figure dataset (Phase 11): tables and a tree, no pipeline run
# ---------------------------------------------------------------------------

FIGURE_GROUPS = [("Alpha one", 12), ("Alpha two", 10), ("Beta three", 14), ("Beta four", 1),
                 ("Gamma", 9)]


def figure_dataset():
    """
    (records, haplotypes, newick) for figure tests: 46 haplotypes H0001-H0046.

    Alpha one (12) and Beta three (14) are clades; Alpha two (10) is split in
    two clades of 5, one nested beside Beta three, so it is not monophyletic;
    Beta four has one haplotype; Gamma (9) is a genus-level group. Alpha
    records lie in 'North Basin', Beta in 'South Basin', Gamma in both; a
    haplotype's fourth member has no coordinates. Years 2005-2020.
    """
    import pandas as pd

    ids, hap_rows = {}, []
    k = 0
    for taxon, n in FIGURE_GROUPS:
        ids[taxon] = []
        for _ in range(n):
            k += 1
            hid = f"H{k:04d}"
            ids[taxon].append(hid)
            hap_rows.append({"haplotype_id": hid, "haplotype_label": f"{taxon} {hid}",
                             "assigned_sp": taxon,
                             "assignment_level": "genus" if " " not in taxon else "species",
                             "n_members": 1 + (k - 1) % 4})
    haplotypes = pd.DataFrame(hap_rows)

    rows, counter = [], 0
    for h in hap_rows:
        for j in range(h["n_members"]):
            counter += 1
            taxon = h["assigned_sp"]
            north = taxon.startswith("Alpha") or (taxon == "Gamma" and counter % 2 == 0)
            has_coords = j != 3
            lat = (10 + counter % 10) if north else (-20 + counter % 10)
            rows.append({
                "processid": f"FIG{counter:04d}", "haplotype_id": h["haplotype_id"],
                "haplotype_label": h["haplotype_label"], "assigned_sp": taxon,
                "lat": float(lat) if has_coords else float("nan"),
                "lon": float(-60 + counter % 30) if has_coords else float("nan"),
                "coord_status": "ok" if has_coords else "missing",
                "ocean_basin": ("North Basin" if north else "South Basin") if has_coords else None,
                "sex": ["Male", "Female", ""][counter % 3],
                "life_stage": ["Adult", "Juvenile"][counter % 2],
                "habitat": ["reef", "pelagic", "estuary"][counter % 3],
                "country/ocean": ["Atlantis", "Lemuria"][counter % 2],
                "collection_date_start": f"{2005 + counter % 16}-06-15",
            })
    records = pd.DataFrame(rows)

    def clade(hids, support, length=0.01):
        return "(" + ",".join(f"{h}:0.004" for h in hids) + f"){support}:{length}"

    a2 = ids["Alpha two"]
    alpha = f"({clade(ids['Alpha one'], 0.95)},{clade(a2[:5], 0.9)})0.8:0.03"
    beta = (f"(({clade(ids['Beta three'], 0.97)},{clade(a2[5:], 0.85)})0.6:0.02,"
            f"{ids['Beta four'][0]}:0.03)0.9:0.03")
    newick = f"({alpha},{beta},{clade(ids['Gamma'], 0.99, 0.1)});"
    return records, haplotypes, newick
