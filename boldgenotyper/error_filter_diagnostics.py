"""Where the singleton error filter's assumption may be false.

The filter treats a singleton haplotype close to another haplotype as a likely
sequencing or PCR error. Two kinds of data break that assumption:

- sources that deposit one sequence per haplotype (common for population
  studies mined from GenBank): every true haplotype appears once, so real rare
  haplotypes look like errors;
- family- or order-level data, where an error-filtered singleton may be the
  only record of a species.

These diagnostics report both; they never change a result.
"""

import pandas as pd

from .taxon_evidence import label_source
from .taxon_labels import label_level

SOURCE_COLUMNS = ["source", "n_records", "n_distinct_sequences", "distinct_fraction",
                  "n_assigned", "n_corrected", "n_excluded_error", "haplotype_only"]
HAPLOTYPE_ONLY_MIN_RECORDS = 10
HAPLOTYPE_ONLY_MIN_DISTINCT = 0.90
CORRECTED_SHARE_WARNING = 0.25

FILTERED_STATUSES = ("assigned", "assigned_corrected", "excluded_error_singleton")
SINGLETON_STATUSES = ("assigned_corrected", "excluded_error_singleton")


def _text(records: pd.DataFrame, col: str) -> pd.Series:
    if col not in records.columns:
        return pd.Series("", index=records.index)
    return records[col].fillna("").astype(str)


def source_diagnostics(records: pd.DataFrame) -> pd.DataFrame:
    """One row per label source (taxon_evidence.label_source), over the records
    the error filter saw (assigned, assigned_corrected, excluded_error_singleton).

    Distinct sequences are counted over the deposited sequences, cleaned
    (whitespace and gaps removed, upper case). haplotype_only marks a source
    with at least HAPLOTYPE_ONLY_MIN_RECORDS records of which at least
    HAPLOTYPE_ONLY_MIN_DISTINCT are distinct sequences.
    """
    rec = records[records["status"].isin(FILTERED_STATUSES)]
    if rec.empty:
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    seq = _text(rec, "nuc").str.replace(r"[\s\-.]", "", regex=True).str.upper()
    src = [label_source(p, i, a) for p, i, a in
           zip(_text(rec, "processid"), _text(rec, "inst"), _text(rec, "insdc_acs"))]
    frame = pd.DataFrame({"source": src, "seq": seq.values, "status": rec["status"].values})
    rows = []
    for source, g in frame.groupby("source"):
        n, distinct = len(g), g["seq"].nunique()
        fraction = round(distinct / n, 4)
        rows.append({
            "source": source, "n_records": n, "n_distinct_sequences": distinct,
            "distinct_fraction": fraction,
            "n_assigned": int((g["status"] == "assigned").sum()),
            "n_corrected": int((g["status"] == "assigned_corrected").sum()),
            "n_excluded_error": int((g["status"] == "excluded_error_singleton").sum()),
            "haplotype_only": bool(n >= HAPLOTYPE_ONLY_MIN_RECORDS
                                   and fraction >= HAPLOTYPE_ONLY_MIN_DISTINCT),
        })
    out = pd.DataFrame(rows, columns=SOURCE_COLUMNS)
    return out.sort_values(["n_records", "source"], ascending=[False, True]).reset_index(drop=True)


def label_in_retained(records: pd.DataFrame) -> pd.Series:
    """For error-filtered singleton records (assigned_corrected,
    excluded_error_singleton) with a species-level label: True when an
    assigned record carries the same label, False when none does. NaN for
    other statuses and for records without a species-level label.
    """
    labels = _text(records, "label")
    retained = set(labels[records["status"] == "assigned"])
    out = pd.Series(float("nan"), index=records.index, dtype=object)
    for idx in records.index[records["status"].isin(SINGLETON_STATUSES)]:
        lab = labels[idx]
        if lab and label_level(lab) == "species":
            out[idx] = lab in retained
    return out


def corrected_share(n_assigned: int, n_corrected: int) -> float:
    """Share of assigned records that are corrected singletons (0 when none are assigned)."""
    total = n_assigned + n_corrected
    return n_corrected / total if total else 0.0
