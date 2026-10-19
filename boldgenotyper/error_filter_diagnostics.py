"""Where the singleton error filter's assumption fails, and the groups it spares.

The filter treats a singleton haplotype close to another haplotype as a likely
sequencing or PCR error. Two kinds of data break that assumption:

- studies that deposit one sequence per haplotype (common for population
  studies mined from GenBank): every true haplotype appears once, so real rare
  haplotypes look like errors. A group of records with one source and one
  label (haplotype_only_groups) is haplotype-only when it has at least
  min_records records of which at least min_distinct are distinct core
  sequences; the filter does not treat its singletons as errors
  (dereplication.identify_haplotypes). The source is the BOLD project
  (processid prefix) or, for records mined from GenBank, the accession prefix;
  judging within one label keeps multi-species barcode-library projects, whose
  sequences differ because their species differ, from being counted.
- family- or order-level data, where an error-filtered singleton may be the
  only record of a species (label_in_retained).
"""

from typing import Optional

import pandas as pd

from .taxon_labels import label_level

GROUP_COLUMNS = ["source", "label", "n_records", "n_distinct_sequences", "distinct_fraction",
                 "haplotype_only"]
SOURCE_COLUMNS = GROUP_COLUMNS + ["n_assigned", "n_corrected", "n_excluded_error"]
HAPLOTYPE_ONLY_MIN_RECORDS = 10          # defaults of config.HaplotypeConfig
HAPLOTYPE_ONLY_MIN_DISTINCT = 0.90
FILTERED_SHARE_WARNING = 0.25

FILTERED_STATUSES = ("assigned", "assigned_corrected", "excluded_error_singleton")
SINGLETON_STATUSES = ("assigned_corrected", "excluded_error_singleton")


def _text(records: pd.DataFrame, col: str) -> pd.Series:
    if col not in records.columns:
        return pd.Series("", index=records.index)
    return records[col].fillna("").astype(str)


def record_groups(records: pd.DataFrame) -> pd.DataFrame:
    """processid, source and label of each record (label: species, else genus)."""
    from .quality_control import record_labels
    from .taxon_evidence import label_source
    labels = record_labels(records) if len(records) else pd.Series(dtype=object)
    return pd.DataFrame({
        "processid": _text(records, "processid").values,
        "source": [label_source(p, i, a) for p, i, a in
                   zip(_text(records, "processid"), _text(records, "inst"), _text(records, "insdc_acs"))],
        "label": labels.fillna("").values if len(records) else [],
    })


def haplotype_only_groups(groups: pd.DataFrame, min_records: int = HAPLOTYPE_ONLY_MIN_RECORDS,
                          min_distinct: float = HAPLOTYPE_ONLY_MIN_DISTINCT) -> pd.DataFrame:
    """One row per source and label (GROUP_COLUMNS).

    groups has processid, source, label and haplotype_key (each record's
    haplotype before the error filter, i.e. its distinct core sequence).
    """
    if groups.empty:
        return pd.DataFrame(columns=GROUP_COLUMNS)
    rows = []
    for (source, label), g in groups.groupby(["source", "label"]):
        n, distinct = len(g), g["haplotype_key"].nunique()
        fraction = round(distinct / n, 4)
        rows.append({"source": source, "label": label, "n_records": n,
                     "n_distinct_sequences": distinct, "distinct_fraction": fraction,
                     "haplotype_only": bool(n >= min_records and fraction >= min_distinct)})
    out = pd.DataFrame(rows, columns=GROUP_COLUMNS)
    return out.sort_values(["haplotype_only", "n_records", "source", "label"],
                           ascending=[False, False, True, True]).reset_index(drop=True)


def sources_table(groups: pd.DataFrame, records: pd.DataFrame) -> pd.DataFrame:
    """sources.csv: the groups with their records' statuses after the filter.

    records has processid, source, label and status.
    """
    if groups.empty:
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    status = records[["source", "label", "status"]]
    counts = (status.assign(n_assigned=status["status"] == "assigned",
                            n_corrected=status["status"] == "assigned_corrected",
                            n_excluded_error=status["status"] == "excluded_error_singleton")
              .groupby(["source", "label"])[["n_assigned", "n_corrected", "n_excluded_error"]].sum()
              .astype(int).reset_index())
    out = groups.merge(counts, on=["source", "label"], how="left")
    for col in ("n_assigned", "n_corrected", "n_excluded_error"):
        out[col] = out[col].fillna(0).astype(int)
    return out[SOURCE_COLUMNS]


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


def filtered_share(n_assigned: int, n_corrected: int, n_excluded: int) -> float:
    """Share of the records the filter saw that it treated as errors (corrected or
    left unassigned); 0 when it saw none."""
    total = n_assigned + n_corrected + n_excluded
    return (n_corrected + n_excluded) / total if total else 0.0


def exempt_processids(groups_by_record: pd.DataFrame, groups: pd.DataFrame) -> set:
    """processids of records in haplotype-only groups (groups_by_record: processid, source, label)."""
    flagged = groups.loc[groups["haplotype_only"], ["source", "label"]]
    if flagged.empty:
        return set()
    hit = groups_by_record.merge(flagged, on=["source", "label"])
    return set(hit["processid"])


def summary_values(groups: Optional[pd.DataFrame], n_exempt: int) -> dict:
    """The summary.csv values about haplotype-only groups."""
    n_groups = int(groups["haplotype_only"].sum()) if groups is not None and len(groups) else 0
    return {"n_haplotype_only_groups": n_groups, "n_singletons_exempt": int(n_exempt)}
