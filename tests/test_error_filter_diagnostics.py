"""Where the singleton error filter's assumption may be false (pre-tag A4)."""

import numpy as np
import pandas as pd

from boldgenotyper import error_filter_diagnostics as efd

BASE = "ACGT" * 50


def _seq(i):
    """Distinct sequences: one changed base per index."""
    s = list(BASE)
    s[i] = "T" if s[i] != "T" else "A"
    return "".join(s)


def _records(source, n, n_distinct, status="assigned"):
    return pd.DataFrame({
        "processid": [f"{source}{i:03d}-19" for i in range(n)],
        "status": status, "label": "Genus alpha", "inst": "Mined from GenBank, NCBI",
        "insdc_acs": [f"{source[2:]}{i:06d}" for i in range(n)],
        "nuc": [_seq(i % n_distinct) for i in range(n)],
    })


def test_haplotype_only_source_flagged():
    rec = pd.concat([_records("GBMNA", 10, 10), _records("GBMNB", 10, 2), _records("GBMNC", 5, 5)],
                    ignore_index=True)
    out = efd.source_diagnostics(rec)
    assert list(out.columns) == efd.SOURCE_COLUMNS
    flags = dict(zip(out["n_records"].astype(str) + "/" + out["n_distinct_sequences"].astype(str),
                     out["haplotype_only"]))
    assert flags == {"10/10": True, "10/2": False, "5/5": False}


def test_source_diagnostics_counts_only_filtered_statuses():
    rec = pd.concat([_records("GBMNA", 10, 10),
                     _records("GBMNX", 4, 4, status="failed_qc")], ignore_index=True)
    out = efd.source_diagnostics(rec)
    assert out["n_records"].sum() == 10


def test_label_in_retained():
    rec = pd.DataFrame({
        "processid": ["a", "b", "c", "d", "e"],
        "status": ["assigned", "assigned_corrected", "assigned_corrected", "assigned",
                   "excluded_error_singleton"],
        "label": ["Genus beta", "Genus beta", "Genus delta", "Genus alpha", "Genus"],
    })
    out = efd.label_in_retained(rec)
    assert out.index.equals(rec.index)
    values = out.tolist()
    assert values[1] is True and values[2] is False
    assert np.isnan(values[0]) and np.isnan(values[3])
    assert np.isnan(values[4])          # genus-only label: cannot tell
