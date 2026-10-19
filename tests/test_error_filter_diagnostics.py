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


def _groups(rows):
    """rows: (source, label, haplotype_key) per record."""
    return pd.DataFrame([{"processid": f"P{i}", "source": src, "label": lab, "haplotype_key": key}
                         for i, (src, lab, key) in enumerate(rows)])


def test_haplotype_only_is_judged_within_one_species():
    rows = [("GenBank:KC835", "Genus alpha", f"k{i}") for i in range(10)]          # 10 of 10 distinct
    rows += [("SOPNG", f"Genus sp{i}", f"m{i}") for i in range(12)]              # 12 species, 1 each
    rows += [("BOLDP", "Genus beta", f"b{i % 2}") for i in range(10)]            # 2 of 10 distinct
    rows += [("GenBank:JN944", "Genus gamma", f"g{i}") for i in range(5)]        # too few records
    out = efd.haplotype_only_groups(_groups(rows), min_records=10, min_distinct=0.90)
    assert list(out.columns) == efd.GROUP_COLUMNS
    flagged = out[out["haplotype_only"]]
    assert list(zip(flagged["source"], flagged["label"])) == [("GenBank:KC835", "Genus alpha")]
    row = out.set_index(["source", "label"]).loc[("BOLDP", "Genus beta")]
    assert (row["n_records"], row["n_distinct_sequences"]) == (10, 2)


def test_record_groups_use_source_and_label():
    rec = pd.DataFrame({"processid": ["GBMNA001-19", "ESHKD057-07"], "species": ["Genus alpha", ""],
                        "genus": ["Genus", "Genus"], "inst": ["Mined from GenBank, NCBI", "Save Our Seas"],
                        "insdc_acs": ["KC835001", ""]})
    out = efd.record_groups(rec).set_index("processid")
    assert out.loc["GBMNA001-19", "source"] == "GenBank:KC835" and out.loc["GBMNA001-19", "label"] == "Genus alpha"
    assert out.loc["ESHKD057-07", "source"] == "ESHKD" and out.loc["ESHKD057-07", "label"] == "Genus"


def test_sources_table_counts_statuses():
    groups = pd.DataFrame({"source": ["S1"], "label": ["Genus alpha"], "n_records": [3],
                           "n_distinct_sequences": [3], "distinct_fraction": [1.0], "haplotype_only": [False]})
    rec = pd.DataFrame({"processid": ["a", "b", "c"], "source": "S1", "label": "Genus alpha",
                        "status": ["assigned", "assigned_corrected", "excluded_error_singleton"]})
    out = efd.sources_table(groups, rec)
    assert list(out.columns) == efd.SOURCE_COLUMNS
    assert out.iloc[0][["n_assigned", "n_corrected", "n_excluded_error"]].tolist() == [1, 1, 1]


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
