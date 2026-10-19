"""Species labels: safe formatting cleanup and likely spelling variants.

clean_label tidies how a label is written (whitespace, underscores, the case of
the genus and epithet) without changing the name. detect_label_variants lists
pairs of labels in one genus that may be the same name written differently:
epithets one or two letters apart, a changed gender ending, or one label that
is the other plus more words. It never merges them: valid species can be two
letters apart and even share haplotypes (e.g. Coregonus hoyi and C. kiyi).
A label map the user confirms (load_label_map / apply_label_map) does the
merging, and the run records it.
"""

import re
from typing import Optional

import pandas as pd

from .taxon_labels import label_level

GENDER_ENDINGS = (("us", "a"), ("us", "um"), ("a", "um"), ("is", "e"), ("er", "ra"), ("er", "rum"))
LABEL_VARIANT_COLUMNS = ["label_a", "label_b", "n_a", "n_b", "edit_distance", "ending_change",
                         "qualifier_form", "n_shared_haplotypes", "in_label_map"]
MAX_EDITS = 2


def clean_label(label) -> str:
    """How a label is written, tidied; the name itself is unchanged.

    Whitespace collapsed, underscores to spaces, the genus capitalized and an
    alphabetic epithet in lower case; other tokens (qualifiers, interim codes)
    as given. "" for a missing label.
    """
    if label is None or (isinstance(label, float) and label != label):
        return ""
    text = re.sub(r"\s+", " ", str(label).replace("_", " ")).strip()
    if text.lower() in ("", "nan", "none"):
        return ""
    tokens = text.split(" ")
    if tokens[0].isalpha():
        tokens[0] = tokens[0][0].upper() + tokens[0][1:].lower()
    if len(tokens) > 1 and tokens[1].isalpha():
        tokens[1] = tokens[1].lower()
    return " ".join(tokens)


def _edit_distance(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def _ending_change(a: str, b: str) -> bool:
    for x, y in GENDER_ENDINGS:
        for e1, e2 in ((a, b), (b, a)):
            if e1.endswith(x) and e2.endswith(y) and e1[:-len(x)] == e2[:-len(y)] and e1[:-len(x)]:
                return True
    return False


def detect_label_variants(records: pd.DataFrame,
                          label_map: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Pairs of labels in one genus that may be one name written differently.

    records has processid, haplotype_id and label, and label_mapped_from when a
    label map was applied: the deposited label is label_mapped_from when it is
    not empty, else label. label_map (load_label_map) sets in_label_map for a
    pair when either label is a `from` in the map. Labels are never changed.
    """
    rec = records.copy()
    deposited = rec["label"].fillna("").astype(str)
    if "label_mapped_from" in rec.columns:
        mapped = rec["label_mapped_from"].fillna("").astype(str)
        deposited = mapped.where(mapped != "", deposited)
    rec = rec.assign(_label=deposited.map(clean_label))
    rec = rec[rec["_label"].map(lambda lab: label_level(lab) == "species")]
    counts = rec["_label"].value_counts()
    haps = (rec[rec["haplotype_id"].notna()].groupby("_label")["haplotype_id"]
            .agg(lambda s: set(s)).to_dict()) if "haplotype_id" in rec.columns else {}
    mapped_from = set(label_map["from"].map(clean_label)) if label_map is not None and len(label_map) else set()

    by_genus = {}
    for lab in counts.index:
        by_genus.setdefault(lab.split(" ")[0], []).append(lab)
    rows = []
    for labels in by_genus.values():
        for i, a in enumerate(labels):
            for b in labels[i + 1:]:
                ea, eb = a.split(" ", 1)[1], b.split(" ", 1)[1]
                single = " " not in ea and " " not in eb
                dist = _edit_distance(ea, eb)
                ending = single and _ending_change(ea, eb)
                qualifier = b.startswith(a + " ") or a.startswith(b + " ")
                if not ((single and dist <= MAX_EDITS) or ending or qualifier):
                    continue
                if (counts[b], a) > (counts[a], b):
                    a, b = b, a
                rows.append({
                    "label_a": a, "label_b": b, "n_a": int(counts[a]), "n_b": int(counts[b]),
                    "edit_distance": dist, "ending_change": bool(ending),
                    "qualifier_form": bool(qualifier),
                    "n_shared_haplotypes": len(haps.get(a, set()) & haps.get(b, set())),
                    "in_label_map": a in mapped_from or b in mapped_from,
                })
    out = pd.DataFrame(rows, columns=LABEL_VARIANT_COLUMNS)
    if len(out):
        out = out.sort_values(["n_shared_haplotypes", "edit_distance", "label_a"],
                              ascending=[False, True, True]).reset_index(drop=True)
    return out
