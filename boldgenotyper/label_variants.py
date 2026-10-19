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
from pathlib import Path
from typing import Optional, Tuple

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


# --- Label map --------------------------------------------------------------------

LABEL_MAP_COLUMNS = ["row", "from", "to", "reason", "source"]
_MAP_REQUIRED = ("from", "to", "reason")


class LabelMapError(ValueError):
    """The label map cannot be used."""


def load_label_map(path) -> pd.DataFrame:
    """Read and check a label map: one row per deposited label to replace.

    CSV, or TSV for .tsv/.txt; columns from, to, reason (required) and source
    (optional); extra columns are ignored. Labels pass through clean_label.
    row is the line in the file. A map with an empty field, a label mapped to
    itself, a repeated `from`, or a `to` that is another row's `from` (a
    chain) is rejected as a whole.
    """
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError) as e:
        raise LabelMapError(f"{path}: cannot be read ({e})")
    df.columns = [str(c).strip().lower() for c in df.columns]
    missing = [c for c in _MAP_REQUIRED if c not in df.columns]
    if missing:
        raise LabelMapError(f"{path}: missing column(s): {', '.join(missing)}")
    if "source" not in df.columns:
        df["source"] = ""
    df["row"] = df.index + 2          # line in the file; the header is line 1
    for c in ("from", "to"):
        df[c] = df[c].map(clean_label)
    for c in ("reason", "source"):
        df[c] = df[c].fillna("").astype(str).str.strip()

    def _lines(mask):
        return ", ".join(map(str, df.loc[mask, "row"]))

    problems = []
    empty = (df["from"] == "") | (df["to"] == "") | (df["reason"] == "")
    if empty.any():
        problems.append(f"line(s) {_lines(empty)} need a from, a to and a reason")
    same = (df["from"] == df["to"]) & ~empty
    if same.any():
        problems.append(f"line(s) {_lines(same)} map a label to itself")
    repeated = df["from"].duplicated(keep=False) & ~empty
    if repeated.any():
        problems.append(f"line(s) {_lines(repeated)} repeat a 'from' label")
    chained = df["to"].isin(set(df["from"])) & ~empty & ~same
    if chained.any():
        problems.append(f"line(s) {_lines(chained)} map to a label that is itself mapped (chain)")
    if problems:
        raise LabelMapError(f"{path}: " + "; ".join(problems))
    return df[LABEL_MAP_COLUMNS].reset_index(drop=True)


def apply_label_map(df: pd.DataFrame, label_map: pd.DataFrame, species_col: str = "species",
                    genus_col: str = "genus") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Replace mapped species labels; keep the deposited ones.

    Returns (a copy of df with species_col replaced where the cleaned deposited
    label is a `from`, genus_col set to the first word of `to` when it differs,
    and species_as_deposited holding the cleaned deposited label; a usage
    table with from, to and n_records).
    """
    out = df.copy()
    given = out[species_col] if species_col in out.columns else pd.Series(pd.NA, index=out.index)
    deposited = given.where(given.isna(), given.map(clean_label))
    out["species_as_deposited"] = deposited
    mapping = dict(zip(label_map["from"], label_map["to"]))
    new = deposited.map(lambda v: mapping.get(v, v) if isinstance(v, str) else v)
    changed = deposited.notna() & (new != deposited)
    out.loc[changed, species_col] = new[changed]
    if genus_col in out.columns:
        new_genus = new[changed].map(lambda v: v.split(" ")[0])
        out.loc[changed, genus_col] = new_genus
    counts = deposited[changed].value_counts()
    usage = pd.DataFrame({"from": label_map["from"], "to": label_map["to"],
                          "n_records": [int(counts.get(f, 0)) for f in label_map["from"]]})
    return out, usage
