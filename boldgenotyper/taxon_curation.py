"""Curated haplotype taxa: accepted taxon calls applied after the label rule.

The label rule (utils.assign_consensus_taxonomy) assigns each haplotype the
species most of its member labels name, or the genus. A curation file records
calls a user has checked against other evidence (see taxon_evidence) and wants
to apply instead, with the reason.

File (CSV, or TSV for .tsv/.txt), one row per curated haplotype:

    seq_md5 or sequence   key of the haplotype (one of the two is required)
    taxon                 the taxon to assign (required)
    reason                why (required)
    curator, date         optional

Matching: a row applies to the haplotype whose seq_md5 equals the row's; else
to the haplotype whose core sequence equals the row's sequence; else to the
haplotype whose core sequence contains or is contained in it (so a call made on
one dataset can apply to another dataset with a different core region). A row
matching no haplotype or several is not applied. Rows naming the same
haplotype are applied once when their taxa agree; otherwise none is applied.
"""

from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

from .taxon_labels import label_level

REQUIRED = ("taxon", "reason")
KEYS = ("seq_md5", "sequence")
OPTIONAL = ("curator", "date")
CURATION_COLUMNS = ["row", "seq_md5", "sequence", "taxon", "reason", "curator", "date"]


class CurationError(ValueError):
    """The curation file cannot be used."""


def _clean_sequence(seq) -> str:
    return "".join(str(seq).split()).upper().replace("-", "").replace(".", "")


def load_curation(path) -> pd.DataFrame:
    """Read and validate a curation file; one row per curated haplotype."""
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError) as e:
        raise CurationError(f"{path}: cannot be read ({e})")
    df.columns = [str(c).strip().lower() for c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if not any(k in df.columns for k in KEYS):
        missing.append("seq_md5 or sequence")
    if missing:
        raise CurationError(f"{path}: missing column(s): {', '.join(missing)}")
    for c in KEYS + OPTIONAL:
        if c not in df.columns:
            df[c] = ""
    for c in KEYS + REQUIRED + OPTIONAL:
        df[c] = df[c].fillna("").astype(str).str.strip()
    df["sequence"] = df["sequence"].map(lambda s: _clean_sequence(s) if s else "")
    df["seq_md5"] = df["seq_md5"].str.lower()
    df["row"] = df.index + 2          # line in the file; the header is line 1
    incomplete = df[(df["taxon"] == "") | (df["reason"] == "")
                    | ((df["seq_md5"] == "") & (df["sequence"] == ""))]
    if len(incomplete):
        raise CurationError(f"{path}: row(s) {', '.join(map(str, incomplete['row']))} need "
                            f"a taxon, a reason and a seq_md5 or sequence")
    return df[CURATION_COLUMNS].reset_index(drop=True)


def match_curation(curation: pd.DataFrame, haplotypes: pd.DataFrame,
                   sequences: Dict[str, str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match curation rows to haplotypes.

    haplotypes has haplotype_id and seq_md5; sequences maps haplotype_id to its
    core sequence. Returns (applied rows with haplotype_id, rows not applied with
    reason_not_applied).
    """
    md5_to_id = dict(zip(haplotypes["seq_md5"].astype(str).str.lower(), haplotypes["haplotype_id"]))
    cores = {h: _clean_sequence(s) for h, s in sequences.items()}
    matched, not_applied = [], []
    for r in curation.to_dict("records"):
        if r["seq_md5"] and r["seq_md5"] in md5_to_id:
            hits = [md5_to_id[r["seq_md5"]]]
        elif r["sequence"]:
            hits = [h for h, s in cores.items() if s == r["sequence"]]
            if not hits:
                hits = [h for h, s in cores.items()
                        if s and (r["sequence"] in s or s in r["sequence"])]
        else:
            hits = []
        if len(hits) == 1:
            matched.append({**r, "haplotype_id": hits[0]})
        else:
            why = ("matches no haplotype" if not hits else
                   f"matches {len(hits)} haplotypes ({';'.join(sorted(hits))})")
            not_applied.append({**r, "reason_not_applied": why})

    applied = []
    by_haplotype: Dict[str, list] = {}
    for r in matched:
        by_haplotype.setdefault(r["haplotype_id"], []).append(r)
    for hid, rows in by_haplotype.items():
        if len({r["taxon"] for r in rows}) == 1:
            applied.append(rows[0])
        else:
            for r in rows:
                not_applied.append({**{k: v for k, v in r.items() if k != "haplotype_id"},
                                    "reason_not_applied": f"rows for {hid} disagree on the taxon"})
    applied_df = pd.DataFrame(applied, columns=CURATION_COLUMNS + ["haplotype_id"])
    not_applied_df = pd.DataFrame(not_applied, columns=CURATION_COLUMNS + ["reason_not_applied"])
    return applied_df, not_applied_df.sort_values("row").reset_index(drop=True)


def apply_curation(taxonomy: pd.DataFrame, applied: pd.DataFrame) -> pd.DataFrame:
    """
    Copy of the taxonomy table (haplotype_id, assigned_sp, assignment_level,
    assignment_notes) with the applied rows' taxa and assignment_source
    ('labels' or 'curated').
    """
    tax = taxonomy.copy()
    tax["assignment_source"] = "labels"
    if "assignment_notes" not in tax.columns:
        tax["assignment_notes"] = ""
    tax["assignment_notes"] = tax["assignment_notes"].astype(object)
    for r in applied.to_dict("records"):
        idx = tax.index[tax["haplotype_id"] == r["haplotype_id"]]
        who = ", ".join(x for x in (r.get("curator", ""), r.get("date", "")) if x)
        note = f"curated: {r['reason']}" + (f" ({who})" if who else "")
        previous = tax.loc[idx, "assignment_notes"].fillna("").astype(str)
        tax.loc[idx, "assigned_sp"] = r["taxon"]
        tax.loc[idx, "assignment_level"] = ("species" if label_level(r["taxon"]) == "species"
                                            else "genus")
        tax.loc[idx, "assignment_source"] = "curated"
        tax.loc[idx, "assignment_notes"] = [note + (f"; labels: {p}" if p else "") for p in previous]
    return tax


def curated_taxa(taxonomy: pd.DataFrame) -> Dict[str, str]:
    """haplotype_id -> assigned_sp for haplotypes whose assignment_source is curated."""
    if taxonomy is None or "assignment_source" not in taxonomy.columns:
        return {}
    cur = taxonomy[taxonomy["assignment_source"] == "curated"]
    return dict(zip(cur["haplotype_id"], cur["assigned_sp"]))
