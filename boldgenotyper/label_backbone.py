"""Suggest a label map by checking a run's labels against a taxonomic backbone.

For each distinct species-level label in a run's records.csv (the deposited
label, before any label map), ask GBIF's species-match service or WoRMS for
the name. A label that matches a different name only fuzzily is suggested as
a spelling fix; a label that is a synonym is suggested with its accepted name.
Exact accepted names are not listed.

The suggestions are a starting point for a label map the user reviews and
confirms (label_variants.load_label_map); nothing here changes a run. This
runs outside the pipeline, so pipeline runs stay offline and reproducible; the
file records the backbone and the date. Network access uses the standard
library only.
"""

import datetime
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict

import pandas as pd

from .label_variants import clean_label
from .taxon_labels import label_level

BACKBONES = {
    "gbif": ("GBIF", "GBIF Backbone Taxonomy"),
    "worms": ("WoRMS", "WoRMS (World Register of Marine Species)"),
}
SUGGESTION_COLUMNS = ["from", "to", "reason", "source", "kind", "backbone_match_type",
                      "backbone_status"]
PAUSE_SECONDS = 0.1
_GBIF_URL = "https://api.gbif.org/v1/species/match?name={}&strict=false"
_WORMS_URL = ("https://www.marinespecies.org/rest/AphiaRecordsByMatchNames"
              "?scientificnames[]={}&marine_only=false")


class BackboneError(RuntimeError):
    """The backbone could not be reached or returned something unreadable."""


def _get_json(url: str, timeout: float):
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            body = response.read()
    except (urllib.error.URLError, OSError) as e:
        raise BackboneError(f"could not reach {url.split('?')[0]}: {e}")
    if not body:
        return None
    try:
        return json.loads(body)
    except ValueError as e:
        raise BackboneError(f"unreadable answer from {url.split('?')[0]}: {e}")


def match_name(name: str, backbone: str = "gbif", timeout: float = 20) -> Dict[str, str]:
    """Match one name: name, backbone, match_type, status, matched_name, accepted_name.

    match_type and status are as the backbone reports them (GBIF: EXACT,
    FUZZY, HIGHERRANK, NONE; ACCEPTED, SYNONYM, ... / WoRMS: exact, near_1,
    phonetic, ...; accepted, unaccepted). Empty strings when nothing matched.
    """
    quoted = urllib.parse.quote(name)
    result = {"name": name, "backbone": backbone, "match_type": "", "status": "",
              "matched_name": "", "accepted_name": ""}
    if backbone == "gbif":
        data = _get_json(_GBIF_URL.format(quoted), timeout) or {}
        status = str(data.get("status", ""))
        matched = str(data.get("canonicalName", ""))
        accepted = str(data.get("species", "")) if "SYNONYM" in status.upper() else matched
        result.update(match_type=str(data.get("matchType", "")), status=status,
                      matched_name=matched, accepted_name=accepted)
    elif backbone == "worms":
        data = _get_json(_WORMS_URL.format(quoted), timeout) or []
        hits = data[0] if data and isinstance(data[0], list) else data
        if hits:
            hit = hits[0]
            result.update(match_type=str(hit.get("match_type", "")),
                          status=str(hit.get("status", "")),
                          matched_name=str(hit.get("scientificname", "")),
                          accepted_name=str(hit.get("valid_name", "") or hit.get("scientificname", "")))
    else:
        raise ValueError(f"unknown backbone {backbone!r} (choose from {', '.join(BACKBONES)})")
    return result


def _is_fuzzy(match_type: str, backbone: str) -> bool:
    mt = match_type.lower()
    return mt == "fuzzy" if backbone == "gbif" else (mt != "" and not mt.startswith("exact"))


def _is_synonym(status: str, backbone: str) -> bool:
    st = status.lower()
    return "synonym" in st if backbone == "gbif" else (st != "" and st != "accepted")


def suggest_label_map(run_dir: Path, backbone: str = "gbif") -> pd.DataFrame:
    """Suggested label-map rows for a run folder's species-level labels."""
    records_csv = Path(run_dir) / "records.csv"
    if not records_csv.exists():
        raise FileNotFoundError(f"{records_csv} not found (a v2 run folder is needed)")
    rec = pd.read_csv(records_csv, dtype=str, keep_default_na=False,
                      usecols=lambda c: c in ("label", "label_mapped_from"))
    labels = rec["label"]
    if "label_mapped_from" in rec.columns:
        labels = rec["label_mapped_from"].where(rec["label_mapped_from"] != "", labels)
    names = sorted({clean_label(x) for x in labels if label_level(clean_label(x)) == "species"})

    short, long_name = BACKBONES[backbone]
    today = datetime.date.today().isoformat()
    rows = []
    for i, name in enumerate(names):
        if i and PAUSE_SECONDS:
            time.sleep(PAUSE_SECONDS)
        m = match_name(name, backbone)
        target = clean_label(m["accepted_name"])
        if not target or target == name:
            continue
        if _is_synonym(m["status"], backbone):
            kind, how = "synonym", "synonym of the accepted name"
        elif _is_fuzzy(m["match_type"], backbone):
            kind, how = "spelling", "fuzzy match"
        else:
            continue
        rows.append({"from": name, "to": target, "reason": f"{short} {how}, {today}",
                     "source": f"{long_name}, accessed {today}", "kind": kind,
                     "backbone_match_type": m["match_type"], "backbone_status": m["status"]})
    return pd.DataFrame(rows, columns=SUGGESTION_COLUMNS)


def write_suggestions(suggestions: pd.DataFrame, path: Path) -> Path:
    """Write the suggestions, replacing path only once the file is complete."""
    path = Path(path)
    tmp = path.with_name(path.name + ".partial")
    suggestions.to_csv(tmp, index=False)
    tmp.replace(path)
    return path
