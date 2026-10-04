"""Key findings: a fixed list of observations about a run, each shown only
when it applies, with counts and a link. None interprets biology."""

import re
from collections import Counter
from dataclasses import dataclass
from typing import List

import pandas as pd

from .data import RunData


@dataclass
class Finding:
    key: str
    text: str
    count: int
    link: str = ""


def _int(value) -> int:
    try:
        return int(value) if value is not None and not pd.isna(value) else 0
    except (TypeError, ValueError):
        return 0


def _pct(n: int, total: int) -> str:
    return f"{100 * n / total:.1f}%" if total else "0%"


def _text(series: pd.Series) -> pd.Series:
    s = series.astype(object).where(series.notna(), "").astype(str).str.strip()
    return s.where(s != "")


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


def key_findings(run: RunData) -> List[Finding]:
    s, haps, rec = run.summary, run.haplotypes, run.records
    out: List[Finding] = []

    n_qc = _int(s.get("n_qc_pass"))
    assigned = _int(s.get("n_assigned")) + _int(s.get("n_assigned_corrected"))
    if n_qc:
        out.append(Finding("assigned", f"{assigned:,} of {n_qc:,} QC-passed records ({_pct(assigned, n_qc)}) "
                                       f"are assigned to {_int(s.get('n_haplotypes')):,} haplotypes.",
                           assigned, "#section-records"))
    exclusions = {"excluded_core_coverage": _int(s.get("n_excluded_core")),
                  "excluded_error_singleton": _int(s.get("n_excluded_error")),
                  "failed_qc": _int(s.get("n_failed_qc")), "no_sequence": _int(s.get("n_no_sequence"))}
    status, n = max(exclusions.items(), key=lambda kv: kv[1])
    if n:
        out.append(Finding("largest_exclusion", f"The most common reason a record is not assigned is "
                                                f"{status}: {n:,} records.", n, "#section-records"))
    corrected = _int(s.get("n_assigned_corrected"))
    tied, far = _int(s.get("n_excluded_error_tied")), _int(s.get("n_excluded_error_far"))
    if corrected or tied or far:
        out.append(Finding("error_singletons",
                           f"{corrected:,} records of error-filtered singletons were assigned to their nearest "
                           f"haplotype; {tied:,} stay unassigned because two haplotypes were equally near and "
                           f"{far:,} because none was within the threshold.", corrected, "#section-records"))

    if not haps.empty:
        not_species = haps[haps.get("assignment_level", pd.Series(dtype=object)).astype(str) != "species"]
        if len(not_species):
            with_suggestion = (int(_text(not_species["suggested_sp"]).notna().sum())
                               if "suggested_sp" in haps.columns else 0)
            out.append(Finding("not_species_level",
                               f"{len(not_species):,} haplotypes are not assigned at species level; "
                               f"{with_suggestion:,} have a suggested taxon (never applied).",
                               len(not_species), "#section-taxonomy"))
        if "n_label_conflicts" in haps.columns:
            conflicts = pd.to_numeric(haps["n_label_conflicts"], errors="coerce").fillna(0)
            n_haps = int((conflicts > 0).sum())
            if n_haps:
                out.append(Finding("label_conflicts",
                                   f"{n_haps:,} haplotypes have members whose label conflicts with the "
                                   f"haplotype's taxon ({int(conflicts.sum()):,} records).",
                                   n_haps, "#section-taxonomy"))
        if "n_species_labels" in haps.columns:
            n_multi = int((pd.to_numeric(haps["n_species_labels"], errors="coerce").fillna(0) > 1).sum())
            if n_multi:
                out.append(Finding("multi_species_labels",
                                   f"{n_multi:,} haplotypes carry more than one species label.",
                                   n_multi, "#section-taxonomy"))
        if "assignment_source" in haps.columns:
            n_cur = int((haps["assignment_source"] == "curated").sum())
            if n_cur:
                out.append(Finding("curated", f"{n_cur:,} haplotypes have a curated taxon.", n_cur,
                                   "#section-taxonomy"))
    n_not_applied = _int((run.applied.get("taxon_curation") or {}).get("n_not_applied"))
    if n_not_applied:
        out.append(Finding("curation_not_applied", f"{n_not_applied:,} curation rows were not applied.",
                           n_not_applied, "#section-taxonomy"))
    if not haps.empty and "is_suspect" in haps.columns:
        suspect = haps[_truthy(haps["is_suspect"])]
        if len(suspect):
            reasons = Counter(re.sub(r"\s*\(.*?\)", "", r).strip()
                              for cell in _text(suspect.get("suspect_reasons", pd.Series(dtype=object))).dropna()
                              for r in str(cell).split(";") if r.strip())
            detail = "; ".join(f"{r} {n}" for r, n in reasons.most_common())
            out.append(Finding("suspect", f"{len(suspect):,} haplotypes are flagged as suspect"
                                          + (f" ({detail})." if detail else "."),
                               len(suspect), "#section-phylogeny"))
    not_mono = run.non_monophyletic()
    if not_mono:
        shown = ", ".join(not_mono[:5]) + (f" and {len(not_mono) - 5} more" if len(not_mono) > 5 else "")
        out.append(Finding("not_monophyletic", f"{len(not_mono):,} taxa do not form a clade in the tree: "
                                               f"{shown}.", len(not_mono), "#section-phylogeny"))

    if not rec.empty and "lat" in rec.columns:
        has_coords = rec["lat"].notna() & rec["lon"].notna()
        status = rec.get("coord_status", pd.Series("", index=rec.index)).astype(str)
        mapped = status == "ok"
        out.append(Finding("coordinates", f"{int(has_coords.sum()):,} of {len(rec):,} records have coordinates.",
                           int(has_coords.sum()), "#section-geography"))
        out.append(Finding("mapped", f"{int(mapped.sum()):,} records have geographic-quality coordinates and "
                                     f"are mapped.", int(mapped.sum()), "#section-geography"))
        not_mapped = status[has_coords & ~mapped].value_counts()
        if len(not_mapped):
            out.append(Finding("not_mapped", f"{int(not_mapped.sum()):,} records with coordinates are not "
                                             f"mapped ({'; '.join(f'{k} {v:,}' for k, v in not_mapped.items())}).",
                               int(not_mapped.sum()), "#section-geography"))
        if run.region_col:
            regions = _text(rec[run.region_col])
            outside = int((mapped & regions.isna()).sum())
            n_regions = int(regions.dropna().nunique())
            out.append(Finding("regions", f"{n_regions:,} {run.region_col} values have records; {outside:,} "
                                          f"mapped records fall outside every region.", n_regions,
                               "#section-geography"))
    if not rec.empty and "collection_date_start" in rec.columns:
        from ..metadata_summary import parse_collection_dates
        dated, _ = parse_collection_dates(rec[["collection_date_start"]].copy(), "collection_date_start")
        years = dated["collection_year"].dropna()
        if len(years):
            out.append(Finding("dates", f"Collection years {int(years.min())}–{int(years.max())}; "
                                        f"{len(years):,} of {len(rec):,} records ({_pct(len(years), len(rec))}) "
                                        f"have a collection date.", len(years), "#section-sampling"))
    return out
