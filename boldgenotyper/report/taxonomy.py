"""Taxonomy and labels: haplotype taxa and what needs attention."""

import pandas as pd

from .data import RunData
from .page import Section, counts, esc, num, table

EVIDENCE_TABLE = ["haplotype_id", "assigned_sp", "n_members", "labels", "nearest_species_haplotype",
                  "nearest_species_distance", "within_taxon_spread", "bin_taxon", "tree_sister_taxon",
                  "tree_sister_support", "suggested_sp", "suggestion_basis"]


def _conflicts_by_label(lc: pd.DataFrame) -> pd.DataFrame:
    if lc.empty:
        return lc
    keep = {c: "first" for c in ("haplotype_taxon", "majority_label", "majority_fraction",
                                 "nearest_haplotype_with_label", "distance_to_nearest_with_label")
            if c in lc.columns}
    g = lc.groupby(["haplotype_id", "label"], dropna=False)
    out = g.agg(keep) if keep else pd.DataFrame(index=g.size().index)
    out["n_records"] = g.size()
    out = out.reset_index()
    sort = ["n_records"] + (["distance_to_nearest_with_label"] if "distance_to_nearest_with_label" in out else [])
    return out.sort_values(sort, ascending=[False] + [True] * (len(sort) - 1))


def build(run: RunData) -> Section:
    haps = run.haplotypes
    if haps.empty:
        return Section("taxonomy", "Taxonomy and labels", '<p class="note">No haplotypes.</p>\n')
    level = haps.get("assignment_level", pd.Series(dtype=object)).astype(str)
    source = haps.get("assignment_source", pd.Series("labels", index=haps.index)).astype(str)
    threshold = (run.config.get("taxonomy") or {}).get("majority_species_threshold")
    html = (f"<p>A haplotype's taxon is the species carried by more than {num(threshold, '{}')} of its "
            "members' counted labels, otherwise the genus; a curation file can replace it "
            "(<code>assignment_source</code> curated). Evidence for each taxon is in "
            "<code>haplotypes.csv</code>.</p>\n"
            + counts([("species level", num(int((level == "species").sum()))),
                      ("genus level or unassigned", num(int((level != "species").sum()))),
                      ("curated", num(int((source == "curated").sum())))]))

    html += "<h3>Needs attention</h3>\n"
    not_species = haps[level != "species"].sort_values("n_members", ascending=False)
    html += ("<h4>Haplotypes not assigned at species level</h4>\n"
             '<p class="note">With the evidence and any suggested taxon (never applied; record an '
             "accepted call in a curation file).</p>\n"
             + table(not_species, EVIDENCE_TABLE, file="haplotypes.csv"))
    html += ("<h4>Label conflicts</h4>\n"
             '<p class="note">Per haplotype and conflicting label, most records first; the evidence for '
             "each record is in <code>label_conflicts.csv</code>.</p>\n"
             + table(_conflicts_by_label(run.label_conflicts),
                     ["haplotype_id", "haplotype_taxon", "label", "n_records", "majority_label",
                      "majority_fraction", "nearest_haplotype_with_label", "distance_to_nearest_with_label"],
                     file="label_conflicts.csv"))
    curated = haps[source == "curated"]
    if len(curated):
        html += "<h4>Curated haplotypes</h4>\n" + table(
            curated, ["haplotype_id", "assigned_sp", "n_members", "assignment_notes"], file="haplotypes.csv")
    not_applied = (run.applied.get("taxon_curation") or {}).get("not_applied") or []
    if not_applied:
        html += ("<h4>Curation rows not applied</h4>\n"
                 + table(pd.DataFrame(not_applied), ["row", "taxon", "reason_not_applied"],
                         file="run_provenance.json"))
    cur = run.applied.get("taxon_curation")
    if cur:
        html += (f'<p class="note">Curation file <code>{esc(cur.get("path"))}</code> (md5 {esc(cur.get("md5"))}): '
                 f'{num(cur.get("n_applied"))} rows applied, {num(cur.get("n_not_applied"))} not applied.</p>\n')
    return Section("taxonomy", "Taxonomy and labels", html)
