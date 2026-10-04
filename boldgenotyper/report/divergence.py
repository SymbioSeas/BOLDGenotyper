"""Divergence (--extended): the label-based gap and the closest species pairs."""

from typing import Optional

import pandas as pd

from .data import RunData
from .page import Section, counts, num, table

GAP = "extended/divergence/label_based_gap.csv"
PAIRS = "extended/divergence/species_pair_divergence.csv"


def applies(run: RunData) -> bool:
    return run.file(GAP).exists() or run.file(PAIRS).exists()


def _csv(run: RunData, rel: str) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(run.file(rel)) if run.file(rel).exists() else None
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return None


def build(run: RunData) -> Section:
    html = ('<p>Uncorrected p-distances between haplotypes; no taxonomic threshold is applied. The '
            "label-based gap uses record labels, so mislabeled records shift it.</p>\n")
    gap = _csv(run, GAP)
    if gap is not None and not gap.empty:
        has_gap = gap["label_based_gap"].astype(str).str.lower() == "true"
        html += ("<h3>Label-based gap</h3>\n"
                 + counts([("labels with a gap", num(int(has_gap.sum()))),
                           ("labels without", num(int((~has_gap).sum())))])
                 + '<p class="note">Labels without a gap first, closest other label first.</p>\n'
                 + table(gap.assign(_g=has_gap).sort_values(["_g", "min_between"]),
                         ["label", "n_records", "n_haplotypes", "max_within", "min_between",
                          "nearest_other_label", "label_based_gap"], file=GAP))
    pairs = _csv(run, PAIRS)
    if pairs is not None and not pairs.empty:
        between = pairs[pairs.get("comparison_type", "between_species") == "between_species"]
        html += ("<h3>Closest species pairs</h3>\n"
                 + table(between.sort_values("min_divergence"),
                         ["species_pair", "min_divergence", "mean_divergence", "max_divergence"], file=PAIRS))
    if run.file("extended/divergence/divergence_heatmap.pdf").exists():
        html += '<p><a href="extended/divergence/divergence_heatmap.pdf">Divergence heatmap (PDF)</a></p>\n'
    return Section("divergence", "Divergence", html)
