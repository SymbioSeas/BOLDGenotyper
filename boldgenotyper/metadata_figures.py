"""Metadata figures (extended/metadata/figures/), by taxon, haplotype and region::

    Per field (cross_tabs):
      <field>                 records per haplotype or taxon, stacked by value
      <field>_by_haplotype    one page per taxon with two or more haplotypes
      <field>_by_region       records per region, stacked by value
    Collection dates (temporal):
      records_by_year         heatmap of records per taxon (or haplotype) and year;
                              pages: records per year stacked by haplotype
      sampling_timeline       first to last collection date per taxon (or haplotype);
                              pages: per haplotype
      records_by_year_by_region  heatmap of records per region and year
    Always:
      metadata_coverage       share of records with a value per field and taxon (or
                              haplotype); pages: per haplotype

"By region" views use the run's region column (the geography chosen for the
run, defined for any shapefile) when it has two or more values; otherwise the
country/ocean field, and the caption says why. Values and records counted are
those of records with a haplotype ("assigned records").
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from . import figure_pages as fp
from .run_figures import _group_labels, _overview_labels, _settings, figure_records, haplotype_table

logger = logging.getLogger(__name__)

COUNTRY = "country/ocean"


def _values(series: pd.Series) -> pd.Series:
    """Text values with empty and missing values as NaN."""
    s = series.astype(object).where(series.notna(), "").astype(str).str.strip()
    return s.where(s != "")


def region_facet(records: pd.DataFrame, region_col: Optional[str]) -> Tuple[str, str]:
    """The column for 'by region' views and the caption with its denominator."""
    rec = figure_records(records)
    n = len(rec)
    if region_col and region_col in rec.columns:
        vals = _values(rec[region_col]).dropna()
        if vals.nunique() >= 2:
            return region_col, (f"{len(vals):,} of {n:,} assigned records have "
                                f"{fp.article(region_col)} {region_col} value.")
        reason = f"fewer than two {region_col} values"
    else:
        reason = "no region analysis in this run"
    vals = _values(rec[COUNTRY]).dropna() if COUNTRY in rec.columns else pd.Series(dtype=object)
    return COUNTRY, (f"By country/ocean ({reason}): {len(vals):,} of {n:,} assigned records "
                     f"have a country/ocean value.")


def _top_rows(counts: pd.Series, n: int) -> List[str]:
    return list(counts.sort_values(ascending=False).index[:n])


def _field_entries(rec, haps, field, facet, settings, out_dir) -> Dict[str, Dict]:
    vals = _values(rec[field])
    has = rec[vals.notna()].assign(_value=vals[vals.notna()])
    if has.empty:
        return {}
    value_counts = has["_value"].value_counts()
    value_labels = fp.cap_items(value_counts, settings.max_colored_items, "values")
    order = list(dict.fromkeys(value_labels[sorted(value_counts.index,
                                                   key=lambda v: (-value_counts[v], v))]))
    order = [v for v in order if not fp.is_other(v)] + [v for v in order if fp.is_other(v)]
    colors = fp.item_colors(order)
    has = has.assign(_vlabel=has["_value"].map(value_labels))
    name = fp.slug(field).lower()
    col, labels, _, _ = _overview_labels(rec, settings)
    what = "haplotype" if col == "haplotype_label" else "taxon"

    def table(df, row_col):
        t = pd.crosstab(df[row_col], df["_vlabel"])
        return t[[v for v in order if v in t.columns]]

    def overview():
        t = table(has, col)
        shown = _top_rows(t.sum(axis=1), settings.max_heatmap_species)
        extra = len(t) - len(shown)
        cap = (f"{len(has):,} of {len(rec):,} assigned records have a {field} value"
               + (f"; the {len(shown)} {what}s with the most are shown ({extra} more in the pages "
                  f"and tabulations)." if extra > 0 else "."))
        return fp.stacked_barh(t.loc[shown], colors, f"{field} by {what}", cap, "records", what)

    def page(group, ids):
        sub = has[has["haplotype_id"].isin(ids)]
        if len(ids) < 2 or sub.empty:
            return None
        hl = _group_labels(haps, ids, settings.max_colored_items)
        t = table(sub.assign(_hap=sub["haplotype_id"].map(hl)), "_hap")
        n_all = int(rec["haplotype_id"].isin(ids).sum())
        return fp.stacked_barh(t, colors, f"{group}: {field} by haplotype",
                               f"{len(sub):,} of {n_all:,} records of {group} have a {field} value.",
                               "records", "haplotype")

    entries = {name: fp.write_paged_figure(name, out_dir, overview, page, fp.species_groups(haps),
                                           pages_name=f"{name}_by_haplotype")}
    facet_col, facet_caption = region_facet(rec, facet)
    if facet_col in has.columns:
        regional = has[_values(has[facet_col]).notna()]
        if not regional.empty:
            t = table(regional, facet_col)
            fig = fp.stacked_barh(t, colors, f"{field} by {facet_col}", facet_caption, "records",
                                  facet_col)
            entries[f"{name}_by_region"] = {
                "overview": [p.name for p in fp.save_figure(fig, Path(out_dir) / f"{name}_by_region")]}
    return entries


def _year_entries(rec, haps, facet, settings, out_dir) -> Dict[str, Dict]:
    from .metadata_summary import parse_collection_dates

    if "collection_date_start" not in rec.columns:
        return {}
    dated, _ = parse_collection_dates(rec, date_col="collection_date_start")
    dated = dated[dated["collection_year"].notna()].copy()
    if dated.empty:
        return {}
    dated["collection_year"] = dated["collection_year"].astype(int)
    years = list(range(dated["collection_year"].min(), dated["collection_year"].max() + 1))
    col, _, _, _ = _overview_labels(rec, settings)
    what = "haplotype" if col == "haplotype_label" else "taxon"
    colors = fp.haplotype_colors(haps, settings.max_colored_items)
    groups = fp.species_groups(haps)
    entries = {}

    def year_overview():
        m = pd.crosstab(dated[col], dated["collection_year"]).reindex(columns=years, fill_value=0)
        shown = _top_rows(m.sum(axis=1), settings.max_heatmap_species)
        extra = len(m) - len(shown)
        cap = (f"{len(dated):,} of {len(rec):,} assigned records have a collection date"
               + (f"; the {len(shown)} {what}s with the most records are shown." if extra > 0 else "."))
        return fp.heatmap(m.loc[shown], f"Records by collection year and {what}", cap, "records")

    def year_page(group, ids):
        sub = dated[dated["haplotype_id"].isin(ids)]
        if len(ids) < 2 or sub.empty:
            return None
        hl = _group_labels(haps, ids, settings.max_colored_items)
        t = pd.crosstab(sub["collection_year"], sub["haplotype_id"].map(hl)).reindex(
            index=range(sub["collection_year"].min(), sub["collection_year"].max() + 1), fill_value=0)
        order = list(dict.fromkeys(hl[h] for h in ids))
        t = t[[c for c in order if c in t.columns]]
        cols = {hl[h]: colors[h] for h in ids}
        n_all = int(rec["haplotype_id"].isin(ids).sum())
        fig = fp.stacked_barh(t.iloc[::-1], cols, f"{group}: records by collection year",
                              f"{len(sub):,} of {n_all:,} records of {group} have a collection date.",
                              "records", "year", sort_rows=False)
        return fig

    entries["records_by_year"] = fp.write_paged_figure("records_by_year", out_dir, year_overview,
                                                       year_page, groups)

    def timeline(df, row_col, row_colors, title, caption):
        import matplotlib.pyplot as plt

        spans = df.groupby(row_col)["parsed_date"].agg(["min", "max", "size"]).sort_values("min")
        fig, ax = plt.subplots(figsize=(10, max(2.5, 0.22 * len(spans) + 1.6)))
        for i, (row, s) in enumerate(spans.iterrows()):
            color = row_colors.get(row, "#4c72b0")
            ax.plot([s["min"], s["max"]], [i, i], color=color, linewidth=3, solid_capstyle="round")
            ax.plot([s["min"], s["max"]], [i, i], "o", color=color, markersize=3)
        ax.set_yticks(range(len(spans)))
        ax.set_yticklabels([f"{r} ({int(n)})" for r, n in zip(spans.index, spans["size"])],
                           fontsize=7 if len(spans) > 20 else 9)
        ax.set_xlabel("collection date")
        ax.set_title(title, fontsize=11, pad=16)
        fp.add_caption(ax, caption)
        return fig

    def timeline_overview():
        return timeline(dated, col, {}, f"First to last collection date by {what}",
                        f"{len(dated):,} of {len(rec):,} assigned records have a collection date; "
                        f"records in parentheses. Dates show sampling, not when a lineage arose.")

    def timeline_page(group, ids):
        sub = dated[dated["haplotype_id"].isin(ids)]
        if len(ids) < 2 or sub.empty:
            return None
        hl = _group_labels(haps, ids, settings.max_colored_items)
        sub = sub.assign(_hap=sub["haplotype_id"].map(hl))
        cols = {hl[h]: colors[h] for h in ids}
        return timeline(sub, "_hap", cols, f"{group}: first to last collection date by haplotype",
                        f"{len(sub):,} records of {group} with a collection date; records in parentheses.")

    entries["sampling_timeline"] = fp.write_paged_figure("sampling_timeline", out_dir,
                                                         timeline_overview, timeline_page, groups)

    facet_col, facet_caption = region_facet(rec, facet)
    if facet_col in dated.columns:
        regional = dated[_values(dated[facet_col]).notna()]
        if not regional.empty:
            m = pd.crosstab(regional[facet_col], regional["collection_year"]).reindex(
                columns=years, fill_value=0)
            fig = fp.heatmap(m, f"Records by collection year and {facet_col}", facet_caption, "records")
            entries["records_by_year_by_region"] = {
                "overview": [p.name for p in fp.save_figure(fig, Path(out_dir) / "records_by_year_by_region")]}
    return entries


def _coverage_entry(rec, haps, fields, settings, out_dir) -> Dict:
    present = [f for f in fields if f in rec.columns]
    if not present:
        return {}
    filled = pd.DataFrame({f: _values(rec[f]).notna() for f in present}, index=rec.index)
    col, _, _, _ = _overview_labels(rec, settings)
    what = "haplotype" if col == "haplotype_label" else "taxon"

    def share(rows: pd.Series) -> pd.DataFrame:
        return (filled.groupby(rows).mean() * 100).round(0)

    def overview():
        m = share(rec[col])
        shown = _top_rows(rec[col].value_counts(), settings.max_heatmap_species)
        cap = (f"Share of each {what}'s records with a value"
               + (f"; the {len(shown)} {what}s with the most records are shown."
                  if len(m) > len(shown) else "."))
        return fp.heatmap(m.loc[shown], f"Metadata coverage by {what}", cap, "% of records", vmax=100)

    def page(group, ids):
        sub = rec["haplotype_id"].isin(ids)
        if len(ids) < 2 or not sub.any():
            return None
        hl = _group_labels(haps, ids, settings.max_colored_items)
        m = share(rec.loc[sub, "haplotype_id"].map(hl))
        return fp.heatmap(m, f"{group}: metadata coverage by haplotype",
                          f"Share of each haplotype's records with a value ({int(sub.sum()):,} records).",
                          "% of records", vmax=100)

    return fp.write_paged_figure("metadata_coverage", out_dir, overview, page, fp.species_groups(haps))


def write_metadata_figures(records: pd.DataFrame, fields: List[str], region_col: Optional[str],
                           out_dir, cross_tabs: bool = True, temporal: bool = True,
                           settings=None) -> Dict[str, Dict]:
    """Metadata figures of a run; returns manifest entries keyed by figure name."""
    settings = _settings(settings)
    rec = figure_records(records)
    haps = haplotype_table(records)
    out_dir = Path(out_dir)
    entries: Dict[str, Dict] = {}
    if cross_tabs:
        for field in fields:
            if field in rec.columns:
                entries.update(_field_entries(rec, haps, field, region_col, settings, out_dir))
    if temporal:
        entries.update(_year_entries(rec, haps, region_col, settings, out_dir))
    entries["metadata_coverage"] = _coverage_entry(rec, haps, fields, settings, out_dir)
    return {k: v for k, v in entries.items() if v}
