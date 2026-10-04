"""A run's figures, drawn from the final per-record table and the tree.

geography/                 distribution_map (overview + one page per taxon)
extended/geography/        region_abundance_relative, region_abundance_total
                           (overview + one page per taxon with two or more haplotypes)
phylogeny/                 tree figures (tree_figures)
extended/metadata/figures/ metadata views (metadata_figures)

Overviews show haplotypes for up to max_overview_haplotypes haplotypes and taxa
above; a panel colors at most max_colored_items items (max_overview_species in
taxon overviews), the rest in gray as "Other (n)". Per-taxon pages color each
haplotype the same way on every page (figure_pages.haplotype_colors).
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import figure_pages as fp

logger = logging.getLogger(__name__)


def _settings(settings):
    if settings is None:
        from .config import VisualizationConfig
        return VisualizationConfig()
    return settings


def figure_records(records: pd.DataFrame) -> pd.DataFrame:
    """Records with a haplotype; taxon is assigned_sp or 'Unassigned'."""
    rec = records[records["haplotype_id"].notna()].copy()
    rec["haplotype_id"] = rec["haplotype_id"].astype(str)
    taxon = rec["assigned_sp"].astype(object).where(rec["assigned_sp"].notna(), "")
    rec["taxon"] = taxon.astype(str).str.strip().replace("", fp.UNASSIGNED)
    if "haplotype_label" not in rec.columns:
        rec["haplotype_label"] = rec["haplotype_id"]
    rec["haplotype_label"] = rec["haplotype_label"].astype(object).where(
        rec["haplotype_label"].notna(), rec["haplotype_id"]).astype(str)
    return rec


def haplotype_table(records: pd.DataFrame) -> pd.DataFrame:
    """haplotype_id, haplotype_label, assigned_sp, n_members from the records."""
    rec = figure_records(records)
    g = rec.groupby("haplotype_id")
    return pd.DataFrame({
        "haplotype_id": list(g.groups),
        "haplotype_label": g["haplotype_label"].first().values,
        "assigned_sp": g["taxon"].first().values,
        "n_members": g.size().values,
    })


def _overview_labels(rec: pd.DataFrame, settings, subset: Optional[pd.DataFrame] = None):
    """(column, item -> display label, colors, noun) for an overview of rec."""
    level = fp.overview_level(rec["haplotype_id"].nunique(), settings.max_overview_haplotypes)
    col, cap, noun = (("haplotype_label", settings.max_colored_items, "haplotypes") if level == "haplotype"
                      else ("taxon", settings.max_overview_species, "species"))
    counts = (subset if subset is not None else rec)[col].value_counts()
    labels = fp.cap_items(counts, cap, noun)
    order = list(dict.fromkeys(labels[sorted(counts.index, key=lambda i: (-counts[i], str(i)))]))
    order = [lab for lab in order if not fp.is_other(lab)] + [lab for lab in order if fp.is_other(lab)]
    return col, labels.to_dict(), fp.item_colors(order), order


def _group_labels(haps: pd.DataFrame, ids: List[str], max_items: int) -> Dict[str, str]:
    """Haplotype ID -> legend label within a group, capped like haplotype_colors."""
    h = haps.set_index("haplotype_id")
    n = pd.to_numeric(h["n_members"], errors="coerce").fillna(0).reindex(ids).fillna(0)
    capped = fp.cap_items(n, max_items, "haplotypes")
    return {hid: (capped[hid] if fp.is_other(capped[hid]) else str(h.loc[hid, "haplotype_label"]))
            for hid in ids}


# ---------------------------------------------------------------------------
# Maps
# ---------------------------------------------------------------------------

def map_extent(lats: pd.Series, lons: pd.Series, buffer_degrees: float = 10.0):
    """(lon_min, lon_max, lat_min, lat_max) and a figure size: no thin strips."""
    from .visualization import calculate_map_extent_with_buffer

    lon_min, lon_max, lat_min, lat_max = calculate_map_extent_with_buffer(
        lats, lons, buffer_degrees=buffer_degrees)
    if lon_max - lon_min < 40:                       # a few sites: keep coastlines in view
        mid = (lon_max + lon_min) / 2
        lon_min, lon_max = max(-180.0, mid - 20), min(180.0, mid + 20)
    min_lat_span = (lon_max - lon_min) / 2.2
    if lat_max - lat_min < min_lat_span:
        mid = (lat_max + lat_min) / 2
        lat_min, lat_max = max(-90.0, mid - min_lat_span / 2), min(90.0, mid + min_lat_span / 2)
    height = min(9.0, max(4.0, 11 * (lat_max - lat_min) / max(lon_max - lon_min, 1e-6)))
    return (lon_min, lon_max, lat_min, lat_max), (11, height)


def _mapped(rec: pd.DataFrame) -> pd.DataFrame:
    ok = rec["coord_status"] == "ok" if "coord_status" in rec.columns else True
    return rec[ok & rec["lat"].notna() & rec["lon"].notna()]


def _draw_map(points: pd.DataFrame, label_col: str, colors: Dict[str, str], order: List[str],
              title: str, caption: str, detail: str):
    """points: lat, lon and a label column; one marker per site and label."""
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from .visualization import _add_cartopy_background, _add_map_gridlines, _point_sizes

    extent, size = map_extent(points["lat"], points["lon"])
    fig = plt.figure(figsize=size)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    _add_cartopy_background(ax, detail_level=detail)
    gl = _add_map_gridlines(ax)
    gl.top_labels = gl.right_labels = False
    sites = points.groupby(["lat", "lon", label_col]).size().reset_index(name="n")
    sites["size"] = _point_sizes(sites["n"], 20, 220)
    totals = points[label_col].value_counts()
    for lab in order:                       # Other last, drawn first underneath
        s = sites[sites[label_col] == lab]
        if s.empty:
            continue
        ax.scatter(s["lon"], s["lat"], s=s["size"], color=colors[lab], alpha=0.85,
                   edgecolors="black", linewidths=0.4, transform=ccrs.PlateCarree(),
                   zorder=3 if fp.is_other(lab) else 4, label=f"{lab} ({totals.get(lab, 0)})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False, fontsize=7,
              title="records in parentheses", title_fontsize=7)
    ax.set_title(title, fontsize=11, pad=16)
    fp.add_caption(ax, caption)
    return fig


def _detail(settings, overview: bool) -> str:
    d = settings.map_background_detail
    return ("full" if overview else "simple") if d == "auto" else d


def distribution_map_entry(records: pd.DataFrame, out_dir, settings=None) -> Dict:
    settings = _settings(settings)
    rec = figure_records(records)
    mapped = _mapped(rec)
    haps = haplotype_table(records)
    colors = fp.haplotype_colors(haps, settings.max_colored_items)
    caption = (f"{len(mapped):,} of {len(rec):,} records with a haplotype have geographic-quality "
               f"coordinates (country or province centroids, [0, 0] and unparseable values are not mapped).")

    def overview():
        if mapped.empty:
            return None
        col, labels, cols, order = _overview_labels(rec, settings, subset=mapped)
        pts = mapped.assign(_label=mapped[col].map(labels))
        what = "haplotype" if col == "haplotype_label" else "taxon"
        return _draw_map(pts, "_label", cols, order, f"Records by {what}", caption,
                         _detail(settings, True))

    def page(group, ids):
        pts = mapped[mapped["haplotype_id"].isin(ids)]
        if pts.empty:
            return None
        labels = _group_labels(haps, ids, settings.max_colored_items)
        pts = pts.assign(_label=pts["haplotype_id"].map(labels))
        order = list(dict.fromkeys(labels[h] for h in ids))
        cols = {labels[h]: colors[h] for h in ids}
        n_all = int(rec["haplotype_id"].isin(ids).sum())
        return _draw_map(pts, "_label", cols, order, f"{group}: records by haplotype",
                         f"{len(pts):,} of {n_all:,} records of {group} have geographic-quality coordinates.",
                         _detail(settings, False))

    return fp.write_paged_figure("distribution_map", out_dir, overview, page, fp.species_groups(haps))


# ---------------------------------------------------------------------------
# Region abundance
# ---------------------------------------------------------------------------

def _draw_region_bars(table: pd.DataFrame, colors: Dict[str, str], relative: bool, title: str,
                      caption: str, region_label: str):
    """table: regions x items (record counts), columns in legend order."""
    import matplotlib.pyplot as plt

    table = table.loc[table.sum(axis=1).sort_values().index]
    values = table.div(table.sum(axis=1), axis=0) * 100 if relative else table
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.45 * len(table) + 1.6)))
    left = np.zeros(len(table))
    for item in table.columns:
        ax.barh(table.index.astype(str), values[item].values, left=left, color=colors[item],
                edgecolor="white", linewidth=0.4, label=f"{item} ({int(table[item].sum())})")
        left += values[item].values
    ax.set_xlabel("share of the region's records (%)" if relative else "records")
    ax.set_ylabel(region_label)
    if relative:
        ax.set_xlim(0, 100)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), ncol=1 if table.shape[1] <= 16 else 2,
              frameon=False, fontsize=7, title="records in parentheses", title_fontsize=7)
    ax.set_title(title, fontsize=11, pad=16)
    fp.add_caption(ax, caption)
    return fig


def region_abundance_entry(records: pd.DataFrame, region_col: str, out_dir, relative: bool,
                           settings=None) -> Dict:
    settings = _settings(settings)
    rec = figure_records(records)
    regioned = rec[rec[region_col].notna() & (rec[region_col].astype(str).str.strip() != "")]
    haps = haplotype_table(records)
    colors = fp.haplotype_colors(haps, settings.max_colored_items)
    name = "region_abundance_relative" if relative else "region_abundance_total"
    caption = (f"{len(regioned):,} of {len(rec):,} records with a haplotype have "
               f"{fp.article(region_col)} {region_col}.")

    def overview():
        if regioned.empty:
            return None
        col, labels, cols, order = _overview_labels(rec, settings, subset=regioned)
        t = pd.crosstab(regioned[region_col], regioned[col].map(labels))
        t = t[[c for c in order if c in t.columns]]
        what = "haplotype" if col == "haplotype_label" else "taxon"
        return _draw_region_bars(t, cols, relative,
                                 f"Records per {region_col} by {what}", caption, region_col)

    def page(group, ids):
        sub = regioned[regioned["haplotype_id"].isin(ids)]
        if len(ids) < 2 or sub.empty:
            return None
        labels = _group_labels(haps, ids, settings.max_colored_items)
        order = list(dict.fromkeys(labels[h] for h in ids))
        cols = {labels[h]: colors[h] for h in ids}
        t = pd.crosstab(sub[region_col], sub["haplotype_id"].map(labels))
        t = t[[c for c in order if c in t.columns]]
        n_all = int(rec["haplotype_id"].isin(ids).sum())
        return _draw_region_bars(t, cols, relative, f"{group}: records per {region_col} by haplotype",
                                 f"{len(sub):,} of {n_all:,} records of {group} have "
                                 f"{fp.article(region_col)} {region_col}.",
                                 region_col)

    return fp.write_paged_figure(name, out_dir, overview, page, fp.species_groups(haps))


def write_geography_figures(records: pd.DataFrame, region_col: Optional[str], geography_dir,
                            extended_dir=None, settings=None) -> Dict[str, Dict]:
    """Manifest entries {'geography': {...}, 'extended/geography': {...}}."""
    out = {"geography": {}, "extended/geography": {}}
    entry = distribution_map_entry(records, geography_dir, settings)
    if entry:
        out["geography"]["distribution_map"] = entry
    if extended_dir is not None and region_col and region_col in records.columns:
        for relative in (True, False):
            name = "region_abundance_relative" if relative else "region_abundance_total"
            entry = region_abundance_entry(records, region_col, extended_dir, relative, settings)
            if entry:
                out["extended/geography"][name] = entry
    return {k: v for k, v in out.items() if v}
