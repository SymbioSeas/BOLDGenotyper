"""Place queries in the reference run: metadata, regions and figures.

The reference is a v2 run folder (run_provenance.json, records.csv,
haplotypes.csv, and phylogeny/tree.nwk when a tree was built). Query metadata
(--query-metadata, CSV or TSV) gives each query's coordinates, collection date,
country and region; a query's region is its `region` value, else its
coordinates placed with the reference run's own shapefile and field.
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from . import metadata, metadata_summary

logger = logging.getLogger(__name__)

QUERY_METADATA_FIELDS = ["query_id", "lat", "lon", "coord", "collection_date", "country/ocean", "region"]
QUERY_REGION_COLUMNS = ["equivalence_set", "species_call", "region", "n_reference_records",
                        "n_queries", "query_ids"]


class QueryMetadataError(ValueError):
    """The query metadata file cannot be used."""


@dataclass
class ReferenceRun:
    run_dir: Path
    records: pd.DataFrame
    haplotypes: pd.DataFrame
    provenance: dict
    region_col: Optional[str]
    tree_path: Optional[Path]


def load_reference(run_dir) -> Optional[ReferenceRun]:
    """The reference run's tables; None unless run_dir is a v2 run folder."""
    if run_dir is None:
        return None
    run_dir = Path(run_dir)
    names = ("run_provenance.json", "records.csv", "haplotypes.csv")
    if not run_dir.is_dir() or not all((run_dir / n).exists() for n in names):
        return None
    provenance = json.loads((run_dir / "run_provenance.json").read_text())
    records = pd.read_csv(run_dir / "records.csv", low_memory=False)
    records, _ = metadata_summary.parse_collection_dates(records, date_col="collection_date_start")
    haplotypes = pd.read_csv(run_dir / "haplotypes.csv")
    applied = provenance.get("applied", {})
    col = applied.get("geo_category")
    region_col = col if applied.get("geographic_analysis") and col in records.columns else None
    tree = run_dir / "phylogeny" / "tree.nwk"
    return ReferenceRun(run_dir, records, haplotypes, provenance, region_col,
                        tree if tree.exists() else None)


def load_query_metadata(path, query_ids: List[str]) -> pd.DataFrame:
    """
    Read --query-metadata: one row per query in query_ids order.

    CSV, or TSV for .tsv/.txt. query_id is required; lat/lon or coord,
    collection_date, country/ocean and region are optional. Rows for IDs that
    are not queries are dropped and listed in attrs["unmatched_ids"]
    (attrs["n_rows"]: rows in the file). Adds
    lat, lon, coord_status, is_geographic_quality and collection_year.
    """
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    except (OSError, UnicodeDecodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        raise QueryMetadataError(f"{path}: cannot be read ({e})")
    df.columns = [str(c).strip().lower() for c in df.columns]
    if "query_id" not in df.columns:
        raise QueryMetadataError(f"{path}: no query_id column")
    df["query_id"] = df["query_id"].str.strip()
    repeated = sorted(df.loc[df["query_id"].duplicated(), "query_id"].unique())
    if repeated:
        raise QueryMetadataError(f"{path}: query_id repeated: {', '.join(repeated)}")
    n_rows = len(df)
    wanted = set(query_ids)
    unmatched = [q for q in df["query_id"] if q not in wanted]
    if unmatched:
        logger.warning(f"Query metadata rows for IDs that are not queries are ignored: "
                       f"{', '.join(unmatched)}")
    df = df[[c for c in QUERY_METADATA_FIELDS if c in df.columns]]
    out = pd.DataFrame({"query_id": list(query_ids)}).merge(df, on="query_id", how="left")
    for col in QUERY_METADATA_FIELDS:
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].fillna("").astype(str).str.strip()
    out = metadata.parse_coord_column(out, coord_col="coord")
    out = metadata.mark_coordinate_quality(out)
    out, _ = metadata_summary.parse_collection_dates(out, date_col="collection_date")
    out.attrs["unmatched_ids"] = unmatched
    out.attrs["n_rows"] = n_rows
    return out


def _regions_from_coordinates(points: pd.DataFrame, reference: ReferenceRun) -> pd.Series:
    """Region per query_id from the reference run's shapefile; empty when unavailable."""
    from . import geographic

    applied = reference.provenance.get("applied", {})
    custom = applied.get("custom_shapefile")
    shapefile = custom or (reference.provenance.get("config", {}).get("geographic", {})
                           .get("goas_shapefile_path"))
    if not shapefile or not Path(shapefile).exists():
        logger.warning(f"Query regions from coordinates skipped: the reference run's shapefile "
                       f"is no longer at {shapefile}")
        return pd.Series(dtype=object)
    if custom:
        assigned = geographic.assign_regions_from_shapefile(
            points, shapefile_path=shapefile, shapefile_field=applied.get("shapefile_field"),
            output_column="region_assigned")
    else:
        assigned = geographic.assign_ocean_basins(
            points, shapefile_path=shapefile, basin_col="region_assigned",
            fallback_to_country_ocean=False)
    regions = assigned.set_index("query_id")["region_assigned"]
    return regions.where(regions.notna() & (regions != "Unknown"), "").astype(str)


def assign_query_regions(meta: pd.DataFrame, reference: Optional[ReferenceRun]) -> pd.DataFrame:
    """
    Add region and region_source (metadata, coordinates or empty).

    A region value wins; otherwise geographic-quality coordinates are placed
    with the reference run's shapefile (custom, else GOaS) when the reference
    ran geographic analysis and the shapefile is still at its path.
    """
    out = meta.copy()
    out["region"] = out["region"].fillna("").astype(str).str.strip()
    out["region_source"] = np.where(out["region"] != "", "metadata", "")
    need = (out["region"] == "") & out["is_geographic_quality"].astype(bool)
    applied = (reference.provenance.get("applied", {}) if reference is not None else {})
    if need.any() and applied.get("geographic_analysis"):
        try:
            regions = _regions_from_coordinates(out.loc[need, ["query_id", "lat", "lon"]], reference)
        except Exception as e:
            logger.warning(f"Query regions from coordinates failed: {e}")
            regions = pd.Series(dtype=object)
        found = out["query_id"].map(regions).fillna("")
        placed = need & (found != "")
        out.loc[placed, "region"] = found[placed]
        out.loc[placed, "region_source"] = "coordinates"
    return out


def query_regions_table(calls: pd.DataFrame, meta: Optional[pd.DataFrame],
                        reference: Optional[ReferenceRun]) -> pd.DataFrame:
    """Reference records and queries per equivalence set and region (QUERY_REGION_COLUMNS)."""
    if reference is None or not reference.region_col:
        return pd.DataFrame(columns=QUERY_REGION_COLUMNS)
    rc = reference.region_col
    rec = reference.records
    rec = rec[rec["haplotype_id"].notna() & rec[rc].notna() & (rec[rc].astype(str).str.strip() != "")]
    query_region = (dict(zip(meta["query_id"], meta["region"])) if meta is not None
                    and "region" in meta.columns else {})
    rows = []
    placed = calls[calls["call_level"] != "none"]
    for eq_set, group in placed.groupby("equivalence_set", sort=False):
        counts = rec[rec["haplotype_id"].isin(eq_set.split(";"))][rc].value_counts()
        by_region = {}
        for q in group["query_id"]:
            if query_region.get(q):
                by_region.setdefault(query_region[q], []).append(q)
        for region in sorted(set(counts.index) | set(by_region)):
            rows.append({
                "equivalence_set": eq_set,
                "species_call": group["species_call"].iloc[0],
                "region": region,
                "n_reference_records": int(counts.get(region, 0)),
                "n_queries": len(by_region.get(region, [])),
                "query_ids": ";".join(by_region.get(region, [])),
            })
    return pd.DataFrame(rows, columns=QUERY_REGION_COLUMNS)


def write_query_regions(table: pd.DataFrame, figures_dir) -> Optional[Path]:
    """query_figures/query_regions.csv; None when the table is empty."""
    if table.empty:
        return None
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)
    path = figures_dir / "query_regions.csv"
    table.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# Figures: a multi-page PDF (query_figures/<name>.pdf) and one SVG per page
# (query_figures/<name>/<page>.svg), text kept as text
# ---------------------------------------------------------------------------

ACCENT = "#c0392b"
REFERENCE_GRAY = "#b8b8b8"
MAX_COLORED = 12


def _slug(text: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_") or "page"


def _save_pages(name: str, pages, out_dir) -> List[Path]:
    """Write pages [(page title, figure), ...]; returns [pdf, svg, ...]."""
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    if not pages:
        return []
    out_dir = Path(out_dir)
    svg_dir = out_dir / name
    svg_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{name}.pdf"
    paths, used = [pdf_path], set()
    with matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
        with PdfPages(pdf_path) as pdf:
            for title, fig in pages:
                slug, n = _slug(title), 2
                while slug in used:
                    slug, n = f"{_slug(title)}_{n}", n + 1
                used.add(slug)
                pdf.savefig(fig, bbox_inches="tight")
                svg = svg_dir / f"{slug}.svg"
                fig.savefig(svg, bbox_inches="tight")
                paths.append(svg)
                plt.close(fig)
    return paths


def _placed(calls: pd.DataFrame) -> pd.DataFrame:
    return calls[calls["call_level"] != "none"]


def _set_haplotypes(eq_sets) -> List[str]:
    return sorted({h for eq in eq_sets for h in str(eq).split(";") if h})


def _caption(fig, lines: List[str]) -> None:
    text = "\n".join(line for line in lines if line)
    if text:
        fig.text(0.01, -0.02, text, ha="left", va="top", fontsize=8, wrap=True)


def plot_query_tree(calls: pd.DataFrame, reference: Optional[ReferenceRun], out_dir) -> List[Path]:
    """
    One page per genus of the called sets' haplotypes: the reference tree
    pruned to that genus, each query written beside its set's tips.
    """
    import copy
    import matplotlib.pyplot as plt
    from Bio import Phylo

    if reference is None or reference.tree_path is None:
        logger.info("Query tree skipped: the reference run has no tree")
        return []
    tree = Phylo.read(str(reference.tree_path), "newick")
    haps = reference.haplotypes.set_index("haplotype_id")

    def genus(hid):
        sp = haps["assigned_sp"].get(hid) if "assigned_sp" in haps.columns else None
        return str(sp).split()[0] if isinstance(sp, str) and sp.strip() else ""

    def label(hid):
        lab = haps["haplotype_label"].get(hid) if "haplotype_label" in haps.columns else None
        return lab if isinstance(lab, str) and lab else hid

    queries_by_hap = {}
    for q, eq in zip(_placed(calls)["query_id"], _placed(calls)["equivalence_set"]):
        for h in str(eq).split(";"):
            queries_by_hap.setdefault(h, []).append(q)
    not_placed = calls.loc[calls["call_level"] == "none", "query_id"].tolist()
    tips = {t.name for t in tree.get_terminals()}
    pages = []
    for g in sorted({genus(h) for h in queries_by_hap if genus(h)}):
        keep = {t for t in tips if genus(t) == g}
        if not keep:
            continue
        sub = copy.deepcopy(tree)
        for t in list(sub.get_terminals()):
            if t.name not in keep:
                sub.prune(t)
        n = sub.count_terminals()
        fig, ax = plt.subplots(figsize=(10, max(3.0, 0.28 * n + 1.5)))

        def label_func(clade):
            if not clade.is_terminal():
                return None
            text = label(clade.name)
            if clade.name in queries_by_hap:
                text += "   <- " + ", ".join(queries_by_hap[clade.name])
            return text

        Phylo.draw(sub, axes=ax, do_show=False, label_func=label_func,
                   label_colors=lambda text: ACCENT if "<-" in text else "black",
                   branch_labels=lambda c: (f"{c.confidence:.2f}" if c.confidence is not None
                                            and c.confidence >= 0.7 else None))
        ax.set_title(f"{g}: queries beside the haplotypes of their equivalence sets", fontsize=11)
        ax.set_ylabel("")
        ax.set_yticks([])
        missing = sorted(h for h in queries_by_hap if genus(h) == g and h not in tips)
        _caption(fig, [
            "Reference tree pruned to this genus; branch values are supports of 0.70 or more.",
            f"Not in the tree (suspect or below --tree-min-members): {', '.join(missing)}" if missing else "",
            f"Not placed: {', '.join(not_placed)}" if not_placed else "",
        ])
        pages.append((g, fig))
    return _save_pages("query_tree", pages, out_dir)


def _haplotype_colors(counts: pd.Series) -> dict:
    """Colors for the most common haplotypes; the rest share REFERENCE_GRAY's darker 'Other'."""
    from .visualization import get_genotype_colors

    top = list(counts.sort_values(ascending=False).index[:MAX_COLORED])
    if len(counts) > MAX_COLORED:
        top = top[:MAX_COLORED - 1]
    colors = dict(zip(top, get_genotype_colors(len(top))))
    return {h: colors.get(h, "#6e6e6e") for h in counts.index}


def plot_query_map(calls: pd.DataFrame, meta: Optional[pd.DataFrame],
                   reference: Optional[ReferenceRun], out_dir) -> List[Path]:
    """
    One page per called taxon: the reference run's geographic-quality records
    of the taxon (gray), of the queries' equivalence sets (colored by
    haplotype), and the queries with coordinates (labeled stars).
    """
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from .run_figures import map_extent
    from .visualization import _add_cartopy_background, _add_map_gridlines

    if reference is None:
        return []
    rec = reference.records
    if "coord_status" in rec.columns:
        rec = rec[rec["coord_status"] == "ok"]
    rec = rec[rec["lat"].notna() & rec["lon"].notna()] if {"lat", "lon"} <= set(rec.columns) else rec.iloc[0:0]
    qmeta = (meta.set_index("query_id") if meta is not None else None)
    placed = _placed(calls).copy()
    placed["page"] = [sp if sp else f"{eq} (no species call)"
                      for sp, eq in zip(placed["species_call"].fillna(""), placed["equivalence_set"])]
    pages = []
    for page, group in placed.groupby("page", sort=True):
        set_haps = _set_haplotypes(group["equivalence_set"])
        set_recs = rec[rec["haplotype_id"].isin(set_haps)]
        species = group["species_call"].iloc[0]
        other = (rec[(rec["assigned_sp"] == species) & ~rec["haplotype_id"].isin(set_haps)]
                 if species else rec.iloc[0:0])
        q_ok, q_missing = [], []
        for q in group["query_id"]:
            if qmeta is not None and q in qmeta.index and bool(qmeta.loc[q, "is_geographic_quality"]):
                q_ok.append((q, qmeta.loc[q, "lat"], qmeta.loc[q, "lon"]))
            else:
                q_missing.append(q)
        lats = pd.concat([set_recs["lat"], other["lat"], pd.Series([x[1] for x in q_ok], dtype=float)])
        lons = pd.concat([set_recs["lon"], other["lon"], pd.Series([x[2] for x in q_ok], dtype=float)])
        if lats.empty:
            continue
        extent, size = map_extent(lats, lons)
        fig = plt.figure(figsize=size)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        _add_cartopy_background(ax, detail_level="simple")
        _add_map_gridlines(ax)
        if len(other):
            ax.scatter(other["lon"], other["lat"], s=14, c=REFERENCE_GRAY, alpha=0.7,
                       transform=ccrs.PlateCarree(), label=f"other {species} haplotypes", zorder=3)
        colors = _haplotype_colors(set_recs["haplotype_id"].value_counts())
        for hid, g in set_recs.groupby("haplotype_id"):
            sites = g.groupby(["lat", "lon"]).size().reset_index(name="n")
            ax.scatter(sites["lon"], sites["lat"], s=20 + 12 * np.sqrt(sites["n"]),
                       c=colors.get(hid, "#6e6e6e"), edgecolors="white", linewidths=0.5,
                       transform=ccrs.PlateCarree(), label=f"{hid} ({len(g)} records)", zorder=4)
        for q, lat, lon in q_ok:
            ax.scatter([lon], [lat], marker="*", s=220, c=ACCENT, edgecolors="black", linewidths=0.6,
                       transform=ccrs.PlateCarree(), zorder=6)
            ax.annotate(q, (lon, lat), xytext=(6, 6), textcoords="offset points", fontsize=8,
                        color=ACCENT, transform=ccrs.PlateCarree(), zorder=7)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles[:MAX_COLORED + 1], labels[:MAX_COLORED + 1], loc="lower left",
                      fontsize=7, framealpha=0.9)
        ax.set_title(f"{page}: reference records and queries (stars)", fontsize=11)
        _caption(fig, [f"Queries without coordinates: {', '.join(q_missing)}" if q_missing else "",
                       "Reference records with geographic-quality coordinates only."])
        pages.append((page, fig))
    return _save_pages("query_map", pages, out_dir)


def plot_query_regions(regions: pd.DataFrame, region_label: str, out_dir) -> List[Path]:
    """One page per equivalence set: reference records per region (bars), queries marked."""
    import matplotlib.pyplot as plt

    pages = []
    for eq_set, t in regions.groupby("equivalence_set", sort=False):
        t = t.sort_values("n_reference_records")
        fig, ax = plt.subplots(figsize=(9, max(2.5, 0.4 * len(t) + 1.2)))
        y = np.arange(len(t))
        ax.barh(y, t["n_reference_records"], color=REFERENCE_GRAY, label="reference records")
        ax.set_yticks(y)
        ax.set_yticklabels(t["region"])
        xmax = max(1, int(t["n_reference_records"].max()))
        for yi, (n_ref, n, ids) in enumerate(zip(t["n_reference_records"], t["n_queries"],
                                                 t["query_ids"].fillna(""))):
            if n:
                x = n_ref + xmax * 0.02
                ax.scatter([x], [yi], marker="*", s=160, c=ACCENT, zorder=3)
                ax.annotate(f"{n} quer{'y' if n == 1 else 'ies'}: {ids}", (x, yi),
                            xytext=(8, -3), textcoords="offset points", fontsize=8, color=ACCENT)
        ax.set_xlim(0, xmax * 1.6)
        ax.set_xlabel("records of the set's haplotypes in the reference run")
        call = t["species_call"].iloc[0]
        ax.set_title(f"{call + ': ' if isinstance(call, str) and call else ''}{eq_set} by {region_label}",
                     fontsize=11)
        pages.append((eq_set, fig))
    return _save_pages("query_regions", pages, out_dir)


def plot_query_timeline(calls: pd.DataFrame, meta: Optional[pd.DataFrame],
                        reference: Optional[ReferenceRun], out_dir) -> List[Path]:
    """One page per equivalence set: reference records by collection year, query years marked."""
    import matplotlib.pyplot as plt

    if reference is None or "collection_year" not in reference.records.columns:
        return []
    rec = reference.records
    years_by_query = (dict(zip(meta["query_id"], meta["collection_year"]))
                      if meta is not None and "collection_year" in meta.columns else {})
    pages = []
    for eq_set, group in _placed(calls).groupby("equivalence_set", sort=False):
        years = rec.loc[rec["haplotype_id"].isin(eq_set.split(";")), "collection_year"].dropna().astype(int)
        q_years = [(q, int(years_by_query[q])) for q in group["query_id"]
                   if pd.notna(years_by_query.get(q, np.nan))]
        if years.empty and not q_years:
            continue
        counts = years.value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(10, 3.5))
        ax.bar(counts.index, counts.values, color=REFERENCE_GRAY, width=0.8, label="reference records")
        top = max(1, int(counts.max()) if len(counts) else 1)
        for i, (q, yr) in enumerate(sorted(q_years, key=lambda x: x[1])):
            ax.axvline(yr, color=ACCENT, linestyle="--", linewidth=1)
            ax.annotate(q, (yr, top * (1.05 + 0.08 * (i % 4))), xytext=(3, 0),
                        textcoords="offset points", fontsize=8, color=ACCENT, ha="left")
        all_years = list(counts.index) + [yr for _, yr in q_years]
        ax.set_xlim(min(all_years) - 1, max(all_years) + 1.5)
        ax.set_ylim(0, top * 1.45)
        ax.set_xlabel("collection year")
        ax.set_ylabel("records")
        call = group["species_call"].iloc[0]
        ax.set_title(f"{call + ': ' if call else ''}{eq_set}: sampling years, queries dashed",
                     fontsize=11)
        undated = [q for q in group["query_id"] if q not in dict(q_years)]
        _caption(fig, [f"Queries without a collection date: {', '.join(undated)}"
                       if meta is not None and undated else ""])
        pages.append((eq_set, fig))
    return _save_pages("query_timeline", pages, out_dir)


def write_query_figures(calls: pd.DataFrame, meta: Optional[pd.DataFrame],
                        reference: Optional[ReferenceRun], out_dir):
    """
    query_regions.csv and the four figures. Returns (figures, warnings):
    figure name -> written paths, and one message per figure that failed.
    """
    out_dir = Path(out_dir)
    figures, warnings = {}, []
    regions = query_regions_table(calls, meta, reference)
    csv = write_query_regions(regions, out_dir)
    if csv is not None:
        figures["query_regions_table"] = [csv]
    jobs = [
        ("query_tree", lambda: plot_query_tree(calls, reference, out_dir)),
        ("query_map", lambda: plot_query_map(calls, meta, reference, out_dir)),
        ("query_regions", lambda: plot_query_regions(
            regions, reference.region_col if reference else "region", out_dir)),
        ("query_timeline", lambda: plot_query_timeline(calls, meta, reference, out_dir)),
    ]
    for name, job in jobs:
        try:
            paths = job()
        except Exception as e:
            logger.warning(f"{name} not drawn: {e}")
            logger.debug("Figure error details:", exc_info=True)
            warnings.append(f"{name}: {e}")
            continue
        if paths:
            figures[name] = paths
    return figures, warnings
