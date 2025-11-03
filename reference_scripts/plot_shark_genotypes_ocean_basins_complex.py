#!/usr/bin/env python3
"""
Cartopy plots for shark genotypes with optional ocean-basin assignment and basin plotting (v11)

What's new in v11
-----------------
- Territory polygons on BOTH simple and basins plots
- Independent controls for:
  * coastal ring width (vs. blob size)
  * clipping mode: coastal | ocean | none
  * polygon simplification to keep files small
- Line width flags for coastlines and point outlines
- Console hint if a species produces no visible territory geometry

All earlier features kept:
- simple/basins modes, strict (precise-only) maps
- bottom-aligned legends (bubble, genotype colors, shapes)
- basin labels (w/out polygons in simple mode), Set2 palette, color-map CSV
- summaries incl. basin long/pivot
"""

import argparse
import os
import re
import textwrap
from typing import Tuple, Literal, List, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader as shpreader
from cartopy.feature import ShapelyFeature
from shapely.geometry import Point, Polygon
from shapely.ops import nearest_points, unary_union
from shapely.geometry.base import BaseGeometry
from shapely.prepared import prep as shapely_prep
from shapely.geometry import MultiPoint
from shapely.geometry import GeometryCollection, Polygon, MultiPolygon
from shapely.ops import polygonize


# ------------------ Coordinate handling ------------------

def parse_coord_raw(val: str) -> Tuple[float, float]:
    if pd.isna(val):
        return float("nan"), float("nan")
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(val))
    if len(nums) >= 2:
        try:
            return float(nums[0]), float(nums[1])
        except ValueError:
            return float("nan"), float("nan")
    return float("nan"), float("nan")


def normalize_lon(lon: float) -> float:
    if pd.isna(lon):
        return lon
    return lon - 360.0 if lon > 180.0 else lon


def infer_coord_order(a: pd.Series, b: pd.Series) -> Literal["lonlat", "latlon"]:
    valid = ~(a.isna() | b.isna())
    if valid.sum() == 0:
        return "lonlat"
    aa = a[valid].abs()
    bb = b[valid].abs()
    cond_latlon = (aa <= 90).mean() > 0.7 and (bb <= 180).mean() > 0.7
    cond_lonlat = (aa <= 180).mean() > 0.7 and (bb <= 90).mean() > 0.7
    if cond_latlon and not cond_lonlat:
        return "latlon"
    if cond_lonlat and not cond_latlon:
        return "lonlat"
    if (bb > 90).mean() > 0.3 and (aa <= 90).mean() > 0.7:
        return "latlon"
    return "lonlat"


def coerce_lonlat(a: pd.Series, b: pd.Series, order: str) -> pd.DataFrame:
    if order == "latlon":
        lat = a
        lon = b.map(normalize_lon)
    else:
        lon = a.map(normalize_lon)
        lat = b
    return pd.DataFrame({"lon": lon, "lat": lat})


# ------------------ Land detection & snapping ------------------

def load_land_geometries(resolution: str = "110m") -> BaseGeometry:
    shp = shpreader.natural_earth(resolution=resolution, category="physical", name="land")
    reader = shpreader.Reader(shp)
    geoms = list(reader.geometries())
    return unary_union(geoms)


def classify_and_adjust_locations(df: pd.DataFrame,
                                  land_union: BaseGeometry,
                                  policy: Literal["mark", "drop", "snap"] = "mark") -> pd.DataFrame:
    def is_on_land(row):
        pt = Point(float(row["lon"]), float(row["lat"]))
        try:
            return land_union.contains(pt)
        except Exception:
            return False

    df = df.copy()
    df["location_type"] = "precise"
    on_land = df.apply(is_on_land, axis=1)
    df.loc[on_land, "location_type"] = "country_centroid"

    if policy == "drop":
        return df[df["location_type"] == "precise"].copy()

    if policy == "snap":
        land_boundary = land_union.boundary
        def snap_point(row):
            if row["location_type"] == "country_centroid":
                p = Point(float(row["lon"]), float(row["lat"]))
                _, nearest_on_land = nearest_points(p, land_boundary)
                return nearest_on_land.x, nearest_on_land.y
            return row["lon"], row["lat"]
        snapped = df.apply(lambda r: snap_point(r), axis=1, result_type="expand")
        snapped.columns = ["lon", "lat"]
        df[["lon", "lat"]] = snapped
        df.loc[df["location_type"] == "country_centroid", "location_type"] = "snapped_to_coast"

    return df


# ------------------ Legend builders ------------------

def build_color_legend_handles(genotype_order: List[str], color_map: Dict[str, str]):
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker='o', linestyle='None',
                      markerfacecolor=color_map[g], markeredgecolor='black',
                      label=str(g)) for g in genotype_order]
    labels = [str(g) for g in genotype_order]
    return handles, labels


def build_bubble_legend_handles(size_min: float, size_scale: float, example_counts: List[int]):
    handles, labels = [], []
    for c in example_counts:
        s = size_min + size_scale * np.sqrt(c)
        h = plt.scatter([], [], s=s, edgecolor="black", facecolor="none", alpha=1.0)
        handles.append(h); labels.append(f"{c}")
    return handles, labels


def build_shape_legend_handles(include_centroid: bool, include_snapped: bool):
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker='o', linestyle='None', markerfacecolor='white',
                      markeredgecolor='black', label='precise')]
    labels = ['precise']
    if include_centroid:
        handles.append(Line2D([], [], marker='^', linestyle='None',
                              markerfacecolor='white', markeredgecolor='black',
                              label='country centroid'))
        labels.append('country centroid')
    if include_snapped:
        handles.append(Line2D([], [], marker='o', linestyle='None',
                              markerfacecolor='white', markeredgecolor='black',
                              label='snapped to coast', markeredgewidth=2))
        labels.append('snapped to coast')
    return handles, labels


# ------------------ Basin handling ------------------

def load_basin_records(shp_path: str, name_field: str):
    """Load polygons; auto-detect name field if needed; log helpful info."""
    try:
        reader = shpreader.Reader(shp_path)
    except Exception as e:
        print(f"[ERROR] Could not open shapefile: {shp_path} ({e})")
        return []
    recs_iter = reader.records()
    try:
        sample = next(recs_iter)
    except StopIteration:
        print(f"[WARN] Shapefile loaded but contains 0 records: {shp_path}")
        return []
    except Exception as e:
        print(f"[ERROR] Failed to iterate shapefile: {e}")
        return []
    attrs = sample.attributes
    if name_field not in attrs:
        candidates = [k for k, v in attrs.items() if isinstance(v, str)]
        picked = candidates[0] if candidates else None
        print(f"[WARN] Name field '{name_field}' not found. "
              f"Available fields: {list(attrs.keys())}. "
              f"{'Using ' + picked if picked else 'No string field found; using \"unknown\"'}")
        name_field = picked if picked else None
    reader = shpreader.Reader(shp_path)
    recs = []
    for rec in reader.records():
        geom = rec.geometry
        if geom is None:
            continue
        name = rec.attributes.get(name_field, "unknown") if name_field else "unknown"
        recs.append((shapely_prep(geom), geom, str(name)))
    print(f"[INFO] Loaded {len(recs)} basin polygon(s) from: {shp_path}")
    return recs


def assign_basin_to_points(df: pd.DataFrame, basin_recs) -> pd.DataFrame:
    def get_basin(lon, lat):
        if np.isnan(lon) or np.isnan(lat):
            return np.nan
        p = Point(float(lon), float(lat))
        for prepped, geom, name in basin_recs:
            try:
                if prepped.contains(p) or geom.contains(p):
                    return name
            except Exception:
                continue
        return np.nan
    df = df.copy()
    df["ocean_basin"] = [get_basin(lon, lat) for lon, lat in zip(df["lon"], df["lat"])]
    return df


def draw_basins(ax, basin_recs, edgecolor="gray"):
    for _, geom, _ in basin_recs:
        ax.add_geometries([geom], crs=ccrs.PlateCarree(),
                          edgecolor=edgecolor, facecolor="none",
                          linewidth=0.4, alpha=0.5)


def label_basins(ax, basin_recs, fontsize=7):
    # Place a label at each polygon representative point (more robust than centroid)
    for _, geom, name in basin_recs:
        try:
            c = geom.representative_point()
            ax.text(c.x, c.y, str(name),
                    transform=ccrs.PlateCarree(),
                    ha="center", va="center",
                    fontsize=fontsize, alpha=0.8,
                    bbox=dict(facecolor="white", alpha=0.5, edgecolor="none", pad=1.5))
        except Exception:
            continue


def add_basin_legend(ax, basin_recs, legend_fontsize=8, legend_title_fontsize=9):
    import matplotlib.patches as mpatches
    names = sorted({name for _, _, name in basin_recs})
    palette = plt.get_cmap("tab20").colors
    name_to_color = {nm: palette[i % len(palette)] for i, nm in enumerate(names)}
    handles = [mpatches.Patch(facecolor=name_to_color[nm], edgecolor="gray", alpha=0.25, label=nm)
               for nm in names]
    leg = ax.legend(handles=handles, title="Ocean basins",
                    loc="upper left", frameon=True,
                    prop={'size': legend_fontsize}, title_fontsize=legend_title_fontsize)
    for lh in leg.legend_handles:
        lh.set_alpha(0.6)
    ax.add_artist(leg)


# ------------------ Territory polygons ------------------

def world_polygon():
    # Rectangle covering lon/lat world extent
    return Polygon([(-180, -90), (-180, 90), (180, 90), (180, -90)])


def build_clip_geometry(land_union: BaseGeometry, mode: str, ring_deg: float):
    """
    Returns the clipping geometry, or None if no clipping.
    - 'coastal': near-coast ocean band (buffer ring)
    - 'ocean'  : world ocean (world minus land)
    - 'none'   : no clipping
    """
    if mode == "none":
        return None
    if mode == "ocean":
        return world_polygon().difference(land_union)
    # coastal ring
    if ring_deg is None or ring_deg <= 0:
        return None
    expanded = land_union.buffer(ring_deg)
    return expanded.difference(land_union)


def territory_from_points(points_lonlat: np.ndarray,
                          method: Literal["buffer", "convex"] = "buffer",
                          buffer_deg: float = 0.5):
    """
    Build a genotype territory from lon/lat points.
    - 'buffer': unary_union of point buffers (radius=buffer_deg), smooth blob
    - 'convex': convex hull of all points; if <3 points, returns small buffer
    """
    pts = [Point(float(lon), float(lat)) for lon, lat in points_lonlat]
    if not pts:
        return None
    if method == "convex":
        if len(pts) < 3:
            return MultiPoint(pts).buffer(max(1e-3, buffer_deg * 0.5))
        return MultiPoint(pts).convex_hull
    # buffer method
    return unary_union([p.buffer(buffer_deg) for p in pts])


def draw_genotype_territories(ax, species_df: pd.DataFrame, color_map: Dict[str, str],
                              land_union: BaseGeometry, territory_params: dict):
    """
    Build & draw genotype territories with optional clipping & simplification.
    territory_params keys:
      method, buffer_deg, ring_deg, clip_mode, simplify_deg, min_area_deg2,
      alpha_fill, edgewidth, overlap_hatch (bool)
    """
    ring_deg = territory_params.get("ring_deg", territory_params["buffer_deg"])
    clip_mode = territory_params.get("clip_mode", "coastal")
    simplify_tol = territory_params.get("simplify_deg", 0.0)
    min_area = territory_params.get("min_area_deg2", 0.0)

    clip_geom = build_clip_geometry(land_union, clip_mode, ring_deg)
    territories = {}
    debug_counts = []

    for g, sdf in species_df.groupby("consensus_group"):
        pts = sdf[["lon", "lat"]].to_numpy()
        poly = territory_from_points(pts, method=territory_params["method"],
                                     buffer_deg=territory_params["buffer_deg"])
        if not poly or poly.is_empty:
            debug_counts.append((g, 0, "no-blob"))
            continue

        # Clip (optional)
        geom = poly if clip_geom is None else poly.intersection(clip_geom)
        if geom.is_empty:
            debug_counts.append((g, 0, "empty-after-clip"))
            continue

        # Keep only polygonal area
        geom = _polygonize_only(geom)
        if not geom or geom.is_empty:
            debug_counts.append((g, 0, "no-polys-after-polygonize"))
            continue

        # Simplify (optional)
        if simplify_tol and simplify_tol > 0:
            try:
                geom = geom.simplify(simplify_tol, preserve_topology=True)
            except Exception:
                pass
            if geom.is_empty:
                debug_counts.append((g, 0, "empty-after-simplify"))
                continue

        # Filter very small crumbs
        geom = _filter_min_area(geom, min_area)
        if not geom or geom.is_empty:
            debug_counts.append((g, 0, "dropped-by-min-area"))
            continue

        # Count parts for debug
        n_parts = len(list(geom.geoms)) if isinstance(geom, MultiPolygon) else 1
        debug_counts.append((g, n_parts, "ok"))
        territories[g] = geom

    if not territories:
        print("[INFO] Territories: no visible polygonal area after clip/polygonize/simplify.")
    else:
        msg = ", ".join([f"{g}:{n}({tag})" for g, n, tag in debug_counts])
        print(f"[INFO] Territories built → {msg}")

    # draw filled polygons (via ShapelyFeature — more robust than add_geometries)
    for g, geom in territories.items():
        feat = ShapelyFeature([geom], ccrs.PlateCarree(),
                              facecolor=color_map[g],
                              edgecolor=color_map[g],
                              linewidth=territory_params["edgewidth"])
        ax.add_feature(feat, zorder=2.0, alpha=territory_params["alpha_fill"])
        
    # visualize overlaps (optional) via pairwise intersections
    if territory_params.get("overlap_hatch", False):
        keys = list(territories.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                inter = territories[keys[i]].intersection(territories[keys[j]])
                inter = _polygonize_only(inter)
                inter = _filter_min_area(inter, min_area)
                if not inter or inter.is_empty:
                    continue
                ax.add_geometries([inter], crs=ccrs.PlateCarree(),
                                  facecolor="none", edgecolor="none",
                                  hatch="////", linewidth=0.0, zorder=2.5, alpha=0.0)
    return list(territories.items())
        
def _polygonize_only(geom):
    """
    Return only polygonal area from an arbitrary geometry:
    - If already Polygon/MultiPolygon, return as-is
    - If GeometryCollection (or lines), polygonize to build area faces
    - Otherwise, return None if no area is present
    """
    if geom.is_empty:
        return None
    if isinstance(geom, (Polygon, MultiPolygon)):
        return geom
    # Try to polygonize any linear parts
    try:
        polys = list(polygonize(geom))
        if not polys:
            return None
        return unary_union(polys)
    except Exception:
        return None


def _filter_min_area(geom, min_area):
    """Drop polygon parts below a minimum area (in degree^2)."""
    if not geom or geom.is_empty:
        return None
    if min_area <= 0:
        return geom
    try:
        if isinstance(geom, Polygon):
            return geom if geom.area >= min_area else None
        if isinstance(geom, MultiPolygon):
            kept = [g for g in geom.geoms if g.area >= min_area]
            if not kept:
                return None
            return unary_union(kept)
    except Exception:
        pass
    return geom


# ------------------ Color mapping ------------------

def parse_legend_counts(arg: str) -> List[int]:
    try:
        vals = [int(x.strip()) for x in arg.split(",") if x.strip()]
        return [v for v in vals if v > 0]
    except Exception:
        return [1, 5, 25]


def get_color_map_for_genotypes(genotypes: List[str], palette: str, colors_arg: str) -> Dict[str, str]:
    if colors_arg:
        provided = [c.strip() for c in colors_arg.split(",") if c.strip()]
        return {g: provided[i % len(provided)] for i, g in enumerate(genotypes)}
    colors = None
    try:
        cmap = plt.get_cmap(palette)
        try:
            base = cmap.colors
        except Exception:
            base = [cmap(i / max(1, len(genotypes) - 1)) for i in range(max(1, len(genotypes)))]
        colors = base
    except Exception:
        try:
            import seaborn as sns
            if palette.lower() == "set2":
                colors = sns.color_palette("Set2", n_colors=max(8, len(genotypes)))
            elif palette.lower() == "dark2":
                colors = sns.color_palette("dark", n_colors=max(8, len(genotypes)))
        except Exception:
            colors = None
    if colors is None:
        try:
            colors = plt.get_cmap("tab20").colors
        except Exception:
            colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                      "#bcbd22", "#17becf"]
    return {g: colors[i % len(colors)] for i, g in enumerate(genotypes)}


def save_color_map_csv(color_map: Dict[str, str], outdir: str):
    import matplotlib.colors as mcolors
    rows = []
    for g, col in color_map.items():
        try:
            hexcol = mcolors.to_hex(col, keep_alpha=False)
        except Exception:
            hexcol = str(col)
        rows.append({"consensus_group": g, "color": hexcol})
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "genotype_color_map.csv"), index=False)


# ------------------ Plotting (per-species) ------------------

def plot_species(species_df: pd.DataFrame, species_name: str, outdir: str,
                 alpha: float, size_min: float, size_scale: float,
                 legend_counts: List[int],
                 color_map: Dict[str, str],
                 color_ncol: Optional[int],
                 legend_fontsize: int, legend_title_fontsize: int,
                 plot_mode: str, basin_recs=None,
                 show_basin_legend=False, label_basin_names=False, basin_label_size=7,
                 coastline_lw: float = 0.6, point_lw: float = 0.6,
                 draw_territories: bool = False, territory_params: dict = None,
                 land_union: BaseGeometry = None):

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(11, 5.5), dpi=300)
    ax = plt.axes(projection=proj)
    ax.set_global()
    ax.coastlines(resolution="110m", linewidth=coastline_lw)
    ax.add_feature(cfeature.LAND, facecolor="lightgray", alpha=0.5)

    if basin_recs:
        if plot_mode == "basins":
            draw_basins(ax, basin_recs)
            if show_basin_legend:
                add_basin_legend(ax, basin_recs, legend_fontsize, legend_title_fontsize)
            if label_basin_names:
                label_basins(ax, basin_recs, fontsize=basin_label_size)
        else:
            if label_basin_names:
                label_basins(ax, basin_recs, fontsize=basin_label_size)

    # Territories (draw beneath points)
    if draw_territories and land_union is not None:
        _ = draw_genotype_territories(ax, species_df, color_map, land_union, territory_params)

    # Aggregate counts by exact lon/lat + genotype + location_type
    agg = (species_df.groupby(["lon", "lat", "consensus_group", "location_type"], as_index=False)
           .size().rename(columns={"size": "n_at_coord"}))

    genotype_order_local = sorted(agg["consensus_group"].dropna().unique())

    # Plot precise/snapped (circles) and centroid (triangles), both colored by genotype
    for loc_type, marker in [("precise", "o"), ("snapped_to_coast", "o"), ("country_centroid", "^")]:
        sub_all = agg[agg["location_type"] == loc_type]
        if sub_all.empty:
            continue
        for g in genotype_order_local:
            sub = sub_all[sub_all["consensus_group"] == g].copy()
            if sub.empty:
                continue
            sub["marker_size"] = size_min + size_scale * np.sqrt(sub["n_at_coord"].astype(float))
            ax.scatter(sub["lon"], sub["lat"], s=sub["marker_size"],
                       alpha=alpha, edgecolor="black", linewidth=point_lw,
                       label=None, c=[color_map[g]], marker=marker, transform=proj, zorder=3.0)

    subtitle = "basins" if plot_mode == "basins" else "coastlines"
    if draw_territories:
        subtitle += " + territories"
    ax.set_title(textwrap.fill(f"{species_name} — global distribution of COI genotypes ({subtitle})", 80))
    plt.subplots_adjust(bottom=0.28)

    # ---- Figure-level legends (bottom row) ----
    color_handles, color_labels = build_color_legend_handles(genotype_order_local, color_map)
    bubble_handles, bubble_labels = build_bubble_legend_handles(size_min, size_scale, legend_counts)
    include_centroid = (agg["location_type"] == "country_centroid").any()
    include_snapped = (agg["location_type"] == "snapped_to_coast").any()
    shape_handles, shape_labels = build_shape_legend_handles(include_centroid, include_snapped)

    fig.legend(bubble_handles, bubble_labels, title="Count at coord", loc="lower left",
               bbox_to_anchor=(0.06, 0.01), frameon=True, ncol=1,
               prop={'size': legend_fontsize}, title_fontsize=legend_title_fontsize)
    ncol_colors = (min(len(color_handles), 6) if color_ncol is None else color_ncol) or 1
    fig.legend(color_handles, color_labels, title="Genotype (consensus_group)", loc="lower center",
               bbox_to_anchor=(0.5, 0.01), frameon=True, ncol=ncol_colors,
               prop={'size': legend_fontsize}, title_fontsize=legend_title_fontsize)
    fig.legend(shape_handles, shape_labels, title="Location type", loc="lower right",
               bbox_to_anchor=(0.94, 0.01), frameon=True, ncol=1,
               prop={'size': legend_fontsize}, title_fontsize=legend_title_fontsize)

    os.makedirs(outdir, exist_ok=True)
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", species_name.strip())
    fig.savefig(os.path.join(outdir, f"{safe_name}_genotypes.svg"), bbox_inches="tight", pad_inches=0.08)
    fig.savefig(os.path.join(outdir, f"{safe_name}_genotypes.pdf"), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


# ------------------ Main ------------------

def main():
    ap = argparse.ArgumentParser(description="Cartopy plots with basins, legends, strict maps, territories, and summaries (v11).")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--coord-order", choices=["auto", "lonlat", "latlon"], default="auto")
    ap.add_argument("--alpha", type=float, default=0.6)
    ap.add_argument("--size-min", type=float, default=12.0)
    ap.add_argument("--size-scale", type=float, default=6.0)
    ap.add_argument("--land-policy", choices=["mark", "drop", "snap"], default="mark")
    ap.add_argument("--legend-counts", type=str, default="1,5,25")
    ap.add_argument("--palette", type=str, default="tab20")
    ap.add_argument("--colors", type=str, default="")
    ap.add_argument("--color-legend-ncol", type=int, default=0)
    ap.add_argument("--legend-fontsize", type=int, default=8)
    ap.add_argument("--legend-title-fontsize", type=int, default=9)
    ap.add_argument("--plot-mode", choices=["simple", "basins"], default="simple")
    ap.add_argument("--basin-shp", type=str, default="")
    ap.add_argument("--basin-name-field", type=str, default="BASIN")
    ap.add_argument("--show-basin-legend", action="store_true")
    ap.add_argument("--label-basins", action="store_true")
    ap.add_argument("--basin-label-size", type=int, default=7)

    # Line widths
    ap.add_argument("--coastline-lw", type=float, default=0.6, help="Coastline line width")
    ap.add_argument("--point-lw", type=float, default=0.6, help="Point outline width")

    # Territories
    ap.add_argument("--draw-territories", action="store_true", help="Draw genotype territories")
    ap.add_argument("--territory-method", choices=["buffer", "convex"], default="buffer",
                    help="Territory geometry method")
    ap.add_argument("--territory-buffer-deg", type=float, default=0.5,
                    help="Buffer radius in degrees for territory blobs (approx distance)")
    ap.add_argument("--territory-alpha-fill", type=float, default=0.25,
                    help="Fill alpha for territory polygons")
    ap.add_argument("--territory-edgewidth", type=float, default=1.0,
                    help="Outline width for territory polygons")
    ap.add_argument("--territory-overlap-hatch", action="store_true",
                    help="Hatch areas where genotype territories overlap")

    # NEW in v11: decoupled ring width, clip mode, simplification
    ap.add_argument("--territory-coastal-ring-deg", type=float, default=None,
                    help="Width of near-coast ring (deg) used for clipping. "
                         "Defaults to --territory-buffer-deg if unset.")
    ap.add_argument("--territory-clip", choices=["coastal", "ocean", "none"], default="coastal",
                    help="Clip territories to 'coastal' ring, all 'ocean' (world minus land), or 'none'.")
    ap.add_argument("--territory-simplify-deg", type=float, default=0.0,
                    help="Simplify tolerance (deg) for territory polygons to shrink file size.")

    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    # Parse coord → lon/lat
    AB = df["coord"].apply(parse_coord_raw).apply(pd.Series).rename(columns={0: "A", 1: "B"})
    df = pd.concat([df, AB], axis=1)
    order = args.coord_order if args.coord_order != "auto" else infer_coord_order(df["A"], df["B"])
    lonlat = coerce_lonlat(df["A"], df["B"], order)
    df = pd.concat([df, lonlat], axis=1)
    df = df.dropna(subset=["lon", "lat", "consensus_group", "COI_common"]).copy()

    # Land detection & policy
    land_union = load_land_geometries(resolution="110m")
    df = classify_and_adjust_locations(df, land_union, policy=args.land_policy)

    # Basins: load whenever shapefile exists (labels in simple mode; assignment & draw in basins mode)
    basin_recs = None
    if args.basin_shp and os.path.exists(args.basin_shp):
        print(f"[INFO] Loading basin shapefile: {args.basin_shp}")
        basin_recs = load_basin_records(args.basin_shp, args.basin_name_field)

    # Basin assignment for analysis tables if polygons available
    if basin_recs:
        df = assign_basin_to_points(df, basin_recs)
    else:
        df["ocean_basin"] = np.nan

    # Colors
    all_genotypes = sorted(df["consensus_group"].dropna().unique())
    color_map = get_color_map_for_genotypes(all_genotypes, palette=args.palette, colors_arg=args.colors)

    # Output dir & color map CSV
    os.makedirs(args.outdir, exist_ok=True)
    save_color_map_csv(color_map, args.outdir)

    # -------- Summaries --------
    summary_species = (
        df.groupby(["COI_common", "consensus_group"])
          .size().reset_index(name="n_samples")
          .sort_values(["COI_common", "n_samples"], ascending=[True, False])
    )
    summary_species.to_csv(os.path.join(args.outdir, "genotype_counts_by_species.csv"), index=False)

    summary_species_coord = (
        df.groupby(["COI_common", "lon", "lat", "consensus_group", "location_type"])
          .size().reset_index(name="n_at_coord")
          .sort_values(["COI_common", "n_at_coord"], ascending=[True, False])
    )
    summary_species_coord.to_csv(os.path.join(args.outdir, "genotype_counts_by_species_coord.csv"), index=False)

    species_totals = (
        df.groupby(["COI_common", "location_type"])
          .size().reset_index(name="n_samples")
          .sort_values(["COI_common", "location_type"], ascending=[True, True])
    )
    species_totals.to_csv(os.path.join(args.outdir, "species_totals_by_location_type.csv"), index=False)

    rec_col_candidates = ["recorded_species", "recorded_common", "Recorded_species", "recorded_name"]
    coi_col_candidates = ["COI_species", "COI_scientific", "COI_name", "COI_Species"]
    recorded_species_col = next((c for c in rec_col_candidates if c in df.columns), "COI_common")
    coi_species_col = next((c for c in coi_col_candidates if c in df.columns), "COI_common")

    basin_long = (
        df.groupby(["consensus_group", recorded_species_col, coi_species_col, "ocean_basin"])
          .size().reset_index(name="n")
          .rename(columns={recorded_species_col: "recorded_species",
                           coi_species_col: "COI_species"})
          .sort_values(["consensus_group", "n"], ascending=[True, False])
    )
    basin_long.to_csv(os.path.join(args.outdir, "genotype_basin_summary_long.csv"), index=False)

    basin_pivot = basin_long.pivot_table(
        index="consensus_group",
        columns=["recorded_species", "COI_species", "ocean_basin"],
        values="n", aggfunc="sum", fill_value=0
    )
    basin_pivot.to_csv(os.path.join(args.outdir, "genotype_basin_summary_pivot.csv"))

    # -------- Plotting --------
    legend_counts = parse_legend_counts(args.legend_counts)
    color_ncol = None if args.color_legend_ncol == 0 else int(args.color_legend_ncol)
    outdir_plots = os.path.join(args.outdir, "basins" if args.plot_mode == "basins" else "simple")
    os.makedirs(outdir_plots, exist_ok=True)

    territory_params = dict(
        method=args.territory_method,
        buffer_deg=args.territory_buffer_deg,
        alpha_fill=args.territory_alpha_fill,
        edgewidth=args.territory_edgewidth,
        overlap_hatch=args.territory_overlap_hatch,
        ring_deg=(args.territory_coastal_ring_deg
                  if args.territory_coastal_ring_deg is not None
                  else args.territory_buffer_deg),
        clip_mode=args.territory_clip,
        simplify_deg=args.territory_simplify_deg,
    )

    # Full set
    for species_name, species_df in df.groupby("COI_common", sort=True):
        plot_species(
            species_df, species_name, outdir_plots,
            alpha=args.alpha, size_min=args.size_min, size_scale=args.size_scale,
            legend_counts=legend_counts,
            color_map=color_map, color_ncol=color_ncol,
            legend_fontsize=args.legend_fontsize, legend_title_fontsize=args.legend_title_fontsize,
            plot_mode=args.plot_mode,
            basin_recs=basin_recs,
            show_basin_legend=(args.show_basin_legend and args.plot_mode == "basins"),
            label_basin_names=args.label_basins, basin_label_size=args.basin_label_size,
            coastline_lw=args.coastline_lw, point_lw=args.point_lw,
            draw_territories=args.draw_territories, territory_params=territory_params,
            land_union=land_union
        )

    # Strict (precise-only)
    strict_subdir = os.path.join(outdir_plots, "strict_precise_only")
    os.makedirs(strict_subdir, exist_ok=True)
    df_strict = df[df["location_type"] == "precise"].copy()
    for species_name, species_df in df_strict.groupby("COI_common", sort=True):
        plot_species(
            species_df, species_name, strict_subdir,
            alpha=args.alpha, size_min=args.size_min, size_scale=args.size_scale,
            legend_counts=legend_counts,
            color_map=color_map, color_ncol=color_ncol,
            legend_fontsize=args.legend_fontsize, legend_title_fontsize=args.legend_title_fontsize,
            plot_mode=args.plot_mode,  # subtitle follows mode; points are precise-only
            basin_recs=basin_recs,
            show_basin_legend=False,  # keep strict plots uncluttered
            label_basin_names=(args.label_basins and bool(basin_recs)),
            basin_label_size=args.basin_label_size,
            coastline_lw=args.coastline_lw, point_lw=args.point_lw,
            draw_territories=args.draw_territories, territory_params=territory_params,
            land_union=land_union
        )


if __name__ == "__main__":
    main()