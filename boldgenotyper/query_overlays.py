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
