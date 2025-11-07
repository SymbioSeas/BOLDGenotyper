"""
Geographic Analysis and Ocean Basin Assignment

This module handles geographic analysis of COI samples using the General
Oceanographic Areas System (GOaS) shapefiles from the United Nations Food
and Agriculture Organization (FAO).

The GOaS system provides standardized ocean basin delineations used for:
- Ocean basin assignment for genetic samples
- Marine vs. terrestrial coordinate validation
- Biogeographic analysis and visualization
- Quality control of geographic metadata

Key Features:
- Loads and validates GOaS shapefiles
- Performs point-in-polygon spatial joins for basin assignment
- Handles edge cases (samples outside basins, on boundaries, missing data)
- Provides fallback mechanisms when GOaS data unavailable
- Generates basin statistics and reports

Spatial Operations:
- Coordinate Reference System (CRS): WGS84 (EPSG:4326)
- Point-in-polygon for basin assignment
- Spatial joins with GeoPandas
- Boundary handling with buffering

Edge Cases Handled:
1. Coordinates outside all defined basins → "Unknown"
2. Coordinates on basin boundaries → First match selected
3. Missing/invalid coordinates → Skipped with warning
4. Terrestrial coordinates → Flagged if validate_marine=True
5. Missing GOaS data → Graceful fallback with instructions

GOaS Data Setup:
If GOaS shapefiles are not available, download from:
https://www.fao.org/geonetwork/srv/en/main.home

Or use the setup script:
    python scripts/setup_goas.py --download

Example Usage:
    >>> from boldgenotyper.geographic import load_goas_data, assign_ocean_basins
    >>> from boldgenotyper.metadata import parse_bold_tsv
    >>>
    >>> # Load data
    >>> df = parse_bold_tsv("sphyrna_lewini.tsv")
    >>> goas = load_goas_data("data/goas/goas_v1.shp")
    >>>
    >>> # Assign ocean basins
    >>> df_with_basins = assign_ocean_basins(
    ...     df, goas,
    ...     lat_col='lat',
    ...     lon_col='lon',
    ...     country_ocean_col='country/ocean'
    ... )
    >>>
    >>> # Get basin statistics
    >>> counts = get_basin_counts(df_with_basins)
    >>> for basin, count in counts.items():
    ...     print(f"{basin}: {count}")

Author: Steph Smith (steph.smith@unc.edu)
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import logging
import warnings
import ast
import re

import pandas as pd
import numpy as np

# Conditional imports for geospatial libraries
try:
    import geopandas as gpd
    from shapely.geometry import Point
    from shapely.ops import nearest_points, unary_union
    from shapely.geometry.base import BaseGeometry
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False
    gpd = None
    Point = None
    nearest_points = None
    unary_union = None
    BaseGeometry = None

# Conditional import for cartopy
try:
    from cartopy.io import shapereader
    CARTOPY_AVAILABLE = True
except ImportError:
    CARTOPY_AVAILABLE = False
    shapereader = None

# Configure logging
logger = logging.getLogger(__name__)


class GOaSDataError(Exception):
    """Raised when GOaS shapefile data cannot be loaded or is invalid."""
    pass


class GeospatialLibraryError(Exception):
    """Raised when required geospatial libraries are not available."""
    pass


def check_geopandas_available() -> None:
    """
    Check if GeoPandas and Shapely are available.

    Raises
    ------
    GeospatialLibraryError
        If geopandas or shapely are not installed

    Notes
    -----
    GeoPandas requires several dependencies:
    - geopandas
    - shapely
    - pyproj
    - fiona

    Install with: pip install geopandas
    """
    if not GEOPANDAS_AVAILABLE:
        raise GeospatialLibraryError(
            "GeoPandas and Shapely are required for geographic operations.\n"
            "Install with: pip install geopandas\n"
            "This will also install shapely, pyproj, and fiona."
        )


def check_cartopy_available() -> None:
    """
    Check if Cartopy is available.

    Raises
    ------
    GeospatialLibraryError
        If cartopy is not installed

    Notes
    -----
    Cartopy is required for land/ocean classification and snap-to-coast features.

    Install with: pip install cartopy
    """
    if not CARTOPY_AVAILABLE:
        raise GeospatialLibraryError(
            "Cartopy is required for land/ocean classification.\n"
            "Install with: pip install cartopy\n"
        )


def load_goas_data(shapefile_path: Union[str, Path]) -> "gpd.GeoDataFrame":
    """
    Load GOaS (General Oceanographic Areas System) shapefile data.

    Parameters
    ----------
    shapefile_path : str or Path
        Path to GOaS shapefile (.shp)

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing ocean basin polygons with standardized CRS (EPSG:4326)

    Raises
    ------
    GeospatialLibraryError
        If geopandas is not installed
    GOaSDataError
        If shapefile cannot be loaded or is invalid
    FileNotFoundError
        If shapefile does not exist

    Notes
    -----
    The GOaS shapefile must contain:
    - Polygon geometries for ocean basins
    - A name field for basin identification

    The loaded data is reprojected to WGS84 (EPSG:4326) if needed.

    GOaS Data Sources:
    - FAO GeoNetwork: https://www.fao.org/geonetwork/srv/en/main.home
    - Search for "General Oceanographic Areas System"

    Examples
    --------
    >>> goas = load_goas_data("data/goas/goas_v1.shp")
    >>> print(f"Loaded {len(goas)} ocean basins")
    >>> print(goas.columns)
    """
    check_geopandas_available()

    shapefile_path = Path(shapefile_path)

    # Check if file exists
    if not shapefile_path.exists():
        raise FileNotFoundError(
            f"GOaS shapefile not found: {shapefile_path}\n\n"
            "Please download GOaS shapefiles from:\n"
            "https://www.fao.org/geonetwork/srv/en/main.home\n\n"
            "Or use the setup script:\n"
            "    python scripts/setup_goas.py --download\n"
        )

    logger.info(f"Loading GOaS shapefile from: {shapefile_path}")

    try:
        # Load shapefile
        goas_gdf = gpd.read_file(shapefile_path)

        # Validate GeoDataFrame
        if len(goas_gdf) == 0:
            raise GOaSDataError("GOaS shapefile is empty")

        if goas_gdf.geometry is None or goas_gdf.geometry.empty:
            raise GOaSDataError("GOaS shapefile has no geometry")

        # Check for valid CRS
        if goas_gdf.crs is None:
            logger.warning("GOaS shapefile has no CRS, assuming WGS84 (EPSG:4326)")
            goas_gdf.set_crs("EPSG:4326", inplace=True)

        # Reproject to WGS84 if needed
        if goas_gdf.crs.to_epsg() != 4326:
            logger.info(f"Reprojecting GOaS data from {goas_gdf.crs} to EPSG:4326")
            goas_gdf = goas_gdf.to_crs("EPSG:4326")

        # Identify basin name column (common variations)
        name_columns = ['name', 'NAME', 'basin', 'BASIN', 'ocean', 'OCEAN',
                       'area_name', 'AREA_NAME', 'region', 'REGION']

        basin_name_col = None
        for col in name_columns:
            if col in goas_gdf.columns:
                basin_name_col = col
                break

        if basin_name_col is None:
            # Use first non-geometry column as basin name
            non_geom_cols = [c for c in goas_gdf.columns if c != 'geometry']
            if non_geom_cols:
                basin_name_col = non_geom_cols[0]
                logger.warning(
                    f"Could not find standard basin name column, using: {basin_name_col}"
                )
            else:
                raise GOaSDataError("GOaS shapefile has no attribute columns")

        # Standardize column name
        if basin_name_col != 'basin_name':
            goas_gdf = goas_gdf.rename(columns={basin_name_col: 'basin_name'})

        logger.info(
            f"Loaded {len(goas_gdf)} ocean basins from GOaS shapefile"
        )

        # Log basin names
        if 'basin_name' in goas_gdf.columns:
            basins = sorted(goas_gdf['basin_name'].unique())
            logger.debug(f"Ocean basins: {', '.join(basins)}")

        return goas_gdf

    except Exception as e:
        if isinstance(e, (GOaSDataError, FileNotFoundError)):
            raise
        raise GOaSDataError(f"Failed to load GOaS shapefile: {e}")


def create_points_geodataframe(
    df: pd.DataFrame,
    lat_col: str = 'lat',
    lon_col: str = 'lon',
    crs: str = "EPSG:4326"
) -> "gpd.GeoDataFrame":
    """
    Create a GeoDataFrame of point geometries from latitude/longitude columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with latitude and longitude columns
    lat_col : str, default='lat'
        Name of latitude column
    lon_col : str, default='lon'
        Name of longitude column
    crs : str, default="EPSG:4326"
        Coordinate reference system (WGS84)

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with Point geometries

    Raises
    ------
    GeospatialLibraryError
        If geopandas is not installed
    ValueError
        If latitude or longitude columns are missing

    Notes
    -----
    - Rows with missing (NaN) coordinates are preserved but have null geometry
    - Invalid coordinates (out of range) are logged as warnings
    - Valid ranges: -90 ≤ lat ≤ 90, -180 ≤ lon ≤ 180

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     'sample_id': ['S1', 'S2'],
    ...     'lat': [23.5, -15.8],
    ...     'lon': [81.2, -47.9]
    ... })
    >>> gdf = create_points_geodataframe(df)
    >>> print(gdf.geometry)
    """
    check_geopandas_available()

    # Validate columns exist
    if lat_col not in df.columns:
        raise ValueError(f"Latitude column '{lat_col}' not found in DataFrame")
    if lon_col not in df.columns:
        raise ValueError(f"Longitude column '{lon_col}' not found in DataFrame")

    # Create copy to avoid modifying original
    df_copy = df.copy()

    # Count valid coordinates
    valid_coords = df_copy[[lat_col, lon_col]].notna().all(axis=1)
    n_valid = valid_coords.sum()
    n_total = len(df_copy)

    logger.info(
        f"Creating point geometries for {n_valid}/{n_total} samples with coordinates"
    )

    # Create Point geometries
    geometry = []
    for idx, row in df_copy.iterrows():
        lat = row[lat_col]
        lon = row[lon_col]

        # Handle missing coordinates
        if pd.isna(lat) or pd.isna(lon):
            geometry.append(None)
            continue

        # Validate coordinate ranges
        if not (-90 <= lat <= 90):
            logger.warning(
                f"Latitude {lat} out of range [-90, 90] for row {idx}, "
                "creating null geometry"
            )
            geometry.append(None)
            continue

        if not (-180 <= lon <= 180):
            logger.warning(
                f"Longitude {lon} out of range [-180, 180] for row {idx}, "
                "creating null geometry"
            )
            geometry.append(None)
            continue

        # Create Point geometry
        geometry.append(Point(lon, lat))  # Note: Point(x, y) = Point(lon, lat)

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(df_copy, geometry=geometry, crs=crs)

    return gdf


def assign_ocean_basins(
    df: pd.DataFrame,
    goas_data: Optional["gpd.GeoDataFrame"] = None,
    shapefile_path: Optional[Union[str, Path]] = None,
    lat_col: str = 'lat',
    lon_col: str = 'lon',
    coord_col: Optional[str] = None,
    country_ocean_col: Optional[str] = None,
    basin_col: str = 'ocean_basin',
    fallback_to_country_ocean: bool = True,
    validate_marine: bool = False
) -> pd.DataFrame:
    """
    Assign ocean basins to samples using GOaS spatial data.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with latitude/longitude coordinates
    goas_data : gpd.GeoDataFrame, optional
        Pre-loaded GOaS GeoDataFrame. If None, will load from shapefile_path
    shapefile_path : str or Path, optional
        Path to GOaS shapefile. Required if goas_data is None
    lat_col : str, default='lat'
        Name of latitude column
    lon_col : str, default='lon'
        Name of longitude column
    country_ocean_col : str, optional
        Name of country/ocean column for fallback
    basin_col : str, default='ocean_basin'
        Name of output basin assignment column
    fallback_to_country_ocean : bool, default=True
        If True, use country/ocean column for samples without coordinates
    validate_marine : bool, default=False
        If True, flag terrestrial coordinates

    Returns
    -------
    pd.DataFrame
        DataFrame with added ocean basin assignment column

    Raises
    ------
    GeospatialLibraryError
        If geopandas is not installed
    ValueError
        If neither goas_data nor shapefile_path is provided

    Notes
    -----
    Basin Assignment Logic:
    1. For samples with coordinates → spatial point-in-polygon
    2. For samples without coordinates and fallback enabled → use country/ocean column
    3. For samples outside all basins → "Unknown"
    4. For samples on boundaries → first match selected
    
    Coordinate handling:
    - If {lat_col, lon_col} are missing but a single coordinate column is present (e.g., "coord" with values like "[-33.8, 151.283]"), set coord_col to that column and it will be parsed to create {lat_col, lon_col}.

    Missing Data Handling:
    - Samples with missing coordinates get basin from country/ocean column (if available)
    - Samples outside all basins are assigned "Unknown"
    - Warning logged for high proportion of "Unknown" assignments

    Examples
    --------
    >>> from boldgenotyper.metadata import parse_bold_tsv
    >>> df = parse_bold_tsv("data.tsv")
    >>> goas = load_goas_data("data/goas/goas_v1.shp")
    >>> df_basins = assign_ocean_basins(df, goas)
    >>> print(df_basins['ocean_basin'].value_counts())
    """
    check_geopandas_available()

    # Load GOaS data if not provided
    if goas_data is None:
        if shapefile_path is None:
            raise ValueError(
                "Either goas_data or shapefile_path must be provided"
            )
        goas_data = load_goas_data(shapefile_path)

    # Create copy to avoid modifying original
    df_copy = df.copy()
    if basin_col not in df_copy.columns:
        df_copy[basin_col] = pd.NA
        
    # Add lat/lon coord columns if needed
    need_latlon = (lat_col not in df_copy.columns) or (lon_col not in df_copy.columns)
    if need_latlon and coord_col is None:
        if 'coord' in df_copy.columns:
            coord_col = 'coord'
        
    if need_latlon and coord_col is not None and coord_col in df_copy.columns:
        # parse strings like "[-33.8, 151.283]" or without brackets
        def _parse_pair(s: str):
            if pd.isna(s):
                return (None, None)
            text = str(s).strip()
            try:
                v = ast.literal_eval(text)
                if isinstance(v, (list, tuple)) and len(v) == 2:
                    return (float(v[0]), float(v[1]))
            except Exception:
                pass
            m = re.search(r'(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)', text)
            if m:
                return (float(m.group(1)), float(m.group(2)))
            return (None, None)

        latlon = df_copy[coord_col].apply(_parse_pair)
        df_copy[lat_col] = latlon.apply(lambda t: t[0])
        df_copy[lon_col] = latlon.apply(lambda t: t[1])

    # Initialize basin column
    df_copy[basin_col] = None

    # Create point geometries
    logger.info("Creating point geometries from coordinates")
    points_gdf = create_points_geodataframe(df_copy, lat_col, lon_col)

    # Filter to samples with valid geometries
    has_geometry = points_gdf.geometry.notna()
    n_with_coords = has_geometry.sum()

    logger.info(
        f"Performing spatial join for {n_with_coords}/{len(df_copy)} samples "
        "with valid coordinates"
    )

    if n_with_coords > 0:
        # Spatial join: point-in-polygon
        # Use 'left' join to preserve all points
        joined = gpd.sjoin(
            points_gdf[has_geometry],
            goas_data[['basin_name', 'geometry']],
            how='left',
            predicate='within'
        )

        # Extract basin assignments
        # Note: sjoin adds index_right column for matched polygon
        if 'basin_name' in joined.columns:
            # Update basin assignments for samples with coordinates
            df_copy.loc[has_geometry, basin_col] = joined['basin_name'].values

        # Log assignment statistics
        assigned = joined['basin_name'].notna().sum()
        unassigned = joined['basin_name'].isna().sum()

        logger.info(
            f"Spatial join results: {assigned} assigned to basins, "
            f"{unassigned} outside all basins"
        )

        if unassigned > 0 and assigned > 0:
            pct_unassigned = (unassigned / n_with_coords) * 100
            if pct_unassigned > 20:
                logger.warning(
                    f"{pct_unassigned:.1f}% of coordinates fall outside GOaS basins. "
                    "This may indicate terrestrial coordinates or incomplete basin coverage."
                )

    # Fallback to country/ocean column for samples without coordinates
    if fallback_to_country_ocean and country_ocean_col is not None:
        if country_ocean_col in df_copy.columns:
            # Identify samples without basin assignment
            no_basin = df_copy[basin_col].isna()
            has_country_ocean = df_copy[country_ocean_col].notna()

            fallback_mask = no_basin & has_country_ocean
            n_fallback = fallback_mask.sum()

            if n_fallback > 0:
                logger.info(
                    f"Using {country_ocean_col} for {n_fallback} samples "
                    "without coordinate-based basin assignment"
                )
                df_copy.loc[fallback_mask, basin_col] = df_copy.loc[
                    fallback_mask, country_ocean_col
                ]

    # Set remaining null basins to "Unknown"
    still_missing = df_copy[basin_col].isna()
    if still_missing.any():
        df_copy.loc[still_missing, basin_col] = "Unknown"
        logger.info(
            f"Assigned {still_missing.sum()} samples to 'Unknown' basin "
            "(no coordinates or outside all basins)"
        )

    # Validate marine coordinates if requested
    if validate_marine:
        logger.info("Validating marine coordinates")
        # Samples assigned to "Unknown" with coordinates might be terrestrial
        unknown_with_coords = (
            (df_copy[basin_col] == "Unknown") &
            df_copy[[lat_col, lon_col]].notna().all(axis=1)
        )

        if unknown_with_coords.any():
            n_potential_terrestrial = unknown_with_coords.sum()
            logger.warning(
                f"Found {n_potential_terrestrial} samples with coordinates "
                "outside ocean basins (potentially terrestrial)"
            )

    return df_copy


def get_basin_counts(
    df: pd.DataFrame,
    basin_col: str = 'ocean_basin',
    sort_by: str = 'count'
) -> Dict[str, int]:
    """
    Get counts of samples per ocean basin.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with basin assignments
    basin_col : str, default='ocean_basin'
        Name of basin column
    sort_by : str, default='count'
        Sort order: 'count' (descending), 'name' (alphabetical), or None

    Returns
    -------
    Dict[str, int]
        Dictionary mapping basin names to sample counts

    Examples
    --------
    >>> counts = get_basin_counts(df)
    >>> for basin, count in counts.items():
    ...     print(f"{basin}: {count}")
    North Atlantic Ocean: 45
    South Pacific Ocean: 32
    Indian Ocean: 18
    """
    if basin_col not in df.columns:
        raise ValueError(f"Basin column '{basin_col}' not found in DataFrame")

    counts = df[basin_col].value_counts().to_dict()

    # Sort if requested
    if sort_by == 'count':
        counts = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))
    elif sort_by == 'name':
        counts = dict(sorted(counts.items(), key=lambda x: x[0]))

    return counts


def get_basin_summary_stats(
    df: pd.DataFrame,
    basin_col: str = 'ocean_basin',
    lat_col: str = 'lat',
    lon_col: str = 'lon'
) -> pd.DataFrame:
    """
    Generate summary statistics for each ocean basin.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with basin assignments and coordinates
    basin_col : str, default='ocean_basin'
        Name of basin column
    lat_col : str, default='lat'
        Name of latitude column
    lon_col : str, default='lon'
        Name of longitude column

    Returns
    -------
    pd.DataFrame
        Summary statistics with columns:
        - basin: Basin name
        - n_samples: Number of samples
        - lat_min, lat_max, lat_mean: Latitude statistics
        - lon_min, lon_max, lon_mean: Longitude statistics

    Examples
    --------
    >>> stats = get_basin_summary_stats(df)
    >>> print(stats)
                        basin  n_samples  lat_min  lat_max  lat_mean  lon_min  lon_max  lon_mean
    0  North Atlantic Ocean         45    10.5     65.2      35.8    -80.1    -10.3      -45.2
    1   South Pacific Ocean         32   -45.2    -5.1     -25.6    -170.5    -70.2     -120.3
    """
    if basin_col not in df.columns:
        raise ValueError(f"Basin column '{basin_col}' not found in DataFrame")

    # Group by basin
    grouped = df.groupby(basin_col)

    stats = []
    for basin_name, group in grouped:
        # Count samples
        n_samples = len(group)

        # Calculate coordinate statistics (only for valid coordinates)
        valid_coords = group[[lat_col, lon_col]].notna().all(axis=1)

        if valid_coords.any():
            lats = group.loc[valid_coords, lat_col]
            lons = group.loc[valid_coords, lon_col]

            basin_stats = {
                'basin': basin_name,
                'n_samples': n_samples,
                'lat_min': lats.min(),
                'lat_max': lats.max(),
                'lat_mean': lats.mean(),
                'lon_min': lons.min(),
                'lon_max': lons.max(),
                'lon_mean': lons.mean(),
            }
        else:
            basin_stats = {
                'basin': basin_name,
                'n_samples': n_samples,
                'lat_min': None,
                'lat_max': None,
                'lat_mean': None,
                'lon_min': None,
                'lon_max': None,
                'lon_mean': None,
            }

        stats.append(basin_stats)

    stats_df = pd.DataFrame(stats)

    # Sort by sample count descending
    stats_df = stats_df.sort_values('n_samples', ascending=False).reset_index(drop=True)

    return stats_df


def validate_marine_coordinates(
    lat: float,
    lon: float,
    goas_data: "gpd.GeoDataFrame",
    buffer_degrees: float = 0.0
) -> bool:
    """
    Validate if coordinates fall within marine areas (ocean basins).

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees
    lon : float
        Longitude in decimal degrees
    goas_data : gpd.GeoDataFrame
        GOaS ocean basin data
    buffer_degrees : float, default=0.0
        Buffer distance in degrees for boundary tolerance

    Returns
    -------
    bool
        True if coordinates fall within ocean basins, False otherwise

    Raises
    ------
    GeospatialLibraryError
        If geopandas is not installed

    Notes
    -----
    This function is useful for quality control to identify potentially
    terrestrial coordinates in marine species datasets.

    Examples
    --------
    >>> goas = load_goas_data("data/goas/goas_v1.shp")
    >>> # Marine coordinate (Gulf of Mexico)
    >>> validate_marine_coordinates(23.5, -81.2, goas)
    True
    >>> # Terrestrial coordinate (Mexico City)
    >>> validate_marine_coordinates(19.4, -99.1, goas)
    False
    """
    check_geopandas_available()

    # Validate input coordinates
    if not (-90 <= lat <= 90):
        raise ValueError(f"Latitude {lat} out of range [-90, 90]")
    if not (-180 <= lon <= 180):
        raise ValueError(f"Longitude {lon} out of range [-180, 180]")

    # Create point geometry
    point = Point(lon, lat)

    # Check if point falls within any ocean basin
    if buffer_degrees > 0:
        # Buffer the point for boundary tolerance
        point_buffered = point.buffer(buffer_degrees)
        is_marine = goas_data.geometry.intersects(point_buffered).any()
    else:
        is_marine = goas_data.geometry.contains(point).any()

    return is_marine


def filter_marine_samples(
    df: pd.DataFrame,
    goas_data: "gpd.GeoDataFrame",
    lat_col: str = 'lat',
    lon_col: str = 'lon',
    remove_terrestrial: bool = True
) -> pd.DataFrame:
    """
    Filter samples to include only marine coordinates.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with coordinates
    goas_data : gpd.GeoDataFrame
        GOaS ocean basin data
    lat_col : str, default='lat'
        Name of latitude column
    lon_col : str, default='lon'
        Name of longitude column
    remove_terrestrial : bool, default=True
        If True, remove terrestrial samples. If False, add flag column.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame (if remove_terrestrial=True) or
        DataFrame with 'is_marine' column (if remove_terrestrial=False)

    Raises
    ------
    GeospatialLibraryError
        If geopandas is not installed

    Examples
    --------
    >>> goas = load_goas_data("data/goas/goas_v1.shp")
    >>> # Remove terrestrial samples
    >>> marine_df = filter_marine_samples(df, goas)
    >>> # Or just flag them
    >>> flagged_df = filter_marine_samples(df, goas, remove_terrestrial=False)
    >>> print(flagged_df[~flagged_df['is_marine']])
    """
    check_geopandas_available()

    df_copy = df.copy()

    # Create point geometries
    points_gdf = create_points_geodataframe(df_copy, lat_col, lon_col)

    # Perform spatial join to identify marine samples
    joined = gpd.sjoin(
        points_gdf,
        goas_data[['geometry']],
        how='left',
        predicate='within'
    )

    # Samples within ocean basins have non-null index_right
    is_marine = joined.index_right.notna()

    n_marine = is_marine.sum()
    n_terrestrial = (~is_marine).sum()

    logger.info(
        f"Marine validation: {n_marine} marine, {n_terrestrial} terrestrial samples"
    )

    if remove_terrestrial:
        df_filtered = df_copy[is_marine].copy()
        return df_filtered
    else:
        df_copy['is_marine'] = is_marine
        return df_copy


def load_land_geometries(resolution: str = "110m") -> "BaseGeometry":
    """
    Load land geometries from Natural Earth data.

    Parameters
    ----------
    resolution : str, default="110m"
        Resolution of Natural Earth data: "110m", "50m", or "10m"
        - 110m: Low resolution, faster loading (recommended for testing)
        - 50m: Medium resolution
        - 10m: High resolution, slower loading

    Returns
    -------
    BaseGeometry
        Unified land geometry (MultiPolygon)

    Raises
    ------
    GeospatialLibraryError
        If geopandas, shapely, or cartopy are not installed

    Notes
    -----
    This function downloads Natural Earth land shapefiles on first use.
    The data is cached by cartopy for subsequent use.

    Natural Earth Data: https://www.naturalearthdata.com/

    Examples
    --------
    >>> land = load_land_geometries(resolution="110m")
    >>> # Check if a point is on land
    >>> from shapely.geometry import Point
    >>> point = Point(-99.1, 19.4)  # Mexico City
    >>> is_land = land.contains(point)
    >>> print(is_land)  # True
    """
    check_geopandas_available()
    check_cartopy_available()

    logger.info(f"Loading Natural Earth land geometries (resolution: {resolution})")

    try:
        # Suppress cartopy download warnings (these are informational only)
        import warnings
        from cartopy.io import DownloadWarning
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DownloadWarning)

            # Download/load Natural Earth land shapefile
            shapefile = shapereader.natural_earth(
                resolution=resolution,
                category="physical",
                name="land"
            )

        # Read geometries
        reader = shapereader.Reader(shapefile)
        geometries = list(reader.geometries())

        logger.info(f"Loaded {len(geometries)} land polygons")

        # Combine all land geometries into a single union
        land_union = unary_union(geometries)

        logger.info("Successfully created unified land geometry")

        return land_union

    except Exception as e:
        raise GeospatialLibraryError(
            f"Failed to load Natural Earth land geometries: {e}\n"
            "Make sure cartopy is properly installed: pip install cartopy"
        )


def classify_and_snap_coordinates(
    df: pd.DataFrame,
    lat_col: str = 'latitude',
    lon_col: str = 'longitude',
    policy: str = 'mark',
    resolution: str = '110m',
    location_type_col: str = 'location_type',
    land_geometry: Optional["BaseGeometry"] = None
) -> pd.DataFrame:
    """
    Classify coordinates as land or ocean and optionally snap land points to coast.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with latitude and longitude columns
    lat_col : str, default='latitude'
        Name of latitude column
    lon_col : str, default='longitude'
        Name of longitude column
    policy : str, default='mark'
        How to handle land-based coordinates:
        - 'mark': Label land points as 'country_centroid', keep coordinates unchanged
        - 'drop': Remove all land-based coordinates
        - 'snap': Snap land points to nearest coastline
    resolution : str, default='110m'
        Natural Earth data resolution: '110m', '50m', or '10m'
        Ignored if land_geometry is provided
    location_type_col : str, default='location_type'
        Name of output column for location type classification
    land_geometry : BaseGeometry, optional
        Pre-loaded land geometry (from load_land_geometries()).
        If provided, this geometry will be used instead of loading from Natural Earth.
        This allows reusing the same geometry across multiple calls for better performance.

    Returns
    -------
    pd.DataFrame
        DataFrame with added location_type column and updated coordinates (if policy='snap')
        Location types:
        - 'precise': Ocean coordinates (original)
        - 'country_centroid': Land coordinates (if policy='mark')
        - 'snapped_to_coast': Land coordinates snapped to coast (if policy='snap')

    Raises
    ------
    GeospatialLibraryError
        If required libraries are not installed
    ValueError
        If invalid policy specified

    Notes
    -----
    This function is useful for marine species datasets where some coordinates
    may incorrectly fall on land due to imprecise GPS or data entry errors.

    The 'snap' policy finds the nearest point on the coastline for each land
    coordinate and updates the latitude/longitude to that coastal location.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'sample': ['Ocean', 'InlandAustralia', 'Coastal'],
    ...     'latitude': [0.0, -25.0, -33.0],
    ...     'longitude': [-150.0, 133.0, 151.0]
    ... })
    >>>
    >>> # Mark land points
    >>> df_marked = classify_and_snap_coordinates(df, policy='mark')
    >>> print(df_marked[['sample', 'location_type']])
    >>>
    >>> # Snap land points to coast
    >>> df_snapped = classify_and_snap_coordinates(df, policy='snap')
    >>> print(df_snapped[['sample', 'latitude', 'longitude', 'location_type']])
    >>>
    >>> # Drop land points
    >>> df_ocean_only = classify_and_snap_coordinates(df, policy='drop')
    >>> print(len(df_ocean_only))  # Only ocean points remain
    >>>
    >>> # Reuse land geometry for better performance with multiple calls
    >>> land = load_land_geometries(resolution="110m")
    >>> df_marked = classify_and_snap_coordinates(df, land_geometry=land, policy='mark')
    >>> df_snapped = classify_and_snap_coordinates(df, land_geometry=land, policy='snap')
    """
    check_geopandas_available()
    check_cartopy_available()

    # Validate policy
    valid_policies = ['mark', 'drop', 'snap']
    if policy not in valid_policies:
        raise ValueError(
            f"Invalid policy '{policy}'. Must be one of: {valid_policies}"
        )

    # Validate columns exist
    if lat_col not in df.columns:
        raise ValueError(f"Latitude column '{lat_col}' not found in DataFrame")
    if lon_col not in df.columns:
        raise ValueError(f"Longitude column '{lon_col}' not found in DataFrame")

    # Create copy to avoid modifying original
    df_copy = df.copy()

    # Load land geometries if not provided
    logger.info(f"Classifying coordinates (policy: {policy})")
    if land_geometry is not None:
        logger.debug("Using pre-loaded land geometry")
        land_union = land_geometry
    else:
        land_union = load_land_geometries(resolution=resolution)

    # Initialize location_type column
    df_copy[location_type_col] = 'precise'

    # Classify each point as land or ocean
    def is_on_land(row):
        """Check if a coordinate point falls on land."""
        try:
            lat = float(row[lat_col])
            lon = float(row[lon_col])

            # Skip invalid coordinates
            if pd.isna(lat) or pd.isna(lon):
                return False

            # Create point geometry
            point = Point(lon, lat)  # Note: Point(x, y) = Point(lon, lat)

            # Check if point is contained in land union
            return land_union.contains(point)

        except Exception as e:
            logger.warning(f"Error checking land status for row: {e}")
            return False

    # Apply land classification
    logger.info("Checking which coordinates fall on land...")
    on_land = df_copy.apply(is_on_land, axis=1)
    n_land = on_land.sum()
    n_ocean = (~on_land).sum()

    logger.info(
        f"Classification complete: {n_ocean} ocean points, {n_land} land points"
    )

    # Mark land points
    df_copy.loc[on_land, location_type_col] = 'country_centroid'

    # Handle based on policy
    if policy == 'drop':
        logger.info(f"Dropping {n_land} land-based coordinates")
        return df_copy[df_copy[location_type_col] == 'precise'].copy()

    elif policy == 'snap':
        if n_land > 0:
            logger.info(f"Snapping {n_land} land points to nearest coastline...")

            # Get land boundary (coastline)
            land_boundary = land_union.boundary

            def snap_point(row):
                """Snap a land point to the nearest coastline."""
                if row[location_type_col] == 'country_centroid':
                    try:
                        lat = float(row[lat_col])
                        lon = float(row[lon_col])
                        point = Point(lon, lat)

                        # Find nearest point on coastline
                        _, nearest_coastal_point = nearest_points(point, land_boundary)

                        # Return new coordinates (lon, lat)
                        return nearest_coastal_point.x, nearest_coastal_point.y

                    except Exception as e:
                        logger.warning(f"Error snapping point: {e}")
                        return row[lon_col], row[lat_col]

                # Return original coordinates for ocean points
                return row[lon_col], row[lat_col]

            # Apply snapping
            snapped_coords = df_copy.apply(
                snap_point, axis=1, result_type='expand'
            )
            snapped_coords.columns = [lon_col, lat_col]

            # Update coordinates
            df_copy[[lon_col, lat_col]] = snapped_coords

            # Update location type for snapped points
            df_copy.loc[
                df_copy[location_type_col] == 'country_centroid',
                location_type_col
            ] = 'snapped_to_coast'

            logger.info("Snapping complete")

    return df_copy
