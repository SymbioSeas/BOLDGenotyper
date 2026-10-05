"""
BOLD TSV Parsing and Coordinate Handling

This module parses BOLD TSV files and handles their coordinates.

Key Responsibilities:

1. Parse BOLD TSV files and extract critical fields:

   - processid: Unique sample identifier (REQUIRED)
   - nuc: COI nucleotide sequence (REQUIRED)
   - coord: Geographic coordinates in format [lat, lon]
   - coord_source: Indicates if coordinates are from centroid
   - country/ocean: Location metadata
   - Other taxonomic and collection metadata

2. Coordinates:

   - Parse coordinate columns and check their ranges
   - Mark coordinate quality (ok, missing, unparseable, zero, out of range, centroid)

3. Data Validation:

   - Check for required fields
   - Handle missing or malformed data gracefully
   - Provide informative warnings for data quality issues

Important Notes:

- Duplicate processids are removed, keeping the first record, with a warning.

Example Usage:
    >>> from boldgenotyper.metadata import parse_bold_tsv
    >>> df = parse_bold_tsv("Sphyrna_lewini.tsv")

Author: Steph Smith (symbioseas@outlook.com)
"""

from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path
import logging
import re
import pandas as pd
import numpy as np
import csv

from . import utils

# Configure logging
logger = logging.getLogger(__name__)


# ============================================================================
# BOLD TSV Parsing
# ============================================================================

def parse_bold_tsv(
    tsv_path: Union[str, Path],
    required_columns: Optional[List[str]] = None,
    encoding: str = 'utf-8',
) -> pd.DataFrame:
    """
    Parse BOLD TSV file and validate required fields.

    Handles various BOLD export formats, encoding issues, and missing columns.
    This validates your BOLD data and ensure it contains the necessary information for genotyping and geographic analysis.

    Parameters
    ----------
    tsv_path : Union[str, Path]
        Path to BOLD TSV file
    required_columns : Optional[List[str]]
        List of required column names. If None, uses default ['processid', 'nuc']
    encoding : str
        File encoding (default: 'utf-8', falls back to 'latin-1' if needed)

    Returns
    -------
    pd.DataFrame
        Parsed metadata with validated fields

    Raises
    ------
    FileNotFoundError
        If TSV file doesn't exist
    ValueError
        If required fields are missing or processid is not unique
    pd.errors.EmptyDataError
        If file is empty

    Examples
    --------
    >>> df = parse_bold_tsv("Sphyrna_lewini.tsv")
    >>> print(df.shape)
    (1500, 86)
    >>> print(df['processid'].nunique())
    1500

    Notes
    -----
    BOLD TSV files typically contain ~86 columns. This function:
    - Handles encoding issues (tries UTF-8, falls back to Latin-1)
    - Validates required columns exist
    - Checks for duplicate processids
    - Converts dtypes where appropriate
    - Logs warnings for data quality issues
    """
    path = Path(tsv_path)

    if not path.exists():
        raise FileNotFoundError(f"BOLD TSV file not found: {path}")

    # Default required columns
    if required_columns is None:
        required_columns = ['processid', 'nuc']

    logger.info(f"Reading BOLD TSV file: {path}")

    def _read_with_fallback(enc: str) -> pd.DataFrame:
        """Read TSV with fallbacks for malformed rows."""
        try:
            return pd.read_csv(
                path,
                sep='\t',
                encoding=enc,
                dtype=str,  # Read all as strings initially
                low_memory=False
            )
        except pd.errors.ParserError as e:
            logger.warning(
                f"Standard TSV parse failed (encoding={enc}): {e}. "
                "Retrying with Python engine and skipping malformed rows."
            )
            return pd.read_csv(
                path,
                sep='\t',
                encoding=enc,
                dtype=str,
                engine='python',
                on_bad_lines='skip',
                quoting=csv.QUOTE_MINIMAL
            )

    # Try reading with specified encoding, fallback to latin-1 if needed
    try:
        df = _read_with_fallback(encoding)
    except UnicodeDecodeError:
        logger.warning("UTF-8 decoding failed; retrying with latin-1 encoding.")
        df = _read_with_fallback('latin-1')

    # Check if file is empty
    if df.empty:
        raise pd.errors.EmptyDataError(f"BOLD TSV file is empty: {path}")

    logger.info(f"Read {len(df)} rows and {len(df.columns)} columns")

    # Validate required columns
    validate_required_columns(df, required_columns)

    # Check for duplicate processids and remove them (keep first occurrence)
    if 'processid' in df.columns:
        duplicates = df['processid'].duplicated()
        if duplicates.any():
            n_duplicates = duplicates.sum()
            dup_ids = df.loc[duplicates, 'processid'].head(5).tolist()
            logger.warning(
                f"Found {n_duplicates} duplicate processids. "
                f"Examples: {dup_ids}. "
                f"Keeping only first occurrence of each duplicate."
            )
            # Remove duplicates, keeping first occurrence
            df = df[~duplicates].copy()
            logger.info(f"After removing duplicates: {len(df)} rows remaining")

    # Clean up whitespace in string columns
    for col in df.columns:
        if utils.is_text_series(df[col]):
            df[col] = df[col].str.strip()

    # Log data quality statistics
    _log_data_quality(df)

    return df


def validate_required_columns(
    df: pd.DataFrame,
    required_columns: List[str]
) -> bool:
    """
    Validate that DataFrame contains required columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate
    required_columns : List[str]
        List of required column names

    Returns
    -------
    bool
        True if all required columns present

    Raises
    ------
    ValueError
        If any required columns are missing
    """
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        available_cols = sorted(df.columns.tolist())
        logger.error(
            f"Missing required columns: {missing_columns}\n"
            f"Available columns: {available_cols[:20]}..."
        )
        raise ValueError(
            f"BOLD TSV is missing required columns: {missing_columns}. "
            f"Found {len(df.columns)} columns total."
        )

    logger.debug(f"All required columns present: {required_columns}")
    return True


def _log_data_quality(df: pd.DataFrame) -> None:
    """Log data quality statistics for BOLD TSV."""
    stats = {
        'total_rows': len(df),
        'total_columns': len(df.columns),
    }

    # Check key columns
    if 'processid' in df.columns:
        stats['unique_processids'] = df['processid'].nunique()

    if 'nuc' in df.columns:
        stats['sequences_present'] = df['nuc'].notna().sum()
        stats['sequences_missing'] = df['nuc'].isna().sum()
        stats['empty_sequences'] = (df['nuc'] == '').sum()

    if 'coord' in df.columns:
        stats['coords_present'] = df['coord'].notna().sum()
        stats['coords_missing'] = df['coord'].isna().sum()

    if 'coord_source' in df.columns:
        centroids = df['coord_source'].str.contains('centroid', case=False, na=False).sum()
        stats['centroid_coords'] = centroids

    logger.info("Data quality summary:")
    for key, value in stats.items():
        logger.info(f"  {key}: {value}")


# ============================================================================
# Coordinate Parsing and Validation
# ============================================================================

def extract_coordinates(
    coord_string: Union[str, float],
    validate_ranges: bool = True
) -> Optional[Tuple[float, float]]:
    """
    Parse BOLD coordinate string from format '[lat, lon]' to floats.

    Handles various coordinate formats from BOLD exports including:
    - '[34.5, -76.2]' (standard format)
    - '34.5, -76.2' (without brackets)
    - '34.5 -76.2' (space-separated)
    - Various whitespace variations

    Parameters
    ----------
    coord_string : Union[str, float]
        Coordinate string from BOLD coord column
    validate_ranges : bool
        If True, validate lat/lon ranges (default: True)

    Returns
    -------
    Optional[Tuple[float, float]]
        (latitude, longitude) tuple, or None if parsing fails

    Examples
    --------
    >>> extract_coordinates('[34.5, -76.2]')
    (34.5, -76.2)
    >>> extract_coordinates('34.5, -76.2')
    (34.5, -76.2)
    >>> extract_coordinates('[0, 0]')
    (0.0, 0.0)
    >>> extract_coordinates('invalid')
    None

    Notes
    -----
    Valid latitude range: -90 to 90
    Valid longitude range: -180 to 180
    Returns None for invalid coordinates without raising an error
    """
    # Handle NaN or None
    if pd.isna(coord_string) or coord_string is None:
        return None

    # Convert to string if needed
    coord_str = str(coord_string).strip()

    # Empty string
    if not coord_str or coord_str.lower() in ['nan', 'none', '']:
        return None

    try:
        # Remove brackets if present
        coord_str = coord_str.strip('[](){} ')

        # Split by comma or whitespace
        if ',' in coord_str:
            parts = coord_str.split(',')
        else:
            parts = coord_str.split()

        if len(parts) != 2:
            logger.debug(f"Coordinate string has wrong number of parts: '{coord_string}'")
            return None

        # Parse values
        val1 = float(parts[0].strip())
        val2 = float(parts[1].strip())

        # Validate ranges if requested
        if validate_ranges:
            if not (-90 <= val1 <= 90):
                logger.debug(f"First value out of latitude range: {val1}")
                return None

            if not (-180 <= val2 <= 180):
                logger.debug(f"Second value out of longitude range: {val2}")
                return None

        return (val1, val2)

    except (ValueError, AttributeError, TypeError) as e:
        logger.debug(f"Failed to parse coordinates '{coord_string}': {e}")
        return None


def infer_coord_order(coord_pairs: List[Tuple[float, float]]) -> str:
    """
    Infer whether coordinates are in lat/lon or lon/lat order.

    Uses heuristics to determine the most likely coordinate order:
    - Latitude must be between -90 and 90
    - Longitude can be between -180 and 180
    - If first values exceed ±90, they must be longitudes

    Parameters
    ----------
    coord_pairs : List[Tuple[float, float]]
        List of coordinate pairs to analyze

    Returns
    -------
    str
        Either 'lat_lon' or 'lon_lat' indicating detected order

    Examples
    --------
    >>> coords = [(34.5, -76.2), (-33.9, 151.2), (41.2, -8.9)]
    >>> infer_coord_order(coords)
    'lat_lon'
    >>> coords = [(-76.2, 34.5), (151.2, -33.9), (-8.9, 41.2)]
    >>> infer_coord_order(coords)
    'lon_lat'

    Notes
    -----
    The inference checks if the first value in pairs exceeds ±90, which
    would indicate it's longitude. If any definitive evidence is found,
    that determines the order. Otherwise defaults to 'lat_lon'.
    """
    if not coord_pairs:
        logger.warning("No coordinate pairs provided for order inference, assuming lat/lon")
        return 'lat_lon'

    # Count definitive evidence for each order
    lon_lat_evidence = 0  # first value > 90 (must be lon)
    lat_lon_evidence = 0  # second value > 90 (must be lon)

    for first, second in coord_pairs:
        # If first value is outside latitude range, it MUST be longitude (lon/lat order)
        if abs(first) > 90:
            lon_lat_evidence += 1
        # If second value is outside latitude range, it MUST be longitude (lat/lon order)
        if abs(second) > 90:
            lat_lon_evidence += 1

    # If we have definitive evidence, use it
    if lon_lat_evidence > 0 and lat_lon_evidence == 0:
        logger.info(f"Inferred coordinate order: lon/lat ({lon_lat_evidence} definitive indicators)")
        return 'lon_lat'
    elif lat_lon_evidence > 0 and lon_lat_evidence == 0:
        logger.info(f"Inferred coordinate order: lat/lon ({lat_lon_evidence} definitive indicators)")
        return 'lat_lon'
    elif lon_lat_evidence > 0 and lat_lon_evidence > 0:
        # Conflicting evidence - this shouldn't happen with consistent data
        logger.warning(
            f"Conflicting coordinate order evidence detected! "
            f"lon/lat indicators: {lon_lat_evidence}, lat/lon indicators: {lat_lon_evidence}. "
            f"Defaulting to lat/lon."
        )
        return 'lat_lon'
    else:
        # No definitive evidence (all values within [-90, 90])
        # Default to lat/lon as it's the more common convention
        logger.info(f"No definitive coordinate order evidence found. Defaulting to lat/lon.")
        return 'lat_lon'


def parse_coordinates_column(
    df: pd.DataFrame,
    coord_column: str = 'coord',
    infer_order: bool = False
) -> pd.DataFrame:
    """
    Parse coordinates column and add latitude/longitude columns.

    Adds two new columns to the DataFrame:
    - latitude: float
    - longitude: float

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with coordinate column
    coord_column : str
        Name of coordinate column (default: 'coord')
    infer_order : bool
        If True, automatically detect and correct lon/lat vs lat/lon order (default: False)

    Returns
    -------
    pd.DataFrame
        DataFrame with added latitude and longitude columns

    Examples
    --------
    >>> df = pd.DataFrame({'coord': ['[34.5, -76.2]', '[0, 0]', None]})
    >>> df_parsed = parse_coordinates_column(df)
    >>> print(df_parsed[['latitude', 'longitude']])

    >>> # With order inference
    >>> df = pd.DataFrame({'coord': ['[-76.2, 34.5]', '[151.2, -33.9]']})
    >>> df_parsed = parse_coordinates_column(df, infer_order=True)
    >>> print(df_parsed[['latitude', 'longitude']])
    """
    if coord_column not in df.columns:
        logger.warning(f"Coordinate column '{coord_column}' not found")
        df['latitude'] = np.nan
        df['longitude'] = np.nan
        return df

    logger.info(f"Parsing coordinates from '{coord_column}' column")

    # If infer_order is True, extract without validation first
    if infer_order:
        # Extract coordinates without range validation
        coords_raw = df[coord_column].apply(
            lambda x: extract_coordinates(x, validate_ranges=False)
        )

        # Get valid coordinate pairs for inference
        valid_coords = [c for c in coords_raw if c is not None]

        if valid_coords:
            detected_order = infer_coord_order(valid_coords)

            if detected_order == 'lon_lat':
                logger.info("Swapping coordinates from lon/lat to lat/lon order")
                # Swap the coordinate pairs
                coords_raw = coords_raw.apply(lambda x: (x[1], x[0]) if x is not None else None)

        # Now validate the (possibly swapped) coordinates
        coords = coords_raw.apply(lambda x: _validate_coordinate_ranges(x) if x is not None else None)
    else:
        # Apply coordinate extraction with validation
        coords = df[coord_column].apply(extract_coordinates)

    # Split into separate columns
    df['latitude'] = coords.apply(lambda x: x[0] if x is not None else np.nan)
    df['longitude'] = coords.apply(lambda x: x[1] if x is not None else np.nan)

    # Log parsing statistics
    n_parsed = coords.notna().sum()
    n_failed = coords.isna().sum()
    logger.info(f"Parsed {n_parsed} coordinates successfully, {n_failed} failed")

    return df


def _validate_coordinate_ranges(coord_pair: Tuple[float, float]) -> Optional[Tuple[float, float]]:
    """
    Validate that coordinate pair has valid latitude/longitude ranges.

    Parameters
    ----------
    coord_pair : Tuple[float, float]
        (latitude, longitude) tuple

    Returns
    -------
    Optional[Tuple[float, float]]
        The coordinate pair if valid, None otherwise
    """
    lat, lon = coord_pair

    # Validate latitude range
    if not (-90 <= lat <= 90):
        logger.debug(f"Latitude out of range: {lat}")
        return None

    # Validate longitude range
    if not (-180 <= lon <= 180):
        logger.debug(f"Longitude out of range: {lon}")
        return None

    return coord_pair


def filter_by_coordinate_quality(
    df: pd.DataFrame,
    exclude_centroids: bool = True,
    exclude_zero_coords: bool = True,
    exclude_missing: bool = True,
    coord_source_column: str = 'coord_source',
) -> pd.DataFrame:
    """
    Filter samples based on coordinate quality.
    
    Uses filtering rules to exclude samples with ambiguous coordinates. This is ESSENTIAL 
    for accurate biogeographic analysis because many BOLD records contain country-level 
    centroids rather than actual collection locations.

    Filtering Rules:
    1. Exclude samples with missing coordinates
    2. Exclude samples with centroid coordinates (coord_source contains 'centroid')
    3. Exclude samples with [0, 0] coordinates (common placeholder)
    4. Exclude samples with invalid coordinate ranges

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed latitude/longitude columns
    exclude_centroids : bool
        Exclude samples where coord_source contains 'centroid' (default: True)
    exclude_zero_coords : bool
        Exclude samples at coordinates [0, 0] (default: True)
    exclude_missing : bool
        Exclude samples with missing coordinates (default: True)
    coord_source_column : str
        Name of coordinate source column (default: 'coord_source')

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only high-quality coordinates

    Examples
    --------
    >>> df = parse_bold_tsv("data.tsv")
    >>> df = parse_coordinates_column(df)
    >>> df_filtered = filter_by_coordinate_quality(df)
    >>> print(f"Retained {len(df_filtered)}/{len(df)} samples")

    Notes
    -----
    This filtering is critical for accurate geography assignment.
    Many BOLD records use country centroids instead of actual 
    collection locations, which can lead to incorrect biogeographic 
    conclusions if not properly filtered.
    """
    n_initial = len(df)
    df_filtered = df.copy()

    # Track filtering reasons
    filter_reasons = {}

    # 1. Exclude missing coordinates
    if exclude_missing and 'latitude' in df_filtered.columns:
        mask_missing = df_filtered['latitude'].isna() | df_filtered['longitude'].isna()
        n_missing = mask_missing.sum()
        df_filtered = df_filtered[~mask_missing]
        filter_reasons['missing_coordinates'] = n_missing
        logger.info(f"Excluded {n_missing} samples with missing coordinates")

    # 2. Exclude centroid coordinates
    if exclude_centroids and coord_source_column in df_filtered.columns:
        mask_centroid = df_filtered[coord_source_column].str.contains(
            'centroid', case=False, na=False
        )
        n_centroid = mask_centroid.sum()
        df_filtered = df_filtered[~mask_centroid]
        filter_reasons['centroid_coordinates'] = n_centroid
        logger.info(f"Excluded {n_centroid} samples with centroid coordinates")

    # 3. Exclude [0, 0] coordinates
    if exclude_zero_coords and 'latitude' in df_filtered.columns:
        mask_zero = (df_filtered['latitude'] == 0) & (df_filtered['longitude'] == 0)
        n_zero = mask_zero.sum()
        df_filtered = df_filtered[~mask_zero]
        filter_reasons['zero_coordinates'] = n_zero
        logger.info(f"Excluded {n_zero} samples with [0, 0] coordinates")

    # Summary
    n_final = len(df_filtered)
    n_excluded = n_initial - n_final
    pct_retained = (n_final / n_initial * 100) if n_initial > 0 else 0

    logger.info(
        f"Coordinate filtering: {n_final}/{n_initial} samples retained "
        f"({pct_retained:.1f}%), {n_excluded} excluded"
    )

    return df_filtered


_COORD_PAIR = re.compile(r"^\s*\[?\s*(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)\s*\]?\s*$")


def parse_coord_column(
    df: pd.DataFrame,
    coord_col: str = "coord",
    lat_col: str = "lat",
    lon_col: str = "lon",
) -> pd.DataFrame:
    """
    Add float lat/lon columns for every record with a parseable coordinate pair.

    Accepts "[lat, lon]" (current BOLD portal exports; integer or decimal
    values, any sign) and "lat, lon" (older exports). Existing non-empty
    lat/lon values are kept (legacy exports with separate columns). Records
    without a parseable pair get NaN.

    Parameters
    ----------
    df : pd.DataFrame
        Records, typically from parse_bold_tsv
    coord_col, lat_col, lon_col : str
        Column names

    Returns
    -------
    pd.DataFrame
        Copy of df with float lat_col and lon_col
    """
    out = df.copy()
    lat = pd.Series(np.nan, index=out.index, dtype=float)
    lon = pd.Series(np.nan, index=out.index, dtype=float)
    if coord_col in out.columns:
        parsed = out[coord_col].fillna("").astype(str).str.extract(_COORD_PAIR)
        lat = pd.to_numeric(parsed[0], errors="coerce")
        lon = pd.to_numeric(parsed[1], errors="coerce")
    for col, values in ((lat_col, lat), (lon_col, lon)):
        if col in out.columns:
            existing = pd.to_numeric(out[col], errors="coerce")
            out[col] = existing.where(existing.notna(), values).astype(float)
        else:
            out[col] = values.astype(float)
    return out


COORD_STATUSES = ("ok", "missing", "unparseable", "centroid", "zero", "out_of_range")


def mark_coordinate_quality(
    df: pd.DataFrame,
    coord_source_column: str = 'coord_source',
    lat_col: str = 'lat',
    lon_col: str = 'lon',
    coord_col: str = 'coord',
) -> pd.DataFrame:
    """
    Mark coordinate quality without removing any records.

    Every record keeps its parsed lat/lon (see parse_coord_column); only
    region assignment and maps use is_geographic_quality.

    Added Columns
    -------------
    coord_status : str
        One of::

            ok            parseable, in range, not [0, 0], not a centroid
            missing       no coordinate given
            unparseable   a coordinate was given but could not be parsed
            centroid      coord_source says the coordinate is a country or province
                          centroid; lat/lon are kept, but the record is not assigned
                          a region or mapped
            zero          exactly [0, 0]
            out_of_range  latitude outside [-90, 90] or longitude outside [-180, 180]
    has_missing_coords, has_centroid_coords, has_zero_coords, has_invalid_coords : bool
        One flag per problem (has_missing_coords covers missing and unparseable)
    is_geographic_quality : bool
        coord_status == "ok"

    Parameters
    ----------
    df : pd.DataFrame
        Records with lat/lon columns (call parse_coord_column first)
    coord_source_column : str
        Name of the coordinate source column (default: 'coord_source')
    lat_col, lon_col, coord_col : str
        Column names

    Returns
    -------
    pd.DataFrame
        Copy of df with the columns above
    """
    df_marked = df.copy()
    idx = df_marked.index

    lat = pd.to_numeric(df_marked[lat_col], errors="coerce") if lat_col in df_marked else pd.Series(np.nan, index=idx)
    lon = pd.to_numeric(df_marked[lon_col], errors="coerce") if lon_col in df_marked else pd.Series(np.nan, index=idx)
    if coord_col in df_marked.columns:
        has_text = df_marked[coord_col].fillna("").astype(str).str.strip() != ""
    else:
        has_text = pd.Series(False, index=idx)
    parsed = lat.notna() & lon.notna()
    if coord_source_column in df_marked.columns:
        centroid = df_marked[coord_source_column].fillna("").astype(str).str.contains(
            "centroid", case=False)
    else:
        centroid = pd.Series(False, index=idx)

    # Precedence (later assignments win): missing/unparseable, out of range, zero, centroid
    status = pd.Series("ok", index=idx, dtype=object)
    status[~parsed & ~has_text] = "missing"
    status[~parsed & has_text] = "unparseable"
    status[parsed & ((lat.abs() > 90) | (lon.abs() > 180))] = "out_of_range"
    status[parsed & (lat == 0) & (lon == 0)] = "zero"
    status[parsed & centroid] = "centroid"

    df_marked['coord_status'] = status
    df_marked['has_missing_coords'] = status.isin(["missing", "unparseable"])
    df_marked['has_centroid_coords'] = status == "centroid"
    df_marked['has_zero_coords'] = status == "zero"
    df_marked['has_invalid_coords'] = status == "out_of_range"
    df_marked['is_geographic_quality'] = status == "ok"

    counts = status.value_counts()
    for s in COORD_STATUSES[1:]:
        if counts.get(s, 0):
            logger.info(f"Marked {counts[s]} samples with coord_status '{s}'")
    n_total = len(df_marked)
    n_geo_quality = int(df_marked['is_geographic_quality'].sum())
    pct_geo_quality = (n_geo_quality / n_total * 100) if n_total > 0 else 0
    logger.info(
        f"Geographic quality: {n_geo_quality}/{n_total} samples ({pct_geo_quality:.1f}%) "
        f"have parseable, non-zero, in-range, non-centroid coordinates"
    )

    return df_marked


def get_coordinate_quality_stats(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Calculate statistics about coordinate quality.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with coordinate columns

    Returns
    -------
    Dict[str, Any]
        Statistics about coordinate quality
    """
    stats = {}

    if 'latitude' in df.columns and 'longitude' in df.columns:
        stats['total_samples'] = len(df)
        stats['with_coordinates'] = df['latitude'].notna().sum()
        stats['missing_coordinates'] = df['latitude'].isna().sum()

        # Zero coordinates
        mask_zero = (df['latitude'] == 0) & (df['longitude'] == 0)
        stats['zero_coordinates'] = mask_zero.sum()

        # Centroid coordinates
        if 'coord_source' in df.columns:
            mask_centroid = df['coord_source'].str.contains(
                'centroid', case=False, na=False
            )
            stats['centroid_coordinates'] = mask_centroid.sum()

        # Coordinate precision (decimal places)
        if df['latitude'].notna().any():
            lat_precision = df['latitude'].dropna().apply(
                lambda x: len(str(x).split('.')[-1]) if '.' in str(x) else 0
            )
            stats['avg_coordinate_precision'] = lat_precision.mean()

    return stats


# ============================================================================
# Validation and Quality Control
# ============================================================================

def validate_metadata_for_analysis(df: pd.DataFrame) -> Tuple[bool, List[str]]:
    """
    Validate metadata DataFrame for downstream analysis.

    Checks for common issues that would prevent analysis from completing
    successfully.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata DataFrame

    Returns
    -------
    Tuple[bool, List[str]]
        (is_valid, list_of_issues)

    Examples
    --------
    >>> df = parse_bold_tsv("data.tsv")
    >>> is_valid, issues = validate_metadata_for_analysis(df)
    >>> if not is_valid:
    ...     for issue in issues:
    ...         print(f"Issue: {issue}")
    """
    issues = []

    # Check required columns
    required = ['processid', 'nuc']
    missing = [col for col in required if col not in df.columns]
    if missing:
        issues.append(f"Missing required columns: {missing}")

    # Check for data
    if len(df) == 0:
        issues.append("DataFrame is empty")
        return False, issues

    # Check processid uniqueness
    if 'processid' in df.columns:
        duplicates = df['processid'].duplicated().sum()
        if duplicates > 0:
            issues.append(f"Found {duplicates} duplicate processids")

    # Check sequences
    if 'nuc' in df.columns:
        empty_seqs = (df['nuc'].isna() | (df['nuc'] == '')).sum()
        if empty_seqs == len(df):
            issues.append("All sequences are empty")
        elif empty_seqs > 0:
            issues.append(f"{empty_seqs} samples have empty sequences")

    # Check coordinates if present
    if 'latitude' in df.columns and 'longitude' in df.columns:
        n_coords = df['latitude'].notna().sum()
        if n_coords == 0:
            issues.append("No valid coordinates found")

    is_valid = len(issues) == 0
    return is_valid, issues
