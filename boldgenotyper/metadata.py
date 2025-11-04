"""
TSV Parsing and Genotype Assignment

This module handles parsing of BOLD TSV files and assignment of samples to
their best-matching consensus genotype groups based on sequence similarity.

Key Responsibilities:
1. Parse BOLD TSV files (~86 columns) extracting critical fields:
   - processid: Unique sample identifier (REQUIRED)
   - nuc: COI nucleotide sequence (REQUIRED)
   - coord: Geographic coordinates in format [lat, lon]
   - coord_source: Indicates if coordinates are from centroid
   - country/ocean: Location metadata
   - Other taxonomic and collection metadata

2. Genotype Assignment:
   - Map each sample to its best-matching consensus sequence
   - Use global edit distance (edlib preferred, fallback to Levenshtein)
   - Require minimum identity threshold (default: 0.90)
   - Add consensus_group column to metadata
   - Generate diagnostics showing identity scores

3. Data Validation:
   - Check for required fields
   - Handle missing or malformed data gracefully
   - Provide informative warnings for data quality issues

Important Notes:
- The genotype assignment uses edit distance rather than alignment-based
  similarity for computational efficiency with large datasets.
- Samples below the minimum identity threshold are flagged but retained
  in the output for manual review.
- The processid field must be unique; duplicates will raise an error.

File Naming Convention:
- Input: {organism}_name.tsv (from BOLD)
- Output: {organism}_with_genotypes.tsv

Example Usage:
    >>> from boldgenotyper.metadata import parse_bold_tsv, assign_genotypes
    >>> df = parse_bold_tsv("Sphyrna_lewini.tsv")
    >>> df_with_genotypes = assign_genotypes(
    ...     tsv_path="Sphyrna_lewini.tsv",
    ...     consensus_fasta="Sphyrna_lewini_consensus.fasta",
    ...     output_path="Sphyrna_lewini_with_genotypes.tsv",
    ...     min_identity=0.90
    ... )

Author: Steph Smith (steph.smith@unc.edu)
"""

from typing import Dict, List, Optional, Tuple, Set, Any, Union
from pathlib import Path
import logging
import re
import pandas as pd
import numpy as np
from functools import partial
import multiprocessing as mp

# Try to import edlib for fast edit distance
try:
    import edlib
    EDLIB_AVAILABLE = True
except ImportError:
    EDLIB_AVAILABLE = False
    logger = logging.getLogger(__name__)
    logger.debug("edlib not available; using pure Python Levenshtein distance")

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
    Performs comprehensive validation and provides informative error messages.

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

    # Try reading with specified encoding
    try:
        df = pd.read_csv(
            path,
            sep='\t',
            encoding=encoding,
            dtype=str,  # Read all as strings initially
            low_memory=False
        )
    except UnicodeDecodeError:
        # Fallback to latin-1 encoding
        logger.warning(f"UTF-8 encoding failed, trying latin-1")
        df = pd.read_csv(
            path,
            sep='\t',
            encoding='latin-1',
            dtype=str,
            low_memory=False
        )

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
    for col in df.select_dtypes(include=['object']).columns:
        df[col] = df[col].str.strip() if df[col].dtype == 'object' else df[col]

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

def extract_coordinates(coord_string: Union[str, float]) -> Optional[Tuple[float, float]]:
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

        # Parse latitude and longitude
        lat = float(parts[0].strip())
        lon = float(parts[1].strip())

        # Validate ranges
        if not (-90 <= lat <= 90):
            logger.debug(f"Latitude out of range: {lat}")
            return None

        if not (-180 <= lon <= 180):
            logger.debug(f"Longitude out of range: {lon}")
            return None

        return (lat, lon)

    except (ValueError, AttributeError, TypeError) as e:
        logger.debug(f"Failed to parse coordinates '{coord_string}': {e}")
        return None


def parse_coordinates_column(df: pd.DataFrame, coord_column: str = 'coord') -> pd.DataFrame:
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

    Returns
    -------
    pd.DataFrame
        DataFrame with added latitude and longitude columns

    Examples
    --------
    >>> df = pd.DataFrame({'coord': ['[34.5, -76.2]', '[0, 0]', None]})
    >>> df_parsed = parse_coordinates_column(df)
    >>> print(df_parsed[['latitude', 'longitude']])
    """
    if coord_column not in df.columns:
        logger.warning(f"Coordinate column '{coord_column}' not found")
        df['latitude'] = np.nan
        df['longitude'] = np.nan
        return df

    logger.info(f"Parsing coordinates from '{coord_column}' column")

    # Apply coordinate extraction
    coords = df[coord_column].apply(extract_coordinates)

    # Split into separate columns
    df['latitude'] = coords.apply(lambda x: x[0] if x is not None else np.nan)
    df['longitude'] = coords.apply(lambda x: x[1] if x is not None else np.nan)

    # Log parsing statistics
    n_parsed = coords.notna().sum()
    n_failed = coords.isna().sum()
    logger.info(f"Parsed {n_parsed} coordinates successfully, {n_failed} failed")

    return df


def filter_by_coordinate_quality(
    df: pd.DataFrame,
    exclude_centroids: bool = True,
    exclude_zero_coords: bool = True,
    exclude_missing: bool = True,
    coord_source_column: str = 'coord_source',
) -> pd.DataFrame:
    """
    Filter samples based on coordinate quality.

    Implements critical filtering rules to exclude samples with low-quality
    or ambiguous coordinates. This is ESSENTIAL for accurate biogeographic
    analysis because many BOLD records contain country-level centroids rather
    than actual collection locations.

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
    This filtering is CRITICAL for accurate ocean basin assignment.
    Many BOLD records use country centroids instead of actual collection
    locations, which can lead to incorrect biogeographic conclusions.
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
# Edit Distance Calculation
# ============================================================================

def levenshtein_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Levenshtein (edit) distance between two sequences.

    Pure Python implementation for fallback when edlib is not available.
    Uses dynamic programming with space optimization.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    int
        Edit distance (number of insertions, deletions, substitutions)

    Examples
    --------
    >>> levenshtein_distance("ATCG", "ATCG")
    0
    >>> levenshtein_distance("ATCG", "ATGG")
    1

    Notes
    -----
    Based on implementation from consensus_group_to_metadata.py
    Time complexity: O(len(seq1) * len(seq2))
    Space complexity: O(min(len(seq1), len(seq2)))
    """
    # Ensure seq1 is shorter for memory efficiency
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1

    previous = list(range(len(seq1) + 1))

    for j, char_b in enumerate(seq2, start=1):
        current = [j]
        for i, char_a in enumerate(seq1, start=1):
            # Cost of operations
            insertion = previous[i] + 1
            deletion = current[i - 1] + 1
            substitution = previous[i - 1] + (char_a != char_b)
            current.append(min(insertion, deletion, substitution))
        previous = current

    return previous[-1]


def calculate_edit_distance(
    seq1: str,
    seq2: str,
    use_edlib: bool = True
) -> Tuple[int, float]:
    """
    Calculate edit distance and percent identity between two sequences.

    Uses edlib for speed if available, falls back to pure Python implementation.
    Calculates global alignment distance (Needleman-Wunsch).

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence
    use_edlib : bool
        Use edlib library if available (default: True)

    Returns
    -------
    Tuple[int, float]
        (edit_distance, percent_identity)
        - edit_distance: number of edits needed to transform seq1 to seq2
        - percent_identity: similarity score from 0-1

    Examples
    --------
    >>> dist, identity = calculate_edit_distance("ATCG", "ATCG")
    >>> print(dist, identity)
    0 1.0
    >>> dist, identity = calculate_edit_distance("ATCG", "ATGG")
    >>> print(dist, identity)
    1 0.75

    Notes
    -----
    Identity calculation: 1 - (edit_distance / max(len(seq1), len(seq2)))
    This differs from aligned sequence identity which only counts matching bases.
    """
    # Convert to uppercase
    s1 = seq1.upper()
    s2 = seq2.upper()

    # Calculate edit distance
    if use_edlib and EDLIB_AVAILABLE:
        # Use edlib for fast global alignment
        result = edlib.align(s1, s2, mode="NW", task="distance")
        distance = result["editDistance"]
    else:
        # Use pure Python implementation
        distance = levenshtein_distance(s1, s2)

    # Calculate percent identity
    max_len = max(len(s1), len(s2))
    if max_len == 0:
        identity = 1.0
    else:
        identity = 1.0 - (distance / max_len)

    return distance, identity


# ============================================================================
# Genotype Assignment
# ============================================================================

def assign_genotypes(
    tsv_path: Union[str, Path],
    consensus_fasta: Union[str, Path],
    output_path: Union[str, Path],
    raw_fasta: Optional[Union[str, Path]] = None,
    min_identity: float = 0.90,
    n_threads: int = 1,
    use_edlib: bool = True,
    diagnostics_path: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Assign each sample to its best-matching consensus genotype.

    Uses global edit distance to find the closest consensus sequence for each
    sample. Samples below the minimum identity threshold are flagged for review.

    Parameters
    ----------
    tsv_path : Union[str, Path]
        Path to BOLD TSV file with metadata
    consensus_fasta : Union[str, Path]
        Path to consensus sequences FASTA
    output_path : Union[str, Path]
        Path for output TSV with genotype assignments
    raw_fasta : Optional[Union[str, Path]]
        Path to raw sequences FASTA. If None, uses 'nuc' column from TSV
    min_identity : float
        Minimum sequence identity for assignment (default: 0.90)
    n_threads : int
        Number of parallel processes (default: 1)
    use_edlib : bool
        Use edlib for fast edit distance if available (default: True)
    diagnostics_path : Optional[Union[str, Path]]
        Path to save diagnostics CSV with identity scores

    Returns
    -------
    pd.DataFrame
        Metadata with added columns:
        - consensus_group: assigned genotype name
        - identity_score: best match identity (0-1)
        - runner_up_group: second-best match (for tie detection)
        - runner_up_identity: second-best identity

    Examples
    --------
    >>> df = assign_genotypes(
    ...     tsv_path="Sphyrna_lewini.tsv",
    ...     consensus_fasta="Sphyrna_lewini_consensus.fasta",
    ...     output_path="Sphyrna_lewini_with_genotypes.tsv",
    ...     min_identity=0.90,
    ...     n_threads=4
    ... )

    Notes
    -----
    Based on consensus_group_to_metadata.py from reference scripts.
    Uses multiprocessing for efficiency with large datasets.
    """
    logger.info("Starting genotype assignment")

    # Load metadata
    logger.info(f"Loading metadata from {tsv_path}")
    metadata = parse_bold_tsv(tsv_path, required_columns=['processid'])

    # Load consensus sequences
    logger.info(f"Loading consensus sequences from {consensus_fasta}")
    consensus_seqs = _load_consensus_sequences(consensus_fasta)
    logger.info(f"Loaded {len(consensus_seqs)} consensus sequences")

    # Get sequences for each processid
    if raw_fasta is not None:
        logger.info(f"Loading raw sequences from {raw_fasta}")
        pid_to_seq = _load_raw_sequences_from_fasta(raw_fasta)
    else:
        logger.info("Using sequences from 'nuc' column in TSV")
        pid_to_seq = _get_sequences_from_tsv(metadata)

    # Prepare tasks for parallel processing
    processids = metadata['processid'].tolist()
    tasks = [(pid, pid_to_seq.get(pid)) for pid in processids]

    logger.info(f"Assigning {len(tasks)} samples to genotypes using {n_threads} threads")

    # Create worker function with fixed parameters
    worker = partial(
        _assign_single_sample,
        consensus_seqs=consensus_seqs,
        min_identity=min_identity,
        use_edlib=use_edlib
    )

    # Run assignment (parallel or serial)
    if n_threads > 1:
        with mp.Pool(processes=n_threads) as pool:
            results = pool.map(worker, tasks)
    else:
        results = [worker(task) for task in tasks]

    # Collect results
    assignments = {r['processid']: r for r in results}

    # Add to metadata
    metadata['consensus_group'] = metadata['processid'].map(
        lambda pid: assignments[pid]['consensus_group']
    )
    metadata['identity_score'] = metadata['processid'].map(
        lambda pid: assignments[pid]['identity_score']
    )
    metadata['runner_up_group'] = metadata['processid'].map(
        lambda pid: assignments[pid]['runner_up_group']
    )
    metadata['runner_up_identity'] = metadata['processid'].map(
        lambda pid: assignments[pid]['runner_up_identity']
    )

    # Log assignment statistics
    _log_assignment_stats(metadata)

    # Save output
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved annotated metadata to {output_path}")

    # Save diagnostics if requested
    if diagnostics_path is not None:
        _save_diagnostics(results, diagnostics_path)

    return metadata


def _load_consensus_sequences(fasta_path: Union[str, Path]) -> List[Tuple[str, str]]:
    """Load consensus sequences from FASTA file."""
    from .utils import read_fasta
    records = read_fasta(fasta_path)
    # Extract group ID (first word of header)
    consensus_seqs = [(hdr.split()[0], seq) for hdr, seq in records]
    return consensus_seqs


def _load_raw_sequences_from_fasta(fasta_path: Union[str, Path]) -> Dict[str, str]:
    """Load raw sequences from FASTA file and extract processids."""
    from .utils import read_fasta, extract_processid_from_header

    records = read_fasta(fasta_path)
    pid_to_seq = {}

    for header, seq in records:
        pid = extract_processid_from_header(header)
        if pid:
            # Store first occurrence only
            pid_to_seq.setdefault(pid, seq)
        else:
            logger.debug(f"Could not extract processid from header: {header}")

    logger.info(f"Extracted {len(pid_to_seq)} sequences with processids from FASTA")
    return pid_to_seq


def _get_sequences_from_tsv(df: pd.DataFrame) -> Dict[str, str]:
    """Get sequences from 'nuc' column of metadata."""
    if 'nuc' not in df.columns:
        raise ValueError("Metadata must have 'nuc' column if raw FASTA not provided")

    pid_to_seq = {}
    for _, row in df.iterrows():
        pid = row['processid']
        seq = row.get('nuc', '')
        if pd.notna(seq) and seq:
            pid_to_seq[pid] = str(seq).upper()

    return pid_to_seq


def _assign_single_sample(
    task: Tuple[str, Optional[str]],
    consensus_seqs: List[Tuple[str, str]],
    min_identity: float,
    use_edlib: bool
) -> Dict[str, Any]:
    """
    Assign a single sample to best-matching consensus.

    Worker function for parallel processing.
    """
    processid, sequence = task

    # Handle missing sequence
    if not sequence or pd.isna(sequence):
        return {
            'processid': processid,
            'consensus_group': None,
            'identity_score': 0.0,
            'runner_up_group': None,
            'runner_up_identity': 0.0,
        }

    # Find best and second-best matches
    best_group = None
    best_identity = -1.0
    runner_up_group = None
    runner_up_identity = -1.0

    for group_id, consensus_seq in consensus_seqs:
        _, identity = calculate_edit_distance(sequence, consensus_seq, use_edlib)

        if identity > best_identity:
            # Shift best to runner-up
            runner_up_group = best_group
            runner_up_identity = best_identity
            # Update best
            best_group = group_id
            best_identity = identity
        elif identity > runner_up_identity:
            # Update runner-up
            runner_up_group = group_id
            runner_up_identity = identity

    # Check minimum identity threshold
    if best_identity < min_identity:
        best_group = None  # Flag as unassigned

    return {
        'processid': processid,
        'consensus_group': best_group,
        'identity_score': round(best_identity, 6),
        'runner_up_group': runner_up_group,
        'runner_up_identity': round(runner_up_identity, 6),
    }


def _log_assignment_stats(df: pd.DataFrame) -> None:
    """Log statistics about genotype assignments."""
    total = len(df)
    assigned = df['consensus_group'].notna().sum()
    unassigned = df['consensus_group'].isna().sum()

    logger.info(f"Assignment summary: {assigned}/{total} assigned, {unassigned} unassigned")

    if 'consensus_group' in df.columns and assigned > 0:
        genotype_counts = df['consensus_group'].value_counts()
        logger.info(f"Number of genotypes: {len(genotype_counts)}")
        logger.info(f"Genotype distribution:")
        for genotype, count in genotype_counts.head(10).items():
            logger.info(f"  {genotype}: {count}")

    # Check for ties (small difference between best and runner-up)
    if 'identity_score' in df.columns and 'runner_up_identity' in df.columns:
        df_assigned = df[df['consensus_group'].notna()].copy()
        if len(df_assigned) > 0:
            df_assigned['identity_diff'] = (
                df_assigned['identity_score'] - df_assigned['runner_up_identity']
            )
            ties = (df_assigned['identity_diff'] < 0.01).sum()
            if ties > 0:
                logger.warning(
                    f"Found {ties} samples with ambiguous assignments "
                    f"(identity difference < 0.01)"
                )


def _save_diagnostics(
    results: List[Dict[str, Any]],
    diagnostics_path: Union[str, Path]
) -> None:
    """Save diagnostics CSV with assignment details."""
    df_diag = pd.DataFrame(results)
    path = Path(diagnostics_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df_diag.to_csv(path, index=False)
    logger.info(f"Saved diagnostics to {path}")


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
