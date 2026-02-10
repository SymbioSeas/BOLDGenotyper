"""
Metadata Analysis Module

This module provides comprehensive analysis of haplotype associations with
non-geographic metadata fields. It addresses two key use cases:
1. Datasets where geographic coordinates are sparse but other metadata is informative
2. Studies incorporating additional metadata variables of interest beyond geography

Categorical Fields Analyzed:
- sex: Specimen sex
- life_stage: Developmental stage
- reproduction: Reproductive status
- country/ocean, country_iso: Geographic jurisdiction
- province/state: Province or state
- realm, biome, ecoregion: Biogeographic classifications
- habitat: Habitat description
- geoid: BOLD geographic identifier

Temporal Field:
- collection_date_start: Collection date (various formats supported)

Philosophy:
Values are presented as-is without standardization. Users know their data best
and should decide how to standardize values. Optional sex normalization is
available via the --normalize-sex flag.

Author: Steph Smith (steph.smith@unc.edu)
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
import warnings

logger = logging.getLogger(__name__)

# Default metadata fields to analyze
DEFAULT_CATEGORICAL_FIELDS = [
    'sex',
    'life_stage',
    'reproduction',
    'country/ocean',
    'country_iso',
    'province/state',
    'realm',
    'biome',
    'ecoregion',
    'habitat',
    'geoid',
]

DATE_FIELD = 'collection_date_start'


def _save_multi_format(fig, output_path: Path, dpi: int = 150):
    """
    Save a figure in PNG, PDF, and SVG formats to their correct directories.

    Routes each format to the correct subdirectory under visualization/metadata/:
    - PDF → visualization/metadata/pdf/
    - SVG → visualization/metadata/svg/
    - PNG → visualization/metadata/png/

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save (unused; saves via plt.savefig for current figure)
    output_path : Path
        Base output path (typically in the pdf/ subdirectory)
    dpi : int
        DPI for PNG output
    """
    output_path = Path(output_path)
    base_name = output_path.stem

    # Determine parent directories for each format
    # output_path is typically in .../metadata/pdf/basename
    # We need to route to .../metadata/svg/ and .../metadata/png/
    parent = output_path.parent
    metadata_base = parent.parent  # .../metadata/

    fmt_dirs = {
        'pdf': parent,  # Keep PDF in the original directory
        'svg': metadata_base / 'svg',
        'png': metadata_base / 'png',
    }

    for fmt, d in fmt_dirs.items():
        d.mkdir(parents=True, exist_ok=True)
        save_path = d / f"{base_name}.{fmt}"
        if fmt == 'png':
            plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        else:
            plt.savefig(save_path, bbox_inches='tight')


# =============================================================================
# Coverage and Variant Detection Functions
# =============================================================================

def analyze_metadata_coverage(
    df: pd.DataFrame,
    fields: List[str]
) -> Dict[str, Any]:
    """
    Calculate coverage statistics for each metadata field.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sample metadata
    fields : List[str]
        List of metadata field names to analyze

    Returns
    -------
    Dict[str, Any]
        Dictionary with coverage statistics for each field including:
        - n_total: Total samples
        - n_with_value: Samples with non-null value
        - pct_coverage: Percentage coverage
        - n_unique: Number of unique values
        - top_values: Top 5 values with counts
        - warning: Any data quality warnings
    """
    total_samples = len(df)
    coverage_stats = {
        'total_samples': total_samples,
        'fields': {}
    }

    for field in fields:
        field_stats = {
            'n_with_value': 0,
            'pct_coverage': 0.0,
            'n_unique': 0,
            'unique_values': [],
            'value_counts': {},
            'top_values': [],
            'warning': None
        }

        if field not in df.columns:
            field_stats['warning'] = f"Field '{field}' not found in data"
            coverage_stats['fields'][field] = field_stats
            continue

        # Get non-null values
        values = df[field].dropna()
        values = values[values != '']  # Also exclude empty strings

        field_stats['n_with_value'] = len(values)
        field_stats['pct_coverage'] = (len(values) / total_samples * 100) if total_samples > 0 else 0.0
        field_stats['n_unique'] = values.nunique()
        field_stats['unique_values'] = sorted(values.unique().tolist())

        # Get value counts
        value_counts = values.value_counts()
        field_stats['value_counts'] = value_counts.to_dict()

        # Top 5 values as formatted strings
        top_5 = value_counts.head(5)
        field_stats['top_values'] = [
            f"{val}: {count}" for val, count in top_5.items()
        ]

        # Check for potential variants
        variants = detect_value_variants(values)
        if variants:
            field_stats['warning'] = f"Multiple value variants detected: {', '.join(variants)}"

        coverage_stats['fields'][field] = field_stats

    return coverage_stats


def detect_value_variants(values: pd.Series) -> List[str]:
    """
    Detect potential variant values (e.g., 'M' vs 'male').

    Checks for common variations like:
    - Case differences (M/Male, F/Female)
    - Abbreviations vs full names
    - Trailing whitespace

    Parameters
    ----------
    values : pd.Series
        Series of string values to check

    Returns
    -------
    List[str]
        List of variant group descriptions (e.g., "M/Male/male")
    """
    if values.empty:
        return []

    # Convert to strings and clean
    str_values = values.astype(str).str.strip()
    unique_vals = str_values.unique()

    variants_found = []
    checked = set()

    for val in unique_vals:
        if val.lower() in checked:
            continue

        # Find all case variants
        val_lower = val.lower()
        case_variants = [v for v in unique_vals if v.lower() == val_lower and v != val]

        if case_variants:
            all_variants = [val] + case_variants
            variants_found.append('/'.join(sorted(all_variants)))
            checked.add(val_lower)

    # Check for common sex variants specifically
    sex_variants = {
        'm': ['M', 'm', 'Male', 'male', 'MALE'],
        'f': ['F', 'f', 'Female', 'female', 'FEMALE']
    }

    for base, variants in sex_variants.items():
        present = [v for v in variants if v in unique_vals]
        if len(present) > 1 and f"{'/'.join(sorted(present))}" not in variants_found:
            variants_found.append('/'.join(sorted(present)))

    return variants_found


def normalize_sex_values(
    df: pd.DataFrame,
    col: str = 'sex'
) -> pd.DataFrame:
    """
    Normalize sex values to 'Male'/'Female' standard.

    Mapping:
    - M, Male, m, male, MALE -> 'Male'
    - F, Female, f, female, FEMALE -> 'Female'
    - Other values preserved as-is

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sex column
    col : str, optional
        Column name containing sex values (default: 'sex')

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized sex values
    """
    if col not in df.columns:
        logger.warning(f"Column '{col}' not found, skipping sex normalization")
        return df

    df = df.copy()

    # Define mappings
    male_values = {'m', 'male', 'M', 'Male', 'MALE'}
    female_values = {'f', 'female', 'F', 'Female', 'FEMALE'}

    def normalize(val):
        if pd.isna(val):
            return val
        val_str = str(val).strip()
        if val_str in male_values:
            return 'Male'
        elif val_str in female_values:
            return 'Female'
        return val_str

    df[col] = df[col].apply(normalize)

    # Log normalization summary
    n_male = (df[col] == 'Male').sum()
    n_female = (df[col] == 'Female').sum()
    n_other = df[col].notna().sum() - n_male - n_female
    logger.info(f"Sex normalization: Male={n_male}, Female={n_female}, Other={n_other}")

    return df


# =============================================================================
# Categorical Analysis Functions
# =============================================================================

def analyze_categorical_field(
    df: pd.DataFrame,
    field: str,
    haplotype_col: str = 'haplotype_sp'
) -> pd.DataFrame:
    """
    Analyze haplotype distribution across categorical values.

    For each haplotype × category combination, calculates:
    - n_samples: Number of samples
    - pct_of_haplotype: Proportion of this haplotype in this category
    - pct_of_value: Proportion of this category that is this haplotype

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments and metadata
    field : str
        Metadata field to analyze
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')

    Returns
    -------
    pd.DataFrame
        Analysis results with columns:
        - haplotype_sp
        - value
        - n_samples
        - pct_of_haplotype
        - pct_of_value
    """
    if field not in df.columns:
        logger.warning(f"Field '{field}' not found in data")
        return pd.DataFrame()

    if haplotype_col not in df.columns:
        logger.warning(f"Haplotype column '{haplotype_col}' not found in data")
        return pd.DataFrame()

    # Filter to rows with both haplotype and field value
    valid_df = df[[haplotype_col, field]].dropna()
    valid_df = valid_df[valid_df[field] != '']

    if valid_df.empty:
        return pd.DataFrame()

    # Calculate cross-tabulation
    cross_tab = pd.crosstab(valid_df[haplotype_col], valid_df[field])

    # Calculate totals
    haplotype_totals = cross_tab.sum(axis=1)
    value_totals = cross_tab.sum(axis=0)

    # Build results
    results = []
    for haplotype in cross_tab.index:
        for value in cross_tab.columns:
            n_samples = cross_tab.loc[haplotype, value]
            if n_samples > 0:
                pct_of_haplotype = (n_samples / haplotype_totals[haplotype] * 100) if haplotype_totals[haplotype] > 0 else 0
                pct_of_value = (n_samples / value_totals[value] * 100) if value_totals[value] > 0 else 0

                results.append({
                    'haplotype_sp': haplotype,
                    'value': value,
                    'n_samples': n_samples,
                    'pct_of_haplotype': round(pct_of_haplotype, 1),
                    'pct_of_value': round(pct_of_value, 1)
                })

    results_df = pd.DataFrame(results)

    # Sort by haplotype and n_samples
    if not results_df.empty:
        results_df = results_df.sort_values(
            ['haplotype_sp', 'n_samples'],
            ascending=[True, False]
        )

    return results_df


def test_haplotype_association(
    df: pd.DataFrame,
    field: str,
    haplotype_col: str = 'haplotype_sp'
) -> Dict[str, Any]:
    """
    Perform chi-square test for haplotype × field association.

    Tests whether haplotype distribution is independent of the
    categorical metadata field.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments and metadata
    field : str
        Metadata field to test
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')

    Returns
    -------
    Dict[str, Any]
        Test results including:
        - field: Field name
        - test_type: 'chi_square' or 'fisher_exact'
        - statistic: Test statistic
        - p_value: P-value
        - degrees_of_freedom: Degrees of freedom (chi-square only)
        - n_samples: Sample size used in test
        - n_haplotypes: Number of haplotypes
        - n_categories: Number of category values
        - warning: Any warnings about test validity
    """
    result = {
        'field': field,
        'test_type': None,
        'statistic': np.nan,
        'p_value': np.nan,
        'degrees_of_freedom': np.nan,
        'n_samples': 0,
        'n_haplotypes': 0,
        'n_categories': 0,
        'warning': None
    }

    if field not in df.columns or haplotype_col not in df.columns:
        result['warning'] = 'Required column(s) not found'
        return result

    # Filter to valid rows
    valid_df = df[[haplotype_col, field]].dropna()
    valid_df = valid_df[valid_df[field] != '']

    if valid_df.empty:
        result['warning'] = 'No valid data for test'
        return result

    # Create contingency table
    contingency_table = pd.crosstab(valid_df[haplotype_col], valid_df[field])

    result['n_samples'] = len(valid_df)
    result['n_haplotypes'] = len(contingency_table.index)
    result['n_categories'] = len(contingency_table.columns)

    # Check minimum requirements
    if result['n_haplotypes'] < 2 or result['n_categories'] < 2:
        result['warning'] = 'Insufficient categories for test (need at least 2×2)'
        return result

    # Perform chi-square test
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            chi2, p_value, dof, expected = stats.chi2_contingency(contingency_table)

        result['test_type'] = 'chi_square'
        result['statistic'] = chi2
        result['p_value'] = p_value
        result['degrees_of_freedom'] = dof

        # Check expected frequencies
        if (expected < 5).sum() > expected.size * 0.2:
            result['warning'] = 'More than 20% of expected frequencies < 5; results may be unreliable'

    except Exception as e:
        result['warning'] = f'Test failed: {str(e)}'

    return result


# =============================================================================
# Temporal Analysis Functions
# =============================================================================

def parse_collection_dates(
    df: pd.DataFrame,
    date_col: str = 'collection_date_start'
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Parse collection dates and add year/month/decade columns.

    Handles multiple date formats including:
    - '12/14/98', '1998-12-14', '14-Dec-1998'
    - ISO formats, various separators

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with collection date column
    date_col : str, optional
        Column containing date strings (default: 'collection_date_start')

    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        - DataFrame with added columns: parsed_date, collection_year,
          collection_month, collection_decade
        - Statistics about parsing success
    """
    df = df.copy()
    parse_stats = {
        'n_total': len(df),
        'n_with_date': 0,
        'n_parsed': 0,
        'pct_parsed': 0.0,
        'date_range': None,
        'n_parse_failures': 0
    }

    if date_col not in df.columns:
        logger.warning(f"Date column '{date_col}' not found")
        df['parsed_date'] = pd.NaT
        df['collection_year'] = np.nan
        df['collection_month'] = np.nan
        df['collection_decade'] = np.nan
        return df, parse_stats

    # Count non-null dates
    non_null = df[date_col].notna() & (df[date_col] != '')
    parse_stats['n_with_date'] = non_null.sum()

    # Parse dates
    df['parsed_date'] = pd.to_datetime(
        df[date_col],
        errors='coerce'
    )

    # Count successful parses
    parsed_mask = df['parsed_date'].notna()
    parse_stats['n_parsed'] = parsed_mask.sum()
    parse_stats['pct_parsed'] = (
        parse_stats['n_parsed'] / parse_stats['n_with_date'] * 100
        if parse_stats['n_with_date'] > 0 else 0.0
    )
    parse_stats['n_parse_failures'] = parse_stats['n_with_date'] - parse_stats['n_parsed']

    # Extract year, month, decade
    df['collection_year'] = df['parsed_date'].dt.year
    df['collection_month'] = df['parsed_date'].dt.month
    df['collection_decade'] = (df['collection_year'] // 10) * 10

    # Get date range
    if parsed_mask.any():
        min_date = df.loc[parsed_mask, 'parsed_date'].min()
        max_date = df.loc[parsed_mask, 'parsed_date'].max()
        parse_stats['date_range'] = {
            'min': str(min_date.date()),
            'max': str(max_date.date())
        }

    logger.info(f"Date parsing: {parse_stats['n_parsed']}/{parse_stats['n_with_date']} dates parsed successfully ({parse_stats['pct_parsed']:.1f}%)")

    return df, parse_stats


def calculate_haplotype_emergence(
    df: pd.DataFrame,
    haplotype_col: str = 'haplotype_sp'
) -> pd.DataFrame:
    """
    Calculate first and last detection dates per haplotype.

    Tracks both emergence (first detection) and potential loss (last detection)
    of haplotypes over time.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed dates and haplotype assignments
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - haplotype_sp
        - first_collection_date
        - first_collection_year
        - first_location (country/ocean if available)
        - last_collection_date
        - last_collection_year
        - last_location (country/ocean if available)
        - years_observed (span from first to last detection)
        - n_total_samples
    """
    if 'parsed_date' not in df.columns:
        logger.warning("No parsed_date column. Run parse_collection_dates first.")
        return pd.DataFrame()

    if haplotype_col not in df.columns:
        logger.warning(f"Haplotype column '{haplotype_col}' not found")
        return pd.DataFrame()

    # Filter to valid dates and haplotypes
    valid_df = df[df['parsed_date'].notna() & df[haplotype_col].notna()].copy()

    if valid_df.empty:
        return pd.DataFrame()

    # Find first and last detection per haplotype
    results = []
    for haplotype, group in valid_df.groupby(haplotype_col):
        first_idx = group['parsed_date'].idxmin()
        first_row = group.loc[first_idx]
        last_idx = group['parsed_date'].idxmax()
        last_row = group.loc[last_idx]

        result = {
            'haplotype_sp': haplotype,
            'first_collection_date': first_row['parsed_date'].strftime('%Y-%m-%d'),
            'first_collection_year': int(first_row['collection_year']),
            'last_collection_date': last_row['parsed_date'].strftime('%Y-%m-%d'),
            'last_collection_year': int(last_row['collection_year']),
            'years_observed': int(last_row['collection_year']) - int(first_row['collection_year']),
            'n_total_samples': len(group)
        }

        # Add first location if available
        if 'country/ocean' in df.columns and pd.notna(first_row.get('country/ocean')):
            result['first_location'] = first_row['country/ocean']
        else:
            result['first_location'] = 'Unknown'

        # Add last location if available
        if 'country/ocean' in df.columns and pd.notna(last_row.get('country/ocean')):
            result['last_location'] = last_row['country/ocean']
        else:
            result['last_location'] = 'Unknown'

        results.append(result)

    emergence_df = pd.DataFrame(results)

    # Sort by first detection date
    emergence_df = emergence_df.sort_values('first_collection_date')

    return emergence_df


def analyze_temporal_distribution(
    df: pd.DataFrame,
    haplotype_col: str = 'haplotype_sp',
    time_unit: str = 'year'
) -> pd.DataFrame:
    """
    Analyze haplotype temporal distribution.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed dates and haplotype assignments
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    time_unit : str, optional
        Time unit for aggregation: 'year', 'month', or 'decade' (default: 'year')

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - time_period
        - haplotype_sp
        - n_samples
        - pct_of_period
    """
    time_col_map = {
        'year': 'collection_year',
        'month': 'collection_month',
        'decade': 'collection_decade'
    }

    if time_unit not in time_col_map:
        raise ValueError(f"time_unit must be one of {list(time_col_map.keys())}")

    time_col = time_col_map[time_unit]

    if time_col not in df.columns:
        logger.warning(f"Time column '{time_col}' not found. Run parse_collection_dates first.")
        return pd.DataFrame()

    if haplotype_col not in df.columns:
        logger.warning(f"Haplotype column '{haplotype_col}' not found")
        return pd.DataFrame()

    # Filter to valid data
    valid_df = df[df[time_col].notna() & df[haplotype_col].notna()].copy()

    if valid_df.empty:
        return pd.DataFrame()

    # Group by time period and haplotype
    counts = valid_df.groupby([time_col, haplotype_col]).size().reset_index(name='n_samples')

    # Calculate period totals
    period_totals = counts.groupby(time_col)['n_samples'].transform('sum')
    counts['pct_of_period'] = (counts['n_samples'] / period_totals * 100).round(1)

    # Rename for clarity
    counts = counts.rename(columns={time_col: 'time_period'})

    # Sort by time period and sample count
    counts = counts.sort_values(['time_period', 'n_samples'], ascending=[True, False])

    return counts


# =============================================================================
# Visualization Functions
# =============================================================================

def get_haplotype_colors(
    df: pd.DataFrame,
    haplotype_col: str = 'haplotype_sp'
) -> Dict[str, str]:
    """
    Get consistent color mapping for haplotypes ordered by abundance.

    Uses the same color scheme as visualization.py for consistency.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')

    Returns
    -------
    Dict[str, str]
        Mapping of haplotype names to hex color codes
    """
    # Reference colors (matches visualization.py)
    REFERENCE_COLORS = ['#8545C1', '#10B3A5', '#FFB031', '#3874F4', '#D975C7', '#132C54']

    if haplotype_col not in df.columns:
        return {}

    # Order by abundance
    haplotypes = df[haplotype_col].dropna().value_counts().index.tolist()
    n_haplotypes = len(haplotypes)

    if n_haplotypes == 0:
        return {}

    # Build color list
    if n_haplotypes <= len(REFERENCE_COLORS):
        colors = REFERENCE_COLORS[:n_haplotypes]
    else:
        # Extend with seaborn husl palette
        try:
            extra_colors = sns.color_palette("husl", n_haplotypes).as_hex()
        except Exception:
            extra_colors = list(mcolors.TABLEAU_COLORS.values())
        colors = REFERENCE_COLORS + [c for c in extra_colors if c not in REFERENCE_COLORS]
        colors = colors[:n_haplotypes]

    return {h: colors[i] for i, h in enumerate(haplotypes)}


def plot_metadata_coverage(
    coverage_stats: Dict[str, Any],
    output_path: Path,
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300
) -> None:
    """
    Plot horizontal bar chart of metadata coverage.

    Parameters
    ----------
    coverage_stats : Dict[str, Any]
        Coverage statistics from analyze_metadata_coverage()
    output_path : Path
        Path to save the plot
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 10x6)
    dpi : int, optional
        Resolution for output (default: 300)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fields = list(coverage_stats['fields'].keys())
    coverages = [coverage_stats['fields'][f]['pct_coverage'] for f in fields]

    # Sort by coverage
    sorted_data = sorted(zip(fields, coverages), key=lambda x: x[1], reverse=True)
    fields, coverages = zip(*sorted_data) if sorted_data else ([], [])

    fig, ax = plt.subplots(figsize=figsize)

    # Color bars by coverage level
    colors = ['#27ae60' if c >= 50 else '#f39c12' if c >= 20 else '#e74c3c' for c in coverages]

    bars = ax.barh(range(len(fields)), coverages, color=colors, edgecolor='black', linewidth=0.5)

    ax.set_yticks(range(len(fields)))
    ax.set_yticklabels(fields)
    ax.set_xlabel('Coverage (%)')
    ax.set_xlim(0, 100)
    ax.set_title('Metadata Field Coverage', fontsize=12, fontweight='bold')

    # Add value labels
    for i, (bar, cov) in enumerate(zip(bars, coverages)):
        ax.text(cov + 1, i, f'{cov:.1f}%', va='center', fontsize=9)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#27ae60', label='High (≥50%)'),
        Patch(facecolor='#f39c12', label='Medium (20-50%)'),
        Patch(facecolor='#e74c3c', label='Low (<20%)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', frameon=False)

    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved metadata coverage plot: {output_path}")


def plot_categorical_by_haplotype(
    df: pd.DataFrame,
    field: str,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp',
    figsize: Tuple[int, int] = (12, 6),
    dpi: int = 300,
    max_categories: int = 20
) -> None:
    """
    Plot stacked bar chart for categorical field by haplotype.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments and metadata
    field : str
        Metadata field to plot
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 12x6)
    dpi : int, optional
        Resolution for output (default: 300)
    max_categories : int, optional
        Maximum categories to display (default: 20)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if field not in df.columns or haplotype_col not in df.columns:
        logger.warning(f"Required columns not found for plotting {field}")
        return

    # Filter to valid data
    valid_df = df[[haplotype_col, field]].dropna()
    valid_df = valid_df[valid_df[field] != '']

    if valid_df.empty:
        logger.warning(f"No valid data for field '{field}'")
        return

    # Get color map if not provided
    if color_map is None:
        color_map = get_haplotype_colors(df, haplotype_col)

    # Create cross-tabulation
    cross_tab = pd.crosstab(valid_df[field], valid_df[haplotype_col], normalize='index')

    # Limit categories if needed
    if len(cross_tab) > max_categories:
        # Keep top categories by total count
        totals = valid_df[field].value_counts()
        top_cats = totals.head(max_categories).index
        cross_tab = cross_tab.loc[cross_tab.index.isin(top_cats)]
        logger.info(f"Limiting to top {max_categories} categories for {field}")

    # Order haplotypes by abundance
    haplotype_order = df[haplotype_col].dropna().value_counts().index.tolist()
    cross_tab = cross_tab[[h for h in haplotype_order if h in cross_tab.columns]]

    fig, ax = plt.subplots(figsize=figsize)

    # Plot stacked bars
    bottom = np.zeros(len(cross_tab))
    x_positions = range(len(cross_tab))

    for haplotype in cross_tab.columns:
        values = cross_tab[haplotype].values
        color = color_map.get(haplotype, '#999999')
        ax.bar(x_positions, values, bottom=bottom, label=haplotype,
               color=color, edgecolor='black', linewidth=0.5)
        bottom += values

    ax.set_xticks(x_positions)
    ax.set_xticklabels(cross_tab.index, rotation=45, ha='right')
    ax.set_ylabel('Relative Abundance')
    ax.set_xlabel(field.replace('_', ' ').title())
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{int(v*100)}%'))

    # Format title
    field_display = field.replace('_', ' ').replace('/', '/\n').title()
    ax.set_title(f'Haplotype Distribution by {field_display}', fontsize=12, fontweight='bold')

    # Legend
    ax.legend(title='Haplotype', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)

    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved categorical plot: {output_path}")


def plot_metadata_heatmap(
    df: pd.DataFrame,
    fields: List[str],
    output_path: Path,
    haplotype_col: str = 'haplotype_sp',
    figsize: Tuple[int, int] = (14, 10),
    dpi: int = 300
) -> None:
    """
    Plot haplotype × metadata heatmap.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments and metadata
    fields : List[str]
        Metadata fields to include
    output_path : Path
        Path to save the plot
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 14x10)
    dpi : int, optional
        Resolution for output (default: 300)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if haplotype_col not in df.columns:
        logger.warning(f"Haplotype column '{haplotype_col}' not found")
        return

    # Filter fields that exist and have data
    available_fields = [f for f in fields if f in df.columns and df[f].notna().any()]

    if not available_fields:
        logger.warning("No fields with data for heatmap")
        return

    # Build matrix of coverage per haplotype
    haplotypes = df[haplotype_col].dropna().value_counts().index.tolist()

    # Calculate coverage matrix
    coverage_matrix = []
    for haplotype in haplotypes:
        hap_df = df[df[haplotype_col] == haplotype]
        row = []
        for field in available_fields:
            coverage = hap_df[field].notna().sum() / len(hap_df) * 100 if len(hap_df) > 0 else 0
            row.append(coverage)
        coverage_matrix.append(row)

    coverage_df = pd.DataFrame(
        coverage_matrix,
        index=haplotypes,
        columns=[f.replace('_', ' ').title() for f in available_fields]
    )

    fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    im = ax.imshow(coverage_df.values, cmap='YlGnBu', aspect='auto')

    # Labels
    ax.set_xticks(range(len(coverage_df.columns)))
    ax.set_xticklabels(coverage_df.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(coverage_df.index)))
    ax.set_yticklabels(coverage_df.index)

    ax.set_xlabel('Metadata Field')
    ax.set_ylabel('Haplotype')
    ax.set_title('Metadata Coverage by Haplotype', fontsize=12, fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Coverage (%)')

    # Add value annotations
    for i in range(len(coverage_df.index)):
        for j in range(len(coverage_df.columns)):
            val = coverage_df.iloc[i, j]
            color = 'white' if val > 50 else 'black'
            ax.text(j, i, f'{val:.0f}', ha='center', va='center',
                   color=color, fontsize=8)

    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved metadata heatmap: {output_path}")


def plot_haplotype_emergence_timeline(
    emergence_df: pd.DataFrame,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (14, 8),
    dpi: int = 300
) -> None:
    """
    Plot timeline of haplotype first and last detections.

    Shows both emergence (first detection) and potential loss (last detection)
    of haplotypes over time, with horizontal bars indicating the observation span.

    Parameters
    ----------
    emergence_df : pd.DataFrame
        Emergence data from calculate_haplotype_emergence()
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 14x8)
    dpi : int, optional
        Resolution for output (default: 300)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if emergence_df.empty:
        logger.warning("No emergence data for timeline plot")
        return

    fig, ax = plt.subplots(figsize=figsize)

    # Convert dates
    emergence_df = emergence_df.copy()
    emergence_df['first_date'] = pd.to_datetime(emergence_df['first_collection_date'])
    if 'last_collection_date' in emergence_df.columns:
        emergence_df['last_date'] = pd.to_datetime(emergence_df['last_collection_date'])
    else:
        emergence_df['last_date'] = emergence_df['first_date']

    emergence_df = emergence_df.sort_values('first_date')

    # Assign y positions (0 to n-1)
    y_positions = list(range(len(emergence_df)))

    # Get colors
    if color_map is None:
        color_map = {}

    # Plot horizontal bars from first to last detection
    for i, (idx, row) in enumerate(emergence_df.iterrows()):
        color = color_map.get(row['haplotype_sp'], '#3874F4')

        # Draw horizontal line from first to last detection
        ax.hlines(
            y=i,
            xmin=row['first_date'],
            xmax=row['last_date'],
            colors=color,
            linewidths=3,
            alpha=0.6
        )

        # First detection marker (circle)
        ax.scatter(
            row['first_date'], i,
            c=[color],
            s=row['n_total_samples'] * 3 + 50,
            marker='o',
            alpha=0.9,
            edgecolors='black',
            linewidths=0.5,
            zorder=10
        )

        # Last detection marker (diamond) - only if different from first
        if row['last_date'] > row['first_date']:
            ax.scatter(
                row['last_date'], i,
                c=[color],
                s=row['n_total_samples'] * 3 + 50,
                marker='D',
                alpha=0.9,
                edgecolors='black',
                linewidths=0.5,
                zorder=10
            )

    # Labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(emergence_df['haplotype_sp'])
    ax.set_xlabel('Collection Date')
    ax.set_ylabel('Haplotype')
    ax.set_title('Haplotype Detection Timeline\n(First ● and Last ◆ Detection)', fontsize=12, fontweight='bold')

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=8, label='First Detection'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='gray',
               markersize=8, label='Last Detection'),
        Line2D([0], [0], color='gray', linewidth=3, alpha=0.6, label='Observation Span'),
    ]

    # Add size reference
    sizes = [10, 50, 100]
    for s in sizes:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgray',
                   markersize=np.sqrt(s*3+50)/2, label=f'n={s}')
        )

    ax.legend(handles=legend_elements, loc='lower right', frameon=True,
              fancybox=True, framealpha=0.9, title='Legend')

    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved emergence timeline: {output_path}")


def plot_temporal_composition(
    temporal_df: pd.DataFrame,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (12, 6),
    dpi: int = 300
) -> None:
    """
    Plot stacked area chart of composition over time.

    Parameters
    ----------
    temporal_df : pd.DataFrame
        Temporal data from analyze_temporal_distribution()
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 12x6)
    dpi : int, optional
        Resolution for output (default: 300)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if temporal_df.empty:
        logger.warning("No temporal data for composition plot")
        return

    # Pivot for stacked area
    pivot_df = temporal_df.pivot(
        index='time_period',
        columns='haplotype_sp',
        values='pct_of_period'
    ).fillna(0)

    # Sort time periods
    pivot_df = pivot_df.sort_index()

    # Get colors
    if color_map is None:
        color_map = {}

    colors = [color_map.get(h, '#999999') for h in pivot_df.columns]

    fig, ax = plt.subplots(figsize=figsize)

    # Plot stacked area
    ax.stackplot(pivot_df.index, pivot_df.values.T, labels=pivot_df.columns,
                 colors=colors, alpha=0.8)

    ax.set_xlabel('Time Period')
    ax.set_ylabel('Relative Abundance (%)')
    ax.set_title('Haplotype Composition Over Time', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 100)

    # Legend
    ax.legend(title='Haplotype', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)

    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved temporal composition: {output_path}")


def plot_temporal_composition_species_faceted(
    df: pd.DataFrame,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp',
    species_col: str = 'assigned_sp',
    figsize_per_facet: Tuple[int, int] = (10, 4),
    dpi: int = 300,
    max_species: int = 12
) -> None:
    """
    Plot species-faceted stacked area charts of haplotype composition over time.

    Creates a separate subplot for each species showing how its haplotype
    composition changes over time.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed dates, haplotype assignments, and species assignments
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    species_col : str, optional
        Column containing species assignments (default: 'assigned_sp')
    figsize_per_facet : Tuple[int, int], optional
        Figure size per facet in inches (default: 10x4)
    dpi : int, optional
        Resolution for output (default: 300)
    max_species : int, optional
        Maximum number of species to plot (default: 12)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if 'collection_year' not in df.columns:
        logger.warning("No collection_year column. Run parse_collection_dates first.")
        return

    if haplotype_col not in df.columns or species_col not in df.columns:
        logger.warning(f"Required columns not found for species-faceted temporal plot")
        return

    # Filter to valid data
    valid_df = df[
        df['collection_year'].notna() &
        df[haplotype_col].notna() &
        df[species_col].notna()
    ].copy()

    if valid_df.empty:
        logger.warning("No valid data for species-faceted temporal plot")
        return

    # Get species ordered by sample count
    species_counts = valid_df[species_col].value_counts()
    species_list = species_counts.head(max_species).index.tolist()

    if len(species_list) == 0:
        return

    # Calculate grid dimensions
    n_species = len(species_list)
    n_cols = min(2, n_species)
    n_rows = (n_species + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_facet[0] * n_cols, figsize_per_facet[1] * n_rows),
        squeeze=False
    )

    # Get colors
    if color_map is None:
        color_map = get_haplotype_colors(valid_df, haplotype_col)

    for idx, species in enumerate(species_list):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        # Filter to this species
        species_df = valid_df[valid_df[species_col] == species]

        # Calculate temporal distribution for this species
        counts = species_df.groupby(['collection_year', haplotype_col]).size().reset_index(name='n_samples')

        if counts.empty:
            ax.set_title(species, fontsize=10, fontweight='bold')
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            continue

        # Calculate percentages
        year_totals = counts.groupby('collection_year')['n_samples'].transform('sum')
        counts['pct'] = (counts['n_samples'] / year_totals * 100)

        # Pivot for stacked area
        pivot_df = counts.pivot(
            index='collection_year',
            columns=haplotype_col,
            values='pct'
        ).fillna(0)
        pivot_df = pivot_df.sort_index()

        # Get colors for this species' haplotypes
        colors = [color_map.get(h, '#999999') for h in pivot_df.columns]

        # Plot stacked area
        ax.stackplot(pivot_df.index, pivot_df.values.T, colors=colors, alpha=0.8)

        ax.set_title(species, fontsize=10, fontweight='bold')
        ax.set_ylim(0, 100)
        ax.set_ylabel('Abundance (%)')
        ax.set_xlabel('Year')

    # Hide unused subplots
    for idx in range(n_species, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    # Add shared legend
    all_haplotypes = valid_df[haplotype_col].unique()
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color_map.get(h, '#999999'), label=h, alpha=0.8)
        for h in sorted(all_haplotypes)
    ]
    fig.legend(
        handles=legend_elements,
        title='Haplotype',
        loc='center right',
        bbox_to_anchor=(1.12, 0.5),
        frameon=False
    )

    fig.suptitle('Haplotype Composition Over Time by Species', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved species-faceted temporal composition: {output_path}")


def plot_emergence_timeline_species_faceted(
    df: pd.DataFrame,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp',
    species_col: str = 'assigned_sp',
    figsize_per_facet: Tuple[int, int] = (10, 5),
    dpi: int = 300,
    max_species: int = 12
) -> None:
    """
    Plot species-faceted haplotype emergence timelines.

    Creates a separate subplot for each species showing when its haplotypes
    were first and last detected.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed dates, haplotype assignments, and species assignments
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    species_col : str, optional
        Column containing species assignments (default: 'assigned_sp')
    figsize_per_facet : Tuple[int, int], optional
        Figure size per facet in inches (default: 10x5)
    dpi : int, optional
        Resolution for output (default: 300)
    max_species : int, optional
        Maximum number of species to plot (default: 12)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if 'parsed_date' not in df.columns:
        logger.warning("No parsed_date column. Run parse_collection_dates first.")
        return

    if haplotype_col not in df.columns or species_col not in df.columns:
        logger.warning(f"Required columns not found for species-faceted emergence plot")
        return

    # Filter to valid data
    valid_df = df[
        df['parsed_date'].notna() &
        df[haplotype_col].notna() &
        df[species_col].notna()
    ].copy()

    if valid_df.empty:
        logger.warning("No valid data for species-faceted emergence plot")
        return

    # Get species ordered by sample count
    species_counts = valid_df[species_col].value_counts()
    species_list = species_counts.head(max_species).index.tolist()

    if len(species_list) == 0:
        return

    # Calculate grid dimensions
    n_species = len(species_list)
    n_cols = min(2, n_species)
    n_rows = (n_species + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_facet[0] * n_cols, figsize_per_facet[1] * n_rows),
        squeeze=False
    )

    # Get colors
    if color_map is None:
        color_map = get_haplotype_colors(valid_df, haplotype_col)

    for idx, species in enumerate(species_list):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        # Filter to this species
        species_df = valid_df[valid_df[species_col] == species]

        # Calculate emergence for this species
        emergence_data = []
        for haplotype, group in species_df.groupby(haplotype_col):
            first_date = group['parsed_date'].min()
            last_date = group['parsed_date'].max()
            emergence_data.append({
                'haplotype': haplotype,
                'first_date': first_date,
                'last_date': last_date,
                'n_samples': len(group)
            })

        if not emergence_data:
            ax.set_title(species, fontsize=10, fontweight='bold')
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            continue

        emergence_df = pd.DataFrame(emergence_data)
        emergence_df = emergence_df.sort_values('first_date')

        # Plot
        y_positions = range(len(emergence_df))

        for i, (_, hap_row) in enumerate(emergence_df.iterrows()):
            color = color_map.get(hap_row['haplotype'], '#3874F4')

            # Horizontal line from first to last
            ax.hlines(
                y=i,
                xmin=hap_row['first_date'],
                xmax=hap_row['last_date'],
                colors=color,
                linewidths=3,
                alpha=0.6
            )

            # First detection marker
            ax.scatter(
                hap_row['first_date'], i,
                c=[color],
                s=hap_row['n_samples'] * 2 + 30,
                marker='o',
                alpha=0.9,
                edgecolors='black',
                linewidths=0.5,
                zorder=10
            )

            # Last detection marker (if different)
            if hap_row['last_date'] > hap_row['first_date']:
                ax.scatter(
                    hap_row['last_date'], i,
                    c=[color],
                    s=hap_row['n_samples'] * 2 + 30,
                    marker='D',
                    alpha=0.9,
                    edgecolors='black',
                    linewidths=0.5,
                    zorder=10
                )

        ax.set_yticks(y_positions)
        ax.set_yticklabels(emergence_df['haplotype'], fontsize=8)
        ax.set_title(species, fontsize=10, fontweight='bold')
        ax.set_xlabel('Date')

    # Hide unused subplots
    for idx in range(n_species, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    # Add shared legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=8, label='First Detection'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='gray',
               markersize=8, label='Last Detection'),
        Line2D([0], [0], color='gray', linewidth=3, alpha=0.6, label='Observation Span'),
    ]
    fig.legend(
        handles=legend_elements,
        loc='center right',
        bbox_to_anchor=(1.08, 0.5),
        frameon=True,
        fancybox=True,
        title='Legend'
    )

    fig.suptitle('Haplotype Detection Timeline by Species\n(First ● and Last ◆ Detection)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved species-faceted emergence timeline: {output_path}")


def plot_temporal_composition_country_faceted(
    df: pd.DataFrame,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp',
    country_col: str = 'country/ocean',
    figsize_per_facet: Tuple[int, int] = (10, 4),
    dpi: int = 300,
    max_countries: int = 16
) -> None:
    """
    Plot country/ocean-faceted stacked area charts of haplotype composition over time.

    Creates a separate subplot for each country/ocean showing how haplotype
    composition changes over time within that geographic region. Useful for
    identifying whether certain haplotypes dominate over time in specific locations.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed dates, haplotype assignments, and country/ocean data
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    country_col : str, optional
        Column containing country/ocean names (default: 'country/ocean')
    figsize_per_facet : Tuple[int, int], optional
        Figure size per facet in inches (default: 10x4)
    dpi : int, optional
        Resolution for output (default: 300)
    max_countries : int, optional
        Maximum number of countries to plot (default: 16)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if 'collection_year' not in df.columns:
        logger.warning("No collection_year column. Run parse_collection_dates first.")
        return

    if haplotype_col not in df.columns or country_col not in df.columns:
        logger.warning(f"Required columns not found for country-faceted temporal plot")
        return

    # Filter to valid data
    valid_df = df[
        df['collection_year'].notna() &
        df[haplotype_col].notna() &
        df[country_col].notna()
    ].copy()

    if valid_df.empty:
        logger.warning("No valid data for country-faceted temporal plot")
        return

    # Get countries ordered by sample count
    country_counts = valid_df[country_col].value_counts()
    country_list = country_counts.head(max_countries).index.tolist()

    if len(country_list) == 0:
        return

    # Calculate grid dimensions
    n_countries = len(country_list)
    n_cols = min(2, n_countries)
    n_rows = (n_countries + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_facet[0] * n_cols, figsize_per_facet[1] * n_rows),
        squeeze=False
    )

    # Get colors
    if color_map is None:
        color_map = get_haplotype_colors(valid_df, haplotype_col)

    for idx, country in enumerate(country_list):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        # Filter to this country
        country_df = valid_df[valid_df[country_col] == country]

        # Calculate temporal distribution for this country
        counts = country_df.groupby(['collection_year', haplotype_col]).size().reset_index(name='n_samples')

        if counts.empty:
            ax.set_title(country, fontsize=10, fontweight='bold')
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            continue

        # Calculate percentages
        year_totals = counts.groupby('collection_year')['n_samples'].transform('sum')
        counts['pct'] = (counts['n_samples'] / year_totals * 100)

        # Pivot for stacked area
        pivot_df = counts.pivot(
            index='collection_year',
            columns=haplotype_col,
            values='pct'
        ).fillna(0)
        pivot_df = pivot_df.sort_index()

        # Get colors for this country's haplotypes
        colors = [color_map.get(h, '#999999') for h in pivot_df.columns]

        # Plot stacked area
        ax.stackplot(pivot_df.index, pivot_df.values.T, colors=colors, alpha=0.8)

        # Add sample count annotation
        total_samples = len(country_df)
        ax.set_title(f'{country}\n(n={total_samples:,})', fontsize=10, fontweight='bold')
        ax.set_ylim(0, 100)
        ax.set_ylabel('Abundance (%)')
        ax.set_xlabel('Year')

    # Hide unused subplots
    for idx in range(n_countries, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    # Add shared legend
    all_haplotypes = valid_df[haplotype_col].unique()
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color_map.get(h, '#999999'), label=h, alpha=0.8)
        for h in sorted(all_haplotypes)
    ]
    fig.legend(
        handles=legend_elements,
        title='Haplotype',
        loc='center right',
        bbox_to_anchor=(1.12, 0.5),
        frameon=False
    )

    fig.suptitle('Haplotype Composition Over Time by Country/Ocean', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved country-faceted temporal composition: {output_path}")


def plot_emergence_timeline_country_faceted(
    df: pd.DataFrame,
    output_path: Path,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp',
    country_col: str = 'country/ocean',
    figsize_per_facet: Tuple[int, int] = (10, 5),
    dpi: int = 300,
    max_countries: int = 16
) -> None:
    """
    Plot country/ocean-faceted haplotype emergence timelines.

    Creates a separate subplot for each country/ocean showing when haplotypes
    were first and last detected in that region.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with parsed dates, haplotype assignments, and country/ocean data
    output_path : Path
        Path to save the plot
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    country_col : str, optional
        Column containing country/ocean names (default: 'country/ocean')
    figsize_per_facet : Tuple[int, int], optional
        Figure size per facet in inches (default: 10x5)
    dpi : int, optional
        Resolution for output (default: 300)
    max_countries : int, optional
        Maximum number of countries to plot (default: 16)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if 'parsed_date' not in df.columns:
        logger.warning("No parsed_date column. Run parse_collection_dates first.")
        return

    if haplotype_col not in df.columns or country_col not in df.columns:
        logger.warning(f"Required columns not found for country-faceted emergence plot")
        return

    # Filter to valid data
    valid_df = df[
        df['parsed_date'].notna() &
        df[haplotype_col].notna() &
        df[country_col].notna()
    ].copy()

    if valid_df.empty:
        logger.warning("No valid data for country-faceted emergence plot")
        return

    # Get countries ordered by sample count
    country_counts = valid_df[country_col].value_counts()
    country_list = country_counts.head(max_countries).index.tolist()

    if len(country_list) == 0:
        return

    # Calculate grid dimensions
    n_countries = len(country_list)
    n_cols = min(2, n_countries)
    n_rows = (n_countries + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_facet[0] * n_cols, figsize_per_facet[1] * n_rows),
        squeeze=False
    )

    # Get colors
    if color_map is None:
        color_map = get_haplotype_colors(valid_df, haplotype_col)

    for idx, country in enumerate(country_list):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        # Filter to this country
        country_df = valid_df[valid_df[country_col] == country]

        # Calculate emergence for this country
        emergence_data = []
        for haplotype, group in country_df.groupby(haplotype_col):
            first_date = group['parsed_date'].min()
            last_date = group['parsed_date'].max()
            emergence_data.append({
                'haplotype': haplotype,
                'first_date': first_date,
                'last_date': last_date,
                'n_samples': len(group)
            })

        if not emergence_data:
            ax.set_title(country, fontsize=10, fontweight='bold')
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            continue

        emergence_df = pd.DataFrame(emergence_data)
        emergence_df = emergence_df.sort_values('first_date')

        # Plot
        y_positions = range(len(emergence_df))

        for i, (_, hap_row) in enumerate(emergence_df.iterrows()):
            color = color_map.get(hap_row['haplotype'], '#3874F4')

            # Horizontal line from first to last
            ax.hlines(
                y=i,
                xmin=hap_row['first_date'],
                xmax=hap_row['last_date'],
                colors=color,
                linewidths=3,
                alpha=0.6
            )

            # First detection marker
            ax.scatter(
                hap_row['first_date'], i,
                c=[color],
                s=hap_row['n_samples'] * 2 + 30,
                marker='o',
                alpha=0.9,
                edgecolors='black',
                linewidths=0.5,
                zorder=10
            )

            # Last detection marker (if different)
            if hap_row['last_date'] > hap_row['first_date']:
                ax.scatter(
                    hap_row['last_date'], i,
                    c=[color],
                    s=hap_row['n_samples'] * 2 + 30,
                    marker='D',
                    alpha=0.9,
                    edgecolors='black',
                    linewidths=0.5,
                    zorder=10
                )

        ax.set_yticks(y_positions)
        ax.set_yticklabels(emergence_df['haplotype'], fontsize=8)
        total_samples = len(country_df)
        ax.set_title(f'{country}\n(n={total_samples:,})', fontsize=10, fontweight='bold')
        ax.set_xlabel('Date')

    # Hide unused subplots
    for idx in range(n_countries, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    # Add shared legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=8, label='First Detection'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='gray',
               markersize=8, label='Last Detection'),
        Line2D([0], [0], color='gray', linewidth=3, alpha=0.6, label='Observation Span'),
    ]
    fig.legend(
        handles=legend_elements,
        loc='center right',
        bbox_to_anchor=(1.08, 0.5),
        frameon=True,
        fancybox=True,
        title='Legend'
    )

    fig.suptitle('Haplotype Detection Timeline by Country/Ocean\n(First ● and Last ◆ Detection)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    # Save in multiple formats, routing each to its correct directory
    _save_multi_format(fig, output_path, dpi=dpi)

    plt.close()
    logger.info(f"Saved country-faceted emergence timeline: {output_path}")


# =============================================================================
# Export Functions
# =============================================================================

def export_metadata_analysis(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    color_map: Optional[Dict[str, str]] = None,
    fields: Optional[List[str]] = None,
    haplotype_col: str = 'haplotype_sp',
    normalize_sex: bool = False,
    temporal_analysis: bool = True
) -> Dict[str, Path]:
    """
    Export all metadata analysis results.

    Runs all analyses and exports results to the specified directory.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments and metadata
    output_dir : Path
        Output directory
    organism : str
        Organism name for file naming
    color_map : Dict[str, str], optional
        Haplotype color mapping
    fields : List[str], optional
        Metadata fields to analyze (default: DEFAULT_CATEGORICAL_FIELDS)
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    normalize_sex : bool, optional
        Whether to normalize sex values (default: False)
    temporal_analysis : bool, optional
        Whether to run temporal analysis (default: True)

    Returns
    -------
    Dict[str, Path]
        Mapping of output names to file paths
    """
    output_dir = Path(output_dir)
    metadata_dir = output_dir / 'metadata_analysis'
    metadata_dir.mkdir(parents=True, exist_ok=True)

    # Visualization directories organized by format
    viz_pdf_dir = output_dir / 'visualization' / 'metadata' / 'pdf'
    viz_svg_dir = output_dir / 'visualization' / 'metadata' / 'svg'
    viz_json_dir = output_dir / 'visualization' / 'metadata' / 'json'
    viz_pdf_dir.mkdir(parents=True, exist_ok=True)
    viz_svg_dir.mkdir(parents=True, exist_ok=True)
    viz_json_dir.mkdir(parents=True, exist_ok=True)
    # Keep backward compatible reference
    viz_dir = viz_pdf_dir

    stats_dir = metadata_dir / 'statistical_tests'
    stats_dir.mkdir(parents=True, exist_ok=True)

    plots_data_dir = output_dir / 'plots' / 'data'
    plots_data_dir.mkdir(parents=True, exist_ok=True)

    if fields is None:
        fields = DEFAULT_CATEGORICAL_FIELDS

    # Get color map if not provided
    if color_map is None:
        color_map = get_haplotype_colors(df, haplotype_col)

    outputs = {}

    # Optional sex normalization
    if normalize_sex and 'sex' in df.columns:
        logger.info("Normalizing sex values...")
        df = normalize_sex_values(df, 'sex')

    # 1. Coverage Analysis
    logger.info("Analyzing metadata coverage...")
    coverage_stats = analyze_metadata_coverage(df, fields)

    # Save coverage JSON
    coverage_json_path = metadata_dir / f'{organism}_metadata_coverage.json'
    with open(coverage_json_path, 'w') as f:
        json.dump(coverage_stats, f, indent=2, default=str)
    outputs['coverage_json'] = coverage_json_path

    # Save coverage summary CSV
    coverage_rows = []
    for field, stats in coverage_stats['fields'].items():
        coverage_rows.append({
            'field': field,
            'n_total': coverage_stats['total_samples'],
            'n_with_value': stats['n_with_value'],
            'pct_coverage': f"{stats['pct_coverage']:.1f}%",
            'n_unique': stats['n_unique'],
            'top_values': '; '.join(stats['top_values'][:5]),
            'warning': stats.get('warning', '')
        })
    coverage_df = pd.DataFrame(coverage_rows)
    coverage_csv_path = metadata_dir / f'{organism}_metadata_summary.csv'
    coverage_df.to_csv(coverage_csv_path, index=False)
    outputs['coverage_csv'] = coverage_csv_path

    # Plot coverage
    coverage_plot_path = viz_dir / f'{organism}_metadata_coverage'
    plot_metadata_coverage(coverage_stats, coverage_plot_path)
    outputs['coverage_plot'] = coverage_plot_path.with_suffix('.png')

    # 2. Categorical Field Analysis
    logger.info("Analyzing categorical fields...")
    association_results = []

    for field in fields:
        if field not in df.columns:
            continue

        # Check if field has sufficient data
        field_stats = coverage_stats['fields'].get(field, {})
        if field_stats.get('n_with_value', 0) < 10:
            logger.debug(f"Skipping {field}: insufficient data")
            continue

        # Analyze field
        analysis_df = analyze_categorical_field(df, field, haplotype_col)
        if not analysis_df.empty:
            analysis_path = metadata_dir / f'{organism}_{field.replace("/", "_")}_analysis.csv'
            analysis_df.to_csv(analysis_path, index=False)
            outputs[f'{field}_analysis'] = analysis_path

        # Statistical test
        test_result = test_haplotype_association(df, field, haplotype_col)
        association_results.append(test_result)

        # Plot if enough categories
        if field_stats.get('n_unique', 0) >= 2:
            plot_path = viz_dir / f'{organism}_{field.replace("/", "_")}_by_haplotype'
            plot_categorical_by_haplotype(df, field, plot_path, color_map, haplotype_col)
            outputs[f'{field}_plot'] = plot_path.with_suffix('.png')

    # Save association test results
    if association_results:
        assoc_df = pd.DataFrame(association_results)
        assoc_path = stats_dir / f'{organism}_association_tests.csv'
        assoc_df.to_csv(assoc_path, index=False)
        outputs['association_tests'] = assoc_path

    # 3. Metadata Heatmap
    available_fields = [f for f in fields if f in df.columns and df[f].notna().any()]
    if available_fields and len(df[haplotype_col].dropna().unique()) > 1:
        heatmap_path = viz_dir / f'{organism}_metadata_heatmap'
        plot_metadata_heatmap(df, available_fields, heatmap_path, haplotype_col)
        outputs['heatmap'] = heatmap_path.with_suffix('.png')

    # 4. Temporal Analysis
    temporal_analysis_performed = False
    if temporal_analysis and DATE_FIELD in df.columns:
        logger.info("Analyzing temporal patterns...")

        # Parse dates
        df, parse_stats = parse_collection_dates(df, DATE_FIELD)

        if parse_stats['n_parsed'] > 0:
            temporal_analysis_performed = True

            # Haplotype emergence (now includes first and last detection)
            emergence_df = calculate_haplotype_emergence(df, haplotype_col)
            if not emergence_df.empty:
                emergence_path = metadata_dir / f'{organism}_haplotype_emergence.csv'
                emergence_df.to_csv(emergence_path, index=False)
                outputs['emergence_csv'] = emergence_path

                # Plot timeline (shows both first and last detection)
                timeline_path = viz_dir / f'{organism}_haplotype_emergence_timeline'
                plot_haplotype_emergence_timeline(emergence_df, timeline_path, color_map)
                outputs['emergence_plot'] = timeline_path.with_suffix('.png')

            # Temporal distribution
            temporal_df = analyze_temporal_distribution(df, haplotype_col, 'year')
            if not temporal_df.empty:
                temporal_path = metadata_dir / f'{organism}_temporal_trends.csv'
                temporal_df.to_csv(temporal_path, index=False)
                outputs['temporal_csv'] = temporal_path

                # Save temporal summary CSV
                temporal_summary_path = metadata_dir / f'{organism}_temporal_summary.csv'
                # Aggregate summary
                summary_data = temporal_df.groupby('time_period').agg({
                    'n_samples': 'sum',
                    'haplotype_sp': 'nunique'
                }).reset_index()
                summary_data.columns = ['year', 'total_samples', 'n_haplotypes']
                summary_data.to_csv(temporal_summary_path, index=False)
                outputs['temporal_summary'] = temporal_summary_path

                # Plot composition over time (all haplotypes)
                composition_path = viz_dir / f'{organism}_temporal_composition'
                plot_temporal_composition(temporal_df, composition_path, color_map)
                outputs['temporal_composition_plot'] = composition_path.with_suffix('.png')

            # 5. Species-Faceted Temporal Plots (if species data available)
            species_col = 'assigned_sp'
            if species_col in df.columns and df[species_col].notna().any():
                n_species = df[species_col].nunique()
                if n_species > 1:
                    logger.info(f"Generating species-faceted temporal plots for {n_species} species...")

                    # Species-faceted temporal composition
                    species_comp_path = viz_dir / f'{organism}_temporal_composition_species_faceted'
                    try:
                        plot_temporal_composition_species_faceted(
                            df, species_comp_path, color_map, haplotype_col, species_col
                        )
                        outputs['temporal_composition_species_faceted'] = species_comp_path.with_suffix('.png')
                    except Exception as e:
                        logger.warning(f"Species-faceted temporal composition failed: {e}")

                    # Species-faceted emergence timeline
                    species_emergence_path = viz_dir / f'{organism}_emergence_timeline_species_faceted'
                    try:
                        plot_emergence_timeline_species_faceted(
                            df, species_emergence_path, color_map, haplotype_col, species_col
                        )
                        outputs['emergence_timeline_species_faceted'] = species_emergence_path.with_suffix('.png')
                    except Exception as e:
                        logger.warning(f"Species-faceted emergence timeline failed: {e}")

            # 6. Country/Ocean-Faceted Temporal Plots (if geographic data available)
            country_col = 'country/ocean'
            if country_col in df.columns and df[country_col].notna().any():
                n_countries = df[country_col].nunique()
                if n_countries > 1:
                    logger.info(f"Generating country/ocean-faceted temporal plots for {n_countries} regions...")

                    # Country-faceted temporal composition
                    country_comp_path = viz_dir / f'{organism}_temporal_composition_country_faceted'
                    try:
                        plot_temporal_composition_country_faceted(
                            df, country_comp_path, color_map, haplotype_col, country_col
                        )
                        outputs['temporal_composition_country_faceted'] = country_comp_path.with_suffix('.png')
                    except Exception as e:
                        logger.warning(f"Country-faceted temporal composition failed: {e}")

                    # Country-faceted emergence timeline
                    country_emergence_path = viz_dir / f'{organism}_emergence_timeline_country_faceted'
                    try:
                        plot_emergence_timeline_country_faceted(
                            df, country_emergence_path, color_map, haplotype_col, country_col
                        )
                        outputs['emergence_timeline_country_faceted'] = country_emergence_path.with_suffix('.png')
                    except Exception as e:
                        logger.warning(f"Country-faceted emergence timeline failed: {e}")

    elif temporal_analysis and DATE_FIELD not in df.columns:
        logger.info(f"Temporal analysis skipped: '{DATE_FIELD}' field not found in data")

    outputs['_temporal_analysis_performed'] = temporal_analysis_performed

    logger.info(f"Metadata analysis complete. {len(outputs)} files exported to {metadata_dir}")

    return outputs


def run_metadata_analysis(
    annotated_df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    fields: Optional[List[str]] = None,
    normalize_sex: bool = False,
    temporal_analysis: bool = True,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp'
) -> Dict[str, Any]:
    """
    Main entry point for running metadata analysis.

    This function is called from the CLI pipeline.

    Parameters
    ----------
    annotated_df : pd.DataFrame
        Annotated DataFrame with haplotype assignments
    output_dir : Path
        Output directory
    organism : str
        Organism name
    fields : List[str], optional
        Metadata fields to analyze
    normalize_sex : bool, optional
        Whether to normalize sex values
    temporal_analysis : bool, optional
        Whether to include temporal analysis
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels

    Returns
    -------
    Dict[str, Any]
        Results including file paths and summary statistics
    """
    logger.info("Starting metadata analysis...")

    # Use default fields if none specified
    if fields is None:
        fields = DEFAULT_CATEGORICAL_FIELDS

    # Track which fields are available and which are missing
    available_fields = [f for f in fields if f in annotated_df.columns]
    missing_fields = [f for f in fields if f not in annotated_df.columns]

    if missing_fields:
        logger.info(f"Metadata fields not found in data (skipping): {', '.join(missing_fields)}")

    if not available_fields:
        logger.warning("No metadata fields found in data. Metadata analysis will be limited.")

    # Check for temporal analysis field
    temporal_field_available = DATE_FIELD in annotated_df.columns
    if temporal_analysis and not temporal_field_available:
        logger.info(f"Temporal analysis field '{DATE_FIELD}' not found. Temporal analysis will be skipped.")

    # Export all analyses (handles missing fields gracefully)
    outputs = export_metadata_analysis(
        df=annotated_df,
        output_dir=output_dir,
        organism=organism,
        color_map=color_map,
        fields=fields,  # Pass all fields; function handles missing ones
        haplotype_col=haplotype_col,
        normalize_sex=normalize_sex,
        temporal_analysis=temporal_analysis
    )

    # Build results summary
    temporal_performed = outputs.pop('_temporal_analysis_performed', False)

    results = {
        'output_files': outputs,
        'n_fields_analyzed': len([k for k in outputs if '_analysis' in k and not k.startswith('_')]),
        'n_visualizations': len([k for k in outputs if '_plot' in k or 'faceted' in k]),
        'temporal_analysis': temporal_performed,
        'available_fields': available_fields,
        'missing_fields': missing_fields,
        'temporal_field_available': temporal_field_available,
    }

    logger.info(f"Metadata analysis complete: {results['n_fields_analyzed']} fields analyzed, "
                f"{results['n_visualizations']} visualizations created")

    if missing_fields:
        logger.info(f"Missing fields: {', '.join(missing_fields)}")

    return results
