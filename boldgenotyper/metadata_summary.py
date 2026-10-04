"""
Metadata Summarization Module

Descriptive summaries of how BOLDGenotyper haplotype assignments distribute
across non-geographic BOLD metadata fields. This module is purely
diagnostic: it produces per-field tabulations (haplotype × value counts)
and temporal coverage (first and last sampling dates, collection-date
trends) to help users assess sampling bias and coverage in a BOLD pull.

It does NOT perform haplotype × metadata hypothesis tests. The historical
chi-square independence tests were removed because (a) the contingency
tables produced by BOLD-scale data are too sparse for the chi-square
approximation to be valid (>20% of expected counts < 5 in every case
study we tested), (b) there is no biologically grounded null hypothesis
being evaluated, and (c) running 10+ tests per dataset with no
multiple-testing correction is statistically indefensible. Users who
require formal tests of population structure should use the population
genetics exports (Arlequin .arp, DnaSP .fas) downstream.

Categorical Fields Tabulated:
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
Values are presented as-is without standardization. Users know their data
best and should decide how to standardize values. Optional sex
normalization is available via the --normalize-sex flag.

Author: Steph Smith (symbioseas@outlook.com)
"""

import logging
import re
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)

# Specimen and sampling fields summarized by default. BOLD's realm, biome and ecoregion
# are derived by BOLD from the coordinates using terrestrial ecoregions (misleading for
# marine and freshwater taxa), geoid is a numeric location code and country_iso repeats
# country/ocean; name them in --metadata-fields to summarize them.
DEFAULT_CATEGORICAL_FIELDS = [
    'sex',
    'life_stage',
    'reproduction',
    'habitat',
    'country/ocean',
    'province/state',
]

DATE_FIELD = 'collection_date_start'


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


# Field-specific codes in BOLD exports (lowercased key -> normalized value; None = missing)
VOCABULARY_CODES = {
    "sex": {"m": "Male", "male": "Male", "f": "Female", "female": "Female",
            "u": None, "unknown": None},
    "life_stage": {"a": "Adult", "i": "Immature", "m": "Mature", "j": "Juvenile"},
    "reproduction": {"s": "Sexual", "a": "Asexual", "u": None},
}
_GIS_SUFFIX = re.compile(r"_simplify-[\d.]+_buffered-[\d.]+")
_WHOLE_NUMBER = re.compile(r"^-?\d+\.0$")


def _normalize_term(field: str, value, map_codes: bool = True):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return np.nan
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    text = str(value).strip()
    if not text:
        return np.nan
    if _WHOLE_NUMBER.match(text):
        text = text[:-2]            # integer codes written as float (BOLD geoid '266.0')
    if not map_codes:
        return text
    codes = VOCABULARY_CODES.get(field, {})
    if text.lower() in codes:
        mapped = codes[text.lower()]
        return np.nan if mapped is None else mapped
    if _GIS_SUFFIX.search(text) or "_" in text:
        parts = [_GIS_SUFFIX.sub("", p.strip()).replace("_", " ").strip() for p in text.split(",")]
        text = ", ".join(p for p in parts if p)
    if text.islower():
        text = text[0].upper() + text[1:]
    return text


def normalize_vocabulary(df: pd.DataFrame, fields: List[str], map_sex_codes: bool = True) -> pd.DataFrame:
    """
    Normalize metadata vocabularies before tabulation (P1-5).

    - Empty and whitespace-only values are missing (not a category).
    - Field codes: sex M/F -> Male/Female and U -> missing; life_stage
      A/I/M/J -> Adult/Immature/Mature/Juvenile; reproduction S/A ->
      Sexual/Asexual and U -> missing.
    - All-lowercase values are capitalized, so 'adult' and 'Adult' merge.
    - BOLD GIS suffixes (e.g. '_simplify-0.001_buffered-0.018') are removed and
      underscores become spaces, per comma-separated item.
    - Integer codes written as float lose the decimal ('266.0' -> '266').

    map_sex_codes=False (--no-normalize-sex) leaves sex values as given, apart
    from blank values becoming missing.
    """
    df = df.copy()
    for field in fields:
        if field in df.columns:
            map_codes = map_sex_codes or field != 'sex'
            df[field] = df[field].map(lambda v, f=field, m=map_codes: _normalize_term(f, v, m))
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


# NOTE: The historical chi-square haplotype × metadata association test was
# removed in the manuscript-prep refactor (see module docstring for rationale).
# Users who require formal hypothesis testing of population structure should
# use the Arlequin / DnaSP exports from popgen_export.py downstream.


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
    First and last collection dates of each haplotype's records (--extended).

    These are sampling dates in BOLD: they reflect when and where studies
    collected specimens, not when a haplotype arose or was lost.

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
        - first_sampled, first_sampled_year
        - first_location (country/ocean if available; empty otherwise)
        - last_sampled, last_sampled_year
        - last_location (country/ocean if available; empty otherwise)
        - years_between_first_and_last_sample
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
            'first_sampled': first_row['parsed_date'].strftime('%Y-%m-%d'),
            'first_sampled_year': int(first_row['collection_year']),
            'last_sampled': last_row['parsed_date'].strftime('%Y-%m-%d'),
            'last_sampled_year': int(last_row['collection_year']),
            'years_between_first_and_last_sample': int(last_row['collection_year']) - int(first_row['collection_year']),
            'n_total_samples': len(group)
        }

        # Add first location if available
        if 'country/ocean' in df.columns and pd.notna(first_row.get('country/ocean')):
            result['first_location'] = first_row['country/ocean']
        else:
            result['first_location'] = np.nan

        # Add last location if available
        if 'country/ocean' in df.columns and pd.notna(last_row.get('country/ocean')):
            result['last_location'] = last_row['country/ocean']
        else:
            result['last_location'] = np.nan

        results.append(result)

    emergence_df = pd.DataFrame(results)

    # Sort by first detection date
    emergence_df = emergence_df.sort_values('first_sampled')

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

# =============================================================================
# Export Functions
# =============================================================================

def export_metadata_summary(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    color_map: Optional[Dict[str, str]] = None,
    fields: Optional[List[str]] = None,
    haplotype_col: str = 'haplotype_sp',
    normalize_sex: bool = True,
    temporal_analysis: bool = True,
    extended: bool = False,
    cross_tabs: bool = True
) -> Dict[str, Path]:
    """
    Export metadata summaries.

    Always writes field coverage (JSON, CSV) and metadata_summary/README.md.
    With extended=True it also writes cross-tabulations on normalized
    vocabularies and the heatmap (if cross_tabs is set) and the temporal and
    first-sampled outputs (if temporal_analysis is set) (P1-3, P1-4, P1-5).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with haplotype assignments and metadata
    output_dir : Path
        Output directory
    organism : str
        Organism name for file naming
    color_map : Dict[str, str], optional
        Unused; kept for callers. Figures are drawn by metadata_figures.
    fields : List[str], optional
        Metadata fields to analyze (default: DEFAULT_CATEGORICAL_FIELDS)
    haplotype_col : str, optional
        Column containing haplotype labels (default: 'haplotype_sp')
    normalize_sex : bool, optional
        Map sex codes (M/F -> Male/Female, U -> missing) (default: True)
    temporal_analysis : bool, optional
        Whether to run temporal analysis with extended=True (default: True)
    extended : bool, optional
        Write the cross-tabulations and temporal outputs (default: False)
    cross_tabs : bool, optional
        With extended=True, write the cross-tabulations and heatmap (default: True)

    Returns
    -------
    Dict[str, Path]
        Mapping of output names to file paths
    """
    output_dir = Path(output_dir)
    metadata_dir = output_dir / 'metadata_summary'
    metadata_dir.mkdir(parents=True, exist_ok=True)

    # Figures are drawn by metadata_figures from records.csv (run_figures)
    if fields is None:
        fields = DEFAULT_CATEGORICAL_FIELDS

    outputs = {}

    # Empty values are missing and codes are expanded before anything is counted
    df = normalize_vocabulary(df, [f for f in fields if f in df.columns], map_sex_codes=normalize_sex)
    if normalize_sex and 'sex' in df.columns:
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

    temporal_analysis_performed = False
    if not extended:
        outputs['readme'] = write_metadata_readme(metadata_dir, organism, outputs)
        outputs['_temporal_analysis_performed'] = temporal_analysis_performed
        logger.info("Metadata coverage written; cross-tabulations and temporal summaries need --extended")
        return outputs

    if cross_tabs:
        # 2. Categorical Field Tabulation (descriptive only — no hypothesis tests)
        logger.info("Tabulating haplotype × metadata distributions...")

        for field in fields:
            if field not in df.columns:
                continue

            # Check if field has sufficient data
            field_stats = coverage_stats['fields'].get(field, {})
            if field_stats.get('n_with_value', 0) < 10:
                logger.debug(f"Skipping {field}: insufficient data")
                continue

            # Per-haplotype × per-value count table
            analysis_df = analyze_categorical_field(df, field, haplotype_col)
            if not analysis_df.empty:
                analysis_path = metadata_dir / f'{organism}_{field.replace("/", "_")}_tabulation.csv'
                analysis_df.to_csv(analysis_path, index=False)
                outputs[f'{field}_tabulation'] = analysis_path


    # 4. Temporal Analysis
    if temporal_analysis and DATE_FIELD in df.columns:
        logger.info("Analyzing temporal patterns...")
        logger.info(f"  Note: {SAMPLING_CAVEAT}")

        # Parse dates
        df, parse_stats = parse_collection_dates(df, DATE_FIELD)

        if parse_stats['n_parsed'] > 0:
            temporal_analysis_performed = True

            # First and last sampling date per haplotype
            first_sampled_df = calculate_haplotype_emergence(df, haplotype_col)
            if not first_sampled_df.empty:
                first_sampled_path = metadata_dir / f'{organism}_haplotype_first_sampled.csv'
                first_sampled_df.to_csv(first_sampled_path, index=False)
                outputs['first_sampled_csv'] = first_sampled_path

            # Temporal distribution
            temporal_df = analyze_temporal_distribution(df, haplotype_col, 'year')
            if not temporal_df.empty:
                temporal_path = metadata_dir / f'{organism}_temporal_trends.csv'
                with open(temporal_path, 'w') as f:
                    f.write(f"# {SAMPLING_CAVEAT}\n")
                    temporal_df.to_csv(f, index=False)
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

    elif temporal_analysis and DATE_FIELD not in df.columns:
        logger.info(f"Temporal analysis skipped: '{DATE_FIELD}' field not found in data")

    outputs['readme'] = write_metadata_readme(metadata_dir, organism, outputs)
    outputs['_temporal_analysis_performed'] = temporal_analysis_performed

    logger.info(f"Metadata summarization complete. {len(outputs)} files exported to {metadata_dir}")

    return outputs


SAMPLING_CAVEAT = ("Counts by collection year reflect when and where specimens were sampled "
                   "and deposited in BOLD, not changes in haplotype frequency.")


def write_metadata_readme(metadata_dir: Path, organism: str, outputs: Dict[str, Path]) -> Path:
    """Describe the files in metadata_summary/ and how to read them."""
    lines = [
        "# Metadata summary",
        "",
        "Descriptive summaries of BOLD record metadata. No hypothesis tests are run.",
        "Values are normalized before counting: empty and whitespace-only values are",
        "missing, field codes are expanded (sex M/F -> Male/Female, U -> missing; life",
        "stage A/I/M/J -> Adult/Immature/Mature/Juvenile; reproduction S/A ->",
        "Sexual/Asexual), case variants are merged and BOLD GIS suffixes are removed.",
        "",
        "## Files",
        "",
        f"- `{organism}_metadata_coverage.json`, `{organism}_metadata_summary.csv`: for each",
        "  field, the number and percentage of records with a value, the number of",
        "  distinct values and the most common values.",
    ]
    if any(k.endswith('_tabulation') for k in outputs):
        lines += [
            f"- `{organism}_<field>_tabulation.csv`: records per haplotype and field value",
            "  (`n_samples`, `pct_of_haplotype`, `pct_of_value`). Records without a value",
            "  are left out.",
        ]
    if 'first_sampled_csv' in outputs:
        lines += [
            f"- `{organism}_haplotype_first_sampled.csv`: the first and last collection date",
            "  among each haplotype's records (`first_sampled`, `last_sampled`), with the",
            "  country/ocean of those records. These are sampling dates, not dates of origin",
            "  or loss.",
        ]
    if 'temporal_csv' in outputs:
        lines += [
            f"- `{organism}_temporal_trends.csv`, `{organism}_temporal_summary.csv`: records",
            "  per haplotype and collection year.",
        ]
    if 'first_sampled_csv' in outputs or 'temporal_csv' in outputs:
        lines += ["", f"Note: {SAMPLING_CAVEAT}"]
    if 'first_sampled_csv' not in outputs and not any(k.endswith('_tabulation') for k in outputs):
        lines += ["", "Cross-tabulations and collection-date summaries are written with `--extended`."]
    path = Path(metadata_dir) / "README.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def run_metadata_summary(
    annotated_df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    fields: Optional[List[str]] = None,
    normalize_sex: bool = True,
    temporal_analysis: bool = True,
    color_map: Optional[Dict[str, str]] = None,
    haplotype_col: str = 'haplotype_sp',
    extended: bool = False,
    cross_tabs: bool = True
) -> Dict[str, Any]:
    """
    Main entry point for running metadata summarization.

    Produces descriptive (non-statistical) summaries of how haplotype
    assignments distribute across BOLD metadata fields, plus temporal
    coverage if a collection-date column is present. Called from the
    CLI pipeline.

    Parameters
    ----------
    annotated_df : pd.DataFrame
        Annotated DataFrame with haplotype assignments
    output_dir : Path
        Output directory
    organism : str
        Organism name
    fields : List[str], optional
        Metadata fields to tabulate
    normalize_sex : bool, optional
        Map sex codes (M/F -> Male/Female, U -> missing) (default: True)
    temporal_analysis : bool, optional
        Whether to include temporal coverage analysis
    color_map : Dict[str, str], optional
        Haplotype color mapping
    haplotype_col : str, optional
        Column containing haplotype labels
    extended : bool, optional
        Write cross-tabulations and temporal outputs (default: False)
    cross_tabs : bool, optional
        With extended=True, write the cross-tabulations (default: True)

    Returns
    -------
    Dict[str, Any]
        Results including file paths and summary statistics
    """
    logger.info("Starting metadata summarization...")

    if fields is None:
        fields = DEFAULT_CATEGORICAL_FIELDS

    available_fields = [f for f in fields if f in annotated_df.columns]
    missing_fields = [f for f in fields if f not in annotated_df.columns]

    if missing_fields:
        logger.info(f"Metadata fields not found in data (skipping): {', '.join(missing_fields)}")

    if not available_fields:
        logger.warning("No metadata fields found in data. Metadata summary will be limited.")

    temporal_field_available = DATE_FIELD in annotated_df.columns
    if extended and temporal_analysis and not temporal_field_available:
        logger.info(f"Temporal field '{DATE_FIELD}' not found. Temporal analysis will be skipped.")

    outputs = export_metadata_summary(
        df=annotated_df,
        output_dir=output_dir,
        organism=organism,
        color_map=color_map,
        fields=fields,
        haplotype_col=haplotype_col,
        normalize_sex=normalize_sex,
        temporal_analysis=temporal_analysis,
        extended=extended,
        cross_tabs=cross_tabs
    )

    temporal_performed = outputs.pop('_temporal_analysis_performed', False)

    results = {
        'output_files': outputs,
        'n_fields_tabulated': len([k for k in outputs if k.endswith('_tabulation')]),
        'n_visualizations': len([k for k in outputs if '_plot' in k or 'faceted' in k]),
        'temporal_analysis': temporal_performed,
        'available_fields': available_fields,
        'missing_fields': missing_fields,
        'temporal_field_available': temporal_field_available,
    }

    logger.info(
        f"Metadata summarization complete: {results['n_fields_tabulated']} fields tabulated, "
        f"{results['n_visualizations']} visualizations created"
    )

    if missing_fields:
        logger.info(f"Missing fields: {', '.join(missing_fields)}")

    return results
