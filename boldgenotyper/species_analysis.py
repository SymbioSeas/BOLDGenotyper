"""
Species-Level Analysis Module

Aggregates haplotype-level data to species level for taxonomically-informed
analyses and visualizations. Provides species-level distribution, diversity,
and divergence metrics.

This module addresses the need for biologically meaningful groupings when
multiple haplotypes belong to the same species, while maintaining haplotype-
level resolution for fine-scale analyses.

Key Functions:
- aggregate_samples_by_species(): Group samples by assigned species
- calculate_species_diversity(): Haplotype diversity metrics per species
- calculate_species_divergence(): Within- and between-species divergence
- generate_species_geographic_summary(): Geographic distribution per species
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
import numpy as np
from collections import defaultdict

logger = logging.getLogger(__name__)


def aggregate_samples_by_species(
    annotated_metadata: pd.DataFrame,
    species_composition: pd.DataFrame,
    min_confidence: float = 0.7,
    output_dir: Optional[Path] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate sample and haplotype data by species.

    Parameters
    ----------
    annotated_metadata : pd.DataFrame
        Metadata with haplotype assignments (processid, haplotype_id, etc.)
    species_composition : pd.DataFrame
        Haplotype species composition (haplotype_id, primary_species, etc.)
    min_confidence : float, optional
        Minimum primary_species_pct to include (default: 0.7 = 70%)
    output_dir : Path, optional
        Directory to save output files

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - species_assignments: Sample-level species assignments
        - species_summary: Species-level summary statistics

    Notes
    -----
    Haplotypes with ambiguous taxonomy (primary_species_pct < min_confidence)
    are excluded from species-level aggregation but retained in haplotype-level
    analyses.
    """
    logger.info(f"Aggregating samples by species (min_confidence={min_confidence:.0%})")

    # Merge metadata with species composition
    # First, ensure primary_species_pct is numeric
    species_composition = species_composition.copy()
    if species_composition['primary_species_pct'].dtype == 'object':
        species_composition['primary_species_pct'] = (
            species_composition['primary_species_pct']
            .str.rstrip('%')
            .astype(float) / 100
        )

    # Merge to get species for each sample
    species_assignments = annotated_metadata.merge(
        species_composition[['haplotype_id', 'primary_species', 'primary_species_pct',
                             'is_ambiguous', 'is_multi_species']],
        on='haplotype_id',
        how='left'
    )

    # Filter to confident species assignments
    # Ensure is_ambiguous is boolean (may contain NaN floats)
    is_not_ambiguous = species_assignments['is_ambiguous'].fillna(False) == False

    confident_assignments = species_assignments[
        (species_assignments['primary_species_pct'] >= min_confidence) &
        is_not_ambiguous
    ].copy()

    n_total = len(species_assignments)
    n_confident = len(confident_assignments)
    n_species = confident_assignments['primary_species'].nunique()

    logger.info(f"  Samples with confident species: {n_confident}/{n_total} ({n_confident/n_total*100:.1f}%)")
    logger.info(f"  Number of species: {n_species}")

    # Calculate species-level summary
    species_summary = confident_assignments.groupby('primary_species').agg({
        'processid': 'count',  # Number of samples
        'haplotype_id': 'nunique',  # Number of haplotypes
    }).reset_index()

    species_summary.columns = ['species', 'n_samples', 'n_haplotypes']
    species_summary = species_summary.sort_values('n_samples', ascending=False)

    # Add haplotype list for each species
    species_haplotypes = (
        confident_assignments.groupby('primary_species')['haplotype_id']
        .apply(lambda x: ', '.join(sorted(x.unique())))
        .reset_index()
    )
    species_haplotypes.columns = ['species', 'haplotype_list']
    species_summary = species_summary.merge(species_haplotypes, on='species')

    # Save outputs if directory provided
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        assignments_file = output_dir / 'species_assignments.csv'
        species_assignments.to_csv(assignments_file, index=False)
        logger.info(f"  ✓ Saved species assignments: {assignments_file}")

        summary_file = output_dir / 'species_summary.csv'
        species_summary.to_csv(summary_file, index=False)
        logger.info(f"  ✓ Saved species summary: {summary_file}")

    return species_assignments, species_summary


def calculate_species_diversity(
    species_assignments: pd.DataFrame,
    output_dir: Optional[Path] = None
) -> pd.DataFrame:
    """
    Calculate haplotype diversity metrics per species.

    Parameters
    ----------
    species_assignments : pd.DataFrame
        Sample-level species assignments
    output_dir : Path, optional
        Directory to save output files

    Returns
    -------
    pd.DataFrame
        Species diversity metrics (nucleotide diversity, haplotype diversity, etc.)

    Notes
    -----
    Metrics calculated:
    - n_haplotypes: Number of unique haplotypes per species
    - n_samples: Number of samples per species
    - haplotype_diversity (h): Probability that two randomly chosen samples
      have different haplotypes
    - simpson_index: 1 - sum(p_i^2) where p_i is frequency of haplotype i
    """
    logger.info("Calculating species-level diversity metrics...")

    # Filter to confident species assignments
    is_not_ambiguous = species_assignments['is_ambiguous'].fillna(False) == False if 'is_ambiguous' in species_assignments.columns else True
    confident = species_assignments[
        is_not_ambiguous &
        (species_assignments['primary_species'].notna())
    ].copy()

    diversity_metrics = []

    for species, group in confident.groupby('primary_species'):
        n_samples = len(group)
        n_haplotypes = group['haplotype_id'].nunique()

        # Calculate haplotype frequencies
        haplotype_counts = group['haplotype_id'].value_counts()
        haplotype_freqs = haplotype_counts / n_samples

        # Haplotype diversity (h) = probability that two random samples differ
        # h = (n/(n-1)) * (1 - sum(p_i^2))
        sum_p_squared = (haplotype_freqs ** 2).sum()
        if n_samples > 1:
            h = (n_samples / (n_samples - 1)) * (1 - sum_p_squared)
        else:
            h = 0.0

        # Simpson's diversity index
        simpson = 1 - sum_p_squared

        diversity_metrics.append({
            'species': species,
            'n_samples': n_samples,
            'n_haplotypes': n_haplotypes,
            'haplotype_diversity': h,
            'simpson_index': simpson,
            'dominant_haplotype': haplotype_counts.index[0],
            'dominant_haplotype_freq': haplotype_freqs.iloc[0]
        })

    diversity_df = pd.DataFrame(diversity_metrics)
    diversity_df = diversity_df.sort_values('n_samples', ascending=False)

    if output_dir:
        output_dir = Path(output_dir)
        diversity_file = output_dir / 'species_diversity.csv'
        diversity_df.to_csv(diversity_file, index=False)
        logger.info(f"  ✓ Saved species diversity metrics: {diversity_file}")

    return diversity_df


def generate_species_geographic_summary(
    species_assignments: pd.DataFrame,
    output_dir: Optional[Path] = None
) -> pd.DataFrame:
    """
    Summarize geographic distribution per species.

    Parameters
    ----------
    species_assignments : pd.DataFrame
        Sample-level species assignments with geographic data
    output_dir : Path, optional
        Directory to save output files

    Returns
    -------
    pd.DataFrame
        Species geographic summary

    Notes
    -----
    Summarizes:
    - Number of samples per species per ocean basin
    - Number of samples per species per country
    - Geographic range extent per species
    """
    logger.info("Generating species-level geographic summary...")

    # Filter to confident species with coordinates
    confident = species_assignments[
        (~species_assignments['is_ambiguous']) &
        (species_assignments['primary_species'].notna())
    ].copy()

    # Ocean basin distribution
    if 'ocean_basin' in confident.columns:
        basin_summary = confident.groupby(['primary_species', 'ocean_basin']).size().reset_index(name='n_samples')
        basin_summary = basin_summary[basin_summary['ocean_basin'] != 'Unknown']

        if output_dir and len(basin_summary) > 0:
            output_dir = Path(output_dir)
            basin_file = output_dir / 'species_basin_distribution.csv'
            basin_summary.to_csv(basin_file, index=False)
            logger.info(f"  ✓ Saved species basin distribution: {basin_file}")

    # Country distribution
    if 'country/ocean' in confident.columns:
        country_summary = confident.groupby(['primary_species', 'country/ocean']).size().reset_index(name='n_samples')
        country_summary = country_summary[country_summary['country/ocean'].notna()]

        if output_dir and len(country_summary) > 0:
            output_dir = Path(output_dir)
            country_file = output_dir / 'species_country_distribution.csv'
            country_summary.to_csv(country_file, index=False)
            logger.info(f"  ✓ Saved species country distribution: {country_file}")

    # Geographic range extent (for species with coordinates)
    if 'coord' in confident.columns:
        range_metrics = []
        for species, group in confident.groupby('primary_species'):
            # Parse coordinates
            coords_with_data = group[group['coord'].notna()].copy()
            if len(coords_with_data) == 0:
                continue

            try:
                # Extract lat/lon from coord string [lat, lon]
                coords = coords_with_data['coord'].str.strip('[]').str.split(', ', expand=True)
                if len(coords.columns) >= 2:
                    lats = pd.to_numeric(coords[0], errors='coerce')
                    lons = pd.to_numeric(coords[1], errors='coerce')

                    valid = ~(lats.isna() | lons.isna())
                    if valid.sum() > 0:
                        lat_range = lats[valid].max() - lats[valid].min()
                        lon_range = lons[valid].max() - lons[valid].min()

                        range_metrics.append({
                            'species': species,
                            'n_samples_with_coords': valid.sum(),
                            'lat_min': lats[valid].min(),
                            'lat_max': lats[valid].max(),
                            'lon_min': lons[valid].min(),
                            'lon_max': lons[valid].max(),
                            'lat_range': lat_range,
                            'lon_range': lon_range
                        })
            except Exception as e:
                logger.debug(f"Could not parse coordinates for {species}: {e}")
                continue

        if range_metrics and output_dir:
            range_df = pd.DataFrame(range_metrics)
            output_dir = Path(output_dir)
            range_file = output_dir / 'species_geographic_range.csv'
            range_df.to_csv(range_file, index=False)
            logger.info(f"  ✓ Saved species geographic range: {range_file}")

    # Return basin summary as primary output
    return basin_summary if 'basin_summary' in locals() else pd.DataFrame()


def get_species_sample_list(
    species_assignments: pd.DataFrame,
    species: str,
    min_confidence: float = 0.7
) -> List[str]:
    """
    Get list of sample processids for a given species.

    Parameters
    ----------
    species_assignments : pd.DataFrame
        Sample-level species assignments
    species : str
        Species name to filter for
    min_confidence : float, optional
        Minimum confidence threshold

    Returns
    -------
    List[str]
        List of processids assigned to this species
    """
    mask = (
        (species_assignments['primary_species'] == species) &
        (species_assignments['primary_species_pct'] >= min_confidence) &
        (~species_assignments['is_ambiguous'])
    )
    return species_assignments[mask]['processid'].tolist()


def get_species_haplotype_list(
    species_assignments: pd.DataFrame,
    species: str,
    min_confidence: float = 0.7
) -> List[str]:
    """
    Get list of haplotype IDs for a given species.

    Parameters
    ----------
    species_assignments : pd.DataFrame
        Sample-level species assignments
    species : str
        Species name to filter for
    min_confidence : float, optional
        Minimum confidence threshold

    Returns
    -------
    List[str]
        List of haplotype IDs assigned to this species
    """
    mask = (
        (species_assignments['primary_species'] == species) &
        (species_assignments['primary_species_pct'] >= min_confidence) &
        (~species_assignments['is_ambiguous'])
    )
    return species_assignments[mask]['haplotype_id'].unique().tolist()


# ============================================================================
# Species-Faceted Haplotype Analysis
# ============================================================================


def generate_species_faceted_haplotype_subsets(
    species_assignments: pd.DataFrame,
    output_dir: Optional[Path] = None,
    min_haplotypes: int = 2
) -> Dict[str, pd.DataFrame]:
    """
    Create separate dataframes for each species containing all haplotypes.

    This enables species-specific analyses of haplotype-level patterns
    (geographic distribution, divergence, etc.) while filtering out
    species with insufficient data.

    Parameters
    ----------
    species_assignments : pd.DataFrame
        Sample-level species assignments with haplotype data
    output_dir : Path, optional
        Directory to save species-specific subset files
    min_haplotypes : int, optional
        Minimum number of haplotypes required to include a species (default: 2)

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping species names to their haplotype-level dataframes

    Notes
    -----
    Filters to confident species assignments (is_ambiguous=False) and
    excludes species with fewer than min_haplotypes unique haplotypes.
    """
    logger.info(f"Generating species-faceted haplotype subsets (min_haplotypes={min_haplotypes})...")

    # Filter to confident species assignments
    confident = species_assignments[
        (~species_assignments['is_ambiguous']) &
        (species_assignments['primary_species'].notna())
    ].copy()

    species_subsets = {}

    for species, group in confident.groupby('primary_species'):
        n_haplotypes = group['haplotype_id'].nunique()

        if n_haplotypes < min_haplotypes:
            logger.debug(f"  Skipping {species}: only {n_haplotypes} haplotype(s)")
            continue

        species_subsets[species] = group.copy()
        logger.info(f"  {species}: {len(group)} samples, {n_haplotypes} haplotypes")

        # Save subset if output directory provided
        if output_dir:
            output_dir = Path(output_dir)
            species_dir = output_dir / 'species_subsets'
            species_dir.mkdir(parents=True, exist_ok=True)

            # Sanitize species name for filename
            safe_name = species.replace(' ', '_').replace('/', '_')
            subset_file = species_dir / f"{safe_name}_haplotypes.csv"
            group.to_csv(subset_file, index=False)
            logger.debug(f"    Saved: {subset_file}")

    logger.info(f"  ✓ Created subsets for {len(species_subsets)} species")

    return species_subsets


def generate_within_species_divergence_matrices(
    divergence_matrix: pd.DataFrame,
    species_summary_csv: Path,
    haplotype_taxonomy_csv: Path,
    output_dir: Path,
    min_haplotypes: int = 2
) -> Dict[str, Any]:
    """
    Create separate divergence matrices for haplotypes within each species.

    Enables visualization of intraspecific divergence patterns, geographic
    structure, and potential cryptic lineages.

    Parameters
    ----------
    divergence_matrix : pd.DataFrame
        Pairwise divergence matrix from haplotype-level analysis
    species_summary_csv : Path
        Species summary CSV from species_analysis module
    haplotype_taxonomy_csv : Path
        Haplotype taxonomy CSV with species assignments
    output_dir : Path
        Output directory for species-specific divergence matrices
    min_haplotypes : int, optional
        Minimum number of haplotypes to create a matrix (default: 2)

    Returns
    -------
    dict
        Results including file paths and summary statistics

    Notes
    -----
    Creates a separate CSV file for each species containing the pairwise
    divergence matrix for all haplotypes assigned to that species.
    """
    output_dir = Path(output_dir)
    species_matrices_dir = output_dir / 'species_divergence_matrices'
    species_matrices_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating within-species divergence matrices...")

    results = {'matrices': {}, 'stats': []}

    # Load data
    species_summary = pd.read_csv(species_summary_csv)
    haplotype_taxonomy = pd.read_csv(haplotype_taxonomy_csv)

    # Create haplotype-to-species mapping
    haplotype_to_species = dict(zip(
        haplotype_taxonomy['haplotype_id'],
        haplotype_taxonomy['assigned_sp']
    ))

    for _, row in species_summary.iterrows():
        species = row['species']
        n_haplotypes = row['n_haplotypes']

        if n_haplotypes < min_haplotypes:
            logger.debug(f"  Skipping {species}: only {n_haplotypes} haplotype(s)")
            continue

        # Get haplotypes for this species
        species_haplotypes = [
            h for h in divergence_matrix.index
            if haplotype_to_species.get(h) == species
        ]

        if len(species_haplotypes) < min_haplotypes:
            continue

        # Extract submatrix
        species_matrix = divergence_matrix.loc[species_haplotypes, species_haplotypes].copy()

        # Calculate summary statistics
        # Get upper triangle values (excluding diagonal)
        mask = np.triu(np.ones(species_matrix.shape), k=1).astype(bool)
        divergence_values = species_matrix.values[mask]
        divergence_values = divergence_values[~np.isnan(divergence_values)]

        if len(divergence_values) > 0:
            stats = {
                'species': species,
                'n_haplotypes': len(species_haplotypes),
                'n_comparisons': len(divergence_values),
                'mean_divergence': np.mean(divergence_values),
                'max_divergence': np.max(divergence_values),
                'min_divergence': np.min(divergence_values),
                'std_divergence': np.std(divergence_values)
            }
            results['stats'].append(stats)

        # Save matrix
        safe_name = species.replace(' ', '_').replace('/', '_')
        matrix_file = species_matrices_dir / f"{safe_name}_divergence_matrix.csv"
        species_matrix.to_csv(matrix_file)
        results['matrices'][species] = matrix_file

        logger.info(f"  ✓ {species}: {len(species_haplotypes)} haplotypes, "
                   f"mean divergence = {stats['mean_divergence']:.4f}")

    # Save summary statistics
    if results['stats']:
        stats_df = pd.DataFrame(results['stats'])
        stats_file = species_matrices_dir / 'within_species_divergence_summary.csv'
        stats_df.to_csv(stats_file, index=False)
        results['stats_file'] = stats_file
        logger.info(f"  ✓ Saved divergence summary: {stats_file}")

    logger.info(f"  ✓ Created matrices for {len(results['matrices'])} species")

    return results
