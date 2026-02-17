#!/usr/bin/env python3
"""
compare_analyses.py — Standalone comparative analysis utility for BOLDGenotyper

Compares species-level and family/genus-level BOLDGenotyper analyses to detect
database contamination and potential misidentifications.

Usage:
    python scripts/compare_analyses.py \\
        --species-level Sphyrna_lewini_output/ \\
        --family-level Sphyrnidae_output/ \\
        --output comparative_analysis/ \\
        --generate-reassignment-table

Requirements:
    pandas, numpy (standard BOLDGenotyper dependencies)

Author: Steph Smith (symbioseas@outlook.com)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from datetime import datetime

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class ComparativeAnalysisError(Exception):
    """Base exception for comparative analysis errors."""
    pass


def _find_diagnostics_path(csv_path: Path, organism: str) -> Optional[Path]:
    """
    Locate the assignment diagnostics file for an analysis.

    Supports legacy and current filenames and minor organism name variants.
    """
    if csv_path is None:
        return None

    diag_dir = csv_path.parent / "genotype_assignments"
    name_variants = {
        organism,
        organism.replace("_", " "),
        organism.replace(" ", "_")
    }

    candidates = []
    for name in name_variants:
        candidates.extend([
            diag_dir / f"{name}_assignment_diagnostics.csv",
            diag_dir / f"{name}_diagnostics.csv"
        ])

    for path in candidates:
        if path.exists():
            return path

    return None


def _load_identity_data(
    df: pd.DataFrame,
    meta: Dict[str, Any]
) -> Tuple[Dict[str, float], pd.DataFrame]:
    """
    Extract identity statistics and per-sample identities.

    Returns a stats dict plus a dataframe with `processid` and `identity`
    (scaled to percentages when source data are proportions).
    """
    csv_path = Path(meta.get('path')) if meta.get('path') else None
    organism = meta.get('organism', '')

    diag_path = _find_diagnostics_path(csv_path, organism) if csv_path else None
    identities_df = pd.DataFrame(columns=['processid', 'identity'])

    if diag_path is not None:
        try:
            diag_df = pd.read_csv(diag_path)
            if 'identity' in diag_df.columns and 'processid' in diag_df.columns:
                identities_df = diag_df[['processid', 'identity']].copy()
        except Exception as e:
            logger.debug(f"Could not read diagnostics at {diag_path}: {e}")

    if identities_df.empty and 'identity' in df.columns:
        identities_df = df[['processid', 'identity']].copy()

    if identities_df.empty:
        return {'mean_identity': np.nan, 'median_identity': np.nan}, identities_df

    # Coerce identity to numeric
    identities_df['identity'] = pd.to_numeric(identities_df['identity'], errors='coerce')
    identities_df = identities_df.dropna(subset=['identity'])

    if identities_df.empty:
        return {'mean_identity': np.nan, 'median_identity': np.nan}, identities_df

    # Normalize to percent if stored as proportion (0-1)
    if identities_df['identity'].max() <= 1.0:
        identities_df = identities_df.copy()
        identities_df['identity'] = identities_df['identity'] * 100

    stats = {
        'mean_identity': identities_df['identity'].mean(),
        'median_identity': identities_df['identity'].median()
    }

    return stats, identities_df


def _load_taxonomy_table(meta: Dict[str, Any]) -> Optional[pd.DataFrame]:
    """
    Load consensus taxonomy table to recover majority species/pct per group.
    """
    csv_path = Path(meta.get('path')) if meta.get('path') else None
    organism = meta.get('organism', '')
    if csv_path is None:
        return None

    tax_dir = csv_path.parent / "taxonomy"
    candidates = [
        tax_dir / f"{organism}_consensus_taxonomy.csv",
        tax_dir / f"{organism.replace('_', ' ')}_consensus_taxonomy.csv",
        tax_dir / f"{organism.replace(' ', '_')}_consensus_taxonomy.csv"
    ]

    tax_path = next((p for p in candidates if p.exists()), None)
    if tax_path is None:
        return None

    try:
        tax_df = pd.read_csv(tax_path)
    except Exception as e:
        logger.debug(f"Could not read taxonomy table {tax_path}: {e}")
        return None

    # Standardize columns
    if 'assigned_sp' in tax_df.columns:
        tax_df = tax_df.rename(columns={'assigned_sp': 'group_majority_species'})

    if 'majority_fraction' in tax_df.columns and 'group_majority_pct' not in tax_df.columns:
        pct = tax_df['majority_fraction']
        tax_df['group_majority_pct'] = pct * 100 if pct.max() <= 1 else pct

    keep_cols = [c for c in ['consensus_group', 'group_majority_species', 'group_majority_pct'] if c in tax_df.columns]
    if not keep_cols:
        return None

    return tax_df[keep_cols]


def load_analysis_results(
    analysis_path: Path
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Load results from a previous BOLDGenotyper analysis.

    Parameters
    ----------
    analysis_path : Path
        Path to analysis output directory or annotated CSV file

    Returns
    -------
    tuple
        (annotated_df, metadata_dict)
    """
    analysis_path = Path(analysis_path)

    # Check if it's a directory or file
    if analysis_path.is_dir():
        # Look for annotated CSV
        csv_files = list(analysis_path.glob("*_annotated.csv"))
        if not csv_files:
            raise ComparativeAnalysisError(f"No annotated CSV found in {analysis_path}")
        csv_path = csv_files[0]
        organism = csv_path.stem.replace("_annotated", "")
    elif analysis_path.is_file() and analysis_path.suffix == '.csv':
        csv_path = analysis_path
        organism = csv_path.stem.replace("_annotated", "")
    else:
        raise ComparativeAnalysisError(f"Invalid analysis path: {analysis_path}")

    logger.info(f"  Loading analysis data from {csv_path}")
    df = pd.read_csv(csv_path)

    # Load parameters if available
    params_path = csv_path.parent / f"{organism}_pipeline_parameters.json"
    if params_path.exists():
        import json
        with open(params_path) as f:
            params = json.load(f)
    else:
        params = {}

    metadata = {
        'organism': organism,
        'path': csv_path,
        'params': params,
        'n_samples': len(df)
    }

    return df, metadata


def generate_comparison_summary(
    species_df: pd.DataFrame,
    family_df: pd.DataFrame,
    species_meta: Dict[str, Any],
    family_meta: Dict[str, Any],
    species_identity_stats: Optional[Dict[str, float]] = None,
    family_identity_stats: Optional[Dict[str, float]] = None
) -> pd.DataFrame:
    """
    Generate high-level comparison summary metrics.

    Parameters
    ----------
    species_df : pd.DataFrame
        Species-level analysis dataframe
    family_df : pd.DataFrame
        Family-level analysis dataframe
    species_meta : dict
        Species-level metadata
    family_meta : dict
        Family-level metadata

    Returns
    -------
    pd.DataFrame
        Summary comparison table
    """
    # Calculate metrics for both analyses
    metrics = []

    # Total samples
    metrics.append({
        'metric': 'total_samples',
        'species_level': len(species_df),
        'family_level': len(family_df)
    })

    # Consensus groups
    species_groups = species_df['consensus_group'].nunique() if 'consensus_group' in species_df.columns else 0
    family_groups = family_df['consensus_group'].nunique() if 'consensus_group' in family_df.columns else 0
    metrics.append({
        'metric': 'consensus_groups',
        'species_level': species_groups,
        'family_level': family_groups
    })

    # Species detected
    species_detected_sp = species_df['species'].nunique() if 'species' in species_df.columns else 0
    family_detected_sp = family_df['species'].nunique() if 'species' in family_df.columns else 0
    metrics.append({
        'metric': 'species_detected',
        'species_level': species_detected_sp,
        'family_level': family_detected_sp
    })

    # Mean identity
    if species_identity_stats is None:
        species_identity_stats, _ = _load_identity_data(species_df, species_meta)
    if family_identity_stats is None:
        family_identity_stats, _ = _load_identity_data(family_df, family_meta)

    species_mean_id = species_identity_stats.get('mean_identity', np.nan)
    family_mean_id = family_identity_stats.get('mean_identity', np.nan)
    metrics.append({
        'metric': 'mean_identity',
        'species_level': species_mean_id,
        'family_level': family_mean_id
    })

    # Groups with 100% majority
    if 'group_majority_pct' in species_df.columns:
        species_100pct = species_df.groupby('consensus_group')['group_majority_pct'].first()
        species_100pct_count = (species_100pct >= 99.9).sum()
    else:
        species_100pct_count = 0

    if 'group_majority_pct' in family_df.columns:
        family_100pct = family_df.groupby('consensus_group')['group_majority_pct'].first()
        family_100pct_count = (family_100pct >= 99.9).sum()
    else:
        family_100pct_count = 0

    metrics.append({
        'metric': 'groups_with_100pct_majority',
        'species_level': species_100pct_count,
        'family_level': family_100pct_count
    })

    # Mixed species groups
    if 'group_species_count' in species_df.columns:
        species_mixed = species_df.groupby('consensus_group')['group_species_count'].first()
        species_mixed_count = (species_mixed > 1).sum()
    else:
        species_mixed_count = 0

    if 'group_species_count' in family_df.columns:
        family_mixed = family_df.groupby('consensus_group')['group_species_count'].first()
        family_mixed_count = (family_mixed > 1).sum()
    else:
        family_mixed_count = 0

    metrics.append({
        'metric': 'mixed_species_groups',
        'species_level': species_mixed_count,
        'family_level': family_mixed_count
    })

    # Potential misidentifications
    if 'potential_misidentification' in species_df.columns:
        species_misid = species_df['potential_misidentification'].sum()
    else:
        species_misid = 0

    if 'potential_misidentification' in family_df.columns:
        family_misid = family_df['potential_misidentification'].sum()
    else:
        family_misid = 0

    metrics.append({
        'metric': 'potential_misidentifications',
        'species_level': species_misid,
        'family_level': family_misid
    })

    # Create dataframe and calculate differences
    summary_df = pd.DataFrame(metrics)

    # Calculate difference and percent change
    summary_df['difference'] = summary_df['family_level'] - summary_df['species_level']

    def calc_pct_change(row):
        if pd.isna(row['species_level']) or row['species_level'] == 0:
            return np.nan
        return (row['difference'] / row['species_level'] * 100)

    summary_df['pct_change'] = summary_df.apply(calc_pct_change, axis=1)

    return summary_df


def generate_genotype_crosswalk(
    species_df: pd.DataFrame,
    family_df: pd.DataFrame,
    species_taxonomy: Optional[pd.DataFrame] = None,
    family_taxonomy: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Track how species-level groups map to family-level groups.

    Parameters
    ----------
    species_df : pd.DataFrame
        Species-level analysis
    family_df : pd.DataFrame
        Family-level analysis

    Returns
    -------
    pd.DataFrame
        Crosswalk table showing group mappings
    """
    species_tax_lookup = {}
    species_pct_lookup = {}
    if species_taxonomy is not None and not species_taxonomy.empty:
        species_tax_lookup = species_taxonomy.set_index('consensus_group')['group_majority_species'].to_dict()
        if 'group_majority_pct' in species_taxonomy.columns:
            species_pct_lookup = species_taxonomy.set_index('consensus_group')['group_majority_pct'].to_dict()

    family_tax_lookup = {}
    family_pct_lookup = {}
    if family_taxonomy is not None and not family_taxonomy.empty:
        family_tax_lookup = family_taxonomy.set_index('consensus_group')['group_majority_species'].to_dict()
        if 'group_majority_pct' in family_taxonomy.columns:
            family_pct_lookup = family_taxonomy.set_index('consensus_group')['group_majority_pct'].to_dict()

    # Merge on processid to connect samples across analyses
    merged = species_df[['processid', 'consensus_group', 'species']].merge(
        family_df[['processid', 'consensus_group']],
        on='processid',
        suffixes=('_species', '_family'),
        how='inner'
    )

    # Count mappings
    crosswalk_data = []

    for sp_group in merged['consensus_group_species'].unique():
        sp_group_samples = merged[merged['consensus_group_species'] == sp_group]
        sp_group_size = len(sp_group_samples)

        # Count how many samples from this species group go to each family group
        family_mapping = sp_group_samples['consensus_group_family'].value_counts()

        for fam_group, n_mapped in family_mapping.items():
            fam_group_size = len(merged[merged['consensus_group_family'] == fam_group])
            pct_of_sp_group = (n_mapped / sp_group_size * 100) if sp_group_size > 0 else 0
            pct_of_fam_group = (n_mapped / fam_group_size * 100) if fam_group_size > 0 else 0

            # Determine mapping type
            if pct_of_sp_group > 50:
                mapping_type = 'majority'
            elif pct_of_sp_group >= 5:
                mapping_type = 'minor'
            elif n_mapped <= 2:
                mapping_type = 'singleton'
            else:
                mapping_type = 'contamination'

            sp_majority_species = species_tax_lookup.get(sp_group)
            fam_majority_species = family_tax_lookup.get(fam_group)

            crosswalk_data.append({
                'species_group': sp_group,
                'species_group_size': sp_group_size,
                'species_group_species': sp_majority_species,
                'species_group_majority_pct': species_pct_lookup.get(sp_group),
                'family_group': fam_group,
                'family_group_size': fam_group_size,
                'family_group_species': fam_majority_species,
                'family_group_majority_pct': family_pct_lookup.get(fam_group),
                'n_samples_mapped': n_mapped,
                'pct_of_species_group': pct_of_sp_group,
                'pct_of_family_group': pct_of_fam_group,
                'mapping_type': mapping_type,
                'species_match': (
                    sp_majority_species == fam_majority_species
                    if sp_majority_species and fam_majority_species
                    else None
                )
            })

    crosswalk_df = pd.DataFrame(crosswalk_data)

    # Sort by species group and mapping count
    crosswalk_df = crosswalk_df.sort_values(
        ['species_group', 'n_samples_mapped'],
        ascending=[True, False]
    )

    return crosswalk_df


def generate_sample_reassignments(
    species_df: pd.DataFrame,
    family_df: pd.DataFrame,
    species_taxonomy: Optional[pd.DataFrame] = None,
    family_taxonomy: Optional[pd.DataFrame] = None,
    family_identity_map: Optional[pd.DataFrame] = None,
    majority_threshold: float = 0.7
) -> pd.DataFrame:
    """
    Generate sample-level reassignment table (Table S4 equivalent).

    Parameters
    ----------
    species_df : pd.DataFrame
        Species-level analysis
    family_df : pd.DataFrame
        Family-level analysis
    majority_threshold : float
        Threshold for flagging misidentifications

    Returns
    -------
    pd.DataFrame
        Sample reassignment table
    """
    species_df = species_df.copy()
    family_df = family_df.copy()

    # Fallback: derive majority pct from majority_fraction if present
    for df in (species_df, family_df):
        if 'majority_fraction' in df.columns and 'group_majority_pct' not in df.columns:
            pct = pd.to_numeric(df['majority_fraction'], errors='coerce')
            df['group_majority_pct'] = pct * 100 if pct.max() <= 1 else pct

    # Attach majority species/pct from taxonomy tables when missing in annotated outputs
    if species_taxonomy is not None and not species_taxonomy.empty:
        species_df = species_df.merge(
            species_taxonomy,
            on='consensus_group',
            how='left'
        )

    if family_taxonomy is not None and not family_taxonomy.empty:
        family_df = family_df.merge(
            family_taxonomy,
            on='consensus_group',
            how='left'
        )

    # Select relevant columns from each analysis
    species_cols = ['processid', 'species', 'assigned_sp', 'consensus_group', 'identity',
                    'group_majority_species', 'group_majority_pct']
    if 'depositor_notes' in species_df.columns:
        species_cols.append('depositor_notes')

    family_cols = ['processid', 'consensus_group', 'group_majority_species',
                   'group_majority_pct', 'identity']

    # Filter to only include columns that exist
    species_cols = [c for c in species_cols if c in species_df.columns]
    family_cols = [c for c in family_cols if c in family_df.columns]

    # Merge
    merged = species_df[species_cols].merge(
        family_df[family_cols],
        on='processid',
        suffixes=('_database', '_family'),
        how='inner'
    )

    # Rename for clarity
    reassignments = merged.copy()
    if 'species_database' in reassignments.columns:
        reassignments['database_species_label'] = reassignments['species_database']
    if 'assigned_sp_database' in reassignments.columns and 'database_species_label' not in reassignments.columns:
        reassignments['database_species_label'] = reassignments['assigned_sp_database']
    if 'consensus_group_family' in reassignments.columns:
        reassignments['family_consensus_group'] = reassignments['consensus_group_family']
    if 'consensus_group_database' in reassignments.columns:
        reassignments['species_consensus_group'] = reassignments['consensus_group_database']

    # Assignment identity (prefer diagnostics-derived values)
    if family_identity_map is not None and not family_identity_map.empty:
        identity_lookup = dict(zip(family_identity_map['processid'], family_identity_map['identity']))
        reassignments['assignment_identity'] = reassignments['processid'].map(identity_lookup)
    if 'identity_family' in reassignments.columns and 'assignment_identity' not in reassignments.columns:
        reassignments['assignment_identity'] = reassignments['identity_family']

    # Assess reassignment confidence
    high_majority_cutoff = max(90, majority_threshold * 100)
    medium_majority_cutoff = max(80, majority_threshold * 100 * 0.9)

    def assess_confidence(row):
        if 'group_majority_pct' not in row or pd.isna(row['group_majority_pct']):
            return 'unknown'

        majority_pct = row['group_majority_pct']
        identity = row.get('assignment_identity', 95)

        # High confidence
        if identity >= 95 and majority_pct >= high_majority_cutoff:
            return 'high'
        elif identity >= 90 or majority_pct >= medium_majority_cutoff:
            return 'medium'
        else:
            return 'low'

    reassignments['reassignment_confidence'] = reassignments.apply(assess_confidence, axis=1)

    # Suggested species (based on family group majority)
    if 'group_majority_species' in reassignments.columns and 'database_species_label' in reassignments.columns:
        reassignments['suggested_species'] = reassignments.apply(
            lambda row: (f"{row['group_majority_species']} - CORRECT"
                        if row['database_species_label'] == row['group_majority_species']
                        else row['group_majority_species']),
            axis=1
        )

    # Select final columns
    final_cols = [
        'processid', 'database_species_label', 'species_consensus_group', 'family_consensus_group',
        'group_majority_species', 'group_majority_pct', 'assignment_identity',
        'reassignment_confidence'
    ]

    if 'depositor_notes' in reassignments.columns:
        final_cols.append('depositor_notes')

    if 'suggested_species' in reassignments.columns:
        final_cols.append('suggested_species')

    final_cols = [c for c in final_cols if c in reassignments.columns]

    return reassignments[final_cols]


def generate_methods_text(
    species_meta: Dict[str, Any],
    family_meta: Dict[str, Any],
    summary_df: pd.DataFrame,
    output_path: Path,
    majority_threshold: float = 0.7
) -> Path:
    """
    Auto-generate methods section text.

    Parameters
    ----------
    species_meta : dict
        Species-level metadata
    family_meta : dict
        Family-level metadata
    summary_df : pd.DataFrame
        Comparison summary
    output_path : Path
        Output path for methods text
    majority_threshold : float
        Threshold used for species-level assignment

    Returns
    -------
    Path
        Path to created file
    """
    # Extract key metrics
    sp_organism = species_meta.get('organism', 'SPECIES')
    fam_organism = family_meta.get('organism', 'FAMILY')

    sp_n = summary_df[summary_df['metric'] == 'total_samples']['species_level'].values[0]
    fam_n = summary_df[summary_df['metric'] == 'total_samples']['family_level'].values[0]

    sp_groups = summary_df[summary_df['metric'] == 'consensus_groups']['species_level'].values[0]
    fam_groups = summary_df[summary_df['metric'] == 'consensus_groups']['family_level'].values[0]

    sp_detected = summary_df[summary_df['metric'] == 'species_detected']['species_level'].values[0]
    fam_detected = summary_df[summary_df['metric'] == 'species_detected']['family_level'].values[0]

    mixed_groups = summary_df[summary_df['metric'] == 'mixed_species_groups']['family_level'].values[0]
    pct_mixed = (mixed_groups / fam_groups * 100) if fam_groups > 0 else 0

    methods_text = f"""### Comparative Taxonomic Analysis

To assess database quality and detect potential misidentifications, we performed
hierarchical taxonomic analysis by comparing species-level and family-level
approaches. Species-level analysis queried BOLD for all samples labeled
"{sp_organism}" (n={sp_n} samples), while family-level analysis retrieved all
{fam_organism} samples (n={fam_n} samples representing {fam_detected} species). Both analyses used
identical BOLDGenotyper parameters (majority threshold={majority_threshold:.2f}).

Species-level analysis generated {sp_groups} consensus haplotypes, with {sp_detected} species detected,
superficially suggesting {"homogeneous, high-quality data" if sp_detected == 1 else "taxonomically diverse samples"}.
However, family-level comparative analysis revealed {mixed_groups} of {fam_groups} consensus
groups ({pct_mixed:.1f}%) contained multiple species, indicating taxonomic heterogeneity
{"invisible to single-species queries" if sp_detected == 1 else "requiring broader context"}.

Samples occurring in mixed-species consensus groups were evaluated for potential
misidentification by comparing database-reported species labels against the
majority species of their genetic cluster. Samples in species-majority groups
with non-matching labels were classified as likely misidentifications. Complete
sample-level reassignments with supporting evidence are documented in the
accompanying data tables.

### Data Quality Implications

The comparative analysis approach revealed that relying solely on database-reported
species labels can mask substantial taxonomic heterogeneity. Approximately {pct_mixed:.1f}%
of genetically distinct groups in our study showed evidence of mixed-species
composition, highlighting the importance of hierarchical taxonomic validation
when mining public barcode databases.
"""

    with open(output_path, 'w') as f:
        f.write(methods_text)

    logger.info(f"  Generated methods text: {output_path}")

    return output_path


def compare_analyses(
    species_path: Path,
    family_path: Path,
    output_dir: Path,
    generate_reassignment_table: bool = True,
    majority_threshold: float = 0.7
) -> Dict[str, Any]:
    """
    Complete comparative analysis workflow.

    Parameters
    ----------
    species_path : Path
        Path to species-level analysis
    family_path : Path
        Path to family-level analysis
    output_dir : Path
        Output directory for comparison results
    generate_reassignment_table : bool
        Whether to create sample reassignment table
    majority_threshold : float
        Threshold for species assignment

    Returns
    -------
    dict
        Summary of results and file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("="*70)
    logger.info("BOLDGenotyper Comparative Analysis")
    logger.info("="*70)

    results = {}

    # Load data
    logger.info("Loading analysis results...")
    species_df, species_meta = load_analysis_results(species_path)
    family_df, family_meta = load_analysis_results(family_path)

    logger.info(f"  Species-level: {species_meta['organism']} ({species_meta['n_samples']} samples)")
    logger.info(f"  Family-level: {family_meta['organism']} ({family_meta['n_samples']} samples)")

    results['species_meta'] = species_meta
    results['family_meta'] = family_meta

    # Supporting lookups
    species_identity_stats, species_identity_map = _load_identity_data(species_df, species_meta)
    family_identity_stats, family_identity_map = _load_identity_data(family_df, family_meta)
    species_taxonomy = _load_taxonomy_table(species_meta)
    family_taxonomy = _load_taxonomy_table(family_meta)

    # Generate comparison summary
    logger.info("Generating comparison summary...")
    summary_df = generate_comparison_summary(
        species_df,
        family_df,
        species_meta,
        family_meta,
        species_identity_stats=species_identity_stats,
        family_identity_stats=family_identity_stats
    )
    summary_path = output_dir / "comparison_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    results['comparison_summary'] = summary_path

    # Generate genotype crosswalk
    logger.info("Generating genotype crosswalk...")
    try:
        crosswalk_df = generate_genotype_crosswalk(
            species_df,
            family_df,
            species_taxonomy=species_taxonomy,
            family_taxonomy=family_taxonomy
        )
        crosswalk_path = output_dir / "genotype_crosswalk.csv"
        crosswalk_df.to_csv(crosswalk_path, index=False)
        results['genotype_crosswalk'] = crosswalk_path
    except Exception as e:
        logger.warning(f"  Could not generate crosswalk: {e}")

    # Generate sample reassignments
    if generate_reassignment_table:
        logger.info("Generating sample reassignment table...")
        try:
            reassignments = generate_sample_reassignments(
                species_df,
                family_df,
                species_taxonomy=species_taxonomy,
                family_taxonomy=family_taxonomy,
                family_identity_map=family_identity_map,
                majority_threshold=majority_threshold
            )
            reassign_path = output_dir / "sample_reassignments.csv"
            reassignments.to_csv(reassign_path, index=False)
            results['sample_reassignments'] = reassign_path
        except Exception as e:
            logger.warning(f"  Could not generate reassignments: {e}")

    # Generate methods text
    logger.info("Generating methods text...")
    try:
        methods_path = output_dir / "methods_text.md"
        generate_methods_text(
            species_meta, family_meta, summary_df, methods_path,
            majority_threshold=majority_threshold
        )
        results['methods_text'] = methods_path
    except Exception as e:
        logger.warning(f"  Could not generate methods text: {e}")

    # Create README
    readme_path = output_dir / "README.md"
    readme_content = f"""# Comparative Analysis Results

## Overview

Comparison of species-level ({species_meta['organism']}) and family-level
({family_meta['organism']}) BOLDGenotyper analyses.

## Key Findings

"""

    # Add key metrics to README
    for _, row in summary_df.iterrows():
        metric = row['metric'].replace('_', ' ').title()
        readme_content += f"- **{metric}**: {row['species_level']} (species) -> {row['family_level']} (family)\n"

    readme_content += f"""

## Files

- `comparison_summary.csv`: High-level comparison metrics
- `genotype_crosswalk.csv`: Mapping of species groups to family groups
- `sample_reassignments.csv`: Sample-level reassignment recommendations
- `methods_text.md`: Ready-to-paste methods section

## Interpretation

See `methods_text.md` for a detailed interpretation and methods description
suitable for manuscript inclusion.

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info("="*70)
    logger.info("Comparative analysis complete")
    logger.info(f"  Output directory: {output_dir}")
    logger.info("="*70)

    return results


def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        prog="compare_analyses.py",
        description="Compare species-level and family-level BOLDGenotyper analyses for contamination detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two completed analyses
  python scripts/compare_analyses.py \\
      --species-level Sphyrna_lewini_output/ \\
      --family-level Sphyrnidae_output/

  # Generate sample reassignment table
  python scripts/compare_analyses.py \\
      --species-level species_analysis/ \\
      --family-level family_analysis/ \\
      --generate-reassignment-table

  # Specify custom output directory
  python scripts/compare_analyses.py \\
      --species-level sp/ --family-level fam/ \\
      --output custom_comparison/
        """
    )

    parser.add_argument(
        '--species-level',
        type=Path,
        required=True,
        help='Path to species-level analysis directory or annotated CSV'
    )

    parser.add_argument(
        '--family-level',
        type=Path,
        required=True,
        help='Path to family-level analysis directory or annotated CSV'
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=Path('comparative_analysis'),
        help='Output directory for comparison results (default: comparative_analysis/)'
    )

    parser.add_argument(
        '--generate-reassignment-table',
        action='store_true',
        help='Generate sample-level reassignment table'
    )

    parser.add_argument(
        '--majority-threshold',
        type=float,
        default=0.7,
        help='Threshold for species-level assignment (default: 0.7)'
    )

    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging verbosity (default: INFO)'
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Validate inputs
    if not args.species_level.exists():
        print(f"Error: Species-level path not found: {args.species_level}", file=sys.stderr)
        return 1

    if not args.family_level.exists():
        print(f"Error: Family-level path not found: {args.family_level}", file=sys.stderr)
        return 1

    # Run comparison
    try:
        results = compare_analyses(
            species_path=args.species_level,
            family_path=args.family_level,
            output_dir=args.output,
            generate_reassignment_table=args.generate_reassignment_table,
            majority_threshold=args.majority_threshold
        )

        print("\nComparative analysis complete!")
        print(f"  Results saved to: {args.output}")
        print("\nGenerated files:")
        for key, path in results.items():
            if isinstance(path, Path):
                print(f"  - {path.name}")

        return 0

    except Exception as e:
        logging.error(f"Comparative analysis failed: {e}", exc_info=True)
        print(f"\nError: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
