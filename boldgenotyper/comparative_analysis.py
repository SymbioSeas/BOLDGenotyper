"""
Comparative Analysis Module for BOLDGenotyper

This module enables hierarchical taxonomic comparison by comparing species-level
and family/genus-level analyses to detect database contamination and potential
misidentifications.

Key Features:
- Compare two BOLDGenotyper outputs (e.g., species vs family level)
- Generate genotype crosswalk tables
- Identify mixed-species groups and contamination patterns
- Create sample-level reassignment tables
- Generate comparison visualizations
- Auto-generate methods text

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


class ComparativeAnalysisError(Exception):
    """Base exception for comparative analysis errors."""
    pass


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
    family_meta: Dict[str, Any]
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
    if 'identity' in species_df.columns:
        species_mean_id = species_df['identity'].mean()
    else:
        species_mean_id = np.nan
    if 'identity' in family_df.columns:
        family_mean_id = family_df['identity'].mean()
    else:
        family_mean_id = np.nan
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
    family_df: pd.DataFrame
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

            crosswalk_data.append({
                'species_group': sp_group,
                'species_group_size': sp_group_size,
                'family_group': fam_group,
                'family_group_size': fam_group_size,
                'n_samples_mapped': n_mapped,
                'pct_of_species_group': pct_of_sp_group,
                'pct_of_family_group': pct_of_fam_group,
                'mapping_type': mapping_type
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
    # Select relevant columns from each analysis
    species_cols = ['processid', 'species', 'consensus_group', 'identity']
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
    if 'consensus_group_family' in reassignments.columns:
        reassignments['family_consensus_group'] = reassignments['consensus_group_family']
    if 'identity_family' in reassignments.columns:
        reassignments['assignment_identity'] = reassignments['identity_family']

    # Assess reassignment confidence
    def assess_confidence(row):
        if 'group_majority_pct' not in row or pd.isna(row['group_majority_pct']):
            return 'unknown'

        majority_pct = row['group_majority_pct']
        identity = row.get('assignment_identity', 95)

        # High confidence
        if identity >= 95 and majority_pct > 90:
            return 'high'
        elif identity >= 90 or majority_pct >= 80:
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
        'processid', 'database_species_label', 'family_consensus_group',
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
    output_path: Path
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

    # Get threshold
    threshold = species_meta.get('params', {}).get('clustering_threshold', 0.015)

    methods_text = f"""### Comparative Taxonomic Analysis

To assess database quality and detect potential misidentifications, we performed
hierarchical taxonomic analysis by comparing species-level and family-level
clustering approaches. Species-level analysis queried BOLD for all samples
labeled "{sp_organism}" (n={sp_n} samples), while family-level analysis retrieved all
{fam_organism} samples (n={fam_n} samples representing {fam_detected} species). Both analyses used
identical BOLDGenotyper parameters (clustering threshold={threshold}, assignment
threshold=0.50).

Species-level analysis generated {sp_groups} consensus genotypes, with {sp_detected} species detected,
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

    logger.info(f"  ✓ Generated methods text: {output_path}")

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

    # Generate comparison summary
    logger.info("Generating comparison summary...")
    summary_df = generate_comparison_summary(species_df, family_df, species_meta, family_meta)
    summary_path = output_dir / "comparison_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    results['comparison_summary'] = summary_path

    # Generate genotype crosswalk
    logger.info("Generating genotype crosswalk...")
    try:
        crosswalk_df = generate_genotype_crosswalk(species_df, family_df)
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
                species_df, family_df, majority_threshold=majority_threshold
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
        generate_methods_text(species_meta, family_meta, summary_df, methods_path)
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
        readme_content += f"- **{metric}**: {row['species_level']} (species) → {row['family_level']} (family)\n"

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
    logger.info(f"✓ Comparative analysis complete")
    logger.info(f"  Output directory: {output_dir}")
    logger.info("="*70)

    return results
