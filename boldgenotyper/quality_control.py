"""
Quality Control Module for BOLDGenotyper

This module provides automated contamination detection and quality control reporting
for genotype assignments. It identifies mixed-species consensus groups, potential
misidentifications, and data quality issues.

Key Features:
- Mixed-species group detection
- Potential misidentification flagging
- Depositor uncertainty note extraction
- Contamination heatmaps and reports
- Purity distribution analysis

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
import re
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


class QualityControlError(Exception):
    """Base exception for quality control errors."""
    pass


# Keywords for depositor uncertainty detection
UNCERTAINTY_KEYWORDS = {
    'cryptic_species': ['cryptic', 'complex', 'cf. gilberti'],
    'uncertain_id': ['uncertain', 'tentative', 'probable', 'likely', 'possibly'],
    'morphology_issues': ['juvenile', 'damaged', 'incomplete', 'poor condition'],
    'taxonomy_qualifier': ['cf.', 'aff.', 'sp.', 'spp.'],
    'mixed_sample': ['mixed', 'contamination', 'multiple']
}


def detect_depositor_uncertainty(
    note: str
) -> Tuple[bool, str, str]:
    """
    Detect uncertainty indicators in depositor notes.

    Parameters
    ----------
    note : str
        Depositor notes field

    Returns
    -------
    tuple
        (has_uncertainty: bool, flag_category: str, pattern_matched: str)
    """
    if pd.isna(note) or not note:
        return False, "", ""

    note_lower = str(note).lower()

    for category, keywords in UNCERTAINTY_KEYWORDS.items():
        for keyword in keywords:
            if keyword.lower() in note_lower:
                return True, category, keyword

    return False, "", ""


def analyze_group_contamination(
    df: pd.DataFrame,
    group_col: str = "consensus_group",
    species_col: str = "species",
    min_samples: int = 3
) -> pd.DataFrame:
    """
    Analyze contamination levels in each consensus group.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    group_col : str
        Column name for consensus groups
    species_col : str
        Column name for species labels
    min_samples : int
        Minimum samples required to analyze a group

    Returns
    -------
    pd.DataFrame
        Summary table with contamination metrics per group
    """
    results = []

    for group_name, group_df in df.groupby(group_col):
        n_samples = len(group_df)

        # Skip very small groups
        if n_samples < min_samples:
            continue

        # Count species
        species_counts = group_df[species_col].value_counts()
        n_species = len(species_counts)

        # Primary species
        primary_species = species_counts.index[0] if len(species_counts) > 0 else "Unknown"
        primary_count = species_counts.iloc[0] if len(species_counts) > 0 else 0
        primary_pct = (primary_count / n_samples * 100) if n_samples > 0 else 0

        # Contaminating species
        if n_species > 1:
            contaminating_species = []
            contamination_count = 0
            contamination_breakdown = []

            for i, (sp, count) in enumerate(species_counts.items()):
                if i > 0:  # Skip primary species
                    contaminating_species.append(sp)
                    contamination_count += count
                    pct = count / n_samples * 100
                    contamination_breakdown.append(f"{sp}: {count} ({pct:.1f}%)")

            contaminating_sp_str = ", ".join(contaminating_species)
            contamination_pct = (contamination_count / n_samples * 100)
            breakdown_str = ", ".join(contamination_breakdown)
        else:
            contaminating_sp_str = ""
            contamination_count = 0
            contamination_pct = 0.0
            breakdown_str = ""

        # Determine flag
        if n_species == 1:
            flag = "CLEAN"
            confidence = "high"
        elif primary_pct >= 90:
            flag = "CONTAMINATED"
            confidence = "high"
        elif primary_pct >= 70:
            flag = "CONTAMINATED"
            confidence = "medium"
        else:
            flag = "AMBIGUOUS"
            confidence = "low"

        results.append({
            'consensus_group': group_name,
            'n_samples': n_samples,
            'n_species': n_species,
            'primary_species': primary_species,
            'primary_pct': primary_pct,
            'contaminating_species': contaminating_sp_str,
            'contamination_count': contamination_count,
            'contamination_pct': contamination_pct,
            'contamination_breakdown': breakdown_str,
            'flags': flag,
            'confidence_level': confidence
        })

    return pd.DataFrame(results)


def add_contamination_columns(
    df: pd.DataFrame,
    group_col: str = "consensus_group",
    species_col: str = "species",
    notes_col: str = "notes",
    majority_threshold: float = 0.7
) -> pd.DataFrame:
    """
    Add contamination analysis columns to the main annotated dataframe.

    This enhances the main output CSV with quality control information:
    - group_majority_species: Most common species in the group
    - group_majority_pct: Percentage of majority species
    - group_species_count: Number of unique species in group
    - matches_group_majority: Does sample match group majority?
    - potential_misidentification: Likely mislabeled sample?
    - depositor_uncertainty_flag: Has notes indicating uncertainty?
    - depositor_notes: Full notes field
    - misidentification_confidence: Confidence level

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    group_col : str
        Column name for consensus groups
    species_col : str
        Column name for species labels
    notes_col : str
        Column name for depositor notes
    majority_threshold : float
        Threshold for flagging misidentifications (default: 0.7)

    Returns
    -------
    pd.DataFrame
        Enhanced dataframe with contamination columns
    """
    df = df.copy()

    # Calculate group-level statistics
    group_stats = df.groupby(group_col).agg({
        species_col: lambda x: x.mode()[0] if len(x.mode()) > 0 else 'Unknown'
    }).rename(columns={species_col: 'group_majority_species'})

    # Calculate majority percentage
    def calc_majority_pct(group_df):
        if len(group_df) == 0:
            return 0
        species_counts = group_df[species_col].value_counts()
        if len(species_counts) == 0:
            return 0
        return (species_counts.iloc[0] / len(group_df)) * 100

    group_majority_pct = df.groupby(group_col, group_keys=False).apply(calc_majority_pct, include_groups=False)
    group_majority_pct.name = 'group_majority_pct'

    # Species count per group
    group_species_count = df.groupby(group_col)[species_col].nunique()
    group_species_count.name = 'group_species_count'

    # Merge back to main dataframe
    df = df.merge(group_stats, left_on=group_col, right_index=True, how='left')
    df = df.merge(group_majority_pct, left_on=group_col, right_index=True, how='left')
    df = df.merge(group_species_count, left_on=group_col, right_index=True, how='left')

    # Calculate matches_group_majority
    df['matches_group_majority'] = df[species_col] == df['group_majority_species']

    # Potential misidentification
    df['potential_misidentification'] = (
        (~df['matches_group_majority']) &
        (df['group_majority_pct'] > (majority_threshold * 100))
    )

    # Depositor uncertainty analysis
    if notes_col in df.columns:
        uncertainty_results = df[notes_col].apply(detect_depositor_uncertainty)
        df['depositor_uncertainty_flag'] = uncertainty_results.apply(lambda x: x[0])
        df['uncertainty_category'] = uncertainty_results.apply(lambda x: x[1])
        df['uncertainty_keyword'] = uncertainty_results.apply(lambda x: x[2])
        df['depositor_notes'] = df[notes_col]
    else:
        df['depositor_uncertainty_flag'] = False
        df['uncertainty_category'] = ""
        df['uncertainty_keyword'] = ""
        df['depositor_notes'] = ""

    # Misidentification confidence
    def assess_confidence(row):
        if not row['potential_misidentification']:
            return 'none'

        # Check if we have identity score
        has_identity = 'identity' in row and not pd.isna(row['identity'])
        identity = float(row['identity']) if has_identity else 95.0

        # High confidence: high identity, clear majority, no uncertainty
        if identity >= 95 and row['group_majority_pct'] > 90 and not row['depositor_uncertainty_flag']:
            return 'high'

        # Medium confidence: moderate identity or moderate majority
        elif identity >= 90 or row['group_majority_pct'] >= 80:
            return 'medium'

        # Low confidence: lower identity or weaker majority
        else:
            return 'low'

    df['misidentification_confidence'] = df.apply(assess_confidence, axis=1)

    return df


def generate_depositor_flags_summary(
    df: pd.DataFrame,
    group_col: str = "consensus_group",
    species_col: str = "species",
    processid_col: str = "processid"
) -> pd.DataFrame:
    """
    Generate summary of samples with depositor uncertainty flags.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with contamination columns
    group_col : str
        Column name for consensus groups
    species_col : str
        Column name for species labels
    processid_col : str
        Column name for process IDs

    Returns
    -------
    pd.DataFrame
        Summary of flagged samples
    """
    flagged = df[df['depositor_uncertainty_flag']].copy()

    if len(flagged) == 0:
        return pd.DataFrame(columns=[
            'processid', 'genotype', 'species_label', 'depositor_note',
            'flag_category', 'pattern_matched'
        ])

    summary = flagged[[
        processid_col, group_col, species_col, 'depositor_notes',
        'uncertainty_category', 'uncertainty_keyword'
    ]].copy()

    summary.columns = [
        'processid', 'genotype', 'species_label', 'depositor_note',
        'flag_category', 'pattern_matched'
    ]

    return summary


def create_contamination_heatmap(
    df: pd.DataFrame,
    output_path: Path,
    group_col: str = "consensus_group",
    species_col: str = "species",
    min_samples: int = 3
) -> None:
    """
    Create a heatmap showing species composition of consensus groups.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    output_path : Path
        Output path for PDF
    group_col : str
        Column name for consensus groups
    species_col : str
        Column name for species labels
    min_samples : int
        Minimum samples to include a group
    """
    # Filter to groups with sufficient samples
    group_sizes = df[group_col].value_counts()
    valid_groups = group_sizes[group_sizes >= min_samples].index
    df_filtered = df[df[group_col].isin(valid_groups)].copy()

    if len(df_filtered) == 0:
        logger.warning("No groups with sufficient samples for heatmap")
        return

    # Create contingency table
    ct = pd.crosstab(df_filtered[group_col], df_filtered[species_col])

    # Sort by group size
    group_order = df_filtered[group_col].value_counts().index
    ct = ct.loc[group_order]

    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, len(ct) * 0.4)))

    # Create heatmap
    sns.heatmap(
        ct,
        annot=True,
        fmt='d',
        cmap='YlOrRd',
        cbar_kws={'label': 'Sample Count'},
        ax=ax,
        linewidths=0.5,
        linecolor='gray'
    )

    ax.set_xlabel('Species', fontsize=12, fontweight='bold')
    ax.set_ylabel('Consensus Group', fontsize=12, fontweight='bold')
    ax.set_title('Species Contamination Heatmap', fontsize=14, fontweight='bold', pad=20)

    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Contamination heatmap saved to {output_path}")


def create_purity_distribution(
    contamination_summary: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Create histogram of majority fractions (group purity).

    Parameters
    ----------
    contamination_summary : pd.DataFrame
        Output from analyze_group_contamination
    output_path : Path
        Output path for PDF
    """
    if len(contamination_summary) == 0:
        logger.warning("No contamination data for purity distribution")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    # Define purity bins and colors
    bins = [0, 70, 90, 99, 100]
    colors = ['#d73027', '#fc8d59', '#fee08b', '#91cf60']
    labels = ['<70% (Ambiguous)', '70-90% (Moderate)', '90-99% (Minor)', '100% (Pure)']

    # Create histogram
    purity = contamination_summary['primary_pct'].values

    # Create histogram with default color first
    counts, bin_edges, patches = ax.hist(purity, bins=bins, edgecolor='black', linewidth=1.2)

    # Color each bar according to its bin
    for patch, color in zip(patches, colors):
        patch.set_facecolor(color)

    # Add threshold line
    ax.axvline(70, color='red', linestyle='--', linewidth=2, label='70% threshold')

    # Labels and styling
    ax.set_xlabel('Majority Fraction (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Groups', fontsize=12, fontweight='bold')
    ax.set_title('Purity Distribution of Consensus Groups', fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper left')
    ax.grid(axis='y', alpha=0.3)

    # Add statistics
    stats_text = f"Mean: {purity.mean():.1f}%\nMedian: {np.median(purity):.1f}%"
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Purity distribution saved to {output_path}")


def generate_quality_control_report(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    group_col: str = "consensus_group",
    species_col: str = "species"
) -> Dict[str, Any]:
    """
    Generate comprehensive quality control report.

    This creates:
    - mixed_species_summary.csv
    - contamination_heatmap.pdf
    - depositor_flags_summary.csv
    - potential_misidentifications.csv
    - purity_distribution.pdf
    - README.md

    Parameters
    ----------
    df : pd.DataFrame
        Enhanced annotated dataframe with contamination columns
    output_dir : Path
        Output directory for QC files
    organism : str
        Organism name
    group_col : str
        Column name for consensus groups
    species_col : str
        Column name for species labels

    Returns
    -------
    dict
        Summary statistics and file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # 1. Analyze contamination
    contamination_summary = analyze_group_contamination(df, group_col, species_col)
    contamination_path = output_dir / "mixed_species_summary.csv"
    contamination_summary.to_csv(contamination_path, index=False)
    results['mixed_species_summary'] = contamination_path

    # 2. Depositor flags
    depositor_flags = generate_depositor_flags_summary(df, group_col, species_col)
    flags_path = output_dir / "depositor_flags_summary.csv"
    depositor_flags.to_csv(flags_path, index=False)
    results['depositor_flags'] = flags_path

    # 3. Potential misidentifications
    misidentified = df[df['potential_misidentification']].copy()
    if len(misidentified) > 0:
        misid_cols = [
            'processid', species_col, group_col, 'group_majority_species',
            'group_majority_pct', 'misidentification_confidence',
            'depositor_uncertainty_flag', 'depositor_notes'
        ]
        misid_cols = [c for c in misid_cols if c in misidentified.columns]
        misid_path = output_dir / "potential_misidentifications.csv"
        misidentified[misid_cols].to_csv(misid_path, index=False)
        results['potential_misidentifications'] = misid_path

    # 4. Contamination heatmap
    heatmap_path = output_dir / "contamination_heatmap.pdf"
    try:
        create_contamination_heatmap(df, heatmap_path, group_col, species_col)
        results['contamination_heatmap'] = heatmap_path
    except Exception as e:
        logger.warning(f"Could not create contamination heatmap: {e}")

    # 5. Purity distribution
    purity_path = output_dir / "purity_distribution.pdf"
    try:
        create_purity_distribution(contamination_summary, purity_path)
        results['purity_distribution'] = purity_path
    except Exception as e:
        logger.warning(f"Could not create purity distribution: {e}")

    # 6. Summary statistics
    n_groups = len(contamination_summary)
    n_mixed = len(contamination_summary[contamination_summary['n_species'] > 1])
    n_pure = len(contamination_summary[contamination_summary['n_species'] == 1])
    pct_contaminated = (n_mixed / n_groups * 100) if n_groups > 0 else 0
    n_flagged_samples = len(df[df['depositor_uncertainty_flag']])
    n_potential_misid = len(df[df['potential_misidentification']])

    results['summary'] = {
        'total_groups': n_groups,
        'pure_groups': n_pure,
        'mixed_groups': n_mixed,
        'pct_contaminated': pct_contaminated,
        'flagged_samples': n_flagged_samples,
        'potential_misidentifications': n_potential_misid
    }

    # 7. README
    readme_path = output_dir / "README.md"
    readme_content = f"""# Quality Control Report for {organism}

## Overview

This directory contains automated quality control analysis identifying potential
contamination and misidentifications in the genotype assignments.

## Files

- **mixed_species_summary.csv**: Summary of contamination in each consensus group
- **contamination_heatmap.pdf**: Visual overview of species × groups
- **depositor_flags_summary.csv**: Samples with uncertainty notes
- **potential_misidentifications.csv**: Flagged samples likely mislabeled
- **purity_distribution.pdf**: Distribution of majority fractions

## Summary Statistics

- **Total consensus groups analyzed**: {n_groups}
- **Pure groups (single species)**: {n_pure} ({n_pure/n_groups*100:.1f}%)
- **Mixed groups (multiple species)**: {n_mixed} ({pct_contaminated:.1f}%)
- **Samples with depositor uncertainty flags**: {n_flagged_samples}
- **Potential misidentifications**: {n_potential_misid}

## Interpretation

### Group Purity Flags

- **CLEAN**: 100% single species - high confidence homogeneous group
- **CONTAMINATED**: <100% but >70% majority - minor contamination present
- **AMBIGUOUS**: <70% majority - unclear species assignment, genus-level appropriate

### Misidentification Confidence

- **high**: >95% identity, clear majority (>90%), no depositor uncertainty
- **medium**: 90-95% identity, moderate majority (70-90%), or has uncertainty notes
- **low**: <90% identity, weak majority (<70%), or flagged as ambiguous

## Recommended Actions

1. Review samples in `potential_misidentifications.csv`
2. Cross-reference with `depositor_flags_summary.csv` for known issues
3. Consider hierarchical (family-level) analysis for contamination context
4. Manually review ambiguous groups with <70% majority

## Methods Integration

Include quality control metrics in your manuscript's results section.
See contamination_summary.csv for group-level statistics suitable for
supplementary tables.
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info(f"Quality control report generated in {output_dir}")

    return results


def print_quality_alert(qc_summary: Dict[str, Any]) -> None:
    """
    Print quality control alerts to console.

    Parameters
    ----------
    qc_summary : dict
        Summary dictionary from generate_quality_control_report
    """
    summary = qc_summary.get('summary', {})

    n_mixed = summary.get('mixed_groups', 0)
    total_groups = summary.get('total_groups', 0)
    n_flagged = summary.get('flagged_samples', 0)
    n_misid = summary.get('potential_misidentifications', 0)

    if n_mixed > 0 or n_flagged > 0 or n_misid > 0:
        print("\n" + "="*70)
        print("⚠️  QUALITY CONTROL ALERTS")
        print("="*70)

        if n_mixed > 0:
            pct = (n_mixed / total_groups * 100) if total_groups > 0 else 0
            print(f"\n🔬 {n_mixed} of {total_groups} consensus groups contain multiple species ({pct:.1f}%)")
            print("   → See: quality_control/mixed_species_summary.csv")

        if n_misid > 0:
            print(f"\n🚨 {n_misid} samples flagged as potential misidentifications")
            print("   → See: quality_control/potential_misidentifications.csv")

        if n_flagged > 0:
            print(f"\n📝 {n_flagged} samples have depositor notes indicating uncertainty")
            print("   → See: quality_control/depositor_flags_summary.csv")

        print("\n💡 TIP: Consider running family-level analysis for quality control")
        print("   Run: boldgenotyper-compare --help")
        print("="*70 + "\n")
    else:
        print("\n✓ Quality control: No significant contamination detected")
