"""
Quality Control Module for BOLDGenotyper

This module provides quality control for COI barcode sequences, including
sequence orientation normalization, ORF validation, dynamic filtering, and contamination
detection.

Key Features:
- COI orientation normalization (forward/reverse complement detection)
- Open reading frame (ORF) validation and quality checking
- Dynamic QC filtering with median-based thresholds
- Mixed-species group detection
- Potential misidentification flagging
- Depositor uncertainty note extraction
- Contamination heatmaps and reports
- Purity distribution analysis

Workflow:
1. Orientation Normalization: Check ORF in both orientations, correct if needed
2. ORF Validation: Check coverage and internal stop codons
3. Dynamic QC: Apply absolute and median-based length/quality filters
4. Contamination Detection: Identify mixed-species groups

Author: Steph Smith (symbioseas@outlook.com)
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
    'cryptic_species': ['cryptic', 'complex', 'novel', 'undescribed'],
    'uncertain_id': ['uncertain', 'tentative', 'probable', 'likely', 'possibly', 'unknown', 'unidentified'],
    'morphology_issues': ['juvenile', 'damaged', 'incomplete', 'poor condition', 'pup', 'decomposed', 'larva', 'nymph', 'egg'],
    'taxonomy_qualifier': ['cf.', 'aff.', 'sp.', 'spp.'],
    'mixed_sample': ['mixed', 'contamination', 'multiple', 'contaminated']
}


# ============================================================================
# COI Orientation Normalization and ORF Validation
# ============================================================================

def apply_orientation_normalization(
    sequences_dict: Dict[str, str],
    genetic_code: int = 2,
    min_orf_coverage: float = 0.7,
    max_internal_stops: int = 2
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """
    Normalize sequence orientation and validate ORF quality.

    Checks each sequence for proper COI orientation by analyzing ORF in both
    forward and reverse complement orientations. Sequences are automatically
    corrected to forward orientation if needed.

    Parameters
    ----------
    sequences_dict : Dict[str, str]
        Dictionary mapping processid -> sequence
    genetic_code : int, optional
        NCBI genetic code table (default: 2 = vertebrate mitochondrial)
    min_orf_coverage : float, optional
        Minimum ORF coverage required (default: 0.7 = 70%)
    max_internal_stops : int, optional
        Maximum internal stop codons allowed (default: 2)

    Returns
    -------
    Tuple[Dict[str, str], pd.DataFrame]
        - corrected_sequences: Dictionary with orientation-corrected sequences
        - orf_stats: DataFrame with ORF validation results per sequence

    Examples
    --------
    >>> sequences = {"sample1": "ATGCCC...", "sample2": "GGGCAT..."}
    >>> corrected, stats = apply_orientation_normalization(sequences)
    >>> stats[['processid', 'orientation', 'orf_valid', 'orf_coverage']]
    """
    from boldgenotyper.utils import check_orf_quality

    corrected_sequences = {}
    orf_results = []

    for processid, sequence in sequences_dict.items():
        if not sequence or pd.isna(sequence):
            # Handle missing sequences
            orf_results.append({
                'processid': processid,
                'orientation': 'unknown',
                'needs_revcomp': False,
                'orf_valid': False,
                'orf_coverage': 0.0,
                'internal_stops': 0,
                'failure_reasons': 'No sequence data'
            })
            corrected_sequences[processid] = sequence
            continue

        # Clean sequence: remove gaps and non-ATGCN characters
        cleaned_sequence = sequence.replace('-', '').replace('.', '').upper()
        # Keep only valid nucleotides
        cleaned_sequence = ''.join(c for c in cleaned_sequence if c in 'ATGCN')

        if not cleaned_sequence:
            # Sequence was all gaps or invalid characters
            orf_results.append({
                'processid': processid,
                'orientation': 'unknown',
                'needs_revcomp': False,
                'orf_valid': False,
                'orf_coverage': 0.0,
                'internal_stops': 0,
                'failure_reasons': 'Sequence contains only gaps or invalid characters'
            })
            corrected_sequences[processid] = cleaned_sequence
            continue

        # Check ORF quality and orientation on cleaned sequence
        orf_check = check_orf_quality(
            cleaned_sequence,
            genetic_code=genetic_code,
            min_coverage=min_orf_coverage,
            max_internal_stops=max_internal_stops
        )

        # Store results (including reading frame)
        orf_results.append({
            'processid': processid,
            'orientation': orf_check['orientation'],
            'frame': orf_check.get('frame', 0),  # Reading frame offset
            'needs_revcomp': orf_check['needs_revcomp'],
            'orf_valid': orf_check['is_valid_orf'],
            'orf_coverage': orf_check['orf_coverage'],
            'internal_stops': orf_check['internal_stops'],
            'failure_reasons': '; '.join(orf_check['failure_reasons']) if orf_check['failure_reasons'] else ''
        })

        # Store corrected sequence
        corrected_sequences[processid] = orf_check['corrected_sequence']

    # Convert to DataFrame
    orf_stats_df = pd.DataFrame(orf_results)

    # Log summary
    n_total = len(orf_results)
    n_reversed = orf_stats_df['needs_revcomp'].sum()
    n_invalid = (~orf_stats_df['orf_valid']).sum()

    # Count frame distribution
    frame_counts = orf_stats_df['frame'].value_counts().sort_index()

    logger.info(f"Orientation normalization complete:")
    logger.info(f"  Total sequences: {n_total}")
    logger.info(f"  Reverse complemented: {n_reversed} ({n_reversed/n_total*100:.1f}%)")
    logger.info(f"  Reading frame distribution:")
    for frame, count in frame_counts.items():
        logger.info(f"    Frame +{frame}: {count} ({count/n_total*100:.1f}%)")
    logger.info(f"  Invalid ORF: {n_invalid} ({n_invalid/n_total*100:.1f}%)")

    if n_invalid > 0:
        logger.warning(
            f"  {n_invalid} sequences failed ORF validation. "
            f"These may be contamination, NUMTs, or sequencing errors."
        )

    return corrected_sequences, orf_stats_df


# ============================================================================
# Dynamic Quality Control Filtering
# ============================================================================

def apply_dynamic_qc_filters(
    df: pd.DataFrame,
    sequences_dict: Dict[str, str],
    orf_stats_df: pd.DataFrame = None,
    min_raw_length_abs: int = 200,
    min_raw_length_frac_of_median: float = 0.7,
    max_raw_N_fraction: float = 0.05,
    require_valid_orf: bool = True,
    processid_col: str = 'processid'
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Apply dynamic quality control filters with median-based thresholds.

    Implements three-tier filtering:
    1. Absolute thresholds (minimum length, max N content)
    2. Median-based relative thresholds (adaptive to dataset)
    3. ORF validation (excludes NUMTs, contamination, pseudogenes)

    Parameters
    ----------
    df : pd.DataFrame
        Metadata dataframe with processid column
    sequences_dict : Dict[str, str]
        Dictionary mapping processid -> sequence
    orf_stats_df : pd.DataFrame, optional
        ORF validation results from apply_orientation_normalization
    min_raw_length_abs : int, optional
        Absolute minimum sequence length (default: 200 bp)
    min_raw_length_frac_of_median : float, optional
        Minimum length as fraction of median (default: 0.7 = 70%)
    max_raw_N_fraction : float, optional
        Maximum fraction of N bases (default: 0.05 = 5%)
    require_valid_orf : bool, optional
        Require valid ORF to pass QC (default: True). Excludes NUMTs,
        contamination, and pseudogenes.
    processid_col : str, optional
        Name of processid column (default: 'processid')

    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        - filtered_df: DataFrame with QC columns and failed samples removed
        - qc_stats: Dictionary with QC statistics

    Notes
    -----
    QC columns added to DataFrame:
    - raw_length: Sequence length in bp
    - raw_N_count: Number of N bases
    - raw_N_fraction: Fraction of N bases
    - qc_pass_abs: Passes absolute thresholds
    - qc_pass_median: Passes median-based thresholds
    - qc_pass_orf: Passes ORF validation (if required)
    - qc_pass: Passes all QC filters
    - qc_fail_reason: Reason for failure (if failed)
    """
    # Calculate sequence metrics
    qc_data = []

    for processid in df[processid_col]:
        sequence = sequences_dict.get(processid, '')

        if not sequence or pd.isna(sequence):
            qc_data.append({
                'processid': processid,
                'raw_length': 0,
                'raw_N_count': 0,
                'raw_N_fraction': 0.0
            })
        else:
            sequence = str(sequence).upper()
            seq_length = len(sequence)
            n_count = sequence.count('N')
            n_fraction = n_count / seq_length if seq_length > 0 else 0.0

            qc_data.append({
                'processid': processid,
                'raw_length': seq_length,
                'raw_N_count': n_count,
                'raw_N_fraction': n_fraction
            })

    qc_df = pd.DataFrame(qc_data)

    # Calculate median length for adaptive threshold
    median_length = qc_df['raw_length'].median()
    min_length_adaptive = median_length * min_raw_length_frac_of_median

    logger.info(f"Dynamic QC thresholds:")
    logger.info(f"  Median sequence length: {median_length:.0f} bp")
    logger.info(f"  Absolute minimum: {min_raw_length_abs} bp")
    logger.info(f"  Adaptive minimum: {min_length_adaptive:.0f} bp ({min_raw_length_frac_of_median:.0%} of median)")
    logger.info(f"  Max N fraction: {max_raw_N_fraction:.1%}")

    # Apply absolute thresholds
    qc_df['qc_pass_abs'] = (
        (qc_df['raw_length'] >= min_raw_length_abs) &
        (qc_df['raw_N_fraction'] <= max_raw_N_fraction)
    )

    # Apply median-based thresholds
    qc_df['qc_pass_median'] = (qc_df['raw_length'] >= min_length_adaptive)

    # Apply ORF validation filter if provided
    if orf_stats_df is not None and require_valid_orf:
        # Merge ORF validation results
        qc_df = qc_df.merge(
            orf_stats_df[['processid', 'orf_valid']],
            on='processid',
            how='left'
        )
        # Default to False if ORF validation missing
        qc_df['orf_valid'] = qc_df['orf_valid'].fillna(False)
        qc_df['qc_pass_orf'] = qc_df['orf_valid']
    else:
        # If ORF validation not required, all pass this criterion
        qc_df['qc_pass_orf'] = True

    # Combined QC pass
    qc_df['qc_pass'] = (
        qc_df['qc_pass_abs'] &
        qc_df['qc_pass_median'] &
        qc_df['qc_pass_orf']
    )

    # Determine failure reasons
    def get_fail_reason(row):
        reasons = []
        if row['raw_length'] < min_raw_length_abs:
            reasons.append(f"Length {row['raw_length']}bp < {min_raw_length_abs}bp")
        if row['raw_length'] < min_length_adaptive:
            reasons.append(f"Length {row['raw_length']}bp < {min_length_adaptive:.0f}bp (adaptive)")
        if row['raw_N_fraction'] > max_raw_N_fraction:
            reasons.append(f"N content {row['raw_N_fraction']:.1%} > {max_raw_N_fraction:.1%}")
        if require_valid_orf and orf_stats_df is not None and not row['qc_pass_orf']:
            reasons.append("Invalid ORF (likely NUMT/contamination/pseudogene)")
        return '; '.join(reasons) if reasons else ''

    qc_df['qc_fail_reason'] = qc_df.apply(get_fail_reason, axis=1)

    # Merge QC columns with input dataframe
    df_with_qc = df.merge(qc_df, on=processid_col, how='left')

    # Calculate statistics
    n_total = len(qc_df)
    n_pass = qc_df['qc_pass'].sum()
    n_fail_abs = (~qc_df['qc_pass_abs']).sum()
    n_fail_median = (~qc_df['qc_pass_median']).sum()
    n_fail_orf = (~qc_df['qc_pass_orf']).sum() if require_valid_orf and orf_stats_df is not None else 0
    n_fail_total = n_total - n_pass

    qc_stats = {
        'n_total': n_total,
        'n_pass': n_pass,
        'n_fail': n_fail_total,
        'n_fail_abs': n_fail_abs,
        'n_fail_median': n_fail_median,
        'n_fail_orf': n_fail_orf,
        'pass_rate': n_pass / n_total if n_total > 0 else 0.0,
        'median_length': median_length,
        'min_length_adaptive': min_length_adaptive
    }

    logger.info(f"QC filtering results:")
    logger.info(f"  Total samples: {n_total}")
    logger.info(f"  Passed QC: {n_pass} ({qc_stats['pass_rate']*100:.1f}%)")
    logger.info(f"  Failed QC: {n_fail_total} ({(1-qc_stats['pass_rate'])*100:.1f}%)")
    logger.info(f"    - Failed absolute thresholds: {n_fail_abs}")
    logger.info(f"    - Failed median-based thresholds: {n_fail_median}")
    if require_valid_orf and orf_stats_df is not None:
        logger.info(f"    - Failed ORF validation: {n_fail_orf}")
        if n_fail_orf > 0:
            logger.warning(
                f"  {n_fail_orf} sequences excluded as likely NUMTs, contamination, or pseudogenes"
            )

    # Filter to only passing samples
    filtered_df = df_with_qc[df_with_qc['qc_pass']].copy()

    logger.info(f"Retained {len(filtered_df)} samples after QC filtering")

    return filtered_df, qc_stats


# ============================================================================
# Depositor Uncertainty and Contamination Detection
# ============================================================================

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
    Analyze contamination levels in each haplotype.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    group_col : str
        Column name for haplotypes
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
    - depositor_notes: Full notes field for reference
    - misidentification_confidence: Confidence level based on contamination and notes

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    group_col : str
        Column name for haplotypes
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
        Column name for haplotypes
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
    Create a heatmap showing species composition of haplotypes.

    Uses mako_r colormap with custom annotations: zeros shown as "×" on white
    background, non-zero values shown with adaptive text color for readability.
    Species names are italicized on x-axis.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    output_path : Path
        Output path for PDF
    group_col : str
        Column name for haplotypes
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

    # Ensure species values exist so groups with missing species are retained
    df_filtered[species_col] = df_filtered[species_col].fillna("Unknown")

    # Create contingency table
    ct = pd.crosstab(df_filtered[group_col], df_filtered[species_col])

    # Sort by group size
    group_order = df_filtered[group_col].value_counts().index
    # Reindex to preserve group order even if some groups had only missing species
    ct = ct.reindex(group_order, fill_value=0)

    # Mask for zero values
    zero_mask = ct.eq(0)

    # Colormap
    cmap = sns.color_palette("mako_r", as_cmap=True)
    vmin, vmax = 0, ct.to_numpy().max()

    # Create figure with adaptive sizing
    fig, ax = plt.subplots(figsize=(max(8, len(ct.columns) * 0.6), max(6, len(ct) * 0.4)))

    # Set background for masked cells (zeros) to white
    ax.set_facecolor("white")

    # Create heatmap with masked zeros
    sns.heatmap(
        ct,
        cmap=cmap,
        ax=ax,
        mask=zero_mask,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="lightgray",
        cbar_kws={'label': 'Sample Count'},
    )

    # Add custom cell annotations with adaptive text color
    def get_text_color(value):
        """Determine text color based on cell background."""
        if vmax > vmin:
            norm = (value - vmin) / (vmax - vmin)
        else:
            norm = 0
        r, g, b, _ = cmap(norm)
        luminance = 0.299 * r + 0.587 * g + 0.114 * b
        return "white" if luminance < 0.5 else "black"

    n_rows, n_cols = ct.shape
    for i in range(n_rows):
        for j in range(n_cols):
            val = ct.iat[i, j]

            if val == 0:
                # Zero cells: white background + × (always black)
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    "×",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                )
            else:
                # Non-zero cells: show value with adaptive color
                text_color = get_text_color(val)
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    f"{int(val)}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=text_color,
                )

    # Axis labels and title
    ax.set_xlabel("BOLD Reported Species")
    ax.set_ylabel("BOLDGenotyper Haplotype")
    ax.set_title(
        "BOLDGenotyper Haplotype vs. BOLD Reported Species",
        fontsize=12,
        pad=12
    )

    # X-axis tick labels: italicized species names
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    for label in ax.get_xticklabels():
        label.set_fontstyle("italic")

    # Y-axis tick labels: haplotype IDs
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    for label in ax.get_yticklabels():
        label.set_fontstyle("italic")

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
    Generate quality control report.

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
        Column name for haplotypes
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

This directory contains quality control analysis identifying potential
contamination and misidentifications in the genotype assignments.

## Files

- **mixed_species_summary.csv**: Summary of contamination in each haplotype
- **contamination_heatmap.pdf**: Visual overview of species × groups
- **depositor_flags_summary.csv**: Samples with depositor uncertainty notes
- **potential_misidentifications.csv**: Flagged samples that were likely misidentified based on group majority and depositor notes
- **purity_distribution.pdf**: Distribution of majority fractions

## Summary Statistics

- **Total haplotypes analyzed**: {n_groups}
- **Pure groups (single species)**: {n_pure} ({n_pure/n_groups*100:.1f}%)
- **Mixed groups (multiple species)**: {n_mixed} ({pct_contaminated:.1f}%)
- **Samples with depositor uncertainty flags**: {n_flagged_samples}
- **Potential misidentifications**: {n_potential_misid}

## Interpretation

### Group Purity Flags

- **CLEAN**: 100% single species - high confidence this haplotype is a homogeneous group
- **CONTAMINATED**: <100% but >70% majority - minor contamination present (some samples were likely originally misidentified or contaminated)
- **AMBIGUOUS**: <70% majority - the haplotype is valid, but default to genus-level assignment because species-level assignment is unclear

### Misidentification Confidence

- **high**: >95% identity, clear majority (>90%), no depositor uncertainty
- **medium**: 90-95% identity, moderate majority (70-90%), or has uncertainty notes
- **low**: <90% identity, weak majority (<70%), or flagged as ambiguous
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
        print("QUALITY CONTROL ALERTS")
        print("="*70)

        if n_mixed > 0:
            pct = (n_mixed / total_groups * 100) if total_groups > 0 else 0
            print(f"\n{n_mixed} of {total_groups} haplotypes contain multiple species ({pct:.1f}%)")
            print("    See: quality_control/mixed_species_summary.csv")

        if n_misid > 0:
            print(f"\n{n_misid} samples flagged as potential misidentifications")
            print("    See: quality_control/potential_misidentifications.csv")

        if n_flagged > 0:
            print(f"\n {n_flagged} samples have depositor notes indicating uncertainty")
            print("    See: quality_control/depositor_flags_summary.csv")

        print("\n TIP: Consider running family-level analysis for quality control")
        print("   Run: python scripts/compare_analyses.py --help")
        print("="*70 + "\n")
    else:
        print("\n Quality control: No significant contamination detected")
