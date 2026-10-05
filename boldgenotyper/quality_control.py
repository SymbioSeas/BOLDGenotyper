"""
Quality Control Module for BOLDGenotyper

This module provides quality control for COI barcode sequences, including
sequence orientation normalization, ORF validation, dynamic filtering, and a
report of label conflicts within haplotypes.

Key Features:

- COI orientation normalization (forward/reverse complement detection)
- Open reading frame (ORF) validation and quality checking
- Dynamic QC filtering with median-based thresholds
- Label conflicts: records whose label differs from their haplotype's majority
  label, with the evidence (counts, majority fraction, distance to the nearest
  haplotype carrying the record's label)
- Depositor-note keyword hits (raw notes plus the keywords matched)

Workflow:

1. Orientation Normalization: Check ORF in both orientations, correct if needed
2. ORF Validation: Check coverage and internal stop codons
3. Dynamic QC: Apply absolute and median-based length/quality filters
4. Label summary: majority labels, missing labels and label conflicts per haplotype

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
    processid_col: str = 'processid',
    return_all: bool = False,
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
    return_all : bool, optional
        Return every row with its QC columns (qc_pass, qc_fail_reason, ...)
        instead of only the passing rows (default: False)

    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        - filtered_df: DataFrame with QC columns and failed samples removed
          (all rows when return_all=True)
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

    if return_all:
        return df_with_qc, qc_stats
    return filtered_df, qc_stats


# ============================================================================
# Record labels, label conflicts and depositor-note keyword hits
# ============================================================================
#
# Record label: the record's species, or its genus when the species is missing
# (a genus-only label); missing when both are empty.
# Majority label of a haplotype: the most common species-level label among its
# members (exact and corrected); genus-only and open-nomenclature labels of the
# majority's genus are compatible; a tie for most common gives no majority.
# Label conflict: a member whose species-level label differs from the majority
# label, any species-level member of a tied haplotype, or a genus-only label of
# another genus. Missing labels are never conflicts; they are counted.

# Keywords searched for (as whole words) in depositor notes. A hit is reported
# with the note itself; it is not an interpretation of the note.
DEPOSITOR_KEYWORDS = [
    'cryptic', 'complex', 'novel', 'undescribed',
    'uncertain', 'tentative', 'probable', 'likely', 'possibly', 'unknown', 'unidentified',
    'juvenile', 'damaged', 'incomplete', 'poor condition', 'pup', 'decomposed', 'larva',
    'nymph', 'egg', 'cf.', 'aff.', 'sp.', 'spp.',
    'mixed', 'contamination', 'multiple', 'contaminated',
]
_KEYWORD_PATTERNS = [(kw, re.compile(r'(?<![A-Za-z])' + re.escape(kw) + r'(?![A-Za-z])', re.IGNORECASE))
                     for kw in DEPOSITOR_KEYWORDS]

LABEL_CONFLICT_COLUMNS = [
    'processid', 'label', 'haplotype_id', 'haplotype_taxon', 'assignment_source',
    'n_members', 'n_labeled', 'majority_label',
    'n_majority', 'majority_fraction', 'nearest_haplotype_with_label',
    'distance_to_nearest_with_label', 'depositor_keyword_hits', 'depositor_note',
]
HAPLOTYPE_LABEL_COLUMNS = [
    'haplotype_id', 'n_members', 'n_labeled', 'n_missing_label', 'n_species_labels',
    'majority_label', 'n_majority', 'majority_fraction', 'labels', 'n_label_conflicts',
]


def find_keyword_hits(note) -> List[str]:
    """Keywords found as whole words in a depositor note, in order of appearance."""
    if note is None or (isinstance(note, float) and pd.isna(note)) or not str(note).strip():
        return []
    text = str(note)
    found = []
    for kw, pattern in _KEYWORD_PATTERNS:
        m = pattern.search(text)
        if m:
            found.append((m.start(), kw))
    return [kw for _, kw in sorted(found)]


def record_labels(df: pd.DataFrame, species_col: str = "species",
                  genus_col: str = "genus") -> pd.Series:
    """Species label, or genus when the species is missing; '' when neither."""
    species = (df[species_col] if species_col in df.columns
               else pd.Series("", index=df.index)).fillna("").astype(str).str.strip()
    genus = (df[genus_col] if genus_col in df.columns
             else pd.Series("", index=df.index)).fillna("").astype(str).str.strip()
    return species.where(species != "", genus)


def _haplotype_label_stats(labels: List[str], curated_taxon: Optional[str] = None) -> Dict[str, Any]:
    """
    Majority label and counts for one haplotype's member labels.

    is_conflict compares a label with the majority label, or with curated_taxon
    when the haplotype's taxon was curated: a species-level label conflicts when
    it differs from a species-level curated taxon; any other label when its
    genus differs from the curated taxon's.
    """
    from .taxon_labels import pool_labels, label_level, genus_of
    pooled = pool_labels(labels)
    species_counts = {}
    for lab in labels:
        if lab and label_level(lab) == "species":
            species_counts[lab] = species_counts.get(lab, 0) + 1
    top = max(species_counts.values()) if species_counts else 0
    tied = [lab for lab, n in species_counts.items() if n == top] if species_counts else []
    tied_genera = {genus_of(lab) for lab in tied}

    def is_conflict(lab):
        if not lab:
            return False
        if curated_taxon:
            if label_level(lab) == "species" and label_level(curated_taxon) == "species":
                return lab != curated_taxon
            return genus_of(lab) != genus_of(curated_taxon)
        if not species_counts:
            return False
        level = label_level(lab)
        if pooled.leading_species is None:          # tie for most common
            return level == "species" or genus_of(lab) not in tied_genera
        if level == "species":
            return lab != pooled.leading_species
        return genus_of(lab) != genus_of(pooled.leading_species)

    return {
        "majority_label": pooled.leading_species,
        "n_majority": pooled.n_leading if pooled.leading_species else top,
        "majority_fraction": round(pooled.support, 4) if pooled.support is not None else np.nan,
        "n_labeled": sum(1 for lab in labels if lab),
        "n_missing_label": pooled.n_missing,
        "n_species_labels": pooled.n_distinct_species,
        "labels_text": pooled.text,
        "is_conflict": is_conflict,
    }


def add_label_columns(
    df: pd.DataFrame,
    group_col: str = "haplotype_id",
    species_col: str = "species",
    genus_col: str = "genus",
    notes_col: str = "notes",
    curated_taxa: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Add record labels, haplotype majority labels and label conflicts.

    Added columns: label, majority_label, n_majority, n_labeled,
    n_missing_label, majority_fraction, label_conflict, depositor_keyword_hits
    (';'-joined), has_depositor_keyword_hit, depositor_notes. Records without a
    haplotype get empty haplotype columns and label_conflict False.
    curated_taxa (haplotype_id -> taxon not from labels: curated or set by
    evidence): members of those haplotypes are compared with that taxon
    instead of the majority label.
    """
    df = df.copy()
    df['label'] = record_labels(df, species_col, genus_col)
    for col in ('majority_label', 'n_majority', 'n_labeled', 'n_missing_label', 'majority_fraction'):
        df[col] = np.nan
    df['majority_label'] = df['majority_label'].astype(object)
    df['label_conflict'] = False

    if group_col in df.columns:
        for hid, idx in df[df[group_col].notna()].groupby(group_col).groups.items():
            labels = df.loc[idx, 'label'].tolist()
            stats = _haplotype_label_stats(labels, (curated_taxa or {}).get(hid))
            df.loc[idx, 'majority_label'] = stats['majority_label']
            df.loc[idx, 'n_majority'] = stats['n_majority']
            df.loc[idx, 'n_labeled'] = stats['n_labeled']
            df.loc[idx, 'n_missing_label'] = stats['n_missing_label']
            df.loc[idx, 'majority_fraction'] = stats['majority_fraction']
            df.loc[idx, 'label_conflict'] = [stats['is_conflict'](lab) for lab in labels]
    df['label_conflict'] = df['label_conflict'].astype(bool)

    notes = df[notes_col] if notes_col in df.columns else pd.Series("", index=df.index)
    hits = notes.map(find_keyword_hits)
    df['depositor_keyword_hits'] = hits.map(lambda h: "; ".join(h))
    df['has_depositor_keyword_hit'] = hits.map(bool)
    df['depositor_notes'] = notes
    return df


def summarize_haplotype_labels(df: pd.DataFrame, group_col: str = "haplotype_id") -> pd.DataFrame:
    """One row per haplotype: members, labeled and missing, majority label, conflicts."""
    rows = []
    for hid, g in df[df[group_col].notna()].groupby(group_col):
        stats = _haplotype_label_stats(g['label'].tolist())
        rows.append({
            'haplotype_id': hid,
            'n_members': len(g),
            'n_labeled': stats['n_labeled'],
            'n_missing_label': stats['n_missing_label'],
            'n_species_labels': stats['n_species_labels'],
            'majority_label': stats['majority_label'],
            'n_majority': stats['n_majority'],
            'majority_fraction': stats['majority_fraction'],
            'labels': stats['labels_text'],
            'n_label_conflicts': int(g['label_conflict'].sum()),
        })
    return pd.DataFrame(rows, columns=HAPLOTYPE_LABEL_COLUMNS)


def build_label_conflicts(df: pd.DataFrame, group_col: str = "haplotype_id",
                          haplotype_sequences: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """One row per label-conflict record, with the evidence.

    haplotype_taxon / assignment_source: the haplotype's assigned_sp and how it
    was assigned (labels or curated), when df carries them.
    nearest_haplotype_with_label / distance_to_nearest_with_label: the nearest
    other haplotype with at least one member carrying the record's label, and
    the p-distance to it (dereplication's definition); empty when no other
    haplotype carries the label or sequences are not given.
    """
    from .dereplication import distances_to_references
    members = df[df[group_col].notna()]
    n_members = members.groupby(group_col).size()
    haps_by_label = members.groupby('label')[group_col].agg(lambda s: sorted(set(s))).to_dict()
    cache = {}
    rows = []
    for _, r in members[members['label_conflict']].iterrows():
        hid, lab = r[group_col], r['label']
        key = (hid, lab)
        if key not in cache:
            cache[key] = (None, np.nan)
            others = [h for h in haps_by_label.get(lab, []) if h != hid]
            if haplotype_sequences and hid in haplotype_sequences and others:
                others = [h for h in others if h in haplotype_sequences]
                if others:
                    d = distances_to_references(haplotype_sequences[hid],
                                                [haplotype_sequences[h] for h in others])
                    k = int(np.argmin(d))
                    cache[key] = (others[k], float(d[k]))
        nearest, dist = cache[key]
        rows.append({
            'processid': r.get('processid'), 'label': lab, 'haplotype_id': hid,
            'haplotype_taxon': r.get('assigned_sp', np.nan),
            'assignment_source': r.get('assignment_source', np.nan),
            'n_members': int(n_members[hid]), 'n_labeled': r['n_labeled'],
            'majority_label': r['majority_label'], 'n_majority': r['n_majority'],
            'majority_fraction': r['majority_fraction'],
            'nearest_haplotype_with_label': nearest, 'distance_to_nearest_with_label': dist,
            'depositor_keyword_hits': r.get('depositor_keyword_hits', ''),
            'depositor_note': r.get('depositor_notes', ''),
        })
    return pd.DataFrame(rows, columns=LABEL_CONFLICT_COLUMNS)


def create_label_composition_heatmap(
    df: pd.DataFrame,
    output_path: Path,
    group_col: str = "haplotype_id",
    label_col: str = "label",
    min_samples: int = 3
) -> None:
    """Heatmap of record counts: haplotypes (rows) by record label (columns).

    Zero cells are shown as "×"; missing labels as "(no species label)".
    """
    group_sizes = df[group_col].value_counts()
    valid_groups = group_sizes[group_sizes >= min_samples].index
    df_filtered = df[df[group_col].isin(valid_groups)].copy()
    if len(df_filtered) == 0:
        logger.warning("No haplotypes with enough members for the label heatmap")
        return
    df_filtered[label_col] = df_filtered[label_col].fillna("").replace("", "(no species label)")
    ct = pd.crosstab(df_filtered[group_col], df_filtered[label_col])
    ct = ct.reindex(df_filtered[group_col].value_counts().index, fill_value=0)
    zero_mask = ct.eq(0)
    cmap = sns.color_palette("mako_r", as_cmap=True)
    vmin, vmax = 0, ct.to_numpy().max()
    fig, ax = plt.subplots(figsize=(max(8, len(ct.columns) * 0.6), max(6, len(ct) * 0.4)))
    ax.set_facecolor("white")
    sns.heatmap(ct, cmap=cmap, ax=ax, mask=zero_mask, vmin=vmin, vmax=vmax,
                linewidths=0.5, linecolor="lightgray", cbar_kws={'label': 'Records'})

    def text_color(value):
        norm = (value - vmin) / (vmax - vmin) if vmax > vmin else 0
        r, g, b, _ = cmap(norm)
        return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else "black"

    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            val = ct.iat[i, j]
            ax.text(j + 0.5, i + 0.5, "×" if val == 0 else f"{int(val)}", ha="center",
                    va="center", fontsize=8, color="black" if val == 0 else text_color(val))
    ax.set_xlabel("Record label (BOLD species, or genus when species is missing)")
    ax.set_ylabel("Haplotype")
    ax.set_title("Record labels by haplotype", fontsize=12, pad=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    for label in ax.get_xticklabels():
        label.set_fontstyle("italic")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Label composition heatmap saved to {output_path}")


def create_majority_fraction_distribution(summary: pd.DataFrame, output_path: Path) -> None:
    """Histogram of each haplotype's majority_fraction (no thresholds drawn)."""
    values = summary['majority_fraction'].dropna().to_numpy() if len(summary) else np.array([])
    if len(values) == 0:
        logger.warning("No majority fractions to plot")
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(values, bins=np.linspace(0, 1, 21), edgecolor='black', linewidth=1.0, color='#4C72B0')
    ax.set_xlabel('Majority fraction (majority label / labeled members)', fontsize=12)
    ax.set_ylabel('Haplotypes', fontsize=12)
    ax.set_title('Majority fraction per haplotype', fontsize=14, pad=20)
    ax.grid(axis='y', alpha=0.3)
    ax.text(0.98, 0.98, f"n = {len(values)}\nmedian = {np.median(values):.2f}",
            transform=ax.transAxes, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Majority fraction distribution saved to {output_path}")


def generate_quality_control_report(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    group_col: str = "haplotype_id",
    species_col: str = "species",
    haplotype_sequences: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Write the label-conflict report.

    Files: label_conflicts.csv, haplotype_label_summary.csv,
    depositor_keyword_hits.csv, label_composition_heatmap.pdf,
    majority_fraction_distribution.pdf, README.md.

    df must carry the columns from add_label_columns. haplotype_sequences
    (haplotype_id -> aligned core sequence) enables the distance evidence.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if 'label_conflict' not in df.columns:
        df = add_label_columns(df, group_col=group_col, species_col=species_col)
    results: Dict[str, Any] = {}

    summary = summarize_haplotype_labels(df, group_col)
    summary.to_csv(output_dir / "haplotype_label_summary.csv", index=False)
    results['haplotype_label_summary'] = output_dir / "haplotype_label_summary.csv"

    conflicts = build_label_conflicts(df, group_col, haplotype_sequences)
    conflicts.to_csv(output_dir / "label_conflicts.csv", index=False)
    results['label_conflicts'] = output_dir / "label_conflicts.csv"

    hits = df[df['has_depositor_keyword_hit']]
    hits_out = pd.DataFrame({
        'processid': hits.get('processid'), 'haplotype_id': hits.get(group_col),
        'label': hits['label'], 'depositor_note': hits['depositor_notes'],
        'keyword_hits': hits['depositor_keyword_hits'],
    })
    hits_out.to_csv(output_dir / "depositor_keyword_hits.csv", index=False)
    results['depositor_keyword_hits'] = output_dir / "depositor_keyword_hits.csv"

    try:
        create_label_composition_heatmap(df, output_dir / "label_composition_heatmap.pdf", group_col)
        results['label_composition_heatmap'] = output_dir / "label_composition_heatmap.pdf"
    except Exception as e:
        logger.warning(f"Could not create label composition heatmap: {e}",
                       extra={"step_failure": "6 post-assignment QC"})
    try:
        create_majority_fraction_distribution(summary, output_dir / "majority_fraction_distribution.pdf")
        results['majority_fraction_distribution'] = output_dir / "majority_fraction_distribution.pdf"
    except Exception as e:
        logger.warning(f"Could not create majority fraction distribution: {e}",
                       extra={"step_failure": "6 post-assignment QC"})

    n_haps = len(summary)
    n_with_conflict = int((summary['n_label_conflicts'] > 0).sum()) if n_haps else 0
    n_conflicts = len(conflicts)
    n_missing = int(summary['n_missing_label'].sum()) if n_haps else 0
    results['summary'] = {
        'total_groups': n_haps,
        'groups_with_label_conflict': n_with_conflict,
        'n_label_conflicts': n_conflicts,
        'n_missing_label': n_missing,
        'n_keyword_hits': len(hits_out),
    }

    readme = f"""# Record labels within haplotypes: {organism}

These tables compare each haplotype's member records by their BOLD labels. They
report counts; they do not decide which label is correct.

## Definitions

- **label**: the record's BOLD species, or its genus when the species is missing
  (a genus-only label); empty when both are missing.
- **majority label**: the most common species-level label among a haplotype's
  members (exact and corrected). Genus-only and open-nomenclature labels
  (`sp.`, `cf.`, `aff.`) of the majority's genus are compatible with it. When
  two or more labels tie for most common, there is no majority label.
- **majority_fraction**: members with the majority label / members with a
  species-level label (plus genus-only labels of other genera). Records without
  a label are not in the denominator; they are counted in `n_missing_label`.
- **label conflict**: a member whose species-level label differs from the
  majority label, any species-level member of a haplotype without a majority
  label, or a genus-only label of another genus.

## Files

- **label_conflicts.csv**: one row per label-conflict record:
  `processid, label, haplotype_id, n_members, n_labeled, majority_label,
  n_majority, majority_fraction, nearest_haplotype_with_label,
  distance_to_nearest_with_label, depositor_keyword_hits, depositor_note`.
  `nearest_haplotype_with_label` is the nearest other haplotype with at least
  one member carrying the record's label, and `distance_to_nearest_with_label`
  the uncorrected p-distance to it (empty when no other haplotype carries it).
- **haplotype_label_summary.csv**: one row per haplotype:
  `haplotype_id, n_members, n_labeled, n_missing_label, n_species_labels,
  majority_label, n_majority, majority_fraction, labels, n_label_conflicts`
  (`labels` lists every label with its count).
- **depositor_keyword_hits.csv**: records whose depositor note contains one of a
  fixed list of keywords (whole words), with the note and the keywords matched.
  A keyword hit is a pointer to read the note, not an interpretation of it.
- **label_composition_heatmap.pdf**: record counts by haplotype and label.
- **majority_fraction_distribution.pdf**: histogram of `majority_fraction`.

## Counts

- Haplotypes: {n_haps}
- Haplotypes with at least one label conflict: {n_with_conflict}
- Label-conflict records: {n_conflicts}
- Members without a label: {n_missing}
- Records with a depositor-note keyword hit: {len(hits_out)}
"""
    (output_dir / "README.md").write_text(readme)
    results['readme'] = output_dir / "README.md"
    logger.info(f"Label report written to {output_dir}")
    return results


def print_quality_alert(qc_summary: Dict[str, Any]) -> None:
    """Print label-conflict counts to the console."""
    s = qc_summary.get('summary', {})
    n_conflicts = s.get('n_label_conflicts', 0)
    if n_conflicts or s.get('n_keyword_hits', 0):
        print("\n" + "=" * 70)
        print("RECORD LABELS")
        print("=" * 70)
        print(f"{s.get('groups_with_label_conflict', 0)} of {s.get('total_groups', 0)} haplotypes "
              f"have label conflicts ({n_conflicts} records); "
              f"{s.get('n_missing_label', 0)} members have no label")
        print("    See: quality_control/label_conflicts.csv and haplotype_label_summary.csv")
        if s.get('n_keyword_hits', 0):
            print(f"{s['n_keyword_hits']} records have depositor-note keyword hits")
            print("    See: quality_control/depositor_keyword_hits.csv")
        print("=" * 70 + "\n")
    else:
        print("\nRecord labels: no label conflicts within haplotypes")
