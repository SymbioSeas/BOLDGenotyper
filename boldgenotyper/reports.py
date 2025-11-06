"""
Quality Control and Analysis Reports

This module generates summary reports and QC metrics for the genotyping pipeline,
including taxonomy conflict reports, assignment summaries, and consensus group
characterizations.
"""

import logging
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def generate_taxonomy_conflicts_report(
    annotated_csv: str,
    diagnostics_csv: str,
    output_csv: str,
) -> None:
    """
    Generate a detailed report of samples with taxonomy conflicts.

    Identifies samples where the reported BOLD species annotation differs from
    the pipeline-assigned species, and provides evidence supporting the correction.

    Parameters
    ----------
    annotated_csv : str
        Path to annotated CSV with assigned_sp and taxonomy_conflict columns
    diagnostics_csv : str
        Path to diagnostics CSV with identity scores and assignment status
    output_csv : str
        Path for output conflict report CSV
    """
    out = Path(output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load data
    df = pd.read_csv(annotated_csv)
    diag = pd.read_csv(diagnostics_csv)

    # Merge
    merged = df.merge(diag, on='processid', how='left')

    # Filter to conflicts with valid assignments
    conflicts = merged[
        (merged['taxonomy_conflict'] == True) &
        (merged['assigned_sp'].notna())
    ].copy()

    if conflicts.empty:
        logger.info("No taxonomy conflicts found with valid species assignments.")
        # Create empty output with expected columns
        empty_df = pd.DataFrame(columns=[
            'processid', 'reported_species', 'assigned_sp', 'consensus_group',
            'identity_to_assigned_group', 'n_samples_in_assigned_group',
            'species_composition_of_assigned_group', 'runner_up_group',
            'runner_up_identity', 'delta_identity', 'assignment_status',
            'evidence_summary'
        ])
        empty_df.to_csv(out, index=False)
        logger.info(f"Empty conflict report saved: {out}")
        return

    # Calculate n_samples_in_assigned_group
    group_sizes = merged.groupby('consensus_group_x')['processid'].count().to_dict()
    conflicts['n_samples_in_assigned_group'] = conflicts['consensus_group_x'].map(group_sizes)

    # Calculate species composition for each consensus group
    species_comp = merged.groupby('consensus_group_x').apply(
        lambda g: ', '.join([
            f"{sp}: {cnt} ({cnt/len(g)*100:.1f}%)"
            for sp, cnt in g['species'].value_counts().items()
        ])
    ).to_dict()
    conflicts['species_composition_of_assigned_group'] = conflicts['consensus_group_x'].map(species_comp)

    # Calculate delta identity
    conflicts['delta_identity'] = conflicts['identity'] - conflicts['runner_up_identity']

    # Generate evidence summary
    def make_evidence_summary(row):
        parts = []
        parts.append(f"Genetic identity to assigned group: {row['identity']:.3f}")
        parts.append(f"Group contains {int(row['n_samples_in_assigned_group'])} samples")

        if pd.notna(row['delta_identity']) and row['delta_identity'] > 0.01:
            parts.append(f"Identity {row['delta_identity']:.3f} higher than runner-up")

        if row['status'] == 'assigned':
            parts.append("High-confidence assignment")
        elif row['status'] == 'low_confidence':
            parts.append("Low-confidence assignment - verify manually")
        elif row['status'] == 'tie':
            parts.append("Tied assignment - verify manually")

        return '; '.join(parts)

    conflicts['evidence_summary'] = conflicts.apply(make_evidence_summary, axis=1)

    # Select and rename columns for output
    report = conflicts[[
        'processid',
        'species',  # reported species from BOLD
        'assigned_sp',
        'consensus_group_x',
        'identity',
        'n_samples_in_assigned_group',
        'species_composition_of_assigned_group',
        'runner_up_group',
        'runner_up_identity',
        'delta_identity',
        'status',
        'evidence_summary'
    ]].copy()

    report.columns = [
        'processid',
        'reported_species',
        'assigned_sp',
        'consensus_group',
        'identity_to_assigned_group',
        'n_samples_in_assigned_group',
        'species_composition_of_assigned_group',
        'runner_up_group',
        'runner_up_identity',
        'delta_identity',
        'assignment_status',
        'evidence_summary'
    ]

    # Sort by identity (highest first)
    report = report.sort_values('identity_to_assigned_group', ascending=False)

    # Save
    report.to_csv(out, index=False)
    logger.info(f"Taxonomy conflict report saved: {out} ({len(report)} conflicts)")


def generate_assignment_summary(
    annotated_csv: str,
    diagnostics_csv: str,
    output_csv: str,
) -> None:
    """
    Generate overall pipeline performance summary statistics.

    Parameters
    ----------
    annotated_csv : str
        Path to annotated CSV with species assignments
    diagnostics_csv : str
        Path to diagnostics CSV with assignment metrics
    output_csv : str
        Path for output summary CSV
    """
    out = Path(output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load data with proper CSV parsing (handle quoted fields with commas)
    df = pd.read_csv(annotated_csv, quoting=1, on_bad_lines='skip')  # quoting=1 is csv.QUOTE_MINIMAL
    diag = pd.read_csv(diagnostics_csv)

    # Merge
    merged = df.merge(diag, on='processid', how='left')

    # Calculate summary metrics
    summary = {}

    # Basic counts
    summary['total_samples'] = len(diag)
    summary['samples_with_sequence'] = len(diag[diag['status'] != 'no_sequence'])

    # Assignment status breakdown
    status_counts = diag['status'].value_counts()
    for status in ['assigned', 'tie', 'low_confidence', 'below_threshold', 'no_sequence']:
        summary[f'n_{status}'] = status_counts.get(status, 0)
        summary[f'pct_{status}'] = (status_counts.get(status, 0) / len(diag) * 100)

    # Identity score statistics (for assigned samples only)
    assigned = diag[diag['status'] == 'assigned']
    if len(assigned) > 0:
        summary['mean_identity_assigned'] = assigned['identity'].mean()
        summary['median_identity_assigned'] = assigned['identity'].median()
        summary['min_identity_assigned'] = assigned['identity'].min()
        summary['max_identity_assigned'] = assigned['identity'].max()
    else:
        summary['mean_identity_assigned'] = np.nan
        summary['median_identity_assigned'] = np.nan
        summary['min_identity_assigned'] = np.nan
        summary['max_identity_assigned'] = np.nan

    # Consensus groups and species
    summary['n_consensus_groups'] = merged['consensus_group_x'].nunique()
    summary['n_species_identified'] = merged['assigned_sp'].nunique()

    # Taxonomy conflicts
    conflicts = merged['taxonomy_conflict'].sum() if 'taxonomy_conflict' in merged.columns else 0
    assigned_with_sp = merged['assigned_sp'].notna().sum()
    summary['n_taxonomy_conflicts'] = conflicts
    summary['pct_taxonomy_conflicts'] = (conflicts / assigned_with_sp * 100) if assigned_with_sp > 0 else 0

    # Convert to DataFrame and transpose for easy reading
    summary_df = pd.DataFrame([summary])

    # Save
    summary_df.to_csv(out, index=False)
    logger.info(f"Assignment summary saved: {out}")


def generate_consensus_characterization(
    annotated_csv: str,
    diagnostics_csv: str,
    output_csv: str,
) -> None:
    """
    Generate detailed characterization of each consensus group.

    Parameters
    ----------
    annotated_csv : str
        Path to annotated CSV with species and geography
    diagnostics_csv : str
        Path to diagnostics CSV with identity scores
    output_csv : str
        Path for output characterization CSV
    """
    out = Path(output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load data
    df = pd.read_csv(annotated_csv)
    diag = pd.read_csv(diagnostics_csv)

    # Merge
    merged = df.merge(diag, on='processid', how='left')

    # Filter to assigned samples with consensus groups
    assigned = merged[
        (merged['consensus_group_x'].notna()) &
        (merged['status'] == 'assigned')
    ].copy()

    if assigned.empty:
        logger.warning("No assigned samples with consensus groups found.")
        empty_df = pd.DataFrame(columns=[
            'consensus_group', 'assigned_sp', 'n_samples', 'mean_identity',
            'min_identity', 'max_identity', 'reported_species_composition',
            'geographic_range', 'n_ocean_basins'
        ])
        empty_df.to_csv(out, index=False)
        return

    # Group by consensus group
    groups = []
    for cg, group_df in assigned.groupby('consensus_group_x'):
        char = {}
        char['consensus_group'] = cg

        # Assigned species (should be consistent within group)
        assigned_species = group_df['assigned_sp'].mode()
        char['assigned_sp'] = assigned_species.iloc[0] if len(assigned_species) > 0 else 'Unknown'

        # Assignment confidence (from assignment_notes if available)
        if 'assignment_notes' in group_df.columns:
            notes = group_df['assignment_notes'].mode()
            char['assignment_confidence'] = notes.iloc[0] if len(notes) > 0 else ''
        else:
            char['assignment_confidence'] = ''

        # Sample count
        char['n_samples'] = len(group_df)

        # Identity statistics
        char['mean_identity'] = group_df['identity'].mean()
        char['min_identity'] = group_df['identity'].min()
        char['max_identity'] = group_df['identity'].max()

        # Reported species composition (original BOLD annotations)
        species_counts = group_df['species'].value_counts()
        species_comp = ', '.join([
            f"{sp}: {cnt} ({cnt/len(group_df)*100:.1f}%)"
            for sp, cnt in species_counts.items()
        ])
        char['reported_species_composition'] = species_comp

        # Geographic range
        if 'ocean_basin' in group_df.columns:
            basins = group_df['ocean_basin'].dropna().unique()
            char['geographic_range'] = ', '.join(sorted([str(b) for b in basins]))
            char['n_ocean_basins'] = len(basins)
        else:
            char['geographic_range'] = ''
            char['n_ocean_basins'] = 0

        groups.append(char)

    # Create DataFrame
    characterization = pd.DataFrame(groups)

    # Sort by sample count (largest first)
    characterization = characterization.sort_values('n_samples', ascending=False)

    # Save
    characterization.to_csv(out, index=False)
    logger.info(f"Consensus characterization saved: {out} ({len(characterization)} groups)")


def generate_sequence_quality_metrics(
    annotated_csv: str,
    output_csv: str,
) -> None:
    """
    Generate sequence quality metrics report.

    Analyzes sequence length distributions, ambiguous base content,
    and other QC metrics.

    Parameters
    ----------
    annotated_csv : str
        Path to annotated CSV with sequence data (nuc column)
    output_csv : str
        Path for output metrics CSV
    """
    out = Path(output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load data
    df = pd.read_csv(annotated_csv)

    # Check if sequence column exists
    if 'nuc' not in df.columns:
        logger.warning("Column 'nuc' not found; skipping sequence quality metrics.")
        empty_df = pd.DataFrame(columns=[
            'metric', 'value'
        ])
        empty_df.to_csv(out, index=False)
        return

    # Filter to samples with sequences
    seqs = df[df['nuc'].notna()].copy()

    if seqs.empty:
        logger.warning("No sequences found; skipping sequence quality metrics.")
        empty_df = pd.DataFrame(columns=['metric', 'value'])
        empty_df.to_csv(out, index=False)
        return

    # Calculate metrics
    metrics = {}

    # Total sequences
    metrics['total_sequences'] = len(df)
    metrics['sequences_with_data'] = len(seqs)
    metrics['sequences_missing'] = len(df) - len(seqs)

    # Sequence length statistics
    seqs['seq_length'] = seqs['nuc'].str.len()
    metrics['mean_seq_length'] = seqs['seq_length'].mean()
    metrics['median_seq_length'] = seqs['seq_length'].median()
    metrics['min_seq_length'] = seqs['seq_length'].min()
    metrics['max_seq_length'] = seqs['seq_length'].max()
    metrics['std_seq_length'] = seqs['seq_length'].std()

    # Length thresholds
    metrics['n_seq_below_500bp'] = (seqs['seq_length'] < 500).sum()
    metrics['pct_seq_below_500bp'] = (seqs['seq_length'] < 500).sum() / len(seqs) * 100
    metrics['n_seq_below_600bp'] = (seqs['seq_length'] < 600).sum()
    metrics['pct_seq_below_600bp'] = (seqs['seq_length'] < 600).sum() / len(seqs) * 100
    metrics['n_seq_above_650bp'] = (seqs['seq_length'] > 650).sum()
    metrics['pct_seq_above_650bp'] = (seqs['seq_length'] > 650).sum() / len(seqs) * 100

    # Ambiguous base content (N's and other IUPAC codes)
    seqs['n_count'] = seqs['nuc'].str.upper().str.count('N')
    seqs['n_pct'] = (seqs['n_count'] / seqs['seq_length'] * 100)

    metrics['mean_n_count'] = seqs['n_count'].mean()
    metrics['median_n_count'] = seqs['n_count'].median()
    metrics['max_n_count'] = seqs['n_count'].max()
    metrics['n_seq_with_Ns'] = (seqs['n_count'] > 0).sum()
    metrics['pct_seq_with_Ns'] = (seqs['n_count'] > 0).sum() / len(seqs) * 100

    # Sequences with high N content (>5% or >1%)
    metrics['n_seq_over_5pct_Ns'] = (seqs['n_pct'] > 5).sum()
    metrics['pct_seq_over_5pct_Ns'] = (seqs['n_pct'] > 5).sum() / len(seqs) * 100
    metrics['n_seq_over_1pct_Ns'] = (seqs['n_pct'] > 1).sum()
    metrics['pct_seq_over_1pct_Ns'] = (seqs['n_pct'] > 1).sum() / len(seqs) * 100

    # Base composition (A, T, G, C percentages)
    for base in ['A', 'T', 'G', 'C']:
        base_counts = seqs['nuc'].str.upper().str.count(base)
        metrics[f'mean_pct_{base}'] = (base_counts / seqs['seq_length'] * 100).mean()

    # GC content
    seqs['gc_count'] = seqs['nuc'].str.upper().str.count('G') + seqs['nuc'].str.upper().str.count('C')
    seqs['gc_pct'] = (seqs['gc_count'] / seqs['seq_length'] * 100)
    metrics['mean_gc_content'] = seqs['gc_pct'].mean()
    metrics['median_gc_content'] = seqs['gc_pct'].median()
    metrics['min_gc_content'] = seqs['gc_pct'].min()
    metrics['max_gc_content'] = seqs['gc_pct'].max()

    # Convert to DataFrame
    metrics_df = pd.DataFrame([
        {'metric': k, 'value': v}
        for k, v in metrics.items()
    ])

    # Save
    metrics_df.to_csv(out, index=False)
    logger.info(f"Sequence quality metrics saved: {out}")
