"""
Quality Control and Analysis Reports

This module generates summary reports and QC metrics for the genotyping pipeline,
including taxonomy conflict reports, assignment summaries, and consensus group
characterizations.
"""

import logging
import json
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


# ============================================================================
# HTML Report Generation
# ============================================================================

# HTML Report CSS
HTML_REPORT_CSS = """
<style>
    :root {
        --primary-color: #2c3e50;
        --secondary-color: #3498db;
        --success-color: #27ae60;
        --warning-color: #f39c12;
        --danger-color: #e74c3c;
        --light-bg: #ecf0f1;
        --border-color: #bdc3c7;
    }

    * {
        margin: 0;
        padding: 0;
        box-sizing: border-box;
    }

    body {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        line-height: 1.6;
        color: #333;
        background-color: #f5f5f5;
        margin: 0;
        padding: 0;
    }

    .container {
        display: flex;
        min-height: 100vh;
        background-color: white;
    }
    
    .geo-map {
        position: relative;
        width: 100%;
        max-width: 100%;
        height: 400px;
        border: 1px solid var(--border-color);
        border-radius: 6px;
        background-color: #f8f9fa;
        margin-top: 15px;
        margin-bottom: 25px;
        overflow: hidden;
    }

    .geo-map-svg {
        width: 100%;
        height: 100%;
    }

    .geo-point {
        fill: var(--secondary-color);
        stroke: white;
        stroke-width: 1;
        opacity: 0.8;
        cursor: pointer;
    }

    .geo-point:hover {
        opacity: 1.0;
    }

    .geo-tooltip {
        position: absolute;
        display: none;
        pointer-events: none;
        background-color: rgba(0, 0, 0, 0.75);
        color: white;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.8em;
        z-index: 10;
        white-space: nowrap;
    }
    
    .header {
        background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
        color: white;
        padding: 30px 40px;
        text-align: center;
    }

    .header h1 {
        font-size: 2em;
        margin-bottom: 8px;
        font-style: italic;
    }

    .header .subtitle {
        font-size: 1em;
        opacity: 0.9;
        margin-bottom: 5px;
    }

    .header .timestamp {
        font-size: 0.85em;
        opacity: 0.8;
    }

    .quick-stats {
        background-color: var(--light-bg);
        padding: 20px 40px;
        display: flex;
        justify-content: space-around;
        flex-wrap: wrap;
        border-bottom: 3px solid var(--secondary-color);
    }

    .quick-stat {
        text-align: center;
        padding: 10px 20px;
    }

    .quick-stat .number {
        font-size: 2em;
        font-weight: bold;
        color: var(--secondary-color);
    }

    .quick-stat .label {
        font-size: 0.9em;
        color: #666;
        margin-top: 5px;
    }

    .content {
        padding: 40px;
        flex: 1;
    }

    .section {
        margin-bottom: 50px;
    }

    .section h2 {
        color: var(--primary-color);
        border-bottom: 2px solid var(--secondary-color);
        padding-bottom: 10px;
        margin-bottom: 20px;
        font-size: 1.8em;
    }

    .section h3 {
        color: var(--primary-color);
        margin-top: 30px;
        margin-bottom: 15px;
        font-size: 1.3em;
    }

    /* Left Sidebar Navigation */
    .sidebar {
        width: 240px;
        background: linear-gradient(180deg, var(--primary-color), #1a252f);
        color: white;
        position: fixed;
        height: 100vh;
        overflow-y: auto;
        padding: 20px 0;
        box-shadow: 2px 0 10px rgba(0,0,0,0.1);
        z-index: 1000;
    }

    .sidebar-header {
        padding: 0 20px 20px 20px;
        border-bottom: 1px solid rgba(255,255,255,0.2);
        margin-bottom: 20px;
    }

    .sidebar-title {
        font-size: 1.1em;
        font-weight: 600;
        margin-bottom: 5px;
        font-style: italic;
    }

    .sidebar-version {
        font-size: 0.75em;
        opacity: 0.7;
    }

    .sidebar-nav {
        list-style: none;
        padding: 0;
        margin: 0;
    }

    .sidebar-nav-item {
        margin: 0;
    }

    .sidebar-nav-link {
        display: block;
        padding: 12px 20px;
        text-decoration: none;
        color: rgba(255,255,255,0.9);
        border-left: 3px solid transparent;
        transition: all 0.2s ease;
        font-size: 0.95em;
    }

    .sidebar-nav-link:hover {
        background-color: rgba(255,255,255,0.1);
        border-left-color: var(--secondary-color);
        color: white;
    }

    .sidebar-nav-link.active {
        background-color: rgba(255,255,255,0.15);
        border-left-color: var(--secondary-color);
        color: white;
        font-weight: 600;
    }

    /* Main Content Area */
    .main-wrapper {
        flex: 1;
        margin-left: 240px;
        display: flex;
        flex-direction: column;
    }
    
    .metric-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
        gap: 20px;
        margin-bottom: 30px;
    }

    .metric-card {
        background-color: white;
        border: 1px solid var(--border-color);
        border-radius: 6px;
        padding: 20px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }

    .metric-card .title {
        font-size: 0.85em;
        color: #666;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-bottom: 10px;
    }

    .metric-card .value {
        font-size: 2em;
        font-weight: bold;
        color: var(--secondary-color);
    }

    .metric-card .description {
        font-size: 0.9em;
        color: #888;
        margin-top: 5px;
    }

    .table-container {
        overflow-x: auto;
        margin: 20px 0;
    }

    table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.9em;
    }

    thead {
        background-color: var(--primary-color);
        color: white;
    }

    th {
        padding: 12px;
        text-align: left;
        font-weight: 600;
    }

    td {
        padding: 10px 12px;
        border-bottom: 1px solid #eee;
    }

    tbody tr:hover {
        background-color: #f8f9fa;
    }

    tbody tr:nth-child(even) {
        background-color: #fafafa;
    }

    .alert {
        padding: 15px;
        border-radius: 6px;
        margin: 20px 0;
        border-left: 4px solid;
    }

    .alert-info {
        background-color: #d1ecf1;
        border-color: #17a2b8;
        color: #0c5460;
    }

    .alert-warning {
        background-color: #fff3cd;
        border-color: #ffc107;
        color: #856404;
    }

    .alert-success {
        background-color: #d4edda;
        border-color: #28a745;
        color: #155724;
    }

    .badge {
        display: inline-block;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.85em;
        font-weight: 600;
    }

    .badge-success {
        background-color: var(--success-color);
        color: white;
    }

    .badge-warning {
        background-color: var(--warning-color);
        color: white;
    }

    .badge-danger {
        background-color: var(--danger-color);
        color: white;
    }

    .badge-info {
        background-color: var(--secondary-color);
        color: white;
    }

    .footer {
        background-color: var(--light-bg);
        padding: 30px 40px;
        text-align: center;
        border-top: 1px solid var(--border-color);
        font-size: 0.9em;
        color: #666;
    }

    .footer a {
        color: var(--secondary-color);
        text-decoration: none;
    }

    .footer a:hover {
        text-decoration: underline;
    }

    @media print {
        body {
            background-color: white;
            padding: 0;
        }
        .container {
            box-shadow: none;
        }
        .section {
            page-break-inside: avoid;
        }
    }

    @media (max-width: 768px) {
        .sidebar {
            width: 200px;
        }
        .main-wrapper {
            margin-left: 200px;
        }
        .header h1 {
            font-size: 1.5em;
        }
        .quick-stats {
            flex-direction: column;
            padding: 15px 20px;
        }
        .metric-grid {
            grid-template-columns: 1fr;
        }
        .content {
            padding: 20px;
        }
    }

    @media (max-width: 576px) {
        .sidebar {
            width: 60px;
            padding: 10px 0;
        }
        .sidebar-header {
            padding: 10px;
        }
        .sidebar-title,
        .sidebar-version {
            display: none;
        }
        .sidebar-nav-link {
            padding: 12px 10px;
            font-size: 0.8em;
            text-align: center;
        }
        .main-wrapper {
            margin-left: 60px;
        }
        .header {
            padding: 20px;
        }
        .content {
            padding: 15px;
        }
    }

    /* Subtab Navigation */
    .subtabs {
        display: flex;
        flex-wrap: wrap;
        gap: 5px;
        margin-bottom: 20px;
        border-bottom: 2px solid var(--border-color);
        padding-bottom: 0;
    }

    .subtab-button {
        padding: 10px 16px;
        background-color: #f5f5f5;
        border: 1px solid var(--border-color);
        border-bottom: none;
        border-radius: 6px 6px 0 0;
        cursor: pointer;
        font-size: 0.9em;
        color: #666;
        transition: all 0.2s ease;
    }

    .subtab-button:hover {
        background-color: #e8e8e8;
        color: var(--primary-color);
    }

    .subtab-button.active {
        background-color: white;
        color: var(--primary-color);
        font-weight: 600;
        border-bottom: 2px solid white;
        margin-bottom: -2px;
    }

    .subtab-content {
        display: none;
        animation: fadeIn 0.3s ease-in;
    }

    .subtab-content.active {
        display: block;
    }

    @keyframes fadeIn {
        from { opacity: 0; }
        to { opacity: 1; }
    }

    /* Visualization Styles */
    .viz-container {
        margin: 30px 0;
        padding: 20px;
        background-color: #f9f9f9;
        border-radius: 8px;
        border: 1px solid var(--border-color);
    }

    .viz-container h3 {
        margin-top: 0;
        color: var(--primary-color);
        font-size: 1.3em;
        margin-bottom: 10px;
    }

    .viz-description {
        color: #666;
        font-size: 0.95em;
        margin-bottom: 15px;
        font-style: italic;
    }

    .viz-image {
        width: 100%;
        max-width: 100%;
        height: auto;
        border: 1px solid #ddd;
        border-radius: 4px;
        background-color: white;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        cursor: zoom-in;
    }

    .viz-image:hover {
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }

    /* Interactive Plot Controls */
    .plot-mode-toggle {
        display: flex;
        gap: 10px;
        margin-bottom: 20px;
        justify-content: center;
    }

    .toggle-btn {
        padding: 10px 24px;
        border: 2px solid var(--primary-color);
        background: white;
        color: var(--primary-color);
        font-weight: 600;
        border-radius: 6px;
        cursor: pointer;
        transition: all 0.3s ease;
        font-size: 14px;
    }

    .toggle-btn:hover {
        background: #f0f8ff;
    }

    .toggle-btn.active {
        background: var(--primary-color);
        color: white;
    }

    .plot-view {
        display: none;
    }

    .plot-view.active {
        display: block;
    }

    .plot-controls {
        background: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
    }

    .controls-section h4 {
        margin-top: 0;
        margin-bottom: 16px;
        color: var(--text-color);
        font-size: 16px;
        font-weight: 600;
    }

    .control-group {
        margin-bottom: 20px;
    }

    .control-group label {
        display: block;
        font-weight: 600;
        margin-bottom: 8px;
        color: var(--text-color);
    }

    .threshold-slider {
        width: 100%;
        height: 6px;
        border-radius: 3px;
        background: #ddd;
        outline: none;
        -webkit-appearance: none;
    }

    .threshold-slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: var(--primary-color);
        cursor: pointer;
        transition: all 0.2s ease;
    }

    .threshold-slider::-webkit-slider-thumb:hover {
        transform: scale(1.2);
        background: #0056b3;
    }

    .threshold-slider::-moz-range-thumb {
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: var(--primary-color);
        cursor: pointer;
        border: none;
    }

    .checkbox-controls {
        display: flex;
        gap: 8px;
        margin-bottom: 12px;
    }

    .genotype-checkboxes {
        max-height: 300px;
        overflow-y: auto;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        padding: 12px;
        background: white;
    }

    .genotype-checkbox-item {
        display: flex;
        align-items: center;
        padding: 6px 0;
        border-bottom: 1px solid #f0f0f0;
    }

    .genotype-checkbox-item:last-child {
        border-bottom: none;
    }

    .genotype-checkbox-item input[type="checkbox"] {
        margin-right: 10px;
        width: 18px;
        height: 18px;
        cursor: pointer;
    }

    .genotype-checkbox-item label {
        margin: 0;
        font-weight: normal;
        cursor: pointer;
        flex: 1;
        display: flex;
        justify-content: space-between;
        align-items: center;
    }

    .genotype-sample-count {
        color: #666;
        font-size: 0.9em;
        margin-left: 8px;
    }

    .btn-small {
        padding: 6px 12px;
        font-size: 13px;
        border: 1px solid #dee2e6;
        background: white;
        color: var(--text-color);
        border-radius: 4px;
        cursor: pointer;
        transition: all 0.2s ease;
    }

    .btn-small:hover {
        background: #e9ecef;
        border-color: #adb5bd;
    }

    .btn-primary {
        padding: 10px 20px;
        background: var(--primary-color);
        color: white;
        border: none;
        border-radius: 6px;
        font-weight: 600;
        cursor: pointer;
        transition: all 0.2s ease;
        font-size: 14px;
        width: 100%;
    }

    .btn-primary:hover {
        background: #0056b3;
        transform: translateY(-1px);
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }

    .plotly-chart {
        min-height: 500px;
        border: 1px solid #ddd;
        border-radius: 4px;
        background: white;
    }

    /* Summary Statistics Panel */
    .stats-panel {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }

    .stats-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 20px;
    }

    .stat-item {
        text-align: center;
        background: rgba(255, 255, 255, 0.15);
        padding: 15px;
        border-radius: 6px;
        backdrop-filter: blur(10px);
    }

    .stat-label {
        font-size: 12px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        color: rgba(255, 255, 255, 0.9);
        margin-bottom: 8px;
    }

    .stat-value {
        font-size: 28px;
        font-weight: 700;
        color: white;
        line-height: 1;
    }

    @media (max-width: 768px) {
        .stats-grid {
            grid-template-columns: 1fr;
            gap: 12px;
        }
    }

    /* Download Options */
    .download-options {
        display: flex;
        gap: 8px;
        margin-bottom: 12px;
        flex-wrap: wrap;
    }

    .btn-download {
        flex: 1;
        min-width: 100px;
        padding: 10px 16px;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 6px;
        font-weight: 600;
        cursor: pointer;
        transition: all 0.3s ease;
        font-size: 13px;
    }

    .btn-download:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(102, 126, 234, 0.4);
    }

    .resolution-control {
        margin-top: 8px;
    }

    .small-label {
        font-size: 13px;
        font-weight: 600;
        display: block;
        margin-bottom: 6px;
        color: var(--text-color);
    }

    .resolution-select {
        width: 100%;
        padding: 8px 12px;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        background: white;
        font-size: 13px;
        cursor: pointer;
        transition: border-color 0.2s ease;
    }

    .resolution-select:hover {
        border-color: var(--primary-color);
    }

    .resolution-select:focus {
        outline: none;
        border-color: var(--primary-color);
        box-shadow: 0 0 0 3px rgba(0, 123, 255, 0.1);
    }

    @media print {
        .viz-container {
            page-break-inside: avoid;
        }
        .subtab-button:not(.active) {
            display: none;
        }
        .subtab-content {
            display: block !important;
        }
        .plot-mode-toggle,
        .plot-controls {
            display: none !important;
        }
        .plot-view {
            display: none !important;
        }
        .plot-view.active {
            display: block !important;
        }
    }
</style>
"""

# HTML Report Template
HTML_REPORT_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="generator" content="BOLDGenotyper">
    <title>{{ organism }} - BOLDGenotyper Report</title>
    {{ css | safe }}
    <!-- Plotly.js for interactive charts -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
</head>
<body>
    <div class="container">
        <!-- Left Sidebar Navigation -->
        <aside class="sidebar">
            <div class="sidebar-header">
                <div class="sidebar-title">{{ organism }}</div>
                <div class="sidebar-version">BOLDGenotyper v{{ version }}</div>
            </div>
            <nav>
                <ul class="sidebar-nav">
                    <li class="sidebar-nav-item">
                        <a href="#section-executive-summary" class="sidebar-nav-link">Summary</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-parameters" class="sidebar-nav-link">Parameters</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-assignment" class="sidebar-nav-link">Assignment</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-taxonomy" class="sidebar-nav-link">Taxonomy</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-geographic" class="sidebar-nav-link">Geography</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-visualizations" class="sidebar-nav-link">Visualizations</a>
                    </li>
                </ul>
            </nav>
        </aside>

        <!-- Main Content Area -->
        <div class="main-wrapper">
            <!-- Header -->
            <div class="header">
                <h1>{{ organism }}</h1>
                <div class="subtitle">BOLDGenotyper Analysis Report</div>
                <div class="timestamp">Generated: {{ timestamp }}</div>
            </div>

            <!-- Quick Stats -->
            {% if quick_stats %}
            <div class="quick-stats">
                {% for stat in quick_stats %}
                <div class="quick-stat">
                    <div class="number">{{ stat.value }}</div>
                    <div class="label">{{ stat.label }}</div>
                </div>
                {% endfor %}
            </div>
            {% endif %}

            <!-- Main Content -->
            <div class="content">
                {{ content | safe }}
            </div>

            <!-- Footer -->
            <div class="footer">
                <p>
                    <strong>BOLDGenotyper</strong> - Automated genotyping pipeline for BOLD barcode data<br>
                    <a href="https://github.com/anthropics/boldgenotyper" target="_blank">GitHub Repository</a>
                </p>
                <p style="margin-top: 15px; font-size: 0.85em;">
                    Report generated on {{ timestamp }}<br>
                    For questions or issues, please visit the GitHub repository.
                </p>
            </div>
        </div>
    </div>

    <script>
        // Switch between subtabs
        function switchSubtab(tabId) {
            // Find the parent section to scope the tab switching
            const tabElement = document.getElementById(tabId);
            if (!tabElement) return;

            const parentSection = tabElement.closest('.section');
            if (!parentSection) return;

            // Hide all subtab contents in this section
            const allContents = parentSection.querySelectorAll('.subtab-content');
            allContents.forEach(content => {
                content.classList.remove('active');
            });

            // Deactivate all subtab buttons in this section
            const allButtons = parentSection.querySelectorAll('.subtab-button');
            allButtons.forEach(button => {
                button.classList.remove('active');
            });

            // Show selected tab content
            tabElement.classList.add('active');

            // Activate corresponding button
            const buttonIndex = Array.from(allContents).indexOf(tabElement);
            if (buttonIndex >= 0 && allButtons[buttonIndex]) {
                allButtons[buttonIndex].classList.add('active');
            }
        }

        // Highlight active section in sidebar based on scroll position
        document.addEventListener('DOMContentLoaded', function() {
            const sections = document.querySelectorAll('.section[id^="section-"]');
            const navLinks = document.querySelectorAll('.sidebar-nav-link');

            function updateActiveLink() {
                let currentSection = '';
                sections.forEach(section => {
                    const sectionTop = section.offsetTop - 100;
                    if (window.pageYOffset >= sectionTop) {
                        currentSection = section.getAttribute('id');
                    }
                });

                navLinks.forEach(link => {
                    link.classList.remove('active');
                    if (link.getAttribute('href') === '#' + currentSection) {
                        link.classList.add('active');
                    }
                });
            }

            window.addEventListener('scroll', updateActiveLink);
            updateActiveLink(); // Initial call

            // Smooth scrolling for anchor links
            navLinks.forEach(link => {
                link.addEventListener('click', function(e) {
                    e.preventDefault();
                    const targetId = this.getAttribute('href').substring(1);
                    const targetElement = document.getElementById(targetId);
                    if (targetElement) {
                        window.scrollTo({
                            top: targetElement.offsetTop - 20,
                            behavior: 'smooth'
                        });
                    }
                });
            });

            // Initialize interactive plots when DOM is loaded
            initializeInteractivePlots();
        });

        // Global state for interactive plots
        const plotState = {};

        // Initialize all interactive plots
        function initializeInteractivePlots() {
            // Find all embedded plot data
            const plotDataElements = document.querySelectorAll('script[type="application/json"][id^="plot-data-"]');

            plotDataElements.forEach(dataElement => {
                const idx = parseInt(dataElement.id.replace('plot-data-', ''));
                const plotData = JSON.parse(dataElement.textContent);

                // Initialize plot state
                plotState[idx] = {
                    data: plotData,
                    threshold: 0,
                    selectedGenotypes: new Set(plotData.genotypes || []),
                    selectedBasins: new Set(plotData.ocean_basins || [])
                };

                // Populate genotype checkboxes
                populateGenotypeCheckboxes(idx);

                // Populate ocean basin checkboxes if available
                if (plotData.ocean_basins && plotData.ocean_basins.length > 0) {
                    populateBasinCheckboxes(idx);
                }
            });
        }

        // Toggle between static and interactive plot views
        function showStaticPlot(plotId, idx) {
            // Update button states
            const buttons = document.querySelectorAll(`button[onclick*="'${plotId}'"]`);
            buttons.forEach(btn => {
                if (btn.textContent.includes('Static')) {
                    btn.classList.add('active');
                } else {
                    btn.classList.remove('active');
                }
            });

            // Show static, hide interactive
            document.getElementById(`${plotId}-static`).style.display = 'block';
            document.getElementById(`${plotId}-interactive`).style.display = 'none';
        }

        function showInteractivePlot(plotId, idx) {
            // Update button states
            const buttons = document.querySelectorAll(`button[onclick*="'${plotId}'"]`);
            buttons.forEach(btn => {
                if (btn.textContent.includes('Interactive')) {
                    btn.classList.add('active');
                } else {
                    btn.classList.remove('active');
                }
            });

            // Show interactive, hide static
            document.getElementById(`${plotId}-static`).style.display = 'none';
            document.getElementById(`${plotId}-interactive`).style.display = 'block';

            // Render the plot if not already rendered
            if (!plotState[idx].rendered) {
                updatePlot(idx);
                plotState[idx].rendered = true;
            }
        }

        // Populate genotype checkboxes
        function populateGenotypeCheckboxes(idx) {
            const container = document.getElementById(`genotype-checkboxes-${idx}`);
            if (!container) return;

            const state = plotState[idx];
            const plotData = state.data;

            container.innerHTML = '';

            const genotypes = plotData.genotypes || [];
            const sampleCounts = plotData.sample_counts || {};

            genotypes.forEach(genotype => {
                const div = document.createElement('div');
                div.className = 'genotype-checkbox-item';

                const checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.id = `genotype-${idx}-${genotype}`;
                checkbox.value = genotype;
                checkbox.checked = state.selectedGenotypes.has(genotype);
                checkbox.onchange = () => {
                    toggleGenotype(idx, genotype);
                    updatePlot(idx);
                };

                const label = document.createElement('label');
                label.htmlFor = checkbox.id;
                label.innerHTML = `
                    <span>${genotype}</span>
                    <span class="genotype-sample-count">n=${sampleCounts[genotype] || 0}</span>
                `;

                div.appendChild(checkbox);
                div.appendChild(label);
                container.appendChild(div);
            });
        }

        // Toggle genotype selection
        function toggleGenotype(idx, genotype) {
            const state = plotState[idx];
            if (state.selectedGenotypes.has(genotype)) {
                state.selectedGenotypes.delete(genotype);
            } else {
                state.selectedGenotypes.add(genotype);
            }
        }

        // Select all genotypes
        function selectAllGenotypes(idx) {
            const state = plotState[idx];
            const genotypes = state.data.genotypes || [];

            genotypes.forEach(g => state.selectedGenotypes.add(g));

            // Update checkbox states
            genotypes.forEach(g => {
                const checkbox = document.getElementById(`genotype-${idx}-${g}`);
                if (checkbox) checkbox.checked = true;
            });

            updatePlot(idx);
        }

        // Deselect all genotypes
        function deselectAllGenotypes(idx) {
            const state = plotState[idx];

            state.selectedGenotypes.clear();

            // Update checkbox states
            const checkboxes = document.querySelectorAll(`input[id^="genotype-${idx}-"]`);
            checkboxes.forEach(cb => cb.checked = false);

            updatePlot(idx);
        }

        // Populate ocean basin checkboxes
        function populateBasinCheckboxes(idx) {
            const container = document.getElementById(`basin-checkboxes-${idx}`);
            if (!container) return;

            const state = plotState[idx];
            const plotData = state.data;
            const basins = plotData.ocean_basins || [];

            container.innerHTML = '';

            basins.forEach(basin => {
                const div = document.createElement('div');
                div.className = 'genotype-checkbox-item';

                const checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.id = `basin-${idx}-${basin}`;
                checkbox.value = basin;
                checkbox.checked = state.selectedBasins.has(basin);
                checkbox.onchange = () => {
                    toggleBasin(idx, basin);
                    updatePlot(idx);
                };

                const label = document.createElement('label');
                label.htmlFor = checkbox.id;
                label.textContent = basin;

                div.appendChild(checkbox);
                div.appendChild(label);
                container.appendChild(div);
            });
        }

        // Toggle ocean basin selection
        function toggleBasin(idx, basin) {
            const state = plotState[idx];
            if (state.selectedBasins.has(basin)) {
                state.selectedBasins.delete(basin);
            } else {
                state.selectedBasins.add(basin);
            }
        }

        // Select all ocean basins
        function selectAllBasins(idx) {
            const state = plotState[idx];
            const basins = state.data.ocean_basins || [];

            basins.forEach(b => state.selectedBasins.add(b));

            // Update checkbox states
            basins.forEach(b => {
                const checkbox = document.getElementById(`basin-${idx}-${b}`);
                if (checkbox) checkbox.checked = true;
            });

            updatePlot(idx);
        }

        // Deselect all ocean basins
        function deselectAllBasins(idx) {
            const state = plotState[idx];

            state.selectedBasins.clear();

            // Update checkbox states
            const checkboxes = document.querySelectorAll(`input[id^="basin-${idx}-"]`);
            checkboxes.forEach(cb => cb.checked = false);

            updatePlot(idx);
        }

        // Update threshold value
        function updateThreshold(idx) {
            const slider = document.getElementById(`threshold-${idx}`);
            const valueDisplay = document.getElementById(`threshold-value-${idx}`);

            if (slider && valueDisplay) {
                plotState[idx].threshold = parseInt(slider.value);
                valueDisplay.textContent = slider.value;
            }
        }

        // Update summary statistics
        function updateStatistics(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // Calculate filtered genotypes
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            // Calculate total samples in filtered data
            let totalSamples = 0;
            filteredGenotypes.forEach(g => {
                totalSamples += plotData.sample_counts[g] || 0;
            });

            // Calculate percentage of total data
            let allSamples = 0;
            (plotData.genotypes || []).forEach(g => {
                allSamples += plotData.sample_counts[g] || 0;
            });
            const percentage = allSamples > 0 ? (totalSamples / allSamples * 100).toFixed(1) : 0;

            // Update DOM elements
            const genotypesStat = document.getElementById(`stat-genotypes-${idx}`);
            const samplesStat = document.getElementById(`stat-samples-${idx}`);
            const percentageStat = document.getElementById(`stat-percentage-${idx}`);

            if (genotypesStat) {
                genotypesStat.textContent = `${filteredGenotypes.length} / ${plotData.genotypes.length}`;
            }
            if (samplesStat) {
                samplesStat.textContent = totalSamples.toLocaleString();
            }
            if (percentageStat) {
                percentageStat.textContent = `${percentage}%`;
            }
        }

        // Update/render the Plotly chart
        function updatePlot(idx) {
            const state = plotState[idx];
            const plotData = state.data;
            const plotType = plotData.plot_type;

            // Update statistics
            updateStatistics(idx);

            if (plotType === 'faceted_abundance') {
                renderFacetedPlot(idx);
            } else if (plotType === 'distribution_map') {
                renderDistributionMap(idx);
            } else {
                renderStackedBarPlot(idx);
            }
        }

        // Render stacked bar plot (for relative and total abundance)
        function renderStackedBarPlot(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // Filter genotypes based on selection and threshold
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            if (filteredGenotypes.length === 0) {
                // Show message if no genotypes selected
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No genotypes selected or meet the minimum sample threshold.</div>';
                }
                return;
            }

            // Prepare Plotly traces (one per genotype)
            const traces = [];
            const basins = plotData.basins || [];
            const basinLabels = plotData.basin_labels || {};
            const colors = plotData.colors || [];
            const counts = plotData.counts || {};

            filteredGenotypes.forEach((genotype, i) => {
                const genotypeCounts = counts[genotype] || {};
                const y = basins.map(basin => genotypeCounts[basin] || 0);

                traces.push({
                    name: genotype,
                    x: basins.map(b => basinLabels[b] || b),
                    y: y,
                    type: 'bar',
                    marker: {
                        color: colors[genotype] || '#888'
                    }
                });
            });

            // Layout configuration
            const isRelative = plotData.plot_type === 'relative_abundance';
            const layout = {
                barmode: 'stack',
                title: '',
                xaxis: {
                    title: 'Ocean Basin',
                    tickangle: 0
                },
                yaxis: {
                    title: isRelative ? 'Relative Abundance' : 'Total Sample Count',
                    rangemode: 'tozero'
                },
                legend: {
                    orientation: 'v',
                    x: 1.02,
                    y: 1,
                    xanchor: 'left'
                },
                margin: {
                    l: 60,
                    r: 200,
                    t: 40,
                    b: 80
                },
                height: 500,
                hovermode: 'closest'
            };

            if (isRelative) {
                layout.barnorm = 'fraction';
                layout.yaxis.tickformat = '.0%';
            }

            // Render plot
            const config = {
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            };

            Plotly.newPlot(`plot-${idx}`, traces, layout, config);
        }

        // Render faceted plot
        function renderFacetedPlot(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // For faceted plots, we show all facets but filter genotypes
            const facets = plotData.facets || {};
            const facetNames = Object.keys(facets);

            if (facetNames.length === 0) {
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No facets available.</div>';
                }
                return;
            }

            // Create traces for all facets
            const traces = [];
            const basins = plotData.basins || [];
            const basinLabels = plotData.basin_labels || {};
            const colors = plotData.colors || [];

            facetNames.forEach(facetName => {
                const facet = facets[facetName];
                const genotypes = facet.genotypes || [];

                genotypes.forEach(genotype => {
                    // Check if this genotype is selected
                    if (!state.selectedGenotypes.has(facetName)) return;

                    const genotypeCounts = facet.counts[genotype] || {};
                    const y = basins.map(basin => genotypeCounts[basin] || 0);

                    traces.push({
                        name: facetName,
                        x: basins.map(b => basinLabels[b] || b),
                        y: y,
                        type: 'bar',
                        marker: {
                            color: colors[genotype] || '#888'
                        }
                    });
                });
            });

            if (traces.length === 0) {
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No genotypes selected.</div>';
                }
                return;
            }

            // Layout configuration
            const layout = {
                barmode: 'stack',
                title: '',
                xaxis: {
                    title: 'Ocean Basin',
                    tickangle: 0
                },
                yaxis: {
                    title: 'Sample Count',
                    rangemode: 'tozero'
                },
                legend: {
                    orientation: 'v',
                    x: 1.02,
                    y: 1,
                    xanchor: 'left'
                },
                margin: {
                    l: 60,
                    r: 200,
                    t: 40,
                    b: 80
                },
                height: 600,
                hovermode: 'closest'
            };

            // Render plot
            const config = {
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            };

            Plotly.newPlot(`plot-${idx}`, traces, layout, config);
        }

        // Render distribution map
        function renderDistributionMap(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // Filter genotypes based on selection and threshold
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            if (filteredGenotypes.length === 0) {
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No genotypes selected or meet the minimum sample threshold.</div>';
                }
                return;
            }

            // Prepare Plotly traces (one per genotype)
            const traces = [];
            const locations = plotData.locations || {};
            const colors = plotData.colors || {};

            filteredGenotypes.forEach(genotype => {
                const genoLocations = locations[genotype];
                if (!genoLocations) return;

                let lats = genoLocations.latitudes || [];
                let lons = genoLocations.longitudes || [];
                let sizes = genoLocations.sizes || [];
                const oceanBasins = genoLocations.ocean_basins || [];

                // Filter by ocean basin if basin data is available and basins are selected
                if (oceanBasins.length > 0 && state.selectedBasins.size > 0) {
                    const filteredIndices = [];
                    oceanBasins.forEach((basin, i) => {
                        if (state.selectedBasins.has(basin)) {
                            filteredIndices.push(i);
                        }
                    });

                    // Filter all arrays by selected basins
                    lats = filteredIndices.map(i => lats[i]);
                    lons = filteredIndices.map(i => lons[i]);
                    sizes = filteredIndices.map(i => sizes[i]);
                }

                // Skip if no points remaining after filtering
                if (lats.length === 0) return;

                // Scale marker sizes for better visibility
                const scaledSizes = sizes.map(s => Math.max(5, Math.min(20, s * 3)));

                traces.push({
                    type: 'scattergeo',
                    mode: 'markers',
                    name: genotype,
                    lat: lats,
                    lon: lons,
                    marker: {
                        size: scaledSizes,
                        color: colors[genotype] || '#888',
                        line: {
                            color: 'black',
                            width: 0.5
                        },
                        opacity: 0.85
                    },
                    hovertemplate: '<b>%{text}</b><br>Lat: %{lat:.2f}<br>Lon: %{lon:.2f}<extra></extra>',
                    text: lats.map(() => genotype)
                });
            });

            // Layout configuration
            const layout = {
                title: '',
                geo: {
                    projection: {
                        type: 'natural earth'
                    },
                    showland: true,
                    landcolor: 'rgb(243, 243, 243)',
                    coastlinecolor: 'rgb(204, 204, 204)',
                    showocean: true,
                    oceancolor: 'rgb(230, 245, 255)',
                    showcountries: true,
                    countrycolor: 'rgb(204, 204, 204)',
                    showlakes: true,
                    lakecolor: 'rgb(230, 245, 255)',
                    resolution: 50
                },
                legend: {
                    orientation: 'v',
                    x: 1.02,
                    y: 1,
                    xanchor: 'left'
                },
                margin: {
                    l: 0,
                    r: 200,
                    t: 40,
                    b: 0
                },
                height: 600,
                hovermode: 'closest'
            };

            // Render plot
            const config = {
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            };

            Plotly.newPlot(`plot-${idx}`, traces, layout, config);
        }

        // Download plot as image (PNG or SVG)
        function downloadPlot(idx, format) {
            const plotDiv = document.getElementById(`plot-${idx}`);
            if (!plotDiv) return;

            // Get selected resolution
            const resolutionSelect = document.getElementById(`resolution-${idx}`);
            const resolution = resolutionSelect ? resolutionSelect.value : '1200x900';
            const [width, height] = resolution.split('x').map(Number);

            // Determine filename based on plot type
            const state = plotState[idx];
            const plotType = state.data.plot_type || 'plot';
            const filename = `${plotType}_filtered_${idx}`;

            // Use Plotly's built-in download functionality
            Plotly.downloadImage(plotDiv, {
                format: format,
                width: width,
                height: height,
                filename: filename
            });
        }

        // Download filtered data as CSV
        function downloadCSV(idx) {
            const state = plotState[idx];
            const plotData = state.data;
            const plotType = plotData.plot_type;

            // Filter genotypes based on current selection and threshold
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            if (filteredGenotypes.length === 0) {
                alert('No data to export with current filters');
                return;
            }

            let csvContent = '';

            if (plotType === 'distribution_map') {
                // CSV for distribution map: genotype, latitude, longitude, sample_count
                csvContent = 'Genotype,Latitude,Longitude,Sample_Count\\n';

                filteredGenotypes.forEach(genotype => {
                    const locations = plotData.locations[genotype];
                    if (!locations) return;

                    const lats = locations.latitudes || [];
                    const lons = locations.longitudes || [];
                    const sizes = locations.sizes || [];

                    for (let i = 0; i < lats.length; i++) {
                        csvContent += `"${genotype}",${lats[i]},${lons[i]},${sizes[i]}\n`;
                    }
                });
            } else {
                // CSV for bar plots: genotype, basin, count, sample_total
                csvContent = 'Genotype,Ocean_Basin,Count,Total_Samples\\n';

                const basins = plotData.basins || [];
                const counts = plotData.counts || {};

                filteredGenotypes.forEach(genotype => {
                    const genotypeCounts = counts[genotype] || {};
                    const totalSamples = plotData.sample_counts[genotype] || 0;

                    basins.forEach(basin => {
                        const count = genotypeCounts[basin] || 0;
                        if (count > 0) {
                            csvContent += `"${genotype}","${basin}",${count},${totalSamples}\n`;
                        }
                    });
                });
            }

            // Create and download CSV file
            const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);

            link.setAttribute('href', url);
            link.setAttribute('download', `${plotType}_filtered_data_${idx}.csv`);
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
    </script>
</body>
</html>
"""


class HTMLReportBuilder:
    """
    Builder class for generating comprehensive HTML reports.

    Assembles pipeline results into a formatted HTML document with
    embedded tables, metrics, and visualizations.
    """

    def __init__(self, organism: str, version: str = "1.0.0"):
        """
        Initialize HTML report builder.

        Parameters
        ----------
        organism : str
            Organism name
        version : str, optional
            BOLDGenotyper version
        """
        self.organism = organism
        self.version = version
        self.sections = []
        self.quick_stats = []

    def add_quick_stat(self, value: str, label: str):
        """Add a quick stat to the header bar."""
        self.quick_stats.append({'value': value, 'label': label})

    def add_section(self, content: str):
        """Add a section to the report."""
        self.sections.append(content)

    def render(self) -> str:
        """
        Render the complete HTML report.

        Returns
        -------
        str
            Complete HTML document
        """
        try:
            from jinja2 import Template
        except ImportError:
            raise ImportError(
                "jinja2 is required for HTML report generation. "
                "Install it with: pip install jinja2"
            )

        from datetime import datetime

        template = Template(HTML_REPORT_TEMPLATE)

        return template.render(
            organism=self.organism,
            version=self.version,
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            css=HTML_REPORT_CSS,
            quick_stats=self.quick_stats,
            content='\n'.join(self.sections)
        )


def _format_number(value: float, decimals: int = 2) -> str:
    """Format a number for display."""
    if pd.isna(value):
        return "N/A"
    if isinstance(value, (int, np.integer)):
        return f"{value:,}"
    return f"{value:,.{decimals}f}"


def _format_percentage(value: float, decimals: int = 1) -> str:
    """Format a percentage for display."""
    if pd.isna(value):
        return "N/A"
    return f"{value:.{decimals}f}%"


def _dataframe_to_html(df: pd.DataFrame, max_rows: int = 100) -> str:
    """
    Convert DataFrame to HTML table.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert
    max_rows : int, optional
        Maximum rows to display (default: 100)

    Returns
    -------
    str
        HTML table string
    """
    if df.empty:
        return '<p class="alert alert-info">No data available</p>'

    # Truncate if too long
    truncated = False
    if len(df) > max_rows:
        df = df.head(max_rows)
        truncated = True

    # Convert to HTML
    html = '<div class="table-container">\n'
    html += df.to_html(index=False, classes='', border=0, escape=True)
    html += '</div>\n'

    if truncated:
        html += f'<p style="color: #888; font-size: 0.9em; margin-top: 10px;">Showing first {max_rows} rows of {len(df)} total</p>\n'

    return html


def _build_executive_summary_section(
    assignment_summary: pd.DataFrame,
    annotated_df: pd.DataFrame
) -> str:
    """Build executive summary section with metric cards."""
    if assignment_summary.empty:
        return (
            '<div class="section" id="section-executive-summary">'
            '<h2>Analysis Summary</h2>'
            '<p class="alert alert-warning">No summary data available</p>'
            '</div>'
        )
    
    summary = assignment_summary.iloc[0] if len(assignment_summary) > 0 else {}

    html = '<div class="section" id="section-executive-summary">\n'
    html += '<h2>Analysis Summary</h2>\n'
    html += '<div class="metric-grid">\n'

    # Sample metrics
    html += '<div class="metric-card">\n'
    html += '  <div class="title">Total Samples</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("total_samples", 0), 0)}</div>\n'
    html += '  <div class="description">Samples analyzed</div>\n'
    html += '</div>\n'

    html += '<div class="metric-card">\n'
    html += '  <div class="title">With Sequence Data</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("samples_with_sequence", 0), 0)}</div>\n'
    pct = (summary.get("samples_with_sequence", 0) / summary.get("total_samples", 1) * 100) if summary.get("total_samples", 0) > 0 else 0
    html += f'  <div class="description">{_format_percentage(pct)} of total</div>\n'
    html += '</div>\n'

    # Genotype metrics
    html += '<div class="metric-card">\n'
    html += '  <div class="title">Consensus Groups</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_consensus_groups", 0), 0)}</div>\n'
    html += '  <div class="description">Unique genotypes identified</div>\n'
    html += '</div>\n'

    html += '<div class="metric-card">\n'
    html += '  <div class="title">Species Identified</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_species_identified", 0), 0)}</div>\n'
    html += '  <div class="description">Distinct species</div>\n'
    html += '</div>\n'

    # Assignment metrics
    html += '<div class="metric-card">\n'
    html += '  <div class="title">Successfully Assigned</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_assigned", 0), 0)}</div>\n'
    html += f'  <div class="description">{_format_percentage(summary.get("pct_assigned", 0))} assignment rate</div>\n'
    html += '</div>\n'

    html += '<div class="metric-card">\n'
    html += '  <div class="title">Mean Identity</div>\n'
    html += f'  <div class="value">{_format_percentage(summary.get("mean_identity_assigned", 0) * 100, 1)}</div>\n'
    html += '  <div class="description">For assigned samples</div>\n'
    html += '</div>\n'

    # Geographic coverage
    if 'ocean_basin' in annotated_df.columns:
        n_basins = annotated_df['ocean_basin'].nunique()
        html += '<div class="metric-card">\n'
        html += '  <div class="title">Ocean Basins</div>\n'
        html += f'  <div class="value">{n_basins}</div>\n'
        html += '  <div class="description">Geographic coverage</div>\n'
        html += '</div>\n'

    # Taxonomy conflicts
    if summary.get("n_taxonomy_conflicts", 0) > 0:
        html += '<div class="metric-card">\n'
        html += '  <div class="title">Taxonomy Conflicts</div>\n'
        html += f'  <div class="value">{_format_number(summary.get("n_taxonomy_conflicts", 0), 0)}</div>\n'
        html += f'  <div class="description">{_format_percentage(summary.get("pct_taxonomy_conflicts", 0))} of assigned</div>\n'
        html += '</div>\n'

    html += '</div>\n'  # Close metric-grid
    html += '</div>\n'  # Close section

    return html


def _build_assignment_section(
    assignment_summary: pd.DataFrame,
    diagnostics_df: pd.DataFrame
) -> str:
    """Build genotype assignment results section."""
    html = '<div class="section" id="section-assignment">\n'
    html += '<h2>Genotype Assignment Results</h2>\n'

    if assignment_summary.empty:
        html += '<p class="alert alert-warning">No assignment summary data available</p>\n'
        html += '</div>\n'
        return html

    summary = assignment_summary.iloc[0]

    # Assignment Status Breakdown
    html += '<h3>Assignment Status Breakdown</h3>\n'
    html += '<div class="table-container">\n'
    html += '<table>\n'
    html += '<thead><tr><th>Status</th><th>Count</th><th>Percentage</th></tr></thead>\n'
    html += '<tbody>\n'

    statuses = [
        ('assigned', 'Successfully Assigned', 'badge-success'),
        ('low_confidence', 'Low Confidence', 'badge-warning'),
        ('tie', 'Tied Assignment', 'badge-warning'),
        ('below_threshold', 'Below Threshold', 'badge-danger'),
        ('no_sequence', 'No Sequence Data', 'badge-info')
    ]

    for status_key, status_label, badge_class in statuses:
        count = summary.get(f'n_{status_key}', 0)
        pct = summary.get(f'pct_{status_key}', 0)
        html += f'<tr><td><span class="badge {badge_class}">{status_label}</span></td>'
        html += f'<td>{_format_number(count, 0)}</td>'
        html += f'<td>{_format_percentage(pct)}</td></tr>\n'

    html += '</tbody>\n</table>\n</div>\n'

    # Identity Score Statistics
    html += '<h3>Identity Score Statistics (Assigned Samples)</h3>\n'
    html += '<div class="table-container">\n'
    html += '<table>\n'
    html += '<thead><tr><th>Metric</th><th>Value</th></tr></thead>\n'
    html += '<tbody>\n'
    html += f'<tr><td>Mean Identity</td><td>{_format_percentage(summary.get("mean_identity_assigned", 0) * 100, 2)}</td></tr>\n'
    html += f'<tr><td>Median Identity</td><td>{_format_percentage(summary.get("median_identity_assigned", 0) * 100, 2)}</td></tr>\n'
    html += f'<tr><td>Minimum Identity</td><td>{_format_percentage(summary.get("min_identity_assigned", 0) * 100, 2)}</td></tr>\n'
    html += f'<tr><td>Maximum Identity</td><td>{_format_percentage(summary.get("max_identity_assigned", 0) * 100, 2)}</td></tr>\n'
    html += '</tbody>\n</table>\n</div>\n'

    html += '</div>\n'  # Close section
    return html


def _build_taxonomy_section(
    taxonomy_dir: Path,
    organism: str
) -> str:
    """Build taxonomy section with consensus taxonomy and species composition."""
    html = '<div class="section" id="section-taxonomy">\n'
    html += '<h2>Taxonomy</h2>\n'

    # Load consensus taxonomy
    consensus_tax_file = taxonomy_dir / f"{organism}_consensus_taxonomy.csv"
    species_by_cons_file = taxonomy_dir / f"{organism}_species_by_consensus.csv"

    if consensus_tax_file.exists():
        html += '<h3>Consensus Group Taxonomy</h3>\n'
        try:
            tax_df = pd.read_csv(consensus_tax_file)
            # Select key columns
            display_cols = [c for c in ['consensus_group', 'assigned_sp', 'assignment_level',
                                        'assignment_notes', 'majority_fraction'] if c in tax_df.columns]
            if display_cols:
                html += _dataframe_to_html(tax_df[display_cols], max_rows=50)
            else:
                html += _dataframe_to_html(tax_df, max_rows=50)
        except Exception as e:
            html += f'<p class="alert alert-warning">Could not load consensus taxonomy: {e}</p>\n'
    else:
        html += '<p class="alert alert-info">Consensus taxonomy file not found</p>\n'

    if species_by_cons_file.exists():
        html += '<h3>Species Composition by Consensus Group</h3>\n'
        try:
            species_df = pd.read_csv(species_by_cons_file)
            html += _dataframe_to_html(species_df, max_rows=50)
        except Exception as e:
            html += f'<p class="alert alert-warning">Could not load species composition: {e}</p>\n'
    else:
        html += '<p class="alert alert-info">Species composition file not found</p>\n'

    html += '</div>\n'  # Close section
    return html


def _build_parameters_section(output_dir: Path, organism: str) -> str:
    """Build pipeline parameters section."""
    html = '<div class="section" id="section-parameters">\n'
    html += '<h2>Pipeline Parameters</h2>\n'

    # Try to load parameters from JSON file
    params_file = output_dir / f"{organism}_pipeline_parameters.json"
    if not params_file.exists():
        html += '<p class="alert alert-info">Parameter information not available</p>\n'
        html += '</div>\n'
        return html

    try:
        import json
        with open(params_file, 'r') as f:
            params = json.load(f)

        html += '<p>The following parameters were used for this analysis:</p>\n'
        html += '<table class="data-table">\n'
        html += '<thead><tr><th>Parameter</th><th>Value</th><th>Description</th></tr></thead>\n'
        html += '<tbody>\n'

        # Clustering threshold
        ct = params.get('clustering_threshold', 'N/A')
        html += '<tr>\n'
        html += f'<td><strong>Clustering Threshold</strong></td>\n'
        html += f'<td><code>{ct}</code>'
        if isinstance(ct, (int, float)):
            html += f' ({(1-ct)*100:.1f}% identity)'
        html += '</td>\n'
        html += '<td>Maximum genetic distance for grouping sequences into consensus genotypes. '
        html += 'Lower values create more groups with tighter genetic similarity.</td>\n'
        html += '</tr>\n'

        # Similarity threshold
        st = params.get('similarity_threshold', 'N/A')
        html += '<tr>\n'
        html += f'<td><strong>Similarity Threshold</strong></td>\n'
        html += f'<td><code>{st}</code>'
        if isinstance(st, (int, float)):
            html += f' ({st*100:.0f}% identity)'
        html += '</td>\n'
        html += '<td>Minimum sequence identity required for assigning samples to genotypes. '
        html += 'Samples below this threshold are marked as unassigned.</td>\n'
        html += '</tr>\n'

        # Tie margin
        tm = params.get('tie_margin', 'N/A')
        html += '<tr>\n'
        html += f'<td><strong>Tie Margin</strong></td>\n'
        html += f'<td><code>{tm}</code>'
        if isinstance(tm, (int, float)):
            html += f' ({tm*100:.1f}% difference)'
        html += '</td>\n'
        html += '<td>Maximum identity difference between top matches to flag as ambiguous. '
        html += 'Samples with (best - runner-up) &lt; tie margin are flagged for manual review.</td>\n'
        html += '</tr>\n'

        # Tie threshold
        tt = params.get('tie_threshold', 'N/A')
        html += '<tr>\n'
        html += f'<td><strong>Tie Threshold</strong></td>\n'
        html += f'<td><code>{tt}</code>'
        if isinstance(tt, (int, float)):
            html += f' ({tt*100:.0f}% identity)'
        html += '</td>\n'
        html += '<td>Minimum best-match identity required to consider tie detection. '
        html += 'Prevents flagging low-quality matches as ties.</td>\n'
        html += '</tr>\n'

        # Additional parameters
        html += '<tr>\n'
        html += f'<td><strong>Threads</strong></td>\n'
        html += f'<td><code>{params.get("threads", "N/A")}</code></td>\n'
        html += '<td>Number of parallel processing threads used.</td>\n'
        html += '</tr>\n'

        html += '<tr>\n'
        html += f'<td><strong>Phylogenetic Tree</strong></td>\n'
        html += f'<td><code>{params.get("build_tree", False)}</code></td>\n'
        html += '<td>Whether phylogenetic tree was constructed.</td>\n'
        html += '</tr>\n'

        html += '</tbody>\n'
        html += '</table>\n'

    except Exception as e:
        html += f'<p class="alert alert-warning">Could not load parameters: {e}</p>\n'

    html += '</div>\n'
    return html


def _build_geographic_section(annotated_df: pd.DataFrame) -> str:
    """Build geographic distribution section."""
    html = '<div class="section" id="section-geographic">\n'
    html += '<h2>Geographic Distribution</h2>\n'

    if 'ocean_basin' not in annotated_df.columns:
        html += '<p class="alert alert-info">No geographic data available</p>\n'
        html += '</div>\n'
        return html

    # Count total samples and samples with Unknown geography
    total_samples = len(annotated_df)
    unknown_geography = (annotated_df['ocean_basin'].fillna('Unknown') == 'Unknown').sum()
    samples_with_geography = total_samples - unknown_geography

    # Add summary information
    html += '<div class="alert alert-info">\n'
    html += f'<strong>Geographic Analysis Summary:</strong> {_format_number(samples_with_geography, 0)} of '
    html += f'{_format_number(total_samples, 0)} samples '
    html += f'({_format_percentage(samples_with_geography / total_samples * 100)}) have defined geographic locations. '
    if unknown_geography > 0:
        html += f'{_format_number(unknown_geography, 0)} samples '
        html += f'({_format_percentage(unknown_geography / total_samples * 100)}) '
        html += 'with unknown or missing geography were excluded from geographic analyses.\n'
    html += '</div>\n'

    # Basin distribution (excluding Unknown)
    html += '<h3>Sample Distribution by Ocean Basin</h3>\n'

    # Filter out Unknown geography for the table
    basin_series = annotated_df['ocean_basin'].fillna('Unknown')
    basin_counts = basin_series[basin_series != 'Unknown'].value_counts().reset_index()

    if len(basin_counts) > 0:
        basin_counts.columns = ['Ocean Basin', 'Sample Count']
        # Calculate percentages based on samples with defined geography
        basin_counts['Percentage'] = (basin_counts['Sample Count'] / samples_with_geography * 100).round(1)
        html += _dataframe_to_html(basin_counts, max_rows=20)
    else:
        html += '<p class="alert alert-warning">No samples with defined ocean basin assignments</p>\n'

    # Genotype × Basin distribution (excluding Unknown)
    if 'consensus_group' in annotated_df.columns:
        html += '<h3>Genotypes per Ocean Basin</h3>\n'

        # Filter out Unknown geography before creating crosstab
        df_with_geography = annotated_df[annotated_df['ocean_basin'].fillna('Unknown') != 'Unknown'].copy()

        if len(df_with_geography) > 0:
            # Create crosstab
            ct = pd.crosstab(
                df_with_geography['consensus_group'].fillna('Unassigned'),
                df_with_geography['ocean_basin']
            )

            # Convert to long format for better display
            ct_reset = ct.reset_index()
            html += _dataframe_to_html(ct_reset, max_rows=30)
        else:
            html += '<p class="alert alert-warning">No genotyped samples with defined ocean basin assignments</p>\n'

    # Missing geography
    if 'lat' in annotated_df.columns and 'lon' in annotated_df.columns:
        missing_coords = annotated_df[annotated_df[['lat', 'lon']].isna().any(axis=1)]
        n_missing = len(missing_coords)
        pct_missing = (n_missing / len(annotated_df) * 100)

        if n_missing > 0:
            html += '<div class="alert alert-warning">\n'
            html += f'<strong>Missing Geography:</strong> {_format_number(n_missing, 0)} samples '
            html += f'({_format_percentage(pct_missing)}) do not have valid coordinate data.\n'
            html += '</div>\n'

    html += '</div>\n'  # Close section
    return html


def _encode_image_to_base64(image_path: Path) -> Optional[str]:
    """
    Encode an image file to base64 for embedding in HTML.

    Parameters
    ----------
    image_path : Path
        Path to image file

    Returns
    -------
    Optional[str]
        Base64-encoded image data URI, or None if encoding failed
    """
    import base64

    if not image_path.exists():
        return None

    try:
        with open(image_path, 'rb') as f:
            image_data = f.read()

        # Determine MIME type from extension
        ext = image_path.suffix.lower()
        mime_types = {
            '.png': 'image/png',
            '.jpg': 'image/jpeg',
            '.jpeg': 'image/jpeg',
            '.gif': 'image/gif',
            '.svg': 'image/svg+xml'
        }
        mime_type = mime_types.get(ext, 'image/png')

        # Encode to base64
        encoded = base64.b64encode(image_data).decode('utf-8')
        return f"data:{mime_type};base64,{encoded}"

    except Exception as e:
        logger.warning(f"Failed to encode image {image_path}: {e}")
        return None


def _build_visualizations_section(
    output_dir: Path,
    organism: str
) -> str:
    """
    Build visualizations section with embedded images.

    Parameters
    ----------
    output_dir : Path
        Base output directory
    organism : str
        Organism name

    Returns
    -------
    str
        HTML content for visualizations section
    """
    html = '<div class="section" id="section-visualizations">\n'
    html += '<h2>Visualizations</h2>\n'

    viz_dir = output_dir / 'visualization'

    if not viz_dir.exists():
        html += '<p class="alert alert-info">No visualizations available</p>\n'
        html += '</div>\n'
        return html

    # Define visualization categories and their display names (in display order)
    viz_categories = [
        {
            'title': 'Identity Distribution',
            'pattern': f'{organism}_identity_distribution.png',
            'description': 'Distribution of sequence identity scores for assigned samples'
        },
        {
            'title': 'Phylogenetic Tree',
            'pattern': f'{organism}_tree.png',
            'description': 'Phylogenetic tree showing relationships between consensus groups'
        },
        {
            'title': 'Relative Abundance by Ocean Basin',
            'pattern': f'{organism}_distribution_bar.png',
            'json_pattern': f'{organism}_distribution_bar_data.json',
            'description': 'Relative abundance of genotypes across ocean basins'
        },
        {
            'title': 'Total Abundance by Ocean Basin',
            'pattern': f'{organism}_totaldistribution_bar.png',
            'json_pattern': f'{organism}_totaldistribution_bar_data.json',
            'description': 'Total sample counts of genotypes across ocean basins'
        },
        {
            'title': 'Total Abundance by Ocean Basin (Faceted)',
            'pattern': f'{organism}_distribution_bar_faceted.png',
            'json_pattern': f'{organism}_distribution_bar_faceted_data.json',
            'description': 'Total sample counts faceted by species or genotype'
        },
        {
            'title': 'Distribution Map',
            'pattern': f'{organism}_distribution_map.png',
            'json_pattern': f'{organism}_distribution_map_data.json',
            'description': 'Geographic distribution of samples'
        },
        {
            'title': 'Distribution Map (Faceted)',
            'pattern': f'{organism}_distribution_map_faceted.png',
            'description': 'Geographic distribution faceted by species or genotype'
        }
    ]

    # Filter to available visualizations
    available_viz = []
    for viz in viz_categories:
        image_path = viz_dir / viz['pattern']
        if image_path.exists():
            encoded_image = _encode_image_to_base64(image_path)
            if encoded_image:
                viz_data = {**viz, 'encoded_image': encoded_image}

                # Load JSON data if available (for interactive plots)
                if 'json_pattern' in viz:
                    json_path = viz_dir / viz['json_pattern']
                    if json_path.exists():
                        try:
                            with open(json_path, 'r') as f:
                                plot_data = json.load(f)
                            viz_data['plot_data'] = plot_data
                            logger.debug(f"Loaded plot data for {viz['title']}: {json_path}")
                        except Exception as e:
                            logger.warning(f"Failed to load plot data from {json_path}: {e}")

                available_viz.append(viz_data)

    if not available_viz:
        html += '<p class="alert alert-info">No visualization images found</p>\n'
        html += '</div>\n'
        return html

    # Create subtabs for visualizations
    html += '<div class="subtabs">\n'
    for idx, viz in enumerate(available_viz):
        active_class = ' active' if idx == 0 else ''
        tab_id = f"viz-tab-{idx}"
        html += f'<button class="subtab-button{active_class}" onclick="switchSubtab(\'{tab_id}\')">{viz["title"]}</button>\n'
    html += '</div>\n'

    # Create subtab content
    for idx, viz in enumerate(available_viz):
        active_class = ' active' if idx == 0 else ''
        tab_id = f"viz-tab-{idx}"
        html += f'<div id="{tab_id}" class="subtab-content{active_class}">\n'
        html += '<div class="viz-container">\n'
        html += f'<h3>{viz["title"]}</h3>\n'
        html += f'<p class="viz-description">{viz["description"]}</p>\n'

        # Add interactive controls for plots with JSON data (except faceted plots)
        if 'plot_data' in viz and viz['plot_data'].get('plot_type') != 'faceted_abundance':
            data_id = f"plot-data-{idx}"
            plot_id = f"plot-{idx}"
            controls_id = f"controls-{idx}"

            # Toggle button between static and interactive
            html += '<div class="plot-mode-toggle">\n'
            html += f'<button class="toggle-btn active" onclick="showStaticPlot(\'{plot_id}\', {idx})">Static Image</button>\n'
            html += f'<button class="toggle-btn" onclick="showInteractivePlot(\'{plot_id}\', {idx})">Interactive Plot</button>\n'
            html += '</div>\n'

            # Static plot (shown by default)
            html += f'<div id="{plot_id}-static" class="plot-view active">\n'
            html += f'<img src="{viz["encoded_image"]}" alt="{viz["title"]}" class="viz-image">\n'
            html += '</div>\n'

            # Interactive plot container (hidden by default)
            html += f'<div id="{plot_id}-interactive" class="plot-view" style="display: none;">\n'

            # Filter controls
            html += f'<div id="{controls_id}" class="plot-controls">\n'

            # Summary Statistics Panel
            html += f'<div id="stats-panel-{idx}" class="stats-panel">\n'
            html += '<div class="stats-grid">\n'
            html += '<div class="stat-item">\n'
            html += '<div class="stat-label">Genotypes Shown</div>\n'
            html += f'<div class="stat-value" id="stat-genotypes-{idx}">--</div>\n'
            html += '</div>\n'
            html += '<div class="stat-item">\n'
            html += '<div class="stat-label">Total Samples</div>\n'
            html += f'<div class="stat-value" id="stat-samples-{idx}">--</div>\n'
            html += '</div>\n'
            html += '<div class="stat-item">\n'
            html += '<div class="stat-label">% of Total Data</div>\n'
            html += f'<div class="stat-value" id="stat-percentage-{idx}">--</div>\n'
            html += '</div>\n'
            html += '</div>\n'
            html += '</div>\n'

            html += '<div class="controls-section">\n'
            html += '<h4>Filter Options</h4>\n'

            # Threshold slider
            html += '<div class="control-group">\n'
            html += '<label for="threshold-{0}">Minimum Sample Count: <span id="threshold-value-{0}">0</span></label>\n'.format(idx)
            html += '<input type="range" id="threshold-{0}" class="threshold-slider" min="0" max="50" value="0" \n'.format(idx)
            html += f'oninput="updateThreshold({idx}); updatePlot({idx})">\n'
            html += '</div>\n'

            # Ocean basin filter (only for distribution maps with ocean basin data)
            if viz['plot_data'].get('plot_type') == 'distribution_map' and viz['plot_data'].get('ocean_basins'):
                html += '<div class="control-group">\n'
                html += '<label>Filter by Ocean Basin:</label>\n'
                html += f'<div class="checkbox-controls">\n'
                html += f'<button class="btn-small" onclick="selectAllBasins({idx})">Select All</button>\n'
                html += f'<button class="btn-small" onclick="deselectAllBasins({idx})">Deselect All</button>\n'
                html += '</div>\n'
                html += f'<div id="basin-checkboxes-{idx}" class="genotype-checkboxes">\n'
                html += '<!-- Basin checkboxes will be populated by JavaScript -->\n'
                html += '</div>\n'
                html += '</div>\n'

            # Genotype checkboxes
            html += '<div class="control-group">\n'
            html += '<label>Select Genotypes:</label>\n'
            html += f'<div class="checkbox-controls">\n'
            html += f'<button class="btn-small" onclick="selectAllGenotypes({idx})">Select All</button>\n'
            html += f'<button class="btn-small" onclick="deselectAllGenotypes({idx})">Deselect All</button>\n'
            html += '</div>\n'
            html += f'<div id="genotype-checkboxes-{idx}" class="genotype-checkboxes">\n'
            html += '<!-- Checkboxes will be populated by JavaScript -->\n'
            html += '</div>\n'
            html += '</div>\n'

            # Download options
            html += '<div class="control-group">\n'
            html += '<label>Export Options:</label>\n'
            html += '<div class="download-options">\n'
            html += f'<button class="btn-download" onclick="downloadPlot({idx}, \'png\')">📥 PNG</button>\n'
            html += f'<button class="btn-download" onclick="downloadPlot({idx}, \'svg\')">📥 SVG</button>\n'
            html += f'<button class="btn-download" onclick="downloadCSV({idx})">📥 CSV Data</button>\n'
            html += '</div>\n'
            html += '<div class="resolution-control">\n'
            html += '<label for="resolution-{0}" class="small-label">Image Resolution:</label>\n'.format(idx)
            html += '<select id="resolution-{0}" class="resolution-select">\n'.format(idx)
            html += '<option value="800x600">Standard (800×600)</option>\n'
            html += '<option value="1200x900" selected>High (1200×900)</option>\n'
            html += '<option value="1600x1200">Print (1600×1200)</option>\n'
            html += '<option value="2400x1800">Publication (2400×1800)</option>\n'
            html += '</select>\n'
            html += '</div>\n'
            html += '</div>\n'

            html += '</div>\n'  # controls-section
            html += '</div>\n'  # plot-controls

            # Plotly chart container
            html += f'<div id="{plot_id}" class="plotly-chart"></div>\n'

            html += '</div>\n'  # plot-view interactive

            # Embed plot data as JSON
            json_str = json.dumps(viz['plot_data'])
            html += f'<script type="application/json" id="{data_id}">\n'
            html += json_str
            html += '\n</script>\n'
        else:
            # No interactive data or faceted plot - just show static image
            html += f'<img src="{viz["encoded_image"]}" alt="{viz["title"]}" class="viz-image">\n'

        html += '</div>\n'  # viz-container
        html += '</div>\n'  # subtab-content

    html += '</div>\n'
    return html


def generate_html_report(
    organism: str,
    output_dir: Path,
    version: str = "1.0.0"
) -> Optional[Path]:
    """
    Generate comprehensive HTML summary report.

    Aggregates all pipeline outputs into a single, formatted HTML document
    with tables, metrics, and analysis summaries.

    Parameters
    ----------
    organism : str
        Organism name
    output_dir : Path
        Base output directory containing reports/, taxonomy/, etc.
    version : str, optional
        BOLDGenotyper version (default: "1.0.0")

    Returns
    -------
    Optional[Path]
        Path to generated HTML report, or None if generation failed

    Examples
    --------
    >>> from pathlib import Path
    >>> report_path = generate_html_report(
    ...     organism="Sphyrna_lewini",
    ...     output_dir=Path("data/Sphyrna_lewini_output")
    ... )
    >>> print(f"Report generated: {report_path}")
    """
    logger.info(f"Generating HTML summary report for {organism}...")

    try:
        # Initialize builder
        builder = HTMLReportBuilder(organism=organism, version=version)

        # Define paths
        reports_dir = output_dir / 'reports'
        taxonomy_dir = output_dir / 'taxonomy'

        assignment_summary_file = reports_dir / f"{organism}_assignment_summary.csv"
        annotated_file = output_dir / f"{organism}_annotated.csv"
        diagnostics_file = output_dir / 'genotype_assignments' / f"{organism}_diagnostics.csv"

        # Load data
        assignment_summary = pd.DataFrame()
        if assignment_summary_file.exists():
            assignment_summary = pd.read_csv(assignment_summary_file)

        annotated_df = pd.DataFrame()
        if annotated_file.exists():
            try:
                annotated_df = pd.read_csv(annotated_file, low_memory=False)
            except Exception as e:
                logger.warning(f"Could not load annotated CSV: {e}")

        diagnostics_df = pd.DataFrame()
        if diagnostics_file.exists():
            diagnostics_df = pd.read_csv(diagnostics_file)

        # Add quick stats
        if not assignment_summary.empty:
            summary = assignment_summary.iloc[0]
            builder.add_quick_stat(
                _format_number(summary.get('total_samples', 0), 0),
                'Total Samples'
            )
            builder.add_quick_stat(
                _format_number(summary.get('n_consensus_groups', 0), 0),
                'Genotypes'
            )
            builder.add_quick_stat(
                _format_percentage(summary.get('pct_assigned', 0)),
                'Assignment Rate'
            )

        # Build sections
        builder.add_section(_build_executive_summary_section(assignment_summary, annotated_df))
        builder.add_section(_build_parameters_section(output_dir, organism))
        builder.add_section(_build_assignment_section(assignment_summary, diagnostics_df))
        builder.add_section(_build_taxonomy_section(taxonomy_dir, organism))
        builder.add_section(_build_geographic_section(annotated_df))
        builder.add_section(_build_visualizations_section(output_dir, organism))

        # Render HTML
        html_content = builder.render()

        # Save report
        output_file = output_dir / f"{organism}_summary_report.html"
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        logger.info(f"HTML report saved: {output_file}")
        return output_file

    except Exception as e:
        logger.error(f"Failed to generate HTML report: {e}", exc_info=True)
        return None
