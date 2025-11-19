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
        padding: 20px;
    }

    .container {
        max-width: 1200px;
        margin: 0 auto;
        background-color: white;
        box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        border-radius: 8px;
        overflow: hidden;
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
        padding: 40px;
        text-align: center;
    }

    .header h1 {
        font-size: 2.5em;
        margin-bottom: 10px;
        font-style: italic;
    }

    .header .subtitle {
        font-size: 1.1em;
        opacity: 0.9;
        margin-bottom: 5px;
    }

    .header .timestamp {
        font-size: 0.9em;
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
    
    .section-tabs {
        margin: 20px 0 10px 0;
        display: flex;
        flex-wrap: wrap;
        border-bottom: 1px solid var(--border-color);
        gap: 8px;
    }

    .section-tab {
        padding: 8px 16px;
        text-decoration: none;
        color: var(--primary-color);
        border: 1px solid var(--border-color);
        border-bottom: none;
        border-radius: 6px 6px 0 0;
        background-color: var(--light-bg);
        font-size: 0.95em;
    }

    .section-tab:hover {
        background-color: white;
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
        .header h1 {
            font-size: 1.8em;
        }
        .quick-stats {
            flex-direction: column;
        }
        .metric-grid {
            grid-template-columns: 1fr;
        }
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

    @media print {
        .viz-container {
            page-break-inside: avoid;
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
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>{{ organism }}</h1>
            <div class="subtitle">BOLDGenotyper Analysis Report</div>
            <div class="timestamp">Generated: {{ timestamp }}</div>
            <div class="timestamp">BOLDGenotyper v{{ version }}</div>
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

        <!-- Section Navigation -->
        <div class="section-tabs">
            <a href="#section-executive-summary" class="section-tab">Executive Summary</a>
            <a href="#section-assignment" class="section-tab">Genotype Assignment</a>
            <a href="#section-taxonomy" class="section-tab">Taxonomy</a>
            <a href="#section-geographic" class="section-tab">Geographic Distribution</a>
            <a href="#section-visualizations" class="section-tab">Visualizations</a>
        </div>

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


def _build_geographic_section(annotated_df: pd.DataFrame) -> str:
    """Build geographic distribution section."""
    html = '<div class="section" id="section-geographic">\n'
    html += '<h2>Geographic Distribution</h2>\n'

    if 'ocean_basin' not in annotated_df.columns:
        html += '<p class="alert alert-info">No geographic data available</p>\n'
        html += '</div>\n'
        return html

    # Basin distribution
    html += '<h3>Sample Distribution by Ocean Basin</h3>\n'

    basin_counts = annotated_df['ocean_basin'].fillna('Unknown').value_counts().reset_index()
    basin_counts.columns = ['Ocean Basin', 'Sample Count']
    basin_counts['Percentage'] = (basin_counts['Sample Count'] / len(annotated_df) * 100).round(1)

    html += _dataframe_to_html(basin_counts, max_rows=20)

    # Genotype × Basin distribution
    if 'consensus_group' in annotated_df.columns:
        html += '<h3>Genotypes per Ocean Basin</h3>\n'

        # Create crosstab
        ct = pd.crosstab(
            annotated_df['consensus_group'].fillna('Unassigned'),
            annotated_df['ocean_basin'].fillna('Unknown')
        )

        # Convert to long format for better display
        ct_reset = ct.reset_index()
        html += _dataframe_to_html(ct_reset, max_rows=30)

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

    # Define visualization categories and their display names
    viz_categories = [
        {
            'title': 'Phylogenetic Tree',
            'pattern': f'{organism}_tree.png',
            'description': 'Phylogenetic tree showing relationships between consensus groups'
        },
        {
            'title': 'Identity Distribution',
            'pattern': f'{organism}_identity_distribution.png',
            'description': 'Distribution of sequence identity scores for assigned samples'
        },
        {
            'title': 'Distribution Map',
            'pattern': f'{organism}_distribution_map.png',
            'description': 'Geographic distribution of samples'
        },
        {
            'title': 'Distribution Map (Faceted)',
            'pattern': f'{organism}_distribution_map_faceted.png',
            'description': 'Geographic distribution faceted by species or genotype'
        },
        {
            'title': 'Relative Abundance by Ocean Basin',
            'pattern': f'{organism}_distribution_bar.png',
            'description': 'Relative abundance of genotypes across ocean basins'
        },
        {
            'title': 'Relative Abundance by Ocean Basin (Faceted)',
            'pattern': f'{organism}_distribution_bar_faceted.png',
            'description': 'Relative abundance faceted by species or genotype'
        },
        {
            'title': 'Total Abundance by Ocean Basin',
            'pattern': f'{organism}_totaldistribution_bar.png',
            'description': 'Total sample counts of genotypes across ocean basins'
        }
    ]

    # Track whether we found any visualizations
    found_any = False

    for viz in viz_categories:
        image_path = viz_dir / viz['pattern']

        if not image_path.exists():
            continue

        found_any = True

        # Encode image
        encoded_image = _encode_image_to_base64(image_path)

        if encoded_image:
            html += '<div class="viz-container">\n'
            html += f'<h3>{viz["title"]}</h3>\n'
            html += f'<p class="viz-description">{viz["description"]}</p>\n'
            html += f'<img src="{encoded_image}" alt="{viz["title"]}" class="viz-image">\n'
            html += '</div>\n'
        else:
            html += '<div class="alert alert-warning">\n'
            html += f'<strong>{viz["title"]}:</strong> Image file found but could not be loaded\n'
            html += '</div>\n'

    if not found_any:
        html += '<p class="alert alert-info">No visualization images found</p>\n'

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
        builder.add_section(_build_assignment_section(assignment_summary, diagnostics_df))
        builder.add_section(_build_taxonomy_section(taxonomy_dir, organism))
        builder.add_section(_build_geographic_section(annotated_df))
        builder.add_section(_build_visualizations_section(output_dir, organism))

        # Render HTML
        html_content = builder.render()

        # Save report
        output_file = reports_dir / f"{organism}_summary_report.html"
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        logger.info(f"HTML report saved: {output_file}")
        return output_file

    except Exception as e:
        logger.error(f"Failed to generate HTML report: {e}", exc_info=True)
        return None
