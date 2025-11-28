"""
Parameter Sweep Module for BOLDGenotyper

This module helps users optimize clustering thresholds by testing multiple values
and visualizing stability. It identifies the "elbow point" where grouping stabilizes
and provides automated recommendations.

Key Features:
- Test multiple threshold values efficiently
- Track sample movement between thresholds
- Generate stability plots and elbow curves
- Provide automated recommendations
- Reuse alignments to save time

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from datetime import datetime
import subprocess
import shutil

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


class ParameterSweepError(Exception):
    """Base exception for parameter sweep errors."""
    pass


def run_single_threshold(
    tsv_path: Path,
    organism: str,
    threshold: float,
    output_dir: Path,
    threads: int = 4,
    keep_intermediates: bool = False
) -> Dict[str, Any]:
    """
    Run BOLDGenotyper pipeline with a single threshold value.

    Parameters
    ----------
    tsv_path : Path
        Input TSV file
    organism : str
        Organism name
    threshold : float
        Clustering threshold to test
    output_dir : Path
        Output directory for this run
    threads : int
        Number of threads
    keep_intermediates : bool
        Whether to keep full output

    Returns
    -------
    dict
        Summary metrics for this threshold
    """
    from . import cli, config

    logger.info(f"  Testing threshold {threshold}...")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load and configure pipeline
    cfg = config.get_default_config()
    cfg = cfg.update(
        dereplication__clustering_threshold=threshold,
        genotype_assignment__min_identity=0.5,
        n_threads=threads,
        output_dir=output_dir,
        log_level='WARNING',  # Reduce verbosity
        phylogenetic__build_tree=False,  # Skip tree for speed
        keep_intermediates=keep_intermediates
    )

    # Run pipeline with minimal output
    try:
        success = cli.run_pipeline(
            tsv_path=tsv_path,
            organism=organism,
            output_dir=output_dir,
            cfg=cfg,
            no_report=True,  # Skip HTML report
            skip_geo=True,   # Skip geographic analysis
            export_plot_data=False
        )

        if not success:
            logger.warning(f"  Pipeline failed for threshold {threshold}")
            return None

    except Exception as e:
        logger.error(f"  Error running threshold {threshold}: {e}")
        return None

    # Extract metrics from output
    metrics = extract_run_metrics(output_dir, organism, threshold)

    # Clean up if not keeping intermediates
    if not keep_intermediates:
        cleanup_intermediate_files(output_dir)

    return metrics


def extract_run_metrics(
    output_dir: Path,
    organism: str,
    threshold: float
) -> Dict[str, Any]:
    """
    Extract summary metrics from a completed run.

    Parameters
    ----------
    output_dir : Path
        Output directory
    organism : str
        Organism name
    threshold : float
        Threshold used

    Returns
    -------
    dict
        Summary metrics
    """
    metrics = {'threshold': threshold}

    # Load annotated data
    annotated_path = output_dir / f"{organism}_annotated.csv"
    if not annotated_path.exists():
        return metrics

    df = pd.read_csv(annotated_path)

    # Basic metrics
    metrics['n_samples'] = len(df)

    if 'consensus_group' in df.columns:
        # Number of groups
        n_groups = df['consensus_group'].nunique()
        metrics['n_groups'] = n_groups

        # Group sizes
        group_sizes = df.groupby('consensus_group').size()
        metrics['n_singletons'] = (group_sizes == 1).sum()
        metrics['mean_group_size'] = group_sizes.mean()
        metrics['median_group_size'] = group_sizes.median()
        metrics['max_group_size'] = group_sizes.max()

    # Identity metrics - read from diagnostics file
    diagnostics_dir = output_dir / "genotype_assignments"
    # Accept both legacy and current diagnostics filenames
    diagnostics_candidates = [
        diagnostics_dir / f"{organism}_assignment_diagnostics.csv",
        diagnostics_dir / f"{organism}_diagnostics.csv"
    ]
    diagnostics_path = next((p for p in diagnostics_candidates if p.exists()), None)

    if diagnostics_path:
        try:
            diag_df = pd.read_csv(diagnostics_path)
            if 'identity' in diag_df.columns:
                # Calculate mean and median identity for assigned samples
                assigned_diag = diag_df[diag_df['identity'].notna()]
                if len(assigned_diag) > 0:
                    identities = assigned_diag['identity']
                    # Identity is stored as proportion (0-1); scale to percent for reporting
                    if identities.max() <= 1.0:
                        identities = identities * 100
                    metrics['mean_identity'] = identities.mean()
                    metrics['median_identity'] = identities.median()
                else:
                    metrics['mean_identity'] = None
                    metrics['median_identity'] = None
            else:
                metrics['mean_identity'] = None
                metrics['median_identity'] = None
        except Exception as e:
            logger.debug(f"Could not read identity from diagnostics: {e}")
            metrics['mean_identity'] = None
            metrics['median_identity'] = None
    else:
        # Fallback: try annotated file
        if 'identity' in df.columns:
            assigned_ann = df[df['identity'].notna()]
            if len(assigned_ann) > 0:
                metrics['mean_identity'] = assigned_ann['identity'].mean()
                metrics['median_identity'] = assigned_ann['identity'].median()
            else:
                metrics['mean_identity'] = None
                metrics['median_identity'] = None
        else:
            metrics['mean_identity'] = None
            metrics['median_identity'] = None

    # Assignment metrics
    if 'consensus_group' in df.columns:
        assigned = df[df['consensus_group'].notna()]
        metrics['pct_assigned'] = (len(assigned) / len(df) * 100) if len(df) > 0 else 0

    # Quality flags
    if 'is_tie' in df.columns:
        ties = df['is_tie'].sum() if df['is_tie'].dtype == bool else (df['is_tie'] == True).sum()
        metrics['pct_tie'] = (ties / len(df) * 100) if len(df) > 0 else 0

    if 'potential_misidentification' in df.columns:
        low_conf = df['potential_misidentification'].sum() if df['potential_misidentification'].dtype == bool else (df['potential_misidentification'] == True).sum()
        metrics['pct_low_confidence'] = (low_conf / len(df) * 100) if len(df) > 0 else 0

    # Species metrics
    if 'species' in df.columns:
        metrics['n_species_detected'] = df['species'].nunique()

    # Mixed groups
    if 'group_species_count' in df.columns:
        mixed = df.groupby('consensus_group')['group_species_count'].first()
        metrics['n_mixed_groups'] = (mixed > 1).sum() if len(mixed) > 0 else 0

    # Group purity
    if 'group_majority_pct' in df.columns:
        purities = df.groupby('consensus_group')['group_majority_pct'].first()
        metrics['avg_group_purity'] = purities.mean() if len(purities) > 0 else 0

    return metrics


def cleanup_intermediate_files(output_dir: Path) -> None:
    """
    Remove intermediate files to save space.

    Parameters
    ----------
    output_dir : Path
        Output directory to clean
    """
    # Keep only essential files
    keep_patterns = [
        '*_annotated.csv',
        '*_consensus.fasta',
        '*_diagnostics.csv',
        '*_consensus_taxonomy.csv'
    ]

    # Remove intermediate directory
    intermediate_dir = output_dir / 'intermediate'
    if intermediate_dir.exists():
        shutil.rmtree(intermediate_dir)

    # Remove other directories
    for subdir in ['visualization', 'phylogenetic', 'reports', 'quality_control', 'divergence_analysis']:
        dir_path = output_dir / subdir
        if dir_path.exists():
            shutil.rmtree(dir_path)


def track_group_membership(
    run_outputs: Dict[float, Path],
    organism: str,
    thresholds: List[float]
) -> pd.DataFrame:
    """
    Track how samples move between groups as threshold varies.

    Uses a biologically meaningful approach: tracks whether samples cluster
    WITH THE SAME PARTNERS, regardless of arbitrary group names.

    Parameters
    ----------
    run_outputs : dict
        Mapping of threshold to output directory
    organism : str
        Organism name
    thresholds : list
        List of thresholds tested

    Returns
    -------
    pd.DataFrame
        Sample-level tracking table with cluster composition changes
    """
    logger.info("  Tracking group membership across thresholds...")

    # Load annotated data for each threshold
    sample_assignments = {}  # threshold -> {sample_id: group_name}
    group_members = {}  # threshold -> {group_name: set(sample_ids)}

    for threshold in sorted(thresholds):
        output_dir = run_outputs[threshold]
        annotated_path = output_dir / f"{organism}_annotated.csv"

        if annotated_path.exists():
            df = pd.read_csv(annotated_path)
            if 'processid' in df.columns and 'consensus_group' in df.columns:
                # Store sample -> group mapping
                assignments = dict(zip(df['processid'], df['consensus_group']))
                sample_assignments[threshold] = assignments

                # Store group -> members mapping
                members = {}
                for sample, group in assignments.items():
                    if pd.notna(group):
                        if group not in members:
                            members[group] = set()
                        members[group].add(sample)
                group_members[threshold] = members

    if not sample_assignments:
        return pd.DataFrame()

    # Get all sample IDs
    all_samples = set()
    for assignments in sample_assignments.values():
        all_samples.update(assignments.keys())

    # Build tracking table
    tracking_data = []

    for sample_id in sorted(all_samples):
        row = {'processid': sample_id}

        # Get assignments at each threshold and track cluster partners
        cluster_partners_by_threshold = []

        for threshold in sorted(thresholds):
            if threshold in sample_assignments:
                group = sample_assignments[threshold].get(sample_id, 'unassigned')
                row[f't_{threshold}'] = group

                # Get the set of samples this sample clusters with
                if pd.notna(group) and group != 'unassigned' and threshold in group_members:
                    partners = group_members[threshold].get(group, set()) - {sample_id}
                    cluster_partners_by_threshold.append((threshold, partners))
                else:
                    cluster_partners_by_threshold.append((threshold, set()))
            else:
                row[f't_{threshold}'] = 'NA'

        # Count biologically meaningful changes:
        # How many times does the cluster composition change significantly?
        n_composition_changes = 0

        for i in range(1, len(cluster_partners_by_threshold)):
            prev_threshold, prev_partners = cluster_partners_by_threshold[i-1]
            curr_threshold, curr_partners = cluster_partners_by_threshold[i]

            # Skip if either set is empty (unassigned)
            if not prev_partners and not curr_partners:
                continue

            # Calculate Jaccard similarity between partner sets
            if prev_partners or curr_partners:
                intersection = len(prev_partners & curr_partners)
                union = len(prev_partners | curr_partners)
                jaccard = intersection / union if union > 0 else 0

                # Consider it a change if Jaccard similarity < 0.7
                # (more than 30% of partners changed)
                if jaccard < 0.7:
                    n_composition_changes += 1

        row['n_changes'] = n_composition_changes

        # Add cluster size information
        cluster_sizes = []
        for threshold, partners in cluster_partners_by_threshold:
            cluster_sizes.append(len(partners) + 1)  # +1 for the sample itself

        if cluster_sizes:
            row['mean_cluster_size'] = np.mean(cluster_sizes)
            row['min_cluster_size'] = min(cluster_sizes)
            row['max_cluster_size'] = max(cluster_sizes)

        # Stability score based on composition changes
        if n_composition_changes == 0:
            stability = 'high'
        elif n_composition_changes <= 2:
            stability = 'medium'
        else:
            stability = 'low'
        row['stability_score'] = stability

        tracking_data.append(row)

    return pd.DataFrame(tracking_data)


def create_threshold_stability_plot(
    summary_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Create multi-panel stability plot.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics for all thresholds
    output_path : Path
        Output path for PDF
    """
    logger.info("  Creating threshold stability plot...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Parameter Sweep Stability Analysis', fontsize=16, fontweight='bold')

    # Panel 1: Number of groups
    ax1 = axes[0, 0]
    ax1.plot(summary_df['threshold'], summary_df['n_groups'], 'o-', linewidth=2, markersize=8)
    ax1.set_xlabel('Clustering Threshold', fontweight='bold')
    ax1.set_ylabel('Number of Groups', fontweight='bold')
    ax1.set_title('Grouping Stability')
    ax1.grid(alpha=0.3)

    # Panel 2: Mean identity (if available)
    ax2 = axes[0, 1]
    if 'mean_identity' in summary_df.columns and summary_df['mean_identity'].notna().any():
        ax2.plot(summary_df['threshold'], summary_df['mean_identity'], 'o-',
                 linewidth=2, markersize=8, color='green')
        ax2.set_xlabel('Clustering Threshold', fontweight='bold')
        ax2.set_ylabel('Mean Identity (%)', fontweight='bold')
        ax2.set_title('Assignment Quality')
        ax2.grid(alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'Identity data not available',
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Assignment Quality (No Data)')
        ax2.axis('off')

    # Panel 3: Assignment rate
    ax3 = axes[1, 0]
    ax3.plot(summary_df['threshold'], summary_df['pct_assigned'], 'o-',
             linewidth=2, markersize=8, color='orange')
    ax3.set_xlabel('Clustering Threshold', fontweight='bold')
    ax3.set_ylabel('% Assigned', fontweight='bold')
    ax3.set_title('Assignment Coverage')
    ax3.grid(alpha=0.3)

    # Panel 4: Singletons
    ax4 = axes[1, 1]
    ax4.plot(summary_df['threshold'], summary_df['n_singletons'], 'o-',
             linewidth=2, markersize=8, color='red')
    ax4.set_xlabel('Clustering Threshold', fontweight='bold')
    ax4.set_ylabel('Number of Singletons', fontweight='bold')
    ax4.set_title('Over-splitting Indicator')
    ax4.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  ✓ Stability plot saved to {output_path}")


def create_elbow_plot(
    summary_df: pd.DataFrame,
    output_path: Path,
    recommended_threshold: Optional[float] = None
) -> None:
    """
    Create elbow plot to identify optimal threshold.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics
    output_path : Path
        Output path for PDF
    recommended_threshold : float, optional
        Recommended threshold to highlight
    """
    logger.info("  Creating elbow plot...")

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Primary axis: Number of groups
    color = 'tab:blue'
    ax1.set_xlabel('Clustering Threshold', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Groups', fontsize=12, fontweight='bold', color=color)
    ax1.plot(summary_df['threshold'], summary_df['n_groups'], 'o-',
             linewidth=2, markersize=10, color=color, label='Number of Groups')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(alpha=0.3)

    # Secondary axis: Mean identity (if available)
    if 'mean_identity' in summary_df.columns and summary_df['mean_identity'].notna().any():
        ax2 = ax1.twinx()
        color = 'tab:red'
        ax2.set_ylabel('Mean Identity (%)', fontsize=12, fontweight='bold', color=color)
        ax2.plot(summary_df['threshold'], summary_df['mean_identity'], 's--',
                 linewidth=2, markersize=8, color=color, label='Mean Identity', alpha=0.7)
        ax2.tick_params(axis='y', labelcolor=color)

    # Highlight recommended threshold
    if recommended_threshold is not None:
        ax1.axvline(recommended_threshold, color='green', linestyle='--',
                    linewidth=2, label=f'Recommended: {recommended_threshold}')

    # Title and legend
    ax1.set_title('Elbow Plot: Optimal Threshold Detection',
                  fontsize=14, fontweight='bold', pad=20)

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  ✓ Elbow plot saved to {output_path}")


def find_elbow_point(summary_df: pd.DataFrame) -> float:
    """
    Find the elbow point in the number of groups curve.

    Uses the "knee" detection algorithm to find where grouping stabilizes.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics

    Returns
    -------
    float
        Recommended threshold at elbow point
    """
    thresholds = summary_df['threshold'].values
    n_groups = summary_df['n_groups'].values

    # Simple elbow detection: find maximum curvature
    # Calculate second derivative (rate of change of slope)
    if len(thresholds) < 3:
        return thresholds[len(thresholds) // 2]

    # Calculate first derivative (slope)
    slopes = np.diff(n_groups) / np.diff(thresholds)

    # Calculate second derivative (curvature)
    curvatures = np.diff(slopes) / np.diff(thresholds[:-1])

    # Find maximum absolute curvature
    elbow_idx = np.argmax(np.abs(curvatures)) + 1

    return thresholds[elbow_idx]


def generate_recommendations(
    summary_df: pd.DataFrame,
    organism: str,
    output_path: Path
) -> float:
    """
    Generate automated recommendations text file.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics
    organism : str
        Organism name
    output_path : Path
        Output path for recommendations file

    Returns
    -------
    float
        Recommended threshold
    """
    logger.info("  Generating recommendations...")

    # Find elbow point
    recommended_threshold = find_elbow_point(summary_df)

    # Get metrics for recommended threshold
    rec_row = summary_df[summary_df['threshold'] == recommended_threshold].iloc[0]

    # Define acceptable range (thresholds with similar grouping)
    rec_n_groups = rec_row['n_groups']

    # Build filter conditions
    conditions = [
        (summary_df['n_groups'] >= rec_n_groups * 0.8),
        (summary_df['n_groups'] <= rec_n_groups * 1.2)
    ]

    # Add adaptive identity filter only if available
    if 'mean_identity' in summary_df.columns and summary_df['mean_identity'].notna().any():
        # Use an adaptive threshold based on the data itself
        # Instead of hardcoded 95%, use the 75th percentile of mean identities
        identity_threshold = summary_df['mean_identity'].quantile(0.75)
        # But don't go below 90% as a safety measure for barcode data
        identity_threshold = max(identity_threshold, 90.0)
        conditions.append(summary_df['mean_identity'] >= identity_threshold)
        logger.debug(f"Using adaptive identity threshold: {identity_threshold:.1f}%")
    else:
        identity_threshold = None

    # Combine conditions
    acceptable = summary_df[pd.concat(conditions, axis=1).all(axis=1)]

    # If no thresholds meet criteria, relax constraints
    if len(acceptable) == 0:
        logger.warning("No thresholds meet all criteria; relaxing constraints")
        # Just use grouping constraint
        acceptable = summary_df[
            (summary_df['n_groups'] >= rec_n_groups * 0.8) &
            (summary_df['n_groups'] <= rec_n_groups * 1.2)
        ]

    if len(acceptable) > 0:
        min_threshold = acceptable['threshold'].min()
        max_threshold = acceptable['threshold'].max()
    else:
        # Fallback: use recommended threshold ±0.01
        min_threshold = max(summary_df['threshold'].min(), recommended_threshold - 0.01)
        max_threshold = min(summary_df['threshold'].max(), recommended_threshold + 0.01)

    # Generate recommendations text
    text = f"""Parameter Sweep Analysis Results
=================================

Dataset: {organism}
Total Samples: {rec_row['n_samples']}
Threshold Range Tested: {summary_df['threshold'].min():.3f} - {summary_df['threshold'].max():.3f}
Number of Thresholds Tested: {len(summary_df)}

SUMMARY:
--------
Tested {len(summary_df)} threshold values across {summary_df['threshold'].min():.3f} to {summary_df['threshold'].max():.3f} range.

KEY FINDINGS:
-------------
1. Grouping stability achieved at threshold >= {recommended_threshold:.3f}
   - At {summary_df['threshold'].min():.3f}: {summary_df.iloc[0]['n_groups']} groups
   - At {recommended_threshold:.3f}: {rec_row['n_groups']} groups (elbow point)
   - At {summary_df['threshold'].max():.3f}: {summary_df.iloc[-1]['n_groups']} groups

2. Assignment quality optimal between {min_threshold:.3f} - {max_threshold:.3f}
   - Assignment rate: {rec_row['pct_assigned']:.1f}%"""

    # Add mean identity only if available
    if 'mean_identity' in rec_row.index and pd.notna(rec_row['mean_identity']):
        text += f"""
   - Mean identity: {rec_row['mean_identity']:.1f}%"""

    text += f"""

3. Singletons and over-splitting:
   - At {summary_df['threshold'].min():.3f}: {summary_df.iloc[0]['n_singletons']} singletons
   - At {recommended_threshold:.3f}: {rec_row['n_singletons']} singletons
   - At {summary_df['threshold'].max():.3f}: {summary_df.iloc[-1]['n_singletons']} singletons

RECOMMENDATIONS:
----------------
PRIMARY: threshold = {recommended_threshold:.3f}
  Rationale:
  - Elbow point in grouping curve"""

    # Add mean identity rationale only if available
    if 'mean_identity' in rec_row.index and pd.notna(rec_row['mean_identity']):
        text += f"""
  - High mean identity ({rec_row['mean_identity']:.1f}%)"""

    text += f"""
  - Balanced singletons ({rec_row['n_singletons']} groups)
  - Assignment rate: {rec_row['pct_assigned']:.1f}%

ACCEPTABLE RANGE: {min_threshold:.3f} - {max_threshold:.3f}
  Rationale:
  - Similar grouping structure (±20%)"""

    # Add adaptive identity threshold information
    if identity_threshold is not None:
        text += f"""
  - Maintains identity ≥{identity_threshold:.1f}% (data-driven threshold)"""

    text += """
  - Similar assignment rates

IMPORTANT NOTES:
----------------
- These recommendations are based on automated elbow detection
- The identity threshold is adaptively calculated from your data
- Always visually inspect the stability plots to confirm recommendations
- Consider biological knowledge and study goals when finalizing threshold
- For different markers or taxonomic levels, these values may need adjustment

AVOID:
  - < {min_threshold:.3f}: May over-split (too many groups)
  - > {max_threshold:.3f}: May under-split (lower identity)

NEXT STEPS:
-----------
1. Review the threshold_stability.pdf plot to visually confirm the elbow point

2. Run full analysis with recommended threshold:
   boldgenotyper {organism}.tsv --clustering-threshold {recommended_threshold:.3f}

3. If results seem over/under-split, adjust within acceptable range

4. Examine group_membership_tracking.csv to see how samples cluster together
   - n_changes: Number of times cluster composition changed significantly
   - stability_score: Overall stability (high/medium/low)
   - Look for samples with high n_changes - these may be problematic

5. Document threshold choice in methods with justification

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(output_path, 'w') as f:
        f.write(text)

    logger.info(f"  ✓ Recommendations saved to {output_path}")

    return recommended_threshold


def run_parameter_sweep(
    tsv_path: Path,
    thresholds: List[float],
    output_dir: Path,
    organism: Optional[str] = None,
    threads: int = 4,
    keep_intermediates: bool = False
) -> Dict[str, Any]:
    """
    Complete parameter sweep workflow.

    Parameters
    ----------
    tsv_path : Path
        Input TSV file
    thresholds : list
        List of thresholds to test
    output_dir : Path
        Output directory
    organism : str, optional
        Organism name
    threads : int
        Number of threads
    keep_intermediates : bool
        Keep full output for each run

    Returns
    -------
    dict
        Results and file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Infer organism if not provided
    if organism is None:
        from .cli import extract_organism_from_path
        organism = extract_organism_from_path(tsv_path)

    logger.info("="*70)
    logger.info(f"BOLDGenotyper Parameter Sweep: {organism}")
    logger.info("="*70)
    logger.info(f"Testing {len(thresholds)} thresholds: {thresholds}")
    logger.info(f"Output: {output_dir}")
    logger.info("="*70)

    results = {'organism': organism, 'thresholds': thresholds}

    # Create runs directory
    runs_dir = output_dir / 'runs'
    runs_dir.mkdir(exist_ok=True)

    # Run each threshold
    run_outputs = {}
    summary_data = []

    for threshold in thresholds:
        threshold_dir = runs_dir / f"threshold_{threshold:.3f}".replace('.', '_')

        metrics = run_single_threshold(
            tsv_path=tsv_path,
            organism=organism,
            threshold=threshold,
            output_dir=threshold_dir,
            threads=threads,
            keep_intermediates=keep_intermediates
        )

        if metrics:
            summary_data.append(metrics)
            run_outputs[threshold] = threshold_dir

    if not summary_data:
        raise ParameterSweepError("No successful runs completed")

    # Create summary table
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('threshold')

    summary_path = output_dir / "sweep_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    results['summary'] = summary_path
    logger.info(f"✓ Summary table saved: {summary_path}")

    # Track group membership
    tracking_df = track_group_membership(run_outputs, organism, thresholds)
    if len(tracking_df) > 0:
        tracking_path = output_dir / "group_membership_tracking.csv"
        tracking_df.to_csv(tracking_path, index=False)
        results['tracking'] = tracking_path
        logger.info(f"✓ Group membership tracking saved: {tracking_path}")

    # Create visualizations
    try:
        stability_path = output_dir / "threshold_stability.pdf"
        create_threshold_stability_plot(summary_df, stability_path)
        results['stability_plot'] = stability_path
    except Exception as e:
        logger.warning(f"Could not create stability plot: {e}")

    # Generate recommendations
    try:
        recommendations_path = output_dir / "recommendations.txt"
        recommended_threshold = generate_recommendations(
            summary_df, organism, recommendations_path
        )
        results['recommendations'] = recommendations_path
        results['recommended_threshold'] = recommended_threshold

        # Create elbow plot
        elbow_path = output_dir / "elbow_plot.pdf"
        create_elbow_plot(summary_df, elbow_path, recommended_threshold)
        results['elbow_plot'] = elbow_path
    except Exception as e:
        logger.warning(f"Could not generate recommendations: {e}")

    # Create README
    readme_path = output_dir / "README.md"
    readme_content = f"""# Parameter Sweep Results for {organism}

## Summary

Tested {len(thresholds)} clustering thresholds to identify optimal parameter values.

**Recommended Threshold:** {results.get('recommended_threshold', 'See recommendations.txt')}

## Files

- `sweep_summary.csv` - Metrics for all tested thresholds
- `threshold_stability.pdf` - Multi-panel stability visualization
- `elbow_plot.pdf` - Optimal threshold detection
- `group_membership_tracking.csv` - Sample clustering stability across thresholds
- `recommendations.txt` - Detailed recommendations and rationale
- `runs/` - Individual outputs for each threshold

## Key Files Explained

### group_membership_tracking.csv

This file tracks how samples cluster together across different thresholds using a **biologically meaningful approach**:

- **Columns t_X.XXX**: Show which consensus group a sample was assigned to at each threshold
- **n_changes**: Number of times the sample's cluster COMPOSITION changed significantly
  - Based on Jaccard similarity of cluster partners (not arbitrary group names)
  - A change is counted when >30% of cluster partners differ between consecutive thresholds
  - Low values (0-2) indicate stable clustering
- **stability_score**: Overall assessment (high/medium/low)
  - high: n_changes = 0 (sample always clusters with same partners)
  - medium: n_changes = 1-2 (minor composition changes)
  - low: n_changes > 2 (unstable clustering - potential problem samples)
- **mean_cluster_size**: Average number of samples in this sample's cluster across thresholds
- **min/max_cluster_size**: Size range of clusters this sample belonged to

**Important:** The group NAMES (e.g., c15, c18) are arbitrary and may change between runs.
What matters is which OTHER samples cluster together, not what the group is called.

## Interpretation

See `recommendations.txt` for detailed analysis and next steps.

The recommended threshold represents the "elbow point" where:
1. Grouping structure stabilizes
2. Identity remains high (data-driven adaptive threshold)
3. Singletons are minimized
4. Assignment coverage is maximized

**Note:** Identity thresholds are calculated adaptively from your data, not hardcoded.
Different markers or taxonomic levels may have different optimal values.

## Next Steps

1. Review `threshold_stability.pdf` to visually confirm the elbow point

2. Examine `group_membership_tracking.csv`:
   - Sort by n_changes descending to find unstable samples
   - Investigate samples with stability_score = 'low'

3. Run full analysis with recommended threshold:
   ```bash
   boldgenotyper {tsv_path.name} --clustering-threshold {results.get('recommended_threshold', 0.015):.3f}
   ```

4. Document your threshold choice and rationale in your methods

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info("="*70)
    logger.info("✓ Parameter sweep complete!")

    # Format recommended threshold properly
    rec_threshold = results.get('recommended_threshold', None)
    if rec_threshold is not None and isinstance(rec_threshold, (int, float)):
        logger.info(f"  Recommended threshold: {rec_threshold:.3f}")
    else:
        logger.info(f"  Recommended threshold: N/A")

    logger.info(f"  Results: {output_dir}")
    logger.info("="*70)

    return results
