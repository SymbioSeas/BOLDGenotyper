"""
Parameter Sweep Module for BOLDGenotyper

Sweeps the singleton error-filtering threshold (min_singleton_distance) across
a range of values to help users choose the setting that best fits their dataset.

The pipeline discovers haplotypes as Exact Sequence Variants (ESVs) — unique
consensus sequences.  Singletons (one-member haplotypes) that are very close
to an existing multi-member haplotype are likely sequencing or PCR errors and
are removed when their nearest-neighbour distance is ≤ min_singleton_distance.

This module runs the pipeline once per threshold value, collects haplotype-level
and sample-level metrics, and produces:
- sweep_summary.csv           per-threshold metric table
- group_membership_tracking.csv  sample-level stability across thresholds
- threshold_stability.pdf     4-panel stability visualisation
- elbow_plot.pdf              optimal-threshold detection
- recommendations.txt         automated interpretation

Author: Steph Smith (symbioseas@outlook.com)
"""

from __future__ import annotations
from typing import Dict, List, Optional, Any
from pathlib import Path
import logging
from datetime import datetime
import shutil

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
    Run BOLDGenotyper pipeline with a single min_singleton_distance value.

    Parameters
    ----------
    tsv_path : Path
        Input TSV file
    organism : str
        Organism name
    threshold : float
        Singleton distance threshold to test.  Singletons with nearest-
        neighbour distance ≤ this value are removed as likely errors.
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
    # max_singleton_distance must be ≥ min_singleton_distance; raise it when
    # the sweep threshold exceeds the default (0.05).
    max_singleton_dist = max(threshold, cfg.haplotype.max_singleton_distance)
    cfg = cfg.update(
        haplotype__min_singleton_distance=threshold,
        haplotype__max_singleton_distance=max_singleton_dist,
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

    Reads two files:
    - haplotypes/{organism}_haplotype_stats.csv  → haplotype counts
    - {organism}_annotated.csv                   → sample-level metrics

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

    # ---------------------------------------------------------------------------
    # Haplotype-level metrics (from haplotype_stats.csv)
    # ---------------------------------------------------------------------------
    stats_path = output_dir / "haplotypes" / f"{organism}_haplotype_stats.csv"
    if stats_path.exists():
        stats_df = pd.read_csv(stats_path)
        metrics['n_haplotypes'] = len(stats_df)
        if 'is_singleton' in stats_df.columns:
            metrics['n_singletons'] = int(stats_df['is_singleton'].sum())
        if 'is_suspect' in stats_df.columns:
            metrics['n_suspect'] = int(stats_df['is_suspect'].sum())

    # ---------------------------------------------------------------------------
    # Sample-level metrics (from annotated CSV)
    # ---------------------------------------------------------------------------
    annotated_path = output_dir / f"{organism}_annotated.csv"
    if not annotated_path.exists():
        return metrics

    df = pd.read_csv(annotated_path)
    metrics['n_samples'] = len(df)

    # Assignment coverage: fraction of samples assigned to a haplotype
    if 'haplotype_id' in df.columns:
        assigned = df[df['haplotype_id'].notna()]
        metrics['pct_assigned'] = (len(assigned) / len(df) * 100) if len(df) > 0 else 0

    # Species richness (from depositor taxonomy — constant across thresholds,
    # included for completeness in the summary table)
    if 'species' in df.columns:
        metrics['n_species_detected'] = df['species'].nunique()

    # Potential misidentification flags — fraction of samples whose assigned
    # haplotype does not match their depositor species.  Can change when
    # singleton removal causes samples to be reassigned.
    if 'potential_misidentification' in df.columns:
        misid = (
            df['potential_misidentification'].sum()
            if df['potential_misidentification'].dtype == bool
            else (df['potential_misidentification'] == True).sum()
        )
        metrics['pct_potential_misid'] = (misid / len(df) * 100) if len(df) > 0 else 0

    # Mixed haplotype groups (haplotypes containing >1 depositor species)
    if 'group_species_count' in df.columns and 'haplotype_id' in df.columns:
        mixed = df.groupby('haplotype_id')['group_species_count'].first()
        metrics['n_mixed_groups'] = int((mixed > 1).sum()) if len(mixed) > 0 else 0

    # Average group purity
    if 'group_majority_pct' in df.columns and 'haplotype_id' in df.columns:
        purities = df.groupby('haplotype_id')['group_majority_pct'].first()
        metrics['avg_group_purity'] = purities.mean() if len(purities) > 0 else 0

    return metrics


def cleanup_intermediate_files(output_dir: Path) -> None:
    """
    Remove per-threshold output directories that are no longer needed after
    metrics have been extracted.  Retains only the annotated CSV and haplotype
    stats that the analysis scripts read.

    Parameters
    ----------
    output_dir : Path
        Output directory to clean
    """
    # These directories are produced by the sweep runs but are not needed
    # once extract_run_metrics has pulled the numbers.
    removable = [
        'haplotype_assignments',
        'species_analysis',
        'taxonomy',
        # Legacy names kept for safety if an older pipeline version wrote them
        'intermediate',
        'visualization',
        'phylogenetic',
        'reports',
        'quality_control',
        'divergence_analysis',
    ]
    for name in removable:
        dir_path = output_dir / name
        if dir_path.exists():
            shutil.rmtree(dir_path)


def track_group_membership(
    run_outputs: Dict[float, Path],
    organism: str,
    thresholds: List[float]
) -> pd.DataFrame:
    """
    Track how samples move between haplotypes as the singleton distance
    threshold varies.

    Uses a biologically meaningful approach: tracks whether samples cluster
    WITH THE SAME PARTNERS, regardless of arbitrary haplotype names.

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
    sample_assignments = {}  # threshold -> {sample_id: haplotype_id}
    group_members = {}       # threshold -> {haplotype_id: set(sample_ids)}

    for threshold in sorted(thresholds):
        output_dir = run_outputs[threshold]
        annotated_path = output_dir / f"{organism}_annotated.csv"

        if annotated_path.exists():
            df = pd.read_csv(annotated_path)
            if 'processid' in df.columns and 'haplotype_id' in df.columns:
                assignments = dict(zip(df['processid'], df['haplotype_id']))
                sample_assignments[threshold] = assignments

                members = {}
                for sample, hap in assignments.items():
                    if pd.notna(hap):
                        members.setdefault(hap, set()).add(sample)
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

        cluster_partners_by_threshold = []

        for threshold in sorted(thresholds):
            if threshold in sample_assignments:
                hap = sample_assignments[threshold].get(sample_id, 'unassigned')
                row[f't_{threshold}'] = hap

                if pd.notna(hap) and hap != 'unassigned' and threshold in group_members:
                    partners = group_members[threshold].get(hap, set()) - {sample_id}
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

            if not prev_partners and not curr_partners:
                continue

            if prev_partners or curr_partners:
                intersection = len(prev_partners & curr_partners)
                union = len(prev_partners | curr_partners)
                jaccard = intersection / union if union > 0 else 0

                # A change is counted when >30% of partners differ
                if jaccard < 0.7:
                    n_composition_changes += 1

        row['n_changes'] = n_composition_changes

        # Cluster size information
        cluster_sizes = [len(partners) + 1 for _, partners in cluster_partners_by_threshold]
        if cluster_sizes:
            row['mean_cluster_size'] = np.mean(cluster_sizes)
            row['min_cluster_size'] = min(cluster_sizes)
            row['max_cluster_size'] = max(cluster_sizes)

        # Stability score
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
    Create 4-panel stability plot showing key metrics vs threshold.

    Panels:
      A. Total haplotypes (primary signal of the sweep)
      B. Singletons retained (expected to drop with rising threshold)
      C. Suspect haplotypes (distant singletons flagged as contamination)
      D. Assignment coverage and misidentification rate

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics for all thresholds
    output_path : Path
        Output path for PDF
    """
    logger.info("  Creating threshold stability plot...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Singleton Distance Threshold — Stability Analysis',
                 fontsize=16, fontweight='bold')

    # Panel A: Number of haplotypes
    ax = axes[0, 0]
    ax.plot(summary_df['threshold'], summary_df['n_haplotypes'],
            'o-', linewidth=2, markersize=8, color='tab:blue')
    ax.set_xlabel('min_singleton_distance')
    ax.set_ylabel('Number of Haplotypes')
    ax.set_title('A. Haplotype Count')
    ax.grid(alpha=0.3)

    # Panel B: Singletons retained
    ax = axes[0, 1]
    if 'n_singletons' in summary_df.columns:
        ax.plot(summary_df['threshold'], summary_df['n_singletons'],
                'o-', linewidth=2, markersize=8, color='coral')
    ax.set_xlabel('min_singleton_distance')
    ax.set_ylabel('Singletons Retained')
    ax.set_title('B. Singleton Count')
    ax.grid(alpha=0.3)

    # Panel C: Suspect haplotypes
    ax = axes[1, 0]
    if 'n_suspect' in summary_df.columns:
        ax.plot(summary_df['threshold'], summary_df['n_suspect'],
                'o-', linewidth=2, markersize=8, color='tab:orange')
    ax.set_xlabel('min_singleton_distance')
    ax.set_ylabel('Suspect Haplotypes')
    ax.set_title('C. Suspect Count')
    ax.grid(alpha=0.3)

    # Panel D: Assignment coverage + misidentification rate (dual y-axis)
    ax = axes[1, 1]
    if 'pct_assigned' in summary_df.columns:
        ax.plot(summary_df['threshold'], summary_df['pct_assigned'],
                'o-', linewidth=2, markersize=8, color='tab:green', label='% Assigned')
        ax.set_ylabel('% Assigned', color='tab:green')
        ax.tick_params(axis='y', labelcolor='tab:green')

    if 'pct_potential_misid' in summary_df.columns:
        ax2 = ax.twinx()
        ax2.plot(summary_df['threshold'], summary_df['pct_potential_misid'],
                 's--', linewidth=2, markersize=6, color='tab:red', alpha=0.7,
                 label='% Misid')
        ax2.set_ylabel('% Potential Misidentification', color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')

    ax.set_xlabel('min_singleton_distance')
    ax.set_title('D. Assignment Quality')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  Stability plot saved to {output_path}")


def create_elbow_plot(
    summary_df: pd.DataFrame,
    output_path: Path,
    recommended_threshold: Optional[float] = None
) -> None:
    """
    Create elbow plot showing haplotype count vs threshold with the
    recommended threshold highlighted.  Secondary axis shows singleton count.

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

    # Primary axis: haplotype count
    color1 = 'tab:blue'
    ax1.set_xlabel('min_singleton_distance', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Haplotypes', fontsize=12, fontweight='bold', color=color1)
    ax1.plot(summary_df['threshold'], summary_df['n_haplotypes'],
             'o-', linewidth=2, markersize=10, color=color1, label='Haplotypes')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(alpha=0.3)

    # Secondary axis: singleton count
    if 'n_singletons' in summary_df.columns and summary_df['n_singletons'].notna().any():
        ax2 = ax1.twinx()
        color2 = 'tab:red'
        ax2.set_ylabel('Singletons Retained', fontsize=12, fontweight='bold', color=color2)
        ax2.plot(summary_df['threshold'], summary_df['n_singletons'],
                 's--', linewidth=2, markersize=8, color=color2, label='Singletons', alpha=0.7)
        ax2.tick_params(axis='y', labelcolor=color2)

    # Highlight recommended threshold
    if recommended_threshold is not None:
        ax1.axvline(recommended_threshold, color='green', linestyle='--',
                    linewidth=2, label=f'Recommended: {recommended_threshold}')

    ax1.set_title('Elbow Plot: Optimal Singleton Distance Threshold',
                  fontsize=14, fontweight='bold', pad=20)

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = (ax2.get_legend_handles_labels()
                        if 'ax2' in dir() else ([], []))
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  Elbow plot saved to {output_path}")


def find_elbow_point(summary_df: pd.DataFrame) -> float:
    """
    Find the elbow point in the haplotype-count curve — the threshold
    beyond which further increases produce diminishing reductions.

    Uses maximum-curvature detection on the second derivative.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics (must contain 'threshold' and 'n_haplotypes')

    Returns
    -------
    float
        Recommended threshold at elbow point
    """
    thresholds = summary_df['threshold'].values
    n_haplotypes = summary_df['n_haplotypes'].values

    if len(thresholds) < 3:
        return thresholds[len(thresholds) // 2]

    # First derivative (slope)
    slopes = np.diff(n_haplotypes) / np.diff(thresholds)

    # Second derivative (curvature)
    curvatures = np.diff(slopes) / np.diff(thresholds[:-1])

    # Maximum absolute curvature = elbow
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

    # Define acceptable range: thresholds with similar haplotype count (±20%)
    rec_n_haplotypes = rec_row['n_haplotypes']
    acceptable = summary_df[
        (summary_df['n_haplotypes'] >= rec_n_haplotypes * 0.8) &
        (summary_df['n_haplotypes'] <= rec_n_haplotypes * 1.2)
    ]

    if len(acceptable) > 0:
        min_threshold = acceptable['threshold'].min()
        max_threshold = acceptable['threshold'].max()
    else:
        min_threshold = max(summary_df['threshold'].min(), recommended_threshold - 0.01)
        max_threshold = min(summary_df['threshold'].max(), recommended_threshold + 0.01)

    # Singleton counts at boundaries
    first_row = summary_df.iloc[0]
    last_row  = summary_df.iloc[-1]

    text = f"""Singleton Distance Threshold — Sweep Results
=============================================

Dataset:                  {organism}
Total Samples:            {rec_row['n_samples']}
Threshold Range Tested:   {summary_df['threshold'].min():.3f} – {summary_df['threshold'].max():.3f}
Number of Thresholds:     {len(summary_df)}

WHAT WAS SWEPT
--------------
min_singleton_distance controls removal of singleton haplotypes (one-member
ESVs) that are very similar to an existing multi-member haplotype.  These
singletons are likely sequencing or PCR errors.  A higher threshold removes
more singletons; a lower threshold retains them.

KEY FINDINGS
------------
1. Haplotype count
   At {summary_df['threshold'].min():.3f}:               {first_row['n_haplotypes']} haplotypes
   At {recommended_threshold:.3f} (elbow):         {rec_row['n_haplotypes']} haplotypes
   At {summary_df['threshold'].max():.3f}:               {last_row['n_haplotypes']} haplotypes

2. Singletons retained
   At {summary_df['threshold'].min():.3f}:               {first_row.get('n_singletons', 'N/A')}
   At {recommended_threshold:.3f} (elbow):         {rec_row.get('n_singletons', 'N/A')}
   At {summary_df['threshold'].max():.3f}:               {last_row.get('n_singletons', 'N/A')}

3. Assignment coverage at elbow:  {rec_row.get('pct_assigned', 'N/A'):.1f}%"""

    if 'pct_potential_misid' in rec_row.index and pd.notna(rec_row.get('pct_potential_misid')):
        text += f"""
   Potential misidentification rate:  {rec_row['pct_potential_misid']:.1f}%"""

    text += f"""

RECOMMENDATIONS
---------------
PRIMARY:  min_singleton_distance = {recommended_threshold:.3f}
  • Elbow point in the haplotype-count curve
  • {rec_row['n_haplotypes']} haplotypes, {rec_row.get('n_singletons', 'N/A')} singletons retained
  • Assignment coverage: {rec_row.get('pct_assigned', 'N/A'):.1f}%

ACCEPTABLE RANGE:  {min_threshold:.3f} – {max_threshold:.3f}
  • Haplotype count within ±20% of the elbow value
  • Use lower end to retain more rare variants
  • Use upper end for a more conservative (error-filtered) set

AVOID
-----
  < {min_threshold:.3f}  — retains singletons that are likely sequencing errors
  > {max_threshold:.3f}  — removes biologically real rare haplotypes

IMPORTANT NOTES
---------------
• These recommendations are based on automated elbow detection.
• Always visually inspect threshold_stability.pdf to confirm.
• Biological knowledge and study goals should inform the final choice.
• group_membership_tracking.csv shows which samples change haplotype
  assignments across thresholds — samples with stability_score = 'low'
  warrant manual review.

NEXT STEPS
----------
1. Review threshold_stability.pdf and elbow_plot.pdf

2. Run the full pipeline with your chosen threshold using the
   --singleton-distance flag:

       boldgenotyper --input your_data.tsv \\
                     --output output/ \\
                     --singleton-distance {recommended_threshold:.3f} \\
                     [other options]

   The default threshold is 0.005 (0.5%).  Always pass the sweep-recommended
   value explicitly so your analysis is reproducible.

3. If the result looks over- or under-split, adjust within the acceptable
   range and re-run.

4. Document the threshold and its rationale in your methods section.

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(output_path, 'w') as f:
        f.write(text)

    logger.info(f"  Recommendations saved to {output_path}")

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

    Runs the pipeline once per threshold value, collects metrics, tracks
    sample-level haplotype membership stability, and generates plots and
    an automated recommendation.

    Parameters
    ----------
    tsv_path : Path
        Input TSV file
    thresholds : list
        min_singleton_distance values to test
    output_dir : Path
        Output directory
    organism : str, optional
        Organism name (inferred from TSV path if omitted)
    threads : int
        Number of threads
    keep_intermediates : bool
        Keep full per-threshold output

    Returns
    -------
    dict
        Results including paths to summary files and recommended threshold
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
    logger.info(f"Sweeping min_singleton_distance across {len(thresholds)} values: {thresholds}")
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
    logger.info(f"Summary table saved: {summary_path}")

    # Track group membership
    tracking_df = track_group_membership(run_outputs, organism, thresholds)
    if len(tracking_df) > 0:
        tracking_path = output_dir / "group_membership_tracking.csv"
        tracking_df.to_csv(tracking_path, index=False)
        results['tracking'] = tracking_path
        logger.info(f"Group membership tracking saved: {tracking_path}")

    # Create visualizations
    try:
        stability_path = output_dir / "threshold_stability.pdf"
        create_threshold_stability_plot(summary_df, stability_path)
        results['stability_plot'] = stability_path
    except Exception as e:
        logger.warning(f"Could not create stability plot: {e}")

    # Generate recommendations and elbow plot
    try:
        recommendations_path = output_dir / "recommendations.txt"
        recommended_threshold = generate_recommendations(
            summary_df, organism, recommendations_path
        )
        results['recommendations'] = recommendations_path
        results['recommended_threshold'] = recommended_threshold

        elbow_path = output_dir / "elbow_plot.pdf"
        create_elbow_plot(summary_df, elbow_path, recommended_threshold)
        results['elbow_plot'] = elbow_path
    except Exception as e:
        logger.warning(f"Could not generate recommendations: {e}")

    # Create README
    readme_path = output_dir / "README.md"
    rec_display = results.get('recommended_threshold', 'see recommendations.txt')
    readme_content = f"""# Parameter Sweep Results — {organism}

## What was swept

`min_singleton_distance` — the divergence cutoff below which singleton
haplotypes (single-member ESVs) are removed as likely sequencing or PCR
errors.  Higher values filter more aggressively; lower values retain more
rare variants.

## Recommended value

**{rec_display}**  (see `recommendations.txt` for full rationale)

## Files

| File | Description |
|------|-------------|
| `sweep_summary.csv` | Per-threshold metrics table |
| `threshold_stability.pdf` | 4-panel stability visualisation |
| `elbow_plot.pdf` | Optimal-threshold detection |
| `group_membership_tracking.csv` | Sample-level haplotype stability |
| `recommendations.txt` | Automated interpretation and next steps |
| `runs/` | Per-threshold pipeline outputs |

## Reading group_membership_tracking.csv

Each row is one sample.  Columns `t_X.XXX` show which haplotype the sample
was assigned to at each tested threshold.  The key derived columns are:

- **n_changes** — number of times the sample's haplotype *composition*
  changed significantly (Jaccard similarity of cluster partners < 0.7).
  Low values (0–2) indicate stable assignment.
- **stability_score** — `high` (0 changes), `medium` (1–2), or `low` (>2).
  Samples scored `low` are worth investigating manually.
- **mean / min / max_cluster_size** — size of the sample's haplotype across
  the sweep.

> Haplotype *names* (e.g. `h1_n177`) are arbitrary and may differ between
> runs.  What matters is which other samples share the same haplotype.

## Interpretation

See `recommendations.txt`.  The recommended threshold is the elbow point
where haplotype count stops dropping sharply — balancing error removal
against retention of biologically real rare variants.

## Generated

{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info("="*70)
    logger.info("Parameter sweep complete!")

    rec_threshold = results.get('recommended_threshold', None)
    if rec_threshold is not None and isinstance(rec_threshold, (int, float)):
        logger.info(f"  Recommended min_singleton_distance: {rec_threshold:.3f}")
    else:
        logger.info(f"  Recommended min_singleton_distance: N/A")

    logger.info(f"  Results: {output_dir}")
    logger.info("="*70)

    return results
