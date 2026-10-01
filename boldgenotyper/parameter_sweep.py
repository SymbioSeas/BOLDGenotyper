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


def _load_run_tables(output_dir: Path, organism: str):
    """
    (haplotype table or None, QC-passed record table or None) of a run folder:
    v2 haplotypes.csv and the qc_pass rows of records.csv, or v1
    haplotypes/<org>_haplotype_stats.csv and <org>_annotated.csv.
    """
    output_dir = Path(output_dir)
    hap_path = output_dir / "haplotypes.csv"
    if not hap_path.exists():
        hap_path = output_dir / "haplotypes" / f"{organism}_haplotype_stats.csv"
    haplotypes = pd.read_csv(hap_path) if hap_path.exists() else None

    records_path = output_dir / "records.csv"
    if records_path.exists():
        records = pd.read_csv(records_path, low_memory=False)
        records = records[records['qc_pass'].astype(str).str.lower() == 'true']
    else:
        annotated_path = output_dir / f"{organism}_annotated.csv"
        records = pd.read_csv(annotated_path, low_memory=False) if annotated_path.exists() else None
    return haplotypes, records


def extract_run_metrics(
    output_dir: Path,
    organism: str,
    threshold: float
) -> Dict[str, Any]:
    """
    Extract summary metrics from a completed run.

    Reads haplotypes.csv (haplotype counts) and the QC-passed rows of
    records.csv (sample-level metrics); v1 run folders: the haplotype stats
    table and <organism>_annotated.csv.

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
    stats_df, df = _load_run_tables(output_dir, organism)
    if stats_df is not None:
        metrics['n_haplotypes'] = len(stats_df)
        if 'is_singleton' in stats_df.columns:
            metrics['n_singletons'] = int(stats_df['is_singleton'].sum())
        if 'is_suspect' in stats_df.columns:
            metrics['n_suspect'] = int(stats_df['is_suspect'].sum())

    # ---------------------------------------------------------------------------
    # Sample-level metrics (from annotated CSV)
    # ---------------------------------------------------------------------------
    if df is None:
        return metrics

    metrics['n_samples'] = len(df)

    # Assignment coverage: fraction of samples assigned to a haplotype
    if 'haplotype_id' in df.columns:
        assigned = df[df['haplotype_id'].notna()]
        metrics['pct_assigned'] = (len(assigned) / len(df) * 100) if len(df) > 0 else 0

    # Species richness (from depositor taxonomy — constant across thresholds,
    # included for completeness in the summary table)
    if 'species' in df.columns:
        metrics['n_species_labels_qc_passed'] = df['species'].nunique()  # P0-9 name

    # Label conflicts: records whose label differs from their haplotype's
    # majority label (quality_control.add_label_columns). Can change when
    # singleton removal causes samples to be reassigned.
    if 'label_conflict' in df.columns:
        conflict = df['label_conflict'].astype(str).str.lower() == 'true'
        metrics['pct_label_conflict'] = (conflict.sum() / len(df) * 100) if len(df) > 0 else 0
        if 'haplotype_id' in df.columns:
            per_hap = conflict.groupby(df['haplotype_id']).any()
            metrics['n_haplotypes_with_label_conflict'] = int(per_hap.sum())

    # Mean majority fraction across haplotypes
    if stats_df is not None and 'majority_fraction' in stats_df.columns:          # v2
        metrics['mean_majority_fraction'] = float(stats_df['majority_fraction'].mean())
    elif 'majority_fraction' in df.columns and 'haplotype_id' in df.columns:
        fractions = df.groupby('haplotype_id')['majority_fraction'].first()
        metrics['mean_majority_fraction'] = float(fractions.mean()) if len(fractions) > 0 else 0

    return metrics


def cleanup_intermediate_files(output_dir: Path) -> None:
    """
    Remove per-threshold output directories that are no longer needed after
    metrics have been extracted.  Retains the record and haplotype tables
    (records.csv, haplotypes.csv; v1: the annotated CSV and haplotype stats).

    Parameters
    ----------
    output_dir : Path
        Output directory to clean
    """
    # These directories are produced by the sweep runs but are not needed
    # once extract_run_metrics has pulled the numbers.
    removable = [
        'divergence', 'phylogeny', 'geography', 'extended', 'exports', 'plots',
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
    sample_md5 = {}          # threshold -> {sample_id: seq_md5 of its haplotype}
    group_members = {}       # threshold -> {haplotype_id: set(sample_ids)}

    for threshold in sorted(thresholds):
        output_dir = run_outputs[threshold]
        haps, df = _load_run_tables(output_dir, organism)

        if df is not None:
            if 'processid' in df.columns and 'haplotype_id' in df.columns:
                assignments = dict(zip(df['processid'], df['haplotype_id']))
                sample_assignments[threshold] = assignments
                # Haplotype IDs are renumbered in every run; the sequence key is comparable
                key = (dict(zip(haps['haplotype_id'], haps['seq_md5']))
                       if haps is not None and 'seq_md5' in haps.columns else {})
                sample_md5[threshold] = {
                    s: (key.get(h, h) if pd.notna(h) else '') for s, h in assignments.items()}

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
                row[f'md5_{threshold}'] = sample_md5[threshold].get(sample_id, 'NA')

                if pd.notna(hap) and hap != 'unassigned' and threshold in group_members:
                    partners = group_members[threshold].get(hap, set()) - {sample_id}
                    cluster_partners_by_threshold.append((threshold, partners))
                else:
                    cluster_partners_by_threshold.append((threshold, set()))
            else:
                row[f't_{threshold}'] = 'NA'
                row[f'md5_{threshold}'] = 'NA'

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

    if 'pct_label_conflict' in summary_df.columns:
        ax2 = ax.twinx()
        ax2.plot(summary_df['threshold'], summary_df['pct_label_conflict'],
                 's--', linewidth=2, markersize=6, color='tab:red', alpha=0.7,
                 label='% Label conflict')
        ax2.set_ylabel('% Records with a label conflict', color='tab:red')
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


def adjacent_change_rates(tracking_df: pd.DataFrame, thresholds: List[float]) -> pd.DataFrame:
    """
    Share of records whose haplotype changes between adjacent thresholds (P3).

    A record's haplotype is compared by sequence (the md5_<threshold> columns
    of group_membership_tracking.csv), so renumbered IDs do not count as
    changes; becoming assigned or unassigned does. Records missing from either
    run ('NA') are left out.
    """
    rows = []
    ordered = sorted(thresholds)
    for lo, hi in zip(ordered, ordered[1:]):
        a, b = tracking_df.get(f'md5_{lo}'), tracking_df.get(f'md5_{hi}')
        if a is None or b is None:
            continue
        a, b = a.fillna(''), b.fillna('')
        both = (a != 'NA') & (b != 'NA')
        n = int(both.sum())
        changed = int((both & (a != b)).sum())
        rows.append({'threshold_from': lo, 'threshold_to': hi, 'n_records': n,
                     'n_changed': changed, 'pct_changed': round(changed / n * 100, 2) if n else 0.0})
    return pd.DataFrame(rows, columns=['threshold_from', 'threshold_to', 'n_records', 'n_changed', 'pct_changed'])


def generate_recommendations(
    summary_df: pd.DataFrame,
    organism: str,
    output_path: Path,
    tracking_df: Optional[pd.DataFrame] = None,
) -> float:
    """
    Write recommendations.txt: the elbow of the haplotype-count curve and, for
    each pair of adjacent thresholds, the share of records whose haplotype
    changes (from group_membership_tracking.csv). No acceptable range is
    given (v1's covered the whole tested range).

    Returns the threshold at the elbow.
    """
    logger.info("  Generating recommendations...")

    recommended_threshold = find_elbow_point(summary_df)
    rows = summary_df.to_dict('records')          # keeps integer counts as integers
    rec_row = next(r for r in rows if r['threshold'] == recommended_threshold)
    first_row, last_row = rows[0], rows[-1]
    lo, hi = summary_df['threshold'].min(), summary_df['threshold'].max()

    def fmt(value, spec):
        if value is None or pd.isna(value):
            return 'N/A'
        return format(int(value) if spec.endswith('d') else value, spec)

    text = f"""Singleton Distance Threshold: Sweep Results
=============================================

Dataset:                  {organism}
QC-passed records:        {fmt(rec_row.get('n_samples'), 'd')}
Threshold range tested:   {lo:.3f} to {hi:.3f}
Number of thresholds:     {len(summary_df)}

WHAT WAS SWEPT
--------------
--singleton-distance (min_singleton_distance): singleton haplotypes within this
p-distance of another haplotype are treated as likely sequencing or PCR errors,
and their records are assigned to the nearest retained haplotype when it is
unique and within the same distance.

HAPLOTYPE COUNT
---------------
   threshold   haplotypes   singletons   assigned (% of QC-passed)
"""
    for r in rows:
        mark = "   <- elbow" if r['threshold'] == recommended_threshold else ""
        text += (f"   {r['threshold']:.3f}       {fmt(r.get('n_haplotypes'), '>6d')}       "
                 f"{fmt(r.get('n_singletons'), '>6d')}       {fmt(r.get('pct_assigned'), '>6.1f')}{mark}\n")

    text += f"""
The elbow (largest change in the slope of the haplotype-count curve) is at
{recommended_threshold:.3f}: {fmt(rec_row['n_haplotypes'], 'd')} haplotypes, against {fmt(first_row['n_haplotypes'], 'd')} at {lo:.3f} and {fmt(last_row['n_haplotypes'], 'd')} at {hi:.3f}.
"""

    rates = adjacent_change_rates(tracking_df, summary_df['threshold'].tolist()) \
        if tracking_df is not None and len(tracking_df) else pd.DataFrame()
    if len(rates):
        text += """
STABILITY
---------
Share of records whose haplotype changes between adjacent thresholds
(compared by sequence; becoming assigned or unassigned counts as a change):

"""
        for _, r in rates.iterrows():
            text += (f"   {r['threshold_from']:.3f} -> {r['threshold_to']:.3f}   "
                     f"{int(r['n_changed']):>6} of {int(r['n_records'])} records   {r['pct_changed']:.1f}%\n")
        text += """
Thresholds where this share is small give nearly the same haplotypes for
nearly every record. Per-record detail: group_membership_tracking.csv.
"""
    text += f"""
NEXT STEPS
----------
1. Review elbow_plot.pdf, threshold_stability.pdf and the stability table above.
   The elbow is a starting point; the choice of threshold is yours.

2. Run the pipeline with the chosen threshold, passed explicitly so the run is
   reproducible (the default is 0.005):

       boldgenotyper your_data.tsv --organism {organism} \\
                     --singleton-distance {recommended_threshold:.3f} [other options]

3. Report the threshold and how it was chosen in your methods.

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
            summary_df, organism, recommendations_path,
            tracking_df=tracking_df if len(tracking_df) > 0 else None,
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
was assigned to at each tested threshold, and `md5_X.XXX` the seq_md5 of that
haplotype (empty when unassigned, NA when the record is not in that run).
recommendations.txt uses the md5 columns to report, for each pair of adjacent
thresholds, the share of records whose haplotype changes.  The key derived
columns are:

- **n_changes** — number of times the sample's haplotype *composition*
  changed significantly (Jaccard similarity of cluster partners < 0.7).
  Low values (0–2) indicate stable assignment.
- **stability_score** — `high` (0 changes), `medium` (1–2), or `low` (>2).
  Samples scored `low` are worth investigating manually.
- **mean / min / max_cluster_size** — size of the sample's haplotype across
  the sweep.

> Haplotype *names* (e.g. `H0001`) are numbered within each run and differ
> between thresholds; compare haplotypes across runs by `seq_md5`.  What
> matters is which other samples share the same haplotype.

## Interpretation

See `recommendations.txt`: the elbow of the haplotype-count curve, and the
share of records whose haplotype changes between adjacent thresholds.

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
