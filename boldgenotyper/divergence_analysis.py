"""
Divergence Analysis Module for BOLDGenotyper

This module provides automatic calculation and reporting of divergence statistics
between and within consensus groups. It helps assess species boundaries, detect
cryptic lineages, and evaluate the barcoding gap.

Key Features:
- Pairwise divergence matrix calculation
- Within-species vs between-species divergence comparison
- Barcoding gap detection and visualization
- Divergence heatmaps
- Species delimitation guidance

Author: Steph Smith (symbioseas@outlook.com)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .taxon_labels import haplotype_species_map

logger = logging.getLogger(__name__)


class DivergenceAnalysisError(Exception):
    """Base exception for divergence analysis errors."""
    pass


def calculate_pairwise_distance(seq1: str, seq2: str) -> float:
    """
    Calculate uncorrected p-distance between two sequences.

    Distance = 1 - (matches / informative_positions)
    Gaps are ignored in the calculation.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    float
        Uncorrected p-distance
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")

    matches = 0
    informative = 0

    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        # Skip positions with gaps or ambiguous bases
        if base1 in '-N' or base2 in '-N':
            continue

        informative += 1
        if base1 == base2:
            matches += 1

    if informative == 0:
        return np.nan

    distance = 1 - (matches / informative)
    return distance


def align_consensus_sequences(
    consensus_fasta: Path,
    output_dir: Path,
    organism: str
) -> Path:
    """
    Align consensus sequences using MAFFT if not already aligned.

    Parameters
    ----------
    consensus_fasta : Path
        Path to consensus sequences FASTA
    output_dir : Path
        Output directory for alignment
    organism : str
        Organism name

    Returns
    -------
    Path
        Path to aligned FASTA file
    """
    import subprocess
    import shutil

    # Check if MAFFT is available
    if not shutil.which('mafft'):
        logger.warning("MAFFT not found. Using unaligned sequences (may affect accuracy).")
        return consensus_fasta

    output_dir.mkdir(parents=True, exist_ok=True)
    aligned_path = output_dir / f"{organism}_consensus_aligned.fasta"

    logger.info("  Aligning consensus sequences with MAFFT...")

    try:
        # Run MAFFT
        cmd = [
            'mafft',
            '--auto',
            '--quiet',
            str(consensus_fasta)
        ]

        with open(aligned_path, 'w') as outf:
            result = subprocess.run(
                cmd,
                stdout=outf,
                stderr=subprocess.PIPE,
                check=True
            )

        logger.info(f"  Aligned sequences saved to {aligned_path}")
        return aligned_path

    except subprocess.CalledProcessError as e:
        logger.warning(f"MAFFT alignment failed: {e}", extra={"step_failure": "9 divergence"})
        logger.warning("Using unaligned sequences (may affect accuracy)")
        return consensus_fasta


def calculate_pairwise_divergence_matrix(
    aligned_fasta: Path,
    taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate pairwise divergence matrix for all consensus sequences.

    Parameters
    ----------
    aligned_fasta : Path
        Path to aligned consensus sequences
    taxonomy_df : pd.DataFrame
        Taxonomy table with consensus_group and assigned_sp columns

    Returns
    -------
    pd.DataFrame
        Square matrix of pairwise divergences
    """
    logger.info("  Calculating pairwise divergence matrix...")

    # Load sequences
    sequences = {}
    for record in SeqIO.parse(aligned_fasta, 'fasta'):
        seq_id = record.id
        sequences[seq_id] = str(record.seq)

    # Build label map: haplotype_id -> haplotype_sp (e.g. "H0001" -> "Salmo salar H0001")
    id_to_sp = {}
    if 'haplotype_id' in taxonomy_df.columns and 'haplotype_sp' in taxonomy_df.columns:
        id_to_sp = dict(zip(taxonomy_df['haplotype_id'], taxonomy_df['haplotype_sp']))

    # Create row labels using haplotype_sp if available, else keep raw seq_id
    labels = []
    seq_ids = list(sequences.keys())

    for seq_id in seq_ids:
        label = id_to_sp.get(seq_id, seq_id)
        labels.append(label)

    # Calculate pairwise distances
    n = len(seq_ids)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            if i == j:
                distance_matrix[i, j] = 0.0
            else:
                dist = calculate_pairwise_distance(
                    sequences[seq_ids[i]],
                    sequences[seq_ids[j]]
                )
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist

    # Create dataframe
    df = pd.DataFrame(
        distance_matrix,
        index=labels,
        columns=labels
    )

    return df


def generate_divergence_summary(
    divergence_matrix: pd.DataFrame,
    taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Generate summary statistics on divergence.

    Parameters
    ----------
    divergence_matrix : pd.DataFrame
        Pairwise divergence matrix
    taxonomy_df : pd.DataFrame
        Taxonomy information

    Returns
    -------
    pd.DataFrame
        Summary table with divergence statistics
    """
    logger.info("  Generating divergence summary...")

    # Extract species information from row labels
    results = []

    labels = divergence_matrix.index.tolist()

    # Build label→species mapping from taxonomy_df
    label_to_species = {}
    if 'haplotype_sp' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        label_to_species = dict(zip(taxonomy_df['haplotype_sp'], taxonomy_df['assigned_sp']))
    # Also map by haplotype_id in case labels are raw IDs
    if 'haplotype_id' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        id_to_sp = dict(zip(taxonomy_df['haplotype_id'], taxonomy_df['assigned_sp']))
        for hid, sp in id_to_sp.items():
            if hid not in label_to_species:
                label_to_species[hid] = sp

    def get_group_species(label):
        """Extract group and species from label using taxonomy mapping."""
        species = label_to_species.get(label, '')
        return label, species

    group_species_map = {label: get_group_species(label) for label in labels}

    # Calculate within-species and between-species divergences
    for i, label1 in enumerate(labels):
        for j, label2 in enumerate(labels):
            if i >= j:  # Skip diagonal and lower triangle
                continue

            group1, sp1 = group_species_map[label1]
            group2, sp2 = group_species_map[label2]

            divergence = divergence_matrix.loc[label1, label2]

            if pd.isna(divergence):
                continue

            # Determine comparison type
            if sp1 == sp2 and sp1:
                comparison_type = 'within_species'
                species_pair = sp1
            else:
                comparison_type = 'between_species'
                species_pair = f"{sp1} vs {sp2}" if sp1 and sp2 else "unknown"

            results.append({
                'comparison_type': comparison_type,
                'group1': group1,
                'group2': group2,
                'species1': sp1,
                'species2': sp2,
                'species_pair': species_pair,
                'divergence': divergence,
            })

    summary_df = pd.DataFrame(results)

    # Calculate aggregate statistics
    if len(summary_df) > 0:
        agg_summary = summary_df.groupby(
            ['comparison_type', 'species_pair']
        )['divergence'].agg(['mean', 'min', 'max'])

        agg_summary = agg_summary.reset_index()
        agg_summary.columns = [
            'comparison_type', 'species_pair', 'mean_divergence',
            'min_divergence', 'max_divergence'
        ]

        return agg_summary
    else:
        return pd.DataFrame()


def analyze_barcoding_gap(
    divergence_matrix: pd.DataFrame,
    taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Deprecated (P1-2): within/between divergence of haplotype majority labels.

    Superseded by compute_label_based_gap (record labels, label-conflict
    records excluded), which the pipeline writes only with --extended.

    Analyze within-species vs between-species divergence (barcoding gap).

    Parameters
    ----------
    divergence_matrix : pd.DataFrame
        Pairwise divergence matrix
    taxonomy_df : pd.DataFrame
        Taxonomy information

    Returns
    -------
    pd.DataFrame
        Barcoding gap analysis per species
    """
    logger.info("  Analyzing barcoding gap...")

    labels = divergence_matrix.index.tolist()

    # Build label→species mapping from taxonomy_df
    label_to_species = {}
    if 'haplotype_sp' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        label_to_species = dict(zip(taxonomy_df['haplotype_sp'], taxonomy_df['assigned_sp']))
    if 'haplotype_id' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        id_to_sp = dict(zip(taxonomy_df['haplotype_id'], taxonomy_df['assigned_sp']))
        for hid, sp in id_to_sp.items():
            if hid not in label_to_species:
                label_to_species[hid] = sp

    group_species_map = {label: (label, label_to_species.get(label, '')) for label in labels}

    # Collect within and between species divergences
    within_species = {}
    between_species = {}

    for i, label1 in enumerate(labels):
        for j, label2 in enumerate(labels):
            if i >= j:
                continue

            group1, sp1 = group_species_map[label1]
            group2, sp2 = group_species_map[label2]

            divergence = divergence_matrix.loc[label1, label2]

            if pd.isna(divergence) or not sp1:
                continue

            if sp1 == sp2:
                if sp1 not in within_species:
                    within_species[sp1] = []
                within_species[sp1].append(divergence)
            else:
                key = tuple(sorted([sp1, sp2]))
                if key not in between_species:
                    between_species[key] = []
                between_species[key].append(divergence)

    # Calculate statistics per species
    results = []

    for species in within_species.keys():
        within = within_species.get(species, [])

        # Get between-species divergences involving this species
        between = []
        for (sp1, sp2), divs in between_species.items():
            if species in (sp1, sp2):
                between.extend(divs)

        if not within:
            continue

        within_mean = np.mean(within)
        within_std = np.std(within)
        within_max = np.max(within)

        if between:
            between_mean = np.mean(between)
            between_std = np.std(between)
            between_min = np.min(between)

            # Check for barcoding gap
            has_gap = between_min > within_max
            gap_size = between_min - within_max if has_gap else 0
        else:
            between_mean = np.nan
            between_std = np.nan
            between_min = np.nan
            has_gap = False
            gap_size = 0

        results.append({
            'species': species,
            'n_within_comparisons': len(within),
            'mean_within': within_mean,
            'sd_within': within_std,
            'max_within': within_max,
            'n_between_comparisons': len(between),
            'mean_between': between_mean,
            'sd_between': between_std,
            'min_between': between_min,
            'barcoding_gap': has_gap,
            'gap_size': gap_size
        })

    return pd.DataFrame(results)


def create_divergence_heatmap(
    divergence_matrix: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Create clustered heatmap of pairwise divergence matrix with dendrograms.

    Generates a hierarchical clustered heatmap with dendrograms showing
    relationships between haplotypes. Uses complete linkage clustering
    and annotates cells with divergence values (0 shown as "×").

    Parameters
    ----------
    divergence_matrix : pd.DataFrame
        Pairwise divergence matrix (will be symmetrized)
    output_path : Path
        Output path for PDF
    """
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform

    logger.info("  Creating divergence heatmap with dendrograms...")

    # Symmetrize matrix (handle any asymmetry from rounding)
    mat = (divergence_matrix + divergence_matrix.T) / 2

    # Build distance array for hierarchical clustering
    dist_array = squareform(mat.values)
    row_linkage = linkage(dist_array, method="complete")
    col_linkage = linkage(dist_array, method="complete")

    # Adaptive sizing based on matrix dimensions
    n = len(mat)
    if n <= 20:
        figsize = (max(9, n * 0.5), max(9, n * 0.5))
        annotation_fontsize = 8
        tick_fontsize = 9
    elif n <= 50:
        figsize = (max(12, n * 0.35), max(12, n * 0.35))
        annotation_fontsize = 6
        tick_fontsize = 7
    elif n <= 100:
        figsize = (max(15, n * 0.25), max(15, n * 0.25))
        annotation_fontsize = 5
        tick_fontsize = 6
    else:
        figsize = (max(20, n * 0.2), max(20, n * 0.2))
        annotation_fontsize = 4
        tick_fontsize = 5

    # Color map settings
    cmap = sns.color_palette("mako_r", as_cmap=True)
    vmin, vmax = 0, mat.values.max()

    # Create clustermap with dendrograms
    cg = sns.clustermap(
        mat,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="lightgray",
        figsize=figsize,
        cbar_kws={"label": "Pairwise divergence"}
    )

    # Style the heatmap
    ax = cg.ax_heatmap
    ax.set_facecolor("white")
    cg.cax.set_visible(False)  # Hide colorbar for cleaner look

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

    # Get the clustered data
    data = cg.data2d
    n_rows, n_cols = data.shape

    for i in range(n_rows):
        for j in range(n_cols):
            val = data.iloc[i, j]
            if val == 0:
                # Diagonal or identical sequences - show ×
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    "×",
                    ha="center",
                    va="center",
                    fontsize=annotation_fontsize,
                    color="black",
                )
            else:
                # Show divergence value with adaptive color
                txt_color = get_text_color(val)
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    f"{val:.3f}",
                    ha="center",
                    va="center",
                    fontsize=annotation_fontsize,
                    color=txt_color,
                )

    # Rotate tick labels for readability
    for tick in ax.get_xticklabels():
        tick.set_fontsize(tick_fontsize)
        tick.set_rotation(45)
        tick.set_ha("right")

    for tick in ax.get_yticklabels():
        tick.set_fontsize(tick_fontsize)

    # Add title
    plt.subplots_adjust(top=0.96)
    cg.fig.suptitle("Haplotype Divergence (Uncorrected P-Distance)", fontsize=11)

    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  Divergence heatmap saved to {output_path}")


def create_barcoding_gap_plot(
    divergence_matrix: pd.DataFrame,
    taxonomy_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Create barcoding gap visualization.

    Parameters
    ----------
    divergence_matrix : pd.DataFrame
        Pairwise divergence matrix
    taxonomy_df : pd.DataFrame
        Taxonomy information
    output_path : Path
        Output path for PDF
    """
    logger.info("  Creating barcoding gap plot...")

    labels = divergence_matrix.index.tolist()

    # Build label→species mapping from taxonomy_df
    label_to_species = {}
    if 'haplotype_sp' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        label_to_species = dict(zip(taxonomy_df['haplotype_sp'], taxonomy_df['assigned_sp']))
    if 'haplotype_id' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        id_to_sp = dict(zip(taxonomy_df['haplotype_id'], taxonomy_df['assigned_sp']))
        for hid, sp in id_to_sp.items():
            if hid not in label_to_species:
                label_to_species[hid] = sp

    group_species_map = {label: (label, label_to_species.get(label, '')) for label in labels}

    # Collect divergences
    within_divergences = []
    between_divergences = []

    for i, label1 in enumerate(labels):
        for j, label2 in enumerate(labels):
            if i >= j:
                continue

            group1, sp1 = group_species_map[label1]
            group2, sp2 = group_species_map[label2]

            divergence = divergence_matrix.loc[label1, label2]

            if pd.isna(divergence):
                continue

            if sp1 == sp2 and sp1:
                within_divergences.append(divergence)
            elif sp1 and sp2:
                between_divergences.append(divergence)

    if not within_divergences and not between_divergences:
        logger.warning("No divergence data for barcoding gap plot")
        return

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot histograms
    bins = np.linspace(0, max(max(within_divergences or [0]), max(between_divergences or [0])) * 1.1, 50)

    if within_divergences:
        ax.hist(within_divergences, bins=bins, alpha=0.6, color='blue',
                label=f'Within-species (n={len(within_divergences)})', edgecolor='black')
        within_mean = np.mean(within_divergences)
        ax.axvline(within_mean, color='blue', linestyle='--', linewidth=2,
                   label=f'Within mean: {within_mean:.3f}')

    if between_divergences:
        ax.hist(between_divergences, bins=bins, alpha=0.6, color='red',
                label=f'Between-species (n={len(between_divergences)})', edgecolor='black')
        between_mean = np.mean(between_divergences)
        ax.axvline(between_mean, color='red', linestyle='--', linewidth=2,
                   label=f'Between mean: {between_mean:.3f}')

    # Check for gap
    if within_divergences and between_divergences:
        within_max = max(within_divergences)
        between_min = min(between_divergences)

        if between_min > within_max:
            # Shade the gap
            ax.axvspan(within_max, between_min, alpha=0.3, color='green',
                      label=f'Barcoding gap: {(between_min - within_max):.3f}')

    ax.set_xlabel('Divergence (p-distance)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title('Barcoding Gap Analysis', fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  Barcoding gap plot saved to {output_path}")


LABEL_GAP_COLUMNS = ["label", "n_records", "n_haplotypes", "max_within", "min_between",
                     "nearest_other_label", "label_based_gap"]


def compute_label_based_gap(distance_by_id: pd.DataFrame, records: pd.DataFrame) -> pd.DataFrame:
    """Within- and between-label divergence from record labels (P1-2, --extended).

    For each species-level record label: the haplotypes of records carrying it
    (label-conflict records excluded when a label_conflict column is present),
    the largest p-distance among those haplotypes (max_within; 0 for one
    haplotype), the smallest p-distance to a haplotype of any other label
    (min_between) and that label, and label_based_gap = min_between >
    max_within. Distances between records are the distances between their
    haplotypes. The result depends on the labels deposited in BOLD and on
    sampling; it is not a species-delimitation test.
    """
    from .taxon_labels import label_level

    rec = records[records["haplotype_id"].notna()].copy()
    rec["label"] = rec["label"].fillna("")
    rec = rec[rec["label"].map(label_level) == "species"]
    if "label_conflict" in rec.columns:
        rec = rec[~rec["label_conflict"].astype(str).str.lower().eq("true")]
    rec = rec[rec["haplotype_id"].isin(distance_by_id.index)]
    haps = rec.groupby("label")["haplotype_id"].agg(lambda s: sorted(set(s)))
    n_records = rec.groupby("label").size()
    rows = []
    for label, h in haps.items():
        within = float(distance_by_id.loc[h, h].to_numpy().max()) if len(h) > 1 else 0.0
        best, best_label = np.inf, None
        for other, g in haps.items():
            if other == label:
                continue
            d = float(distance_by_id.loc[h, g].to_numpy().min())
            if d < best:
                best, best_label = d, other
        min_between = best if np.isfinite(best) else np.nan
        rows.append({
            "label": label, "n_records": int(n_records[label]), "n_haplotypes": len(h),
            "max_within": within, "min_between": min_between, "nearest_other_label": best_label,
            "label_based_gap": (min_between > within) if np.isfinite(best) else np.nan,
        })
    return pd.DataFrame(rows, columns=LABEL_GAP_COLUMNS)


def create_label_gap_plot(gap_df: pd.DataFrame, output_path: Path) -> None:
    """Scatter of max_within against min_between per label, with the 1:1 line."""
    import matplotlib.pyplot as plt
    data = gap_df.dropna(subset=["min_between"])
    if data.empty:
        logger.info("  No labels with a between-label comparison; label gap plot skipped")
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(data["max_within"], data["min_between"], s=30, alpha=0.8, color="#4C72B0",
               edgecolors="black", linewidths=0.4)
    top = float(max(data["max_within"].max(), data["min_between"].max(), 0.01)) * 1.05
    ax.plot([0, top], [0, top], color="gray", linestyle="-", linewidth=1)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel("Largest p-distance among haplotypes with the label (max_within)")
    ax.set_ylabel("Smallest p-distance to another label (min_between)")
    ax.set_title("Label-based within vs between divergence (one point per label)")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def generate_divergence_analysis(
    consensus_fasta: Path,
    taxonomy_csv: Path,
    output_dir: Path,
    organism: str,
    extended: bool = False,
    records: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    """
    Complete divergence analysis workflow.

    With extended=True and records (the annotated table with haplotype_id,
    label and label_conflict), also writes label_based_gap.csv and its plot.

    Parameters
    ----------
    consensus_fasta : Path
        Path to consensus sequences FASTA
    taxonomy_csv : Path
        Path to consensus taxonomy CSV
    output_dir : Path
        Output directory for divergence analysis
    organism : str
        Organism name

    Returns
    -------
    dict
        Results and file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating divergence analysis...")

    results = {}

    # Load taxonomy
    taxonomy_df = pd.read_csv(taxonomy_csv)

    # Align sequences if needed
    aligned_fasta = align_consensus_sequences(
        consensus_fasta,
        output_dir,
        organism
    )

    # Calculate pairwise divergence matrix
    divergence_matrix = calculate_pairwise_divergence_matrix(
        aligned_fasta,
        taxonomy_df
    )

    # Save matrix. The label-keyed copy stays in staging for the species
    # analyses (labels carry the taxon); the run folder gets the copy keyed by
    # haplotype_id, which does not change when a taxon assignment does.
    matrix_path = output_dir / "pairwise_divergence_matrix.csv"
    divergence_matrix.to_csv(matrix_path)
    results['divergence_matrix'] = matrix_path
    to_id = {}
    if 'haplotype_sp' in taxonomy_df.columns and 'haplotype_id' in taxonomy_df.columns:
        to_id = dict(zip(taxonomy_df['haplotype_sp'], taxonomy_df['haplotype_id']))
    distance_by_id = divergence_matrix.rename(index=lambda x: to_id.get(x, x),
                                              columns=lambda x: to_id.get(x, x))
    distance_by_id.index.name = 'haplotype_id'
    by_id_path = output_dir / "pairwise_p_distance_by_id.csv"
    distance_by_id.to_csv(by_id_path)
    results['divergence_matrix_by_id'] = by_id_path
    logger.info(f"  Saved divergence matrix: {by_id_path}")

    # Generate divergence summary
    summary_df = generate_divergence_summary(divergence_matrix, taxonomy_df)
    if len(summary_df) > 0:
        summary_path = output_dir / "divergence_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        results['divergence_summary'] = summary_path
        logger.info(f"  Saved divergence summary: {summary_path}")

    # Label-based barcoding gap: only with --extended (P1-2)
    if extended and records is not None:
        gap_df = compute_label_based_gap(distance_by_id, records)
        gap_path = output_dir / "label_based_gap.csv"
        gap_df.to_csv(gap_path, index=False)
        results['label_based_gap'] = gap_path
        results['gap_stats'] = gap_df.to_dict('records')
        logger.info(f"  Saved label-based gap table: {gap_path}")

    # Create visualizations
    try:
        heatmap_path = output_dir / "divergence_heatmap.pdf"
        create_divergence_heatmap(divergence_matrix, heatmap_path)
        results['heatmap'] = heatmap_path
    except Exception as e:
        logger.warning(f"Could not create divergence heatmap: {e}", extra={"step_failure": "9 divergence"})

    if 'label_based_gap' in results:
        try:
            gap_plot_path = output_dir / "label_based_gap.pdf"
            create_label_gap_plot(gap_df, gap_plot_path)
            if gap_plot_path.exists():
                results['gap_plot'] = gap_plot_path
        except Exception as e:
            logger.warning(f"Could not create label gap plot: {e}", extra={"step_failure": "9 divergence"})

    # Create README
    gap_files_text = ""
    if 'label_based_gap' in results:
        gap_files_text = (
            "- **label_based_gap.csv** (with --extended): per species-level record label, the\n"
            "  largest p-distance among haplotypes carrying it (`max_within`), the smallest\n"
            "  p-distance to a haplotype of another label (`min_between`, `nearest_other_label`)\n"
            "  and `label_based_gap` (min_between > max_within). Label-conflict records are\n"
            "  excluded. It reflects the labels deposited in BOLD and sampling; it is not a\n"
            "  species-delimitation test.\n"
            "- **label_based_gap.pdf** (with --extended): max_within against min_between per label\n"
        )
    readme_path = output_dir / "README.md"
    readme_content = f"""# Divergence Analysis for {organism}

## Overview

This directory contains pairwise divergence analysis results for haplotype sequences.

## Files

- **pairwise_divergence_matrix.csv**: Square matrix of uncorrected p-distances
- **divergence_summary.csv**: Summary statistics by comparison type (mean, min, max divergence)
- **divergence_heatmap.pdf**: Visual heatmap of pairwise divergences
{gap_files_text}

## Haplotype Label Format

Haplotype IDs (`H0001`, `H0002`, ...) carry no counts. They are numbered by
member count (largest first), then by `seq_md5`. Member counts are in the
`n_members` column of the haplotype table; `seq_md5` (the md5 of the
uppercase, ungapped core sequence) is the stable key for comparing haplotypes
across runs.

## FASTA Files

- **haplotypes.fasta** (in haplotypes/): ungapped haplotype sequences for downstream use
- **{organism}_consensus_aligned.fasta** (this folder): the aligned sequences used for the
  divergence calculation

## Methodology

Divergences are uncorrected p-distances between haplotype sequences:
```
distance = 1 - (matches / informative_positions)
```

Gaps and ambiguous bases (N) are excluded from calculations. Species groupings
come from each haplotype's majority label; the distances are reported without
taxonomic thresholds, because within- and between-species COI divergence varies
across taxa.
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info(f"  Divergence analysis complete")

    return results


def print_divergence_summary(results: Dict[str, Any]) -> None:
    """
    Print divergence analysis summary to console.

    Parameters
    ----------
    results : dict
        Results from generate_divergence_analysis
    """
    if 'gap_stats' not in results:
        return

    print("\n" + "="*70)
    print("DIVERGENCE ANALYSIS")
    print("="*70)

    n_gap = sum(1 for g in results['gap_stats'] if g.get('label_based_gap') is True)
    n_compared = sum(1 for g in results['gap_stats'] if not pd.isna(g.get('min_between')))
    print(f"Label-based gap (min_between > max_within): {n_gap} of {n_compared} labels "
          f"with a between-label comparison")
    print("  See divergence_analysis/label_based_gap.csv (record labels; depends on sampling)")
    print("="*70 + "\n")


# ============================================================================
# Species-Level Divergence Analysis
# ============================================================================


def generate_species_divergence_analysis(
    divergence_matrix: pd.DataFrame,
    species_summary_csv: Path,
    haplotype_taxonomy_csv: Path,
    output_dir: Path
) -> Dict[str, Any]:
    """
    Analyze genetic divergence at the species level.

    Uses existing haplotype-level divergence matrix and species assignments
    to calculate within-species and between-species divergence.

    Parameters
    ----------
    divergence_matrix : pd.DataFrame
        Pairwise divergence matrix from haplotype-level analysis
    species_summary_csv : Path
        Species summary CSV from species_analysis module
    haplotype_taxonomy_csv : Path
        Haplotype taxonomy CSV with species assignments
    output_dir : Path
        Output directory for species-level divergence analysis

    Returns
    -------
    dict
        Results and file paths

    Notes
    -----
    Calculates:
    - Mean/max/min divergence within each species (across all haplotypes)
    - Mean/max/min divergence between species pairs
    (The barcoding gap is in divergence_analysis/label_based_gap.csv, from
    record labels; this summary uses haplotype majority labels and reports
    distances only.)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating species-level divergence analysis...")

    results = {}

    # Load species summary and haplotype taxonomy
    species_summary = pd.read_csv(species_summary_csv)
    haplotype_taxonomy = pd.read_csv(haplotype_taxonomy_csv)

    # Create a mapping from haplotype_id to species
    # The divergence matrix is indexed by display labels (haplotype_sp)
    haplotype_to_species = haplotype_species_map(haplotype_taxonomy)

    # Filter divergence matrix to only include haplotypes with species assignments
    valid_haplotypes = [h for h in divergence_matrix.index if h in haplotype_to_species]
    div_matrix_filtered = divergence_matrix.loc[valid_haplotypes, valid_haplotypes].copy()

    # Calculate species-level statistics
    species_divergence = []

    for species in species_summary['species']:
        # Get all haplotypes for this species
        species_haplotypes = [
            h for h in div_matrix_filtered.index
            if haplotype_to_species.get(h) == species
        ]

        if len(species_haplotypes) < 2:
            # Skip species with only one haplotype
            continue

        # Within-species divergence (all pairwise comparisons within species)
        within_values = []
        for i, h1 in enumerate(species_haplotypes):
            for h2 in species_haplotypes[i+1:]:
                div_val = div_matrix_filtered.loc[h1, h2]
                if not pd.isna(div_val):
                    within_values.append(div_val)

        # Between-species divergence (comparisons to all other species)
        between_values = []
        other_haplotypes = [
            h for h in div_matrix_filtered.index
            if h not in species_haplotypes
        ]

        for h1 in species_haplotypes:
            for h2 in other_haplotypes:
                div_val = div_matrix_filtered.loc[h1, h2]
                if not pd.isna(div_val):
                    between_values.append(div_val)

        # Calculate statistics
        if within_values:
            mean_within = np.mean(within_values)
            max_within = np.max(within_values)
            min_within = np.min(within_values)
        else:
            mean_within = max_within = min_within = np.nan

        if between_values:
            mean_between = np.mean(between_values)
            max_between = np.max(between_values)
            min_between = np.min(between_values)
        else:
            mean_between = max_between = min_between = np.nan

        species_divergence.append({
            'species': species,
            'n_haplotypes': len(species_haplotypes),
            'mean_within_species': mean_within,
            'max_within_species': max_within,
            'min_within_species': min_within,
            'mean_between_species': mean_between,
            'max_between_species': max_between,
            'min_between_species': min_between,
        })

    # Create dataframe and save
    if species_divergence:
        species_div_df = pd.DataFrame(species_divergence)
        species_div_df = species_div_df.sort_values('species')

        div_summary_path = output_dir / "species_divergence_summary.csv"
        species_div_df.to_csv(div_summary_path, index=False)
        results['species_divergence_summary'] = div_summary_path
        logger.info(f"  Saved species divergence summary: {div_summary_path}")

        # Store for console output
        results['species_div_stats'] = species_div_df.to_dict('records')
    else:
        logger.warning("  No species with multiple haplotypes found for divergence analysis")

    return results
