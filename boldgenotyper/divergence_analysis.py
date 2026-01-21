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

Author: Steph Smith (steph.smith@unc.edu)
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

        logger.info(f"  ✓ Aligned sequences saved to {aligned_path}")
        return aligned_path

    except subprocess.CalledProcessError as e:
        logger.warning(f"MAFFT alignment failed: {e}")
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

    # Create labels mapping consensus_group to species
    group_to_species = {}
    if 'consensus_group' in taxonomy_df.columns and 'assigned_sp' in taxonomy_df.columns:
        group_to_species = dict(zip(
            taxonomy_df['consensus_group'],
            taxonomy_df['assigned_sp']
        ))

    # Create row labels (consensus_group_species format)
    labels = []
    seq_ids = list(sequences.keys())

    for seq_id in seq_ids:
        # Extract consensus group from sequence ID
        # Assuming format like "consensus_c1" or just "c1"
        group = seq_id.replace('consensus_', '')
        species = group_to_species.get(group, group_to_species.get(seq_id, ''))

        if species:
            label = f"{group}_{species.replace(' ', '_')}"
        else:
            label = group

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

    # Parse labels to extract group and species
    def parse_label(label):
        """Extract consensus group and species from label."""
        parts = label.split('_')
        if len(parts) >= 2:
            group = parts[0]
            species = '_'.join(parts[1:])
        else:
            group = label
            species = ''
        return group, species

    group_species_map = {label: parse_label(label) for label in labels}

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

            # Categorize divergence
            if comparison_type == 'within_species':
                if divergence < 0.01:
                    interpretation = 'low_divergence_within_species'
                elif divergence < 0.02:
                    interpretation = 'moderate_divergence_within_species'
                else:
                    interpretation = 'high_divergence_within_species'
            else:
                if divergence < 0.03:
                    interpretation = 'low_interspecific'
                elif divergence < 0.05:
                    interpretation = 'typical_interspecific'
                else:
                    interpretation = 'high_interspecific'

            results.append({
                'comparison_type': comparison_type,
                'group1': group1,
                'group2': group2,
                'species1': sp1,
                'species2': sp2,
                'species_pair': species_pair,
                'divergence': divergence,
                'interpretation': interpretation
            })

    summary_df = pd.DataFrame(results)

    # Calculate aggregate statistics
    if len(summary_df) > 0:
        agg_summary = summary_df.groupby(
            ['comparison_type', 'species_pair']
        )['divergence'].agg(['mean', 'std', 'min', 'max', 'count'])

        agg_summary = agg_summary.reset_index()
        agg_summary.columns = [
            'comparison_type', 'species_pair', 'mean_divergence',
            'sd_divergence', 'min_divergence', 'max_divergence', 'n_comparisons'
        ]

        return agg_summary
    else:
        return pd.DataFrame()


def analyze_barcoding_gap(
    divergence_matrix: pd.DataFrame,
    taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """
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

    # Parse labels
    def parse_label(label):
        parts = label.split('_')
        if len(parts) >= 2:
            group = parts[0]
            species = '_'.join(parts[1:])
        else:
            group = label
            species = ''
        return group, species

    group_species_map = {label: parse_label(label) for label in labels}

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

    # Determine figure size based on matrix size
    n = len(mat)
    figsize = (max(9, n * 0.4), max(9, n * 0.4))

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
                    fontsize=8,
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
                    fontsize=7.5,
                    color=txt_color,
                )

    # Rotate tick labels for readability
    for tick in ax.get_xticklabels():
        tick.set_fontsize(9)
        tick.set_rotation(45)
        tick.set_ha("right")

    for tick in ax.get_yticklabels():
        tick.set_fontsize(9)

    # Add title
    plt.subplots_adjust(top=0.96)
    cg.fig.suptitle("Haplotype Divergence (Uncorrected P-Distance)", fontsize=11)

    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"  ✓ Divergence heatmap saved to {output_path}")


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

    # Parse labels
    def parse_label(label):
        parts = label.split('_')
        if len(parts) >= 2:
            group = parts[0]
            species = '_'.join(parts[1:])
        else:
            group = label
            species = ''
        return group, species

    group_species_map = {label: parse_label(label) for label in labels}

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

    logger.info(f"  ✓ Barcoding gap plot saved to {output_path}")


def generate_divergence_analysis(
    consensus_fasta: Path,
    taxonomy_csv: Path,
    output_dir: Path,
    organism: str
) -> Dict[str, Any]:
    """
    Complete divergence analysis workflow.

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

    # Save matrix
    matrix_path = output_dir / "pairwise_divergence_matrix.csv"
    divergence_matrix.to_csv(matrix_path)
    results['divergence_matrix'] = matrix_path
    logger.info(f"  ✓ Saved divergence matrix: {matrix_path}")

    # Generate divergence summary
    summary_df = generate_divergence_summary(divergence_matrix, taxonomy_df)
    if len(summary_df) > 0:
        summary_path = output_dir / "divergence_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        results['divergence_summary'] = summary_path
        logger.info(f"  ✓ Saved divergence summary: {summary_path}")

    # Analyze barcoding gap
    gap_df = analyze_barcoding_gap(divergence_matrix, taxonomy_df)
    if len(gap_df) > 0:
        gap_path = output_dir / "within_vs_between_species.csv"
        gap_df.to_csv(gap_path, index=False)
        results['barcoding_gap'] = gap_path
        logger.info(f"  ✓ Saved barcoding gap analysis: {gap_path}")

        # Store summary stats for console output
        results['gap_stats'] = gap_df.to_dict('records')

    # Create visualizations
    try:
        heatmap_path = output_dir / "divergence_heatmap.pdf"
        create_divergence_heatmap(divergence_matrix, heatmap_path)
        results['heatmap'] = heatmap_path
    except Exception as e:
        logger.warning(f"Could not create divergence heatmap: {e}")

    try:
        gap_plot_path = output_dir / "barcoding_gap.pdf"
        create_barcoding_gap_plot(divergence_matrix, taxonomy_df, gap_plot_path)
        results['gap_plot'] = gap_plot_path
    except Exception as e:
        logger.warning(f"Could not create barcoding gap plot: {e}")

    # Create README
    readme_path = output_dir / "README.md"
    readme_content = f"""# Divergence Analysis for {organism}

## Overview

This directory contains pairwise divergence analysis results for consensus sequences.

## Files

- **pairwise_divergence_matrix.csv**: Square matrix of uncorrected p-distances
- **divergence_summary.csv**: Summary statistics by comparison type
- **within_vs_between_species.csv**: Barcoding gap analysis per species
- **divergence_heatmap.pdf**: Visual heatmap of pairwise divergences
- **barcoding_gap.pdf**: Within vs between species divergence distributions

## Interpretation

### Divergence Categories

**Within-Species:**
- Low (<1%): Typical for most COI barcodes
- Moderate (1-2%): Structured populations or recent divergence
- High (>2%): Potential cryptic species or misidentification

**Between-Species:**
- Low (<3%): Recently diverged or hybridizing species
- Typical (3-5%): Standard interspecific divergence for COI
- High (>5%): Deeply divergent lineages

### Barcoding Gap

A "barcoding gap" exists when the minimum between-species divergence exceeds
the maximum within-species divergence. This indicates clear species boundaries
and supports species delimitation.

## Methodology

Divergences calculated as uncorrected p-distance:
```
distance = 1 - (matches / informative_positions)
```

Gaps and ambiguous bases (N) are excluded from calculations.
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info(f"  ✓ Divergence analysis complete")

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

    for species_data in results['gap_stats']:
        species = species_data['species']
        mean_within = species_data['mean_within']
        max_within = species_data['max_within']
        mean_between = species_data.get('mean_between', np.nan)
        min_between = species_data.get('min_between', np.nan)
        has_gap = species_data['barcoding_gap']

        print(f"\n{species}:")
        print(f"  Within-species divergence: {mean_within:.3f} (max: {max_within:.3f})")

        if not pd.isna(mean_between):
            print(f"  Between-species divergence: {mean_between:.3f} (min: {min_between:.3f})")

        # Warnings for high within-species divergence
        if mean_within > 0.02:
            print(f"  ⚠️  WARNING: Within-species divergence >{2:.1f}% suggests potential cryptic lineages")
            print("      Consider:")
            print("      1. Nuclear marker validation")
            print("      2. Morphological examination")
            print("      3. Biogeographic context")

        # Barcoding gap status
        if has_gap:
            gap_size = species_data['gap_size']
            print(f"  ✓ Barcoding gap present: {gap_size:.3f}")
        else:
            if not pd.isna(min_between):
                overlap = max_within - min_between
                print(f"  ✗ Barcoding gap absent (overlap: {overlap:.3f})")

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
    - Barcoding gap assessment for each species
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating species-level divergence analysis...")

    results = {}

    # Load species summary and haplotype taxonomy
    species_summary = pd.read_csv(species_summary_csv)
    haplotype_taxonomy = pd.read_csv(haplotype_taxonomy_csv)

    # Create a mapping from haplotype_id to species
    haplotype_to_species = dict(zip(
        haplotype_taxonomy['haplotype_id'],
        haplotype_taxonomy['assigned_sp']
    ))

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

        # Barcoding gap
        has_gap = False
        gap_size = np.nan
        if not pd.isna(max_within) and not pd.isna(min_between):
            has_gap = min_between > max_within
            gap_size = min_between - max_within if has_gap else np.nan

        species_divergence.append({
            'species': species,
            'n_haplotypes': len(species_haplotypes),
            'mean_within_species': mean_within,
            'max_within_species': max_within,
            'min_within_species': min_within,
            'mean_between_species': mean_between,
            'max_between_species': max_between,
            'min_between_species': min_between,
            'barcoding_gap_present': has_gap,
            'gap_size': gap_size
        })

    # Create dataframe and save
    if species_divergence:
        species_div_df = pd.DataFrame(species_divergence)
        species_div_df = species_div_df.sort_values('species')

        div_summary_path = output_dir / "species_divergence_summary.csv"
        species_div_df.to_csv(div_summary_path, index=False)
        results['species_divergence_summary'] = div_summary_path
        logger.info(f"  ✓ Saved species divergence summary: {div_summary_path}")

        # Store for console output
        results['species_div_stats'] = species_div_df.to_dict('records')
    else:
        logger.warning("  ⚠ No species with multiple haplotypes found for divergence analysis")

    return results
