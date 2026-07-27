"""
Haplotype Discovery and Dereplication Module

This module implements a haplotype-first workflow for COI barcode analysis,
identifying Exact Sequence Variants (ESVs) from aligned regions and
flagging potentially suspect haplotypes based on distance and ORF quality.

The haplotype-first workflow:
1. Extract COI sequences from BOLD TSV file
2. Perform multiple sequence alignment using MAFFT
3. Mask gappy columns (>50% gaps) to remove insertion artifacts
4. Extract core shared region (positions covered by ≥80% of sequences)
5. Identify unique haplotypes (exact sequence variants) from core region
6. Assign samples to haplotypes via direct mapping (no re-calculation)
7. Calculate pairwise distances between haplotypes
8. Flag suspect haplotypes (singletons distant from nearest neighbor, ORF failures)
9. Generate haplotype FASTA with diagnostic information

Note: Step 6 uses direct mapping from ESV discovery rather than identity-based
re-assignment. This avoids false "tie" detections that occur when ESV haplotypes
differ by <0.3% (typical tie threshold), which is common with 1-2 bp differences.

Some Definitions for Context:
- **Haplotype (ESV)**: Exact Sequence Variant - unique sequence in core region
  Unlike clustering-based approaches, ESVs preserve all biological variation
  without arbitrary distance thresholds. Porter & Hajibabaei (2020) recommend ESVs for high-resolution COI barcoding.

- **Core Region**: The portion of alignment covered by most sequences (default: 80% of sequences have valid bases at that position)
  This addresses variable COI sequence lengths in BOLD (150-1550 bp) by focusing
  on the universally sequenced region. Gappy columns (>50% gaps) are masked before
  core region extraction to prevent insertion artifacts from fragmenting coverage.

- **Suspect Haplotype Flagging**: Singletons are flagged as suspect if:
  * Distance to nearest neighbor > threshold (default: 5%)
  * OR ORF validation failed (low coverage, internal stops)
  * OR suspect COI sequence (contamination indicator)

- **Distance Calculation**: p-distance on ungapped core region alignment
  Ignores gaps and N bases to focus on actual sequence differences.


Dependencies:
- MAFFT v7+ for multiple sequence alignment
- numpy for distance calculations
- Biopython for sequence handling

Example Usage:
    >>> from boldgenotyper.dereplication import identify_haplotypes
    >>> haplotypes = identify_haplotypes(
    ...     tsv_path="Sphyrna_lewini.tsv",
    ...     fasta_path="Sphyrna_lewini.fasta",
    ...     output_dir="results/",
    ...     min_core_coverage=0.8
    ... )

References:
    Porter & Hajibabaei (2020). Over 2.5 million COI sequences in GenBank
    and growing. PLOS ONE 15(9): e0238765.

Author: Steph Smith (symbioseas@outlook.com)
"""

from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
import logging
import subprocess
import shutil
import os
import re
from collections import Counter

import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# Configure logging
logger = logging.getLogger(__name__)


class DereplicationError(Exception):
    """Base exception for dereplication errors."""
    pass


class AlignmentError(DereplicationError):
    """Error during sequence alignment."""
    pass



def run_mafft_alignment(
    input_fasta: str,
    output_fasta: str,
    mafft_options: Optional[List[str]] = None
) -> None:
    """
    Perform multiple sequence alignment using MAFFT.

    Parameters
    ----------
    input_fasta : str
        Path to input FASTA file with unaligned sequences
    output_fasta : str
        Path to output aligned FASTA file
    mafft_options : List[str], optional
        Additional MAFFT command line options (default: ["--auto"])

    Raises
    ------
    AlignmentError
        If MAFFT is not found or alignment fails

    Examples
    --------
    >>> run_mafft_alignment("sequences.fasta", "aligned.fasta")
    >>> run_mafft_alignment("sequences.fasta", "aligned.fasta",
    ...                     mafft_options=["--maxiterate", "1000"])
    """
    # Check if MAFFT is available
    if shutil.which('mafft') is None:
        raise AlignmentError(
            "MAFFT not found in PATH. Please install MAFFT:\n"
            "  - macOS: brew install mafft\n"
            "  - Ubuntu/Debian: sudo apt-get install mafft\n"
            "  - conda: conda install -c bioconda mafft"
        )

    # Set default options
    if mafft_options is None:
        mafft_options = ["--auto"]

    # Build command
    cmd = ["mafft"] + mafft_options + [input_fasta]

    logger.info(f"Running MAFFT alignment: {' '.join(cmd)}")

    try:
        with open(output_fasta, 'w') as out_handle:
            result = subprocess.run(
                cmd,
                stdout=out_handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        logger.info(f"MAFFT alignment completed: {output_fasta}")

    except subprocess.CalledProcessError as e:
        error_msg = f"MAFFT alignment failed:\n{e.stderr}"
        logger.error(error_msg)
        raise AlignmentError(error_msg) from e
    except OSError as e:
        error_msg = f"Failed to write alignment output: {e}"
        logger.error(error_msg)
        raise AlignmentError(error_msg) from e


# ============================================================================
# Haplotype-First Workflow Functions
# ============================================================================

def identify_unique_haplotypes(
    core_alignment: List[str],
    headers: List[str],
    orf_validation_results: Optional[pd.DataFrame] = None
) -> Tuple[List[SeqRecord], pd.DataFrame]:
    """
    Identify unique haplotypes (exact sequence variants) from core region alignment.

    Parameters
    ----------
    core_alignment : List[str]
        Aligned core-region sequences (all same length)
    headers : List[str]
        Sequence headers (processids)
    orf_validation_results : Optional[pd.DataFrame]
        DataFrame with ORF validation results (processid, orf_valid, etc.)

    Returns
    -------
    Tuple[List[SeqRecord], pd.DataFrame]
        - haplotype_sequences: List of unique haplotype SeqRecords
        - haplotype_table: DataFrame mapping samples to haplotypes

    Notes
    -----
    Unlike clustering, this identifies exact sequence variants (ESVs) without
    any distance threshold. Each unique sequence in the core region defines
    a separate haplotype.
    """
    from collections import defaultdict

    # Group sequences by unique core sequence
    # Store both ungapped (for grouping) and aligned (for distance calculation)
    haplotype_groups = defaultdict(list)
    aligned_seqs = {}  # ungapped_seq -> aligned_seq

    for header, seq in zip(headers, core_alignment):
        # Ungapped sequence for haplotype ID
        ungapped = seq.replace('-', '')
        haplotype_groups[ungapped].append(header)
        # Store aligned version (first occurrence)
        if ungapped not in aligned_seqs:
            aligned_seqs[ungapped] = seq

    # Create haplotype records
    haplotype_records = []
    haplotype_mapping = []

    haplotype_id = 1
    for ungapped_seq, members in sorted(haplotype_groups.items(), key=lambda x: -len(x[1])):
        n_members = len(members)

        # Create haplotype ID (simplified format: h1_n8)
        hap_id = f"h{haplotype_id}_n{n_members}"

        # Create SeqRecord with ALIGNED sequence (for distance calculation)
        # The aligned version is needed for pairwise distance matrix
        seq_record = SeqRecord(
            Seq(aligned_seqs[ungapped_seq]),
            id=hap_id,
            description=""
        )
        haplotype_records.append(seq_record)

        # Record mapping
        for processid in members:
            haplotype_mapping.append({
                'processid': processid,
                'haplotype_id': hap_id,
                'haplotype_number': haplotype_id,
                'n_members': n_members,
                'is_singleton': n_members == 1
            })

        haplotype_id += 1

    # Convert mapping to DataFrame
    haplotype_df = pd.DataFrame(haplotype_mapping)

    # Merge with ORF validation if provided and contains non-null data
    if orf_validation_results is not None and len(orf_validation_results) > 0:
        orf_cols = ['processid', 'orf_valid', 'orf_coverage', 'internal_stops']
        available_cols = [c for c in orf_cols if c in orf_validation_results.columns]
        # Only merge if we have ORF data columns beyond processid
        if len(available_cols) > 1:
            orf_data = orf_validation_results[available_cols]
            # Only include columns that have non-null data
            cols_with_data = ['processid']
            for col in available_cols[1:]:
                if orf_data[col].notna().any():
                    cols_with_data.append(col)
            if len(cols_with_data) > 1:
                haplotype_df = haplotype_df.merge(
                    orf_data[cols_with_data],
                    on='processid',
                    how='left'
                )

    logger.info(f"Identified {len(haplotype_records)} unique haplotypes from {len(headers)} sequences")
    n_singletons = sum(1 for h in haplotype_records if h.id.endswith("_n1"))
    logger.info(f"  Singletons: {n_singletons} ({n_singletons/len(haplotype_records)*100:.1f}%)")

    return haplotype_records, haplotype_df


def _compute_distance(seq1: str, seq2: str) -> float:
    """
    Compute distance between two aligned sequences.

    Distance = 1 - (matches / valid_sites)
    where valid_sites are positions with A/C/G/T in both sequences.

    Parameters
    ----------
    seq1 : str
        First aligned sequence (uppercase)
    seq2 : str
        Second aligned sequence (uppercase)

    Returns
    -------
    float
        Distance value between 0.0 and 1.0
        Returns 1.0 if no valid comparison sites exist
    """
    matches = 0
    valid_sites = 0

    for base1, base2 in zip(seq1, seq2):
        if base1 in "ACGT" and base2 in "ACGT":
            valid_sites += 1
            if base1 == base2:
                matches += 1

    if valid_sites == 0:
        return 1.0

    return 1.0 - (matches / valid_sites)


def calculate_pairwise_distances(alignment: List[SeqRecord]) -> np.ndarray:
    """
    Calculate pairwise sequence distances ignoring gaps and Ns.

    Computes a condensed distance matrix suitable for
    scipy.cluster.hierarchy.linkage(). Distances are calculated as:
    distance = 1 - (matching_bases / total_valid_bases)

    Only positions with valid bases (A, C, G, T) in both sequences are
    considered. Gaps ('-') and ambiguous bases ('N') are ignored.

    Parameters
    ----------
    alignment : List[SeqRecord]
        Multiple sequence alignment (all sequences must have same length)

    Returns
    -------
    np.ndarray
        Condensed distance matrix (1D array of length n*(n-1)/2)
        Suitable for scipy.cluster.hierarchy.linkage()

    Raises
    ------
    ValueError
        If alignment is empty or sequences have different lengths
    """
    if not alignment:
        raise ValueError("Alignment is empty")

    n_seqs = len(alignment)
    if n_seqs < 2:
        raise ValueError("Need at least 2 sequences for distance calculation")

    aln_length = len(alignment[0].seq)
    if not all(len(rec.seq) == aln_length for rec in alignment):
        raise ValueError("All sequences in alignment must have the same length")

    logger.info(f"Calculating pairwise distances for {n_seqs} sequences")

    seqs = [str(rec.seq).upper() for rec in alignment]

    dist_vec = []
    total_comparisons = n_seqs * (n_seqs - 1) // 2
    comparison_count = 0

    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            distance = _compute_distance(seqs[i], seqs[j])
            dist_vec.append(distance)

            comparison_count += 1
            if comparison_count % 10000 == 0:
                logger.debug(
                    f"Distance calculation progress: {comparison_count}/{total_comparisons} "
                    f"({100*comparison_count/total_comparisons:.1f}%)"
                )

    logger.info(
        f"Distance calculation complete: {len(dist_vec)} pairwise distances"
    )

    return np.array(dist_vec)


def filter_error_singletons(
    haplotype_records: List[SeqRecord],
    haplotype_mapping: pd.DataFrame,
    min_singleton_distance: float = 0.005
) -> Tuple[List[SeqRecord], pd.DataFrame, Dict[str, int]]:
    """
    Filter likely sequencing/PCR error-derived singleton haplotypes.

    Removes singleton haplotypes that differ by ≤min_singleton_distance from their
    nearest neighbor, as these are statistically consistent with sequencing/PCR errors
    rather than biological variation.

    Parameters
    ----------
    haplotype_records : List[SeqRecord]
        Unique haplotype sequences (aligned core regions)
    haplotype_mapping : pd.DataFrame
        DataFrame mapping samples to haplotypes
    min_singleton_distance : float, optional
        Minimum divergence threshold for retaining singletons (default: 0.005 = 0.5%)
        Set to 0.0 to disable filtering (retain all singletons)

    Returns
    -------
    Tuple[List[SeqRecord], pd.DataFrame, Dict[str, int]]
        - Filtered haplotype records
        - Filtered haplotype mapping (with removed samples excluded)
        - Statistics dictionary with filtering counts

    Notes
    -----
    **Rationale for Error Filtering:**

    Sequencing and PCR errors create false singleton haplotypes that differ by 1-3 bp
    from true haplotypes:

    - **Sequencing errors:**
      - Illumina: ~0.1-0.3% error rate (Schirmer et al. 2015)
      - Sanger: ~0.1% error rate
      - In 650 bp sequence: ~1-2 errors per read

    - **PCR errors:**
      - Taq polymerase: ~0.01-0.1% error rate per bp
      - Accumulated across PCR cycles

    A singleton differing by 1-2 bp (0.15-0.3%) from its nearest neighbor is
    statistically more likely to be an error than a genuine rare variant.

    **Empirical validation:**
    Analysis of BOLD Sphyrnidae data showed ~85% of singletons fell within
    0.3% divergence, with median divergence of 0.0%, strongly suggesting
    error-driven inflation.

    **Default threshold (0.5%):**
    - Removes ~85% of singletons (likely errors)
    - Preserves singletons >0.5% divergent (rare variants, geographic structure)
    - Conservative balance between error removal and diversity preservation

    References
    ----------
    Schirmer et al. (2015). Insight into biases and sequencing errors for amplicon
    sequencing with the Illumina MiSeq platform. Nucleic Acids Research, 43(6), e37.

    Edgar (2013). UPARSE: highly accurate OTU sequences from microbial amplicon reads.
    Nature Methods, 10(10), 996-998.
    """
    if min_singleton_distance <= 0.0:
        # Filtering disabled
        return haplotype_records, haplotype_mapping, {
            'n_total': len(haplotype_records),
            'n_singletons': sum(1 for r in haplotype_records if '_n1' in r.id),
            'n_filtered': 0,
            'n_retained': len(haplotype_records)
        }

    # Calculate pairwise distances
    distance_matrix_condensed = calculate_pairwise_distances(haplotype_records)
    from scipy.spatial.distance import squareform
    distance_matrix = squareform(distance_matrix_condensed)

    # Identify singletons to keep
    records_to_keep = []
    haplotype_ids_to_keep = []

    n_singletons = 0
    n_filtered = 0

    for i, record in enumerate(haplotype_records):
        hap_id = record.id
        n_members = int(re.search(r'_n(\d+)', hap_id).group(1))
        is_singleton = n_members == 1

        if is_singleton:
            n_singletons += 1

            # Find nearest neighbor distance
            distances = distance_matrix[i, :]
            distances_to_others = [d for j, d in enumerate(distances) if j != i]

            if len(distances_to_others) > 0:
                nearest_neighbor_dist = min(distances_to_others)
            else:
                nearest_neighbor_dist = np.nan

            # Filter if below threshold
            if nearest_neighbor_dist <= min_singleton_distance:
                n_filtered += 1
                continue  # Don't keep this singleton

        # Keep this haplotype (multi-member or distant singleton)
        records_to_keep.append(record)
        haplotype_ids_to_keep.append(hap_id)

    # Filter mapping to exclude samples from removed singletons
    mapping_filtered = haplotype_mapping[
        haplotype_mapping['haplotype_id'].isin(haplotype_ids_to_keep)
    ].copy()

    # Collect statistics
    stats = {
        'n_total': len(haplotype_records),
        'n_singletons': n_singletons,
        'n_filtered': n_filtered,
        'n_retained': len(records_to_keep),
        'filter_rate': n_filtered / n_singletons if n_singletons > 0 else 0.0
    }

    logger.info(f"Singleton error filtering (threshold: >{min_singleton_distance*100:.1f}%):")
    logger.info(f"  Total haplotypes: {stats['n_total']}")
    logger.info(f"  Singletons: {stats['n_singletons']} ({stats['n_singletons']/stats['n_total']*100:.1f}%)")
    logger.info(f"  Filtered (likely errors): {stats['n_filtered']} ({stats['filter_rate']*100:.1f}% of singletons)")
    logger.info(f"  Retained haplotypes: {stats['n_retained']} ({stats['n_retained']/stats['n_total']*100:.1f}%)")

    return records_to_keep, mapping_filtered, stats


def flag_suspect_haplotypes(
    haplotype_records: List[SeqRecord],
    haplotype_df: pd.DataFrame,
    max_singleton_distance: float = 0.05,
    distance_outlier_min_abs: float = 0.10,
    distance_outlier_quantile: float = 0.95
) -> pd.DataFrame:
    """
    Flag suspect haplotypes based on distance and ORF quality.

    NOTE: This function flags potentially problematic haplotypes but does NOT filter them.
    Error-derived singletons should be filtered using filter_error_singletons() BEFORE
    calling this function.

    Suspect haplotypes are flagged if they are:
    1. Singletons with distance to nearest neighbor > threshold
    2. Any haplotype with distance > outlier threshold
    3. Haplotypes with failed ORF validation (if available)

    Parameters
    ----------
    haplotype_records : List[SeqRecord]
        Unique haplotype sequences
    haplotype_df : pd.DataFrame
        DataFrame mapping samples to haplotypes (with ORF validation if available)
    max_singleton_distance : float, optional
        Distance threshold for flagging singleton haplotypes (default: 0.05 = 5%)
    distance_outlier_min_abs : float, optional
        Absolute minimum distance for outlier flagging (default: 0.10 = 10%)
    distance_outlier_quantile : float, optional
        Quantile threshold for outlier detection (default: 0.95)

    Returns
    -------
    pd.DataFrame
        Haplotype-level statistics with suspect flags

    Notes
    -----
    This implements quality control recommended by Porter & Hajibabaei (2020)
    to flag likely sequencing errors, contamination, or other artifacts while
    preserving genuine biological variation.

    This flagging step is SEPARATE from and AFTER error singleton filtering:
    - filter_error_singletons(): Removes ≤0.5% divergent singletons (errors)
    - flag_suspect_haplotypes(): Flags >5% divergent singletons (contamination)
    """
    # Calculate pairwise distances between haplotypes
    distance_matrix_condensed = calculate_pairwise_distances(haplotype_records)

    # Convert condensed distance matrix to square form for easier indexing
    from scipy.spatial.distance import squareform
    distance_matrix = squareform(distance_matrix_condensed)

    # For each haplotype, find nearest neighbor distance
    haplotype_stats = []

    for i, record in enumerate(haplotype_records):
        # Get distances to all other haplotypes (exclude self)
        distances = distance_matrix[i, :]
        distances_to_others = [d for j, d in enumerate(distances) if j != i]

        if len(distances_to_others) > 0:
            nearest_neighbor_dist = min(distances_to_others)
        else:
            nearest_neighbor_dist = np.nan

        # Extract haplotype info from description
        hap_id = record.id
        n_members = int(re.search(r'_n(\d+)', hap_id).group(1))
        is_singleton = n_members == 1

        haplotype_stats.append({
            'haplotype_id': hap_id,
            'n_members': n_members,
            'is_singleton': is_singleton,
            'nearest_neighbor_distance': nearest_neighbor_dist
        })

    haplotype_stats_df = pd.DataFrame(haplotype_stats)

    # Calculate outlier threshold
    if len(haplotype_stats_df) > 0 and not haplotype_stats_df['nearest_neighbor_distance'].isna().all():
        distance_outlier_threshold = haplotype_stats_df['nearest_neighbor_distance'].quantile(
            distance_outlier_quantile
        )
    else:
        distance_outlier_threshold = distance_outlier_min_abs

    logger.info(f"Distance-based flagging thresholds:")
    logger.info(f"  Singleton distance threshold: {max_singleton_distance:.3f}")
    logger.info(f"  Outlier absolute threshold: {distance_outlier_min_abs:.3f}")
    logger.info(f"  Outlier quantile threshold ({distance_outlier_quantile:.0%}): {distance_outlier_threshold:.3f}")

    # Flag suspect haplotypes
    def flag_haplotype(row):
        reasons = []

        # Singleton with high distance
        if row['is_singleton'] and row['nearest_neighbor_distance'] > max_singleton_distance:
            reasons.append(f"Singleton distant from nearest neighbor ({row['nearest_neighbor_distance']:.3f})")

        # Outlier (any haplotype with very high distance)
        if (row['nearest_neighbor_distance'] >= distance_outlier_min_abs and
            row['nearest_neighbor_distance'] >= distance_outlier_threshold):
            reasons.append(f"Distance outlier ({row['nearest_neighbor_distance']:.3f})")

        return '; '.join(reasons) if reasons else ''

    haplotype_stats_df['suspect_distance_reason'] = haplotype_stats_df.apply(flag_haplotype, axis=1)
    haplotype_stats_df['suspect_distance'] = haplotype_stats_df['suspect_distance_reason'] != ''

    # Check for ORF failures (if ORF validation data available)
    if 'orf_valid' in haplotype_df.columns:
        # For each haplotype, check if ANY member has ORF failure
        orf_failures_by_haplotype = haplotype_df.groupby('haplotype_id').agg({
            'orf_valid': lambda x: (x == False).any() if x.notna().any() else False,  # True if any member has ORF failure
            'orf_coverage': 'min',  # Worst ORF coverage
            'internal_stops': 'max'  # Worst internal stops
        }).reset_index()
        orf_failures_by_haplotype.columns = ['haplotype_id', 'has_orf_failure', 'min_orf_coverage', 'max_internal_stops']

        haplotype_stats_df = haplotype_stats_df.merge(orf_failures_by_haplotype, on='haplotype_id', how='left')

        # Add ORF failure to suspect flags
        haplotype_stats_df['suspect_orf_reason'] = haplotype_stats_df.apply(
            lambda row: f"ORF failure (coverage={row['min_orf_coverage']:.2f}, stops={row['max_internal_stops']})"
            if row.get('has_orf_failure', False) else '',
            axis=1
        )
    else:
        haplotype_stats_df['has_orf_failure'] = False
        haplotype_stats_df['suspect_orf_reason'] = ''

    # Combined suspect flag
    haplotype_stats_df['is_suspect'] = (
        haplotype_stats_df['suspect_distance'] |
        haplotype_stats_df.get('has_orf_failure', False)
    )

    # Combined reasons
    def combine_reasons(row):
        reasons = []
        if row['suspect_distance_reason']:
            reasons.append(row['suspect_distance_reason'])
        if row.get('suspect_orf_reason', ''):
            reasons.append(row['suspect_orf_reason'])
        return '; '.join(reasons)

    haplotype_stats_df['suspect_reasons'] = haplotype_stats_df.apply(combine_reasons, axis=1)

    # Log summary
    n_suspect = haplotype_stats_df['is_suspect'].sum()
    n_suspect_distance = haplotype_stats_df['suspect_distance'].sum()
    n_suspect_orf = haplotype_stats_df.get('has_orf_failure', pd.Series([False])).sum()

    logger.info(f"Suspect haplotype flagging:")
    logger.info(f"  Total haplotypes: {len(haplotype_stats_df)}")
    logger.info(f"  Suspect haplotypes: {n_suspect} ({n_suspect/len(haplotype_stats_df)*100:.1f}%)")
    logger.info(f"    - Distance-based: {n_suspect_distance}")
    if 'orf_valid' in haplotype_df.columns:
        logger.info(f"    - ORF-based: {n_suspect_orf}")

    return haplotype_stats_df


def merge_zero_divergence_haplotypes(
    haplotype_records: List[SeqRecord],
    haplotype_mapping: pd.DataFrame,
    distance_threshold: float = 0.0
) -> Tuple[List[SeqRecord], pd.DataFrame, Dict[str, int]]:
    """
    Merge haplotypes with zero or near-zero divergence.

    This function addresses cases where sequences are identical in overlapping
    regions but have different 5'/3' endpoints or gap patterns, resulting in
    separate haplotypes during discovery but 0% divergence in pairwise comparison.

    For DNA barcoding applications, sequences with 0% divergence in their overlap
    should be considered the same haplotype regardless of length differences, as
    they represent the same biological sequence with different coverage.

    Parameters
    ----------
    haplotype_records : List[SeqRecord]
        Haplotype sequences (aligned core regions)
    haplotype_mapping : pd.DataFrame
        DataFrame mapping samples to haplotypes
    distance_threshold : float, optional
        Maximum divergence for merging (default: 0.0 = only identical sequences)

    Returns
    -------
    Tuple[List[SeqRecord], pd.DataFrame, Dict[str, int]]
        - Merged haplotype records
        - Updated haplotype mapping
        - Statistics dictionary

    Notes
    -----
    **Merging Strategy:**

    When multiple haplotypes have 0 divergence:
    1. Group them into connected components (transitively identical)
    2. Keep the haplotype with the most members as the representative
    3. Reassign all samples from merged haplotypes to the representative
    4. Update haplotype IDs and member counts

    **Example:**

    Before merging:
    - haplotype_h1_n100: 100 samples, 591 bp
    - haplotype_h5_n10: 10 samples, 587 bp (identical in 587 bp overlap)

    After merging:
    - haplotype_h1_n110: 110 samples (combined)

    **Rationale:**

    In COI barcoding, primer variation and sequence quality result in variable
    5'/3' endpoints. Sequences that are 100% identical in their overlap represent
    the same biological haplotype and should not be split based solely on
    coverage differences. This is especially important for:
    - Low-frequency haplotypes (avoid artificial splitting)
    - Geographic structure analysis (avoid artificial diversity inflation)
    - Species delimitation (avoid over-splitting)

    Examples
    --------
    >>> records, mapping, stats = merge_zero_divergence_haplotypes(
    ...     haplotype_records, haplotype_mapping
    ... )
    >>> print(f"Merged {stats['n_merged']} haplotypes")
    """
    if len(haplotype_records) <= 1:
        return haplotype_records, haplotype_mapping, {
            'n_input': len(haplotype_records),
            'n_merged': 0,
            'n_output': len(haplotype_records)
        }

    logger.info(f"Checking for zero-divergence haplotypes to merge...")

    # Calculate pairwise distances
    distance_matrix_condensed = calculate_pairwise_distances(haplotype_records)
    from scipy.spatial.distance import squareform
    distance_matrix = squareform(distance_matrix_condensed)

    # Find pairs with divergence <= threshold
    n = len(haplotype_records)
    haplotype_ids = [rec.id for rec in haplotype_records]

    # Build merge groups using union-find
    from collections import defaultdict
    parent = {hap_id: hap_id for hap_id in haplotype_ids}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Find all pairs with divergence <= threshold
    merge_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i, j] <= distance_threshold:
                merge_pairs.append((haplotype_ids[i], haplotype_ids[j]))
                union(haplotype_ids[i], haplotype_ids[j])

    if not merge_pairs:
        logger.info(f"  No zero-divergence haplotypes found - no merging needed")
        return haplotype_records, haplotype_mapping, {
            'n_input': len(haplotype_records),
            'n_merged': 0,
            'n_output': len(haplotype_records)
        }

    # Group haplotypes by their root (connected components)
    groups = defaultdict(list)
    for hap_id in haplotype_ids:
        root = find(hap_id)
        groups[root].append(hap_id)

    # For each group, select representative (largest member count)
    haplotype_merges = {}  # old_id -> new_id
    representatives = {}   # new_id -> SeqRecord

    for root, group_members in groups.items():
        if len(group_members) == 1:
            # No merging needed for this group
            hap_id = group_members[0]
            record = next(r for r in haplotype_records if r.id == hap_id)
            representatives[hap_id] = record
        else:
            # Select representative: largest member count, or first if tied
            member_counts = {}
            for hap_id in group_members:
                n_members = int(re.search(r'_n(\d+)', hap_id).group(1))
                member_counts[hap_id] = n_members

            representative_id = max(group_members, key=lambda x: member_counts[x])
            representative_record = next(r for r in haplotype_records if r.id == representative_id)

            # Map all group members to representative
            for hap_id in group_members:
                if hap_id != representative_id:
                    haplotype_merges[hap_id] = representative_id

            representatives[representative_id] = representative_record

            logger.info(f"  Merging {len(group_members)} haplotypes into {representative_id}:")
            for hap_id in sorted(group_members):
                if hap_id == representative_id:
                    logger.info(f"    ✓ {hap_id} (representative, n={member_counts[hap_id]})")
                else:
                    logger.info(f"    → {hap_id} (n={member_counts[hap_id]}) merged into {representative_id}")

    # Update mapping
    if haplotype_merges:
        haplotype_mapping = haplotype_mapping.copy()
        haplotype_mapping['haplotype_id'] = haplotype_mapping['haplotype_id'].replace(haplotype_merges)

    # Recalculate member counts and update records
    member_counts = haplotype_mapping.groupby('haplotype_id').size()

    merged_records = []
    for hap_id, record in representatives.items():
        n_members = member_counts.get(hap_id, 0)
        hap_number = int(re.search(r'h(\d+)', hap_id).group(1))

        updated_record = SeqRecord(
            record.seq,
            id=hap_id,
            description=f"Haplotype {hap_number} with {n_members} sequences"
        )
        merged_records.append(updated_record)

    logger.info(f"  Merged {len(haplotype_merges)} haplotypes with zero divergence")
    logger.info(f"  Output: {len(merged_records)} unique haplotypes")

    stats = {
        'n_input': len(haplotype_records),
        'n_merged': len(haplotype_merges),
        'n_output': len(merged_records),
        'n_groups': len([g for g in groups.values() if len(g) > 1])
    }

    return merged_records, haplotype_mapping, stats


def identify_haplotypes(
    tsv_path: Union[str, Path],
    fasta_path: Union[str, Path],
    output_dir: Union[str, Path],
    min_core_coverage: float = 0.8,
    gap_mask_threshold: float = 0.5,
    min_core_length: int = 200,
    min_singleton_distance: float = 0.005,
    max_singleton_distance: float = 0.05,
    orf_validation_df: Optional[pd.DataFrame] = None,
    cleanup_intermediates: bool = False
) -> Tuple[List[SeqRecord], pd.DataFrame, pd.DataFrame]:
    """
    Main haplotype-first workflow: identify ESVs from core region.
    Includes two-stage singleton quality control to remove sequencing/PCR errors
    while preserving biological variation.

    Parameters
    ----------
    tsv_path : Union[str, Path]
        Path to BOLD TSV file
    fasta_path : Union[str, Path]
        Path to input FASTA file (orientation-normalized sequences)
    output_dir : Union[str, Path]
        Output directory for results
    min_core_coverage : float, optional
        Minimum sequence coverage for core region (default: 0.8 = 80%)
    gap_mask_threshold : float, optional
        Gap threshold for masking columns before core region extraction.
        Columns with >gap_mask_threshold gaps are masked (default: 0.5 = 50%)
    min_core_length : int, optional
        Minimum core region length (default: 200 bp)
    min_singleton_distance : float, optional
        Minimum divergence for retaining singletons (default: 0.005 = 0.5%)
        Singletons ≤0.5% divergent are filtered as likely errors.
        Set to 0.0 to disable filtering.
    max_singleton_distance : float, optional
        Distance threshold for flagging singletons (default: 0.05 = 5%)
        Singletons >5% divergent are flagged as potential contamination.
    orf_validation_df : Optional[pd.DataFrame]
        ORF validation results from quality control
    cleanup_intermediates : bool, optional
        Remove intermediate files (default: False)

    Returns
    -------
    Tuple[List[SeqRecord], pd.DataFrame, pd.DataFrame]
        - haplotype_records: Unique haplotype sequences (after filtering)
        - haplotype_mapping: Sample-to-haplotype mapping (after filtering)
        - haplotype_stats: Haplotype-level statistics with flags

    Notes
    -----
    **Two-Stage Singleton Quality Control:**

    1. **Error filtering (Step 4.5):**
       - Removes singletons ≤0.5% divergent (likely sequencing/PCR errors)
       - Applied BEFORE flagging
       - Default: ENABLED (min_singleton_distance=0.005)

    2. **Suspect flagging (Step 5):**
       - Flags singletons >5% divergent (likely contamination/misID)
       - Applied AFTER error filtering
       - Default: depends on HaplotypeConfig.flag_suspect_haplotypes

    **Workflow:**
    1. Align sequences with MAFFT
    2. Mask gappy columns (>gap_mask_threshold gaps) to remove insertion artifacts
    3. Extract core shared region (covered by ≥min_core_coverage of sequences)
    4. Identify unique haplotypes (ESVs)
    4.5. Filter error-derived singletons (NEW)
    5. Calculate distances and flag suspects
    6. Write output files
    """
    from boldgenotyper.utils import extract_core_region

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("Starting haplotype-first workflow")
    logger.info("=" * 70)

    # Step 1: Align sequences
    logger.info("Step 1: Aligning sequences with MAFFT...")
    aligned_fasta = output_path / f"{Path(fasta_path).stem}_aligned.fasta"
    run_mafft_alignment(str(fasta_path), str(aligned_fasta))

    # Step 2: Load alignment
    logger.info("Step 2: Loading alignment...")
    # Guard against empty FASTA (all sequences failed QC)
    if aligned_fasta.stat().st_size == 0 or not any(
        line.startswith(b">") for line in aligned_fasta.open("rb")
    ):
        raise ValueError(
            "No sequences remain after quality control — cannot proceed with haplotype discovery. "
            "All sequences failed ORF validation. If your organism is an invertebrate, try "
            "--genetic-code 5 (invertebrate mitochondrial). Check the ORF validation report in "
            "the intermediate/quality_control/ directory for details."
        )
    alignment = list(AlignIO.read(str(aligned_fasta), "fasta"))
    aligned_seqs = [str(record.seq) for record in alignment]
    headers = [record.id for record in alignment]

    logger.info(f"Loaded alignment: {len(aligned_seqs)} sequences, {len(aligned_seqs[0])} positions")

    # Step 3: Extract core region
    logger.info("Step 3: Extracting core shared region...")
    core_result = extract_core_region(
        aligned_seqs,
        headers,
        min_coverage=min_core_coverage,
        gap_threshold=gap_mask_threshold,
        min_core_length=min_core_length
    )

    if core_result is None:
        raise DereplicationError(
            "Core region extraction failed - region < 100 bp absolute minimum. "
            "This dataset may have insufficient overlap in sequenced regions. "
            "Consider checking sequence quality and alignment coverage."
        )

    core_seqs, core_headers = core_result
    logger.info(f"Core region: {len(core_seqs)} sequences, {len(core_seqs[0])} positions")

    # Step 4: Identify unique haplotypes
    logger.info("Step 4: Identifying unique haplotypes...")
    haplotype_records, haplotype_mapping = identify_unique_haplotypes(
        core_seqs,
        core_headers,
        orf_validation_results=orf_validation_df
    )

    # Step 4.5: Filter error-derived singletons
    logger.info("Step 4.5: Filtering error-derived singletons...")
    haplotype_records, haplotype_mapping, filter_stats = filter_error_singletons(
        haplotype_records,
        haplotype_mapping,
        min_singleton_distance=min_singleton_distance
    )

    # Step 5: Flag suspect haplotypes
    logger.info("Step 5: Flagging suspect haplotypes...")
    haplotype_stats = flag_suspect_haplotypes(
        haplotype_records,
        haplotype_mapping,
        max_singleton_distance=max_singleton_distance
    )

    # Step 5.5: Merge zero-divergence haplotypes
    logger.info("Step 5.5: Merging zero-divergence haplotypes...")
    haplotype_records, haplotype_mapping, merge_stats = merge_zero_divergence_haplotypes(
        haplotype_records,
        haplotype_mapping,
        distance_threshold=0.0
    )

    # Update haplotype stats after merging
    if merge_stats['n_merged'] > 0:
        # Recalculate stats for merged haplotypes
        logger.info("  Recalculating haplotype statistics after merging...")
        haplotype_stats = flag_suspect_haplotypes(
            haplotype_records,
            haplotype_mapping,
            max_singleton_distance=max_singleton_distance
        )

    # Step 6: Write haplotype FASTA (ungapped for downstream use)
    haplotype_fasta = output_path / f"{Path(tsv_path).stem}_haplotypes.fasta"

    # Create ungapped versions for writing to FASTA
    ungapped_records = []
    for rec in haplotype_records:
        ungapped_seq = str(rec.seq).replace('-', '')
        ungapped_rec = SeqRecord(
            Seq(ungapped_seq),
            id=rec.id,
            description=rec.description
        )
        ungapped_records.append(ungapped_rec)

    SeqIO.write(ungapped_records, haplotype_fasta, "fasta")
    logger.info(f"Wrote {len(ungapped_records)} haplotypes to {haplotype_fasta}")

    # Step 7: Write haplotype tables
    haplotype_mapping_file = output_path / f"{Path(tsv_path).stem}_haplotype_mapping.csv"
    haplotype_mapping.to_csv(haplotype_mapping_file, index=False)
    logger.info(f"Wrote haplotype mapping to {haplotype_mapping_file}")

    haplotype_stats_file = output_path / f"{Path(tsv_path).stem}_haplotype_stats.csv"
    # Remove ORF columns if they are all empty/null
    orf_cols = ['min_orf_coverage', 'max_internal_stops', 'suspect_orf_reason', 'has_orf_failure']
    cols_to_drop = [c for c in orf_cols if c in haplotype_stats.columns and haplotype_stats[c].isna().all()]
    # Also drop if all values are empty strings or False
    for c in orf_cols:
        if c in haplotype_stats.columns and c not in cols_to_drop:
            if haplotype_stats[c].dtype == 'object' and (haplotype_stats[c] == '').all():
                cols_to_drop.append(c)
            elif haplotype_stats[c].dtype == 'bool' and not haplotype_stats[c].any():
                cols_to_drop.append(c)
    haplotype_stats_clean = haplotype_stats.drop(columns=cols_to_drop, errors='ignore')
    haplotype_stats_clean.to_csv(haplotype_stats_file, index=False)
    logger.info(f"Wrote haplotype statistics to {haplotype_stats_file}")

    # Cleanup
    if cleanup_intermediates and aligned_fasta.exists():
        aligned_fasta.unlink()
        logger.debug(f"Removed intermediate file: {aligned_fasta}")

    logger.info("=" * 70)
    logger.info("Haplotype-first workflow completed successfully")
    logger.info("=" * 70)

    return haplotype_records, haplotype_mapping, haplotype_stats
