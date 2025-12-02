"""
Haplotype Discovery and Dereplication Module

This module implements a haplotype-first workflow for COI barcode analysis,
identifying Exact Sequence Variants (ESVs) from aligned core regions and
flagging potentially suspect haplotypes based on distance and ORF quality.

The haplotype-first workflow:
1. Extract COI sequences from BOLD TSV file
2. Perform multiple sequence alignment using MAFFT
3. Extract core shared region (positions covered by ≥80% of sequences)
4. Apply gap masking to remove low-information columns
5. Identify unique haplotypes (exact sequence variants) from core region
6. Calculate pairwise distances between haplotypes
7. Flag suspect haplotypes (singletons distant from nearest neighbor, ORF failures)
8. Generate haplotype FASTA with diagnostic information

Key Concepts:
- **Haplotype (ESV)**: Exact Sequence Variant - unique sequence in core region
  Unlike clustering-based approaches, ESVs preserve all biological variation
  without arbitrary distance thresholds.

- **Core Region**: The portion of alignment covered by most sequences (default: 80%)
  This addresses variable COI sequence lengths in BOLD (150-1550 bp) by focusing
  on the universally sequenced region.

- **Suspect Haplotype Flagging**: Singletons are flagged as suspect if:
  * Distance to nearest neighbor > threshold (default: 5%)
  * OR ORF validation failed (low coverage, internal stops)
  * OR suspect COI sequence (contamination indicator)

- **Distance Calculation**: p-distance on ungapped core region alignment
  Ignores gaps and N bases to focus on actual sequence differences.

Workflow Comparison:
- OLD (clustering-first): align → cluster → generate consensus
- NEW (haplotype-first): align → extract core → identify ESVs → flag suspects

This implements the ESV approach recommended by Porter & Hajibabaei (2020) for
high-resolution COI barcoding while maintaining data quality through suspect flagging.

Dependencies:
- MAFFT v7+ for multiple sequence alignment
- trimAl for alignment trimming (optional)
- scipy for distance calculations
- Biopython for sequence handling

Example Usage:
    >>> from boldgenotyper.dereplication import identify_haplotypes
    >>> haplotypes = identify_haplotypes(
    ...     tsv_path="Sphyrna_lewini.tsv",
    ...     output_dir="results/",
    ...     min_core_coverage=0.8
    ... )

Legacy Support:
    The clustering-based workflow is still available via cluster_sequences()
    and generate_consensus() for backward compatibility.

References:
    Porter & Hajibabaei (2020). Over 2.5 million COI sequences in GenBank
    and growing. PLOS ONE 15(9): e0238765.

Author: Steph Smith (steph.smith@unc.edu)
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
from scipy.cluster.hierarchy import linkage, fcluster

# Configure logging
logger = logging.getLogger(__name__)


class DereplicationError(Exception):
    """Base exception for dereplication errors."""
    pass


class AlignmentError(DereplicationError):
    """Error during sequence alignment."""
    pass


class TrimmingError(DereplicationError):
    """Error during alignment trimming."""
    pass


class ClusteringError(DereplicationError):
    """Error during sequence clustering."""
    pass


def check_external_tools() -> Dict[str, bool]:
    """
    Check if required external tools (MAFFT, trimAl) are available.

    Returns
    -------
    Dict[str, bool]
        Dictionary mapping tool names to availability status

    Examples
    --------
    >>> tools = check_external_tools()
    >>> if not tools['mafft']:
    ...     print("MAFFT not found, please install it")
    """
    tools = {
        'mafft': shutil.which('mafft') is not None,
        'trimal': shutil.which('trimal') is not None
    }
    return tools


def filter_by_ungapped_length(
    alignment: List[SeqRecord],
    min_length: int
) -> List[SeqRecord]:
    """
    Filter alignment by ungapped sequence length.

    Removes sequences where the ungapped (non-gap) length is below the threshold.
    This prevents short sequence fragments from being included in clustering.

    Parameters
    ----------
    alignment : List[SeqRecord]
        List of aligned sequence records
    min_length : int
        Minimum ungapped sequence length to retain

    Returns
    -------
    List[SeqRecord]
        Filtered alignment with only sequences >= min_length (ungapped)

    Examples
    --------
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> seqs = [
    ...     SeqRecord(Seq("ATCG--ATCG"), id="seq1"),  # 8 bp ungapped
    ...     SeqRecord(Seq("AT----ATCG"), id="seq2"),  # 6 bp ungapped
    ... ]
    >>> filtered = filter_by_ungapped_length(seqs, min_length=7)
    >>> len(filtered)
    1
    """
    filtered = []
    for record in alignment:
        # Calculate ungapped length (remove gaps and ambiguous positions)
        ungapped_seq = str(record.seq).replace('-', '').replace('N', '')
        ungapped_length = len(ungapped_seq)

        if ungapped_length >= min_length:
            filtered.append(record)
        else:
            logger.debug(
                f"  Filtering {record.id}: ungapped length {ungapped_length} < {min_length}"
            )

    return filtered


def filter_consensus_by_length(
    consensus_records: Dict[str, SeqRecord],
    min_length_ratio: float
) -> Dict[str, SeqRecord]:
    """
    Filter consensus sequences by relative length.

    Removes consensus sequences that are suspiciously short compared to the median,
    as these likely represent sequence fragments rather than true genotypes.

    Parameters
    ----------
    consensus_records : Dict[str, SeqRecord]
        Dictionary of consensus sequence records (id -> SeqRecord)
    min_length_ratio : float
        Minimum length as fraction of median (e.g., 0.75 = 75% of median)

    Returns
    -------
    Dict[str, SeqRecord]
        Filtered consensus records

    Examples
    --------
    >>> consensus = {
    ...     'c1': SeqRecord(Seq("A" * 500), id="c1"),  # Full length
    ...     'c2': SeqRecord(Seq("A" * 490), id="c2"),  # Full length
    ...     'c3': SeqRecord(Seq("A" * 200), id="c3"),  # Fragment
    ... }
    >>> filtered = filter_consensus_by_length(consensus, min_length_ratio=0.75)
    >>> len(filtered)  # c3 removed (200 < 0.75 * 495)
    2
    """
    if len(consensus_records) <= 1:
        return consensus_records

    # Calculate ungapped lengths
    lengths = {}
    for cons_id, record in consensus_records.items():
        ungapped_seq = str(record.seq).replace('-', '')
        lengths[cons_id] = len(ungapped_seq)

    # Calculate median length
    median_length = np.median(list(lengths.values()))
    min_allowed_length = median_length * min_length_ratio

    logger.debug(f"  Consensus length stats: median={median_length:.0f} bp, "
                 f"minimum={min_allowed_length:.0f} bp (ratio={min_length_ratio})")

    # Filter by length
    filtered = {}
    for cons_id, record in consensus_records.items():
        if lengths[cons_id] >= min_allowed_length:
            filtered[cons_id] = record
        else:
            logger.info(
                f"  Filtering short consensus {cons_id}: "
                f"{lengths[cons_id]} bp < {min_allowed_length:.0f} bp "
                f"({min_length_ratio:.0%} of median)"
            )

    return filtered


def dereplicate_from_fasta(
    input_fasta: Union[str, Path],
    output_dir: Union[str, Path],
    threshold: float = 0.02,
    frequency_cutoff: float = 0.7,
    mafft_options: Optional[List[str]] = None,
    trimal_options: Optional[List[str]] = None,
    cleanup_intermediates: bool = False,
    organism_name: Optional[str] = None,
    min_post_trim_length: int = 300,
    min_consensus_length_ratio: float = 0.75
    ) -> Dict[str, SeqRecord]:
    input_fasta = Path(input_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tools = check_external_tools()
    if not tools['mafft']:
        raise DereplicationError("MAFFT not found. Please install MAFFT.")
    if not tools['trimal']:
        raise DereplicationError("trimAl not found. Please install trimAl.")

    if organism_name is None:
        organism_name = input_fasta.stem

    aligned_fasta   = output_dir / f"{organism_name}_aligned.fasta"
    trimmed_fasta   = output_dir / f"{organism_name}_trimmed.fasta"
    consensus_fasta = output_dir / f"{organism_name}_consensus.fasta"

    # Align
    run_mafft_alignment(str(input_fasta), str(aligned_fasta), mafft_options=mafft_options)
    # Trim
    run_trimal_trimming(str(aligned_fasta), str(trimmed_fasta), trimal_options=trimal_options)
    # Load trimmed alignment
    alignment = list(AlignIO.read(str(trimmed_fasta), "fasta"))

    # Stage 2: Post-trimming length filtering
    # Remove sequences that became too short after trimming
    if min_post_trim_length > 0:
        original_count = len(alignment)
        alignment = filter_by_ungapped_length(alignment, min_post_trim_length)
        filtered_count = original_count - len(alignment)
        if filtered_count > 0:
            logger.info(f"  Filtered {filtered_count} sequences shorter than {min_post_trim_length} bp after trimming")
            if len(alignment) == 0:
                raise DereplicationError(
                    f"All sequences were filtered out after trimming. "
                    f"Consider lowering min_post_trim_length (currently {min_post_trim_length})"
                )

    # Distances
    distances = calculate_pairwise_distances(alignment)
    # Cluster
    labels = cluster_sequences(distances, threshold=threshold)
    # Consensus
    clusters = {}
    for rec, cid in zip(alignment, labels):
        clusters.setdefault(cid, []).append(rec)
    consensus_records = {}
    consensus_metadata = []  # Track reference counts for metadata file
    for cid, seqs in sorted(clusters.items()):
        cons = generate_consensus(seqs, cid, frequency_cutoff=frequency_cutoff)
        consensus_records[cons.id] = cons
        # Store metadata: consensus_group, cluster_id, n_reference
        consensus_metadata.append({
            'consensus_group': cons.id,
            'cluster_id': cid,
            'n_reference': len(seqs)
        })

    # Stage 3: Consensus length filtering
    # Remove consensus sequences that are suspiciously short compared to median
    if min_consensus_length_ratio > 0 and len(consensus_records) > 1:
        original_count = len(consensus_records)
        consensus_records = filter_consensus_by_length(
            consensus_records,
            min_consensus_length_ratio
        )
        filtered_count = original_count - len(consensus_records)
        if filtered_count > 0:
            logger.info(f"  Filtered {filtered_count} consensus sequences with length <{min_consensus_length_ratio:.0%} of median")
        # Also filter metadata to match
        filtered_ids = set(consensus_records.keys())
        consensus_metadata = [m for m in consensus_metadata if m['consensus_group'] in filtered_ids]

    SeqIO.write(consensus_records.values(), str(consensus_fasta), "fasta")

    # Write consensus metadata CSV
    consensus_metadata_csv = output_dir / f"{organism_name}_consensus_metadata.csv"
    import pandas as pd
    pd.DataFrame(consensus_metadata).to_csv(consensus_metadata_csv, index=False)
    logger.info(f"Wrote consensus metadata to {consensus_metadata_csv}")

    if cleanup_intermediates:
        for p in (aligned_fasta, trimmed_fasta):
            if p.exists():
                p.unlink()
    return consensus_records
    
    
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


def run_trimal_trimming(
    input_fasta: str,
    output_fasta: str,
    trimal_options: Optional[List[str]] = None
) -> None:
    """
    Trim multiple sequence alignment using trimAl.

    Parameters
    ----------
    input_fasta : str
        Path to input aligned FASTA file
    output_fasta : str
        Path to output trimmed FASTA file
    trimal_options : List[str], optional
        Additional trimAl command line options (default: ["-automated1"])

    Raises
    ------
    TrimmingError
        If trimAl is not found or trimming fails

    Notes
    -----
    The -automated1 option automatically selects the best trimming method
    based on alignment properties. It balances gap removal with preservation
    of phylogenetic signal.

    Examples
    --------
    >>> run_trimal_trimming("aligned.fasta", "trimmed.fasta")
    >>> run_trimal_trimming("aligned.fasta", "trimmed.fasta",
    ...                     trimal_options=["-gt", "0.8", "-st", "0.001"])
    """
    # Check if trimAl is available
    if shutil.which('trimal') is None:
        raise TrimmingError(
            "trimAl not found in PATH. Please install trimAl:\n"
            "  - macOS: brew install trimal\n"
            "  - Ubuntu/Debian: sudo apt-get install trimal\n"
            "  - conda: conda install -c bioconda trimal"
        )

    # Set default options
    if trimal_options is None:
        trimal_options = ["-automated1"]

    # Build command
    cmd = ["trimal", "-in", input_fasta, "-out", output_fasta] + trimal_options

    logger.info(f"Running trimAl: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        logger.info(f"trimAl completed: {output_fasta}")

    except subprocess.CalledProcessError as e:
        error_msg = f"trimAl trimming failed:\n{e.stderr}"
        logger.error(error_msg)
        raise TrimmingError(error_msg) from e


def calculate_pairwise_distances(alignment: List[SeqRecord]) -> np.ndarray:
    """
    Calculate pairwise sequence distances ignoring gaps and Ns.

    This function computes a condensed distance matrix suitable for
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

    Notes
    -----
    The distance calculation ignores positions where either sequence has
    a gap or ambiguous base. This focuses on actual sequence differences
    rather than missing data artifacts.

    For large alignments (>1000 sequences), this function may consume
    significant memory. The distance matrix size is O(n²) where n is the
    number of sequences.

    Examples
    --------
    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("trimmed.fasta", "fasta")
    >>> distances = calculate_pairwise_distances(list(alignment))
    >>> print(f"Distance matrix shape: {distances.shape}")
    Distance matrix shape: (4950,)  # For 100 sequences: 100*99/2 = 4950
    """
    if not alignment:
        raise ValueError("Alignment is empty")

    n_seqs = len(alignment)
    if n_seqs < 2:
        raise ValueError("Need at least 2 sequences for distance calculation")

    # Validate alignment length
    aln_length = len(alignment[0].seq)
    if not all(len(rec.seq) == aln_length for rec in alignment):
        raise ValueError("All sequences in alignment must have the same length")

    logger.info(f"Calculating pairwise distances for {n_seqs} sequences")

    # Convert sequences to uppercase strings for faster comparison
    seqs = [str(rec.seq).upper() for rec in alignment]

    # Calculate condensed distance matrix
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

    # Return maximum distance if no valid sites
    if valid_sites == 0:
        return 1.0

    return 1.0 - (matches / valid_sites)


def cluster_sequences(
    distance_matrix: np.ndarray,
    threshold: float = 0.01,
    method: str = "average"
) -> np.ndarray:
    """
    Perform hierarchical clustering on sequence distance matrix.

    Parameters
    ----------
    distance_matrix : np.ndarray
        Condensed distance matrix from calculate_pairwise_distances()
    threshold : float, optional
        Distance threshold for cluster formation (default: 0.01 = 99% identity)
        Sequences with distance ≤ threshold are grouped together
    method : str, optional
        Linkage method for hierarchical clustering (default: "average")
        Options: "single", "complete", "average", "weighted", "ward"

    Returns
    -------
    np.ndarray
        Cluster labels (1-indexed integers) for each sequence

    Raises
    ------
    ClusteringError
        If clustering fails or invalid parameters provided

    Notes
    -----
    - "average" linkage (UPGMA) is recommended for molecular sequences
    - Threshold of 0.01 corresponds to 99% sequence identity
    - Cluster labels are 1-indexed (not 0-indexed)

    Examples
    --------
    >>> distances = calculate_pairwise_distances(alignment)
    >>> labels = cluster_sequences(distances, threshold=0.01)
    >>> print(f"Found {labels.max()} clusters")
    Found 15 clusters
    """
    if len(distance_matrix) == 0:
        raise ClusteringError("Distance matrix is empty")

    if threshold < 0 or threshold > 1:
        raise ClusteringError(
            f"Invalid threshold: {threshold}. Must be between 0 and 1"
        )

    valid_methods = ["single", "complete", "average", "weighted", "ward"]
    if method not in valid_methods:
        raise ClusteringError(
            f"Invalid linkage method: {method}. "
            f"Must be one of {valid_methods}"
        )

    logger.info(
        f"Clustering sequences with {method} linkage, threshold={threshold}"
    )

    try:
        # Perform hierarchical clustering
        linkage_matrix = linkage(distance_matrix, method=method)

        # Cut dendrogram at threshold to get clusters
        labels = fcluster(linkage_matrix, t=threshold, criterion='distance')

        n_clusters = labels.max()
        logger.info(f"Clustering complete: {n_clusters} clusters formed")

        # Log cluster size distribution
        cluster_sizes = Counter(labels)
        logger.debug(f"Cluster size distribution: {dict(cluster_sizes)}")

        return labels

    except Exception as e:
        error_msg = f"Clustering failed: {e}"
        logger.error(error_msg)
        raise ClusteringError(error_msg) from e


def generate_consensus(
    sequences: List[SeqRecord],
    cluster_id: int,
    frequency_cutoff: float = 0.7
) -> SeqRecord:
    """
    Generate consensus sequence using majority rule.

    At each position in the alignment, the most common base is selected
    if it meets the frequency cutoff. Gaps and ambiguous bases (N) are
    ignored when counting. If no base meets the cutoff or no valid bases
    exist, 'N' is assigned.

    Parameters
    ----------
    sequences : List[SeqRecord]
        Aligned sequences in cluster (all must have same length)
    cluster_id : int
        Cluster identifier for naming the consensus
    frequency_cutoff : float, optional
        Minimum fraction (0-1) required to call consensus base (default: 0.7)
        If the most common base occurs in <70% of sequences, 'N' is used

    Returns
    -------
    SeqRecord
        Consensus sequence with ID: consensus_c{cluster_id}_n{sample_count}

    Raises
    ------
    ValueError
        If sequences are empty, have different lengths, or invalid cutoff

    Notes
    -----
    - Only A, C, G, T bases are considered for consensus calling
    - Gaps ('-') and ambiguous bases ('N') are ignored
    - The frequency cutoff helps identify ambiguous positions
    - Lower cutoff (e.g., 0.5) is more permissive, higher (e.g., 0.9) more stringent

    Examples
    --------
    >>> cluster_seqs = [alignment[i] for i in range(10)]
    >>> consensus = generate_consensus(cluster_seqs, cluster_id=1)
    >>> print(consensus.id)
    consensus_c1_n10
    """
    if not sequences:
        raise ValueError("Cannot generate consensus from empty sequence list")

    if frequency_cutoff < 0 or frequency_cutoff > 1:
        raise ValueError(
            f"Invalid frequency_cutoff: {frequency_cutoff}. Must be between 0 and 1"
        )

    # Validate alignment length
    aln_length = len(sequences[0].seq)
    if not all(len(rec.seq) == aln_length for rec in sequences):
        raise ValueError(
            "All sequences must have the same length for consensus generation"
        )

    n_seqs = len(sequences)
    logger.debug(
        f"Generating consensus for cluster {cluster_id} "
        f"({n_seqs} sequences, length={aln_length})"
    )

    consensus_bases = []

    # Process each position in the alignment
    for pos in range(aln_length):
        # Extract column at this position
        column = [rec.seq[pos].upper() for rec in sequences]

        # Filter to valid bases only (A, C, G, T)
        valid_bases = [base for base in column if base in "ACGT"]

        if not valid_bases:
            # No valid bases at this position
            consensus_bases.append('N')
        else:
            # Count base frequencies
            base_counts = Counter(valid_bases)
            most_common_base, count = base_counts.most_common(1)[0]

            # Check if most common base meets frequency cutoff
            frequency = count / len(valid_bases)
            if frequency >= frequency_cutoff:
                consensus_bases.append(most_common_base)
            else:
                # No clear majority
                consensus_bases.append('N')

    # Create consensus sequence
    consensus_seq = "".join(consensus_bases)

    # Degap the consensus sequence (remove alignment gaps and padding)
    # Remove '-' (gaps) and 'N' (ambiguous/padding positions)
    # This converts the aligned consensus (~7500 bp with padding) to raw consensus (~650 bp)
    consensus_seq_degapped = consensus_seq.replace('-', '').replace('N', '')

    consensus_id = f"consensus_c{cluster_id}"

    consensus_record = SeqRecord(
        Seq(consensus_seq_degapped),
        id=consensus_id,
        description=f"Consensus sequence for cluster {cluster_id} ({n_seqs} reference sequences)"
    )

    logger.debug(f"Generated consensus: {consensus_id} (degapped from {len(consensus_seq)} to {len(consensus_seq_degapped)} bp, {n_seqs} reference sequences)")

    return consensus_record


def dereplicate_sequences(
    tsv_path: str,
    output_dir: str,
    threshold: float = 0.01,
    frequency_cutoff: float = 0.7,
    organism_name: Optional[str] = None,
    mafft_options: Optional[List[str]] = None,
    trimal_options: Optional[List[str]] = None,
    cleanup_intermediates: bool = False
) -> Dict[str, SeqRecord]:
    """
    Dereplicate COI sequences through alignment, clustering, and consensus generation.

    This is the main orchestration function that performs the complete dereplication
    workflow:
    1. Extract sequences from BOLD TSV file
    2. Write unaligned FASTA
    3. Run MAFFT alignment
    4. Run trimAl trimming
    5. Calculate pairwise distances
    6. Cluster sequences by distance threshold
    7. Generate consensus for each cluster
    8. Write consensus sequences to FASTA

    Parameters
    ----------
    tsv_path : str
        Path to BOLD TSV file containing sequence data
    output_dir : str
        Directory for output files (created if doesn't exist)
    threshold : float, optional
        Distance threshold for clustering (default: 0.01 = 99% identity)
    frequency_cutoff : float, optional
        Minimum fraction to call consensus base (default: 0.7)
    organism_name : str, optional
        Organism name for file naming (extracted from TSV if None)
    mafft_options : List[str], optional
        Custom MAFFT options (default: ["--auto"])
    trimal_options : List[str], optional
        Custom trimAl options (default: ["-automated1"])
    cleanup_intermediates : bool, optional
        If True, remove intermediate alignment files (default: False)

    Returns
    -------
    Dict[str, SeqRecord]
        Dictionary mapping consensus IDs to SeqRecord objects

    Raises
    ------
    DereplicationError
        If any step in the dereplication workflow fails
    FileNotFoundError
        If TSV file doesn't exist
    ValueError
        If TSV file has no sequences or invalid parameters

    Notes
    -----
    Output files created in output_dir:
    - {organism}_sequences.fasta - Unaligned input sequences
    - {organism}_aligned.fasta - MAFFT alignment
    - {organism}_trimmed.fasta - Trimmed alignment
    - {organism}_consensus.fasta - Final consensus sequences

    If cleanup_intermediates=True, only the consensus FASTA is retained.

    Examples
    --------
    >>> consensus_seqs = dereplicate_sequences(
    ...     tsv_path="Sphyrna_lewini.tsv",
    ...     output_dir="results/",
    ...     threshold=0.01
    ... )
    >>> print(f"Generated {len(consensus_seqs)} consensus sequences")
    Generated 15 consensus sequences

    >>> # With custom options
    >>> consensus_seqs = dereplicate_sequences(
    ...     tsv_path="data.tsv",
    ...     output_dir="results/",
    ...     threshold=0.005,  # 99.5% identity
    ...     frequency_cutoff=0.8,  # Stricter consensus
    ...     mafft_options=["--maxiterate", "1000", "--localpair"],
    ...     cleanup_intermediates=True
    ... )
    """
    logger.info("=" * 70)
    logger.info("Starting sequence dereplication workflow")
    logger.info("=" * 70)

    # Validate inputs
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check external tools
    tools = check_external_tools()
    if not tools['mafft']:
        raise DereplicationError("MAFFT not found. Please install MAFFT.")
    if not tools['trimal']:
        raise DereplicationError("trimAl not found. Please install trimAl.")

    # Extract organism name from TSV filename if not provided
    if organism_name is None:
        organism_name = tsv_path.stem

    # Define output file paths
    sequences_fasta = output_dir / f"{organism_name}_sequences.fasta"
    aligned_fasta = output_dir / f"{organism_name}_aligned.fasta"
    trimmed_fasta = output_dir / f"{organism_name}_trimmed.fasta"
    consensus_fasta = output_dir / f"{organism_name}_consensus.fasta"

    # Step 1: Extract sequences from TSV
    logger.info(f"Step 1/7: Reading sequences from {tsv_path}")
    try:
        df = pd.read_csv(tsv_path, sep='\t')

        if 'nuc' not in df.columns:
            raise ValueError("TSV file missing 'nuc' column")

        # Filter out empty sequences
        df = df[df['nuc'].notna() & (df['nuc'] != '')]

        if len(df) == 0:
            raise ValueError("No valid sequences found in TSV file")

        logger.info(f"Found {len(df)} sequences")

    except Exception as e:
        raise DereplicationError(f"Failed to read TSV file: {e}") from e
        
    # Step 2: Write sequences to FASTA
    logger.info(f"Step 2/7: Writing unaligned sequences to {sequences_fasta}")
    try:
        records = []
        
        # allowed IUPAC DNA + gap
        ALLOWED = set("ACGTURYKMSWBDHVN-")
        
        for idx, row in df.iterrows():
            # Create sequence ID from processid or use index
            seq_id = row.get('processid', f"seq_{idx}")
            
            # Sanitize sequence
            seq_raw = str(row['nuc'])
            seq = re.sub(r'\s+', '', seq_raw).upper().replace('U', 'T')
            
            invalid = set(re.findall(r'[^ACGTURYKMSWBDHVN-]', seq))
            if invalid:
                logger.warning(f"Skipping {seq_id}: invalid characters {invalid}")
                continue
                
            if len(seq) < 100:
                logger.warning(f"Skipping {seq_id}: Sequence too short ({len(seq)} < 100)")
                continue
                
            gap_frac = seq.count('-') / len(seq) if len(seq) else 1.0
            if gap_frac > 0.8:
                logger.warning(f"Skipping {seq_id}: Excessive gaps ({gap_frac:.1%})")
                continue
            
            records.append(SeqRecord(Seq(seq), id=seq_id, description=""))
            
        if not records:
            raise ValueError("No sequences remained after sanitation/filters")

        SeqIO.write(records, sequences_fasta, "fasta")
        logger.info(f"Wrote {len(records)} sequences to FASTA")

    except Exception as e:
        raise DereplicationError(f"Failed to write FASTA file: {e}") from e

    # Step 3: Run MAFFT alignment
    logger.info(f"Step 3/7: Running MAFFT alignment")
    try:
        run_mafft_alignment(
            str(sequences_fasta),
            str(aligned_fasta),
            mafft_options=mafft_options
        )
    except AlignmentError as e:
        raise DereplicationError(f"MAFFT alignment failed: {e}") from e

    # Step 4: Run trimAl trimming
    logger.info(f"Step 4/7: Running trimAl trimming")
    try:
        run_trimal_trimming(
            str(aligned_fasta),
            str(trimmed_fasta),
            trimal_options=trimal_options
        )
    except TrimmingError as e:
        raise DereplicationError(f"trimAl trimming failed: {e}") from e

    # Step 5: Calculate pairwise distances
    logger.info(f"Step 5/7: Calculating pairwise distances")
    try:
        alignment = list(AlignIO.read(trimmed_fasta, "fasta"))
        distances = calculate_pairwise_distances(alignment)
    except Exception as e:
        raise DereplicationError(f"Distance calculation failed: {e}") from e

    # Step 6: Cluster sequences
    logger.info(f"Step 6/7: Clustering sequences (threshold={threshold})")
    try:
        cluster_labels = cluster_sequences(distances, threshold=threshold)
    except ClusteringError as e:
        raise DereplicationError(f"Clustering failed: {e}") from e

    # Step 7: Generate consensus sequences
    logger.info(f"Step 7/7: Generating consensus sequences")
    try:
        # Group sequences by cluster
        clusters = {}
        for seq_record, cluster_id in zip(alignment, cluster_labels):
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(seq_record)

        # Generate consensus for each cluster
        consensus_records = {}
        consensus_metadata = []  # Track reference counts for metadata file
        for cluster_id, cluster_seqs in sorted(clusters.items()):
            consensus = generate_consensus(
                cluster_seqs,
                cluster_id,
                frequency_cutoff=frequency_cutoff
            )
            consensus_records[consensus.id] = consensus
            # Store metadata: consensus_group, cluster_id, n_reference
            consensus_metadata.append({
                'consensus_group': consensus.id,
                'cluster_id': cluster_id,
                'n_reference': len(cluster_seqs)
            })

        # Write consensus sequences to file
        SeqIO.write(consensus_records.values(), consensus_fasta, "fasta")

        # Write consensus metadata CSV
        consensus_metadata_csv = output_dir / f"{organism_name}_consensus_metadata.csv"
        pd.DataFrame(consensus_metadata).to_csv(consensus_metadata_csv, index=False)
        logger.info(f"Wrote consensus metadata to {consensus_metadata_csv}")

        logger.info(f"Generated {len(consensus_records)} consensus sequences")
        logger.info(f"Consensus sequences written to: {consensus_fasta}")

        # Log cluster statistics
        cluster_sizes = [len(seqs) for seqs in clusters.values()]
        logger.info(
            f"Cluster size statistics: "
            f"min={min(cluster_sizes)}, "
            f"max={max(cluster_sizes)}, "
            f"mean={np.mean(cluster_sizes):.1f}, "
            f"median={np.median(cluster_sizes):.1f}"
        )

    except Exception as e:
        raise DereplicationError(f"Consensus generation failed: {e}") from e

    # Cleanup intermediate files if requested
    if cleanup_intermediates:
        logger.info("Cleaning up intermediate files")
        for filepath in [sequences_fasta, aligned_fasta, trimmed_fasta]:
            if filepath.exists():
                filepath.unlink()
                logger.debug(f"Removed: {filepath}")

    logger.info("=" * 70)
    logger.info("Dereplication workflow completed successfully")
    logger.info("=" * 70)

    return consensus_records


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

        # Create haplotype ID
        hap_id = f"haplotype_h{haplotype_id}_n{n_members}"

        # Create SeqRecord with ALIGNED sequence (for distance calculation)
        # The aligned version is needed for pairwise distance matrix
        seq_record = SeqRecord(
            Seq(aligned_seqs[ungapped_seq]),
            id=hap_id,
            description=f"Haplotype {haplotype_id} with {n_members} sequences"
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

    # Merge with ORF validation if provided
    if orf_validation_results is not None:
        haplotype_df = haplotype_df.merge(
            orf_validation_results[['processid', 'orf_valid', 'orf_coverage', 'internal_stops']],
            on='processid',
            how='left'
        )

    logger.info(f"Identified {len(haplotype_records)} unique haplotypes from {len(headers)} sequences")
    n_singletons = sum(1 for h in haplotype_records if h.description.endswith("1 sequences"))
    logger.info(f"  Singletons: {n_singletons} ({n_singletons/len(haplotype_records)*100:.1f}%)")

    return haplotype_records, haplotype_df


def flag_suspect_haplotypes(
    haplotype_records: List[SeqRecord],
    haplotype_df: pd.DataFrame,
    max_singleton_distance: float = 0.05,
    distance_outlier_min_abs: float = 0.10,
    distance_outlier_quantile: float = 0.95
) -> pd.DataFrame:
    """
    Flag suspect haplotypes based on distance and ORF quality.

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


def identify_haplotypes(
    tsv_path: Union[str, Path],
    fasta_path: Union[str, Path],
    output_dir: Union[str, Path],
    min_core_coverage: float = 0.8,
    gap_mask_threshold: float = 0.5,
    min_core_length: int = 200,
    max_singleton_distance: float = 0.05,
    orf_validation_df: Optional[pd.DataFrame] = None,
    cleanup_intermediates: bool = False
) -> Tuple[List[SeqRecord], pd.DataFrame, pd.DataFrame]:
    """
    Main haplotype-first workflow: identify ESVs from core region.

    This is the new recommended workflow replacing the clustering-based approach.

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
        Gap threshold for column masking (default: 0.5)
    min_core_length : int, optional
        Minimum core region length (default: 200 bp)
    max_singleton_distance : float, optional
        Distance threshold for flagging singletons (default: 0.05)
    orf_validation_df : Optional[pd.DataFrame]
        ORF validation results from quality control
    cleanup_intermediates : bool, optional
        Remove intermediate files (default: False)

    Returns
    -------
    Tuple[List[SeqRecord], pd.DataFrame, pd.DataFrame]
        - haplotype_records: Unique haplotype sequences
        - haplotype_mapping: Sample-to-haplotype mapping
        - haplotype_stats: Haplotype-level statistics with flags

    Notes
    -----
    Workflow:
    1. Align sequences with MAFFT
    2. Extract core shared region (covered by ≥min_core_coverage of sequences)
    3. Apply gap masking
    4. Identify unique haplotypes (ESVs)
    5. Calculate distances and flag suspects
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

    # Step 5: Flag suspect haplotypes
    logger.info("Step 5: Flagging suspect haplotypes...")
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
    haplotype_stats.to_csv(haplotype_stats_file, index=False)
    logger.info(f"Wrote haplotype statistics to {haplotype_stats_file}")

    # Cleanup
    if cleanup_intermediates and aligned_fasta.exists():
        aligned_fasta.unlink()
        logger.debug(f"Removed intermediate file: {aligned_fasta}")

    logger.info("=" * 70)
    logger.info("Haplotype-first workflow completed successfully")
    logger.info("=" * 70)

    return haplotype_records, haplotype_mapping, haplotype_stats
