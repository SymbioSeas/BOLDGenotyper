"""
Sequence Clustering and Consensus Generation

This module handles the dereplication of COI sequences from BOLD data through
multiple sequence alignment, distance calculation, hierarchical clustering, and
consensus sequence generation.

The dereplication workflow:
1. Extract COI sequences from BOLD TSV file
2. Create FASTA file with formatted headers
3. Perform multiple sequence alignment using MAFFT
4. Trim alignment using trimAl (--automated1)
5. Calculate pairwise sequence distances (ignoring gaps and Ns)
6. Perform hierarchical clustering based on distance threshold
7. Generate consensus sequences for each cluster using majority rule

Key Concepts:
- Default clustering threshold: 0.02 (98% sequence identity)
  Rationale: This is the standard threshold for COI-based species delimitation
  and captures meaningful intraspecific genetic variation while distinguishing
  genotypes that may represent distinct evolutionary lineages.

- Distance calculation: Pairwise comparisons ignore gaps and ambiguous (N) bases
  to focus on actual sequence differences rather than missing data artifacts.

- Consensus generation: Uses majority rule voting at each position. Ambiguous
  bases are assigned when no clear majority exists (e.g., 50/50 split).

- Output naming: consensus_c{cluster_id}_n{sample_count}
  Example: consensus_c1_n45 represents cluster 1 with 45 samples

Dependencies:
- MAFFT v7+ for multiple sequence alignment
- trimAl for alignment trimming
- scipy for hierarchical clustering
- Biopython for sequence handling

Example Usage:
    >>> from boldgenotyper.dereplication import dereplicate_sequences
    >>> consensus_seqs = dereplicate_sequences(
    ...     tsv_path="Sphyrna_lewini.tsv",
    ...     output_dir="results/",
    ...     threshold=0.02
    ... )

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
                logger.warning(f"Skipping {seq_id}: Excessive gaps ({gap_frac:..1%})")
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
