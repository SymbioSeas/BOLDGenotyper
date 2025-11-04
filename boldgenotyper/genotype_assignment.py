"""
Genotype Assignment Module

This module assigns individual sequences to consensus genotype groups based on
sequence similarity. Each sample (identified by processid) is matched to the
closest consensus sequence using global edit distance (Levenshtein distance).

The assignment workflow:
1. Load consensus sequences (genotype representatives)
2. Load raw sequences from FASTA file
3. Extract processids from FASTA headers
4. For each sequence, calculate edit distance to all consensus sequences
5. Assign to best-matching consensus group (if above identity threshold)
6. Generate diagnostic reports for quality control

Key Concepts:
- Edit Distance: Minimum number of insertions, deletions, and substitutions
  needed to transform one sequence into another (global alignment)

- Identity Score: 1 - (edit_distance / max_length), represents similarity

- Minimum Identity Threshold: Default 0.90 (90% similarity required for assignment)
  Sequences below this threshold are marked as unassigned

- Tie Detection: Flags cases where multiple consensus groups have nearly
  identical scores, indicating ambiguous assignment

- Runner-up Tracking: Records second-best match for diagnostic purposes

Dependencies:
- edlib (optional): Fast C-based edit distance calculation
- Falls back to pure Python implementation if edlib unavailable

Example Usage:
    >>> from boldgenotyper.genotype_assignment import assign_genotypes
    >>> results = assign_genotypes(
    ...     metadata_path="data.tsv",
    ...     fasta_path="sequences.fasta",
    ...     consensus_path="consensus.fasta",
    ...     output_path="data_with_genotypes.tsv",
    ...     min_identity=0.90,
    ...     n_processes=4
    ... )

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
import re
import csv
from collections import namedtuple
from functools import partial
import multiprocessing as mp

import pandas as pd
import numpy as np

# import taxonomy dataclass
try:
    from .config import TaxonomyConfig
except Exception:
    from config import TaxonomyConfig

# Configure logging
logger = logging.getLogger(__name__)

# Try to import edlib for fast edit distance
try:
    import edlib
    EDLIB_AVAILABLE = True
except ImportError:
    EDLIB_AVAILABLE = False
    edlib = None

# Named tuple for worker tasks
AssignmentTask = namedtuple("AssignmentTask", ["processid", "sequence"])


class GenotypeAssignmentError(Exception):
    """Base exception for genotype assignment errors."""
    pass


class ProcessIDExtractionError(GenotypeAssignmentError):
    """Error extracting processid from FASTA header."""
    pass


# Regular expression to extract processid from FASTA headers
# Matches the last underscore-delimited field before a dot or end-of-line
# Example: "Sphyrna_lewini_ANGBF11456-15.COI-5P" -> "ANGBF11456-15"
PROCESSID_REGEX = re.compile(r"_(?P<pid>[^.\s_]+)(?:[.\s]|$)")


def check_edlib_available() -> bool:
    """
    Check if edlib library is available for fast edit distance.

    Returns
    -------
    bool
        True if edlib is available, False otherwise

    Examples
    --------
    >>> if check_edlib_available():
    ...     print("Using edlib for fast alignment")
    ... else:
    ...     print("Using pure Python implementation")
    """
    return EDLIB_AVAILABLE


def levenshtein_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Levenshtein distance (edit distance) between two sequences.

    This pure Python implementation uses dynamic programming with O(n*m) time
    complexity and O(min(n,m)) space complexity (optimized to use only two rows).

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    int
        Edit distance (minimum number of insertions, deletions, substitutions)

    Notes
    -----
    The algorithm swaps seq1 and seq2 if seq1 is longer, to minimize memory usage.
    This implementation is used as a fallback when edlib is not available.

    Examples
    --------
    >>> levenshtein_distance("ACGT", "ACGT")
    0
    >>> levenshtein_distance("ACGT", "AGGT")
    1
    >>> levenshtein_distance("ACGT", "TGCA")
    4
    """
    # Ensure seq1 is the shorter sequence for memory efficiency
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1

    # Initialize first row (0, 1, 2, ..., len(seq1))
    previous_row = list(range(len(seq1) + 1))

    # Dynamic programming: build up distance matrix row by row
    for j, char2 in enumerate(seq2, start=1):
        current_row = [j]  # First column value
        for i, char1 in enumerate(seq1, start=1):
            # Calculate costs for three operations
            insertion_cost = previous_row[i] + 1
            deletion_cost = current_row[i - 1] + 1
            substitution_cost = previous_row[i - 1] + (char1 != char2)

            # Take minimum cost
            current_row.append(min(insertion_cost, deletion_cost, substitution_cost))

        previous_row = current_row

    return previous_row[-1]


def calculate_edit_distance(seq1: str, seq2: str, use_edlib: bool = True) -> int:
    """
    Calculate edit distance between two sequences using edlib or fallback.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence
    use_edlib : bool, optional
        If True and edlib available, use edlib (default: True)

    Returns
    -------
    int
        Edit distance

    Examples
    --------
    >>> calculate_edit_distance("ACGT", "ACGT")
    0
    >>> calculate_edit_distance("ACGT", "AGGT", use_edlib=False)
    1
    """
    if use_edlib and EDLIB_AVAILABLE:
        result = edlib.align(seq1, seq2, mode="NW", task="distance")
        return result["editDistance"]
    else:
        return levenshtein_distance(seq1, seq2)


def calculate_identity(edit_distance: int, len1: int, len2: int) -> float:
    """
    Calculate sequence identity from edit distance.

    Identity is calculated as: 1 - (edit_distance / max_length)
    This represents the fraction of positions that match.

    Parameters
    ----------
    edit_distance : int
        Edit distance between sequences
    len1 : int
        Length of first sequence
    len2 : int
        Length of second sequence

    Returns
    -------
    float
        Identity score between 0.0 and 1.0

    Examples
    --------
    >>> calculate_identity(0, 100, 100)
    1.0
    >>> calculate_identity(10, 100, 100)
    0.9
    >>> calculate_identity(50, 100, 100)
    0.5
    """
    max_length = max(len1, len2)
    if max_length == 0:
        return 1.0
    return 1.0 - (edit_distance / max_length)


def extract_processid_from_header(header: str) -> Optional[str]:
    """
    Extract processid from FASTA header.

    The processid is extracted as the last underscore-delimited field
    before a dot or end-of-line.

    Parameters
    ----------
    header : str
        FASTA header line (with or without leading >)

    Returns
    -------
    str or None
        Extracted processid, or None if pattern not found

    Examples
    --------
    >>> extract_processid_from_header("Sphyrna_lewini_ANGBF11456-15.COI-5P")
    'ANGBF11456-15'
    >>> extract_processid_from_header("species_BOLD12345")
    'BOLD12345'
    >>> extract_processid_from_header(">Genus_species_ABC123.marker")
    'ABC123'
    """
    match = PROCESSID_REGEX.search(header)
    if match:
        return match.group("pid")
    return None


def read_fasta_simple(fasta_path: str) -> List[Tuple[str, str]]:
    """
    Read FASTA file without BioPython dependency.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file

    Returns
    -------
    List[Tuple[str, str]]
        List of (header, sequence) tuples

    Raises
    ------
    FileNotFoundError
        If FASTA file doesn't exist
    GenotypeAssignmentError
        If FASTA file is malformed

    Examples
    --------
    >>> records = read_fasta_simple("sequences.fasta")
    >>> for header, seq in records:
    ...     print(f"{header}: {len(seq)} bp")
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    records = []
    current_header = None
    current_seq_lines = []

    try:
        with open(fasta_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip('\n')

                if not line:  # Skip empty lines
                    continue

                if line.startswith('>'):
                    # Save previous record
                    if current_header is not None:
                        sequence = "".join(current_seq_lines).upper()
                        records.append((current_header, sequence))

                    # Start new record
                    current_header = line[1:].strip()
                    current_seq_lines = []
                else:
                    current_seq_lines.append(line.strip())

            # Save last record
            if current_header is not None:
                sequence = "".join(current_seq_lines).upper()
                records.append((current_header, sequence))

    except Exception as e:
        raise GenotypeAssignmentError(
            f"Error reading FASTA file {fasta_path}: {e}"
        ) from e

    if not records:
        raise GenotypeAssignmentError(f"No sequences found in FASTA file: {fasta_path}")

    logger.info(f"Read {len(records)} sequences from {fasta_path}")
    return records


def find_best_consensus_match(
    sequence: str,
    consensus_groups: List[Tuple[str, str]],
    min_identity: float = 0.90,
    use_edlib: bool = True
) -> Dict[str, Any]:
    """
    Find best matching consensus group for a sequence.

    Parameters
    ----------
    sequence : str
        Query sequence to assign
    consensus_groups : List[Tuple[str, str]]
        List of (group_id, consensus_sequence) tuples
    min_identity : float, optional
        Minimum identity required for assignment (default: 0.90)
    use_edlib : bool, optional
        Use edlib if available (default: True)

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:
        - 'best_group': Best matching group ID (or None if below threshold)
        - 'best_identity': Identity to best match
        - 'runner_up_group': Second best group ID
        - 'runner_up_identity': Identity to runner-up
        - 'is_tie': Boolean, True if best and runner-up are very close (diff < 0.01)
        - 'is_low_confidence': Boolean, True if best identity < min_identity + 0.05

    Examples
    --------
    >>> consensus = [("c1_n10", "ACGTACGT"), ("c2_n5", "TGCATGCA")]
    >>> result = find_best_consensus_match("ACGTACGT", consensus)
    >>> result['best_group']
    'c1_n10'
    >>> result['best_identity']
    1.0
    """
    best_group = None
    best_identity = -1.0
    runner_up_group = None
    runner_up_identity = -1.0

    # Compare to all consensus sequences
    for group_id, consensus_seq in consensus_groups:
        # Calculate edit distance
        edit_dist = calculate_edit_distance(sequence, consensus_seq, use_edlib=use_edlib)

        # Convert to identity score
        identity = calculate_identity(edit_dist, len(sequence), len(consensus_seq))

        # Update best and runner-up
        if identity > best_identity:
            # Current best becomes runner-up
            runner_up_group = best_group
            runner_up_identity = best_identity
            # New best
            best_group = group_id
            best_identity = identity
        elif identity > runner_up_identity:
            # Update runner-up
            runner_up_group = group_id
            runner_up_identity = identity

    # Determine if assignment is below threshold
    if best_identity < min_identity:
        best_group = None  # Unassigned

    # Detect ties (best and runner-up are very close)
    is_tie = False
    if runner_up_identity > 0 and (best_identity - runner_up_identity) < 0.01:
        is_tie = True

    # Detect low-confidence assignments (barely above threshold)
    is_low_confidence = False
    if best_group is not None and best_identity < (min_identity + 0.05):
        is_low_confidence = True

    return {
        'best_group': best_group,
        'best_identity': best_identity,
        'runner_up_group': runner_up_group,
        'runner_up_identity': runner_up_identity,
        'is_tie': is_tie,
        'is_low_confidence': is_low_confidence
    }


def _assignment_worker(
    task: Tuple[str, Optional[str]],
    consensus_groups: List[Tuple[str, str]],
    min_identity: float,
    use_edlib: bool
) -> Dict[str, Any]:
    """
    Worker function for parallel genotype assignment.

    Parameters
    ----------
    task : Tuple[str, Optional[str]]
        (processid, sequence) tuple
    consensus_groups : List[Tuple[str, str]]
        List of (group_id, consensus_sequence)
    min_identity : float
        Minimum identity threshold
    use_edlib : bool
        Whether to use edlib

    Returns
    -------
    Dict[str, Any]
        Assignment results for this processid
    """
    processid, sequence = task

    # Handle missing sequence
    if sequence is None or sequence == "":
        return {
            'processid': processid,
            'consensus_group': None,
            'identity': 0.0,
            'runner_up_group': None,
            'runner_up_identity': 0.0,
            'is_tie': False,
            'is_low_confidence': False,
            'status': 'no_sequence'
        }

    # Find best match
    result = find_best_consensus_match(
        sequence, consensus_groups, min_identity, use_edlib
    )

    # Add processid and status
    result['processid'] = processid
    if result['best_group'] is None:
        result['status'] = 'below_threshold'
    elif result['is_tie']:
        result['status'] = 'tie'
    elif result['is_low_confidence']:
        result['status'] = 'low_confidence'
    else:
        result['status'] = 'assigned'

    # Rename for output consistency
    result['consensus_group'] = result.pop('best_group')
    result['identity'] = result.pop('best_identity')

    return result


def assign_species_to_sample(
    query_fasta: str,
    db_path: str,
    config: TaxonomyConfig
) -> List[Dict]:
    """
    Run BLASTn/VSEARCH for each query sequence and return per-sample 
    seq-based assignment with level, identity, qcov, ties, etc.
    """
    # 1) run search (subprocess to blastn -task megablast -perc_identity; or vsearch --usearch_global)
    # 2) parse hits, compute coverage and identity; sort, compute top2 delta, ties
    # 3) apply thresholds & LCA to species/genus; set fields:
    #    seq_sp, seq_level, seq_best_identity, seq_qcov, seq_top2_delta, n_top_ties, low_confidence_flag
    # 4) return records keyed by processid
    
def assign_species_to_consensus(
    consensus_fasta: str,
    db_path: str,
    config: TaxonomyConfig
) -> pd.DataFrame:
    """
    Classify each consensus sequence (untrimmed from Step 6 FASTA)
    with same thresholds; return DataFrame:
      consensus_group, cluster_seq_sp, cluster_seq_level, cluster_seq_best_identity, cluster_seq_qcov, cluster_seq_top2_delta, n_top_ties
    """

def assign_genotypes(
    metadata_path: str,
    fasta_path: str,
    consensus_path: str,
    output_path: str,
    min_identity: float = 0.90,
    n_processes: int = 1,
    diagnostics_path: Optional[str] = None
) -> Dict[str, Any]:
    """
    Assign genotype groups to sequences based on consensus matching.

    This is the main orchestration function that:
    1. Loads metadata, raw sequences, and consensus sequences
    2. Extracts processids from FASTA headers
    3. Assigns each sequence to best-matching consensus group
    4. Generates diagnostic reports
    5. Writes updated metadata with genotype assignments

    Parameters
    ----------
    metadata_path : str
        Path to metadata TSV file (must have 'processid' column)
    fasta_path : str
        Path to raw sequences FASTA file
    consensus_path : str
        Path to consensus sequences FASTA file
    output_path : str
        Path for output TSV with 'consensus_group' column added
    min_identity : float, optional
        Minimum identity threshold for assignment (default: 0.90)
    n_processes : int, optional
        Number of parallel processes (default: 1)
    diagnostics_path : str, optional
        Path for diagnostics CSV output (default: None, no diagnostics)

    Returns
    -------
    Dict[str, Any]
        Summary statistics:
        - 'total': Total number of samples
        - 'assigned': Number successfully assigned
        - 'unassigned': Number not assigned
        - 'no_sequence': Number with missing sequences
        - 'below_threshold': Number below identity threshold
        - 'ties': Number with ambiguous ties
        - 'low_confidence': Number with low-confidence assignments

    Raises
    ------
    FileNotFoundError
        If input files don't exist
    GenotypeAssignmentError
        If metadata is missing 'processid' column or other errors
    ValueError
        If parameters are invalid

    Examples
    --------
    >>> results = assign_genotypes(
    ...     metadata_path="data.tsv",
    ...     fasta_path="sequences.fasta",
    ...     consensus_path="consensus.fasta",
    ...     output_path="data_with_genotypes.tsv",
    ...     min_identity=0.90,
    ...     n_processes=4,
    ...     diagnostics_path="diagnostics.csv"
    ... )
    >>> print(f"Assigned {results['assigned']}/{results['total']} sequences")
    """
    logger.info("=" * 70)
    logger.info("Starting genotype assignment workflow")
    logger.info("=" * 70)

    # Validate inputs
    metadata_path = Path(metadata_path)
    fasta_path = Path(fasta_path)
    consensus_path = Path(consensus_path)
    output_path = Path(output_path)

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not consensus_path.exists():
        raise FileNotFoundError(f"Consensus file not found: {consensus_path}")

    if min_identity < 0 or min_identity > 1:
        raise ValueError(f"min_identity must be between 0 and 1, got {min_identity}")

    if n_processes < 1:
        raise ValueError(f"n_processes must be >= 1, got {n_processes}")

    # Check edlib availability
    use_edlib = EDLIB_AVAILABLE
    if use_edlib:
        logger.info("Using edlib for fast edit distance calculation")
    else:
        logger.info("edlib not available, using pure Python Levenshtein distance")
        logger.warning("Consider installing edlib for faster performance: pip install edlib")

    # Step 1: Load metadata
    logger.info(f"Step 1/6: Loading metadata from {metadata_path}")
    try:
        metadata_df = pd.read_csv(metadata_path, sep='\t', dtype=str)
        if 'processid' not in metadata_df.columns:
            raise GenotypeAssignmentError(
                "Metadata TSV missing required 'processid' column"
            )
        logger.info(f"Loaded {len(metadata_df)} samples from metadata")
    except Exception as e:
        raise GenotypeAssignmentError(f"Failed to load metadata: {e}") from e

    # Step 2: Load consensus sequences
    logger.info(f"Step 2/6: Loading consensus sequences from {consensus_path}")
    consensus_records = read_fasta_simple(str(consensus_path))
    consensus_groups = []
    for header, seq in consensus_records:
        # Extract consensus group ID (e.g., "consensus_c7_n381")
        group_id = header.split()[0]  # Take first word
        consensus_groups.append((group_id, seq))
    logger.info(f"Loaded {len(consensus_groups)} consensus groups")

    # Step 3: Load raw sequences and extract processids
    logger.info(f"Step 3/6: Loading raw sequences from {fasta_path}")
    raw_records = read_fasta_simple(str(fasta_path))

    processid_to_seq = {}
    headers_without_processid = []

    for header, seq in raw_records:
        processid = extract_processid_from_header(header)
        if processid:
            # Store first occurrence if duplicates exist
            processid_to_seq.setdefault(processid, seq)
        else:
            headers_without_processid.append(header)

    logger.info(f"Extracted {len(processid_to_seq)} unique processids from FASTA")

    if headers_without_processid:
        logger.warning(
            f"Could not extract processid from {len(headers_without_processid)} headers"
        )
        if len(headers_without_processid) <= 5:
            for header in headers_without_processid:
                logger.warning(f"  Missing processid: {header}")

    # Step 4: Prepare tasks for parallel processing
    logger.info(f"Step 4/6: Preparing assignment tasks")
    processids = metadata_df['processid'].astype(str).tolist()

    tasks = []
    for pid in processids:
        sequence = processid_to_seq.get(pid)
        tasks.append((pid, sequence))

    logger.info(f"Prepared {len(tasks)} assignment tasks")

    # Step 5: Perform parallel assignment
    logger.info(f"Step 5/6: Assigning genotypes (using {n_processes} processes)")

    worker_func = partial(
        _assignment_worker,
        consensus_groups=consensus_groups,
        min_identity=min_identity,
        use_edlib=use_edlib
    )

    if n_processes > 1:
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(worker_func, tasks)
    else:
        results = list(map(worker_func, tasks))

    logger.info(f"Completed genotype assignment for {len(results)} samples")

    # Step 6: Collect results and generate outputs
    logger.info(f"Step 6/6: Writing outputs")

    # Create mapping for metadata
    processid_to_group = {}
    for result in results:
        processid_to_group[result['processid']] = result['consensus_group']

    # Add consensus_group column to metadata
    metadata_df['consensus_group'] = metadata_df['processid'].map(processid_to_group)

    # Write updated metadata
    metadata_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote updated metadata to {output_path}")

    # Write diagnostics if requested
    if diagnostics_path:
        diagnostics_path = Path(diagnostics_path)
        with open(diagnostics_path, 'w', newline='') as f:
            fieldnames = [
                'processid', 'consensus_group', 'identity',
                'runner_up_group', 'runner_up_identity',
                'is_tie', 'is_low_confidence', 'status'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for result in results:
                writer.writerow({
                    'processid': result['processid'],
                    'consensus_group': result['consensus_group'] or '',
                    'identity': round(result['identity'], 6),
                    'runner_up_group': result['runner_up_group'] or '',
                    'runner_up_identity': round(result['runner_up_identity'], 6),
                    'is_tie': result['is_tie'],
                    'is_low_confidence': result['is_low_confidence'],
                    'status': result['status']
                })

        logger.info(f"Wrote diagnostics to {diagnostics_path}")

    # Calculate summary statistics
    stats = {
        'total': len(results),
        'assigned': sum(1 for r in results if r['consensus_group'] is not None),
        'unassigned': sum(1 for r in results if r['consensus_group'] is None),
        'no_sequence': sum(1 for r in results if r['status'] == 'no_sequence'),
        'below_threshold': sum(1 for r in results if r['status'] == 'below_threshold'),
        'ties': sum(1 for r in results if r['is_tie']),
        'low_confidence': sum(1 for r in results if r['is_low_confidence'])
    }

    # Log summary
    logger.info("=" * 70)
    logger.info("Genotype assignment summary:")
    logger.info(f"  Total samples: {stats['total']}")
    logger.info(f"  Successfully assigned: {stats['assigned']} ({100*stats['assigned']/stats['total']:.1f}%)")
    logger.info(f"  Unassigned: {stats['unassigned']} ({100*stats['unassigned']/stats['total']:.1f}%)")
    logger.info(f"    - No sequence in FASTA: {stats['no_sequence']}")
    logger.info(f"    - Below identity threshold: {stats['below_threshold']}")
    logger.info(f"  Diagnostic flags:")
    logger.info(f"    - Ties (ambiguous): {stats['ties']}")
    logger.info(f"    - Low confidence: {stats['low_confidence']}")
    logger.info("=" * 70)

    return stats
