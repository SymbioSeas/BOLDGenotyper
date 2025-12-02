"""
Haplotype Assignment Module

This module assigns individual sequences to haplotype groups based on
sequence similarity. Each sample (identified by processid) is matched to the
closest haplotype sequence using global edit distance (Levenshtein distance).

The assignment workflow:
1. Load haplotype sequences (unique ESVs from core region)
2. Load raw sequences from FASTA file
3. Extract processids from FASTA headers
4. For each sequence, calculate edit distance to all haplotype sequences
5. Assign to best-matching haplotype (if above identity threshold)
6. Generate diagnostic reports for quality control

Key Concepts:
- Edit Distance: Minimum number of insertions, deletions, and substitutions
  needed to transform one sequence into another (global alignment)

- Identity Score: 1 - (edit_distance / max_length), represents similarity

- Minimum Identity Threshold: Default 0.90 (90% similarity required for assignment)
  Sequences below this threshold are marked as unassigned

- Tie Detection: Flags cases where multiple haplotypes have nearly
  identical scores, indicating ambiguous assignment

- Runner-up Tracking: Records second-best match for diagnostic purposes

Terminology:
- **Haplotype**: Unique sequence variant (ESV) identified from core region
- **Haplotype ID**: Format haplotype_hX_nY where X = haplotype number, Y = sample count
- Previous "genotype" and "haplotype_id" terminology updated to "haplotype"

Dependencies:
- edlib (optional): Fast C-based edit distance calculation
- Falls back to pure Python implementation if edlib unavailable

Example Usage:
    >>> from boldgenotyper.haplotype_assignment import assign_haplotypes
    >>> results = assign_haplotypes(
    ...     metadata_path="data.tsv",
    ...     fasta_path="sequences.fasta",
    ...     haplotype_path="haplotypes.fasta",
    ...     output_path="data_with_haplotypes.tsv",
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


class HaplotypeAssignmentError(Exception):
    """Base exception for haplotype assignment errors."""
    pass


class ProcessIDExtractionError(HaplotypeAssignmentError):
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


def parse_cigar(cigar: str) -> Dict[str, int]:
    """
    Parse CIGAR string and return operation counts.

    CIGAR format uses extended CIGAR notation:
    - '=' : Match (bases are identical)
    - 'X' : Mismatch (bases differ)
    - 'I' : Insertion to target (extra bases in query)
    - 'D' : Deletion from target (missing bases in query)

    Parameters
    ----------
    cigar : str
        CIGAR string from edlib alignment (e.g., "5=1X4=2I")

    Returns
    -------
    Dict[str, int]
        Dictionary with keys '=', 'X', 'I', 'D' and their counts

    Examples
    --------
    >>> parse_cigar("10=")
    {'=': 10, 'X': 0, 'I': 0, 'D': 0}
    >>> parse_cigar("5=1X4=")
    {'=': 9, 'X': 1, 'I': 0, 'D': 0}
    >>> parse_cigar("10=5I")
    {'=': 10, 'X': 0, 'I': 5, 'D': 0}
    """
    operations = {'=': 0, 'X': 0, 'I': 0, 'D': 0}

    # Parse CIGAR: numbers followed by operation symbols
    # Pattern matches: one or more digits followed by operation character
    pattern = r'(\d+)([=XIDM])'
    for count_str, op in re.findall(pattern, cigar):
        count = int(count_str)
        if op in operations:
            operations[op] += count
        elif op == 'M':
            # 'M' can mean match or mismatch in standard CIGAR
            # In edlib extended CIGAR, this shouldn't occur, but handle it
            logger.warning(f"Found 'M' operation in CIGAR (ambiguous). CIGAR: {cigar}")
            operations['='] += count

    return operations


def calculate_identity_with_cigar(
    seq1: str,
    seq2: str,
    use_edlib: bool = True
) -> Dict[str, Any]:
    """
    Calculate identity using alignment path (CIGAR string).

    This function uses edlib's full alignment path to compute a more accurate
    identity metric that is robust to length differences. It calculates:

    - **target_identity**: matches / len(target)
      This asks: "What fraction of the consensus is represented in the sample?"
      Recommended for haplotype assignment where consensus is the reference.

    - **classic_identity**: 1 - (edit_distance / max_length)
      The original method (kept for comparison and backwards compatibility).

    The target_identity metric is more lenient when samples have extra bases
    (e.g., noisy 5'/3' ends) but still penalizes missing bases appropriately.

    Parameters
    ----------
    seq1 : str
        Query sequence (typically the sample)
    seq2 : str
        Target sequence (typically the consensus)
    use_edlib : bool, optional
        If True and edlib available, use edlib with task='path' (default: True)
        If False or edlib unavailable, falls back to classic calculation

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - edit_distance : int - Total edit operations
        - matches : int - Number of matching bases (from CIGAR '=')
        - mismatches : int - Number of mismatching bases (from CIGAR 'X')
        - insertions : int - Insertions to target (from CIGAR 'I')
        - deletions : int - Deletions from target (from CIGAR 'D')
        - target_length : int - Length of target sequence
        - query_length : int - Length of query sequence
        - target_identity : float - matches / target_length
        - classic_identity : float - 1 - (edit_dist / max_length)
        - alignment_identity : float - matches / (matches + mismatches + gaps)
        - cigar : str - Full CIGAR string (or None if not using edlib)
        - method : str - "edlib_cigar" or "classic_fallback"

    Raises
    ------
    ValueError
        If both sequences are empty

    Notes
    -----
    When edlib is not available or use_edlib=False, the function falls back
    to the classic identity calculation. In this case, CIGAR-derived fields
    (matches, mismatches, insertions, deletions, cigar) will be None or 0.

    Examples
    --------
    >>> result = calculate_identity_with_cigar("ACTGACTG", "ACTGACTG")
    >>> result['target_identity']
    1.0
    >>> result['matches']
    8

    >>> # Sample with noisy 3' end
    >>> result = calculate_identity_with_cigar("ACTGNNNN", "ACTG")
    >>> result['target_identity']  # 100% of consensus is represented
    1.0
    >>> result['classic_identity']  # Penalized for extra bases
    0.5
    """
    # Handle empty sequences
    if len(seq1) == 0 and len(seq2) == 0:
        raise ValueError("Cannot calculate identity for two empty sequences")

    if len(seq1) == 0 or len(seq2) == 0:
        alignment_length = len(seq1) + len(seq2)
        return {
            'edit_distance': max(len(seq1), len(seq2)),
            'matches': 0,
            'mismatches': 0,
            'insertions': 0,
            'deletions': 0,
            'target_length': len(seq2),
            'query_length': len(seq1),
            'target_identity': 0.0,
            'classic_identity': 0.0,
            'alignment_length': alignment_length,
            'alignment_identity': 0.0,
            'cigar': None,
            'method': 'empty_sequence'
        }

    # Fast path for identical sequences
    if seq1 == seq2:
        alignment_length = len(seq1)
        return {
            'edit_distance': 0,
            'matches': len(seq1),
            'mismatches': 0,
            'insertions': 0,
            'deletions': 0,
            'target_length': len(seq2),
            'query_length': len(seq1),
            'target_identity': 1.0,
            'classic_identity': 1.0,
            'alignment_length': alignment_length,
            'alignment_identity': 1.0,
            'cigar': f"{len(seq1)}=",
            'method': 'identical'
        }

    # Use edlib with CIGAR if available
    if use_edlib and EDLIB_AVAILABLE:
        try:
            # Request full alignment path
            result = edlib.align(seq1, seq2, mode="NW", task="path")

            if result['editDistance'] == -1:
                # Alignment failed (shouldn't happen with mode="NW")
                logger.warning(f"Edlib alignment failed for sequences of length {len(seq1)}, {len(seq2)}")
                raise RuntimeError("Edlib alignment failed")

            # Parse CIGAR to get operation counts
            cigar = result['cigar']
            ops = parse_cigar(cigar)

            # Calculate identities
            target_length = len(seq2)
            max_length = max(len(seq1), len(seq2))
            alignment_length = sum(ops.values())

            target_identity = ops['='] / target_length if target_length > 0 else 0.0
            classic_identity = 1.0 - (result['editDistance'] / max_length) if max_length > 0 else 0.0
            alignment_identity = ops['='] / alignment_length if alignment_length > 0 else 0.0

            return {
                'edit_distance': result['editDistance'],
                'matches': ops['='],
                'mismatches': ops['X'],
                'insertions': ops['I'],
                'deletions': ops['D'],
                'target_length': target_length,
                'query_length': len(seq1),
                'target_identity': target_identity,
                'classic_identity': classic_identity,
                'alignment_length': alignment_length,
                'alignment_identity': alignment_identity,
                'cigar': cigar,
                'method': 'edlib_cigar'
            }

        except Exception as e:
            logger.warning(f"Error using edlib CIGAR calculation: {e}. Falling back to classic method.")
            # Fall through to classic method

    # Fallback: classic identity calculation without CIGAR
    edit_dist = calculate_edit_distance(seq1, seq2, use_edlib=False)
    classic_identity = calculate_identity(edit_dist, len(seq1), len(seq2))

    # Without CIGAR we approximate match counts using the max length
    max_length = max(len(seq1), len(seq2))
    target_length = len(seq2)
    approx_matches = max_length - edit_dist
    alignment_length = max_length if max_length > 0 else 0
    alignment_identity = classic_identity if alignment_length > 0 else 0.0

    return {
        'edit_distance': edit_dist,
        'matches': approx_matches,
        'mismatches': None,  # Unknown without CIGAR
        'insertions': None,  # Unknown without CIGAR
        'deletions': None,  # Unknown without CIGAR
        'target_length': target_length,
        'query_length': len(seq1),
        'target_identity': approx_matches / target_length if target_length > 0 else 0.0,
        'classic_identity': classic_identity,
        'alignment_length': alignment_length,
        'alignment_identity': alignment_identity,
        'cigar': None,
        'method': 'classic_fallback'
    }


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
    HaplotypeAssignmentError
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
        raise HaplotypeAssignmentError(
            f"Error reading FASTA file {fasta_path}: {e}"
        ) from e

    if not records:
        raise HaplotypeAssignmentError(f"No sequences found in FASTA file: {fasta_path}")

    logger.info(f"Read {len(records)} sequences from {fasta_path}")
    return records


def find_best_haplotype_match(
    sequence: str,
    haplotype_groups: List[Tuple[str, str]],
    min_identity: float = 0.90,
    use_edlib: bool = True,
    identity_method: str = "target_based",
    tie_margin: float = 0.003,
    tie_min_identity: float = 0.95,
) -> Dict[str, Any]:
    """
    Find best matching haplotype for a sequence.

    Supports two identity calculation methods:
    - "target_based" (default): matches / haplotype_length
      More robust to length differences and noisy 5'/3' ends
    - "classic": 1 - (edit_distance / max_length)
      Original method for backwards compatibility

    Parameters
    ----------
    sequence : str
        Query sequence to assign
    haplotype_groups : List[Tuple[str, str]]
        List of (haplotype_id, haplotype_sequence) tuples
    min_identity : float, optional
        Minimum identity required for assignment (default: 0.90)
    use_edlib : bool, optional
        Use edlib if available (default: True)
    identity_method : str, optional
        Identity calculation method: "target_based" or "classic" (default: "target_based")
    tie_margin : float, optional
        Maximum allowed difference between best and runner-up identity to be considered a tie (default: 0.003 = 0.3%).
    tie_min_identity : float, optional
        Minimum best-identity required before we even consider calling a tie (default: 0.95).

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:
        - 'best_haplotype': Best matching haplotype ID (or None if below threshold)
        - 'best_identity': Identity to best match (using selected method)
        - 'classic_identity': Identity using classic metric (for comparison)
        - 'target_identity': Identity using target-based metric (for comparison)
        - 'matches': Number of matching bases
        - 'mismatches': Number of mismatching bases
        - 'insertions': Number of insertions
        - 'deletions': Number of deletions
        - 'edit_distance': Total edit distance
        - 'cigar': CIGAR string for best match
        - 'length_discrepancy': abs(query_length - target_length)
        - 'identity_method': Which method was used for best_identity
        - 'runner_up_haplotype': Second best haplotype ID
        - 'runner_up_identity': Identity to runner-up
        - 'is_tie': Boolean, True if best and runner-up are very close (diff < 0.01)
        - 'is_low_confidence': Boolean, True if best identity < min_identity + 0.05

    Examples
    --------
    >>> haplotypes = [("h1_n10", "ACGTACGT"), ("h2_n5", "TGCATGCA")]
    >>> result = find_best_haplotype_match("ACGTACGT", haplotypes)
    >>> result['best_haplotype']
    'h1_n10'
    >>> result['best_identity']
    1.0
    """
    best_haplotype = None
    best_identity = -1.0
    best_details = None  # Store detailed alignment info for best match
    runner_up_haplotype = None
    runner_up_identity = -1.0
    runner_up_details = None

    # Validate identity_method
    if identity_method not in ["target_based", "classic"]:
        raise ValueError(
            f"identity_method must be 'target_based' or 'classic', "
            f"got '{identity_method}'"
        )

    # Compare to all haplotype sequences
    for haplotype_id, haplotype_seq in haplotype_groups:
        # Calculate identity using CIGAR-based method
        result = calculate_identity_with_cigar(sequence, haplotype_seq, use_edlib=use_edlib)

        # Select identity metric based on method
        if identity_method == "target_based":
            identity = result['target_identity']
            secondary_identity = result.get('alignment_identity', result['classic_identity'])
            tertiary_identity = result['classic_identity']
        else:  # classic
            identity = result['classic_identity']
            secondary_identity = result.get('alignment_identity', result['target_identity'])
            tertiary_identity = result['target_identity']

        candidate_score = (
            identity,
            secondary_identity,
            tertiary_identity,
            -result['edit_distance'],
        )
        best_score = (
            best_identity,
            best_details['alignment_identity'] if best_details else -1.0,
            best_details['classic_identity'] if best_details else -1.0,
            -best_details['edit_distance'] if best_details else float('-inf'),
        )
        runner_score = (
            runner_up_identity,
            runner_up_details['alignment_identity'] if runner_up_details else -1.0,
            runner_up_details['classic_identity'] if runner_up_details else -1.0,
            -runner_up_details['edit_distance'] if runner_up_details else float('-inf'),
        )

        # Update best and runner-up
        if candidate_score > best_score:
            # Current best becomes runner-up
            runner_up_haplotype = best_haplotype
            runner_up_identity = best_identity
            runner_up_details = best_details
            # New best
            best_haplotype = haplotype_id
            best_identity = identity
            best_details = result  # Save full details for best match
        elif candidate_score > runner_score:
            # Update runner-up
            runner_up_haplotype = haplotype_id
            runner_up_identity = identity
            runner_up_details = result

    # Determine if assignment is below threshold
    if best_identity < min_identity:
        best_haplotype = None  # Unassigned

    # Detect ties (best and runner-up are very close)
    is_tie = False
    best_align_id = best_details['alignment_identity'] if best_details else 0.0
    runner_align_id = runner_up_details['alignment_identity'] if runner_up_details else 0.0
    primary_delta = best_identity - runner_up_identity
    align_delta = best_align_id - runner_align_id
    if (
        best_haplotype is not None
        and runner_up_identity > 0
        and best_identity >= tie_min_identity
        and primary_delta < tie_margin
        and abs(align_delta) < tie_margin
    ):
        is_tie = True

    # Detect low-confidence assignments (barely above threshold)
    is_low_confidence = False
    if best_haplotype is not None and best_identity < (min_identity + 0.03):
        is_low_confidence = True

    # Calculate length discrepancy
    length_discrepancy = abs(best_details['query_length'] - best_details['target_length']) if best_details else 0

    return {
        'best_haplotype': best_haplotype,
        'best_identity': best_identity,
        'target_identity': best_details['target_identity'] if best_details else 0.0,
        'classic_identity': best_details['classic_identity'] if best_details else 0.0,
        'alignment_identity': best_details['alignment_identity'] if best_details else 0.0,
        'identity_method': identity_method,
        'matches': best_details['matches'] if best_details else 0,
        'mismatches': best_details['mismatches'] if best_details else None,
        'insertions': best_details['insertions'] if best_details else None,
        'deletions': best_details['deletions'] if best_details else None,
        'edit_distance': best_details['edit_distance'] if best_details else 0,
        'cigar': best_details['cigar'] if best_details else None,
        'length_discrepancy': length_discrepancy,
        'runner_up_haplotype': runner_up_haplotype,
        'runner_up_identity': runner_up_identity,
        'is_tie': is_tie,
        'is_low_confidence': is_low_confidence
    }


def _assignment_worker(
    task: Tuple[str, Optional[str]],
    haplotype_groups: List[Tuple[str, str]],
    min_identity: float,
    use_edlib: bool,
    identity_method: str = "target_based",
    tie_margin: float = 0.003,
    tie_min_identity: float = 0.95,
) -> Dict[str, Any]:
    """
    Worker function for parallel haplotype assignment.

    Parameters
    ----------
    task : Tuple[str, Optional[str]]
        (processid, sequence) tuple
    haplotype_groups : List[Tuple[str, str]]
        List of (group_id, haplotype_sequence)
    min_identity : float
        Minimum identity threshold
    use_edlib : bool
        Whether to use edlib
    identity_method : str, optional
        Identity calculation method (default: "target_based")
    tie_margin : float, optional
        Maximum allowed difference between best and runner-up identity to be considered a tie.
    tie_min_identity : float, optional
        Minimum best-identity required before we consider calling a tie.


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
            'haplotype_id': None,
            'identity': 0.0,
            'target_identity': 0.0,
            'classic_identity': 0.0,
            'alignment_identity': 0.0,
            'identity_method': identity_method,
            'matches': 0,
            'mismatches': None,
            'insertions': None,
            'deletions': None,
            'edit_distance': 0,
            'cigar': None,
            'length_discrepancy': 0,
            'runner_up_haplotype': None,
            'runner_up_identity': 0.0,
            'is_tie': False,
            'is_low_confidence': False,
            'status': 'no_sequence'
        }

    # Find best match
    result = find_best_haplotype_match(
        sequence=sequence,
        haplotype_groups=haplotype_groups,
        min_identity=min_identity,
        use_edlib=use_edlib,
        identity_method=identity_method,
        tie_margin=tie_margin,
        tie_min_identity=tie_min_identity,
    )

    # Add processid and status
    result['processid'] = processid
    if result['best_haplotype'] is None:
        result['status'] = 'below_threshold'
    elif result['is_tie']:
        result['status'] = 'tie'
    elif result['is_low_confidence']:
        result['status'] = 'low_confidence'
    else:
        result['status'] = 'assigned'

    # Rename for output consistency
    result['haplotype_id'] = result.pop('best_haplotype')
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
      haplotype_id, cluster_seq_sp, cluster_seq_level, cluster_seq_best_identity, cluster_seq_qcov, cluster_seq_top2_delta, n_top_ties
    """


def compute_species_composition(
    assignments_with_metadata: pd.DataFrame,
    species_column: str = 'species_name',
    min_abundance_pct: float = 1.0,
    logger: Optional[logging.Logger] = None
) -> pd.DataFrame:
    """
    Compute species composition for each haplotype.

    Args:
        assignments_with_metadata: DataFrame with haplotype_id and species columns
        species_column: Name of the column containing species information
        min_abundance_pct: Only report species above this percentage (default 1%)
        logger: Optional logger for warnings

    Returns:
        DataFrame with columns:
            - haplotype_id: Name of the haplotype
            - n_total_samples: Total number of samples in this haplotype
            - n_species: Number of distinct species (above min_abundance_pct)
            - species_list: Comma-separated list of all significant species
            - primary_species: Most abundant species
            - primary_species_pct: Percentage of primary species (as string with %)
            - species_composition: Detailed breakdown with counts and percentages
            - is_multi_species: True if >1 species with >min_abundance_pct
            - is_ambiguous: True if no species has >70% abundance
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Check if species column exists
    if species_column not in assignments_with_metadata.columns:
        logger.warning(f"Species column '{species_column}' not found in metadata. Returning empty composition.")
        return pd.DataFrame()

    # Filter to only assigned samples (where haplotype_id is not None/NaN)
    assigned = assignments_with_metadata[
        assignments_with_metadata['haplotype_id'].notna()
    ].copy()

    if len(assigned) == 0:
        logger.warning("No assigned samples found for species composition analysis.")
        return pd.DataFrame()

    # Group by haplotype_id and compute composition
    composition_data = []

    for group, group_df in assigned.groupby('haplotype_id'):
        # Get species counts, handling NaN values
        species_series = group_df[species_column].fillna('Unknown')
        species_counts = species_series.value_counts()
        total_samples = len(group_df)

        # Calculate percentages
        species_pcts = (species_counts / total_samples * 100)

        # Filter species by min_abundance
        significant_species = species_pcts[species_pcts >= min_abundance_pct]

        # Build detailed composition string
        composition_parts = [
            f"{sp}: {species_counts[sp]} ({species_pcts[sp]:.1f}%)"
            for sp in significant_species.index
        ]
        composition_str = "; ".join(composition_parts)

        # Compute flags
        is_multi_species = len(significant_species) > 1
        primary_species = species_counts.index[0]
        primary_pct = species_pcts.iloc[0]
        is_ambiguous = primary_pct < 70.0

        composition_data.append({
            'haplotype_id': group,
            'n_total_samples': total_samples,
            'n_species': len(significant_species),
            'species_list': ", ".join(significant_species.index),
            'primary_species': primary_species,
            'primary_species_pct': f"{primary_pct:.1f}%",
            'species_composition': composition_str,
            'is_multi_species': is_multi_species,
            'is_ambiguous': is_ambiguous
        })

    result_df = pd.DataFrame(composition_data).sort_values('n_total_samples', ascending=False)

    # Log summary
    if len(result_df) > 0:
        multi_count = result_df['is_multi_species'].sum()
        ambiguous_count = result_df['is_ambiguous'].sum()
        logger.info(f"Species composition computed for {len(result_df)} haplotypes")
        if multi_count > 0:
            logger.info(f"  - {multi_count} groups contain multiple species")
        if ambiguous_count > 0:
            logger.info(f"  - {ambiguous_count} groups are taxonomically ambiguous (<70% primary species)")

    return result_df


def assign_haplotypes(
    metadata_path: str,
    fasta_path: str,
    consensus_path: str,
    output_path: str,
    min_identity: float = 0.90,
    n_processes: int = 1,
    diagnostics_path: Optional[str] = None,
    identity_method: str = "target_based",
    tie_margin: float = 0.003,
    tie_min_identity: float = 0.95,
) -> Dict[str, Any]:
    """
    Assign genotype groups to sequences based on consensus matching.

    This is the main orchestration function that:
    1. Loads metadata, raw sequences, and consensus sequences
    2. Extracts processids from FASTA headers
    3. Assigns each sequence to best-matching haplotype
    4. Generates diagnostic reports
    5. Writes updated metadata with haplotype assignments

    Parameters
    ----------
    metadata_path : str
        Path to metadata TSV file (must have 'processid' column)
    fasta_path : str
        Path to raw sequences FASTA file
    consensus_path : str
        Path to consensus sequences FASTA file
    output_path : str
        Path for output TSV with 'haplotype_id' column added
    min_identity : float, optional
        Minimum identity threshold for assignment (default: 0.90)
    n_processes : int, optional
        Number of parallel processes (default: 1)
    diagnostics_path : str, optional
        Path for diagnostics CSV output (default: None, no diagnostics)
    identity_method : str, optional
        Identity calculation method: "target_based" or "classic" (default: "target_based")
        - "target_based": matches / consensus_length (robust to length differences)
        - "classic": 1 - (edit_distance / max_length) (backwards compatibility)
    tie_margin : float, optional
        Maximum allowed difference between best and runner-up identity to be considered a tie (default: 0.003 = 0.3%).
    tie_min_identity : float, optional
        Minimum best-identity required before we consider calling a tie (default: 0.95).


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
    HaplotypeAssignmentError
        If metadata is missing 'processid' column or other errors
    ValueError
        If parameters are invalid

    Examples
    --------
    >>> results = assign_haplotypes(
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
    logger.info("Starting haplotype assignment workflow")
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

    if identity_method not in ["target_based", "classic"]:
        raise ValueError(
            f"identity_method must be 'target_based' or 'classic', "
            f"got '{identity_method}'"
        )

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
            raise HaplotypeAssignmentError(
                "Metadata TSV missing required 'processid' column"
            )
        logger.info(f"Loaded {len(metadata_df)} samples from metadata")
    except Exception as e:
        raise HaplotypeAssignmentError(f"Failed to load metadata: {e}") from e

    # Step 2: Load haplotype sequences
    logger.info(f"Step 2/6: Loading haplotype sequences from {consensus_path}")
    haplotype_records = read_fasta_simple(str(consensus_path))
    haplotype_groups = []
    for header, seq in haplotype_records:
        # Extract haplotype ID (e.g., "haplotype_h7_n381")
        haplotype_id_str = header.split()[0]  # Take first word
        haplotype_groups.append((haplotype_id_str, seq))
    logger.info(f"Loaded {len(haplotype_groups)} haplotypes")

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
    logger.info(f"Step 5/6: Assigning haplotypes (using {n_processes} processes)")
    logger.info(f"Identity calculation method: {identity_method}")
    logger.info(f"Tie margin: {tie_margin}, tie_min_identity: {tie_min_identity}")

    worker_func = partial(
        _assignment_worker,
        haplotype_groups=haplotype_groups,
        min_identity=min_identity,
        use_edlib=use_edlib,
        identity_method=identity_method,
        tie_margin=tie_margin,
        tie_min_identity=tie_min_identity,
    )

    if n_processes > 1:
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(worker_func, tasks)
    else:
        results = list(map(worker_func, tasks))

    logger.info(f"Completed haplotype assignment for {len(results)} samples")

    # Step 6: Collect results and generate outputs
    logger.info(f"Step 6/6: Writing outputs")

    # Create mapping for metadata
    processid_to_group = {}
    for result in results:
        processid_to_group[result['processid']] = result['haplotype_id']

    # Add haplotype_id column to metadata
    # NOTE:
    # - haplotype_id IDs come directly from the consensus FASTA headers
    #   (e.g., "consensus_c1_n84") produced by the dereplication step.
    # - These IDs must remain stable and unmodified so that downstream
    #   taxonomy (Phase 4) and phylogenetic relabeling (Phase 5) can
    #   match them exactly to the consensus_taxonomy.csv table and the
    #   tree tip labels.
    metadata_df['haplotype_id'] = metadata_df['processid'].map(processid_to_group)

    # Write updated metadata
    metadata_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote updated metadata to {output_path}")
        
    # Count samples assigned to each haplotype
    # Filter out None values (unassigned samples)
    assigned_groups = metadata_df[metadata_df['haplotype_id'].notna()]['haplotype_id']
    group_counts = assigned_groups.value_counts().to_dict()

    # Create mapping from old names (consensus_cX) to new names (consensus_cX_nZ)
    # where Z = number of samples assigned to that group
    old_to_new_names = {}
    for group, count in group_counts.items():
        # group is like "consensus_c7" and count is number of samples
        new_name = f"{group}_n{count}"
        old_to_new_names[group] = new_name
        logger.debug(f"Renaming {group} to {new_name} ({count} samples assigned)")

    # Update haplotype names in metadata
    metadata_df['haplotype_id'] = metadata_df['haplotype_id'].map(
        lambda x: old_to_new_names.get(x, x) if pd.notna(x) else x
    )

    # Update haplotype names in results list for diagnostics
    for result in results:
        if result['haplotype_id'] is not None:
            result['haplotype_id'] = old_to_new_names.get(result['haplotype_id'], result['haplotype_id'])
        if result['runner_up_haplotype'] is not None:
            result['runner_up_haplotype'] = old_to_new_names.get(result['runner_up_haplotype'], result['runner_up_haplotype'])

    logger.info(f"Renamed {len(old_to_new_names)} haplotypes with sample counts (_nZ suffix)")

    # Write updated metadata
    metadata_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote updated metadata to {output_path}")

    # Compute and write species composition
    # Try to find species column - check common variations
    species_column = None
    for col_name in ['species_name', 'species', 'species name']:
        if col_name in metadata_df.columns:
            species_column = col_name
            break

    if species_column:
        logger.info("Computing species composition for haplotypes...")
        species_comp_df = compute_species_composition(
            metadata_df,
            species_column=species_column,
            min_abundance_pct=1.0,
            logger=logger
        )

        if not species_comp_df.empty:
            # Determine output path for species composition
            # Write to same directory as diagnostics file if provided, otherwise use output_path parent
            if diagnostics_path:
                diagnostics_path_obj = Path(diagnostics_path)
                organism_name = diagnostics_path_obj.stem.replace('_diagnostics', '')
                composition_path = diagnostics_path_obj.parent / f"{organism_name}_species_composition.csv"
            else:
                output_path_obj = Path(output_path)
                organism_name = output_path_obj.stem.replace('_with_genotypes', '')
                composition_path = output_path_obj.parent / f"{organism_name}_species_composition.csv"

            species_comp_df.to_csv(composition_path, index=False)
            logger.info(f"Wrote species composition to {composition_path}")

            # Log warnings for multi-species or ambiguous groups
            multi_species_groups = species_comp_df[species_comp_df['is_multi_species']]
            if not multi_species_groups.empty:
                logger.warning(
                    f"Found {len(multi_species_groups)} haplotypes with multiple species. "
                    f"See {composition_path} for details."
                )
        else:
            logger.info("No species composition data to write (no assigned samples with species information)")
    else:
        logger.info("Species column not found in metadata. Skipping species composition analysis.")

    # Write diagnostics if requested
    if diagnostics_path:
        diagnostics_path = Path(diagnostics_path)
        with open(diagnostics_path, 'w', newline='') as f:
            fieldnames = [
                'processid', 'haplotype_id', 'n_assigned_to_consensus',
                'identity', 'target_identity', 'classic_identity',
                'identity_method', 'matches', 'mismatches', 'insertions', 'deletions',
                'edit_distance', 'length_discrepancy',
                'runner_up_haplotype', 'runner_up_identity',
                'is_tie', 'is_low_confidence', 'status'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for result in results:
                writer.writerow({
                    'processid': result['processid'],
                    'haplotype_id': result['haplotype_id'] or '',
                    'n_assigned_to_consensus': result.get('n_assigned_to_consensus', ''),
                    'identity': round(result['identity'], 6),
                    'target_identity': round(result['target_identity'], 6),
                    'classic_identity': round(result['classic_identity'], 6),
                    'identity_method': result['identity_method'],
                    'matches': result['matches'],
                    'mismatches': result['mismatches'] if result['mismatches'] is not None else '',
                    'insertions': result['insertions'] if result['insertions'] is not None else '',
                    'deletions': result['deletions'] if result['deletions'] is not None else '',
                    'edit_distance': result['edit_distance'],
                    'length_discrepancy': result['length_discrepancy'],
                    'runner_up_haplotype': result['runner_up_haplotype'] or '',
                    'runner_up_identity': round(result['runner_up_identity'], 6),
                    'is_tie': result['is_tie'],
                    'is_low_confidence': result['is_low_confidence'],
                    'status': result['status']
                })

        logger.info(f"Wrote diagnostics to {diagnostics_path}")

    # Calculate summary statistics
    stats = {
        'total': len(results),
        'assigned': sum(1 for r in results if r['haplotype_id'] is not None),
        'unassigned': sum(1 for r in results if r['haplotype_id'] is None),
        'no_sequence': sum(1 for r in results if r['status'] == 'no_sequence'),
        'below_threshold': sum(1 for r in results if r['status'] == 'below_threshold'),
        'ties': sum(1 for r in results if r['is_tie']),
        'low_confidence': sum(1 for r in results if r['is_low_confidence'])
    }

    # Log summary
    logger.info("=" * 70)
    logger.info("Haplotype assignment summary:")
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
