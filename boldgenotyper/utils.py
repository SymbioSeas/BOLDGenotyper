"""
Helper Functions and Utilities

This module provides common utility functions used throughout the BOLDGenotyper
package, including file I/O operations, external tool verification, logging
configuration, and general-purpose helper functions.

Key Utilities:
1. File Operations
   - Cross-platform path handling using pathlib
   - Safe file reading and writing
   - Automatic directory creation
   - Organism name extraction from filenames

2. External Tool Management
   - Check for required bioinformatics tools (MAFFT, trimAl, PhyML)
   - Version checking
   - Helpful error messages with installation instructions
   - Path discovery in system PATH

3. Logging Configuration
   - Centralized logging setup
   - Progress tracking for long-running operations
   - Structured error reporting
   - User-friendly status messages

4. Data Validation
   - Input file format validation
   - Required field checking
   - Type conversion with error handling
   - Edge case detection

5. General Helpers
   - String manipulation (organism names with special characters)
   - Numeric utilities (sequence identity calculation)
   - Time and memory tracking
   - Configuration file parsing

Design Principles:
- Functions are small, focused, and well-documented
- All file operations use pathlib for cross-platform compatibility
- Error messages are informative and actionable
- Graceful degradation when optional dependencies are missing

Example Usage:
    >>> from boldgenotyper.utils import check_external_tool, extract_organism_name
    >>> if check_external_tool("mafft", min_version="7.0"):
    ...     print("MAFFT is available")
    >>> organism = extract_organism_name("Sphyrna_lewini_BOLD_data.tsv")
    >>> print(organism)  # Output: Sphyrna_lewini

Author: Steph Smith (steph.smith@unc.edu)
"""

from typing import Optional, Dict, Any, List, Tuple, Union
from pathlib import Path
import logging
import subprocess
import re
import sys
import shutil
import os
from datetime import datetime
from collections import Counter

# Configure module logger
logger = logging.getLogger(__name__)


# ============================================================================
# Logging Configuration
# ============================================================================

def setup_logging(
    log_level: str = "INFO",
    log_file: Optional[str] = None,
    format_string: Optional[str] = None,
) -> logging.Logger:
    """
    Configure logging for BOLDGenotyper.

    Sets up a logger with console and optional file output. The logger uses
    a clear, informative format suitable for both interactive use and log
    file analysis.

    Parameters
    ----------
    log_level : str, optional
        Logging level: DEBUG, INFO, WARNING, ERROR, or CRITICAL (default: INFO)
    log_file : str, optional
        Path to log file. If None, logs only to console (default: None)
    format_string : str, optional
        Custom format string for log messages. If None, uses default format

    Returns
    -------
    logging.Logger
        Configured logger instance

    Examples
    --------
    >>> logger = setup_logging(log_level="DEBUG", log_file="analysis.log")
    >>> logger.info("Starting analysis")

    Notes
    -----
    The default format includes timestamp, level, and message:
    [2025-11-03 10:30:45] INFO: Starting analysis
    """
    # Get root logger for the package
    package_logger = logging.getLogger("boldgenotyper")
    package_logger.setLevel(getattr(logging, log_level.upper()))

    # Remove existing handlers to avoid duplicates
    package_logger.handlers.clear()

    # Default format
    if format_string is None:
        format_string = "[%(asctime)s] %(levelname)s: %(message)s"

    formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level.upper()))
    console_handler.setFormatter(formatter)
    package_logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file is not None:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(getattr(logging, log_level.upper()))
        file_handler.setFormatter(formatter)
        package_logger.addHandler(file_handler)

        package_logger.info(f"Logging to file: {log_file}")

    return package_logger


def log_function_call(func_name: str, **kwargs) -> None:
    """
    Log a function call with its parameters.

    Useful for tracking pipeline execution and debugging.

    Parameters
    ----------
    func_name : str
        Name of the function being called
    **kwargs
        Function parameters to log

    Examples
    --------
    >>> log_function_call("dereplicate_sequences", threshold=0.01, n_sequences=150)
    """
    params = ", ".join(f"{k}={v}" for k, v in kwargs.items())
    logger.debug(f"Calling {func_name}({params})")


# ============================================================================
# External Tool Management
# ============================================================================

def check_external_tool(
    tool_name: str,
    min_version: Optional[str] = None,
) -> bool:
    """
    Check if external bioinformatics tool is available and meets version requirements.

    Searches system PATH for the tool and optionally verifies version compatibility.
    Provides helpful error messages with installation instructions if tool is missing.

    Parameters
    ----------
    tool_name : str
        Name of tool to check (e.g., 'mafft', 'trimal', 'phyml')
    min_version : str, optional
        Minimum required version (e.g., '7.0', '1.4.1')

    Returns
    -------
    bool
        True if tool is available and meets version requirement, False otherwise

    Examples
    --------
    >>> if check_external_tool("mafft", min_version="7.0"):
    ...     print("MAFFT is ready")
    >>> if not check_external_tool("phyml"):
    ...     print("Please install PhyML")

    Notes
    -----
    - Tool must be in system PATH
    - Version checking uses tool-specific commands (--version, -v, etc.)
    - Returns False (with warning) if tool not found or version too old
    """
    # Check if tool exists in PATH
    tool_path = shutil.which(tool_name)

    if tool_path is None:
        logger.warning(f"Tool '{tool_name}' not found in PATH")
        logger.info(get_tool_installation_instructions(tool_name))
        return False

    logger.debug(f"Found {tool_name} at: {tool_path}")

    # If no version requirement, we're done
    if min_version is None:
        return True

    # Get tool version
    current_version = get_tool_version(tool_name)

    if current_version is None:
        logger.warning(f"Could not determine version of {tool_name}")
        return True  # Assume it's OK if we can't check version

    # Compare versions
    if compare_versions(current_version, min_version) < 0:
        logger.warning(
            f"{tool_name} version {current_version} is older than "
            f"required version {min_version}"
        )
        return False

    logger.debug(f"{tool_name} version {current_version} meets requirement (>= {min_version})")
    return True


def get_tool_version(tool_name: str) -> Optional[str]:
    """
    Get version string for an external tool.

    Tries common version flags (--version, -v, -version, version) and parses output.

    Parameters
    ----------
    tool_name : str
        Name of tool

    Returns
    -------
    Optional[str]
        Version string if found, None otherwise

    Examples
    --------
    >>> version = get_tool_version("mafft")
    >>> print(version)  # e.g., "7.490"
    """
    # Common version flags to try
    version_flags = ["--version", "-version", "-v", "version"]

    for flag in version_flags:
        try:
            result = subprocess.run(
                [tool_name, flag],
                capture_output=True,
                text=True,
                timeout=5
            )

            # Combine stdout and stderr (some tools output to stderr)
            output = result.stdout + result.stderr

            # Try to find version number
            version = parse_version_string(output)
            if version:
                return version

        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            continue

    return None


def parse_version_string(text: str) -> Optional[str]:
    """
    Extract version number from tool output.

    Looks for patterns like "v7.490", "version 1.4.1", etc.

    Parameters
    ----------
    text : str
        Tool output text

    Returns
    -------
    Optional[str]
        Version string if found, None otherwise

    Examples
    --------
    >>> parse_version_string("MAFFT v7.490 (2021/Oct/30)")
    '7.490'
    """
    # Common version patterns
    patterns = [
        r'v?(\d+\.\d+\.?\d*)',  # v7.490, 1.4.1, etc.
        r'version\s+(\d+\.\d+\.?\d*)',  # version 7.490
        r'(\d+\.\d+\.?\d*)',  # bare number
    ]

    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return match.group(1)

    return None


def compare_versions(version1: str, version2: str) -> int:
    """
    Compare two version strings.

    Parameters
    ----------
    version1 : str
        First version (e.g., "7.490")
    version2 : str
        Second version (e.g., "7.0")

    Returns
    -------
    int
        Negative if version1 < version2, 0 if equal, positive if version1 > version2

    Examples
    --------
    >>> compare_versions("7.490", "7.0")
    1
    >>> compare_versions("1.2.3", "1.2.4")
    -1
    """
    # Split into components and convert to integers
    def version_tuple(v):
        return tuple(int(x) for x in re.findall(r'\d+', v))

    v1 = version_tuple(version1)
    v2 = version_tuple(version2)

    if v1 < v2:
        return -1
    elif v1 > v2:
        return 1
    else:
        return 0


def get_tool_installation_instructions(tool_name: str) -> str:
    """
    Get installation instructions for missing external tools.

    Provides conda and alternative installation methods for common bioinformatics tools.

    Parameters
    ----------
    tool_name : str
        Name of tool

    Returns
    -------
    str
        Installation instructions

    Examples
    --------
    >>> print(get_tool_installation_instructions("mafft"))
    """
    instructions = {
        "mafft": """
MAFFT Installation:
  Via conda: conda install -c bioconda mafft
  Via apt:   sudo apt-get install mafft
  Via brew:  brew install mafft
  Website:   https://mafft.cbrc.jp/alignment/software/
""",
        "trimal": """
trimAl Installation:
  Via conda: conda install -c bioconda trimal
  Via apt:   sudo apt-get install trimal
  Website:   http://trimal.cgenomics.org/
""",
        "phyml": """
PhyML Installation:
  Via conda: conda install -c bioconda phyml
  Via apt:   sudo apt-get install phyml
  Website:   http://www.atgc-montpellier.fr/phyml/
""",
    }

    return instructions.get(
        tool_name.lower(),
        f"Please install {tool_name} and ensure it is in your system PATH"
    )


# ============================================================================
# File I/O and Path Handling
# ============================================================================

def create_output_directory(output_dir: Union[str, Path]) -> Path:
    """
    Create output directory if it doesn't exist.

    Handles path creation with proper error handling and logging.

    Parameters
    ----------
    output_dir : Union[str, Path]
        Path to output directory

    Returns
    -------
    Path
        Path object for output directory

    Raises
    ------
    OSError
        If directory cannot be created due to permissions or other issues

    Examples
    --------
    >>> output_path = create_output_directory("results/analysis_2025")
    >>> print(output_path.exists())
    True
    """
    path = Path(output_dir)

    try:
        path.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created/verified output directory: {path}")
        return path
    except OSError as e:
        logger.error(f"Failed to create directory {path}: {e}")
        raise


def safe_file_path(
    base_dir: Union[str, Path],
    filename: str,
    extension: Optional[str] = None
) -> Path:
    """
    Create a safe file path with optional extension.

    Ensures directory exists and handles filename sanitization.

    Parameters
    ----------
    base_dir : Union[str, Path]
        Base directory
    filename : str
        Filename (will be sanitized)
    extension : str, optional
        File extension (with or without leading dot)

    Returns
    -------
    Path
        Safe file path

    Examples
    --------
    >>> path = safe_file_path("results", "Sphyrna lewini", ".tsv")
    >>> print(path)
    results/Sphyrna_lewini.tsv
    """
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)

    # Sanitize filename
    safe_name = sanitize_filename(filename)

    # Add extension if provided
    if extension:
        if not extension.startswith('.'):
            extension = '.' + extension
        safe_name += extension

    return base / safe_name


def sanitize_filename(filename: str) -> str:
    """
    Sanitize filename for cross-platform compatibility.

    Replaces spaces with underscores and removes problematic characters.

    Parameters
    ----------
    filename : str
        Original filename

    Returns
    -------
    str
        Sanitized filename

    Examples
    --------
    >>> sanitize_filename("Great White Shark (2025)")
    'Great_White_Shark_2025'
    """
    # Replace spaces with underscores
    safe = filename.replace(' ', '_')

    # Remove problematic characters
    safe = re.sub(r'[^\w\-.]', '_', safe)

    # Remove multiple underscores
    safe = re.sub(r'_+', '_', safe)

    # Remove leading/trailing underscores
    safe = safe.strip('_')

    return safe


def extract_organism_name(file_path: Union[str, Path]) -> str:
    """
    Extract organism name from BOLD TSV filename.

    Handles various naming conventions and special characters. Attempts to
    intelligently extract the organism name from common BOLD filename patterns.

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to BOLD TSV file

    Returns
    -------
    str
        Extracted organism name (cleaned for use in output filenames)

    Examples
    --------
    >>> extract_organism_name("Sphyrna_lewini.tsv")
    'Sphyrna_lewini'
    >>> extract_organism_name("/path/to/Great White Shark_BOLD.tsv")
    'Great_White_Shark'
    >>> extract_organism_name("Carcharodon_carcharias_whiteshark_data.tsv")
    'Carcharodon_carcharias'

    Notes
    -----
    Removes common suffixes: _BOLD, _data, _sequences, etc.
    """
    path = Path(file_path)
    basename = path.stem  # Filename without extension

    # Common suffixes to remove
    suffixes_to_remove = [
        '_BOLD', '_bold',
        '_data', '_Data',
        '_sequences', '_Sequences',
        '_barcode', '_Barcode',
        '_COI', '_coi',
        '_download', '_Download',
    ]

    # Remove suffixes
    cleaned = basename
    for suffix in suffixes_to_remove:
        if cleaned.endswith(suffix):
            cleaned = cleaned[:-len(suffix)]

    # Sanitize for output filenames
    cleaned = sanitize_filename(cleaned)

    # If we ended up with something too short, use original
    if len(cleaned) < 3:
        cleaned = sanitize_filename(basename)

    return cleaned


def read_fasta(fasta_path: Union[str, Path]) -> List[Tuple[str, str]]:
    """
    Read FASTA file and return list of (header, sequence) tuples.

    Pure Python implementation, no BioPython required. Handles multi-line
    sequences and various FASTA formats.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file

    Returns
    -------
    List[Tuple[str, str]]
        List of (header, sequence) tuples

    Raises
    ------
    FileNotFoundError
        If FASTA file doesn't exist
    ValueError
        If file is not valid FASTA format

    Examples
    --------
    >>> records = read_fasta("sequences.fasta")
    >>> for header, seq in records:
    ...     print(f"{header}: {len(seq)} bp")

    Notes
    -----
    - Headers do not include the '>' character
    - Sequences are converted to uppercase
    - Empty lines are ignored
    """
    path = Path(fasta_path)

    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records = []
    header = None
    seq_lines = []

    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')

            if not line:  # Skip empty lines
                continue

            if line.startswith('>'):
                # Save previous record if exists
                if header is not None:
                    sequence = ''.join(seq_lines).upper()
                    records.append((header, sequence))

                # Start new record
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())

        # Save last record
        if header is not None:
            sequence = ''.join(seq_lines).upper()
            records.append((header, sequence))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    logger.debug(f"Read {len(records)} sequences from {path}")
    return records


def write_fasta(
    records: List[Tuple[str, str]],
    output_path: Union[str, Path],
    wrap_width: int = 80
) -> None:
    """
    Write sequences to FASTA file.

    Parameters
    ----------
    records : List[Tuple[str, str]]
        List of (header, sequence) tuples
    output_path : Union[str, Path]
        Output FASTA file path
    wrap_width : int, optional
        Line width for sequence wrapping (default: 80, 0 for no wrapping)

    Examples
    --------
    >>> records = [("seq1", "ATCG"), ("seq2", "GCTA")]
    >>> write_fasta(records, "output.fasta")
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, 'w') as fh:
        for header, sequence in records:
            fh.write(f">{header}\n")

            if wrap_width > 0:
                # Wrap sequence to specified width
                for i in range(0, len(sequence), wrap_width):
                    fh.write(sequence[i:i+wrap_width] + '\n')
            else:
                fh.write(sequence + '\n')

    logger.debug(f"Wrote {len(records)} sequences to {path}")


# ============================================================================
# Validation Functions
# ============================================================================

def validate_fasta_file(fasta_path: Union[str, Path]) -> bool:
    """
    Validate that file is properly formatted FASTA.

    Checks for basic FASTA format compliance without being overly strict.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file

    Returns
    -------
    bool
        True if valid FASTA format, False otherwise

    Examples
    --------
    >>> if validate_fasta_file("sequences.fasta"):
    ...     print("Valid FASTA")
    """
    try:
        records = read_fasta(fasta_path)

        # Check that we have at least one record
        if not records:
            logger.warning(f"No sequences found in {fasta_path}")
            return False

        # Check that sequences are not empty
        for header, seq in records:
            if not seq:
                logger.warning(f"Empty sequence for header: {header}")
                return False

        return True

    except Exception as e:
        logger.warning(f"FASTA validation failed for {fasta_path}: {e}")
        return False


def validate_tsv_file(
    tsv_path: Union[str, Path],
    required_columns: Optional[List[str]] = None
) -> bool:
    """
    Validate that TSV file exists and contains required columns.

    Parameters
    ----------
    tsv_path : Union[str, Path]
        Path to TSV file
    required_columns : List[str], optional
        List of required column names

    Returns
    -------
    bool
        True if valid, False otherwise

    Examples
    --------
    >>> if validate_tsv_file("metadata.tsv", ["processid", "nuc"]):
    ...     print("Valid metadata file")
    """
    path = Path(tsv_path)

    if not path.exists():
        logger.error(f"TSV file not found: {path}")
        return False

    if required_columns is None:
        return True

    # Read first line to check headers
    try:
        with open(path, 'r') as fh:
            header_line = fh.readline().strip()
            headers = header_line.split('\t')

            missing = [col for col in required_columns if col not in headers]

            if missing:
                logger.error(
                    f"TSV file {path} is missing required columns: {missing}"
                )
                return False

            return True

    except Exception as e:
        logger.error(f"Error validating TSV file {path}: {e}")
        return False


def validate_sequence(sequence: str, min_length: int = 100) -> Tuple[bool, str]:
    """
    Validate DNA sequence quality.

    Checks for minimum length, valid characters, and excessive ambiguous bases.

    Parameters
    ----------
    sequence : str
        DNA sequence
    min_length : int, optional
        Minimum sequence length (default: 100)

    Returns
    -------
    Tuple[bool, str]
        (is_valid, reason) - True if valid, False with reason if not

    Examples
    --------
    >>> is_valid, reason = validate_sequence("ATCGATCG", min_length=5)
    >>> if not is_valid:
    ...     print(f"Invalid: {reason}")
    """
    seq = sequence.upper().strip()

    # Check length
    if len(seq) < min_length:
        return False, f"Sequence too short ({len(seq)} < {min_length})"

    # Check for valid characters
    valid_chars = set('ACGTN-')
    invalid = set(seq) - valid_chars
    if invalid:
        return False, f"Invalid characters: {invalid}"

    # Check for excessive Ns
    n_count = seq.count('N')
    n_percent = (n_count / len(seq)) * 100
    if n_percent > 50:
        return False, f"Excessive ambiguous bases ({n_percent:.1f}% N)"

    # Check for excessive gaps (for aligned sequences)
    gap_count = seq.count('-')
    gap_percent = (gap_count / len(seq)) * 100
    if gap_percent > 80:
        return False, f"Excessive gaps ({gap_percent:.1f}%)"

    return True, "Valid"


# ============================================================================
# Sequence Utilities
# ============================================================================

def calculate_sequence_identity(
    seq1: str,
    seq2: str,
    ignore_gaps: bool = True
) -> float:
    """
    Calculate percent identity between two sequences.

    For aligned sequences, calculates identity as the fraction of matching
    positions. Optionally ignores positions with gaps.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence
    ignore_gaps : bool, optional
        Ignore positions with gaps (default: True)

    Returns
    -------
    float
        Percent identity (0.0-1.0)

    Examples
    --------
    >>> identity = calculate_sequence_identity("ATCG", "ATCG")
    >>> print(identity)
    1.0
    >>> identity = calculate_sequence_identity("ATCG", "ATGG")
    >>> print(identity)
    0.75

    Notes
    -----
    Based on the compute_distance function from msa_to_consensus.py
    This is the aligned sequence version (not edit distance)
    """
    s1 = seq1.upper()
    s2 = seq2.upper()

    if len(s1) != len(s2):
        raise ValueError("Sequences must be same length (aligned)")

    matches = 0
    valid_sites = 0

    for a, b in zip(s1, s2):
        if ignore_gaps:
            # Only count positions where both are valid bases
            if a in 'ACGT' and b in 'ACGT':
                valid_sites += 1
                if a == b:
                    matches += 1
        else:
            # Count all positions
            valid_sites += 1
            if a == b:
                matches += 1

    if valid_sites == 0:
        return 0.0

    return matches / valid_sites


def get_sequence_stats(sequence: str) -> Dict[str, Any]:
    """
    Calculate statistics for a DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence

    Returns
    -------
    Dict[str, Any]
        Dictionary with length, GC content, base composition, etc.

    Examples
    --------
    >>> stats = get_sequence_stats("ATCGATCG")
    >>> print(stats['length'])
    8
    >>> print(stats['gc_content'])
    0.5
    """
    seq = sequence.upper().strip()
    length = len(seq)

    # Count bases
    base_counts = Counter(seq)

    # Calculate GC content (excluding N and gaps)
    valid_bases = sum(base_counts.get(b, 0) for b in 'ACGT')
    gc_count = base_counts.get('G', 0) + base_counts.get('C', 0)
    gc_content = gc_count / valid_bases if valid_bases > 0 else 0.0

    return {
        'length': length,
        'gc_content': gc_content,
        'base_counts': dict(base_counts),
        'n_count': base_counts.get('N', 0),
        'n_percent': (base_counts.get('N', 0) / length * 100) if length > 0 else 0,
        'gap_count': base_counts.get('-', 0),
    }


# ============================================================================
# Time and Formatting Utilities
# ============================================================================

def format_elapsed_time(seconds: float) -> str:
    """
    Format elapsed time in human-readable format.

    Parameters
    ----------
    seconds : float
        Elapsed time in seconds

    Returns
    -------
    str
        Formatted time string (e.g., "2h 15m 30s")

    Examples
    --------
    >>> format_elapsed_time(3661)
    '1h 1m 1s'
    >>> format_elapsed_time(45)
    '45s'
    """
    if seconds < 60:
        return f"{seconds:.0f}s"

    minutes = seconds / 60
    if minutes < 60:
        return f"{minutes:.1f}m"

    hours = minutes / 60
    minutes_remainder = minutes % 60

    if hours < 24:
        return f"{int(hours)}h {int(minutes_remainder)}m"

    days = hours / 24
    hours_remainder = hours % 24
    return f"{int(days)}d {int(hours_remainder)}h"


def format_file_size(size_bytes: int) -> str:
    """
    Format file size in human-readable format.

    Parameters
    ----------
    size_bytes : int
        Size in bytes

    Returns
    -------
    str
        Formatted size (e.g., "1.5 MB")

    Examples
    --------
    >>> format_file_size(1536)
    '1.5 KB'
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"


def get_timestamp() -> str:
    """
    Get current timestamp string.

    Returns
    -------
    str
        Timestamp in ISO format

    Examples
    --------
    >>> timestamp = get_timestamp()
    >>> print(timestamp)
    2025-11-03T10:30:45
    """
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


# ============================================================================
# Process ID Extraction (BOLD-specific)
# ============================================================================

def extract_processid_from_header(header: str) -> Optional[str]:
    """
    Extract BOLD processid from FASTA header.

    Uses regex pattern to find processid in common BOLD header formats.
    Based on the regex from consensus_group_to_metadata.py.

    Parameters
    ----------
    header : str
        FASTA header (without '>')

    Returns
    -------
    Optional[str]
        Extracted processid or None if not found

    Examples
    --------
    >>> extract_processid_from_header("Sphyrna_lewini_ANGBF11456-15.COI-5P")
    'ANGBF11456-15'
    >>> extract_processid_from_header("Species_GBMIN12345-20")
    'GBMIN12345-20'

    Notes
    -----
    Pattern: underscore + processid + (dot/space/end)
    Processid is typically the last underscore-separated field before a delimiter
    """
    # Pattern from reference script
    pid_pattern = re.compile(r"_(?P<pid>[^.\s_]+)(?:[.\s]|$)")

    match = pid_pattern.search(header)
    if match:
        return match.group("pid")

    return None


# ============================================================================
# Progress Tracking
# ============================================================================

class ProgressTracker:
    """
    Simple progress tracker for long-running operations.

    Examples
    --------
    >>> tracker = ProgressTracker(total=100, description="Processing")
    >>> for i in range(100):
    ...     tracker.update()
    >>> tracker.finish()
    """

    def __init__(self, total: int, description: str = "Progress"):
        """
        Initialize progress tracker.

        Parameters
        ----------
        total : int
            Total number of items to process
        description : str
            Description of the operation
        """
        self.total = total
        self.description = description
        self.current = 0
        self.start_time = datetime.now()
        self.last_log_percent = 0

    def update(self, n: int = 1) -> None:
        """
        Update progress by n items.

        Parameters
        ----------
        n : int
            Number of items completed
        """
        self.current += n
        percent = (self.current / self.total) * 100

        # Log at 10% intervals
        if percent - self.last_log_percent >= 10:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            rate = self.current / elapsed if elapsed > 0 else 0
            eta = (self.total - self.current) / rate if rate > 0 else 0

            logger.info(
                f"{self.description}: {self.current}/{self.total} "
                f"({percent:.1f}%) - ETA: {format_elapsed_time(eta)}"
            )
            self.last_log_percent = int(percent / 10) * 10

    def finish(self) -> None:
        """Log completion."""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        logger.info(
            f"{self.description} complete: {self.total} items "
            f"in {format_elapsed_time(elapsed)}"
        )
        
# ============================================================================
# Taxonomy Assignment Helper
# ============================================================================
import re
import pandas as pd
from typing import Tuple, Dict

def _norm_species(s: str) -> str:
    if pd.isna(s):
        return ""
    s = str(s).strip()
    s = re.sub(r"\s+", " ", s)
    return s
    
def _to_genus(species: str) -> str:
    species = _norm_species(species)
    return species.split(" ")[0] if species else ""
    
def assign_consensus_taxonomy(
    df: pd.DataFrame,
    group_col: str = "consensus_group",
    species_col: str = "species",
    genus_col: str = "genus",
    majority_threshold: float = 0.5
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each consensus group, tally 'species', pick a dominant species if any has
    > majority_threshold of samples; otherwise fall back to genus.

    Returns:
      - assign_table: one row per consensus_group with chosen assigned_sp (+metadata)
      - species_counts: long table of counts per (consensus_group, species)
    """
    df = df.copy()

    # Normalize species/genus text
    if species_col not in df.columns:
        raise ValueError(f"Column '{species_col}' not found")
    if genus_col not in df.columns:
        # derive genus from species if genus column absent
        df[genus_col] = df[species_col].map(_to_genus)

    df[species_col] = df[species_col].map(_norm_species)
    df[genus_col] = df[genus_col].map(_norm_species)

    # Species tallies per group
    counts = (
        df.groupby([group_col, species_col], dropna=False)
          .size()
          .rename("n")
          .reset_index()
    )
    # n in group
    n_by_group = counts.groupby(group_col)["n"].sum().rename("n_in_group")
    counts = counts.merge(n_by_group, on=group_col, how="left")
    counts["frac"] = counts["n"] / counts["n_in_group"]
    
    # majority fraction per group (max frac among species)
    maj = (
        counts.groupby(group_col, as_index=False)["frac"]
            .max()
            .rename(columns={"frac": "majority_fraction"})
    )
    
        
    # pick winner per group (majority species; else fall back to genus
    def _choose(group_df: pd.DataFrame) -> Dict[str, str]:
        gname = group_df[group_col].iloc[0]
        
        # sort by count desc, then species name for determinism
        gsorted = group_df.sort_values(["n", species_col], ascending=[False, True])
        top = gsorted.iloc[0]
        tie = (gsorted["n"].values == top["n"]).sum() > 1
        
        # majority species?
        if (not tie) and (top["frac"] > majority_threshold) and isinstance(top[species_col], str) and top[species_col]:
            assigned_sp = top[species_col]
            level = "species"
            notes = f"majority {top['frac']:.2f}"
        else:
            # fallback: mode genus among members of this consensus group
            members = (
                df.loc[df[group_col] == gname, genus_col]
                    .dropna().astype(str).str.strip()
            )
            assigned_sp = members.value_counts().idxmax() if not members.empty else ""
            level = "genus" if assigned_sp else "unassigned"
            notes = "tie or not majority; fell back to genus" if assigned_sp else "no genus available"
            
        return {
            "consensus_group": gname,
            "assigned_sp": assigned_sp,
            "assignment_level": level,
            "assignment_notes": notes,
        }

    # pick winner per group (kept as nested function so it can see df/*_col/threshold)
    def _choose_series(g: pd.DataFrame) -> pd.Series:
        d = _choose(g)  # your existing function that returns a dict
        return pd.Series(d, dtype="object")
    
    assign = (
        counts.groupby(group_col, group_keys=False)   # <- prevents group col duplication
              .apply(_choose_series)
              .reset_index(drop=True)
    )
    
    # Merge in majority_fraction
    assign = assign.merge(
        maj.rename(columns={group_col: "consensus_group"}),
        on="consensus_group",
        how="left"
    )
    
    # Ensure column order and types
    expected_cols = ["consensus_group", "assigned_sp", "assignment_level", "assignment_notes", "majority_fraction"]
    missing = [c for c in expected_cols if c not in assign.columns]
    if missing:
        raise RuntimeError(f"assign_consensus_taxonomy: missing columns in result: {missing}")
    
    # Build the species-by-group table (long summary)
    species_counts = counts[[group_col, species_col, "n", "frac", "n_in_group"]].rename(
        columns={species_col: "reported_species"}
    )

    return assign, species_counts
    

def pick_final_group_taxon(
    cluster_sp: str,
    cluster_level: str,
    cluster_id: float,
    cluster_qcov: float,
    majority_sp: str,
    majority_level: str,
    majority_frac: float,
    cfg_taxonomy,
):
    """
    Decide final consensus-group taxon based on sequence (cluster) vs metadata majority.
    
    Returns: (final_sp, final_level, provenance)
        provenance = {"cluster_seq","majority_species","cluster_genus","majority_genus","none"}
    """
    # Threshold checks
    id_ok = (cluster_id or 0) >= cfg_taxonomy.min_identity_pct
    cov_ok = (cluster_qcov or 0) >= cfg_taxonomy.min_query_cov_pct
    maj_ok = (majority_frac or 0) >=cfg_taxonomy.majority_species_threshold
    
    # Prefer high-confidence sequence-based call
    if(cluster_level == "species") and id_ok and cov_ok and (cluster_sp or "").strip():
        return cluster_sp, "species", "cluster_seq"
        
    # Fallback: strong metadata majority at species level
    if (majority_level == "species") and maj_ok and (majority_sp or "").strip():
        return majority_sp, "species", "majority_species"
        
    # Genus-level fallback
    if (cluster_level == "genus") and (cluster_sp or "").strip():
        return cluster_sp, "genus", "cluster_genus"
        
    if (majority_level == "genus") and (majority_sp or "").strip():
        return majority_sp, "genus", "majority_genus"
        
    return "", "unassigned", "none"