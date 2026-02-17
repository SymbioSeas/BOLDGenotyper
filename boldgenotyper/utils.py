"""
Helper Functions and Utilities

This module provides common utility functions used throughout the BOLDGenotyper
package, including file I/O operations, external tool verification, logging
configuration, sequence validation, and taxonomy assignment helpers.

Key Utilities:
1. Logging Configuration
   - Centralized logging setup with console and file output

2. External Tool Management
   - Check for required bioinformatics tools (MAFFT, trimAl, FastTree)
   - Version checking and installation instructions

3. File Operations
   - FASTA reading and writing (pure Python, no BioPython required)
   - Filename sanitization and organism name extraction from BOLD filenames

4. Sequence Validation
   - DNA sequence quality checks (length, N-content, gaps)
   - COI open reading frame validation and orientation detection
   - Core alignment region extraction for variable-length barcodes

5. Taxonomy Assignment
   - Consensus taxonomy for haplotype groups (majority-rule species assignment)
   - Multi-source taxonomy resolution (sequence-based vs metadata-based)

Example Usage:
    >>> from boldgenotyper.utils import check_external_tool, extract_organism_name
    >>> if check_external_tool("mafft", min_version="7.0"):
    ...     print("MAFFT is available")
    >>> organism = extract_organism_name("Sphyrna_lewini_BOLD_data.tsv")
    >>> print(organism)  # Output: Sphyrna_lewini

Author: Steph Smith (symbioseas@outlook.com)
"""

from typing import Optional, Dict, Any, List, Tuple, Union, Iterable
from pathlib import Path
import logging
import subprocess
import re
import sys
import shutil
import os

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
        Name of tool to check (e.g., 'mafft', 'trimal', 'fasttree')
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
    >>> if not check_external_tool("fasttree"):
    ...     print("Please install FastTree")

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
        "fasttree": """
FastTree Installation:
  Via conda: conda install -c bioconda fasttree
  Via apt:   sudo apt-get install fasttree
  Website:   http://www.microbesonline.org/fasttree/
""",
    }

    return instructions.get(
        tool_name.lower(),
        f"Please install {tool_name} and ensure it is in your system PATH"
    )


# ============================================================================
# File I/O and Path Handling
# ============================================================================

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
# Formatting Utilities
# ============================================================================

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


# ============================================================================
# Process ID Extraction (BOLD-specific)
# ============================================================================

def extract_processid_from_header(header: str) -> Optional[str]:
    """
    Extract BOLD processid from FASTA header.

    Uses regex pattern to find processid in common BOLD header formats.

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
# Taxonomy Assignment Helper
# ============================================================================
import pandas as pd

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
    majority_threshold: float = 0.5,
    all_groups: Optional[Iterable[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each consensus group, tally 'species', pick a dominant species if any has
    > majority_threshold of samples; otherwise fall back to genus.

    Returns:
      - assign_table: one row per consensus_group with chosen assigned_sp (+metadata)
      - species_counts: long table of counts per (consensus_group, species)

    Parameters
    ----------
    all_groups : Optional[Iterable[str]]
        Optional complete list of consensus/haplotype IDs that should appear in
        the output even if no samples are assigned (e.g., orphan haplotypes used
        only in phylogenetic trees). Missing groups are added with unassigned
        taxonomy so downstream relabeling can still find a tip label.
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
    def _choose(gname: str, group_df: pd.DataFrame) -> Dict[str, str]:
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
            group_col: gname,  # Use parameter instead of hardcoded name
            "assigned_sp": assigned_sp,
            "assignment_level": level,
            "assignment_notes": notes,
        }

    # Build assignment list by iterating through groups explicitly
    assign_list = []
    for gname, group_df in counts.groupby(group_col):
        result = _choose(gname, group_df)
        assign_list.append(result)

    # Include groups that had zero assigned samples (e.g., orphan haplotypes kept for phylogeny)
    if all_groups is not None:
        observed = set(counts[group_col].unique())
        # Helper: consider variants with trailing _n\d+ as equivalent to base
        def _is_observed_variant(name: str) -> bool:
            if name in observed:
                return True
            return any(obs.startswith(f"{name}_n") for obs in observed)

        for gname in all_groups:
            if _is_observed_variant(gname):
                continue
            assign_list.append({
                group_col: gname,
                "assigned_sp": "",
                "assignment_level": "unassigned",
                "assignment_notes": "no samples assigned; added to retain label",
                "majority_fraction": 0.0,
            })

    assign = pd.DataFrame(assign_list)

    # Merge in majority_fraction (keep original column name), coalescing if it already exists
    assign = assign.merge(
        maj,
        on=group_col,
        how="left",
        suffixes=("", "_maj")
    )
    if "majority_fraction_maj" in assign.columns:
        if "majority_fraction" in assign.columns:
            assign["majority_fraction"] = assign["majority_fraction"].fillna(assign["majority_fraction_maj"])
        else:
            assign["majority_fraction"] = assign["majority_fraction_maj"]
        assign = assign.drop(columns=["majority_fraction_maj"])

    # Ensure column order and types
    expected_cols = [group_col, "assigned_sp", "assignment_level", "assignment_notes", "majority_fraction"]
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


# ============================================================================
# COI Validation and Translation
# ============================================================================

def translate_dna(
    sequence: str,
    genetic_code: int = 2,
    to_stop: bool = False
) -> str:
    """
    Translate DNA sequence to protein using specified genetic code.

    Parameters
    ----------
    sequence : str
        DNA sequence to translate
    genetic_code : int, optional
        NCBI genetic code table number (default: 2 = vertebrate mitochondrial)
        See: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    to_stop : bool, optional
        If True, stop translation at first stop codon (default: False)

    Returns
    -------
    str
        Translated protein sequence

    Examples
    --------
    >>> translate_dna("ATGCGATAA", genetic_code=2)
    'MR*'
    >>> translate_dna("ATGCGATAA", genetic_code=2, to_stop=True)
    'MR'

    Notes
    -----
    Uses Bio.Seq if available for accurate genetic code support.
    Falls back to standard genetic code if BioPython not installed.

    References
    ----------
    NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """
    try:
        from Bio.Seq import Seq
        from Bio.Data import CodonTable

        # Clean sequence to ensure no gaps or invalid characters
        cleaned_seq = sequence.upper().replace('-', '').replace('.', '')
        cleaned_seq = ''.join(c for c in cleaned_seq if c in 'ATGCN')

        if not cleaned_seq:
            return ""

        seq_obj = Seq(cleaned_seq)
        protein = seq_obj.translate(table=genetic_code, to_stop=to_stop)
        return str(protein)

    except ImportError:
        logger.warning("BioPython not available. Using standard genetic code for translation.")
        # Fallback: standard genetic code (not mitochondrial)
        standard_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        sequence = sequence.upper()
        protein = []

        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                aa = standard_code.get(codon, 'X')
                if to_stop and aa == '*':
                    break
                protein.append(aa)

        return ''.join(protein)


def reverse_complement(sequence: str) -> str:
    """
    Compute reverse complement of DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence

    Returns
    -------
    str
        Reverse complement sequence

    Examples
    --------
    >>> reverse_complement("ATGC")
    'GCAT'
    """
    complement = str.maketrans('ATGCN', 'TACGN')
    return sequence.upper().translate(complement)[::-1]


def find_best_orf(
    sequence: str,
    genetic_code: int = 2,
    min_orf_length: int = 150
) -> Dict[str, Any]:
    """
    Find best open reading frame in both orientations and all three frames.

    Tests both forward and reverse complement orientations, and all three
    reading frames (0, 1, 2), returning the combination with the fewest
    internal stops and longest ORF.

    Parameters
    ----------
    sequence : str
        DNA sequence to analyze
    genetic_code : int, optional
        NCBI genetic code table (default: 2 = vertebrate mitochondrial)
    min_orf_length : int, optional
        Minimum ORF length in nucleotides (default: 150)

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
        - 'orientation': 'forward' or 'reverse'
        - 'frame': Reading frame offset (0, 1, or 2)
        - 'orf_length_nt': ORF length in nucleotides
        - 'orf_length_aa': ORF length in amino acids
        - 'orf_coverage': Fraction of sequence covered by ORF
        - 'internal_stops': Number of internal stop codons
        - 'protein': Translated protein sequence
        - 'corrected_sequence': Sequence in correct orientation and frame

    Examples
    --------
    >>> result = find_best_orf("ATGAAATAA", genetic_code=2)
    >>> result['orientation']
    'forward'
    >>> result['internal_stops']
    0

    Notes
    -----
    Internal stops are stop codons before the final codon.
    Mitochondrial genetic code (2) differs from standard code,
    particularly in stop codon usage (TGA codes for Trp, not stop).

    This function tries all 6 possible combinations (2 orientations × 3 frames)
    and selects the best one based on:
    1. Fewest internal stop codons (primary criterion)
    2. Longest ORF length (secondary criterion for ties)

    This is critical for COI barcoding because primer-amplified sequences
    may start at different positions relative to the gene start, requiring
    frame-shift correction.
    """
    results = []

    for orientation, seq in [('forward', sequence), ('reverse', reverse_complement(sequence))]:
        # Try all three reading frames (0, 1, 2)
        for frame in [0, 1, 2]:
            # Apply frame offset
            frame_seq = seq[frame:]

            # Trim to codon boundary to avoid partial codon warnings
            trimmed_seq = frame_seq[:len(frame_seq) - (len(frame_seq) % 3)]

            if not trimmed_seq or len(trimmed_seq) < 3:
                # Sequence too short to contain even one codon
                results.append({
                    'orientation': orientation,
                    'frame': frame,
                    'orf_length_nt': 0,
                    'orf_length_aa': 0,
                    'orf_coverage': 0.0,
                    'internal_stops': 999,  # High penalty for invalid
                    'protein': '',
                    'corrected_sequence': seq
                })
                continue

            # Translate without stopping at first stop
            protein = translate_dna(trimmed_seq, genetic_code=genetic_code, to_stop=False)

            # Count internal stops (stops before the last codon)
            internal_stops = protein[:-1].count('*') if len(protein) > 0 else 0

            # Calculate ORF metrics based on trimmed sequence
            orf_length_aa = len(protein)
            orf_length_nt = len(trimmed_seq)  # Use actual trimmed length
            orf_coverage = orf_length_nt / len(sequence) if len(sequence) > 0 else 0.0

            results.append({
                'orientation': orientation,
                'frame': frame,
                'orf_length_nt': orf_length_nt,
                'orf_length_aa': orf_length_aa,
                'orf_coverage': orf_coverage,
                'internal_stops': internal_stops,
                'protein': protein,
                'corrected_sequence': seq  # Return full sequence in correct orientation
            })

    # Choose best combination: fewest internal stops, then longest ORF
    best = min(results, key=lambda x: (x['internal_stops'], -x['orf_length_nt']))

    return best


def check_orf_quality(
    sequence: str,
    genetic_code: int = 2,
    min_coverage: float = 0.7,
    max_internal_stops: int = 2
) -> Dict[str, Any]:
    """
    Check ORF quality and orientation for COI sequence validation.

    Determines if sequence is likely a valid COI barcode by checking:
    - Open reading frame coverage
    - Number of internal stop codons
    - Correct orientation (forward vs reverse complement)

    Parameters
    ----------
    sequence : str
        DNA sequence to validate
    genetic_code : int, optional
        NCBI genetic code table (default: 2 = vertebrate mitochondrial)
    min_coverage : float, optional
        Minimum ORF coverage fraction (default: 0.7 = 70%)
    max_internal_stops : int, optional
        Maximum internal stop codons allowed (default: 2)

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
        - 'is_valid_orf': Boolean, True if passes quality checks
        - 'orientation': 'forward' or 'reverse'
        - 'needs_revcomp': Boolean, True if should be reverse complemented
        - 'orf_coverage': Fraction of sequence in ORF
        - 'internal_stops': Number of internal stop codons
        - 'corrected_sequence': Sequence in correct orientation
        - 'protein': Translated protein
        - 'failure_reasons': List of reasons if validation failed

    Examples
    --------
    >>> result = check_orf_quality("ATGAAAGCGTAA", genetic_code=2)
    >>> result['is_valid_orf']
    True
    >>> result['needs_revcomp']
    False

    Notes
    -----
    This implements the COI validation approach from Porter & Hajibabaei (2020).
    Sequences that fail validation may be:
    - Non-COI contamination
    - Pseudogenes (NUMTs)
    - Sequencing errors
    - Chimeric sequences

    References
    ----------
    Porter & Hajibabaei (2020). Over 2.5 million COI sequences in GenBank
    and growing. PLOS ONE 15(9): e0238765.
    """
    # Find best ORF
    orf_result = find_best_orf(sequence, genetic_code=genetic_code)

    # Check quality criteria
    failure_reasons = []

    if orf_result['orf_coverage'] < min_coverage:
        failure_reasons.append(
            f"Low ORF coverage: {orf_result['orf_coverage']:.2f} < {min_coverage}"
        )

    if orf_result['internal_stops'] > max_internal_stops:
        failure_reasons.append(
            f"Too many internal stops: {orf_result['internal_stops']} > {max_internal_stops}"
        )

    is_valid_orf = len(failure_reasons) == 0
    needs_revcomp = orf_result['orientation'] == 'reverse'

    return {
        'is_valid_orf': is_valid_orf,
        'orientation': orf_result['orientation'],
        'frame': orf_result.get('frame', 0),
        'needs_revcomp': needs_revcomp,
        'orf_coverage': orf_result['orf_coverage'],
        'internal_stops': orf_result['internal_stops'],
        'corrected_sequence': orf_result['corrected_sequence'],
        'protein': orf_result['protein'],
        'failure_reasons': failure_reasons
    }


# ============================================================================
# Core Region Extraction and Masking
# ============================================================================

def compute_core_region_coverage(
    alignment: List[str],
    min_coverage: float = 0.8
) -> Dict[str, Any]:
    """
    Compute per-position coverage statistics for aligned sequences.

    Identifies which alignment positions are covered by a sufficient
    fraction of sequences, enabling extraction of a "core" region.

    Parameters
    ----------
    alignment : List[str]
        List of aligned sequences (all same length)
    min_coverage : float, optional
        Minimum fraction of sequences required (default: 0.8 = 80%)

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
        - 'coverage_per_position': List of coverage fractions [0.0-1.0]
        - 'core_start': Start position of core region
        - 'core_end': End position of core region (exclusive)
        - 'core_length': Length of core region
        - 'alignment_length': Total alignment length
        - 'n_sequences': Number of sequences

    Examples
    --------
    >>> alignment = ["ATGC--", "AT-CGA", "ATGCGA"]
    >>> result = compute_core_region_coverage(alignment, min_coverage=0.67)
    >>> result['core_length']
    4

    Notes
    -----
    Coverage is defined as the fraction of sequences with a non-gap
    character at each position. Positions with coverage >= min_coverage
    define the core region.

    This addresses the variable-length COI sequence problem in BOLD,
    where sequences range from 150-1550 bp.
    """
    if not alignment or len(alignment) == 0:
        return {
            'coverage_per_position': [],
            'core_start': 0,
            'core_end': 0,
            'core_length': 0,
            'alignment_length': 0,
            'n_sequences': 0
        }

    n_sequences = len(alignment)
    alignment_length = len(alignment[0])

    # Compute coverage per position
    coverage_per_position = []
    for pos in range(alignment_length):
        # Count non-gap characters at this position
        non_gaps = sum(1 for seq in alignment if pos < len(seq) and seq[pos] != '-')
        coverage = non_gaps / n_sequences if n_sequences > 0 else 0.0
        coverage_per_position.append(coverage)

    # Find core region (continuous stretch with >= min_coverage)
    # Use the longest such stretch
    core_start = 0
    core_end = 0
    current_start = None
    longest_length = 0

    for pos, cov in enumerate(coverage_per_position):
        if cov >= min_coverage:
            if current_start is None:
                current_start = pos
        else:
            if current_start is not None:
                # End of a core stretch
                length = pos - current_start
                if length > longest_length:
                    longest_length = length
                    core_start = current_start
                    core_end = pos
                current_start = None

    # Check if core extends to end
    if current_start is not None:
        length = alignment_length - current_start
        if length > longest_length:
            core_start = current_start
            core_end = alignment_length

    return {
        'coverage_per_position': coverage_per_position,
        'core_start': core_start,
        'core_end': core_end,
        'core_length': core_end - core_start,
        'alignment_length': alignment_length,
        'n_sequences': n_sequences
    }


def mask_alignment_gaps(
    alignment: List[str],
    gap_threshold: float = 0.5
) -> List[str]:
    """
    Mask alignment columns with excessive gaps.

    Replaces columns with >gap_threshold gaps with all gaps,
    removing low-information positions from the alignment.

    Parameters
    ----------
    alignment : List[str]
        List of aligned sequences (all same length)
    gap_threshold : float, optional
        Maximum gap fraction before masking (default: 0.5 = 50%)

    Returns
    -------
    List[str]
        Masked alignment with gappy columns replaced by all gaps

    Examples
    --------
    >>> alignment = ["ATGC", "AT--", "ATGC"]
    >>> masked = mask_alignment_gaps(alignment, gap_threshold=0.5)
    >>> masked
    ['ATG-', 'AT--', 'ATG-']

    Notes
    -----
    This is a simple gap-based masking approach. Columns with >50% gaps
    typically provide little phylogenetic information and may introduce noise.
    """
    if not alignment or len(alignment) == 0:
        return alignment

    n_sequences = len(alignment)
    alignment_length = len(alignment[0])

    # Identify columns to mask
    columns_to_mask = []
    for pos in range(alignment_length):
        # Count gaps at this position
        n_gaps = sum(1 for seq in alignment if pos < len(seq) and seq[pos] == '-')
        gap_fraction = n_gaps / n_sequences if n_sequences > 0 else 0.0

        if gap_fraction > gap_threshold:
            columns_to_mask.append(pos)

    # Mask columns
    masked_alignment = []
    for seq in alignment:
        masked_seq = list(seq)
        for pos in columns_to_mask:
            if pos < len(masked_seq):
                masked_seq[pos] = '-'
        masked_alignment.append(''.join(masked_seq))

    logger.debug(f"Masked {len(columns_to_mask)}/{alignment_length} columns with >{gap_threshold:.1%} gaps")

    return masked_alignment


def extract_core_region(
    alignment: List[str],
    headers: List[str],
    min_coverage: float = 0.8,
    gap_threshold: float = 0.5,
    min_core_length: int = 200
) -> Optional[Tuple[List[str], List[str]]]:
    """
    Extract core shared region from variable-length alignment.

    Masks gappy columns (insertion artifacts) before identifying the
    region covered by most sequences. This prevents low-coverage insertion
    columns from fragmenting the core region in DNA barcoding datasets.

    Parameters
    ----------
    alignment : List[str]
        List of aligned sequences (all same length)
    headers : List[str]
        Sequence headers (same order as alignment)
    min_coverage : float, optional
        Minimum sequence coverage for core region (default: 0.8 = 80%)
    gap_threshold : float, optional
        Gap threshold for masking columns before coverage computation.
        Columns with >gap_threshold gaps are masked (default: 0.5 = 50%).
        Set to 1.0 to disable gap masking.
    min_core_length : int, optional
        Minimum core region length (default: 200 bp)

    Returns
    -------
    Optional[Tuple[List[str], List[str]]]
        (core_sequences, headers) if core region found, None otherwise

    Examples
    --------
    >>> alignment = ["ATGCATGC", "AT--ATGC", "ATGCATGC"]
    >>> headers = ["seq1", "seq2", "seq3"]
    >>> core_seqs, core_headers = extract_core_region(alignment, headers)

    Notes
    -----
    Workflow:
    1. Mask gappy columns in full alignment (removes insertion artifacts)
    2. Compute coverage per position on masked alignment
    3. Identify longest core region with >= min_coverage
    4. Extract core region from all sequences
    5. Remove all-gap sequences from core

    This addresses the challenge of COI sequences with variable 5'/3' coverage
    in BOLD datasets. Gap masking before coverage computation prevents
    insertion artifacts (columns with >50% gaps) from fragmenting the core
    region, which is critical for DNA barcoding where true indels are rare.
    """
    # Step 1: Mask gappy columns in full alignment BEFORE computing coverage
    # This prevents insertion artifacts from fragmenting the core region
    masked_alignment = alignment
    if gap_threshold < 1.0:
        masked_alignment = mask_alignment_gaps(alignment, gap_threshold=gap_threshold)
        logger.debug(f"Applied gap masking with threshold {gap_threshold:.1%}")

        # Remove all-gap columns created by masking to prevent fragmentation
        # Identify all-gap columns
        if masked_alignment:
            alignment_length = len(masked_alignment[0])
            all_gap_columns = []
            for pos in range(alignment_length):
                if all(seq[pos] == '-' for seq in masked_alignment):
                    all_gap_columns.append(pos)

            # Remove all-gap columns
            if all_gap_columns:
                filtered_alignment = []
                for seq in masked_alignment:
                    filtered_seq = ''.join(seq[i] for i in range(len(seq)) if i not in all_gap_columns)
                    filtered_alignment.append(filtered_seq)
                masked_alignment = filtered_alignment
                logger.debug(f"Removed {len(all_gap_columns)} all-gap columns after masking")

    # Step 2: Compute coverage on masked alignment
    coverage_result = compute_core_region_coverage(masked_alignment, min_coverage=min_coverage)

    # Check if core region meets minimum length with fallback
    core_length = coverage_result['core_length']
    min_fallback_length = 100  # Absolute minimum

    if core_length < min_fallback_length:
        # Core region too short even for fallback
        logger.error(
            f"Core region too short: {core_length} bp < {min_fallback_length} bp "
            f"absolute minimum. Cannot proceed with haplotype discovery."
        )
        return None
    elif core_length < min_core_length:
        # Between fallback and recommended: warn but proceed
        logger.warning(
            f"Core region shorter than recommended: {core_length} bp "
            f"< {min_core_length} bp. Proceeding with fallback threshold. "
            f"Results may have reduced phylogenetic resolution."
        )

    logger.info(
        f"Identified core region: positions {coverage_result['core_start']}-"
        f"{coverage_result['core_end']} ({core_length} bp)"
    )

    # Step 3: Extract core region from masked alignment
    core_start = coverage_result['core_start']
    core_end = coverage_result['core_end']

    core_sequences = [seq[core_start:core_end] for seq in masked_alignment]

    # Step 4: Remove sequences with insufficient core region coverage
    #
    # This filtering step addresses a critical issue in COI barcoding where primer
    # variation and sequence quality result in variable 5'/3' endpoints. Without
    # this filter, sequences that are identical in their overlap but have different
    # coverage would be treated as separate haplotypes, artificially inflating
    # diversity estimates.
    #
    # **Problem**: In test data (Sphyrnidae), this resulted in 44 zero-divergence
    # haplotype pairs before filtering - sequences that were 100% identical in
    # their overlap but differed only in length.
    #
    # **Solution**: Filter sequences that don't meet 95% of the reference coverage.
    # Combined with the zero-divergence merge step in dereplication.py, this
    # completely eliminates artificial haplotype splitting (44 → 0 problematic pairs).
    #
    # See also: merge_zero_divergence_haplotypes() in dereplication.py
    #
    filtered_core = []
    filtered_headers = []

    # Calculate expected ungapped core length (use 95th percentile as reference)
    ungapped_lengths = []
    for seq in core_sequences:
        ungapped_seq = seq.replace('-', '')
        if len(ungapped_seq) > 0:
            ungapped_lengths.append(len(ungapped_seq))

    if not ungapped_lengths:
        logger.error("No sequences with data in core region")
        return None

    # Use 95th percentile as reference (handles outliers better than max)
    import numpy as np
    reference_length = np.percentile(ungapped_lengths, 95)
    min_coverage_threshold = 0.95  # Require 95% of reference length
    min_required_length = reference_length * min_coverage_threshold

    logger.info(
        f"Core region coverage filter: reference length = {reference_length:.0f} bp, "
        f"minimum required = {min_required_length:.0f} bp ({min_coverage_threshold:.0%})"
    )

    removed_count = 0
    for seq, header in zip(core_sequences, headers):
        ungapped_seq = seq.replace('-', '')
        ungapped_length = len(ungapped_seq)

        if ungapped_length == 0:
            # All-gap sequence
            logger.debug(f"Removed all-gap sequence from core: {header}")
            removed_count += 1
        elif ungapped_length < min_required_length:
            # Insufficient coverage - this prevents identical sequences with different
            # lengths from being treated as unique haplotypes
            logger.debug(
                f"Removed low-coverage sequence from core: {header} "
                f"({ungapped_length:.0f} bp < {min_required_length:.0f} bp minimum)"
            )
            removed_count += 1
        else:
            filtered_core.append(seq)
            filtered_headers.append(header)

    if removed_count > 0:
        logger.info(f"Filtered {removed_count} sequences with insufficient core coverage")

    logger.info(
        f"Core region extracted: {len(filtered_core)} sequences, "
        f"{len(filtered_core[0])} positions"
    )

    return filtered_core, filtered_headers
