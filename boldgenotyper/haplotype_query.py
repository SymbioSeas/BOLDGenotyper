"""
Haplotype Query Module

This module enables querying of new COI sequences against previously identified
haplotypes from a completed BOLD analysis. It allows users to determine whether
new sequences (e.g., reference genomes, museum specimens, field samples) match
any known haplotypes without re-running the full pipeline.

Key Features:
- Multi-FASTA support for batch processing
- Local alignment to handle length differences
- Multiple output formats (CSV, JSON, detailed text report)
- Metadata enrichment from previous analysis
- Match quality classification (perfect, high, good, moderate, low)
- No identity threshold - reports all matches ranked by quality

Use Case Example:
-----------------
After completing a BOLD analysis on Sphyrnidae, a user obtains a S. gilberti
reference genome and wants to determine if the COX1 sequence matches any
identified haplotypes. This module performs local alignment against all
haplotypes and reports that it's a 100% match to haplotype H0015 (containing
5 S. lewini samples).

Dependencies:
- BioPython (for sequence alignment)
- pandas (for data management)
- pathlib (for file operations)

Author: Steph Smith (symbioseas@outlook.com)
"""

from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
import json
from datetime import datetime
from dataclasses import dataclass, field
import warnings

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import query_assignment

# Configure logging
logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class AlignmentResult:
    """
    Result of aligning one query to one haplotype (overlap alignment).

    Coordinates are 0-based, end-exclusive, on the cleaned query and the
    ungapped haplotype. The overlap is the span from the first to the last
    column where both sequences have a base; end gaps outside it are free.

    Attributes
    ----------
    query_id, haplotype_id : str
    score : float
        Alignment score
    identity_pct : float
        Matches / columns with a base in both sequences x 100 (gaps excluded)
    matches : int
    aligned_length : int
        Columns with a base in both sequences
    query_length, haplotype_length : int
        Ungapped lengths
    divergence : float
        1 - identity_pct / 100
    aligned_query, aligned_haplotype : str
        Gapped overlap, for display
    n_mismatches : int
    n_gaps : int
        Gap columns inside the overlap (either sequence)
    gap_runs : list of int
        Lengths of internal gap runs
    indel_in_query : bool
        Any internal gap run whose length is not a multiple of 3 (COI is
        protein-coding; such an indel suggests a sequencing error or pseudogene)
    query_start, query_end, haplotype_start, haplotype_end : int
        Overlap coordinates
    query_coverage : float
        Query bases in the overlap / query_length
    overlap_coverage : float
        Query bases in the overlap / min(query_length, haplotype_length); the
        --min-query-coverage threshold applies to this value
    differences : int
        n_mismatches + n_gaps
    haplotype_region : str
        Ungapped haplotype bases in the overlap (used to test whether
        haplotypes are identical over the query's region)
    """
    query_id: str
    haplotype_id: str
    score: float
    identity_pct: float
    matches: int
    aligned_length: int
    query_length: int
    haplotype_length: int
    divergence: float
    aligned_query: str
    aligned_haplotype: str
    n_mismatches: int = 0
    n_gaps: int = 0
    gap_runs: List[int] = field(default_factory=list)
    indel_in_query: bool = False
    query_start: int = 0
    query_end: int = 0
    haplotype_start: int = 0
    haplotype_end: int = 0
    query_coverage: float = 0.0
    overlap_coverage: float = 0.0
    differences: int = 0
    haplotype_region: str = ""


# ============================================================================
# Core Alignment Functions
# ============================================================================

# Overlap alignment: global, with free end gaps on both sequences, so a query
# shorter or longer than the reference core is compared over their whole
# overlap. A gap costs more than several mismatches (6 > 2k for k <= 2), so a
# mismatch near an end is not traded for a gap plus a shifted match.
ALIGNMENT_SCORING = {
    "mode": "overlap (global, free end gaps)",
    "match": 1, "mismatch": -1, "gap_open": -6, "gap_extend": -1, "end_gap": 0,
}


def make_aligner() -> PairwiseAligner:
    """PairwiseAligner configured for query-to-haplotype overlap alignment."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = ALIGNMENT_SCORING["match"]
    aligner.mismatch_score = ALIGNMENT_SCORING["mismatch"]
    aligner.open_gap_score = ALIGNMENT_SCORING["gap_open"]
    aligner.extend_gap_score = ALIGNMENT_SCORING["gap_extend"]
    aligner.end_gap_score = ALIGNMENT_SCORING["end_gap"]   # after the internal scores
    return aligner


def clean_query_sequence(seq: str) -> Tuple[str, bool]:
    """Uppercase a query and remove whitespace and gap characters ('-', '.').

    Returns (clean_sequence, had_gap_characters).
    """
    text = "".join(str(seq).split()).upper()
    had_gaps = "-" in text or "." in text
    return text.replace("-", "").replace(".", ""), had_gaps


def load_query_sequences(query_fasta: Path, min_length: int = 100,
                         max_length: int = 2000) -> List[SeqRecord]:
    """Read, clean and length-filter query sequences (lengths are ungapped)."""
    records = list(SeqIO.parse(str(query_fasta), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {query_fasta}")
    logger.info(f"Loaded {len(records)} query sequence(s)")

    cleaned, n_gapped = [], 0
    for record in records:
        seq, had_gaps = clean_query_sequence(str(record.seq))
        n_gapped += had_gaps
        if len(seq) < min_length:
            logger.warning(f"Query {record.id} too short ({len(seq)} bp < {min_length} bp), skipping")
            continue
        if len(seq) > max_length:
            logger.warning(f"Query {record.id} too long ({len(seq)} bp > {max_length} bp), skipping")
            continue
        cleaned.append(SeqRecord(Seq(seq), id=record.id, description=""))
    if n_gapped:
        logger.warning(
            f"Query input looks aligned: {n_gapped} of {len(records)} sequences contained "
            f"'-' or '.'; these were removed before matching, and lengths are ungapped"
        )
    return cleaned


def align_query_to_haplotype(
    query_seq: str,
    query_id: str,
    haplotype_seq: str,
    haplotype_id: str,
    aligner: Optional[PairwiseAligner] = None,
) -> AlignmentResult:
    """
    Align a query to a haplotype (overlap alignment) and measure the result.

    Every base where the two sequences overlap is compared, so mismatches at
    the ends of a query are counted (v1's local alignment dropped them).

    Examples
    --------
    >>> r = align_query_to_haplotype("ATCGATCG", "query1", "GGATCGATCGGG", "H0001")
    >>> r.identity_pct, r.differences
    (100.0, 0)
    """
    query_clean, _ = clean_query_sequence(query_seq)
    haplotype_clean, _ = clean_query_sequence(haplotype_seq)
    aligner = aligner or make_aligner()
    qlen, hlen = len(query_clean), len(haplotype_clean)

    empty = AlignmentResult(
        query_id=query_id, haplotype_id=haplotype_id, score=0.0, identity_pct=0.0,
        matches=0, aligned_length=0, query_length=qlen, haplotype_length=hlen,
        divergence=1.0, aligned_query="", aligned_haplotype="",
    )
    if not qlen or not hlen:
        return empty
    alignment = aligner.align(query_clean, haplotype_clean)[0]
    coords = alignment.coordinates
    segments = [(int(coords[0][k]), int(coords[1][k]),
                 int(coords[0][k + 1] - coords[0][k]), int(coords[1][k + 1] - coords[1][k]))
                for k in range(coords.shape[1] - 1)]
    blocks = [i for i, (_, _, dq, dh) in enumerate(segments) if dq and dh]
    if not blocks:
        return empty
    first, last = blocks[0], blocks[-1]

    matches = n_mismatches = 0
    gap_runs = []
    aligned_q, aligned_h = [], []
    for i in range(first, last + 1):
        qs, hs, dq, dh = segments[i]
        if dq and dh:
            qpart, hpart = query_clean[qs:qs + dq], haplotype_clean[hs:hs + dh]
            m = sum(a == b for a, b in zip(qpart, hpart))
            matches += m
            n_mismatches += dq - m
            aligned_q.append(qpart)
            aligned_h.append(hpart)
        elif dq:      # query bases opposite a gap in the haplotype (insertion)
            gap_runs.append(dq)
            aligned_q.append(query_clean[qs:qs + dq])
            aligned_h.append("-" * dq)
        elif dh:      # haplotype bases opposite a gap in the query (deletion)
            gap_runs.append(dh)
            aligned_q.append("-" * dh)
            aligned_h.append(haplotype_clean[hs:hs + dh])

    q_start, h_start = segments[first][0], segments[first][1]
    q_end = segments[last][0] + segments[last][2]
    h_end = segments[last][1] + segments[last][3]
    aligned_length = matches + n_mismatches
    identity_pct = matches / aligned_length * 100 if aligned_length else 0.0
    n_gaps = sum(gap_runs)
    query_bases = q_end - q_start

    return AlignmentResult(
        query_id=query_id,
        haplotype_id=haplotype_id,
        score=float(alignment.score),
        identity_pct=identity_pct,
        matches=matches,
        aligned_length=aligned_length,
        query_length=qlen,
        haplotype_length=hlen,
        divergence=1.0 - identity_pct / 100.0,
        aligned_query="".join(aligned_q),
        aligned_haplotype="".join(aligned_h),
        n_mismatches=n_mismatches,
        n_gaps=n_gaps,
        gap_runs=gap_runs,
        indel_in_query=any(run % 3 for run in gap_runs),
        query_start=q_start,
        query_end=q_end,
        haplotype_start=h_start,
        haplotype_end=h_end,
        query_coverage=query_bases / qlen,
        overlap_coverage=min(query_bases / min(qlen, hlen), 1.0),
        differences=n_mismatches + n_gaps,
        haplotype_region=haplotype_clean[h_start:h_end],
    )


def classify_match_quality(identity_pct: float) -> str:
    """
    Deprecated: interpretive identity labels (removed from outputs in v2.0.0).

    The labels ("same species" at 97%, "same genus" at 95%) do not hold across
    taxa; v2 reports numbers only (identity, differences, coverage).
    """
    warnings.warn(
        "classify_match_quality is deprecated and no longer used in outputs; "
        "it will be removed in a future release.",
        DeprecationWarning, stacklevel=2,
    )
    if identity_pct >= 100.0:
        return "perfect"
    elif identity_pct >= 99.5:
        return "high"
    elif identity_pct >= 97.0:
        return "good"
    elif identity_pct >= 95.0:
        return "moderate"
    else:
        return "low"


# ============================================================================
# Metadata Loading
# ============================================================================

def load_haplotype_metadata(analysis_dir: Path) -> Optional[pd.DataFrame]:
    """
    Load haplotype metadata from previous analysis.

    Attempts to load and merge:
    - Haplotype taxonomy (assigned_sp, assignment_level, haplotype_sp)
    - Haplotype stats (n_members, is_singleton, is_suspect)

    Parameters
    ----------
    analysis_dir : Path
        Path to previous analysis directory

    Returns
    -------
    Optional[pd.DataFrame]
        Merged metadata DataFrame with haplotype_id as key, or None if not found

    Examples
    --------
    >>> metadata = load_haplotype_metadata(Path("analysis/Carcharhiniformes_full"))
    >>> metadata.columns
    Index(['haplotype_id', 'assigned_sp', 'assignment_level', 'haplotype_sp', 'n_members', ...])
    """
    analysis_path = Path(analysis_dir)

    stats_file = None
    taxonomy_file = None

    for pattern in ["haplotypes/*_haplotype_stats.csv", "*_haplotype_stats.csv"]:
        matches = list(analysis_path.glob(pattern))
        if matches:
            stats_file = matches[0]
            break

    for pattern in ["taxonomy/*_haplotype_taxonomy.csv", "*_haplotype_taxonomy.csv"]:
        matches = list(analysis_path.glob(pattern))
        if matches:
            taxonomy_file = matches[0]
            break

    if not stats_file and not taxonomy_file:
        logger.warning(
            f"Could not find haplotype stats or taxonomy files in {analysis_dir}. "
            "Metadata enrichment will be skipped."
        )
        return None

    try:
        metadata = None

        if stats_file:
            metadata = pd.read_csv(stats_file)
            logger.info(f"Loaded haplotype stats from {stats_file}")

        if taxonomy_file:
            taxonomy_df = pd.read_csv(taxonomy_file)
            logger.info(f"Loaded haplotype taxonomy from {taxonomy_file}")
            tax_cols = [c for c in ['haplotype_id', 'assigned_sp', 'assignment_level', 'haplotype_sp']
                        if c in taxonomy_df.columns]
            if metadata is not None:
                metadata = metadata.merge(taxonomy_df[tax_cols], on='haplotype_id', how='left')
            else:
                metadata = taxonomy_df[tax_cols]

        return metadata

    except Exception as e:
        logger.error(f"Error loading metadata: {e}")
        return None


# ============================================================================
# Main Query Function
# ============================================================================

def query_against_haplotypes(
    query_fasta: Path,
    haplotype_fasta: Path,
    analysis_dir: Optional[Path] = None,
    top_n: int = 10,
    min_length: int = 100,
    max_length: int = 2000,
    min_identity: float = 0.97,
    tie_margin: Optional[float] = None,
    min_query_coverage: float = 0.95,
    max_mismatch_delta: int = 0,
    min_species_support: float = 0.70,
) -> Tuple[List[AlignmentResult], Optional[pd.DataFrame], List["query_assignment.QueryVerdict"]]:
    """
    Query sequences against haplotype consensuses.

    Main orchestration function that:
    1. Loads query sequences (validates length)
    2. Loads haplotype consensuses
    3. Performs pairwise alignments
    4. Loads metadata (if analysis_dir provided)
    5. Returns ranked results

    Parameters
    ----------
    query_fasta : Path
        Path to query FASTA file (single or multi-FASTA)
    haplotype_fasta : Path
        Path to haplotype consensus FASTA from previous analysis
    analysis_dir : Optional[Path]
        Path to previous analysis directory for metadata enrichment
    top_n : int, optional
        Number of top matches to keep per query (default: 10)
    min_length : int, optional
        Minimum query sequence length (default: 100)
    max_length : int, optional
        Maximum query sequence length (default: 2000)

    Returns
    -------
    Tuple[List[AlignmentResult], Optional[pd.DataFrame]]
        (alignment_results, metadata_df)

    Raises
    ------
    FileNotFoundError
        If input files don't exist
    ValueError
        If no valid sequences found

    Examples
    --------
    >>> results, metadata, verdicts = query_against_haplotypes(
    ...     Path("query.fasta"),
    ...     Path("analysis/haplotypes/Organism_haplotypes.fasta"),
    ...     analysis_dir=Path("analysis/")
    ... )
    >>> len(results)
    50  # 5 queries × 10 top matches each
    """
    # Validate input files
    if not query_fasta.exists():
        raise FileNotFoundError(f"Query file not found: {query_fasta}")
    if not haplotype_fasta.exists():
        raise FileNotFoundError(f"Haplotype file not found: {haplotype_fasta}")

    # Load query sequences: cleaned (gaps removed, uppercase), length-filtered (Q1)
    logger.info(f"Loading query sequences from {query_fasta}")
    valid_queries = load_query_sequences(query_fasta, min_length=min_length, max_length=max_length)
    if not valid_queries:
        raise ValueError("No valid query sequences after length filtering")
    logger.info(f"{len(valid_queries)} query sequence(s) passed length validation")

    # Load haplotype consensuses
    logger.info(f"Loading haplotype consensuses from {haplotype_fasta}")
    haplotype_records = list(SeqIO.parse(haplotype_fasta, "fasta"))

    if not haplotype_records:
        raise ValueError(f"No haplotypes found in {haplotype_fasta}")

    logger.info(f"Loaded {len(haplotype_records)} haplotype(s)")

    # Load metadata if analysis directory provided
    metadata = None
    if analysis_dir:
        metadata = load_haplotype_metadata(analysis_dir)

    # Member labels for species support (Q6): exact and corrected members of
    # each haplotype in a v2 run folder; None otherwise (Q-v1 decision).
    member_labels, label_source = query_assignment.load_member_labels(analysis_dir)
    if member_labels is None:
        logger.warning(f"Species support not computed: {label_source}")
    else:
        logger.info(f"Species support uses {label_source} "
                    f"(members with status assigned or assigned_corrected)")

    # Perform alignments
    logger.info(f"Performing alignments: {len(valid_queries)} queries × {len(haplotype_records)} haplotypes")
    all_results = []
    verdicts = []
    aligner = make_aligner()

    for query_record in valid_queries:
        query_seq = str(query_record.seq)
        query_id = query_record.id

        # Align against all haplotypes
        query_results = []
        for haplotype_record in haplotype_records:
            haplotype_seq = str(haplotype_record.seq)
            haplotype_id = haplotype_record.id

            result = align_query_to_haplotype(
                query_seq, query_id,
                haplotype_seq, haplotype_id,
                aligner=aligner,
            )
            query_results.append(result)

        # Rank (eligible first, then fewer differences, ...) and compute the
        # verdict from the full list before truncating to top N for display.
        query_results = query_assignment.rank_matches(
            query_results, min_identity=min_identity, min_query_coverage=min_query_coverage)
        verdicts.append(query_assignment.assign_query(
            query_id, query_results, member_labels,
            min_identity=min_identity, min_query_coverage=min_query_coverage,
            max_mismatch_delta=max_mismatch_delta, min_species_support=min_species_support,
            tie_margin=tie_margin,
        ))
        all_results.extend(query_results[:top_n])

        # Log best match
        best = query_results[0]
        logger.info(
            f"  {query_id}: Best match = {best.haplotype_id} "
            f"({best.identity_pct:.2f}%, {best.differences} differences, "
            f"coverage {best.overlap_coverage:.2f})"
        )

    logger.info(f"Alignment complete: {len(all_results)} results generated")

    return all_results, metadata, verdicts


# ============================================================================
# Output Formatting
# ============================================================================

def results_to_dataframe(
    results: List[AlignmentResult],
    metadata: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Convert alignment results to pandas DataFrame.

    Parameters
    ----------
    results : List[AlignmentResult]
        List of alignment results
    metadata : Optional[pd.DataFrame]
        Haplotype metadata to merge

    Returns
    -------
    pd.DataFrame
        Results table with one row per alignment
    """
    # Convert to DataFrame
    records = []
    for result in results:
        records.append({
            'query_id': result.query_id,
            'haplotype_id': result.haplotype_id,
            'identity_pct': result.identity_pct,
            'matches': result.matches,
            'aligned_length': result.aligned_length,
            'query_length': result.query_length,
            'haplotype_length': result.haplotype_length,
            'divergence': result.divergence,
            'alignment_score': result.score,
            'n_mismatches': result.n_mismatches,
            'n_gaps': result.n_gaps,
            'differences': result.differences,
            'indel_in_query': result.indel_in_query,
            'query_coverage': round(result.query_coverage, 4),
            'overlap_coverage': round(result.overlap_coverage, 4),
        })

    df = pd.DataFrame(records)

    # Rank within each query, in the order given (rank_matches ordering)
    df['rank'] = df.groupby('query_id').cumcount() + 1

    # Merge with metadata if provided
    if metadata is not None:
        # Merge on haplotype_id
        df = df.merge(
            metadata,
            on='haplotype_id',
            how='left'
        )

    # Reorder columns — species diagnostics immediately after haplotype_id
    core_cols = ['query_id', 'rank', 'haplotype_id', 'assigned_sp', 'haplotype_sp',
                 'assignment_level', 'differences', 'n_mismatches', 'n_gaps',
                 'indel_in_query', 'identity_pct', 'query_coverage', 'overlap_coverage',
                 'matches', 'aligned_length', 'query_length', 'haplotype_length',
                 'divergence', 'alignment_score']
    col_order = [c for c in core_cols if c in df.columns]
    extra_cols = [c for c in df.columns if c not in col_order]
    col_order.extend(extra_cols)

    return df[col_order]


def format_alignment_display(result: AlignmentResult, width: int = 120) -> str:
    """
    Format alignment for human-readable display.

    Parameters
    ----------
    result : AlignmentResult
        Alignment result to format
    width : int, optional
        Maximum line width for wrapping (default: 120)

    Returns
    -------
    str
        Formatted alignment string

    Examples
    --------
    >>> print(format_alignment_display(result))
    Query:     ATCGATCG
    Haplotype: ATCGATCG
    Match:     ||||||||
    """
    aligned_query = result.aligned_query
    aligned_haplotype = result.aligned_haplotype

    # Create match string
    match_str = ''.join(
        '|' if q == h else ' '
        for q, h in zip(aligned_query, aligned_haplotype)
    )

    # Format in blocks
    lines = []
    for i in range(0, len(aligned_query), width):
        query_block = aligned_query[i:i+width]
        haplotype_block = aligned_haplotype[i:i+width]
        match_block = match_str[i:i+width]

        lines.append(f"Query:     {query_block}")
        lines.append(f"Haplotype: {haplotype_block}")
        lines.append(f"Match:     {match_block}")
        lines.append("")

    return '\n'.join(lines)


def write_results(
    results: List[AlignmentResult],
    output_dir: Path,
    metadata: Optional[pd.DataFrame] = None,
    haplotype_file: Path = None,
    analysis_dir: Path = None,
    verdicts=None,
    parameters: Optional[Dict[str, Any]] = None,
    reference_run: Optional[Dict[str, Any]] = None,
    record_matches: Optional[Dict[str, int]] = None,
) -> None:
    """
    Write query results in all output formats.

    Generates:
    - CSV table (query_results.csv)
    - JSON structured data (query_results.json)
    - Detailed text report (query_results_detailed.txt)

    Parameters
    ----------
    results : List[AlignmentResult]
        Alignment results
    output_dir : Path
        Output directory
    metadata : Optional[pd.DataFrame]
        Haplotype metadata
    haplotype_file : Path, optional
        Path to haplotype file (for report header)
    analysis_dir : Path, optional
        Path to analysis directory (for report header)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert to DataFrame
    df = results_to_dataframe(results, metadata)

    # Write CSV
    csv_file = output_dir / "query_results.csv"
    df.to_csv(csv_file, index=False)
    logger.info(f"Wrote CSV results to {csv_file}")

    # Write JSON
    json_file = output_dir / "query_results.json"
    write_json_results(results, json_file, metadata)
    logger.info(f"Wrote JSON results to {json_file}")

    # Write detailed text report
    txt_file = output_dir / "query_results_detailed.txt"
    write_text_report(results, txt_file, metadata, haplotype_file, analysis_dir, verdicts)
    logger.info(f"Wrote detailed report to {txt_file}")

    # Write per-query verdict summary
    if verdicts is not None:
        write_query_summary(verdicts, output_dir, parameters=parameters,
                            reference_run=reference_run, record_matches=record_matches)


QUERY_SUMMARY_COLUMNS = [
    "query_id", "query_length", "aligned_length", "query_coverage", "n_mismatches",
    "n_gaps", "best_identity", "call_level", "equivalence_set", "n_equivalent",
    "set_identical_over_region", "pooled_member_labels", "species_call",
    "species_support", "flags", "n_exact_record_matches", "reference_run_id",
]

POOLED_LABELS_NOTE = ("pooled_member_labels and species_support use the labels of every record "
                      "assigned to the equivalence set's haplotypes, including assigned_corrected "
                      "records; genus-only labels of the leading genus are compatible")


def query_summary_header(parameters: Optional[Dict[str, Any]] = None,
                         reference_run: Optional[Dict[str, Any]] = None) -> List[str]:
    """'#' comment lines written above the query_summary.csv header."""
    lines = [f"# BOLDGenotyper query summary; read with pandas.read_csv(path, comment='#')"]
    if reference_run:
        for key, value in reference_run.items():
            lines.append(f"# reference_run.{key}: {value}")
    else:
        lines.append("# reference_run: unknown (no run_provenance.json)")
    for key, value in (parameters or {}).items():
        lines.append(f"# parameter.{key}: {value}")
    lines.append("# alignment: " + ", ".join(f"{k}={v}" for k, v in ALIGNMENT_SCORING.items()))
    lines.append(f"# note: {POOLED_LABELS_NOTE}")
    return lines


def write_query_summary(verdicts, output_dir, parameters=None, reference_run=None,
                        record_matches=None):
    """Write one row per query to query_summary.csv, below '#' header lines."""
    run_id = (reference_run or {}).get("provenance_md5", "")
    rows = []
    for v in verdicts:
        rows.append({
            "query_id": v.query_id,
            "query_length": v.query_length,
            "aligned_length": v.aligned_length,
            "query_coverage": round(v.query_coverage, 4),
            "n_mismatches": v.n_mismatches,
            "n_gaps": v.n_gaps,
            "best_identity": round(v.best_identity, 4),
            "call_level": v.call_level,
            "equivalence_set": ";".join(v.equivalence_set),
            "n_equivalent": v.n_equivalent,
            "set_identical_over_region": v.set_identical_over_region,
            "pooled_member_labels": v.pooled_member_labels,
            "species_call": v.species_call,
            "species_support": v.species_support,
            "flags": ";".join(v.flags),
            "n_exact_record_matches": (record_matches or {}).get(v.query_id) if record_matches is not None else None,
            "reference_run_id": run_id,
        })
    summary_file = Path(output_dir) / "query_summary.csv"
    with open(summary_file, "w", newline="") as f:
        f.write("\n".join(query_summary_header(parameters, reference_run)) + "\n")
        pd.DataFrame(rows, columns=QUERY_SUMMARY_COLUMNS).to_csv(f, index=False)
    logger.info(f"Wrote query verdict summary: {summary_file}")
    return summary_file


def write_json_results(
    results: List[AlignmentResult],
    output_file: Path,
    metadata: Optional[pd.DataFrame] = None
) -> None:
    """Write results in JSON format."""
    # Group results by query
    queries_dict = {}
    for result in results:
        if result.query_id not in queries_dict:
            queries_dict[result.query_id] = {
                'query_id': result.query_id,
                'query_length': result.query_length,
                'matches': []
            }

        # Get metadata for this haplotype
        hap_metadata = {}
        if metadata is not None and 'haplotype_id' in metadata.columns:
            hap_meta_row = metadata[metadata['haplotype_id'] == result.haplotype_id]
            if not hap_meta_row.empty:
                hap_metadata = hap_meta_row.iloc[0].to_dict()

        match_entry = {
            'rank': len(queries_dict[result.query_id]['matches']) + 1,
            'haplotype_id': result.haplotype_id,
            'identity_pct': result.identity_pct,
            'differences': result.differences,
            'n_mismatches': result.n_mismatches,
            'n_gaps': result.n_gaps,
            'indel_in_query': result.indel_in_query,
            'query_coverage': result.query_coverage,
            'overlap_coverage': result.overlap_coverage,
            'alignment': {
                'matches': result.matches,
                'length': result.aligned_length,
                'score': result.score,
                'query_start': result.query_start,
                'query_end': result.query_end,
                'haplotype_start': result.haplotype_start,
                'haplotype_end': result.haplotype_end,
            },
            'metadata': hap_metadata
        }

        queries_dict[result.query_id]['matches'].append(match_entry)

    # Create output structure
    output = {
        'query_metadata': {
            'timestamp': datetime.now().isoformat(),
            'n_queries': len(queries_dict),
            'n_haplotypes': len(set(r.haplotype_id for r in results))
        },
        'queries': list(queries_dict.values())
    }

    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)


def write_text_report(
    results: List[AlignmentResult],
    output_file: Path,
    metadata: Optional[pd.DataFrame] = None,
    haplotype_file: Path = None,
    analysis_dir: Path = None,
    verdicts=None
) -> None:
    """Write detailed text report with alignments."""
    # Group by query
    queries_dict = {}
    for result in results:
        if result.query_id not in queries_dict:
            queries_dict[result.query_id] = []
        queries_dict[result.query_id].append(result)

    with open(output_file, 'w') as f:
        # Header
        f.write("=" * 80 + "\n")
        f.write("BOLDGenotyper Haplotype Query Results\n")
        f.write("=" * 80 + "\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Queries Processed: {len(queries_dict)}\n")
        f.write(f"Haplotypes Searched: {len(set(r.haplotype_id for r in results))}\n")
        if haplotype_file:
            f.write(f"Haplotype File: {haplotype_file}\n")
        if analysis_dir:
            f.write(f"Reference Analysis: {analysis_dir}\n")
        f.write("\n")

        # Summary
        f.write("Query Summary:\n")
        f.write("-" * 80 + "\n")
        verdict_by_id = {v.query_id: v for v in (verdicts or [])}
        for query_id, query_results in queries_dict.items():
            v = verdict_by_id.get(query_id)
            if v is not None:
                if v.call_level == "none":
                    hap = f"no eligible match (best {v.best_identity:.2f}%, coverage {v.query_coverage:.2f})"
                else:
                    hap = (f"{v.call_level} {';'.join(v.equivalence_set)} "
                           f"({v.differences} difference{'' if v.differences == 1 else 's'}, "
                           f"{v.best_identity:.2f}% identity, "
                           f"coverage {v.query_coverage:.2f})")
                sp = (f"species {v.species_call} (support {v.species_support:.2f})" if v.species_call
                      else f"no species call (support {v.species_support:.2f})" if v.species_support is not None
                      else "no species call")
                f.write(f"{query_id} : {hap}\n")
                f.write(f"{' ' * len(query_id)}   {sp}; pooled labels: {v.pooled_member_labels or '-'}\n")
                if v.flags:
                    f.write(f"{' ' * len(query_id)}   flags: {', '.join(v.flags)}\n")
                continue

            best = query_results[0]
            f.write(
                f"{query_id:30s}: best {best.haplotype_id} "
                f"({best.identity_pct:.2f}%, {best.differences} differences)\n"
            )

        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("Detailed Results\n")
        f.write("=" * 80 + "\n\n")

        # Detailed results for each query
        for query_id, query_results in queries_dict.items():
            best = query_results[0]

            f.write(f"Query: {query_id}\n")
            f.write(f"Best Match: {best.haplotype_id} (Identity: {best.identity_pct:.2f}%)\n")

            # Add metadata if available
            if metadata is not None and 'haplotype_id' in metadata.columns:
                hap_meta = metadata[metadata['haplotype_id'] == best.haplotype_id]
                if not hap_meta.empty:
                    meta_row = hap_meta.iloc[0]
                    species = meta_row.get('assigned_sp', meta_row.get('primary_species', ''))
                    if species:
                        n_members = meta_row.get('n_members', meta_row.get('n_total_samples', 'unknown'))
                        level = meta_row.get('assignment_level', '')
                        level_str = f", {level}-level assignment" if level and level != 'species' else ''
                        f.write(f"Species: {species} ({n_members} samples{level_str})\n")

            f.write("\n")
            f.write(f"Alignment ({best.aligned_length} bp):\n")
            f.write(format_alignment_display(best))
            f.write(f"Matches: {best.matches}/{best.aligned_length} bp ({best.identity_pct:.2f}%)\n")
            f.write(f"Mismatches: {best.n_mismatches}; gap columns: {best.n_gaps}"
                    f"{' (indel not a multiple of 3)' if best.indel_in_query else ''}; "
                    f"query coverage {best.query_coverage:.3f}\n")

            # Show other top matches
            if len(query_results) > 1:
                f.write(f"\nOther Top Matches:\n")
                for i, result in enumerate(query_results[1:], 2):
                    sp_str = ''
                    if metadata is not None and 'haplotype_id' in metadata.columns:
                        hap_meta = metadata[metadata['haplotype_id'] == result.haplotype_id]
                        if not hap_meta.empty:
                            sp = hap_meta.iloc[0].get('assigned_sp', hap_meta.iloc[0].get('primary_species', ''))
                            if sp:
                                sp_str = f" — {sp}"
                    f.write(
                        f"  {i}. {result.haplotype_id}: {result.identity_pct:.2f}%, "
                        f"{result.differences} differences{sp_str}\n"
                    )

            f.write("\n" + "-" * 80 + "\n\n")
