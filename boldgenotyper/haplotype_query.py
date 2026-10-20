"""
Haplotype Query Module

This module enables querying of new COI sequences against previously identified
haplotypes from a completed BOLD analysis. It allows users to determine whether
new sequences (e.g., reference genomes, museum specimens, field samples) match
any known haplotypes without re-running the full pipeline.

Key Features:

- Multi-FASTA support for batch processing; aligned input is cleaned (gaps removed)
- Overlap alignment (global, free end gaps): every base where query and
  haplotype overlap is compared, whatever their lengths
- Numbers, not labels: differences (mismatches + gap columns), identity,
  coverage, indel flags
- Equivalence sets: the haplotypes a query cannot tell apart
- Species support from the pooled labels of every member record of the set
- The reference run is identified by the md5 of its run_provenance.json
- Optional search of every QC-passed reference record (--search-records)
- Output files are written by query_outputs (calls, candidates, alignments,
  provenance, README), query_overlays (figures) and query_report (HTML)

Use Case Example:
-----------------
After completing a BOLD analysis on Sphyrnidae, a user obtains a S. gilberti
reference genome and wants to determine if the COX1 sequence matches any
identified haplotypes. This module aligns it against all
haplotypes and reports that it's a 100% match to haplotype H0015 (containing
5 S. lewini samples).

Dependencies:

- BioPython (for sequence alignment)
- pandas (for data management)
- pathlib (for file operations)

Author: Steph Smith (symbioseas@outlook.com)
"""

from typing import Dict, List, Tuple, Optional
from pathlib import Path
import logging
from dataclasses import dataclass, field

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
        Query bases in the overlap / min(query_length, haplotype_length)
    query_overhang : int
        Query bases beyond the haplotype's own ends: before the overlap where
        the haplotype has no bases before it, or after it where the haplotype
        has none after it (e.g. a 5' mini-barcode that starts before the
        haplotype core). They cannot be compared.
    compared_coverage : float
        Query bases in the overlap / min(query_length - query_overhang,
        haplotype_length): the share of what the haplotype could cover that the
        alignment covers. The --min-query-coverage threshold applies to this
        value (2.1.0); it equals overlap_coverage when nothing hangs over
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
    query_overhang: int = 0
    compared_coverage: float = 0.0
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
    # Query bases beyond the haplotype's ends cannot be compared and do not count against
    # compared_coverage; unaligned query bases opposite haplotype bases do (2.1.0).
    overhang = (q_start if h_start == 0 else 0) + (qlen - q_end if h_end == hlen else 0)

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
        query_overhang=overhang,
        compared_coverage=min(query_bases / min(qlen - overhang, hlen), 1.0),
        differences=n_mismatches + n_gaps,
        haplotype_region=haplotype_clean[h_start:h_end],
    )


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

    # v2 layout: one haplotype table (haplotype_label is v1's haplotype_sp)
    v2_table = analysis_path / "haplotypes.csv"
    if v2_table.exists():
        try:
            metadata = pd.read_csv(v2_table)
            metadata['haplotype_sp'] = metadata.get('haplotype_label')
            logger.info(f"Loaded haplotype table from {v2_table}")
            return metadata
        except Exception as e:
            logger.error(f"Error loading {v2_table}: {e}")
            return None

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
    min_query_coverage: float = 0.95,
    max_mismatch_delta: int = 0,
    min_species_support: float = 0.70,
    query_records: Optional[List[SeqRecord]] = None,
    haplotype_taxa: Optional[Dict[str, "query_assignment.HaplotypeTaxon"]] = None,
    min_overlap_length: int = query_assignment.DEFAULT_MIN_OVERLAP_LENGTH,
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
    haplotype_taxa : dict, optional
        haplotype_id -> HaplotypeTaxon for species calls; None loads them from
        analysis_dir's haplotypes.csv (query_assignment.load_haplotype_taxa)

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
    if query_records is not None:
        valid_queries = list(query_records)
    else:
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
    if haplotype_taxa is None:
        haplotype_taxa = query_assignment.load_haplotype_taxa(analysis_dir)
    n_curated = sum(t.source == "curated" for t in (haplotype_taxa or {}).values())
    if n_curated:
        logger.info(f"{n_curated} haplotype(s) have curated taxa; sets containing one are "
                    f"called from the haplotypes' taxa")

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
            query_results, min_identity=min_identity, min_query_coverage=min_query_coverage,
            min_overlap_length=min_overlap_length)
        verdicts.append(query_assignment.assign_query(
            query_id, query_results, member_labels,
            min_identity=min_identity, min_query_coverage=min_query_coverage,
            max_mismatch_delta=max_mismatch_delta, min_species_support=min_species_support,
            haplotype_taxa=haplotype_taxa, min_overlap_length=min_overlap_length,
        ))
        all_results.extend(query_results[:top_n])

        # Log best match
        best = query_results[0]
        logger.info(
            f"  {query_id}: Best match = {best.haplotype_id} "
            f"({best.identity_pct:.2f}%, {best.differences} differences, "
            f"coverage {best.compared_coverage:.2f})"
        )

    logger.info(f"Alignment complete: {len(all_results)} results generated")

    return all_results, metadata, verdicts


# ============================================================================
# Output Formatting
# ============================================================================

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


RECORD_MATCH_COLUMNS = ["query_id", "processid", "match_type", "label", "haplotype_id",
                        "status", "record_length"]


def search_records(queries: List[SeqRecord], records_fasta: Path,
                   record_table: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Exact-match records for each query (Q10, --search-records).

    A record matches when the query occurs in it unchanged
    (query_within_record) or it occurs unchanged in the query
    (record_within_query). Searches every QC-passed record, including records
    that are not assigned to any haplotype.
    """
    records = [(r.id, str(r.seq).upper().replace("-", ""))
               for r in SeqIO.parse(str(records_fasta), "fasta")]
    info = (record_table.drop_duplicates("processid").set_index("processid")
            if record_table is not None else None)
    rows = []
    for q in queries:
        qs = str(q.seq).upper()
        for pid, rs in records:
            if qs in rs:
                match_type = "query_within_record"
            elif rs in qs:
                match_type = "record_within_query"
            else:
                continue
            meta = info.loc[pid] if info is not None and pid in info.index else None
            rows.append({
                "query_id": q.id, "processid": pid, "match_type": match_type,
                "label": meta["label"] if meta is not None else None,
                "haplotype_id": meta["haplotype_id"] if meta is not None else None,
                "status": meta["status"] if meta is not None else None,
                "record_length": len(rs),
            })
    return pd.DataFrame(rows, columns=RECORD_MATCH_COLUMNS)


