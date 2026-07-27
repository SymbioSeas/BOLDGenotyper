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
haplotypes and reports that it's a 100% match to haplotype_h15_n5 (containing
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
from dataclasses import dataclass

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class AlignmentResult:
    """Result of pairwise alignment between query and haplotype.

    Attributes
    ----------
    query_id : str
        Query sequence identifier
    haplotype_id : str
        Haplotype identifier
    score : float
        Raw alignment score
    identity_pct : float
        Percent identity in aligned region (0-100)
    matches : int
        Number of matching bases
    aligned_length : int
        Length of aligned region
    query_length : int
        Length of query sequence (ungapped)
    haplotype_length : int
        Length of haplotype consensus (ungapped)
    divergence : float
        Distance (1 - identity/100)
    match_quality : str
        Quality classification: perfect, high, good, moderate, low
    aligned_query : str
        Query sequence in alignment
    aligned_haplotype : str
        Haplotype sequence in alignment
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
    match_quality: str
    aligned_query: str
    aligned_haplotype: str


# ============================================================================
# Core Alignment Functions
# ============================================================================

def align_query_to_haplotype(
    query_seq: str,
    query_id: str,
    haplotype_seq: str,
    haplotype_id: str
) -> AlignmentResult:
    """
    Perform local alignment between query and haplotype sequences.

    Uses local alignment (not global) to handle length differences where
    query sequences (e.g., 1557 bp full COX1) are often much longer than
    haplotype consensuses (e.g., 184 bp core region).

    Parameters
    ----------
    query_seq : str
        Query DNA sequence
    query_id : str
        Query sequence identifier
    haplotype_seq : str
        Haplotype consensus sequence
    haplotype_id : str
        Haplotype identifier

    Returns
    -------
    AlignmentResult
        Alignment results with identity metrics

    Notes
    -----
    Uses BioPython's Bio.Align.PairwiseAligner (local mode) with scoring:
    - Match: +2
    - Mismatch: -1
    - Gap open: -2
    - Gap extension: -0.5

    Identity is calculated only in the aligned region, ignoring gaps.

    Examples
    --------
    >>> result = align_query_to_haplotype(
    ...     "ATCGATCG", "query1",
    ...     "ATCG", "haplotype_h1"
    ... )
    >>> result.identity_pct
    100.0
    """
    # Clean sequences (uppercase, remove whitespace)
    query_clean = query_seq.upper().replace(' ', '').replace('\n', '')
    haplotype_clean = haplotype_seq.upper().replace(' ', '').replace('\n', '')

    # Perform local alignment
    _aligner = PairwiseAligner()
    _aligner.mode = 'local'
    _aligner.match_score = 2
    _aligner.mismatch_score = -1
    _aligner.open_gap_score = -2
    _aligner.extend_gap_score = -0.5
    alignments = _aligner.align(query_clean, haplotype_clean)

    if not alignments:
        # No alignment found (should be rare)
        logger.warning(
            f"No alignment found between {query_id} and {haplotype_id}"
        )
        return AlignmentResult(
            query_id=query_id,
            haplotype_id=haplotype_id,
            score=0.0,
            identity_pct=0.0,
            matches=0,
            aligned_length=0,
            query_length=len(query_clean),
            haplotype_length=len(haplotype_clean),
            divergence=1.0,
            match_quality="no_match",
            aligned_query="",
            aligned_haplotype=""
        )

    # Take best alignment
    best_alignment = alignments[0]
    aligned_query = best_alignment[0]
    aligned_haplotype = best_alignment[1]
    score = best_alignment.score

    # Calculate identity in aligned region
    # Count matches where both sequences have a base (not gap)
    matches = 0
    aligned_positions = 0

    for q_base, h_base in zip(aligned_query, aligned_haplotype):
        if q_base != '-' and h_base != '-':
            aligned_positions += 1
            if q_base == h_base:
                matches += 1

    # Calculate metrics
    identity_pct = (matches / aligned_positions * 100) if aligned_positions > 0 else 0.0
    divergence = 1.0 - (identity_pct / 100.0)

    # Classify match quality
    match_quality = classify_match_quality(identity_pct)

    return AlignmentResult(
        query_id=query_id,
        haplotype_id=haplotype_id,
        score=score,
        identity_pct=identity_pct,
        matches=matches,
        aligned_length=aligned_positions,
        query_length=len(query_clean),
        haplotype_length=len(haplotype_clean),
        divergence=divergence,
        match_quality=match_quality,
        aligned_query=aligned_query,
        aligned_haplotype=aligned_haplotype
    )


def classify_match_quality(identity_pct: float) -> str:
    """
    Classify alignment quality based on DNA barcoding standards.

    Parameters
    ----------
    identity_pct : float
        Percent identity (0-100)

    Returns
    -------
    str
        Quality classification

    Notes
    -----
    Classification thresholds:
    - perfect: 100% - exact haplotype match
    - high: ≥99.5% - likely same haplotype, minor sequencing variation
    - good: ≥97% - same species, possibly different haplotype
    - moderate: ≥95% - same genus, divergent haplotype
    - low: <95% - different species or contamination
    """
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
    max_length: int = 2000
) -> Tuple[List[AlignmentResult], Optional[pd.DataFrame]]:
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
    >>> results, metadata = query_against_haplotypes(
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

    # Load query sequences
    logger.info(f"Loading query sequences from {query_fasta}")
    query_records = list(SeqIO.parse(query_fasta, "fasta"))

    if not query_records:
        raise ValueError(f"No sequences found in {query_fasta}")

    logger.info(f"Loaded {len(query_records)} query sequence(s)")

    # Validate and filter query sequences by length
    valid_queries = []
    for record in query_records:
        seq_len = len(str(record.seq))

        if seq_len < min_length:
            logger.warning(
                f"Query {record.id} too short ({seq_len} bp < {min_length} bp), skipping"
            )
            continue

        if seq_len > max_length:
            logger.warning(
                f"Query {record.id} too long ({seq_len} bp > {max_length} bp), skipping"
            )
            continue

        valid_queries.append(record)

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

    # Perform alignments
    logger.info(f"Performing alignments: {len(valid_queries)} queries × {len(haplotype_records)} haplotypes")
    all_results = []

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
                haplotype_seq, haplotype_id
            )
            query_results.append(result)

        # Sort by identity (descending) and keep top N
        query_results.sort(key=lambda x: x.identity_pct, reverse=True)
        all_results.extend(query_results[:top_n])

        # Log best match
        best = query_results[0]
        logger.info(
            f"  {query_id}: Best match = {best.haplotype_id} "
            f"({best.identity_pct:.2f}%, {best.match_quality})"
        )

    logger.info(f"Alignment complete: {len(all_results)} results generated")

    return all_results, metadata


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
            'match_quality': result.match_quality,
        })

    df = pd.DataFrame(records)

    # Add rank within each query
    df['rank'] = df.groupby('query_id')['identity_pct'].rank(ascending=False, method='first').astype(int)

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
                 'assignment_level', 'identity_pct', 'match_quality', 'matches',
                 'aligned_length', 'query_length', 'haplotype_length', 'divergence',
                 'alignment_score']
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
    analysis_dir: Path = None
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
    write_text_report(results, txt_file, metadata, haplotype_file, analysis_dir)
    logger.info(f"Wrote detailed report to {txt_file}")


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
            'match_quality': result.match_quality,
            'alignment': {
                'matches': result.matches,
                'length': result.aligned_length,
                'score': result.score
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
    analysis_dir: Path = None
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
        for query_id, query_results in queries_dict.items():
            best = query_results[0]
            quality_desc = {
                'perfect': 'Perfect match',
                'high': 'High confidence match',
                'good': 'Good match',
                'moderate': 'Moderate match',
                'low': 'Low quality match'
            }.get(best.match_quality, 'Match')

            f.write(
                f"{query_id:30s}: {quality_desc} to {best.haplotype_id} "
                f"({best.identity_pct:.2f}%)\n"
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
            f.write(f"Match Quality: {best.match_quality}\n")

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
                        f"  {i}. {result.haplotype_id}: {result.identity_pct:.2f}% "
                        f"({result.match_quality}){sp_str}\n"
                    )

            f.write("\n" + "-" * 80 + "\n\n")
