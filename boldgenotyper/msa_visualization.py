"""
Multiple Sequence Alignment (MSA) Visualization Module

Generates publication-quality MSA visualizations with phylogenetic ordering,
consensus sequences, and conservation metrics. Designed for integration with
the BOLDGenotyper phylogenetic workflow.

Features:
- Phylogeny-ordered sequences (matches tree topology)
- Chunked display for long alignments
- Nucleotide coloring (Chemistry_NT scheme)
- Consensus sequence and sequence logo
- Conservation bar charts
- Adaptive sizing for different dataset sizes

The module integrates with the phylogenetic workflow to create MSA plots
ordered by tree topology, making it easy to identify conserved regions
and lineage-specific mutations.

Author: BOLDGenotyper Development Team
"""

from pathlib import Path
from typing import List, Optional, Tuple
import logging

from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def create_phylo_ordered_msa(
    alignment_path: Path,
    tree_path: Path,
    output_dir: Path,
    organism: str,
    chunk_size: int = 100,
    max_sequences: int = 50,
    color_scheme: str = "Chemistry_NT",
    show_consensus: bool = True,
    show_logo: bool = True,
    output_formats: List[str] = None
) -> List[Path]:
    """
    Create phylogeny-ordered MSA visualization with chunking.

    Generates multiple MSA plots for long alignments, with sequences
    ordered by phylogenetic tree topology.

    Parameters
    ----------
    alignment_path : Path
        Path to aligned FASTA file
    tree_path : Path
        Path to phylogenetic tree (Newick)
    output_dir : Path
        Output directory for plots
    organism : str
        Organism name for file naming
    chunk_size : int, optional
        Alignment positions per chunk (default: 100)
    max_sequences : int, optional
        Maximum sequences to display (default: 50)
    color_scheme : str, optional
        Color scheme for nucleotides (default: "Nucleotide")
        Options: "Nucleotide", "Purine/Pyrimidine", "Clustal", "Taylor"
    show_consensus : bool, optional
        Show consensus bar (default: True)
    show_logo : bool, optional
        Show sequence logo (default: True)
    output_formats : List[str], optional
        Output formats (default: ["pdf", "png"])

    Returns
    -------
    List[Path]
        Paths to generated MSA plot files

    Notes
    -----
    For alignments >100 sequences, the function will downsample to
    max_sequences phylogenetically diverse representatives.

    For alignments >chunk_size bp, multiple plots are generated:
    - {organism}_MSA_chunk1.pdf (positions 1-100)
    - {organism}_MSA_chunk2.pdf (positions 101-200)
    - etc.

    Examples
    --------
    >>> plots = create_phylo_ordered_msa(
    ...     alignment_path="aligned.fasta",
    ...     tree_path="tree.nwk",
    ...     output_dir="output/",
    ...     organism="Sphyrnidae"
    ... )
    >>> len(plots)
    7  # For 650 bp alignment with 100 bp chunks
    """
    try:
        from pymsaviz import MsaViz
    except ImportError:
        logger.error(
            "pymsaviz is not installed. MSA visualization requires pymsaviz.\n"
            "Install with: pip install pymsaviz>=0.4.0"
        )
        return []

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if output_formats is None:
        output_formats = ["pdf", "png"]

    logger.info(f"Generating MSA visualization for {organism}")

    # Load alignment
    try:
        alignment = AlignIO.read(alignment_path, "fasta")
        logger.info(f"  Loaded alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} bp")
    except Exception as e:
        logger.error(f"Failed to load alignment from {alignment_path}: {e}")
        return []

    # Order by tree
    try:
        ordered_alignment = order_by_tree(alignment, tree_path)
    except Exception as e:
        logger.warning(f"Failed to order by tree: {e}. Using original order.")
        ordered_alignment = alignment

    # Downsample if needed
    if len(ordered_alignment) > max_sequences:
        ordered_alignment = downsample_alignment(
            ordered_alignment,
            max_sequences,
            tree_path
        )
        logger.info(f"  Downsampled to {len(ordered_alignment)} sequences")

    # Create chunks
    chunks = chunk_alignment(ordered_alignment, chunk_size)
    logger.info(f"  Creating {len(chunks)} chunks")

    # Generate plots for each chunk
    output_files = []
    for i, (start, end, chunk_records) in enumerate(chunks, 1):
        for fmt in output_formats:
            output_file = output_dir / f"{organism}_MSA_chunk{i}.{fmt}"

            try:
                # pymsaviz expects alignment records, not a list
                # Convert list of SeqRecords to MultipleSeqAlignment if needed
                if isinstance(chunk_records, list):
                    chunk_aln = MultipleSeqAlignment(chunk_records)
                else:
                    chunk_aln = chunk_records

                # Create MSA plot
                # Note: show_consensus and show_logo parameters control display but
                # pymsaviz v0.5.0 doesn't have separate add_consensus_bar/add_seq_logo methods
                mv = MsaViz(
                    chunk_aln,
                    color_scheme=color_scheme,
                    show_label=True,
                    show_seq_char=True,
                    show_consensus=show_consensus,
                    wrap_length=None,
                    show_count=False
                )

                # Save
                dpi_value = 300 if fmt == "png" else None
                mv.savefig(output_file, dpi=dpi_value)
                output_files.append(output_file)
                logger.info(f"  Saved chunk {i} ({start}-{end} bp): {output_file}")

            except Exception as e:
                logger.error(f"Failed to generate MSA chunk {i}: {e}")
                continue

    return output_files


def order_by_tree(
    alignment: MultipleSeqAlignment,
    tree_path: Path
) -> MultipleSeqAlignment:
    """
    Reorder alignment sequences to match phylogenetic tree tip order.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Input alignment
    tree_path : Path
        Path to phylogenetic tree

    Returns
    -------
    MultipleSeqAlignment
        Reordered alignment

    Raises
    ------
    FileNotFoundError
        If tree file doesn't exist
    ValueError
        If tree and alignment have no common sequences
    """
    if not Path(tree_path).exists():
        raise FileNotFoundError(f"Tree file not found: {tree_path}")

    # Read tree
    tree = Phylo.read(tree_path, "newick")

    # Get tip order (leaf nodes in tree traversal order)
    tip_order = [tip.name.replace("'", "") for tip in tree.get_terminals()]

    # Create mapping
    seq_dict = {rec.id: rec for rec in alignment}

    # Reorder
    ordered_records = []
    for tip in tip_order:
        if tip in seq_dict:
            ordered_records.append(seq_dict[tip])

    # Add any sequences not in tree (shouldn't happen but safe)
    remaining = [rec for rec in alignment if rec.id not in tip_order]
    ordered_records.extend(remaining)

    if not ordered_records:
        raise ValueError("No sequences matched between tree and alignment")

    return MultipleSeqAlignment(ordered_records)


def chunk_alignment(
    alignment: MultipleSeqAlignment,
    chunk_size: int
) -> List[Tuple[int, int, List[SeqRecord]]]:
    """
    Split alignment into chunks for visualization.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Full alignment
    chunk_size : int
        Number of positions per chunk

    Returns
    -------
    List[Tuple[int, int, List[SeqRecord]]]
        List of (start, end, records) tuples

    Examples
    --------
    >>> chunks = chunk_alignment(alignment, 100)
    >>> for start, end, records in chunks:
    ...     print(f"Chunk: {start}-{end} ({len(records)} sequences)")
    """
    aln_length = alignment.get_alignment_length()
    chunks = []

    for start in range(0, aln_length, chunk_size):
        end = min(start + chunk_size, aln_length)

        # Extract chunk from each sequence
        chunk_records = [rec[start:end] for rec in alignment]

        chunks.append((start + 1, end, chunk_records))  # 1-indexed

    return chunks


def downsample_alignment(
    alignment: MultipleSeqAlignment,
    max_sequences: int,
    tree_path: Path
) -> MultipleSeqAlignment:
    """
    Downsample alignment to max_sequences phylogenetically diverse representatives.

    Uses tree to select phylogenetically dispersed sequences rather than
    random sampling. Currently implements simple truncation; phylogenetically-
    informed sampling is planned for future enhancement.

    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Full alignment
    max_sequences : int
        Target number of sequences
    tree_path : Path
        Phylogenetic tree for guided sampling

    Returns
    -------
    MultipleSeqAlignment
        Downsampled alignment

    Notes
    -----
    Current implementation uses simple truncation (first N sequences).
    Future versions will implement phylogenetically-informed downsampling
    to maximize phylogenetic diversity in the sampled set.
    """
    # TODO: Implement phylogenetically-informed downsampling
    # For now, just take first max_sequences
    logger.warning(
        f"Downsampling from {len(alignment)} to {max_sequences} sequences. "
        "Using simple truncation for now."
    )
    return alignment[:max_sequences]


def validate_msa_inputs(
    alignment_path: Path,
    tree_path: Path,
    output_dir: Path
) -> Tuple[bool, Optional[str]]:
    """
    Validate inputs for MSA visualization.

    Parameters
    ----------
    alignment_path : Path
        Path to alignment file
    tree_path : Path
        Path to tree file
    output_dir : Path
        Output directory

    Returns
    -------
    Tuple[bool, Optional[str]]
        (is_valid, error_message)

    Examples
    --------
    >>> valid, error = validate_msa_inputs(aln_path, tree_path, out_dir)
    >>> if not valid:
    ...     print(f"Error: {error}")
    """
    # Check alignment file
    if not alignment_path.exists():
        return False, f"Alignment file not found: {alignment_path}"

    # Check tree file
    if not tree_path.exists():
        return False, f"Tree file not found: {tree_path}"

    # Check if output directory can be created
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        return False, f"Cannot create output directory {output_dir}: {e}"

    return True, None
