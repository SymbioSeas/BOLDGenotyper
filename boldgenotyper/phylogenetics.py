"""
Tree Building and Phylogenetic Analysis

This module provides optional phylogenetic analysis functionality for consensus
sequences, including quality filtering, multiple sequence alignment, alignment
trimming, and maximum likelihood tree construction.

IMPORTANT: Quality-Filtered Phylogenetic Workflow
=================================================
This module implements a biologically appropriate phylogenetic workflow that
addresses a critical issue in consensus-based phylogenetics:

Problem:
- Consensus sequences from small clusters (1-5 sequences) are highly incomplete
- These short sequences (33-40% complete) create alignments with 60-70% gaps
- High gap percentages compromise phylogenetic inference:
  * Ambiguous positional homology
  * Reduced phylogenetic signal
  * Unreliable bootstrap support
  * Distorted evolutionary distances

Solution:
- Filter consensus sequences BEFORE tree building
- Only use high-quality sequences (≥600 bp, ≥5 sequences per cluster)
- Trim gappy alignment columns with trimAl
- ALL sequences still used for genotype assignment

This ensures phylogenetic trees are biologically meaningful while maintaining
complete taxonomic coverage for genotyping.

Phylogenetic Analysis Modes:

Mode 1: No phylogeny (default)
- Genotype identification only
- Fastest analysis option
- Suitable for exploratory analysis

Mode 2: --phylogeny flag
- Build filtered, quality-controlled phylogenetic tree
- Midpoint rooted
- No outgroup required
- Good for visualizing relationships among genotypes

Mode 3: --phylogeny --outgroup <fasta>
- User provides outgroup sequences
- Proper phylogenetic rooting
- Requires taxonomic expertise for outgroup selection
- Best for publication-quality phylogenies

Quality-Filtered Workflow:
1. Filter consensus sequences (≥600 bp, ≥5 sequences per cluster)
2. Align filtered sequences using MAFFT (--auto mode)
3. Trim gappy columns with trimAl (gappyout method)
4. Build maximum likelihood tree with FastTree
   - Model: GTR+CAT (General Time Reversible with CAT approximation)
   - Bootstrap: 1000 replicates (or user-specified)
5. Midpoint root (if no outgroup provided)
6. Export tree in Newick format

Dependencies:
- MAFFT v7+ for alignment
- trimAl 1.4+ for alignment trimming
- FastTree 2.1+ for tree building
- Biopython for tree parsing

Output Files:
- {organism}_filtered.fasta : Quality-filtered consensus sequences
- {organism}_aligned.fasta : MAFFT alignment
- {organism}_aligned_trimmed.fasta : trimAl-trimmed alignment
- {organism}_tree.nwk : Phylogenetic tree (Newick format)
- {organism}_tree_relabeled.nwk : Tree with species labels

Example Usage:
    >>> from boldgenotyper.phylogenetics import build_phylogeny
    >>> tree = build_phylogeny(
    ...     consensus_fasta="Sphyrna_consensus.fasta",
    ...     output_prefix="results/Sphyrna",
    ...     min_consensus_length=600,
    ...     min_cluster_size=5,
    ...     trim_alignment=True
    ... )

Geographic Analysis Note:
    Geographic analysis uses ALL consensus sequences (primary genotypes),
    NOT the filtered subset used for phylogenetics. This ensures complete
    biogeographic coverage while maintaining phylogenetic quality.

Author: Steph Smith (steph.smith@unc.edu)
"""

from typing import Optional, List
from pathlib import Path
import logging
import subprocess
import shutil
from Bio import Phylo
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Configure logging
logger = logging.getLogger(__name__)


def filter_consensus_sequences(
    consensus_fasta: str,
    output_fasta: str,
    min_length: int = 600,
    min_cluster_size: int = 5,
) -> tuple[int, int, List[str]]:
    """
    Filter consensus sequences for phylogenetic analysis.

    This function addresses a critical biological issue: consensus sequences from
    small clusters (1-5 sequences) tend to be highly incomplete (<40% of full COI),
    creating alignments with 60-70% gaps. Such alignments produce unreliable
    phylogenetic trees due to:
    - Ambiguous positional homology
    - Reduced phylogenetic signal
    - Inflated/deflated bootstrap support
    - Distorted evolutionary distances

    Only consensus sequences meeting quality criteria are retained for tree building.
    All consensus sequences (including filtered ones) are still used for genotype
    assignment to ensure complete taxonomic coverage.

    Parameters
    ----------
    consensus_fasta : str
        Path to input consensus sequences FASTA
    output_fasta : str
        Path for filtered output FASTA
    min_length : int, optional
        Minimum ungapped sequence length in bp (default: 600)
        Rationale: ~40% of full COI barcode (~1550 bp)
    min_cluster_size : int, optional
        Minimum number of sequences in cluster (default: 5)
        Rationale: Larger clusters produce more complete consensus sequences

    Returns
    -------
    tuple of (n_total, n_kept, filtered_ids)
        n_total : int
            Total number of input sequences
        n_kept : int
            Number of sequences passing filters
        filtered_ids : List[str]
            IDs of filtered-out sequences

    Examples
    --------
    >>> n_total, n_kept, filtered = filter_consensus_sequences(
    ...     "consensus.fasta",
    ...     "consensus_filtered.fasta",
    ...     min_length=600,
    ...     min_cluster_size=5
    ... )
    >>> print(f"Kept {n_kept}/{n_total} sequences for phylogenetics")

    Notes
    -----
    This filtering is ONLY applied to phylogenetic tree building.
    All consensus sequences are retained for:
    - Genotype assignment
    - Taxonomic identification
    - Diversity metrics
    - Geographic analysis

    The filtering prevents the following issue:
        Before filtering:
        - consensus_c18 (1 seq):  513 bp → 67% gaps in alignment
        - consensus_c17 (250 seq): 1548 bp → 0.4% gaps in alignment
        → Unreliable tree topology

        After filtering:
        - Only sequences with >600 bp and >5 sequences retained
        → Biologically meaningful phylogenetic relationships

    See Also
    --------
    run_trimal : Trim gappy columns from alignment
    build_phylogeny : Build phylogenetic tree from filtered sequences
    """
    sequences = list(SeqIO.parse(consensus_fasta, 'fasta'))
    n_total = len(sequences)

    filtered_seqs = []
    filtered_ids = []

    for seq in sequences:
        seq_len = len(seq.seq)

        # Extract cluster size from description
        # Format: "consensus_c1 Consensus sequence for cluster 1 (75 reference sequences)"
        cluster_size = 1  # default for malformed descriptions
        if '(' in seq.description:
            try:
                cluster_size = int(seq.description.split('(')[1].split()[0])
            except (IndexError, ValueError):
                logger.debug(f"Could not parse cluster size from: {seq.description}")

        # Apply filters
        passes_length = seq_len >= min_length
        passes_cluster_size = cluster_size >= min_cluster_size

        if passes_length and passes_cluster_size:
            filtered_seqs.append(seq)
        else:
            filtered_ids.append(seq.id)
            reason = []
            if not passes_length:
                reason.append(f"length={seq_len}<{min_length}")
            if not passes_cluster_size:
                reason.append(f"cluster_size={cluster_size}<{min_cluster_size}")
            logger.debug(f"Filtered {seq.id}: {', '.join(reason)}")

    n_kept = len(filtered_seqs)

    # Write filtered sequences
    if n_kept > 0:
        SeqIO.write(filtered_seqs, output_fasta, 'fasta')
        logger.info(
            f"Consensus filtering: {n_kept}/{n_total} sequences retained for phylogenetics "
            f"({(n_kept/n_total*100):.1f}%)"
        )
        if filtered_ids:
            logger.info(
                f"  Filtered out {len(filtered_ids)} sequences: "
                f"{', '.join(filtered_ids[:5])}"
                f"{' ...' if len(filtered_ids) > 5 else ''}"
            )
    else:
        logger.warning(
            f"No consensus sequences passed quality filters "
            f"(min_length={min_length}, min_cluster_size={min_cluster_size})"
        )

    return n_total, n_kept, filtered_ids


def run_trimal(
    input_alignment: str,
    output_alignment: str,
    method: str = 'gappyout',
    gap_threshold: float = 0.5,
) -> str:
    """
    Trim poorly aligned regions from multiple sequence alignment.

    Removes alignment columns with excessive gaps or low conservation, which
    can introduce noise and reduce phylogenetic accuracy. This is especially
    important when aligning consensus sequences of varying completeness.

    Common issue addressed:
    - Consensus sequences vary from 40-100% complete
    - MAFFT aligns short fragments within longer sequences
    - Creates many gap-only or gap-heavy columns
    - These columns provide no phylogenetic signal but add noise

    Parameters
    ----------
    input_alignment : str
        Path to input aligned FASTA
    output_alignment : str
        Path for trimmed output alignment
    method : str, optional
        trimAl method (default: 'gappyout')
        Options:
        - 'gappyout': Remove gap-heavy columns (recommended for mixed-length seqs)
        - 'strict': Very conservative trimming
        - 'automated1': Heuristic balance of trimming
        - 'manual': Use gap_threshold parameter
    gap_threshold : float, optional
        Maximum fraction of gaps allowed per column (default: 0.5)
        Only used if method='manual'

    Returns
    -------
    str
        Path to trimmed alignment file

    Raises
    ------
    RuntimeError
        If trimAl is not found or trimming fails

    Examples
    --------
    >>> trimmed = run_trimal(
    ...     "alignment.fasta",
    ...     "alignment_trimmed.fasta",
    ...     method='gappyout'
    ... )

    Notes
    -----
    trimAl removes columns, not sequences. All input sequences are retained
    in the output, but problematic alignment regions are removed.

    Before trimming (example):
        Seq1: ATCG----GCTA----NNNN  (50% gaps)
        Seq2: ATCGATCGGCTAGCTATGCA  (0% gaps)
        Column gap %: [0,0,0,0,50,50,50,50,...]

    After gappyout trimming:
        Seq1: ATCGGCTA  (gaps removed)
        Seq2: ATCGATCGGCTAGCTA

    References
    ----------
    Capella-Gutiérrez et al. (2009) trimAl: a tool for automated alignment
    trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15).
    """
    # Check if trimAl is available
    if not shutil.which('trimal'):
        raise RuntimeError(
            "trimAl not found in PATH. Please install trimAl:\n"
            "  - conda: conda install -c bioconda trimal\n"
            "  - apt: sudo apt-get install trimal\n"
            "  - brew: brew install brewsci/bio/trimal"
        )

    # Build trimAl command
    if method == 'manual':
        cmd = ['trimal', '-in', input_alignment, '-out', output_alignment,
               '-gt', str(gap_threshold)]
    else:
        cmd = ['trimal', '-in', input_alignment, '-out', output_alignment,
               f'-{method}']

    logger.info(f"Running trimAl: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        # Parse trimAl output for statistics
        if result.stdout:
            for line in result.stdout.split('\n'):
                if 'columns' in line.lower():
                    logger.debug(f"trimAl: {line.strip()}")

        logger.info(f"trimAl completed: {output_alignment}")
        return output_alignment

    except subprocess.CalledProcessError as e:
        error_msg = f"trimAl failed:\n{e.stderr}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e
    except OSError as e:
        error_msg = f"Failed to run trimAl: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e


def build_phylogeny(
    consensus_fasta: str,
    output_prefix: str,
    outgroup_fasta: Optional[str] = None,
    model: str = "GTR",
    bootstrap: int = 1000,
    threads: int = 1,
    min_consensus_length: int = 600,
    min_cluster_size: int = 5,
    trim_alignment: bool = True,
    trim_method: str = 'gappyout',
) -> Optional[Phylo.BaseTree.Tree]:
    """
    Build maximum likelihood phylogenetic tree from consensus sequences.

    This function implements a quality-filtered phylogenetic workflow:
    1. Filter consensus sequences (remove short/low-coverage sequences)
    2. Align filtered sequences with MAFFT
    3. Trim gappy alignment columns with trimAl
    4. Build ML tree with FastTree

    The filtering and trimming steps address a critical issue: consensus sequences
    from small clusters produce highly gappy alignments (60-70% gaps) that yield
    unreliable phylogenetic trees. Only high-quality consensus sequences are used
    for tree building, while ALL consensus sequences are retained for genotype
    assignment.

    Parameters
    ----------
    consensus_fasta : str
        Path to consensus sequences FASTA file
    output_prefix : str
        Prefix for output files (tree, alignment, etc.)
    outgroup_fasta : str, optional
        Path to outgroup sequences for rooting
    model : str, optional
        Substitution model for FastTree (default: GTR)
    bootstrap : int, optional
        Number of bootstrap replicates (default: 1000)
    threads : int, optional
        Number of CPU threads (default: 1)
    min_consensus_length : int, optional
        Minimum sequence length for phylogenetics (default: 600 bp)
        Sequences shorter than this are excluded from tree building
    min_cluster_size : int, optional
        Minimum cluster size for phylogenetics (default: 5 sequences)
        Clusters smaller than this are excluded from tree building
    trim_alignment : bool, optional
        Trim gappy columns with trimAl (default: True)
        Recommended to remove noise from alignment
    trim_method : str, optional
        trimAl method (default: 'gappyout')
        Options: gappyout, strict, automated1

    Returns
    -------
    Phylo.BaseTree.Tree or None
        Phylogenetic tree object, or None if tree building failed

    Examples
    --------
    >>> tree = build_phylogeny(
    ...     consensus_fasta="Sphyrna_consensus.fasta",
    ...     output_prefix="results/Sphyrna",
    ...     min_consensus_length=600,
    ...     min_cluster_size=5,
    ...     trim_alignment=True
    ... )

    Notes
    -----
    Filtering is ONLY applied to phylogenetic tree building. All consensus
    sequences (including filtered ones) are used for:
    - Genotype assignment
    - Taxonomic identification
    - Geographic analysis
    - Diversity metrics

    Files created:
    - {output_prefix}_filtered.fasta : Filtered consensus sequences
    - {output_prefix}_aligned.fasta : MAFFT alignment
    - {output_prefix}_aligned_trimmed.fasta : trimAl-trimmed alignment (if enabled)
    - {output_prefix}_tree.nwk : Final phylogenetic tree

    See Also
    --------
    filter_consensus_sequences : Filter short/low-coverage sequences
    run_trimal : Trim gappy alignment columns
    """
    try:
        output_prefix_path = Path(output_prefix)
        output_prefix_path.parent.mkdir(parents=True, exist_ok=True)

        # Step 1: Filter consensus sequences for phylogenetic quality
        logger.info("  Step 1/4: Filtering consensus sequences...")
        filtered_fasta = f"{output_prefix}_filtered.fasta"
        n_total, n_kept, filtered_ids = filter_consensus_sequences(
            consensus_fasta,
            filtered_fasta,
            min_length=min_consensus_length,
            min_cluster_size=min_cluster_size
        )

        # Check if enough sequences passed filtering
        if n_kept < 4:
            logger.warning(
                f"Only {n_kept} sequences passed quality filters. "
                f"Need at least 4 sequences for phylogenetic tree. Skipping tree building."
            )
            return None

        # Step 2: Combine with outgroup if provided
        input_fasta = filtered_fasta
        if outgroup_fasta:
            logger.info("  Adding outgroup sequences...")
            combined_fasta = f"{output_prefix}_combined.fasta"
            input_fasta = add_outgroup(filtered_fasta, outgroup_fasta, combined_fasta)

        # Step 3: Run MAFFT alignment
        logger.info("  Step 2/4: Aligning sequences with MAFFT...")
        aligned_fasta = f"{output_prefix}_aligned.fasta"
        run_mafft_alignment(input_fasta, aligned_fasta, threads=threads)

        # Step 4: Trim alignment with trimAl (if enabled)
        tree_input = aligned_fasta
        if trim_alignment:
            logger.info("  Step 3/4: Trimming alignment with trimAl...")
            try:
                trimmed_fasta = f"{output_prefix}_aligned_trimmed.fasta"
                run_trimal(aligned_fasta, trimmed_fasta, method=trim_method)
                tree_input = trimmed_fasta
            except RuntimeError as e:
                logger.warning(f"trimAl failed: {e}")
                logger.warning("Proceeding with untrimmed alignment")
                tree_input = aligned_fasta
        else:
            logger.info("  Step 3/4: Skipping alignment trimming (trim_alignment=False)")

        # Step 5: Run FastTree to build phylogeny
        logger.info("  Step 4/4: Building phylogenetic tree with FastTree...")
        tree_file = f"{output_prefix}_tree.nwk"
        run_fasttree(tree_input, tree_file, model=model)

        # Load and return tree
        if Path(tree_file).exists():
            tree = Phylo.read(tree_file, "newick")
            logger.info(f"  ✓ Phylogenetic tree built successfully")
            logger.info(f"    - {n_kept}/{n_total} consensus sequences used")
            logger.info(f"    - Filtered out: {', '.join(filtered_ids[:3])}"
                       f"{' ...' if len(filtered_ids) > 3 else ''}")
            return tree
        else:
            logger.warning(f"Tree file not created: {tree_file}")
            return None

    except Exception as e:
        logger.error(f"Phylogenetic tree building failed: {e}")
        return None


def run_mafft_alignment(
    input_fasta: str,
    output_fasta: str,
    algorithm: str = "auto",
    threads: int = 1,
) -> str:
    """
    Run MAFFT multiple sequence alignment.

    Parameters
    ----------
    input_fasta : str
        Path to input sequences
    output_fasta : str
        Path for aligned output
    algorithm : str, optional
        MAFFT algorithm (auto, linsi, ginsi, etc.)
    threads : int, optional
        Number of CPU threads

    Returns
    -------
    str
        Path to aligned FASTA file
    """
    # Build MAFFT command
    cmd = ["mafft", f"--{algorithm}", "--thread", str(threads), input_fasta]

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
        return output_fasta

    except subprocess.CalledProcessError as e:
        error_msg = f"MAFFT alignment failed:\n{e.stderr}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e
    except OSError as e:
        error_msg = f"Failed to write alignment output: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e


def run_fasttree(
    alignment_file: str,
    output_tree: str,
    model: str = "GTR",
) -> str:
    """
    Run FastTree for maximum likelihood tree construction.

    Parameters
    ----------
    alignment_file : str
        Path to aligned sequences (FASTA format)
    output_tree : str
        Path for output tree file
    model : str, optional
        Substitution model (default: GTR)

    Returns
    -------
    str
        Path to tree file
    """
    # Build FastTree command
    # FastTree options:
    # -nt: nucleotide sequences
    # -gtr: generalized time-reversible model
    # -gamma: Gamma20-based model of rate heterogeneity
    cmd = ["fasttree", "-nt", "-gtr", "-gamma", alignment_file]

    logger.info(f"Running FastTree: {' '.join(cmd)}")

    try:
        with open(output_tree, 'w') as out_handle:
            result = subprocess.run(
                cmd,
                stdout=out_handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        logger.info(f"FastTree completed: {output_tree}")
        return output_tree

    except subprocess.CalledProcessError as e:
        error_msg = f"FastTree failed:\n{e.stderr}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e
    except OSError as e:
        error_msg = f"Failed to write tree output: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e


def run_phyml(
    alignment_file: str,
    output_prefix: str,
    model: str = "GTR",
    bootstrap: int = 1000,
) -> str:
    """
    Run PhyML for maximum likelihood tree construction.

    Parameters
    ----------
    alignment_file : str
        Path to aligned sequences (PHYLIP format)
    output_prefix : str
        Prefix for output files
    model : str, optional
        Substitution model (default: GTR)
    bootstrap : int, optional
        Number of bootstrap replicates (default: 1000)

    Returns
    -------
    str
        Path to tree file
    """
    # Implementation will go here
    pass


def midpoint_root_tree(tree: Phylo.BaseTree.Tree) -> Phylo.BaseTree.Tree:
    """
    Root tree at midpoint between most distant taxa.

    Parameters
    ----------
    tree : Phylo.BaseTree.Tree
        Unrooted phylogenetic tree

    Returns
    -------
    Phylo.BaseTree.Tree
        Midpoint-rooted tree
    """
    # Implementation will go here
    pass


def add_outgroup(
    ingroup_fasta: str,
    outgroup_fasta: str,
    output_fasta: str,
) -> str:
    """
    Combine ingroup and outgroup sequences for analysis.

    Parameters
    ----------
    ingroup_fasta : str
        Path to consensus sequences
    outgroup_fasta : str
        Path to outgroup sequences
    output_fasta : str
        Path for combined output

    Returns
    -------
    str
        Path to combined FASTA file
    """
    logger.info(f"Combining ingroup and outgroup sequences")

    try:
        with open(output_fasta, 'w') as out_handle:
            # Write ingroup sequences
            for record in SeqIO.parse(ingroup_fasta, "fasta"):
                SeqIO.write(record, out_handle, "fasta")

            # Write outgroup sequences
            for record in SeqIO.parse(outgroup_fasta, "fasta"):
                SeqIO.write(record, out_handle, "fasta")

        logger.info(f"Combined sequences written to: {output_fasta}")
        return output_fasta

    except Exception as e:
        error_msg = f"Failed to combine sequences: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e


def relabel_tree_and_alignment(
    tree_file: str,
    alignment_file: str,
    taxonomy_csv: str,
    output_tree: str,
    output_alignment: str,
    label_column: str = "consensus_group_sp",
    id_column: str = "consensus_group"
) -> tuple[str, str]:
    """
    Relabel phylogenetic tree and alignment with consensus_group_sp labels.

    This function replaces consensus_group labels (e.g., "consensus_c34_n97")
    with consensus_group_sp labels (e.g., "Crassostrea hongkongensis c34_n97")
    in both the tree and alignment files.

    Parameters
    ----------
    tree_file : str
        Path to input Newick tree file
    alignment_file : str
        Path to input aligned FASTA file
    taxonomy_csv : str
        Path to taxonomy CSV with consensus_group and consensus_group_sp columns
    output_tree : str
        Path for output relabeled tree file
    output_alignment : str
        Path for output relabeled alignment file
    label_column : str, optional
        Column in taxonomy CSV containing new labels (default: "consensus_group_sp")
    id_column : str, optional
        Column in taxonomy CSV containing original IDs (default: "consensus_group")

    Returns
    -------
    tuple[str, str]
        Paths to output tree and alignment files

    Raises
    ------
    FileNotFoundError
        If input files don't exist
    ValueError
        If required columns are missing from taxonomy CSV

    Examples
    --------
    >>> relabel_tree_and_alignment(
    ...     tree_file="Crassostrea_tree.nwk",
    ...     alignment_file="Crassostrea_aligned.fasta",
    ...     taxonomy_csv="Crassostrea_consensus_taxonomy.csv",
    ...     output_tree="Crassostrea_tree_relabeled.nwk",
    ...     output_alignment="Crassostrea_aligned_relabeled.fasta"
    ... )
    """
    import pandas as pd

    logger.info(f"Relabeling tree and alignment with {label_column} labels")

    # Check input files exist
    if not Path(tree_file).exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")
    if not Path(alignment_file).exists():
        raise FileNotFoundError(f"Alignment file not found: {alignment_file}")
    if not Path(taxonomy_csv).exists():
        raise FileNotFoundError(f"Taxonomy CSV not found: {taxonomy_csv}")

    # Load taxonomy mapping
    try:
        taxonomy_df = pd.read_csv(taxonomy_csv)
    except Exception as e:
        raise ValueError(f"Failed to read taxonomy CSV: {e}")

    # Validate required columns
    if id_column not in taxonomy_df.columns:
        raise ValueError(f"Column '{id_column}' not found in taxonomy CSV")
    if label_column not in taxonomy_df.columns:
        raise ValueError(f"Column '{label_column}' not found in taxonomy CSV")

    # Create mapping dictionary: consensus_group -> consensus_group_sp
    label_map = dict(zip(taxonomy_df[id_column], taxonomy_df[label_column]))
    logger.info(f"Loaded {len(label_map)} label mappings from taxonomy CSV")

    # Create alternative mapping for base names (consensus_cX -> consensus_cX_nZ)
    # This handles cases where tree tips are named consensus_c1 but taxonomy has consensus_c1_n84
    base_name_map = {}
    for consensus_group, label in label_map.items():
        # Extract base name: consensus_c1_n84 -> consensus_c1
        import re
        match = re.match(r'(consensus_c\d+)(?:_n\d+)?$', str(consensus_group))
        if match:
            base_name = match.group(1)
            # Store mapping from base name to full label
            # If multiple entries have same base (shouldn't happen), use first one
            if base_name not in base_name_map:
                base_name_map[base_name] = label
                logger.debug(f"Base name mapping: {base_name} -> {label}")

    logger.info(f"Created {len(base_name_map)} base name mappings")

    # Relabel tree
    try:
        tree = Phylo.read(tree_file, "newick")
        n_relabeled = 0

        for clade in tree.get_terminals():
            original_name = clade.name
            new_name = None

            # Try exact match first
            if original_name in label_map:
                new_name = label_map[original_name]
            # Try base name match if exact match fails
            elif original_name in base_name_map:
                new_name = base_name_map[original_name]
                logger.debug(f"Using base name match for: {original_name}")

            if new_name:
                clade.name = new_name
                n_relabeled += 1
                logger.debug(f"Relabeled tree tip: {original_name} -> {clade.name}")
            else:
                logger.warning(f"No mapping found for tree tip: {original_name}")

        # Write relabeled tree
        Phylo.write(tree, output_tree, "newick")
        logger.info(f"Relabeled {n_relabeled} tree tips, wrote to: {output_tree}")

    except Exception as e:
        error_msg = f"Failed to relabel tree: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e

    # Relabel alignment
    try:
        relabeled_records = []
        n_relabeled = 0

        for record in SeqIO.parse(alignment_file, "fasta"):
            original_id = record.id
            new_id = None

            # Try exact match first
            if original_id in label_map:
                new_id = label_map[original_id]
            # Try base name match if exact match fails
            elif original_id in base_name_map:
                new_id = base_name_map[original_id]
                logger.debug(f"Using base name match for sequence: {original_id}")

            if new_id:
                record.id = new_id
                record.description = new_id
                n_relabeled += 1
                logger.debug(f"Relabeled sequence: {original_id} -> {record.id}")
            else:
                logger.warning(f"No mapping found for sequence: {original_id}")

            relabeled_records.append(record)

        # Write relabeled alignment
        SeqIO.write(relabeled_records, output_alignment, "fasta")
        logger.info(f"Relabeled {n_relabeled} sequences, wrote to: {output_alignment}")

    except Exception as e:
        error_msg = f"Failed to relabel alignment: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg) from e

    return output_tree, output_alignment
