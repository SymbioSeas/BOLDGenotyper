"""
Tree Building and Phylogenetic Analysis

This module provides optional phylogenetic analysis functionality for consensus
sequences, including multiple sequence alignment, maximum likelihood tree
construction, and bootstrap support calculation.

Phylogenetic Analysis Modes:

Mode 1: No phylogeny (default)
- Genotype identification only
- Fastest analysis option
- Suitable for exploratory analysis

Mode 2: --phylogeny flag
- Build tree from consensus sequences
- Midpoint rooted
- No outgroup required
- Good for visualizing relationships among genotypes

Mode 3: --phylogeny --outgroup <fasta>
- User provides outgroup sequences
- Proper phylogenetic rooting
- Requires taxonomic expertise for outgroup selection
- Best for publication-quality phylogenies

Why Phylogeny is Optional:
- Outgroup selection requires taxonomic knowledge
- Not all users need phylogenetic trees
- Makes pipeline more flexible for exploratory analysis
- Reduces computational time for basic genotyping

Workflow:
1. Align consensus sequences using MAFFT (--auto mode)
2. Optional: Trim alignment with trimAl
3. Build maximum likelihood tree with PhyML
   - Model: GTR+G+I (General Time Reversible with gamma and invariant sites)
   - Bootstrap: 1000 replicates (or user-specified)
4. Midpoint root (if no outgroup provided)
5. Export tree in Newick and visualization formats

Dependencies:
- MAFFT v7+ for alignment
- PhyML 3.0+ for tree building
- Biopython for tree parsing
- ete3 or DendroPy for tree manipulation

Tree Output:
- Newick format: {organism}_phylogeny.tree
- PDF visualization: {organism}_phylogeny.pdf
- Bootstrap support values included as node labels

Example Usage:
    >>> from boldgenotyper.phylogenetics import build_phylogeny
    >>> tree = build_phylogeny(
    ...     consensus_fasta="Sphyrna_lewini_consensus.fasta",
    ...     output_prefix="results/Sphyrna_lewini",
    ...     outgroup_fasta="outgroup_sequences.fasta",
    ...     bootstrap=1000
    ... )

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


def build_phylogeny(
    consensus_fasta: str,
    output_prefix: str,
    outgroup_fasta: Optional[str] = None,
    model: str = "GTR",
    bootstrap: int = 1000,
    threads: int = 1,
) -> Optional[Phylo.BaseTree.Tree]:
    """
    Build maximum likelihood phylogenetic tree from consensus sequences.

    Parameters
    ----------
    consensus_fasta : str
        Path to consensus sequences FASTA file
    output_prefix : str
        Prefix for output files (tree, alignment, etc.)
    outgroup_fasta : str, optional
        Path to outgroup sequences for rooting
    model : str, optional
        Substitution model for PhyML (default: GTR)
    bootstrap : int, optional
        Number of bootstrap replicates (default: 1000)
    threads : int, optional
        Number of CPU threads (default: 1)

    Returns
    -------
    Phylo.BaseTree.Tree or None
        Phylogenetic tree object, or None if tree building failed
    """
    try:
        output_prefix_path = Path(output_prefix)
        output_prefix_path.parent.mkdir(parents=True, exist_ok=True)

        # Combine with outgroup if provided
        input_fasta = consensus_fasta
        if outgroup_fasta:
            combined_fasta = f"{output_prefix}_combined.fasta"
            input_fasta = add_outgroup(consensus_fasta, outgroup_fasta, combined_fasta)

        # Run MAFFT alignment
        aligned_fasta = f"{output_prefix}_aligned.fasta"
        run_mafft_alignment(input_fasta, aligned_fasta, threads=threads)

        # Run FastTree to build phylogeny
        tree_file = f"{output_prefix}_tree.nwk"
        run_fasttree(aligned_fasta, tree_file, model=model)

        # Load and return tree
        if Path(tree_file).exists():
            tree = Phylo.read(tree_file, "newick")
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
