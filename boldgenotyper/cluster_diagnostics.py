#!/usr/bin/env python

"""
Cluster diagnostics for dereplication + genotype assignment.

This script summarizes, per consensus cluster:

- n_reference: number of sequences used to build the consensus
               (parsed from consensus ID: consensus_c{cluster_id}_n{n})
- n_assigned: number of samples ultimately assigned to that consensus group
- n_cluster_members_fail: n_reference - n_assigned
  (approximate, assuming dereplication and assignment used the same dataset)
- Counts of tie / low_confidence / below_threshold / no_sequence statuses
- Identity distribution summaries (min/mean/median/max) per cluster
- Optional intra-cluster distance stats from the trimmed alignment
- A global dendrogram of all sequences to visually assess cluster tightness

ASSUMPTIONS
-----------
- Consensus IDs follow the pattern: "consensus_c{cluster_id}_n{n_reference}"
- The diagnostics file is the TSV written by genotype_assignment.assign_genotypes()
  via the diagnostics_path argument, with columns:
  ['processid', 'consensus_group', 'identity', 'target_identity',
   'classic_identity', 'identity_method', 'matches', 'mismatches',
   'insertions', 'deletions', 'edit_distance', 'length_discrepancy',
   'runner_up_group', 'runner_up_identity', 'is_tie', 'is_low_confidence',
   'status']
- If dereplication and assignment use exactly the same underlying dataset,
  then n_cluster_members_fail ≈ (n_reference - n_assigned). If derep used
  additional sequences (e.g. BOLD refs) that are not in your assignment run,
  this difference will include “members not present in assignment” as well.

DEPENDENCIES
------------
- Python 3.8+
- Biopython
- NumPy
- pandas
- SciPy
- matplotlib

USAGE
-----
Example:

    python cluster_diagnostics.py \
        --consensus consensus_Sphyrnidae.fasta \
        --diagnostics Sphyrnidae_assignment_diagnostics.tsv \
        --output-dir diagnostics_Sphyrnidae \
        --alignment Sphyrnidae_trimmed_alignment.fasta \
        --threshold 0.01
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

from Bio import SeqIO

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Import distance calculation from your dereplication module
from .dereplication import calculate_pairwise_distances, cluster_sequences

logger = logging.getLogger(__name__)


def parse_consensus_fasta(consensus_path: Path) -> pd.DataFrame:
    """
    Parse consensus FASTA and metadata to extract cluster information.

    Reads:
    - consensus FASTA: consensus_group and consensus_length
    - consensus metadata CSV: cluster_id and n_reference

    Expected ID pattern: "consensus_c{cluster_id}"
    e.g., "consensus_c24"
    """
    # Read consensus sequences to get consensus_group and length
    records = list(SeqIO.parse(str(consensus_path), "fasta"))
    fasta_rows = []

    for rec in records:
        header = rec.id
        token = header.split()[0]
        if not token.startswith("consensus_c"):
            logger.warning(f"Unexpected consensus ID format: {token}")
            continue

        fasta_rows.append({
            "consensus_group": token,
            "consensus_length": len(rec.seq),
        })

    fasta_df = pd.DataFrame(fasta_rows)
    if fasta_df.empty:
        raise ValueError(
            f"No valid consensus IDs parsed from {consensus_path}. "
            "Check naming convention."
        )

    # Read consensus metadata CSV to get cluster_id and n_reference
    metadata_path = consensus_path.parent / f"{consensus_path.stem}_metadata.csv"
    if not metadata_path.exists():
        raise ValueError(
            f"Consensus metadata CSV not found: {metadata_path}. "
            "This file should be generated during dereplication."
        )

    metadata_df = pd.read_csv(metadata_path)
    required_cols = ['consensus_group', 'cluster_id', 'n_reference']
    if not all(col in metadata_df.columns for col in required_cols):
        raise ValueError(
            f"Consensus metadata CSV missing required columns. "
            f"Expected: {required_cols}, Found: {list(metadata_df.columns)}"
        )

    # Merge FASTA data with metadata
    df = fasta_df.merge(metadata_df[required_cols], on='consensus_group', how='left')

    if df['cluster_id'].isna().any() or df['n_reference'].isna().any():
        raise ValueError(
            "Some consensus groups in FASTA not found in metadata CSV. "
            "Check that files match."
        )

    return df


def load_diagnostics(diagnostics_path: Path) -> pd.DataFrame:
    """
    Load the genotype assignment diagnostics file.

    Supports both TSV and CSV:
    - If the extension is .tsv or .txt, uses tab as a separator.
    - Otherwise, lets pandas infer the separator (CSV by default).
    """
    suffix = diagnostics_path.suffix.lower()
    if suffix in [".tsv", ".txt"]:
        df = pd.read_csv(diagnostics_path, sep="\t")
    else:
        # Assume CSV or let pandas infer
        df = pd.read_csv(diagnostics_path)

    # Normalize types a bit
    if "consensus_group" in df.columns:
        df["consensus_group"] = df["consensus_group"].fillna("").astype(str)
    if "status" in df.columns:
        df["status"] = df["status"].astype(str)
    return df
    

def compute_cluster_assignment_stats(
    consensus_df: pd.DataFrame,
    diag_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    For each consensus_group, compute:

    - n_reference
    - n_assigned (status == 'assigned')
    - n_tie
    - n_low_confidence
    - n_below_threshold
    - n_no_sequence
    - n_cluster_members_fail ≈ n_reference - n_assigned
    - identity summary stats (identity column)

    Note: consensus_df has base names (consensus_cX) from dereplication,
    while diag_df may have names with sample count suffix (consensus_cX_nZ)
    from genotype assignment. We extract the base name to match them.
    """
    # Ensure expected columns
    required_diag_cols = {
        "processid",
        "consensus_group",
        "identity",
        "status",
    }
    missing = required_diag_cols - set(diag_df.columns)
    if missing:
        raise ValueError(
            f"Diagnostics file is missing required columns: {missing}"
        )

    # Extract base consensus group name (strip _nZ suffix if present)
    # consensus_c7_n125 -> consensus_c7
    # consensus_c7 -> consensus_c7 (unchanged if no suffix)
    def extract_base_name(name):
        if pd.isna(name) or name == '':
            return name
        # Match pattern: consensus_cX_nZ -> extract consensus_cX
        # Use regex to remove _nZ suffix where Z is a number
        import re
        match = re.match(r'(consensus_c\d+)(?:_n\d+)?$', str(name))
        if match:
            return match.group(1)
        return name

    diag_df = diag_df.copy()
    diag_df['base_consensus_group'] = diag_df['consensus_group'].apply(extract_base_name)

    # Pre-aggregate diagnostics per base consensus group
    group = diag_df.groupby("base_consensus_group", dropna=False)

    stats_rows = []

    for _, row in consensus_df.iterrows():
        group_name = row["consensus_group"]  # This is the base name (consensus_cX)
        cluster_id = row["cluster_id"]
        n_reference = row["n_reference"]
        consensus_length = row["consensus_length"]

        if group_name in group.groups:
            sub = group.get_group(group_name)
            # Get the full consensus group name with _nZ suffix for the output
            # (use the most common one if there are multiple)
            full_group_names = sub['consensus_group'].value_counts()
            if len(full_group_names) > 0:
                full_group_name = full_group_names.index[0]
            else:
                full_group_name = group_name
        else:
            # No samples ever assigned or considered for this group
            sub = diag_df.iloc[0:0]  # empty
            full_group_name = group_name

        n_total_samples = len(sub)
        n_assigned = (sub["status"] == "assigned").sum()
        n_tie = (sub["status"] == "tie").sum()
        n_low_conf = (sub["status"] == "low_confidence").sum()
        n_below = (sub["status"] == "below_threshold").sum()
        n_no_seq = (sub["status"] == "no_sequence").sum()

        # Approximate: how many derep cluster members failed to get assigned
        n_cluster_members_fail = max(0, n_reference - n_assigned)

        # Identity stats for all samples whose best group is this consensus
        if n_total_samples > 0:
            identities = sub["identity"].values.astype(float)
            identity_min = float(np.min(identities))
            identity_mean = float(np.mean(identities))
            identity_median = float(np.median(identities))
            identity_max = float(np.max(identities))
        else:
            identity_min = identity_mean = identity_median = identity_max = np.nan

        stats_rows.append(
            {
                "consensus_group": full_group_name,  # Use full name with _nZ suffix
                "cluster_id": cluster_id,
                "n_reference": n_reference,
                "n_assigned": int(n_assigned),
                "n_cluster_members_fail": int(n_cluster_members_fail),
                "n_tie": int(n_tie),
                "n_low_confidence": int(n_low_conf),
                "n_below_threshold": int(n_below),
                "n_no_sequence": int(n_no_seq),
                "n_samples_considered": int(n_total_samples),
                "consensus_length": int(consensus_length),
                "identity_min": identity_min,
                "identity_mean": identity_mean,
                "identity_median": identity_median,
                "identity_max": identity_max,
            }
        )

    stats_df = pd.DataFrame(stats_rows).sort_values("cluster_id")
    return stats_df


def compute_intracluster_distance_stats(
    alignment_path: Path,
    threshold: float,
) -> Tuple[pd.DataFrame, np.ndarray, List[str], np.ndarray]:
    """
    Compute pairwise distances from the trimmed alignment and summarize
    intra-cluster distances.

    Returns:
    - cluster_distance_df: per-cluster distance summary stats
    - dist_condensed: condensed distance matrix (1D)
    - seq_ids: list of sequence IDs in same order as distance matrix
    - cluster_labels: array of cluster labels per sequence (1..n_clusters)
    """
    # Load alignment
    aln_records = list(SeqIO.parse(str(alignment_path), "fasta"))
    if len(aln_records) == 0:
        raise ValueError(f"No sequences found in alignment: {alignment_path}")

    seq_ids = [rec.id for rec in aln_records]

    # Compute pairwise distances (condensed)
    logger.info("Computing pairwise distances for alignment")
    dist_condensed = calculate_pairwise_distances(aln_records)

    # Cluster to recover cluster labels (threshold must match dereplication)
    logger.info(f"Clustering sequences at threshold={threshold}")
    cluster_labels = cluster_sequences(dist_condensed, threshold=threshold)

    # Convert to square matrix to easily extract intra-cluster distances
    dist_matrix = squareform(dist_condensed)

    # Summarize per cluster
    rows = []
    unique_clusters = np.unique(cluster_labels)
    for cid in unique_clusters:
        idxs = np.where(cluster_labels == cid)[0]
        if len(idxs) <= 1:
            # Singletons: no intra-cluster distances
            rows.append(
                {
                    "cluster_id": int(cid),
                    "cluster_size": int(len(idxs)),
                    "intra_min": np.nan,
                    "intra_mean": np.nan,
                    "intra_median": np.nan,
                    "intra_max": np.nan,
                    "intra_p95": np.nan,
                }
            )
            continue

        submat = dist_matrix[np.ix_(idxs, idxs)]
        # Take upper triangle without diagonal
        iu = np.triu_indices_from(submat, k=1)
        dists = submat[iu]

        rows.append(
            {
                "cluster_id": int(cid),
                "cluster_size": int(len(idxs)),
                "intra_min": float(np.min(dists)),
                "intra_mean": float(np.mean(dists)),
                "intra_median": float(np.median(dists)),
                "intra_max": float(np.max(dists)),
                "intra_p95": float(np.quantile(dists, 0.95)),
            }
        )

    cluster_distance_df = pd.DataFrame(rows).sort_values("cluster_id")
    return cluster_distance_df, dist_condensed, seq_ids, cluster_labels


def plot_global_dendrogram(
    dist_condensed: np.ndarray,
    seq_ids: List[str],
    output_dir: Path,
    title: str = "Global sequence dendrogram",
) -> None:
    """
    Plot a global dendrogram of all sequences using the condensed distance matrix.

    Saves:
    - cluster_dendrogram.pdf
    - cluster_dendrogram.svg
    """
    logger.info("Plotting global dendrogram")

    # Convert to linkage matrix
    Z = linkage(dist_condensed, method="average")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(6.0, len(seq_ids) * 0.01)))
    dendrogram(
        Z,
        labels=seq_ids,
        leaf_rotation=90,
        leaf_font_size=6,
        ax=ax,
    )
    ax.set_ylabel("Distance")
    ax.set_title(title)

    pdf_path = output_dir / "cluster_dendrogram.pdf"
    svg_path = output_dir / "cluster_dendrogram.svg"

    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)

    logger.info(f"Saved dendrogram to {pdf_path} and {svg_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Diagnostics for dereplication clusters and genotype assignment"
    )
    parser.add_argument(
        "--consensus",
        required=True,
        help="Path to consensus FASTA (e.g., consensus_Sphyrnidae.fasta)",
    )
    parser.add_argument(
        "--diagnostics",
        required=True,
        help="Path to genotype assignment diagnostics TSV",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for diagnostics",
    )
    parser.add_argument(
        "--alignment",
        required=False,
        help=(
            "Path to trimmed alignment FASTA used in dereplication. "
            "If provided, script will compute intra-cluster distance stats "
            "and plot a global dendrogram."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help=(
            "Distance threshold used for dereplication clustering "
            "(must match the dereplication threshold if alignment is provided). "
            "Default: 0.01"
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO)",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    consensus_path = Path(args.consensus)
    diagnostics_path = Path(args.diagnostics)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== Cluster diagnostics ===")
    logger.info(f"Consensus FASTA: {consensus_path}")
    logger.info(f"Diagnostics TSV: {diagnostics_path}")
    logger.info(f"Output directory: {output_dir}")

    # Load inputs
    consensus_df = parse_consensus_fasta(consensus_path)
    diag_df = load_diagnostics(diagnostics_path)

    # Compute assignment stats per consensus group
    cluster_stats_df = compute_cluster_assignment_stats(consensus_df, diag_df)

    # Optionally compute intra-cluster distances and dendrogram
    if args.alignment:
        alignment_path = Path(args.alignment)
        (
            intra_df,
            dist_condensed,
            seq_ids,
            cluster_labels,
        ) = compute_intracluster_distance_stats(
            alignment_path, threshold=args.threshold
        )

        # Merge intra-cluster stats into the main stats table
        cluster_stats_df = cluster_stats_df.merge(
            intra_df,
            on="cluster_id",
            how="left",
            suffixes=("", "_intra"),
        )

        # Plot global dendrogram
        plot_global_dendrogram(
            dist_condensed=dist_condensed,
            seq_ids=seq_ids,
            output_dir=output_dir,
            title="Global sequence dendrogram (trimmed alignment)",
        )

    # Write summary table
    out_tsv = output_dir / "cluster_diagnostics.tsv"
    cluster_stats_df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Wrote cluster diagnostics to {out_tsv}")

    logger.info("Diagnostics completed.")


if __name__ == "__main__":
    main()