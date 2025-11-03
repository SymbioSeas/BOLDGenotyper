#!/usr/bin/env python3
"""
msa_to_consensus.py

For each FASTA in a directory:
 1) Align with MAFFT → alignments/<base>.aln.fasta
 2) Trim the MSA with trimAl (-automated1) → trimmed/<base>.trimmed.fasta
 3) Compute distance matrix (1 - pct_identity), ignoring gaps ('-') and Ns
 4) Cluster (hierarchical; default cutoff = 0.01)
 5) Build majority-rule consensus per cluster, ignoring gaps and Ns
 6) Write consensus/<base>_consensus.fasta

Requires:
  - mafft in your PATH
  - trimAl in your PATH
  - Python3 + biopython + scipy:
      pip install biopython scipy
"""

import os
import sys
import argparse
import subprocess
import shutil
import textwrap
from collections import Counter

from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, fcluster


def compute_distance(seq1: str, seq2: str) -> float:
    """1 - identity between two aligned sequences, ignoring gaps/Ns."""
    matches = sites = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a in "ACGT" and b in "ACGT":
            sites += 1
            if a == b:
                matches += 1
    return 1.0 if sites == 0 else 1.0 - (matches / sites)


def process_file(fasta_path, dirs, threshold, freq_cutoff, wrap_width):
    base = os.path.splitext(os.path.basename(fasta_path))[0]
    aln_path     = os.path.join(dirs["align"],    f"{base}.aln.fasta")
    trimmed_path = os.path.join(dirs["trimmed"],  f"{base}.trimmed.fasta")
    out_path     = os.path.join(dirs["consensus"],f"{base}_consensus.fasta")

    print(f"\n==> {base}")

    # 1) MAFFT
    try:
        with open(aln_path, "w") as h:
            subprocess.run(
                ["mafft", "--auto", fasta_path],
                stdout=h, stderr=subprocess.PIPE,
                text=True, check=True
            )
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"MAFFT failed on {base}:\n{e.stderr}\n")
        return

    # 2) trimAl
    try:
        subprocess.run(
            ["trimal", "-in", aln_path, "-out", trimmed_path, "-automated1"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"trimAl failed on {base}:\n{e.stderr}\n")
        return

    # 3) Read trimmed MSA
    alignment = AlignIO.read(trimmed_path, "fasta")
    seqs = [str(rec.seq) for rec in alignment]
    n = len(seqs)

    # 4) Distance vector
    dist_vec = []
    for i in range(n):
        for j in range(i+1, n):
            dist_vec.append(compute_distance(seqs[i], seqs[j]))

    # 5) Clustering
    Z = linkage(dist_vec, method="average")
    labels = fcluster(Z, t=threshold, criterion="distance")

    # 6) Group into clusters
    clusters = {}
    for rec, cid in zip(alignment, labels):
        clusters.setdefault(cid, []).append(rec)

    # 7) Consensus
    L = alignment.get_alignment_length()
    with open(out_path, "w") as out_h:
        for cid, members in clusters.items():
            cons = []
            for i in range(L):
                col = [r.seq[i] for r in members]
                legit = [b.upper() for b in col if b.upper() in "ACGT"]
                if not legit:
                    cons.append("N")
                else:
                    base, cnt = Counter(legit).most_common(1)[0]
                    cons.append(base if cnt/len(legit) >= freq_cutoff else "N")
            seq = "".join(cons)

            out_h.write(f">consensus_c{cid}_n{len(members)}\n")
            if wrap_width > 0:
                out_h.write(textwrap.fill(seq, width=wrap_width) + "\n")
            else:
                out_h.write(seq + "\n")

    print(f"   → consensus: {os.path.basename(out_path)}")


def main():
    p = argparse.ArgumentParser(
        description="MSA→trim→cluster→consensus over all FASTAs in a directory"
    )
    p.add_argument("indir",
                   help="Directory containing FASTA files (*.fa, *.fasta, *.fna)")
    p.add_argument("--threshold", type=float, default=0.01,
                   help="Distance cutoff (default 0.01 = 99%% id)")
    p.add_argument("--freq-cutoff", type=float, default=0.7,
                   help="Min fraction to call consensus base (default 0.7)")
    p.add_argument("--wrap", type=int, default=80,
                   help="Wrap width for consensus sequences (default 80; 0 to disable)")
    args = p.parse_args()

    # check tools & input
    if not os.path.isdir(args.indir):
        sys.exit(f"ERROR: not a directory: {args.indir}")
    for tool in ("mafft", "trimal"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: '{tool}' not found in PATH")

    # make subdirs
    dirs = {d: os.path.join(args.indir, d) for d in ("align","trimmed","consensus")}
    for path in dirs.values():
        os.makedirs(path, exist_ok=True)

    # gather FASTAs
    exts = (".fa", ".fasta", ".fna")
    fastas = sorted(f for f in os.listdir(args.indir)
                    if f.lower().endswith(exts))
    if not fastas:
        sys.exit("No FASTA files found.")

    # run
    for fn in fastas:
        process_file(os.path.join(args.indir, fn),
                     dirs, args.threshold, args.freq_cutoff, args.wrap)


if __name__ == "__main__":
    main()