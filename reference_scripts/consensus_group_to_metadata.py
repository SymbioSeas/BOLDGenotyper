#!/usr/bin/env python3
"""
Map each processid in metadata to a consensus sequence group.

Inputs
------
1) Metadata TSV with column 'processid'
2) Raw FASTA where headers contain the processid as a substring after the last
   underscore and before the first dot or end-of-line, e.g.:
     >Sphyrna_lewini_ANGBF11456-15.COI-5P
   -> processid = ANGBF11456-15
3) Consensus FASTA (e.g., ..._consensus.fasta) with headers like:
     >consensus_c7_n381

Output
------
Updated TSV with new column 'consensus_group' mapping each processid to the
best-matching consensus header.

Notes
-----
- Because consensus sequences are built from trimmed MSAs and majority-rule
  inference, they generally won't be exact copies of any single input read.
  We therefore assign by best sequence similarity (global edit distance).
- If edlib is installed, we use it for speed. Otherwise we use a pure-Python
  Levenshtein distance implementation.
- A diagnostics CSV can be emitted to inspect identities and ties.

Usage
-----
source ~/venvs/depredation/bin/activate
python consensus_group_to_metadata.py \
  --metadata Sphyrna_lewini_scallopedhammerhead.tsv \
  --fasta Sphyrna_lewini_scallopedhammerhead.fasta \
  --consensus Sphyrna_lewini_scallopedhammerhead_consensus.fasta \
  --out Sphyrna_lewini_scallopedhammerhead_with_consensus.tsv

Optional:
  --min-identity 0.90 \
  --threads 8 \
  --diag sphyrna_mapping_diagnostics.csv
"""

import argparse
import csv
import math
import multiprocessing as mp
import os
import re
import sys
from collections import namedtuple
Task = namedtuple("Task", ["pid", "seq"])
from functools import partial

try:
    import pandas as pd
except ImportError:
    pd = None
    # We'll fall back to csv.DictReader/writer if pandas not available.

# Optional fast aligner
try:
    import edlib
    HAVE_EDLIB = True
except Exception:
    HAVE_EDLIB = False


def read_fasta(path):
    """Return list of (header, seq) for a FASTA file (no Biopython required)."""
    records = []
    header = None
    seq_lines = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines).upper()))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_lines).upper()))
    return records


PID_REGEX = re.compile(r"_(?P<pid>[^.\s_]+)(?:[.\s]|$)")

def extract_pid_from_header(header):
    """
    Extract processid from raw FASTA header.
    Example:
      'Sphyrna_lewini_ANGBF11456-15.COI-5P' -> 'ANGBF11456-15'
    """
    m = PID_REGEX.search(header)
    if m:
        return m.group("pid")
    return None


def levenshtein_distance(a, b):
    """Pure-Python Levenshtein distance (edit distance)."""
    # Ensure a is the shorter for memory efficiency
    if len(a) > len(b):
        a, b = b, a
    previous = list(range(len(a) + 1))
    for j, ch_b in enumerate(b, start=1):
        current = [j]
        for i, ch_a in enumerate(a, start=1):
            ins = previous[i] + 1
            dele = current[i - 1] + 1
            sub = previous[i - 1] + (ch_a != ch_b)
            current.append(min(ins, dele, sub))
        previous = current
    return previous[-1]


def identity_from_distance(d, len_a, len_b):
    """Compute global identity-like score from edit distance."""
    denom = max(len_a, len_b)
    if denom == 0:
        return 1.0
    return 1.0 - (d / denom)


def best_match_for_sequence(seq, consensus_list, min_identity, use_edlib):
    """
    Given a sequence, return (best_group, best_identity, second_group, second_identity).
    consensus_list: list of (group_id, consensus_seq)
    """
    best_gid = None
    best_iden = -1.0
    second_gid = None
    second_iden = -1.0

    for gid, cseq in consensus_list:
        if use_edlib:
            # Global alignment distance
            res = edlib.align(seq, cseq, mode="NW", task="distance")
            d = res["editDistance"]
        else:
            d = levenshtein_distance(seq, cseq)
        iden = identity_from_distance(d, len(seq), len(cseq))

        if iden > best_iden:
            second_gid, second_iden = best_gid, best_iden
            best_gid, best_iden = gid, iden
        elif iden > second_iden:
            second_gid, second_iden = gid, iden

    if best_iden < min_identity:
        return (None, best_iden, second_gid, second_iden)
    return (best_gid, best_iden, second_gid, second_iden)


def load_metadata(path):
    if pd is not None:
        df = pd.read_csv(path, sep="\t", dtype=str)
        if "processid" not in df.columns:
            sys.exit("ERROR: metadata TSV is missing required column 'processid'.")
        return df
    # Fallback csv
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "processid" not in reader.fieldnames:
            sys.exit("ERROR: metadata TSV is missing required column 'processid'.")
        for row in reader:
            rows.append({k: (v if v is not None else "") for k, v in row.items()})
    return rows


def write_metadata(df_or_rows, out_path):
    if pd is not None and isinstance(df_or_rows, pd.DataFrame):
        df_or_rows.to_csv(out_path, sep="\t", index=False)
        return
    # Fallback csv writer for list-of-dicts
    rows = df_or_rows
    fieldnames = list(rows[0].keys()) if rows else ["processid", "consensus_group"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Append consensus group to metadata by aligning raw reads to consensus sequences."
    )
    ap.add_argument("--metadata", required=True, help="Input metadata .tsv with column 'processid'")
    ap.add_argument("--fasta", required=True, help="Raw FASTA with headers containing processid")
    ap.add_argument("--consensus", required=True, help="Consensus FASTA (headers like >consensus_c7_n381)")
    ap.add_argument("--out", required=True, help="Output metadata TSV with 'consensus_group' column added")
    ap.add_argument("--min-identity", type=float, default=0.90,
                    help="Minimum identity to accept an assignment (default 0.90)")
    ap.add_argument("--threads", type=int, default=1, help="Number of processes (parallelism)")
    ap.add_argument("--diag", default=None, help="Optional diagnostics CSV with identities and runner-up")
    return ap.parse_args()


def main():
    args = parse_args()

    use_edlib = HAVE_EDLIB
    if use_edlib:
        print("[info] Using edlib for fast global edit distance.")
    else:
        print("[info] edlib not found; using pure-Python Levenshtein distance (slower).")

    # Load inputs
    metadata = load_metadata(args.metadata)
    raw_records = read_fasta(args.fasta)
    cons_records = read_fasta(args.consensus)

    # Build consensus list (group_id -> seq)
    consensus_list = []
    for hdr, seq in cons_records:
        gid = hdr.split()[0]  # keep 'consensus_cX_nY'
        consensus_list.append((gid, seq))

    # Map processid -> sequence from raw FASTA
    pid_to_seq = {}
    missing_pid_headers = []
    for hdr, seq in raw_records:
        pid = extract_pid_from_header(hdr)
        if pid:
            # only record first occurrence if duplicated headers exist
            pid_to_seq.setdefault(pid, seq)
        else:
            missing_pid_headers.append(hdr)

    # Prepare list of processids to annotate (preserve order from metadata)
    if pd is not None and isinstance(metadata, pd.DataFrame):
        pids = metadata["processid"].astype(str).tolist()
    else:
        pids = [row["processid"] for row in metadata]

    # Work function for parallel mapping      
    tasks = [(pid,  pid_to_seq.get(pid)) for pid in pids]

    worker = partial(_worker_map, consensus_list=consensus_list,
                     min_identity=args.min_identity, use_edlib=use_edlib)

    if args.threads > 1:
        with mp.Pool(processes=args.threads) as pool:
            results = pool.map(worker, tasks)
    else:
        results = list(map(worker, tasks))

    # Collect outputs
    pid_to_group = {}
    diag_rows = []
    for r in results:
        pid_to_group[r["pid"]] = r["group"]
        if args.diag:
            diag_rows.append(r)

    # Append to metadata and write
    if pd is not None and isinstance(metadata, pd.DataFrame):
        metadata["consensus_group"] = metadata["processid"].map(pid_to_group)
        write_metadata(metadata, args.out)
    else:
        # list-of-dicts
        for row in metadata:
            row["consensus_group"] = pid_to_group.get(row["processid"])
        write_metadata(metadata, args.out)

    if args.diag:
        with open(args.diag, "w", newline="") as fh:
            fieldnames = ["pid", "group", "identity", "runner_up_group", "runner_up_identity"]
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for row in diag_rows:
                writer.writerow(row)

    # Summary
    n_total = len(pids)
    n_assigned = sum(1 for g in pid_to_group.values() if g)
    n_unassigned = n_total - n_assigned
    print(f"[done] Wrote: {args.out}")
    if args.diag:
        print(f"[info] Diagnostics: {args.diag}")
    print(f"[summary] total={n_total} assigned={n_assigned} unassigned={n_unassigned}")


def _worker_map(task, consensus_list, min_identity, use_edlib):
    pid, seq = task
    if not seq:
        return {"pid": pid, "group": None, "identity": 0.0,
                "runner_up_group": None, "runner_up_identity": 0.0}
    best_gid, best_iden, second_gid, second_iden = best_match_for_sequence(
        seq, consensus_list, min_identity, use_edlib
    )
    return {"pid": pid, "group": best_gid, "identity": round(best_iden, 6),
            "runner_up_group": second_gid, "runner_up_identity": round(second_iden, 6)}


if __name__ == "__main__":
    main()