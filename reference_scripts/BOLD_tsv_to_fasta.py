#!/usr/bin/env python3
"""
tsv_to_fasta.py

Given a tab-delimited file with columns "record_id", "species", and "nuc",
writes a FASTA file named `<basename>.fasta` in the same directory as the input.
Headers are formatted as: >species_record_id
"""

import argparse
import csv
import os
import sys
import textwrap

def parse_args():
    p = argparse.ArgumentParser(description="Convert TSV of sequences to FASTA")
    p.add_argument("tsv", help="Input TSV file (must have 'record_id', 'species', 'nuc' columns)")
    p.add_argument(
        "--wrap", type=int, default=80,
        help="Wrap nucleotide sequences to this width (default: 80; set 0 to disable wrapping)"
    )
    return p.parse_args()

def main():
    args = parse_args()

    # Verify input file exists
    if not os.path.isfile(args.tsv):
        sys.exit(f"Error: file not found: {args.tsv}")

    # Build output path
    base, _ = os.path.splitext(os.path.basename(args.tsv))
    out_dir = os.path.dirname(os.path.abspath(args.tsv))
    fasta_path = os.path.join(out_dir, f"{base}.fasta")

    # Read TSV and write FASTA
    with open(args.tsv, newline='') as infile, open(fasta_path, 'w') as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        # Ensure required columns are present
        for col in ("record_id", "species", "nuc"):
            if col not in reader.fieldnames:
                sys.exit(f"Error: column '{col}' not found in header")

        for row in reader:
            rid = row["record_id"].strip()
            species = row["species"].strip().replace(" ", "_")
            seq = row["nuc"].strip()

            if not seq:
                # Skip empty sequences (optional)
                continue

            header = f">{species}_{rid}"
            outfile.write(header + "\n")

            if args.wrap > 0:
                # wrap sequence for readability
                wrapped = textwrap.fill(seq, width=args.wrap)
                outfile.write(wrapped + "\n")
            else:
                # no wrapping
                outfile.write(seq + "\n")

    print(f"Wrote {fasta_path}")

if __name__ == "__main__":
    main()