#!/usr/bin/env python3
"""
Extract per-run pipeline metrics and per-phase timing from benchmark log files.

Produces: benchmarking/runtime/benchmark_details.csv

Columns:
  run_name, organism, n_sequences, n_qc_passed,
  n_haplotypes_raw, n_haplotypes_filtered, n_haplotypes_phylo, n_species,
  phase_1_sec, phase_2_sec, phase_3_sec, phase_4_sec, phase_5_sec,
  phase_6_sec, phase_7_sec, phase_8_sec, phase_9_sec, phase_9_5_sec,
  phase_10_sec, total_time_sec
"""

import re
import pandas as pd
from datetime import datetime
from pathlib import Path


def parse_log(log_path: Path) -> dict:
    text = log_path.read_text()
    row = {'run_name': log_path.stem}

    # --- Scalar metrics (single values extracted by regex) ---

    m = re.search(r'Read (\d+) rows', text)
    row['n_sequences'] = int(m.group(1)) if m else None

    m = re.search(r'Retained (\d+) samples after QC', text)
    row['n_qc_passed'] = int(m.group(1)) if m else None

    # "Identified 359 unique haplotypes from 782 sequences"
    m = re.search(r'Identified (\d+) unique haplotypes from (\d+) sequences', text)
    row['n_haplotypes_raw'] = int(m.group(1)) if m else None

    # "Retained haplotypes: 187 (52.1%)"
    m = re.search(r'Retained haplotypes:\s*(\d+)', text)
    row['n_haplotypes_filtered'] = int(m.group(1)) if m else None

    # "Consensus filtering: 51/173 sequences retained for phylogenetics"
    m = re.search(r'Consensus filtering:\s*(\d+)/(\d+)\s*sequences retained', text)
    row['n_haplotypes_phylo'] = int(m.group(1)) if m else None

    # "Number of species: 106"
    m = re.search(r'Number of species:\s*(\d+)', text)
    row['n_species'] = int(m.group(1)) if m else None

    # --- Per-phase timing ---
    # Each phase line looks like:
    #   [2026-01-31 09:59:36] INFO: PHASE 1: Data Loading ...
    # We collect (timestamp, phase_label) pairs and compute durations.

    phase_pattern = re.compile(
        r'\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]\s+INFO:\s+PHASE\s+([\d.]+):'
    )
    phases = []
    for m in phase_pattern.finditer(text):
        ts = datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S')
        label = m.group(2)
        phases.append((ts, label))

    # Also grab the pipeline completion timestamp for the final boundary
    m = re.search(
        r'\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]\s+INFO:.*Pipeline completed',
        text
    )
    end_ts = datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S') if m else None

    # Deduplicate phases: keep first occurrence of each label
    # (PHASE 6.5 appears twice in the log — once for species aggregation,
    # once for plot export; we want only the first)
    seen = set()
    unique_phases = []
    for ts, label in phases:
        if label not in seen:
            seen.add(label)
            unique_phases.append((ts, label))
    phases = unique_phases

    # Map phase labels to canonical output column names
    label_to_col = {
        '1': 'phase_1_sec',       # Data Loading & QC
        '2': 'phase_2_sec',       # Pre-processing QC
        '3': 'phase_3_sec',       # Haplotype Discovery
        '4': 'phase_4_sec',       # Haplotype Assignment
        '5': 'phase_5_sec',       # Taxonomy
        '6': 'phase_6_sec',       # Post-assignment QC
        '6.5': 'phase_6_5_sec',   # Species-Level Aggregation
        '7': 'phase_7_sec',       # Geographic Enhancement
        '8': 'phase_8_sec',       # Phylogenetic Analysis
        '9': 'phase_9_sec',       # Divergence Analysis
        '9.5': 'phase_9_5_sec',   # Metadata Analysis
        '10': 'phase_10_sec',     # Visualization
    }

    # Compute duration of each phase as (next_phase_start - this_phase_start)
    # Last phase duration = (pipeline_end - last_phase_start)
    for i, (ts, label) in enumerate(phases):
        col = label_to_col.get(label)
        if col is None:
            continue  # skip unlabelled phases (e.g. the duplicate 6.5)

        if i + 1 < len(phases):
            next_ts = phases[i + 1][0]
        elif end_ts:
            next_ts = end_ts
        else:
            next_ts = None

        if next_ts:
            row[col] = (next_ts - ts).total_seconds()
        else:
            row[col] = None

    # Total time from first phase to pipeline completion
    if phases and end_ts:
        row['total_time_sec'] = (end_ts - phases[0][0]).total_seconds()

    return row


def main():
    results_dir = Path('benchmarking/runtime/results')
    log_files = sorted(results_dir.glob('*.log'))

    rows = []
    for log_path in log_files:
        row = parse_log(log_path)
        # Extract organism from run_name
        m = re.match(r'^(.+?)_n\d+_rep\d+$', row['run_name'])
        row['organism'] = m.group(1) if m else None
        rows.append(row)

    df = pd.DataFrame(rows)

    # Drop the one malformed entry (no organism match)
    df = df[df['organism'].notna()].copy()

    # Only keep runs that have the metrics we need
    required = ['n_sequences', 'n_haplotypes_filtered', 'total_time_sec']
    before = len(df)
    df = df.dropna(subset=required)
    dropped = before - len(df)
    if dropped:
        print(f"Dropped {dropped} runs missing required fields")

    output = Path('benchmarking/runtime/benchmark_details.csv')
    df.to_csv(output, index=False)
    print(f"Wrote {len(df)} runs to {output}")
    print(f"\nColumns: {list(df.columns)}")
    print(f"\nSample (Carcharhiniformes_n1000_rep1):")
    sample = df[df['run_name'] == 'Carcharhiniformes_n1000_rep1']
    if len(sample):
        print(sample.T.to_string())


if __name__ == '__main__':
    main()
