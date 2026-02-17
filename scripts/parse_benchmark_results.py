#!/usr/bin/env python3
"""
Parse runtime benchmark outputs from run_runtime_benchmarks.sh into
the CSV format expected by analyze_runtime_benchmarks.py.

Reads two file types per benchmark run:
  - {name}_time.txt : /usr/bin/time -l output (real/user/sys times, memory)
  - {name}.log      : pipeline stdout (sequence count, success/failure)

Produces:
  benchmarking/runtime/runtime_benchmark_results.csv
"""

import re
import pandas as pd
from pathlib import Path


def parse_time_file(time_path: Path) -> dict:
    """
    Parse a /usr/bin/time -l output file.

    The file may contain Python warnings on stderr before the actual
    time output.  We locate fields by their label, not line number.
    """
    text = time_path.read_text()

    result = {
        'real_seconds': None,
        'user_seconds': None,
        'sys_seconds': None,
        'peak_memory_bytes': None,
    }

    # "real" timing line:  "      291.92 real       284.59 user         5.50 sys"
    timing_match = re.search(
        r'([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys', text
    )
    if timing_match:
        result['real_seconds'] = float(timing_match.group(1))
        result['user_seconds'] = float(timing_match.group(2))
        result['sys_seconds'] = float(timing_match.group(3))

    # "peak memory footprint" (macOS 10.13+, in bytes)
    mem_match = re.search(r'([\d]+)\s+peak memory footprint', text)
    if mem_match:
        result['peak_memory_bytes'] = int(mem_match.group(1))
    else:
        # Fall back to "maximum resident set size"
        mem_match = re.search(r'([\d]+)\s+maximum resident set size', text)
        if mem_match:
            result['peak_memory_bytes'] = int(mem_match.group(1))

    return result


def parse_log_file(log_path: Path) -> dict:
    """
    Parse a pipeline .log file for sequence count and success status.
    """
    text = log_path.read_text()

    result = {
        'n_sequences_read': None,
        'success': False,
        'error': None,
    }

    # Sequence count: "Read 989 rows and 91 columns"
    seq_match = re.search(r'Read (\d+) rows and \d+ columns', text)
    if seq_match:
        result['n_sequences_read'] = int(seq_match.group(1))

    # Success: "Pipeline completed successfully"
    if 'Pipeline completed successfully' in text:
        result['success'] = True
    else:
        # Capture the first ERROR line as the failure reason
        error_match = re.search(r'ERROR: (.+)', text)
        if error_match:
            result['error'] = error_match.group(1).strip()
        else:
            result['error'] = 'Unknown failure (no ERROR line found)'

    return result


def extract_n_from_filename(name: str) -> int:
    """
    Pull the intended sample size from a benchmark filename.
    e.g. "Carcharhiniformes_n1000_rep2" → 1000
    """
    match = re.search(r'_n(\d+)_', name)
    if match:
        return int(match.group(1))
    # Fallback: try end of string
    match = re.search(r'_n(\d+)$', name)
    return int(match.group(1)) if match else None


def main():
    results_dir = Path('benchmarking/runtime/results')
    if not results_dir.exists():
        raise FileNotFoundError(
            f"Expected results directory not found: {results_dir}\n"
            "Run scripts/run_runtime_benchmarks.sh first."
        )

    # Collect all _time.txt files (each has a matching .log)
    time_files = sorted(results_dir.glob('*_time.txt'))
    if not time_files:
        raise FileNotFoundError(
            f"No *_time.txt files found in {results_dir}\n"
            "Run scripts/run_runtime_benchmarks.sh first."
        )

    rows = []
    for time_path in time_files:
        # Derive companion log path: remove "_time" suffix
        # e.g. Carcharhiniformes_n1000_rep1_time.txt → Carcharhiniformes_n1000_rep1.log
        stem = time_path.stem  # "Carcharhiniformes_n1000_rep1_time"
        base_name = stem.replace('_time', '')  # "Carcharhiniformes_n1000_rep1"
        log_path = time_path.parent / f"{base_name}.log"

        if not log_path.exists():
            print(f"WARNING: No log file for {time_path.name}, skipping")
            continue

        # Parse both files
        time_data = parse_time_file(time_path)
        log_data = parse_log_file(log_path)

        # If /usr/bin/time produced no output at all (process crashed before
        # time could record), the fields will be None — mark as failed.
        if time_data['real_seconds'] is None:
            log_data['success'] = False
            if log_data['error'] is None:
                log_data['error'] = 'No timing data captured (process may have crashed)'

        # Compute derived metrics
        total_time = time_data['real_seconds']
        peak_memory_mb = (
            time_data['peak_memory_bytes'] / (1024 * 1024)
            if time_data['peak_memory_bytes'] else None
        )

        # CPU % ≈ (user + sys) / real × 100
        if all(v is not None for v in [
            time_data['user_seconds'],
            time_data['sys_seconds'],
            time_data['real_seconds'],
        ]) and time_data['real_seconds'] > 0:
            avg_cpu_percent = (
                (time_data['user_seconds'] + time_data['sys_seconds'])
                / time_data['real_seconds'] * 100
            )
        else:
            avg_cpu_percent = None

        rows.append({
            'run_name': base_name,
            'input_file': f"benchmarking/runtime/{base_name}.tsv",
            'n_sequences': extract_n_from_filename(base_name),
            'n_sequences_read': log_data['n_sequences_read'],
            'total_time': total_time,
            'peak_memory_mb': round(peak_memory_mb, 1) if peak_memory_mb else None,
            'avg_cpu_percent': round(avg_cpu_percent, 1) if avg_cpu_percent else None,
            'success': log_data['success'],
            'error': log_data['error'],
        })

    df = pd.DataFrame(rows)

    # Print summary
    n_total = len(df)
    n_success = df['success'].sum()
    n_failed = n_total - n_success
    print(f"\nParsed {n_total} benchmark runs")
    print(f"  Successful: {n_success}")
    print(f"  Failed:     {n_failed}")

    if n_failed > 0:
        print("\nFailed runs:")
        for _, row in df[~df['success']].iterrows():
            print(f"  {row['run_name']}: {row['error']}")

    # Print summary table for successful runs
    if n_success > 0:
        success_summary = (
            df[df['success']]
            .groupby('n_sequences')
            .agg(
                runs=('total_time', 'count'),
                mean_time_sec=('total_time', 'mean'),
                mean_memory_mb=('peak_memory_mb', 'mean'),
                mean_cpu_pct=('avg_cpu_percent', 'mean'),
            )
            .round(1)
        )
        print("\nSummary by dataset size:")
        print(success_summary.to_string())

    # Write CSV
    output_path = Path('benchmarking/runtime/runtime_benchmark_results.csv')
    df.to_csv(output_path, index=False)
    print(f"\nWrote: {output_path}")


if __name__ == '__main__':
    main()
