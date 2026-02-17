# scripts/run_benchmark_with_timing.py
import time
import subprocess
import psutil
import pandas as pd
import sys
from pathlib import Path

def run_with_monitoring(input_file, output_dir):
    """
    Run BOLDGenotyper with detailed timing and resource monitoring.
    """
    results = {
        'input_file': str(input_file),
        'n_sequences': 0,
        'total_time': 0,
        'peak_memory_mb': 0,
        'avg_cpu_percent': 0,
        'success': False,
        'error': None
    }

    # Count sequences in input
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        results['n_sequences'] = len(df)
    except Exception as e:
        results['error'] = f"Failed to read input: {e}"
        return results

    # Start monitoring
    start_time = time.time()
    cpu_samples = []
    memory_samples = []

    try:
        # Launch process
        process = subprocess.Popen(
            ['boldgenotyper', str(input_file), '--output-dir', str(output_dir), '--build-tree'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Monitor while running
        ps_process = psutil.Process(process.pid)

        while process.poll() is None:
            try:
                cpu_samples.append(ps_process.cpu_percent(interval=1))
                memory_samples.append(ps_process.memory_info().rss / 1024 / 1024)  # MB
                time.sleep(1)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break

        # Wait for completion
        stdout, stderr = process.communicate()
        end_time = time.time()

        # Check success
        if process.returncode == 0:
            results['success'] = True
        else:
            results['error'] = stderr.decode('utf-8')

        # Calculate metrics
        results['total_time'] = end_time - start_time
        results['peak_memory_mb'] = max(memory_samples) if memory_samples else 0
        results['avg_cpu_percent'] = sum(cpu_samples) / len(cpu_samples) if cpu_samples else 0

    except Exception as e:
        results['error'] = str(e)

    return results

def run_all_benchmarks(benchmark_dir, output_csv):
    """
    Run all benchmark datasets and compile results.
    """
    results_list = []

    # Find all benchmark TSV files
    benchmark_files = sorted(Path(benchmark_dir).rglob("*.tsv"))

    print(f"Found {len(benchmark_files)} benchmark files")

    for i, input_file in enumerate(benchmark_files, 1):
        print(f"\nRunning benchmark {i}/{len(benchmark_files)}: {input_file.name}")

        # Create output directory
        output_dir = Path("benchmarking/runtime/output") / input_file.stem
        output_dir.mkdir(parents=True, exist_ok=True)

        # Run benchmark
        result = run_with_monitoring(input_file, output_dir)

        # Print summary
        if result['success']:
            print(f"  Success! Time: {result['total_time']:.1f}s, Memory: {result['peak_memory_mb']:.0f} MB")
        else:
            print(f"  FAILED: {result['error']}")

        results_list.append(result)

    # Save results
    results_df = pd.DataFrame(results_list)
    results_df.to_csv(output_csv, index=False)
    print(f"\nResults saved to: {output_csv}")

    return results_df

if __name__ == "__main__":
    results = run_all_benchmarks(
        benchmark_dir="benchmarking/runtime",
        output_csv="benchmarking/runtime/runtime_benchmark_results.csv"
    )

    print("\nBenchmark Summary:")
    print(results.groupby('n_sequences').agg({
        'total_time': ['mean', 'std'],
        'peak_memory_mb': ['mean', 'max']
    }))