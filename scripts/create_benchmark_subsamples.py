# scripts/create_benchmark_subsamples.py
import pandas as pd
import random
import os

def subsample_dataset(input_file, output_dir, sizes, n_replicates=3, seed=42):
    """
    Create subsampled datasets for benchmarking.

    Parameters:
    -----------
    input_file : str
        Path to full BOLD TSV file
    output_dir : str
        Directory to save subsampled files
    sizes : list of int
        Sample sizes to create
    n_replicates : int
        Number of replicates per size
    seed : int
        Random seed for reproducibility
    """
    random.seed(seed)

    # Read full dataset with lenient parsing for files with variable column counts
    try:
        # Try normal read first
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except pd.errors.ParserError as e:
        # If that fails, use lenient parsing that skips problematic rows
        # This handles files with inconsistent column counts (extra or embedded delimiters)
        print(f"  Note: File has parsing issues ({str(e)[:50]}...), using lenient parsing")
        df = pd.read_csv(input_file, sep='\t', on_bad_lines='skip')
    total_samples = len(df)

    print(f"Input: {input_file}")
    print(f"Total samples: {total_samples}")

    os.makedirs(output_dir, exist_ok=True)

    for size in sizes:
        if size > total_samples:
            print(f"Warning: Requested size {size} exceeds total samples {total_samples}, skipping")
            continue

        for rep in range(1, n_replicates + 1):
            # Random sample
            sample_df = df.sample(n=size, random_state=seed + rep)

            # Output filename
            basename = os.path.splitext(os.path.basename(input_file))[0]
            output_file = f"{output_dir}/{basename}_n{size}_rep{rep}.tsv"

            # Save
            sample_df.to_csv(output_file, sep='\t', index=False)
            print(f"Created: {output_file}")

# Usage example
if __name__ == "__main__":
    # Sphyrnidae (small)
    subsample_dataset(
        input_file="data/Sphyrnidae.tsv",
        output_dir="benchmarking/runtime/Sphyrnidae",
        sizes=[100, 250, 500, 750, 1000]
    )
    
    # Panulirus (medium)
    subsample_dataset(
        input_file="data/Panulirus.tsv",
        output_dir="benchmarking/runtime/Panulirus",
        sizes=[100, 500, 750, 1000, 2000]
    )
    
    # Salmonidae (large)
    subsample_dataset(
        input_file="data/Salmonidae.tsv",
        output_dir="benchmarking/runtime/Salmonidae",
        sizes=[100, 500, 1000, 2000, 5000]
    )
    
    # Carcharhiniformes (very large)
    subsample_dataset(
        input_file="data/Carcharhiniformes.tsv",
        output_dir="benchmarking/runtime/Carcharhiniformes",
        sizes=[100, 500, 1000, 2500, 5000, 7500, 10000]
    )
    
    # Pieridae (largest)
    subsample_dataset(
        input_file="data/Pieridae.tsv",
        output_dir="benchmarking/runtime/Pieridae",
        sizes=[100, 500, 1000, 2500, 5000, 10000, 20000]
    )

