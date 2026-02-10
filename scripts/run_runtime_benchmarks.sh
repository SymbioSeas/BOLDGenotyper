#!/bin/bash
# scripts/run_runtime_benchmarks.sh

# Create output directory
mkdir -p benchmarking/runtime/results

# Function to run benchmarked analysis
run_benchmark() {
    input_file=$1
    output_prefix=$2

    # Extract size and replicate from filename
    size=$(echo $input_file | sed -E 's/.*_n([0-9]+).*/\1/')
    rep=$(echo $input_file | sed -E 's/.*_rep([0-9]+).*/\1/')

    # Output files
    log_file="benchmarking/runtime/results/${output_prefix}_n${size}_rep${rep}.log"
    time_file="benchmarking/runtime/results/${output_prefix}_n${size}_rep${rep}_time.txt"
    mem_file="benchmarking/runtime/results/${output_prefix}_n${size}_rep${rep}_mem.txt"

    echo "Running: $input_file (n=$size, rep=$rep)"

    # Run with time and memory tracking
    /usr/bin/time -l boldgenotyper $input_file \
        --output-dir benchmarking/runtime/output/${output_prefix}_n${size}_rep${rep} \
        --build-tree \
        2> $time_file \
        > $log_file

    # Extract timing from log
    grep "Pipeline completed" $log_file || echo "ERROR: Pipeline failed" >> $log_file

    echo "Completed: n=$size, rep=$rep"
}

export -f run_benchmark

# Run all Panulirus benchmarks
find benchmarking/runtime/Panulirus -name "*.tsv" | while read file; do
    run_benchmark "$file" "Panulirus"
done

# Run all Sphyrnidae benchmarks
find benchmarking/runtime/Sphyrnidae -name "*.tsv" | while read file; do
    run_benchmark "$file" "Sphyrnidae"
done

# Run all Carcharhinidae benchmarks (large datasets - may take hours)
find benchmarking/runtime/Carcharhiniformes -name "*.tsv" | while read file; do
    run_benchmark "$file" "Carcharhiniformes"
done

# Run all Salmonidae benchmarks (very large - may take many hours)
find benchmarking/runtime/Salmonidae -name "*.tsv" | while read file; do
    run_benchmark "$file" "Salmonidae"
done

# Run all Pieridae benchmarks (largest - may take many hours)
find benchmarking/runtime/Pieridae -name "*.tsv" | while read file; do
    run_benchmark "$file" "Pieridae"
done

echo "All runtime benchmarks complete!"