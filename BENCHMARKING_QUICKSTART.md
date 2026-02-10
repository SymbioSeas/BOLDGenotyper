# BOLDGenotyper Benchmarking Quick Start Guide

**Purpose**: Get started with benchmarking immediately
**Time to read**: 5 minutes
**Date**: 2026-01-28

---

## Step 1: Download Datasets from BOLD (Today)

Navigate to https://boldsystems.org/ and download these 5 datasets:

| # | Search Term | Download As | Expected Size | Save As |
|---|-------------|-------------|---------------|---------|
| 1 | `Sphyrnidae` (family, COI-5P) | Combined TSV | ~1,400 sequences | `data/Sphyrnidae.tsv` |
| 2 | `Carcharhinidae` (family, COI-5P) | Combined TSV | ~8,000-10,000 sequences | `data/Carcharhinidae.tsv` |
| 3 | `Panulirus` (genus, COI-5P) | Combined TSV | ~1,500 sequences | `data/Panulirus.tsv` |
| 4 | `Salmonidae` (family, COI-5P) | Combined TSV | ~10,000+ sequences | `data/Salmonidae.tsv` |
| 5 | `Pieridae` (family, COI-5P) | Combined TSV | ~3,000-5,000 sequences | `data/Pieridae.tsv` |

**After download**, verify file sizes:
```bash
cd /Users/stesmith/Documents/depredation/boldgenotyper/data
ls -lh *.tsv
wc -l *.tsv  # Check line counts
```

---

## Step 2: Quick Test Runs (This Week)

Before full benchmarking, run quick tests to verify pipeline compatibility:

```bash
# Test each dataset with default parameters (no tree, fast)
boldgenotyper data/Sphyrnidae.tsv --no-geo
boldgenotyper data/Carcharhinidae.tsv --no-geo
boldgenotyper data/Panulirus.tsv --no-geo
boldgenotyper data/Salmonidae.tsv --no-geo
boldgenotyper data/Pieridae.tsv --no-geo
```

**Verify**:
- [ ] All datasets parse correctly
- [ ] No immediate errors
- [ ] Output files generated
- [ ] Reasonable runtime (<30 min for most)

---

## Step 3: Create Benchmarking Directory Structure

```bash
cd /Users/stesmith/Documents/depredation/boldgenotyper

# Create directories
mkdir -p benchmarking/runtime/{Panulirus,Sphyrnidae,Carcharhinidae,Salmonidae}/
mkdir -p benchmarking/runtime/{output,results}/
mkdir -p benchmarking/accuracy/
mkdir -p benchmarking/parameter_sweep/{Sphyrnidae,Panulirus,Salmonidae,Pieridae}/
mkdir -p benchmarking/comparative_analysis/
mkdir -p scripts/
```

---

## Step 4: Create Subsampling Script

Save this as `scripts/create_benchmark_subsamples.py`:

```python
import pandas as pd
import random
import os
from pathlib import Path

def subsample_dataset(input_file, output_dir, sizes, n_replicates=3, seed=42):
    """Create subsampled datasets for benchmarking."""
    random.seed(seed)

    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    total_samples = len(df)

    print(f"Input: {input_file}")
    print(f"Total samples: {total_samples}")

    os.makedirs(output_dir, exist_ok=True)

    for size in sizes:
        if size > total_samples:
            print(f"Warning: Size {size} exceeds total {total_samples}, skipping")
            continue

        for rep in range(1, n_replicates + 1):
            sample_df = df.sample(n=size, random_state=seed + rep)

            basename = Path(input_file).stem
            output_file = f"{output_dir}/{basename}_n{size}_rep{rep}.tsv"

            sample_df.to_csv(output_file, sep='\t', index=False)
            print(f"  Created: {output_file}")

if __name__ == "__main__":
    # Panulirus (small)
    subsample_dataset(
        "data/Panulirus.tsv",
        "benchmarking/runtime/Panulirus",
        [100, 250, 500, 750, 1000]
    )

    # Sphyrnidae (medium)
    subsample_dataset(
        "data/Sphyrnidae.tsv",
        "benchmarking/runtime/Sphyrnidae",
        [100, 500, 1000, 1500]
    )

    # Carcharhinidae (large)
    subsample_dataset(
        "data/Carcharhinidae.tsv",
        "benchmarking/runtime/Carcharhinidae",
        [100, 500, 1000, 2500, 5000, 7500, 10000]
    )

    # Salmonidae (very large)
    subsample_dataset(
        "data/Salmonidae.tsv",
        "benchmarking/runtime/Salmonidae",
        [100, 500, 1000, 5000, 10000]
    )

    print("\nAll subsamples created!")
```

Run it:
```bash
python scripts/create_benchmark_subsamples.py
```

---

## Step 5: Run Runtime Benchmarks (Allow 1-3 Days)

This will take substantial time. Consider running overnight or over a weekend.

Save as `scripts/run_all_benchmarks.sh`:

```bash
#!/bin/bash

echo "Starting runtime benchmarks..."
echo "This will take several hours to days depending on dataset sizes"
echo "Start time: $(date)"

# Function to run and time a single analysis
run_benchmark() {
    input_file=$1
    echo ""
    echo "=========================================="
    echo "Running: $input_file"
    echo "Time: $(date)"
    echo "=========================================="

    # Extract info from filename
    basename=$(basename "$input_file" .tsv)
    output_dir="benchmarking/runtime/output/${basename}"

    # Run with time tracking
    time boldgenotyper "$input_file" \
        --output-dir "$output_dir" \
        --build-tree \
        2>&1 | tee "benchmarking/runtime/results/${basename}.log"

    echo "Completed: $basename at $(date)"
}

export -f run_benchmark

# Run all Panulirus (should be quick)
echo "=== Starting Panulirus benchmarks ==="
find benchmarking/runtime/Panulirus -name "*.tsv" -type f | sort | while read f; do
    run_benchmark "$f"
done

# Run all Sphyrnidae
echo "=== Starting Sphyrnidae benchmarks ==="
find benchmarking/runtime/Sphyrnidae -name "*.tsv" -type f | sort | while read f; do
    run_benchmark "$f"
done

# Run all Carcharhinidae (will take longest)
echo "=== Starting Carcharhinidae benchmarks ==="
find benchmarking/runtime/Carcharhinidae -name "*.tsv" -type f | sort | while read f; do
    run_benchmark "$f"
done

# Run all Salmonidae
echo "=== Starting Salmonidae benchmarks ==="
find benchmarking/runtime/Salmonidae -name "*.tsv" -type f | sort | while read f; do
    run_benchmark "$f"
done

echo ""
echo "=========================================="
echo "All benchmarks complete!"
echo "End time: $(date)"
echo "=========================================="
```

Run it:
```bash
chmod +x scripts/run_all_benchmarks.sh
nohup ./scripts/run_all_benchmarks.sh > benchmarking/runtime/benchmark_run.log 2>&1 &
```

**Monitor progress**:
```bash
tail -f benchmarking/runtime/benchmark_run.log
```

---

## Step 6: Run Parameter Sweeps (While Benchmarks Run)

These are independent and can run in parallel:

```bash
# Run in separate terminal or after runtime benchmarks
boldgenotyper-sweep data/Sphyrnidae.tsv \
    --output-dir benchmarking/parameter_sweep/Sphyrnidae \
    --min-threshold 0.01 --max-threshold 0.10 --step 0.005

boldgenotyper-sweep data/Panulirus.tsv \
    --output-dir benchmarking/parameter_sweep/Panulirus \
    --min-threshold 0.01 --max-threshold 0.10 --step 0.005

boldgenotyper-sweep data/Salmonidae.tsv \
    --output-dir benchmarking/parameter_sweep/Salmonidae \
    --min-threshold 0.01 --max-threshold 0.08 --step 0.005

boldgenotyper-sweep data/Pieridae.tsv \
    --output-dir benchmarking/parameter_sweep/Pieridae \
    --min-threshold 0.01 --max-threshold 0.10 --step 0.005
```

---

## Step 7: Prepare Accuracy Validation Data

For Sphyrnidae, you need to prepare ground truth from your manual curation:

```python
# scripts/prepare_sphyrnidae_ground_truth.py
import pandas as pd

# Load your manually curated data (from PLOS ONE manuscript)
# Adjust path/columns as needed
manual = pd.read_csv("path/to/your/manual_curation.csv")

# Should have columns: processid, manual_haplotype_id
# Format for comparison
ground_truth = manual[['processid', 'manual_haplotype_id']].copy()
ground_truth.columns = ['processid', 'manual_haplotype']

# Save
ground_truth.to_csv("benchmarking/accuracy/sphyrnidae_manual_curation.csv", index=False)
print(f"Saved ground truth for {len(ground_truth)} samples")
```

---

## Step 8: Run Comparative Analysis

```bash
# First run species-level (if not already done)
boldgenotyper data/Sphyrna_lewini.tsv \
    --output-dir data/Sphyrna_lewini_output \
    --clustering-threshold 0.025

# Then family-level
boldgenotyper data/Sphyrnidae.tsv \
    --output-dir data/Sphyrnidae_output \
    --clustering-threshold 0.030

# Finally comparative analysis
boldgenotyper-compare \
    --species-level data/Sphyrna_lewini_output \
    --family-level data/Sphyrnidae_output \
    --output-dir benchmarking/comparative_analysis/Sphyrnidae_comparison
```

---

## Step 9: Analyze Results

Once benchmarks complete, run analysis scripts from `BENCHMARKING_PROTOCOL.md`:

1. `scripts/analyze_runtime_benchmarks.py` → Figures + tables
2. `scripts/calculate_accuracy_metrics.py` → Accuracy stats
3. `scripts/analyze_parameter_sweeps.py` → Sweep validation
4. `scripts/analyze_comparative_results.py` → QC metrics

---

## Timeline Estimate

| Phase | Time Required | Can Run in Parallel? |
|-------|---------------|---------------------|
| Download datasets | 1-2 hours | No (sequential downloads) |
| Quick test runs | 2-3 hours | Yes (all 5 in parallel) |
| Create subsamples | 10 minutes | No |
| Runtime benchmarks | 1-3 days | No (sequential for accuracy) |
| Parameter sweeps | 4-12 hours | Yes (all 4 in parallel) |
| Accuracy validation | 2-4 hours | No (manual curation prep) |
| Comparative analysis | 1-2 hours | No |
| Data analysis | 4-6 hours | No |
| **Total** | **~1-2 weeks** | **with parallelization** |

---

## Quick Troubleshooting

**Problem**: Downloads too slow from BOLD
- **Solution**: Download during off-peak hours (evening/weekend)
- **Alternative**: Contact BOLD Systems for bulk data access

**Problem**: Benchmarks taking too long
- **Solution**: Reduce replicate count from 3 to 1
- **Alternative**: Skip largest size categories (e.g., 10,000 sequences)

**Problem**: Running out of memory
- **Solution**: Close other applications, use swap space
- **Alternative**: Run largest benchmarks on HPC cluster if available

**Problem**: Can't find manual curation data for Sphyrnidae
- **Solution**: Extract from your PLOS ONE manuscript supplementary data
- **Alternative**: Use one of your output runs as "manual" baseline

---

## Checklist

Before starting:
- [ ] All 5 datasets downloaded from BOLD
- [ ] Quick test runs completed successfully
- [ ] Directory structure created
- [ ] Subsampling script ready
- [ ] Enough disk space (recommend 50+ GB free)
- [ ] Time allocated for long-running benchmarks

After benchmarking:
- [ ] Runtime data collected for all size categories
- [ ] Parameter sweep completed for all datasets
- [ ] Accuracy validation data prepared
- [ ] Comparative analysis results generated
- [ ] All analysis scripts run successfully
- [ ] Figures and tables generated

---

## Next Steps After Benchmarking

Once all benchmarking is complete:

1. **Review results** - Check if metrics meet expectations
2. **Generate manuscript figures** - Use analysis scripts
3. **Write Methods section** - Use templates from BENCHMARKING_PROTOCOL.md
4. **Write Results section** - Report key findings
5. **Create supplementary materials** - Tables and extended results

---

## Questions Before You Start?

Review these documents for details:
- `CASE_STUDY_DATASET_RECOMMENDATIONS.md` - Dataset selection rationale
- `BENCHMARKING_PROTOCOL.md` - Complete protocols and analysis
- `MANUSCRIPT_PREPARATION_PLAN.md` - Overall publication strategy

**Ready to begin!** Start with Step 1: Download datasets from BOLD.

---

**Document created**: 2026-01-28
**Estimated time to complete all benchmarking**: 1-2 weeks
**Questions**: Review full protocol documents or ask for clarification
