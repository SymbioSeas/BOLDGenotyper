# scripts/analyze_parameter_sweeps.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

datasets = ['Sphyrnidae', 'Panulirus', 'Salmonidae', 'Pieridae', 'Carcharhiniformes']

# ---------------------------------------------------------------------------
# Load sweep summaries and enrich with per-run haplotype counts
# ---------------------------------------------------------------------------
all_data = {}

for dataset in datasets:
    summary_path = Path(f"benchmarking/parameter_sweep/{dataset}/sweep_summary.csv")
    if not summary_path.exists():
        print(f"Warning: {summary_path} not found")
        continue

    df = pd.read_csv(summary_path)

    # Derive haplotype counts and singleton counts from per-run stats files
    n_haplotypes = []
    n_singletons = []
    for _, row in df.iterrows():
        threshold_label = f"{row['threshold']:.3f}".replace('.', '_')
        stats_path = (
            Path(f"benchmarking/parameter_sweep/{dataset}/runs")
            / f"threshold_{threshold_label}"
            / "haplotypes"
            / f"{dataset}_haplotype_stats.csv"
        )
        if stats_path.exists():
            stats = pd.read_csv(stats_path)
            n_haplotypes.append(len(stats))
            n_singletons.append(stats['is_singleton'].sum() if 'is_singleton' in stats.columns else None)
        else:
            n_haplotypes.append(None)
            n_singletons.append(None)

    df['n_haplotypes'] = n_haplotypes
    df['n_singletons'] = n_singletons
    all_data[dataset] = df

if not all_data:
    print("ERROR: No sweep results found. Run scripts/run_parameter_sweeps.sh first.")
    exit(1)

# ---------------------------------------------------------------------------
# Figure 1: Per-organism haplotype count vs threshold (2x3 grid)
#
# X-axis is standardized to 0.01–0.10 across all panels so elbow positions
# are visually comparable.
# ---------------------------------------------------------------------------

# Compute shared x-axis limits across all datasets
all_thresholds = [d['threshold'] for d in all_data.values()]
x_min = min(t.min() for t in all_thresholds)
x_max = 0.10   # standardize upper bound for cross-dataset comparison

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()

summary_data = []

for i, dataset in enumerate(datasets):
    if dataset not in all_data:
        continue

    df = all_data[dataset]
    ax = axes[i]

    ax.plot(df['threshold'], df['n_haplotypes'], 'o-', linewidth=2, markersize=6, color='steelblue')
    ax.plot(df['threshold'], df['n_singletons'], 's--', linewidth=1.5, markersize=5, color='coral', alpha=0.7)

    ax.set_xlim(x_min - 0.002, x_max + 0.002)
    ax.set_xlabel('min_singleton_distance', fontsize=10)
    ax.set_ylabel('Count', fontsize=10)
    ax.set_title(f'{dataset}', fontsize=12, fontweight='bold')
    ax.legend(['Haplotypes', 'Singletons'], fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)

    summary_data.append({
        'dataset': dataset,
        'n_sequences': df['n_samples'].iloc[0],
        'threshold_range_tested': f"{df['threshold'].min()}-{df['threshold'].max()}",
        'min_haplotypes': int(df['n_haplotypes'].min()),
        'max_haplotypes': int(df['n_haplotypes'].max()),
    })

# Hide empty subplot
axes[len(datasets)].axis('off')

plt.suptitle('Haplotype Count vs Singleton Distance Threshold', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('benchmarking/parameter_sweep/Figure_Parameter_Sweep_All_Datasets.png',
            dpi=300, bbox_inches='tight')
plt.savefig('benchmarking/parameter_sweep/Figure_Parameter_Sweep_All_Datasets.pdf',
            bbox_inches='tight')
print("Saved: Figure_Parameter_Sweep_All_Datasets.png/pdf")

# ---------------------------------------------------------------------------
# Figure 2: Cross-organism comparison
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Absolute haplotype counts
ax = axes[0]
for dataset, df in all_data.items():
    ax.plot(df['threshold'], df['n_haplotypes'], 'o-', linewidth=2, markersize=5, label=dataset, alpha=0.8)
ax.set_xlabel('min_singleton_distance', fontsize=11)
ax.set_ylabel('Number of Haplotypes', fontsize=11)
ax.set_title('A. Absolute Haplotype Counts', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, title='Dataset', title_fontsize=9)
ax.grid(True, alpha=0.3)

# Panel B: Normalized to max (threshold sensitivity)
ax = axes[1]
for dataset, df in all_data.items():
    max_h = df['n_haplotypes'].max()
    ax.plot(df['threshold'], df['n_haplotypes'] / max_h, 'o-', linewidth=2, markersize=5, label=dataset, alpha=0.8)
ax.set_xlabel('min_singleton_distance', fontsize=11)
ax.set_ylabel('Haplotypes (normalized to max)', fontsize=11)
ax.set_title('B. Threshold Sensitivity', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, title='Dataset', title_fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.05)

plt.tight_layout()
plt.savefig('benchmarking/parameter_sweep/Figure_Parameter_Sweep_Summary.png',
            dpi=300, bbox_inches='tight')
plt.savefig('benchmarking/parameter_sweep/Figure_Parameter_Sweep_Summary.pdf',
            bbox_inches='tight')
print("Saved: Figure_Parameter_Sweep_Summary.png/pdf")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv("benchmarking/parameter_sweep/sweep_summary_all_datasets.csv", index=False)

print("\nParameter Sweep Summary:")
print(summary_df.to_string(index=False))
