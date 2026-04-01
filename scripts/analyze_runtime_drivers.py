#!/usr/bin/env python3
"""
Identify what actually drives BOLDGenotyper runtime and produce
publication-ready figures demonstrating the relationship.

Input:  benchmarking/runtime/benchmark_details.csv
        (produced by scripts/parse_benchmark_details.py)

Output figures:
  Figure_Runtime_Drivers.png/pdf     – runtime vs haplotype count vs input size
  Figure_Phase_Breakdown.png/pdf     – where does the time actually go?
  Figure_Input_vs_Haplotypes.png/pdf – why some organisms scale cleanly

Output table:
  benchmarking/runtime/runtime_driver_correlations.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')

df = pd.read_csv('benchmarking/runtime/benchmark_details.csv')
organisms = sorted(df['organism'].unique())
palette = dict(zip(organisms, sns.color_palette("Set2", n_colors=len(organisms))))

# ---------------------------------------------------------------------------
# Figure 1: What drives runtime?  n_sequences vs n_haplotypes_filtered
#
# Three panels, same y-axis (total runtime), different x-axes.
# Left:  runtime vs input sequences  (what we plotted before)
# Middle: runtime vs filtered haplotypes (the real driver)
# Right: runtime vs phylo subset size (driver of the high-variance tail)
#
# Each organism gets its own colour and a linear trend on the group means.
# ---------------------------------------------------------------------------

means = (
    df.groupby(['organism', 'n_sequences'])
    .agg(
        mean_time=('total_time_sec', 'mean'),
        mean_h_filt=('n_haplotypes_filtered', 'mean'),
        mean_h_phylo=('n_haplotypes_phylo', 'mean'),
    )
    .reset_index()
)

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

x_configs = [
    ('n_sequences',          'Input Sequences',           'mean_time'),
    ('n_haplotypes_filtered','Filtered Haplotypes (h)',   'mean_time'),
    ('n_haplotypes_phylo',   'Phylogenetic Subset Size',  'mean_time'),
]

for ax, (xcol, xlabel, _) in zip(axes, x_configs):
    # Raw points
    sns.scatterplot(
        data=df, x=xcol, y='total_time_sec',
        hue='organism', style='organism',
        alpha=0.45, s=50, ax=ax, legend=(ax is axes[0])
    )
    # Per-organism linear trend on group means
    for org in organisms:
        grp = df[df['organism'] == org].dropna(subset=[xcol, 'total_time_sec'])
        if len(grp) < 3:
            continue
        X = grp[xcol].values.reshape(-1, 1)
        y = grp['total_time_sec'].values
        model = LinearRegression().fit(X, y)

        x_line = np.linspace(X.min(), X.max(), 80).reshape(-1, 1)
        y_line = np.clip(model.predict(x_line), 0, None)
        ax.plot(x_line, y_line, color=palette[org], linewidth=1.8)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.grid(True, alpha=0.25)

axes[0].set_ylabel('Total Runtime (seconds)', fontsize=11)
axes[0].set_title('A. vs Input Size', fontsize=12, fontweight='bold')
axes[1].set_title('B. vs Haplotype Count', fontsize=12, fontweight='bold')
axes[2].set_title('C. vs Phylogenetic Subset', fontsize=12, fontweight='bold')
axes[0].legend(title='Organism', fontsize=7, title_fontsize=8, loc='upper left')

plt.tight_layout()
plt.savefig('benchmarking/runtime/Figure_Runtime_Drivers.png', dpi=300, bbox_inches='tight')
plt.savefig('benchmarking/runtime/Figure_Runtime_Drivers.pdf', bbox_inches='tight')
print("Saved: Figure_Runtime_Drivers.png/pdf")

# ---------------------------------------------------------------------------
# Figure 2: Where does the time go? (phase breakdown)
#
# Stacked horizontal bar chart: one bar per organism, showing mean seconds
# in each phase.  Phases grouped into three logical buckets for clarity.
# ---------------------------------------------------------------------------

phase_cols = {
    'Data & QC':                ['phase_1_sec', 'phase_2_sec'],
    'Alignment & Dereplication':['phase_3_sec', 'phase_4_sec', 'phase_5_sec', 'phase_6_sec', 'phase_6_5_sec'],
    'Geography':                ['phase_7_sec'],
    'Phylogenetics':            ['phase_8_sec'],
    'Divergence':               ['phase_9_sec'],
    'Metadata':                 ['phase_9_5_sec'],
    'Visualization':            ['phase_10_sec'],
}

# Aggregate: sum the raw phase columns, then take mean per organism
phase_sums = pd.DataFrame()
for bucket, cols in phase_cols.items():
    existing = [c for c in cols if c in df.columns]
    if existing:
        phase_sums[bucket] = df[existing].sum(axis=1)

phase_sums['organism'] = df['organism'].values
phase_means = phase_sums.groupby('organism').mean()

# Sort organisms by total mean runtime (ascending) for a cleaner chart
phase_means['_total'] = phase_means.sum(axis=1)
phase_means = phase_means.sort_values('_total').drop(columns='_total')

bucket_colors = {
    'Data & QC':                 '#4e79a7',
    'Alignment & Dereplication': '#f28e2b',
    'Geography':                 '#e15759',
    'Phylogenetics':             '#76b7b2',
    'Divergence':                '#59a14f',
    'Metadata':                  '#edc948',
    'Visualization':             '#b07aa1',
}

fig, ax = plt.subplots(figsize=(10, 5))
phase_means.plot(
    kind='barh', stacked=True, ax=ax,
    color=[bucket_colors[b] for b in phase_means.columns],
    edgecolor='white', linewidth=0.5
)
ax.set_xlabel('Mean Runtime (seconds)', fontsize=11)
ax.set_ylabel('')
ax.set_title('Runtime by Pipeline Phase', fontsize=13, fontweight='bold')
ax.legend(title='Phase', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
ax.grid(True, alpha=0.25, axis='x')

plt.tight_layout()
plt.savefig('benchmarking/runtime/Figure_Phase_Breakdown.png', dpi=300, bbox_inches='tight')
plt.savefig('benchmarking/runtime/Figure_Phase_Breakdown.pdf', bbox_inches='tight')
print("Saved: Figure_Phase_Breakdown.png/pdf")

# ---------------------------------------------------------------------------
# Figure 3: Why do some organisms scale cleanly?
#
# Show input size vs haplotype count per organism.  Where h grows linearly
# with n, runtime is predictable.  Where it doesn't, runtime is noisy.
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: n_sequences vs n_haplotypes_filtered
ax = axes[0]
sns.scatterplot(
    data=df, x='n_sequences', y='n_haplotypes_filtered',
    hue='organism', style='organism',
    alpha=0.55, s=55, ax=ax
)
for org in organisms:
    grp = df[df['organism'] == org].dropna(subset=['n_sequences', 'n_haplotypes_filtered'])
    if len(grp) < 3:
        continue
    X = grp['n_sequences'].values.reshape(-1, 1)
    y = grp['n_haplotypes_filtered'].values
    model = LinearRegression().fit(X, y)
    x_line = np.linspace(X.min(), X.max(), 80).reshape(-1, 1)
    ax.plot(x_line, np.clip(model.predict(x_line), 0, None),
            color=palette[org], linewidth=1.8)

ax.set_xlabel('Input Sequences', fontsize=11)
ax.set_ylabel('Filtered Haplotypes', fontsize=11)
ax.set_title('A. Input Size → Haplotype Count', fontsize=12, fontweight='bold')
ax.legend(title='Organism', fontsize=7, title_fontsize=8)
ax.grid(True, alpha=0.25)

# Panel B: n_haplotypes_filtered vs n_haplotypes_phylo
ax = axes[1]
plot_df = df.dropna(subset=['n_haplotypes_filtered', 'n_haplotypes_phylo'])
sns.scatterplot(
    data=plot_df, x='n_haplotypes_filtered', y='n_haplotypes_phylo',
    hue='organism', style='organism',
    alpha=0.55, s=55, ax=ax, legend=False
)
for org in organisms:
    grp = plot_df[plot_df['organism'] == org]
    if len(grp) < 3:
        continue
    X = grp['n_haplotypes_filtered'].values.reshape(-1, 1)
    y = grp['n_haplotypes_phylo'].values
    model = LinearRegression().fit(X, y)
    x_line = np.linspace(X.min(), X.max(), 80).reshape(-1, 1)
    ax.plot(x_line, np.clip(model.predict(x_line), 0, None),
            color=palette[org], linewidth=1.8)

ax.set_xlabel('Filtered Haplotypes', fontsize=11)
ax.set_ylabel('Phylogenetic Subset Size', fontsize=11)
ax.set_title('B. Haplotype Count → Tree Size', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig('benchmarking/runtime/Figure_Input_vs_Haplotypes.png', dpi=300, bbox_inches='tight')
plt.savefig('benchmarking/runtime/Figure_Input_vs_Haplotypes.pdf', bbox_inches='tight')
print("Saved: Figure_Input_vs_Haplotypes.png/pdf")

# ---------------------------------------------------------------------------
# Statistical summary: R² comparison table
#
# For each organism (and pooled), fit runtime ~ n_sequences and
# runtime ~ n_haplotypes_filtered.  Print and save the comparison.
# ---------------------------------------------------------------------------

predictor_pairs = [
    ('n_sequences',          'Input Sequences'),
    ('n_haplotypes_filtered','Filtered Haplotypes'),
    ('n_haplotypes_phylo',   'Phylo Subset Size'),
]

rows = []
for org in ['[Pooled]'] + organisms:
    subset = df if org == '[Pooled]' else df[df['organism'] == org]
    for xcol, label in predictor_pairs:
        clean = subset.dropna(subset=[xcol, 'total_time_sec'])
        if len(clean) < 3:
            continue
        X = clean[xcol].values.reshape(-1, 1)
        y = clean['total_time_sec'].values
        model = LinearRegression().fit(X, y)
        r2 = model.score(X, y)
        rows.append({
            'organism': org,
            'predictor': label,
            'R²': round(r2, 3),
            'n_runs': len(clean),
        })

corr_df = pd.DataFrame(rows)
# Pivot for readability
pivot = corr_df.pivot_table(
    index='organism', columns='predictor', values='R²', aggfunc='first'
)[['Input Sequences', 'Filtered Haplotypes', 'Phylo Subset Size']]

print("\n--- R² of total_time ~ predictor (linear fit) ---")
print(pivot.to_string())

corr_df.to_csv('benchmarking/runtime/runtime_driver_correlations.csv', index=False)
print("\nSaved: runtime_driver_correlations.csv")

# ---------------------------------------------------------------------------
# Print the key takeaway numbers for the manuscript
# ---------------------------------------------------------------------------
print("\n--- Key numbers for manuscript ---")
pooled = corr_df[corr_df['organism'] == '[Pooled]'].set_index('predictor')['R²']
print(f"Pooled R² (runtime ~ input sequences):   {pooled['Input Sequences']}")
print(f"Pooled R² (runtime ~ filtered haplotypes): {pooled['Filtered Haplotypes']}")
print(f"Pooled R² (runtime ~ phylo subset):       {pooled['Phylo Subset Size']}")

# Best-fitting organism per predictor
for pred in ['Input Sequences', 'Filtered Haplotypes']:
    best = (
        corr_df[(corr_df['predictor'] == pred) & (corr_df['organism'] != '[Pooled]')]
        .sort_values('R²', ascending=False)
        .iloc[0]
    )
    print(f"\nBest R² for '{pred}': {best['organism']} (R²={best['R²']})")

# ---------------------------------------------------------------------------
# Salmonidae vs Carcharhiniformes phase-breakdown comparison
# (manuscript discussion point: dataset structure drives runtime, not just size)
# ---------------------------------------------------------------------------
print("\n--- Salmonidae vs Carcharhiniformes phase breakdown (manuscript) ---")
for org in ['Salmonidae', 'Carcharhiniformes']:
    if org not in phase_means.index:
        continue
    row = phase_means.loc[org]
    total = row.sum()
    aln_pct   = row.get('Alignment & Dereplication', 0) / total * 100
    meta_pct  = row.get('Metadata', 0) / total * 100
    viz_pct   = row.get('Visualization', 0) / total * 100
    geo_pct   = row.get('Geography', 0) / total * 100
    print(f"\n  {org} (mean total across all subsample sizes: {total:.0f}s)")
    print(f"    Alignment & Dereplication: {aln_pct:.1f}%")
    print(f"    Geography:                 {geo_pct:.1f}%")
    print(f"    Metadata:                  {meta_pct:.1f}%")
    print(f"    Visualization:             {viz_pct:.1f}%")

# Per-haplotype singleton rate: singletons removed / total haplotypes before filtering
for org in ['Salmonidae', 'Carcharhiniformes']:
    grp = df[df['organism'] == org]
    raw   = grp['n_haplotypes_raw'].mean()
    filt  = grp['n_haplotypes_filtered'].mean()
    if raw > 0:
        pct_removed = (raw - filt) / raw * 100
        print(f"\n  {org}: mean raw haplotypes={raw:.0f}, filtered={filt:.0f}, "
              f"singleton removal={pct_removed:.1f}% of raw")
