#!/usr/bin/env python3
"""
Generate Figure 5 for the BOLDGenotyper MER manuscript.

Runtime performance benchmarks across 5 datasets.

Panel A — Runtime scaling: total pipeline time vs. number of input sequences.
           Per-dataset lines with mean ± SD shading. Datasets coloured by
           environment (marine / freshwater / terrestrial); marine datasets
           distinguished by linestyle. Reference annotation for the largest
           tested N value per dataset.

Panel B — Phase breakdown: mean time per pipeline phase, stacked horizontal
           bars, one bar per dataset. Phases grouped into 7 logical buckets.
           Shows that Alignment & Dereplication dominates for high-haplotype-
           count datasets (Carcharhiniformes); Geography / Metadata / Vis
           phases dominate for highly sub-divided datasets (Salmonidae).

Panel C — Runtime efficiency: seconds per input sequence at the largest
           tested dataset size, coloured by environment. Bars annotated with
           the rate and the benchmark N used.

Output:
    benchmarking/runtime/Figure5.png
    benchmarking/runtime/Figure5.pdf
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RT_DIR  = Path('benchmarking/runtime')
RT_CSV  = RT_DIR / 'runtime_benchmark_results.csv'
DET_CSV = RT_DIR / 'benchmark_details.csv'
SUM_CSV = RT_DIR / 'runtime_summary_table.csv'

# ---------------------------------------------------------------------------
# Colors — consistent with Figure 4
# ---------------------------------------------------------------------------
ENV_COLORS = {
    'Marine':      '#3182bd',
    'Freshwater':  '#31a354',
    'Terrestrial': '#a63603',
}

DATASET_ENV = {
    'Carcharhiniformes': 'Marine',
    'Sphyrnidae':        'Marine',
    'Panulirus':         'Marine',
    'Salmonidae':        'Freshwater',
    'Pieridae':          'Terrestrial',
}

# Short labels for in-plot annotations
SHORT_LABEL = {
    'Carcharhiniformes': 'Carcharhiniformes',
    'Sphyrnidae':        'Sphyrnidae',
    'Panulirus':         'Panulirus',
    'Salmonidae':        'Salmonidae',
    'Pieridae':          'Pieridae',
}

# Marine datasets distinguished by linestyle
MARINE_STYLE = {
    'Carcharhiniformes': '-',
    'Sphyrnidae':        '--',
    'Panulirus':         ':',
}

# Phase groupings and colours for Panel B
PHASE_COLS = {
    'Data & QC':                ['phase_1_sec', 'phase_2_sec'],
    'Alignment & Dereplication':['phase_3_sec', 'phase_4_sec', 'phase_5_sec',
                                 'phase_6_sec', 'phase_6_5_sec'],
    'Geography':                ['phase_7_sec'],
    'Phylogenetics':            ['phase_8_sec'],
    'Divergence':               ['phase_9_sec'],
    'Metadata':                 ['phase_9_5_sec'],
    'Visualization':            ['phase_10_sec'],
}

PHASE_COLORS = {
    'Data & QC':                '#90CAF9',
    'Alignment & Dereplication':'#E65100',
    'Geography':                '#43A047',
    'Phylogenetics':            '#8E24AA',
    'Divergence':               '#FFB300',
    'Metadata':                 '#26A69A',
    'Visualization':            '#BDBDBD',
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
rt  = pd.read_csv(RT_CSV)
det = pd.read_csv(DET_CSV)
smry= pd.read_csv(SUM_CSV)

organisms = sorted(rt['run_name'].str.split('_n').str[0].unique())

# Add organism column to raw results
rt['organism'] = rt['run_name'].str.split('_n').str[0]

# Group means + SD per organism × n_sequences
grp = (
    rt[rt['success'] == True]
    .groupby(['organism', 'n_sequences'])
    .agg(
        mean_time=('total_time', 'mean'),
        sd_time=('total_time', 'std'),
    )
    .reset_index()
)

# ---------------------------------------------------------------------------
# Phase aggregation for Panel B
# ---------------------------------------------------------------------------
phase_sums = pd.DataFrame()
for bucket, cols in PHASE_COLS.items():
    existing = [c for c in cols if c in det.columns]
    if existing:
        phase_sums[bucket] = det[existing].sum(axis=1)

phase_sums['organism'] = det['organism'].values

# Use rows at max tested N per organism for phase breakdown (most representative)
max_n_per_org = det.groupby('organism')['n_sequences'].transform('max')
mask = det['n_sequences'] == max_n_per_org
phase_at_max = phase_sums[mask.values].copy()
phase_at_max['organism'] = det.loc[mask, 'organism'].values

phase_means = phase_at_max.groupby('organism')[list(PHASE_COLS.keys())].mean()
# Sort by total runtime ascending (smallest dataset at top of bar chart)
phase_means['_total'] = phase_means.sum(axis=1)
phase_means = phase_means.sort_values('_total').drop(columns='_total')

# ---------------------------------------------------------------------------
# Panel C data: seconds/sequence at largest tested N
# ---------------------------------------------------------------------------
smry['organism'] = smry['organism']
max_n_smry = smry.loc[smry.groupby('organism')['n_sequences'].idxmax()].copy()
max_n_smry['sec_per_seq'] = max_n_smry['mean_time_sec'] / max_n_smry['n_sequences']
max_n_smry['environment'] = max_n_smry['organism'].map(DATASET_ENV)
# Sort by environment then sec_per_seq descending
env_order = {'Marine': 0, 'Freshwater': 1, 'Terrestrial': 2}
max_n_smry['_eo'] = max_n_smry['environment'].map(env_order)
max_n_smry = max_n_smry.sort_values(['_eo', 'sec_per_seq'], ascending=[True, False])

# ---------------------------------------------------------------------------
# Figure layout — 3 columns, unequal widths
# ---------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')
fig = plt.figure(figsize=(15, 5.5))
gs  = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.4, 1.0], wspace=0.28)

ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2])

# ===========================================================================
# PANEL A — Runtime scaling
# ===========================================================================
org_order = ['Carcharhiniformes', 'Salmonidae', 'Pieridae', 'Sphyrnidae', 'Panulirus']

for org in org_order:
    env = DATASET_ENV[org]
    col = ENV_COLORS[env]
    ls  = MARINE_STYLE.get(org, '-')
    sub = grp[grp['organism'] == org].sort_values('n_sequences')
    if sub.empty:
        continue

    # Shaded SD band
    ax_a.fill_between(
        sub['n_sequences'],
        (sub['mean_time'] - sub['sd_time']).clip(lower=0),
        sub['mean_time'] + sub['sd_time'],
        color=col, alpha=0.10,
    )
    # Mean line
    ax_a.plot(
        sub['n_sequences'], sub['mean_time'],
        color=col, linewidth=2.0, linestyle=ls,
        marker='o', markersize=4.5, zorder=3,
    )
    # Only label the two clearly-separated lines inline; rest identified via legend
    last = sub.iloc[-1]
    if org == 'Carcharhiniformes':
        mid = sub.iloc[len(sub) // 2]
        ax_a.text(mid['n_sequences'] + 150, mid['mean_time'] - 50,
                  'Carcharhiniformes', fontsize=7.5, color=col, va='bottom')
    elif org == 'Salmonidae':
        ax_a.text(last['n_sequences'] - 300, last['mean_time'] + 20,
                  'Salmonidae', fontsize=7.5, color=col, va='bottom', ha='right')

ax_a.set_xlabel('Input sequences (N)', fontsize=10)
ax_a.set_ylabel('Total pipeline runtime (seconds)', fontsize=10)
ax_a.set_title('A.  Runtime scaling with dataset size', fontsize=11, fontweight='bold')
ax_a.grid(True, alpha=0.3)

# Legend: environment patches + marine linestyle legend
env_patches = [
    mpatches.Patch(facecolor=ENV_COLORS['Marine'],      label='Marine'),
    mpatches.Patch(facecolor=ENV_COLORS['Freshwater'],  label='Freshwater'),
    mpatches.Patch(facecolor=ENV_COLORS['Terrestrial'], label='Terrestrial'),
]
marine_lines = [
    mlines.Line2D([], [], color='#3182bd', linestyle='-',  linewidth=1.5, label='Carcharhiniformes'),
    mlines.Line2D([], [], color='#3182bd', linestyle='--', linewidth=1.5, label='Sphyrnidae'),
    mlines.Line2D([], [], color='#3182bd', linestyle=':',  linewidth=1.5, label='Panulirus'),
]
ax_a.legend(
    handles=env_patches + marine_lines,
    fontsize=7.5, title='Dataset', title_fontsize=8,
    framealpha=0.9, loc='upper left',
)

# ===========================================================================
# PANEL B — Phase breakdown (at max tested N per dataset)
# ===========================================================================
y_labels = phase_means.index.tolist()
y_pos    = np.arange(len(y_labels))
left     = np.zeros(len(y_labels))

for phase in PHASE_COLS:
    if phase not in phase_means.columns:
        continue
    vals = phase_means[phase].values
    ax_b.barh(
        y_pos, vals, left=left,
        color=PHASE_COLORS[phase],
        height=0.6, edgecolor='white', linewidth=0.5,
        label=phase,
    )
    # Label segments > 25 s to avoid clutter
    for i, (v, l) in enumerate(zip(vals, left)):
        if v > 25:
            ax_b.text(
                l + v / 2, i, f'{v:.0f}s',
                ha='center', va='center',
                fontsize=6.5, color='white', fontweight='bold',
            )
    left = left + vals

# Annotate total runtime at end of each bar
for i, org in enumerate(y_labels):
    total_s = phase_means.loc[org].sum()
    ax_b.text(
        left[i] + 3, i,
        f'{total_s:.0f} s ({total_s/60:.1f} min)',
        va='center', ha='left', fontsize=7.5, color='#333333',
    )

# Colour y-tick labels by environment
for label_obj, org in zip(ax_b.get_yticklabels(), y_labels):
    label_obj.set_color(ENV_COLORS[DATASET_ENV[org]])

ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(y_labels, fontsize=9)
ax_b.set_xlabel('Mean runtime (seconds)', fontsize=10)
ax_b.set_xlim(0, left.max() * 1.30)
ax_b.set_title('B.  Runtime by pipeline phase\n(at largest tested N per dataset)',
               fontsize=11, fontweight='bold')
ax_b.legend(
    loc='lower right', fontsize=7, title='Phase', title_fontsize=7.5,
    framealpha=0.9, ncol=1,
)
ax_b.grid(True, axis='x', alpha=0.3)
ax_b.set_axisbelow(True)

# ===========================================================================
# PANEL C — Runtime efficiency (s/seq at max N)
# ===========================================================================
y_pos_c  = np.arange(len(max_n_smry))
colors_c = [ENV_COLORS[e] for e in max_n_smry['environment']]

ax_c.barh(
    y_pos_c, max_n_smry['sec_per_seq'],
    color=colors_c, height=0.55,
    edgecolor='white', linewidth=0.8, zorder=3,
)

for i, row in enumerate(max_n_smry.itertuples()):
    ax_c.text(
        row.sec_per_seq + 0.002, i,
        f'{row.sec_per_seq:.2f} s/seq\n(n={row.n_sequences:,})',
        va='center', ha='left', fontsize=7.5, color='#222222',
    )

ax_c.set_yticks(y_pos_c)
ax_c.set_yticklabels(max_n_smry['organism'], fontsize=9)
ax_c.set_xlabel('Seconds per input sequence', fontsize=10)
ax_c.set_xlim(0, max_n_smry['sec_per_seq'].max() * 1.55)
ax_c.set_title('C.  Runtime efficiency\n(at largest tested N)',
               fontsize=11, fontweight='bold')
ax_c.grid(True, axis='x', alpha=0.3)
ax_c.set_axisbelow(True)
ax_c.invert_yaxis()

# Colour y-tick labels by environment
for label_obj, org in zip(ax_c.get_yticklabels(), max_n_smry['organism'].tolist()):
    label_obj.set_color(ENV_COLORS[DATASET_ENV[org]])

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
plt.tight_layout()

out_png = RT_DIR / 'Figure5.png'
out_pdf = RT_DIR / 'Figure5.pdf'
fig.savefig(str(out_png), dpi=300, bbox_inches='tight')
fig.savefig(str(out_pdf), bbox_inches='tight')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
