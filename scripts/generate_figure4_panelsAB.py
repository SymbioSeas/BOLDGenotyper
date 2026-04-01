#!/usr/bin/env python3
"""
Generate Figure 4 Panels A and B for the BOLDGenotyper MER manuscript.

Panel A — 5-panel elbow plot grid (one sub-panel per dataset):
    Haplotype count vs. min_singleton_distance threshold.
    Elbow point annotated with dashed vertical line and label.
    Sub-panels colored by environment type.
    Standardised x-axis across all sub-panels for direct comparison.

Panel B — Normalised threshold sensitivity plot (all 5 datasets overlaid):
    Haplotype count as a fraction of the count at the minimum threshold tested.
    Datasets coloured by environment; Pieridae clearly diverges.
    Shows WHY per-dataset optimisation is necessary.

Output:
    benchmarking/parameter_sweep/Figure4_PanelA.png
    benchmarking/parameter_sweep/Figure4_PanelA.pdf
    benchmarking/parameter_sweep/Figure4_PanelB.png
    benchmarking/parameter_sweep/Figure4_PanelB.pdf
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Dataset metadata
# ---------------------------------------------------------------------------
DATASETS = [
    dict(name='Sphyrnidae',       label='Sphyrnidae\n(family, marine)',       environment='Marine'),
    dict(name='Panulirus',         label='Panulirus\n(genus, marine)',         environment='Marine'),
    dict(name='Carcharhiniformes', label='Carcharhiniformes\n(order, marine)', environment='Marine'),
    dict(name='Salmonidae',        label='Salmonidae\n(family, freshwater)',   environment='Freshwater'),
    dict(name='Pieridae',          label='Pieridae\n(family, terrestrial)',    environment='Terrestrial'),
]

ENV_COLORS = {
    'Marine':      '#3182bd',
    'Freshwater':  '#31a354',
    'Terrestrial': '#a63603',
}

ENV_LIGHT = {
    'Marine':      '#deebf7',
    'Freshwater':  '#e5f5e0',
    'Terrestrial': '#fee6ce',
}

SWEEP_ROOT = Path('benchmarking/parameter_sweep')
OUT_ROOT   = SWEEP_ROOT

# ---------------------------------------------------------------------------
# Load sweep data and parse elbow threshold for each dataset
# ---------------------------------------------------------------------------
def load_sweep(name: str) -> pd.DataFrame:
    return pd.read_csv(SWEEP_ROOT / name / 'sweep_summary.csv')


def parse_elbow(name: str) -> float:
    text = (SWEEP_ROOT / name / 'recommendations.txt').read_text()
    return float(re.search(r'PRIMARY:\s+min_singleton_distance\s*=\s*([0-9.]+)', text).group(1))


for ds in DATASETS:
    ds['df']    = load_sweep(ds['name'])
    ds['elbow'] = parse_elbow(ds['name'])

# ---------------------------------------------------------------------------
# PANEL A — 5-panel elbow plot grid
# ---------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')

# 3 columns × 2 rows; 5th dataset goes in row 2 col 2 (centre)
fig_a, axes = plt.subplots(2, 3, figsize=(12, 7))
axes_flat   = axes.flatten()

# Hide the 6th (empty) subplot
axes_flat[5].set_visible(False)
# Move 4th and 5th datasets into cols 1 and 2 of row 2 (indices 3 and 4)
# Default layout: [0]=Sphyrnidae, [1]=Panulirus, [2]=Carcharhiniformes,
#                 [3]=Salmonidae, [4]=Pieridae, [5]=empty
# That's already fine — Salmonidae and Pieridae in row 2 cols 0 and 1,
# empty in col 2.

for i, ds in enumerate(DATASETS):
    ax   = axes_flat[i]
    df   = ds['df']
    name = ds['name']
    env  = ds['environment']
    col  = ENV_COLORS[env]
    bg   = ENV_LIGHT[env]
    elbow_t = ds['elbow']

    thresholds   = df['threshold'].values
    n_haplotypes = df['n_haplotypes'].values
    n_singletons = df['n_singletons'].values

    # Elbow index
    elbow_idx = np.where(thresholds == elbow_t)[0]
    elbow_n   = n_haplotypes[elbow_idx[0]] if len(elbow_idx) > 0 else None

    # Background shading by environment
    ax.set_facecolor(bg)

    # Haplotype count line
    ax.plot(thresholds, n_haplotypes, color=col, linewidth=2.0,
            marker='o', markersize=4.5, zorder=3, label='Haplotypes')

    # Singleton count line (secondary, thinner)
    ax.plot(thresholds, n_singletons, color=col, linewidth=1.0,
            marker='s', markersize=3, alpha=0.55, linestyle='--',
            zorder=2, label='Singletons retained')

    # Elbow annotation
    if elbow_n is not None:
        ax.axvline(elbow_t, color='#444444', linewidth=1.1, linestyle=':', zorder=4)
        ax.scatter([elbow_t], [elbow_n], color='white', edgecolors=col,
                   s=60, zorder=5, linewidths=2.0)
        ax.text(elbow_t + (thresholds[-1] - thresholds[0]) * 0.025,
                elbow_n,
                f'  {elbow_t:.3f}\n  ({elbow_n:,} hap.)',
                fontsize=7.5, color='#222222', va='center', zorder=6)

    ax.set_title(ds['label'], fontsize=9.5, fontweight='bold', pad=4)
    ax.set_xlabel('min_singleton_distance', fontsize=8.5)
    ax.set_ylabel('Count', fontsize=8.5)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    # Environment legend (small patch in corner)
    env_patch = mpatches.Patch(facecolor=col, edgecolor='white', label=env)
    ax.legend(handles=[env_patch], loc='upper right', fontsize=7,
              framealpha=0.85, handlelength=1.2)

fig_a.suptitle(
    'A. Haplotype count vs. singleton-distance threshold\n'
    '(elbow point = optimal min_singleton_distance)',
    fontsize=11, fontweight='bold', y=1.01
)
plt.tight_layout()

out_a_png = OUT_ROOT / 'Figure4_PanelA.png'
out_a_pdf = OUT_ROOT / 'Figure4_PanelA.pdf'
fig_a.savefig(str(out_a_png), dpi=300, bbox_inches='tight')
fig_a.savefig(str(out_a_pdf), bbox_inches='tight')
print(f"Saved: {out_a_png}")
print(f"Saved: {out_a_pdf}")

# ---------------------------------------------------------------------------
# PANEL B — Normalised threshold sensitivity (all datasets overlaid)
# ---------------------------------------------------------------------------
fig_b, ax_b = plt.subplots(figsize=(7, 4.5))

for ds in DATASETS:
    df   = ds['df']
    env  = ds['environment']
    col  = ENV_COLORS[env]

    thresholds   = df['threshold'].values
    n_haplotypes = df['n_haplotypes'].values
    norm         = n_haplotypes / n_haplotypes[0]   # normalise to min threshold

    elbow_t = ds['elbow']
    elbow_idx = np.where(thresholds == elbow_t)[0]
    elbow_norm = norm[elbow_idx[0]] if len(elbow_idx) > 0 else None

    linestyle = '--' if ds['name'] == 'Pieridae' else '-'
    ax_b.plot(thresholds, norm,
              color=col, linewidth=2.0, linestyle=linestyle,
              marker='o', markersize=4,
              label=f"{ds['name']} ({env})", zorder=3)

    # Mark elbow
    if elbow_norm is not None:
        ax_b.scatter([elbow_t], [elbow_norm],
                     color='white', edgecolors=col,
                     s=55, zorder=5, linewidths=2.0)

# Reference line at the vertebrate consensus threshold 0.015
ax_b.axvline(0.015, color='#888888', linewidth=0.9, linestyle=':', alpha=0.7, zorder=2)
ax_b.text(0.015 + 0.001, 1.003, 'consensus\n(0.015)',
          fontsize=7.5, color='#888888', va='bottom')

# Pieridae annotation
pie_df = next(d for d in DATASETS if d['name'] == 'Pieridae')
pie_thresholds = pie_df['df']['threshold'].values
pie_norm = pie_df['df']['n_haplotypes'].values / pie_df['df']['n_haplotypes'].values[0]
ax_b.annotate(
    'Pieridae\n(terrestrial)',
    xy=(pie_thresholds[-1], pie_norm[-1]),
    xytext=(pie_thresholds[-1] - 0.025, pie_norm[-1] - 0.06),
    fontsize=8, color=ENV_COLORS['Terrestrial'],
    arrowprops=dict(arrowstyle='->', color=ENV_COLORS['Terrestrial'], lw=1.1),
)

ax_b.set_xlabel('min_singleton_distance', fontsize=10)
ax_b.set_ylabel('Haplotype count (fraction of minimum threshold)', fontsize=10)
ax_b.set_title(
    'B. Normalised threshold sensitivity across datasets\n'
    '(open circles = elbow / optimal threshold; dashed = Pieridae)',
    fontsize=11, fontweight='bold'
)
ax_b.set_ylim(top=ax_b.get_ylim()[1] * 1.06)  # headroom for annotations
ax_b.grid(True, alpha=0.3)

# Compact legend: one entry per environment
legend_patches = [
    mpatches.Patch(facecolor=ENV_COLORS['Marine'],      edgecolor='white', label='Marine'),
    mpatches.Patch(facecolor=ENV_COLORS['Freshwater'],  edgecolor='white', label='Freshwater'),
    mpatches.Patch(facecolor=ENV_COLORS['Terrestrial'], edgecolor='white', label='Terrestrial'),
]
ax_b.legend(handles=legend_patches + [
    plt.Line2D([0], [0], color='#888888', linewidth=1.5, linestyle=':',
               label='consensus threshold (0.015)'),
],
    loc='lower left', fontsize=8.5, framealpha=0.9,
    title='Environment', title_fontsize=9,
)

plt.tight_layout()

out_b_png = OUT_ROOT / 'Figure4_PanelB.png'
out_b_pdf = OUT_ROOT / 'Figure4_PanelB.pdf'
fig_b.savefig(str(out_b_png), dpi=300, bbox_inches='tight')
fig_b.savefig(str(out_b_pdf), bbox_inches='tight')
print(f"Saved: {out_b_png}")
print(f"Saved: {out_b_pdf}")
