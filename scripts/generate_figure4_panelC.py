#!/usr/bin/env python3
"""
Generate Figure 4 Panel C for the BOLDGenotyper MER manuscript.

Panel C: Summary bar chart of optimal singleton-distance threshold across all
five benchmark datasets, colour-coded by environment type.

Key biological message: marine and freshwater datasets converge on a common
optimal threshold (~0.015), while the terrestrial Pieridae dataset requires a
substantially higher value (0.085), justifying per-dataset parameter optimisation.

Output:
    benchmarking/parameter_sweep/Figure4_PanelC.png
    benchmarking/parameter_sweep/Figure4_PanelC.pdf
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Dataset metadata
# ---------------------------------------------------------------------------
DATASETS = [
    dict(
        name='Carcharhiniformes',
        label='Carcharhiniformes\n(order, marine)',
        environment='Marine',
        n_sequences=9244,
        taxonomic_level='Order',
    ),
    dict(
        name='Sphyrnidae',
        label='Sphyrnidae\n(family, marine)',
        environment='Marine',
        n_sequences=1087,
        taxonomic_level='Family',
    ),
    dict(
        name='Panulirus',
        label='Panulirus\n(genus, marine)',
        environment='Marine',
        n_sequences=755,
        taxonomic_level='Genus',
    ),
    dict(
        name='Salmonidae',
        label='Salmonidae\n(family, freshwater)',
        environment='Freshwater',
        n_sequences=5546,
        taxonomic_level='Family',
    ),
    dict(
        name='Pieridae',
        label='Pieridae\n(family, terrestrial)',
        environment='Terrestrial',
        n_sequences=3484,
        taxonomic_level='Family',
    ),
]

ENV_COLORS = {
    'Marine':      '#3182bd',
    'Freshwater':  '#31a354',
    'Terrestrial': '#a63603',
}

SWEEP_ROOT = Path('benchmarking/parameter_sweep')

# ---------------------------------------------------------------------------
# Parse elbow threshold and haplotype count from recommendations.txt
# ---------------------------------------------------------------------------
def parse_recommendations(dataset_name: str) -> dict:
    path = SWEEP_ROOT / dataset_name / 'recommendations.txt'
    text = path.read_text()

    # Extract elbow threshold: "At 0.015 (elbow):" or "PRIMARY:  min_singleton_distance = 0.085"
    threshold_match = re.search(
        r'PRIMARY:\s+min_singleton_distance\s*=\s*([0-9.]+)', text
    )
    if not threshold_match:
        raise ValueError(f"Could not parse threshold from {path}")
    threshold = float(threshold_match.group(1))

    # Extract haplotype count at elbow
    haplotype_match = re.search(
        r'At [0-9.]+ \(elbow\):\s+([0-9.]+)\s+haplotypes', text
    )
    n_haplotypes = int(float(haplotype_match.group(1))) if haplotype_match else None

    # Extract acceptable range low/high
    range_match = re.search(
        r'ACCEPTABLE RANGE:\s+([0-9.]+)\s+.+?\s+([0-9.]+)', text
    )
    range_low  = float(range_match.group(1)) if range_match else threshold
    range_high = float(range_match.group(2)) if range_match else threshold

    return dict(threshold=threshold, n_haplotypes=n_haplotypes,
                range_low=range_low, range_high=range_high)

# ---------------------------------------------------------------------------
# Assemble data
# ---------------------------------------------------------------------------
for ds in DATASETS:
    ds.update(parse_recommendations(ds['name']))

# Sort: ascending threshold, then alphabetical within same threshold
DATASETS.sort(key=lambda d: (d['threshold'], d['name']))

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')

fig, ax = plt.subplots(figsize=(8, 5))

y_pos  = np.arange(len(DATASETS))
labels = [d['label'] for d in DATASETS]
thresholds = [d['threshold'] for d in DATASETS]
colors = [ENV_COLORS[d['environment']] for d in DATASETS]

bars = ax.barh(
    y_pos, thresholds,
    color=colors,
    height=0.55,
    edgecolor='white', linewidth=0.8,
    zorder=3,
)

# Threshold value annotation at end of each bar
for i, ds in enumerate(DATASETS):
    ax.text(
        ds['threshold'] + 0.0015,
        i,
        f"{ds['threshold']:.3f}",
        va='center', ha='left',
        fontsize=9, fontweight='bold',
        color='#222222',
    )
    # Haplotype count annotation inside/below bar
    if ds['n_haplotypes'] is not None:
        ax.text(
            0.001,
            i - 0.27,
            f"n = {ds['n_haplotypes']:,} haplotypes at elbow",
            va='center', ha='left',
            fontsize=7, color='white' if ds['threshold'] > 0.03 else '#444444',
            style='italic',
        )

# Reference line at 0.015 (consensus marine/freshwater value)
ax.axvline(0.015, color='#555555', linewidth=1.0, linestyle='--', alpha=0.6, zorder=2)
ax.text(0.015 + 0.0008, len(DATASETS) - 0.15,
        'consensus\n(0.015)', fontsize=7.5, color='#555555', va='top')

# Axes formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel('Optimal singleton-distance threshold (min_singleton_distance)', fontsize=10)
ax.set_title('C. Optimal Threshold by Dataset and Environment', fontsize=12, fontweight='bold')
ax.set_xlim(0, 0.105)
ax.grid(True, axis='x', alpha=0.3, zorder=1)
ax.set_axisbelow(True)
ax.invert_yaxis()  # highest threshold at bottom so Pieridae stands out

# Environment legend
legend_patches = [
    mpatches.Patch(facecolor=ENV_COLORS['Marine'],      edgecolor='white', label='Marine'),
    mpatches.Patch(facecolor=ENV_COLORS['Freshwater'],  edgecolor='white', label='Freshwater'),
    mpatches.Patch(facecolor=ENV_COLORS['Terrestrial'], edgecolor='white', label='Terrestrial'),
]
ax.legend(handles=legend_patches, loc='lower right', fontsize=9, framealpha=0.9,
          title='Environment', title_fontsize=9)

plt.tight_layout()

out_base = SWEEP_ROOT / 'Figure4_PanelC'
plt.savefig(str(out_base) + '.png', dpi=300, bbox_inches='tight')
plt.savefig(str(out_base) + '.pdf', bbox_inches='tight')
print(f"Saved: {out_base}.png")
print(f"Saved: {out_base}.pdf")
