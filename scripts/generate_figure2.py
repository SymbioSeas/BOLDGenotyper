#!/usr/bin/env python3
"""
Generate Figure 2 for the BOLDGenotyper MER manuscript.

Flagship case study: Sphyrnidae family-level analysis.

Panel A — World map: all coordinate-bearing samples plotted by species,
           circle size proportional to co-located sample count.
           GOaS ocean basin boundaries overlaid.
Panel B — Stacked bar chart: species composition by ocean basin
           (basin-assigned samples only, Unknown excluded).

Color palette is shared between both panels.

Output:
    scripts/Figure2.png
    scripts/Figure2.pdf
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
OUT_DIR   = Path('output/Sphyrnidae_full')
SHP_PATH  = Path('shapefiles/GOaS_v1_20211214/goas_v01.shp')
MAP_DATA  = OUT_DIR / 'plots/data/distribution_map.csv'
BAR_DATA  = OUT_DIR / 'plots/data/distribution_bar_relative.csv'

# ---------------------------------------------------------------------------
# Species color palette  (consistent across both panels and Figure 3)
# ---------------------------------------------------------------------------
SPECIES_COLORS = {
    'Sphyrna lewini':    '#1565C0',   # dark blue  — focal species
    'Sphyrna zygaena':  '#E65100',   # deep orange
    'Sphyrna mokarran': '#2E7D32',   # dark green
    'Sphyrna tiburo':   '#C62828',   # dark red
    'Sphyrna tudes':    '#6A1B9A',   # purple
    'Sphyrna corona':   '#4E342E',   # brown
    'Eusphyra blochii': '#00838F',   # teal
    'Sphyrna gilberti': '#9E9E9E',   # grey
}
FALLBACK_COLOR = '#BDBDBD'

# Shortened basin labels for Panel B
BASIN_LABELS = {
    'Indian Ocean':                                'Indian Ocean',
    'North Atlantic Ocean':                        'N. Atlantic',
    'South Atlantic Ocean':                        'S. Atlantic',
    'North Pacific Ocean':                         'N. Pacific',
    'South Pacific Ocean':                         'S. Pacific',
    'South China and Easter Archipelagic Seas':    'S. China Seas',
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
map_df = pd.read_csv(MAP_DATA)
bar_df = pd.read_csv(BAR_DATA)

# Map: only samples with real coordinates
map_plot = map_df[map_df['lat'].notna() & map_df['lon'].notna()].copy()

# Bar: exclude Unknown basin; aggregate to species level
bar_plot = (
    bar_df[bar_df['ocean_basin'] != 'Unknown']
    .groupby(['ocean_basin', 'species'])['n_samples']
    .sum()
    .reset_index()
)
bar_plot['basin_label'] = bar_plot['ocean_basin'].map(BASIN_LABELS).fillna(bar_plot['ocean_basin'])

# Pivot for stacked bars; order basins by total S. lewini count descending
pivot = bar_plot.pivot_table(index='basin_label', columns='species',
                              values='n_samples', aggfunc='sum').fillna(0)
# Sort by S. lewini count (ascending so largest is at top after invert_yaxis)
# Keep all basins even if S. lewini count = 0
lewini_col = pivot.get('Sphyrna lewini', pd.Series(0, index=pivot.index))
pivot = pivot.reindex(lewini_col.sort_values(ascending=True).index)

# Ordered species list for consistent stacking
species_order = [s for s in SPECIES_COLORS if s in pivot.columns]

# ---------------------------------------------------------------------------
# Load GOaS shapefile for basin boundaries
# ---------------------------------------------------------------------------
goas = None
if SHP_PATH.exists():
    goas = gpd.read_file(SHP_PATH)

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')
fig = plt.figure(figsize=(11, 10))
gs  = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.55, 1.0], hspace=0.12)

# ===========================================================================
# PANEL A — World map
# ===========================================================================
ax_map = fig.add_subplot(gs[0], projection=ccrs.Robinson())
ax_map.set_global()
ax_map.add_feature(cfeature.LAND,  facecolor='#F5F5F5', zorder=1)
ax_map.add_feature(cfeature.OCEAN, facecolor='#E3F2FD', zorder=0)
ax_map.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='#9E9E9E', zorder=2)
ax_map.gridlines(draw_labels=False, linewidth=0.3, color='#BDBDBD',
                 alpha=0.5, linestyle='--')

# GOaS basin boundaries
if goas is not None:
    for _, row in goas.iterrows():
        ax_map.add_geometries(
            [row.geometry], crs=ccrs.PlateCarree(),
            facecolor='none', edgecolor='#1565C0',
            linewidth=0.7, linestyle='--', alpha=0.45, zorder=3
        )

# Plot samples, one species at a time so legend entries are ordered
for species in species_order:
    grp = map_plot[map_plot['species'] == species]
    if grp.empty:
        continue
    # Use n_at_location for sizing; if column missing, default size 30
    sizes = grp['n_at_location'].clip(upper=30) * 5 if 'n_at_location' in grp.columns else 30
    genus, epithet = species.split(' ', 1)
    ax_map.scatter(
        grp['lon'], grp['lat'],
        transform=ccrs.PlateCarree(),
        s=sizes,
        color=SPECIES_COLORS[species],
        alpha=0.75, linewidths=0.3, edgecolors='white',
        zorder=5, label=rf'$\it{{{genus}\ {epithet}}}$'
    )

# Size legend (sample count reference)
for ref_n, label in [(1, 'n=1'), (5, 'n=5'), (20, 'n=20')]:
    ax_map.scatter([], [], s=min(ref_n, 30) * 5, color='#757575',
                   alpha=0.7, label=label, edgecolors='white', linewidths=0.3)

handles, labels = ax_map.get_legend_handles_labels()
# Separate species handles from size handles
sp_handles  = [h for h, l in zip(handles, labels) if l.startswith('$')]
sp_labels   = [l for l in labels if l.startswith('$')]
sz_handles  = [h for h, l in zip(handles, labels) if l.startswith('n=')]
sz_labels   = [l for l in labels if l.startswith('n=')]

leg_sp = ax_map.legend(sp_handles, sp_labels, loc='lower left',
                        fontsize=7.5, title='Species', title_fontsize=8.5,
                        framealpha=0.9, ncol=1, markerscale=1.2,
                        handletextpad=0.4)
ax_map.add_artist(leg_sp)
ax_map.legend(sz_handles, sz_labels, loc='lower right',
              fontsize=7.5, title='Sample count', title_fontsize=8,
              framealpha=0.9, handletextpad=0.4)

ax_map.set_title('A.  Geographic distribution of Sphyrnidae COI haplotypes',
                 fontsize=11, fontweight='bold', pad=6)

# ===========================================================================
# PANEL B — Stacked bar chart (species by ocean basin)
# ===========================================================================
ax_bar = fig.add_subplot(gs[1])

# Convert to percentages per basin
row_totals = pivot[species_order].sum(axis=1)
pct_pivot  = pivot[species_order].div(row_totals, axis=0) * 100

basins = pct_pivot.index.tolist()
y_pos  = np.arange(len(basins))
left   = np.zeros(len(basins))

for species in species_order:
    if species not in pct_pivot.columns:
        continue
    vals = pct_pivot[species].values
    genus, epithet = species.split(' ', 1)
    bars = ax_bar.barh(
        y_pos, vals, left=left,
        color=SPECIES_COLORS[species],
        height=0.6, edgecolor='white', linewidth=0.5,
        label=rf'$\it{{{genus}\ {epithet}}}$',
    )
    # Label segments > 8% to avoid clutter
    for i, (v, l) in enumerate(zip(vals, left)):
        if v > 8:
            ax_bar.text(l + v / 2, i, f'{v:.0f}%',
                       ha='center', va='center', fontsize=6.5,
                       color='white', fontweight='bold')
    left = left + vals

# Annotate with total n per basin
for i, basin in enumerate(basins):
    n_total = int(row_totals[basin])
    ax_bar.text(101, i, f'n={n_total}', va='center', ha='left',
                fontsize=7.5, color='#444444')

ax_bar.set_yticks(y_pos)
ax_bar.set_yticklabels(basins, fontsize=9)
ax_bar.set_xlabel('Relative abundance (%)', fontsize=10)
ax_bar.set_xlim(0, 108)
ax_bar.set_title('B.  Species composition by ocean basin (coordinate-assigned samples)',
                 fontsize=11, fontweight='bold')
ax_bar.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0),
              fontsize=7.5, title='Species', title_fontsize=8.5,
              framealpha=0.9, ncol=1)
ax_bar.grid(True, axis='x', alpha=0.3, zorder=0)
ax_bar.set_axisbelow(True)
ax_bar.invert_yaxis()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_png = Path('scripts/Figure2.png')
out_pdf = Path('scripts/Figure2.pdf')
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
