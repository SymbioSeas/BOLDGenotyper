#!/usr/bin/env python3
"""
Generate Figure 3 for the BOLDGenotyper MER manuscript.

Flagship case study: Sphyrnidae phylogenetic tree.

Tips are colored by species. The two S. lewini phylogeographic groups
(consistent with Smith & Black 2026) are highlighted with shaded backgrounds.
Support values shown at internal nodes (FastTree local support ≥ 0.70).

Uses the same species color palette as Figure 2.

Output:
    scripts/Figure3.png
    scripts/Figure3.pdf
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio import Phylo

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TREE_PATH = Path('output/Sphyrnidae_full/phylogenetic/Sphyrnidae_tree_relabeled.nwk')

# ---------------------------------------------------------------------------
# Species color palette — identical to generate_figure2.py
# ---------------------------------------------------------------------------
SPECIES_COLORS = {
    'Sphyrna lewini':    '#1565C0',
    'Sphyrna zygaena':  '#E65100',
    'Sphyrna mokarran': '#2E7D32',
    'Sphyrna tiburo':   '#C62828',
    'Sphyrna tudes':    '#6A1B9A',
    'Sphyrna corona':   '#4E342E',
    'Eusphyra blochii': '#00838F',
    'Sphyrna gilberti': '#9E9E9E',
}
FALLBACK_COLOR = '#333333'

# S. lewini phylogeographic groups (parsed from tree topology)
# Group A: predominates Indian Ocean and Atlantic basins
# Group B: predominates western Pacific / S. China Seas
LEWINI_GROUP_A = {'h1', 'h3', 'h9', 'h12', 'h20', 'h21'}
LEWINI_GROUP_B = {'h6', 'h22'}

GROUP_A_COLOR = '#1565C0'   # dark blue
GROUP_B_COLOR = '#64B5F6'   # light blue — same species, distinct clade

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def parse_species(tip_name: str) -> str:
    parts = tip_name.strip().split()
    if len(parts) >= 2 and parts[0][0].isupper():
        return f"{parts[0]} {parts[1]}"
    return tip_name

def parse_haplotype_id(tip_name: str) -> str:
    for part in tip_name.strip().split():
        if part.startswith('h') and '_' in part:
            return part.split('_')[0]
    return ''

def tip_color(tip_name: str) -> str:
    sp  = parse_species(tip_name)
    hid = parse_haplotype_id(tip_name)
    if sp == 'Sphyrna lewini':
        if hid in LEWINI_GROUP_A:
            return GROUP_A_COLOR
        if hid in LEWINI_GROUP_B:
            return GROUP_B_COLOR
    return SPECIES_COLORS.get(sp, FALLBACK_COLOR)

def clean_label(tip_name: str) -> str:
    """'Sphyrna lewini h1_n179' → 'S. lewini  h1  (n=179)'"""
    parts = tip_name.strip().split()
    if len(parts) >= 3 and parts[0][0].isupper():
        abbr    = parts[0][0] + '.'
        epithet = parts[1]
        hap     = parts[2]
        if '_n' in hap:
            hid, n = hap.split('_n', 1)
            return f"{abbr} {epithet}  {hid}  (n={n})"
        return f"{abbr} {epithet}  {hap}"
    return tip_name

# ---------------------------------------------------------------------------
# Load tree — keep original clade names so label_func can map them
# ---------------------------------------------------------------------------
tree = Phylo.read(str(TREE_PATH), 'newick')
tree.ladderize()

# Map original name → clean display label and color
terminals      = tree.get_terminals()
orig_names     = [c.name or '' for c in terminals]
clean_labels   = {n: clean_label(n) for n in orig_names}
color_map      = {n: tip_color(n) for n in orig_names}

# Compute y-positions directly from terminal order.
# Biopython Phylo.draw() assigns y = 1 (bottom) to n_tips (top),
# stepping through get_terminals() from first to last.
n_tips = len(terminals)
tip_y  = {orig_names[i]: (i + 1) for i in range(n_tips)}

group_a_ys = [tip_y[n] for n in orig_names if parse_haplotype_id(n) in LEWINI_GROUP_A]
group_b_ys = [tip_y[n] for n in orig_names if parse_haplotype_id(n) in LEWINI_GROUP_B]

# ---------------------------------------------------------------------------
# Draw tree with clean labels; recolor text objects afterward
# ---------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')
fig, ax = plt.subplots(figsize=(10, 8))

Phylo.draw(
    tree,
    axes=ax,
    label_func=lambda c: clean_labels.get(c.name, '') if c.is_terminal() else '',
    do_show=False,
)

# Recolor tip text objects by matching their text content to our clean labels
clean_to_color = {v: color_map[k] for k, v in clean_labels.items()}
for txt in ax.texts:
    color = clean_to_color.get(txt.get_text())
    if color:
        txt.set_color(color)
        txt.set_fontsize(8.5)

# ---------------------------------------------------------------------------
# Highlight S. lewini phylogeographic groups using computed y-positions
# ---------------------------------------------------------------------------
x_max     = ax.get_xlim()[1]
bracket_x = x_max * 0.68   # shade up to ~2/3 of x-range (before tip labels)

if group_a_ys:
    y0, y1 = min(group_a_ys) - 0.45, max(group_a_ys) + 0.45
    ax.fill_betweenx([y0, y1], 0, bracket_x,
                     color=GROUP_A_COLOR, alpha=0.08, zorder=0)
    ax.text(0.002, (y0 + y1) / 2, 'Group A',
            fontsize=8, color=GROUP_A_COLOR, fontweight='bold',
            va='center', ha='left', zorder=5)

if group_b_ys:
    y0, y1 = min(group_b_ys) - 0.45, max(group_b_ys) + 0.45
    ax.fill_betweenx([y0, y1], 0, bracket_x,
                     color=GROUP_B_COLOR, alpha=0.16, zorder=0)
    ax.text(0.002, (y0 + y1) / 2, 'Group B',
            fontsize=8, color=GROUP_B_COLOR, fontweight='bold',
            va='center', ha='left', zorder=5)

# ---------------------------------------------------------------------------
# Axis formatting
# ---------------------------------------------------------------------------
ax.set_xlabel('Branch length (substitutions per site)', fontsize=10)
ax.set_ylabel('')
ax.set_yticks([])
ax.set_title(
    'Sphyrnidae COI maximum-likelihood phylogeny\n'
    '(FastTree, GTR model; local support values at nodes)',
    fontsize=11, fontweight='bold'
)

# ---------------------------------------------------------------------------
# Legend — species present in tree, plus S. lewini group entries
# ---------------------------------------------------------------------------
species_in_tree = []
seen = set()
for orig in orig_names:
    sp = parse_species(orig)
    if sp not in seen:
        seen.add(sp)
        species_in_tree.append(sp)

legend_patches = []
for sp in species_in_tree:
    g, e   = sp.split(' ', 1)
    color  = SPECIES_COLORS.get(sp, FALLBACK_COLOR)
    # S. lewini gets two entries (Group A / B) — skip the generic entry
    if sp == 'Sphyrna lewini':
        continue
    legend_patches.append(
        mpatches.Patch(facecolor=color, label=rf'$\it{{{g}\ {e}}}$')
    )

# Add S. lewini group entries at top
legend_patches = [
    mpatches.Patch(facecolor=GROUP_A_COLOR,
                   label=rf'$\it{{S.\ lewini}}$ Group A (Indian Ocean / Atlantic)'),
    mpatches.Patch(facecolor=GROUP_B_COLOR,
                   label=rf'$\it{{S.\ lewini}}$ Group B (W. Pacific)'),
] + legend_patches

# Place legend outside the axes to the right so it never covers tip labels
ax.legend(handles=legend_patches,
          loc='upper left',
          bbox_to_anchor=(1.02, 1.0),
          bbox_transform=ax.transAxes,
          fontsize=7.5, title='Species', title_fontsize=8.5,
          framealpha=0.92, borderpad=0.8)

plt.tight_layout()

out_png = Path('scripts/Figure3.png')
out_pdf = Path('scripts/Figure3.pdf')
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
