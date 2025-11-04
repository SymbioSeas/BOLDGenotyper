"""
Publication-Ready Visualization Generation

This module provides functions for creating high-quality figures to visualize
genetic partitioning patterns and biogeographic distributions. All figures are
designed to be publication-ready with appropriate resolution, styling, and
formatting.

Figure Types:
1. Global Distribution Map
   - Points colored by genotype
   - Point size proportional to sample count at location
   - World map background with country boundaries
   - Legend with genotype colors
   - Scale bar and coordinates

2. Ocean Basin Abundance Plot
   - Stacked bar chart showing genotype proportions per basin
   - Bars arranged by ocean basin
   - Colors consistent with distribution map
   - Y-axis shows relative abundance (0-100%)
   - Sample counts annotated

3. Phylogenetic Tree (optional)
   - Maximum likelihood tree with bootstrap support
   - Branch tips colored by genotype
   - Bootstrap values shown for key nodes
   - Scale bar for genetic distance

4. Sequence Alignment Viewer (optional)
   - Pairwise alignment highlighting differences
   - Color-coded nucleotides
   - Consensus sequence shown

Design Specifications:
- Color palette: Colorblind-friendly (seaborn 'colorblind' or 'tab10')
  For consistency with reference analysis, first 3 colors approximate:
  * Purple (#9D7ABE)
  * Teal (#5AB4AC)
  * Yellow (#F2CC8F)
- Output formats: PNG (300 DPI) and PDF (vector)
- Figure size: Publication-ready (e.g., 10x6 inches for maps)
- Font: Arial or Helvetica, 10-12pt
- File naming: {organism}_distribution_map.png, {organism}_ocean_basins.pdf

The visualization module handles color palette generation dynamically based on
the number of genotypes detected (scales gracefully to 20+ genotypes).

Example Usage:
    >>> from boldgenotyper.visualization import plot_distribution_map, plot_ocean_basin_abundance
    >>> plot_distribution_map(
    ...     df=metadata_df,
    ...     output_path="Sphyrna_lewini_distribution_map.png",
    ...     genotype_column="consensus_group"
    ... )
    >>> plot_ocean_basin_abundance(
    ...     df=metadata_df,
    ...     output_path="Sphyrna_lewini_ocean_basins.pdf"
    ... )

Author: Steph Smith (steph.smith@unc.edu)
"""

from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Configure logging
logger = logging.getLogger(__name__)

# Define reference color palette
REFERENCE_COLORS = ['#9D7ABE', '#5AB4AC', '#F2CC8F']


def get_genotype_colors(n_genotypes: int) -> List[str]:
    """
    Generate colorblind-friendly color palette for genotypes.

    Uses reference colors for first 3 genotypes, then extends with
    seaborn colorblind palette.

    Parameters
    ----------
    n_genotypes : int
        Number of genotypes to assign colors

    Returns
    -------
    List[str]
        List of hex color codes
    """
    if n_genotypes <= 0:
        return []
    base = REFERENCE_COLORS[:]
    if n_genotypes <= len(base):
        return base[:n_genotypes]
    try:
        pal = sns.color_palette("colorblind", max(n_genotypes, 10)).as_hex()
    except Exception:
        # fallback to matplotlib tab colors
        pal = list(mcolors.TABLEAU_COLORS.values())
        while len(pal) < n_genotypes:
            pal += pal  # repeat if needed
    # ensure reference colors lead
    remainder = [c for c in pal if c not in base]
    out = (base + remainder)[:n_genotypes]
    return out
    
    pass


def plot_distribution_map(
    df: pd.DataFrame,
    output_path: str,
    genotype_column: str = "consensus_group",
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300,
) -> None:
    """
    Create global distribution map with points colored by genotype.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata with coordinates and genotype assignments
    output_path : str
        Path for output figure (PNG or PDF)
    genotype_column : str, optional
        Column containing genotype assignments
    latitude_col : str, optional
        Column containing latitude values
    longitude_col : str, optional
        Column containing longitude values
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 10x6)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # drop rows without coords or genotype
    need = [genotype_column, latitude_col, longitude_col]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found in dataframe")
    d = df.dropna(subset=[genotype_column, latitude_col, longitude_col]).copy()
    if d.empty:
        logger.warning("No rows with valid coordinates; skipping map.")
        return

    genos = sorted(d[genotype_column].dropna().unique())
    colors = get_genotype_colors(len(genos))
    color_map = {g: colors[i] for i, g in enumerate(genos)}
    d["_color"] = d[genotype_column].map(color_map)

    use_cartopy = True
    try:
        _ = ccrs.PlateCarree()
    except Exception:
        use_cartopy = False
        logger.warning("Cartopy unavailable; drawing scatter without basemap.")

    plt.figure(figsize=figsize)
    if use_cartopy:
        ax = plt.axes(projection=ccrs.Robinson())
        ax.add_feature(cfeature.LAND, zorder=0, edgecolor="black", linewidth=0.2, facecolor="#f2f2f2")
        ax.add_feature(cfeature.OCEAN, zorder=0, facecolor="#d9edf7")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.3)
        ax.gridlines(draw_labels=False, linewidth=0.2, color="gray", alpha=0.5)
        for g in genos:
            sub = d[d[genotype_column] == g]
            ax.scatter(
                sub[longitude_col], sub[latitude_col],
                transform=ccrs.PlateCarree(), s=20, alpha=0.8, label=str(g),
                color=color_map[g], edgecolors="black", linewidths=0.2,
            )
    else:
        ax = plt.gca()
        for g in genos:
            sub = d[d[genotype_column] == g]
            ax.scatter(
                sub[longitude_col], sub[latitude_col],
                s=20, alpha=0.8, label=str(g),
                color=color_map[g], edgecolors="black", linewidths=0.2,
            )
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
        ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
        ax.grid(True, linestyle="--", linewidth=0.3, alpha=0.5)

    ax.legend(title="Genotype", loc="lower left", bbox_to_anchor=(1.02, 0.0), frameon=False)
    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()
    
    pass


def plot_ocean_basin_abundance(
    df: pd.DataFrame,
    output_path: str,
    genotype_column: str = "consensus_group",
    basin_column: str = "ocean_basin",
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300,
) -> None:
    """
    Create stacked bar chart showing genotype proportions by ocean basin.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata with genotype and ocean basin assignments
    output_path : str
        Path for output figure (PNG or PDF)
    genotype_column : str, optional
        Column containing genotype assignments
    basin_column : str, optional
        Column containing ocean basin names
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 10x6)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    need = [genotype_column, basin_column]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found in dataframe")
    d = df.copy()
    d[basin_column] = d[basin_column].fillna("Unknown")
    d[genotype_column] = d[genotype_column].fillna("Unassigned")

    # proportions by basin
    counts = (
        d.groupby([basin_column, genotype_column]).size().rename("n").reset_index()
    )
    totals = counts.groupby(basin_column)["n"].transform("sum")
    counts["prop"] = counts["n"] / totals

    basins = counts[basin_column].unique().tolist()
    genos  = counts[genotype_column].unique().tolist()
    colors = get_genotype_colors(len(genos))
    color_map = {g: colors[i] for i, g in enumerate(genos)}
    label_map = {}
    if "consensus_group" in df.columns and "consensus_group_sp" in df.columns:
        tmp = df[["consensus_group", "consensus_group_sp"]].dropna().drop_duplicates()
        label_map = dict(zip(tmp["consensus_group"], tmp["consensus_group_sp"]))

    # pivot to stacked proportions
    wide = counts.pivot(index=basin_column, columns=genotype_column, values="prop").fillna(0.0)
    wide = wide.reindex(index=sorted(wide.index))  # sort basins
    wide = wide[[g for g in genos if g in wide.columns]]

    plt.figure(figsize=figsize)
    ax = plt.gca()
    bottom = None
    for g in wide.columns:
        vals = wide[g].values
        ax.bar(wide.index, vals, bottom=bottom, label=str(g), color=color_map[g])
        bottom = vals if bottom is None else bottom + vals

    ax.set_ylabel("Relative abundance")
    ax.set_xlabel("Ocean basin")
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f"{int(v*100)}%"))
    ax.legend(title="Genotype", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()

    pass


def plot_phylogenetic_tree(
    tree_file: str,
    output_path: str,
    genotype_colors: Optional[Dict[str, str]] = None,
    show_bootstrap: bool = True,
    bootstrap_threshold: int = 70,
    figsize: Tuple[int, int] = (8, 10),
    dpi: int = 300,
    label_map: Optional[Dict[str, str]] = None,
) -> None:
    """
    Visualize phylogenetic tree with genotype colors and bootstrap support.

    Parameters
    ----------
    tree_file : str
        Path to Newick format tree file
    output_path : str
        Path for output figure
    genotype_colors : Dict[str, str], optional
        Mapping of genotype names to colors
    show_bootstrap : bool, optional
        Display bootstrap values (default: True)
    bootstrap_threshold : int, optional
        Minimum bootstrap value to display (default: 70)
    figsize : Tuple[int, int], optional
        Figure size in inches
    dpi : int, optional
        Resolution for PNG output
    """
    from Bio import Phylo

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    tree = Phylo.read(str(tree_file), "newick")
    
    # Apply label_map
    if label_map:
        for clade in tree.get_terminals():
            if clade.name in label_map:
                clade.name = label_map[clade.name]

    # Map tip colors
    tip_colors = {}
    if genotype_colors:
        tip_colors = {}
        for clade in tree.get_terminals():
            g = clade.name
            if g in genotype_colors:
                tip_colors[clade] = genotype_colors[g]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(
        tree, do_show=False, axes=ax, label_colors=tip_colors if tip_colors else None
    )

    # Optionally annotate bootstrap values (internal node confidences)
    if show_bootstrap:
        for clade in tree.get_nonterminals():
            if clade.confidence is not None:
                val = float(clade.confidence)
                if val >= bootstrap_threshold:
                    x = clade.branch_length if clade.branch_length else 0.0
                    # Bio.Phylo's draw already places text; skipping extra labels to avoid clutter

    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()
    
    pass
