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

    Point sizes are scaled by the number of samples at each location.

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

    # Count samples with coordinates and genotypes separately for diagnostic message
    total_samples = len(df)
    samples_with_coords = df[[latitude_col, longitude_col]].notna().all(axis=1).sum()
    samples_with_genotype = df[genotype_column].notna().sum()

    d = df.dropna(subset=[genotype_column, latitude_col, longitude_col]).copy()
    samples_with_both = len(d)

    if d.empty or samples_with_both < 5:
        msg = (
            f"Insufficient data for distribution map (need at least 5 samples with both coordinates and genotypes). "
            f"Total samples: {total_samples}, "
            f"with coordinates: {samples_with_coords} ({samples_with_coords/total_samples*100:.1f}%), "
            f"with genotypes: {samples_with_genotype} ({samples_with_genotype/total_samples*100:.1f}%), "
            f"with BOTH: {samples_with_both} ({samples_with_both/total_samples*100:.1f}%). "
            f"Skipping map generation for: {output_path}"
        )
        logger.warning(msg)
        print(f"⚠ WARNING: {msg}")
        return

    # Count samples at each location for sizing
    d['_count'] = d.groupby([latitude_col, longitude_col, genotype_column])[genotype_column].transform('count')
    # Scale point sizes: min 20, max 200, proportional to count
    min_size, max_size = 20, 200
    if d['_count'].max() > 1:
        d['_size'] = min_size + (d['_count'] - d['_count'].min()) / (d['_count'].max() - d['_count'].min()) * (max_size - min_size)
    else:
        d['_size'] = 50  # default size if all counts are 1

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
        ax = plt.axes(projection=ccrs.PlateCarree())
        # Suppress cartopy download warnings (these are informational only)
        import warnings
        from cartopy.io import DownloadWarning
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DownloadWarning)
            ax.add_feature(cfeature.LAND, zorder=0, edgecolor="black", linewidth=0.2, facecolor="#f2f2f2")
            ax.add_feature(cfeature.OCEAN, zorder=0, facecolor="#d9edf7")
            ax.add_feature(cfeature.COASTLINE, linewidth=0.3)
        # Add gridlines with labels (PlateCarree projection supports proper gridline labels)
        gl = ax.gridlines(draw_labels=True, linewidth=0.2, color="gray", alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 10}
        for g in genos:
            sub = d[d[genotype_column] == g]
            ax.scatter(
                sub[longitude_col], sub[latitude_col],
                transform=ccrs.PlateCarree(), s=sub['_size'], alpha=0.8, label=str(g),
                color=color_map[g], edgecolors="black", linewidths=0.2,
            )
    else:
        ax = plt.gca()
        for g in genos:
            sub = d[d[genotype_column] == g]
            ax.scatter(
                sub[longitude_col], sub[latitude_col],
                s=sub['_size'], alpha=0.8, label=str(g),
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


def plot_distribution_map_faceted(
    df: pd.DataFrame,
    output_path: str,
    genotype_column: str = "consensus_group",
    species_column: str = "assigned_sp",
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
    width: int = 10,
    height_per_species: int = 5,
    dpi: int = 300,
) -> None:
    """
    Create distribution map with separate facets for each species.

    Each facet shows genotypes for a single species, stacked vertically
    in a single column. Species labels are italicized. Point sizes are
    scaled by the number of samples at each location.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata with coordinates, genotype, and species assignments
    output_path : str
        Path for output figure (PNG or PDF)
    genotype_column : str, optional
        Column containing genotype assignments
    species_column : str, optional
        Column containing species assignments for faceting
    latitude_col : str, optional
        Column containing latitude values
    longitude_col : str, optional
        Column containing longitude values
    width : int, optional
        Figure width in inches (default: 10)
    height_per_species : int, optional
        Height per species facet in inches (default: 5)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Validate required columns
    need = [genotype_column, species_column, latitude_col, longitude_col]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found in dataframe")

    # Count samples with different attributes for diagnostic message
    total_samples = len(df)
    samples_with_coords = df[[latitude_col, longitude_col]].notna().all(axis=1).sum()
    samples_with_genotype = df[genotype_column].notna().sum()
    samples_with_species = df[species_column].notna().sum()

    # Filter to valid rows
    d = df.dropna(subset=[genotype_column, species_column, latitude_col, longitude_col]).copy()
    samples_with_all = len(d)

    if d.empty or samples_with_all < 5:
        msg = (
            f"Insufficient data for faceted distribution map (need at least 5 samples with coordinates, genotypes, and species). "
            f"Total samples: {total_samples}, "
            f"with coordinates: {samples_with_coords} ({samples_with_coords/total_samples*100:.1f}%), "
            f"with genotypes: {samples_with_genotype} ({samples_with_genotype/total_samples*100:.1f}%), "
            f"with species: {samples_with_species} ({samples_with_species/total_samples*100:.1f}%), "
            f"with ALL: {samples_with_all} ({samples_with_all/total_samples*100:.1f}%). "
            f"Skipping faceted map generation for: {output_path}"
        )
        logger.warning(msg)
        print(f"⚠ WARNING: {msg}")
        return

    # Count samples at each location for sizing
    d['_count'] = d.groupby([latitude_col, longitude_col, genotype_column, species_column])[genotype_column].transform('count')

    # Get unique species (sorted)
    species_list = sorted(d[species_column].unique())
    n_species = len(species_list)

    if n_species == 0:
        logger.warning("No species found; skipping faceted map.")
        return

    # Get all genotypes and assign consistent colors
    all_genotypes = sorted(d[genotype_column].unique())
    all_colors = get_genotype_colors(len(all_genotypes))
    color_map = {g: all_colors[i] for i, g in enumerate(all_genotypes)}

    # Check if cartopy is available
    use_cartopy = True
    try:
        _ = ccrs.PlateCarree()
    except Exception:
        use_cartopy = False
        logger.warning("Cartopy unavailable; drawing faceted maps without basemap.")

    # Create figure with cartopy projection for each subplot
    logger.info(f"Creating faceted distribution map with {'cartopy' if use_cartopy else 'simple'} backgrounds")
    fig_height = height_per_species * n_species
    fig = plt.figure(figsize=(width, fig_height))

    axes = []
    for i in range(n_species):
        if use_cartopy:
            ax = fig.add_subplot(n_species, 1, i + 1, projection=ccrs.PlateCarree())
        else:
            ax = fig.add_subplot(n_species, 1, i + 1)
        axes.append(ax)

    for idx, species in enumerate(species_list):
        ax = axes[idx]
        species_data = d[d[species_column] == species].copy()
        species_genotypes = sorted(species_data[genotype_column].unique())

        # Scale point sizes for this species
        min_size, max_size = 20, 200
        if species_data['_count'].max() > 1:
            species_data['_size'] = min_size + (species_data['_count'] - species_data['_count'].min()) / (species_data['_count'].max() - species_data['_count'].min()) * (max_size - min_size)
        else:
            species_data['_size'] = 50

        if use_cartopy:
            # Add cartopy map features
            import warnings
            from cartopy.io import DownloadWarning
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=DownloadWarning)
                ax.add_feature(cfeature.LAND, zorder=0, edgecolor="black", linewidth=0.2, facecolor="#f2f2f2")
                ax.add_feature(cfeature.OCEAN, zorder=0, facecolor="#d9edf7")
                ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

            # Add gridlines with labels
            gl = ax.gridlines(draw_labels=True, linewidth=0.2, color="gray", alpha=0.5, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 10}
            gl.ylabel_style = {'size': 10}

            # Plot scatter points for each genotype with transform
            for g in species_genotypes:
                sub = species_data[species_data[genotype_column] == g]
                if len(sub) > 0:
                    ax.scatter(
                        sub[longitude_col], sub[latitude_col],
                        transform=ccrs.PlateCarree(),
                        s=sub['_size'], alpha=0.8,
                        color=color_map[g], edgecolors="black", linewidths=0.2,
                    )
        else:
            # Fallback to simple scatter without cartopy
            for g in species_genotypes:
                sub = species_data[species_data[genotype_column] == g]
                if len(sub) > 0:
                    ax.scatter(
                        sub[longitude_col], sub[latitude_col],
                        s=sub['_size'], alpha=0.8,
                        color=color_map[g], edgecolors="black", linewidths=0.2,
                    )
            ax.set_xlabel("Longitude", fontsize=10)
            ax.set_ylabel("Latitude", fontsize=10)
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            ax.grid(True, linestyle="--", linewidth=0.3, alpha=0.5)

        # Add species name as title (genotype info is in the species name itself)
        ax.set_title(species, fontsize=11, fontweight='bold', fontstyle='italic', loc='left', pad=10)

    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved faceted distribution map: {out}")


def plot_ocean_basin_abundance_faceted(
    df: pd.DataFrame,
    output_path: str,
    genotype_column: str = "consensus_group",
    species_column: str = "assigned_sp",
    basin_column: str = "ocean_basin",
    width: int = 9,
    height_per_species: int = 5,
    dpi: int = 300,
) -> None:
    """
    Create bar chart with separate facets for each genotype.

    Each facet shows sample counts by ocean basin for a single genotype,
    with n-values annotated above each bar. Genotype labels use consensus_group_sp
    when available for clearer species identification.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata with genotype, species, and ocean basin assignments
    output_path : str
        Path for output figure (PNG or PDF)
    genotype_column : str, optional
        Column containing genotype assignments
    species_column : str, optional
        Column containing species assignments for faceting
    basin_column : str, optional
        Column containing ocean basin names
    width : int, optional
        Figure width in inches (default: 9)
    height_per_species : int, optional
        Height per species facet in inches (default: 5)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Validate required columns
    need = [genotype_column, species_column, basin_column]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found in dataframe")

    # Prepare data
    d = df.copy()
    d[basin_column] = d[basin_column].fillna("Unknown")
    d[genotype_column] = d[genotype_column].fillna("Unassigned")

    # Filter to valid species
    d = d.dropna(subset=[species_column])
    if d.empty:
        logger.warning("No rows with valid species; skipping faceted basin plot.")
        return

    # Create label map for genotypes (consensus_group -> consensus_group_sp)
    label_map = {}
    if "consensus_group_sp" in d.columns:
        tmp = d[[genotype_column, "consensus_group_sp"]].dropna().drop_duplicates()
        label_map = dict(zip(tmp[genotype_column], tmp["consensus_group_sp"]))

    # Get unique genotypes (sorted) - these will be our facets
    genotype_list = sorted(d[genotype_column].unique())
    n_genotypes = len(genotype_list)

    if n_genotypes == 0:
        logger.warning("No genotypes found; skipping faceted basin plot.")
        return

    # Get all genotypes and assign consistent colors
    all_colors = get_genotype_colors(n_genotypes)
    color_map = {g: all_colors[i] for i, g in enumerate(genotype_list)}

    # Get ALL ocean basins across all genotypes (so we show all basins even if 0)
    all_basins = sorted(d[basin_column].unique())

    # Create figure with facets (one per genotype)
    fig_height = height_per_species * n_genotypes
    fig, axes = plt.subplots(n_genotypes, 1, figsize=(width, fig_height))

    # Handle single genotype case
    if n_genotypes == 1:
        axes = [axes]

    for idx, genotype in enumerate(genotype_list):
        ax = axes[idx]
        genotype_data = d[d[genotype_column] == genotype]

        # Count samples by basin for this genotype
        counts = (
            genotype_data.groupby(basin_column)
            .size().reset_index(name="n")
        )

        # Create series indexed by basin for easier plotting
        count_series = counts.set_index(basin_column)["n"]
        # Reindex to include ALL basins (even if this genotype has 0 in some basins)
        count_series = count_series.reindex(index=all_basins, fill_value=0)

        # Plot bars
        bars = ax.bar(count_series.index, count_series.values, color=color_map[genotype])

        # Add n-value annotations above each bar
        for i, (basin, count) in enumerate(count_series.items()):
            if count > 0:  # Only annotate bars with data
                ax.text(i, count, f"n={int(count)}",
                       ha='center', va='bottom', fontsize=9)

        ax.set_ylabel("Sample count")
        ax.set_xlabel("Ocean basin")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha="right")

        # Set y-axis to start at 0 and add a bit of padding for annotations
        max_count = count_series.max()
        ax.set_ylim(0, max_count * 1.15 if max_count > 0 else 1)

    # Adjust margins: more space on left for genotype labels
    plt.subplots_adjust(left=0.15, right=0.95, hspace=0.85)

    # Now add genotype labels positioned based on actual axes positions
    for idx, (genotype, ax) in enumerate(zip(genotype_list, axes)):
        # Get axes position in figure coordinates
        pos = ax.get_position()
        y_center = (pos.y0 + pos.y1) / 2  # center of this axes in figure coords

        # Use consensus_group_sp label if available, otherwise use genotype ID
        display_label = label_map.get(genotype, str(genotype))

        # Add genotype label (outside left margin, centered on facet y-axis)
        # Use italic style for species names (consensus_group_sp format)
        fig.text(0.02, y_center, display_label,
                fontsize=12, fontweight='bold', fontstyle='italic',
                verticalalignment='center', rotation=90, transform=fig.transFigure)

    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved faceted basin abundance plot: {out}")


def plot_phylogenetic_tree(
    tree_file: str,
    output_path: str,
    genotype_colors: Optional[Dict[str, str]] = None,
    show_bootstrap: bool = True,
    bootstrap_threshold: int = 70,
    figsize: Optional[Tuple[int, int]] = None,
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
        Figure size in inches. If None, automatically scales based on number of tips.
        Minimum: (8, 10), scales as: height = max(10, n_tips * 0.3)
    dpi : int, optional
        Resolution for PNG output
    label_map : Dict[str, str], optional
        Mapping to rename tip labels (e.g., consensus_c1 -> species name)
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

    # Calculate figure size based on number of tips if not specified
    if figsize is None:
        n_tips = len(list(tree.get_terminals()))
        # Scale height with number of tips: 0.3 inches per tip, minimum 10, maximum 50
        height = max(10, min(50, n_tips * 0.3))
        # Width scales more slowly: base 8 + extra for large trees
        width = 8 if n_tips <= 30 else min(14, 8 + (n_tips - 30) * 0.1)
        figsize = (width, height)
        logger.info(f"Auto-scaled tree figure size to {figsize} for {n_tips} tips")

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


def plot_identity_distribution(
    diagnostics_csv: str,
    output_path: str,
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300,
) -> None:
    """
    Plot distribution of identity scores for assigned samples.

    Shows histogram with density overlay to visualize the quality and
    confidence of genotype assignments.

    Parameters
    ----------
    diagnostics_csv : str
        Path to diagnostics CSV with identity scores
    output_path : str
        Path for output figure (PNG or PDF)
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 10x6)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load diagnostics
    df = pd.read_csv(diagnostics_csv)

    # Filter to samples with identity scores (exclude no_sequence)
    df = df[df['identity'] > 0].copy()

    if df.empty:
        logger.warning("No samples with identity scores; skipping identity distribution plot.")
        return

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram with density
    ax.hist(df['identity'], bins=50, alpha=0.6, color='#5AB4AC', edgecolor='black', linewidth=0.5, density=True, label='Distribution')

    # Add density curve
    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(df['identity'])
        x_range = np.linspace(df['identity'].min(), df['identity'].max(), 200)
        ax.plot(x_range, kde(x_range), color='#9D7ABE', linewidth=2, label='Density')
    except:
        pass  # Skip KDE if it fails

    # Add vertical lines for mean and median
    mean_val = df['identity'].mean()
    median_val = df['identity'].median()
    ax.axvline(mean_val, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color='orange', linestyle='--', linewidth=1.5, label=f'Median: {median_val:.3f}')

    # Labels and formatting
    ax.set_xlabel('Identity Score', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Distribution of Genotype Assignment Identity Scores', fontsize=14, fontweight='bold')
    ax.legend(frameon=False, fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved identity distribution plot: {out}")


def plot_identity_by_status(
    diagnostics_csv: str,
    output_path: str,
    figsize: Tuple[int, int] = (12, 6),
    dpi: int = 300,
) -> None:
    """
    Plot identity scores faceted by assignment status.

    Shows box plots comparing identity distributions across different
    assignment statuses (assigned, low_confidence, tie, etc.).

    Parameters
    ----------
    diagnostics_csv : str
        Path to diagnostics CSV with identity and status
    output_path : str
        Path for output figure (PNG or PDF)
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 12x6)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load diagnostics
    df = pd.read_csv(diagnostics_csv)

    # Filter to samples with identity scores
    df = df[df['identity'] > 0].copy()

    if df.empty:
        logger.warning("No samples with identity scores; skipping identity by status plot.")
        return

    # Get status order (by frequency)
    status_order = df['status'].value_counts().index.tolist()

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create box plot
    positions = range(len(status_order))
    bp_data = [df[df['status'] == status]['identity'].values for status in status_order]

    bp = ax.boxplot(bp_data, positions=positions, widths=0.6, patch_artist=True,
                    showmeans=True, meanline=True,
                    boxprops=dict(facecolor='#5AB4AC', alpha=0.7),
                    medianprops=dict(color='red', linewidth=2),
                    meanprops=dict(color='orange', linewidth=2, linestyle='--'),
                    whiskerprops=dict(color='black', linewidth=1),
                    capprops=dict(color='black', linewidth=1))

    # Add sample counts
    for i, status in enumerate(status_order):
        n = len(df[df['status'] == status])
        ax.text(i, 0.4, f'n={n}', ha='center', va='top', fontsize=9, style='italic')

    # Labels and formatting
    ax.set_xticks(positions)
    ax.set_xticklabels(status_order, rotation=20, ha='right')
    ax.set_xlabel('Assignment Status', fontsize=12)
    ax.set_ylabel('Identity Score', fontsize=12)
    ax.set_title('Identity Scores by Assignment Status', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.grid(True, axis='y', linestyle='--', alpha=0.3)

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, label='Median'),
        Line2D([0], [0], color='orange', linewidth=2, linestyle='--', label='Mean')
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=10)

    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved identity by status plot: {out}")


def plot_assignment_status(
    diagnostics_csv: str,
    output_path: str,
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300,
) -> None:
    """
    Plot assignment status breakdown as a stacked bar chart.

    Shows the proportion of samples in each assignment status category.

    Parameters
    ----------
    diagnostics_csv : str
        Path to diagnostics CSV with assignment status
    output_path : str
        Path for output figure (PNG or PDF)
    figsize : Tuple[int, int], optional
        Figure size in inches (default: 10x6)
    dpi : int, optional
        Resolution for PNG output (default: 300)
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load diagnostics
    df = pd.read_csv(diagnostics_csv)

    # Count by status
    status_counts = df['status'].value_counts()
    total = len(df)

    # Define colors for each status
    status_colors = {
        'assigned': '#5AB4AC',
        'low_confidence': '#F2CC8F',
        'tie': '#E07A7A',
        'below_threshold': '#D3D3D3',
        'no_sequence': '#A9A9A9'
    }

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left panel: Pie chart
    colors = [status_colors.get(status, '#CCCCCC') for status in status_counts.index]
    wedges, texts, autotexts = ax1.pie(
        status_counts.values,
        labels=status_counts.index,
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 10}
    )
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    ax1.set_title('Assignment Status Distribution', fontsize=12, fontweight='bold')

    # Right panel: Bar chart with counts
    ax2.bar(range(len(status_counts)), status_counts.values,
            color=[status_colors.get(status, '#CCCCCC') for status in status_counts.index],
            alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.set_xticks(range(len(status_counts)))
    ax2.set_xticklabels(status_counts.index, rotation=20, ha='right')
    ax2.set_ylabel('Number of Samples', fontsize=12)
    ax2.set_xlabel('Assignment Status', fontsize=12)
    ax2.set_title('Sample Counts by Status', fontsize=12, fontweight='bold')
    ax2.grid(True, axis='y', linestyle='--', alpha=0.3)

    # Add count labels on bars
    for i, (status, count) in enumerate(status_counts.items()):
        ax2.text(i, count + total*0.01, f'{count}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    if out.suffix.lower() == ".png":
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        plt.savefig(out, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved assignment status plot: {out}")
