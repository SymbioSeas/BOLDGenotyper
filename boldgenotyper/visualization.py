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
    # Implementation will go here
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
    # Implementation will go here
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
    # Implementation will go here
    pass


def plot_phylogenetic_tree(
    tree_file: str,
    output_path: str,
    genotype_colors: Optional[Dict[str, str]] = None,
    show_bootstrap: bool = True,
    bootstrap_threshold: int = 70,
    figsize: Tuple[int, int] = (8, 10),
    dpi: int = 300,
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
    # Implementation will go here
    pass
