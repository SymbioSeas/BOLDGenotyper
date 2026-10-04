"""
Shared plotting helpers: the haplotype color palette and map drawing.

The run's figures are drawn by tree_figures, run_figures and
metadata_figures (organized by species, with manifests); the query figures
by query_overlays. This module keeps what they share:

- get_genotype_colors: a colorblind-friendly palette, reference colors first
- calculate_map_extent_with_buffer: a map extent around the points
- _add_cartopy_background, _add_map_gridlines: cartopy map backgrounds
- _point_sizes: marker sizes scaled by records per site

Author: Steph Smith (symbioseas@outlook.com)
"""

from typing import List, Optional, Tuple
import logging
import pandas as pd
import matplotlib.colors as mcolors
import seaborn as sns
import cartopy.feature as cfeature

# Configure logging
logger = logging.getLogger(__name__)

# Define reference color palette
REFERENCE_COLORS = ['#8545C1', '#10B3A5', '#FFB031', '#3874F4', '#D975C7', "#132C54"]


def calculate_map_extent_with_buffer(
    latitudes: pd.Series,
    longitudes: pd.Series,
    buffer_degrees: float = 20.0
) -> Tuple[float, float, float, float]:
    """
    Calculate map extent with buffer, constrained to valid lat/lon ranges.

    Adds buffer around data extent while ensuring bounds don't exceed
    valid geographic coordinates.

    Parameters
    ----------
    latitudes : pd.Series
        Latitude values
    longitudes : pd.Series
        Longitude values
    buffer_degrees : float, optional
        Buffer margin in degrees (default: 20.0)

    Returns
    -------
    lon_min, lon_max, lat_min, lat_max : floats
        Map extent constrained to [-180, 180] × [-90, 90]

    Examples
    --------
    >>> lats = pd.Series([10, 20, 30])
    >>> lons = pd.Series([-100, -90, -80])
    >>> calculate_map_extent_with_buffer(lats, lons, buffer_degrees=10)
    (-110.0, -70.0, 0.0, 40.0)
    """
    lat_min, lat_max = latitudes.min(), latitudes.max()
    lon_min, lon_max = longitudes.min(), longitudes.max()

    # Add buffer
    lat_min = max(-90, lat_min - buffer_degrees)
    lat_max = min(90, lat_max + buffer_degrees)
    lon_min = max(-180, lon_min - buffer_degrees)
    lon_max = min(180, lon_max + buffer_degrees)

    return lon_min, lon_max, lat_min, lat_max


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
        pal = sns.color_palette("husl", max(n_genotypes, 10)).as_hex()
    except Exception:
        # fallback to matplotlib tab colors
        pal = list(mcolors.TABLEAU_COLORS.values())
        while len(pal) < n_genotypes:
            pal += pal  # repeat if needed
    # ensure reference colors lead
    remainder = [c for c in pal if c not in base]
    out = (base + remainder)[:n_genotypes]
    return out


def _add_map_gridlines(ax):
    """Draw coordinate labels on a map without visible grid lines.

    Width-0 lines hide the grid while keeping labels. The line style must be
    solid: matplotlib scales dash patterns by line width, and a dashed style at
    width 0 fails with "At least one value in the dash list must be positive".
    """
    return ax.gridlines(draw_labels=True, linewidth=0, color="gray", alpha=0.5, linestyle="-")


def _add_cartopy_background(ax, detail_level: str = "full"):
    """
    Add cartopy map background features with configurable detail level.

    This function provides performance optimization for cartopy-based maps by
    allowing simplified backgrounds that render faster and produce smaller files.

    Parameters
    ----------
    ax : cartopy.mpl.geoaxes.GeoAxes
        Cartopy axes object with PlateCarree projection
    detail_level : str
        Level of detail for map background:
        - "full": Detailed land/ocean/coastline polygons (publication quality)
        - "simple": Coastlines only with solid colors (70-90% faster rendering)

    Notes
    -----
    File size impact for SVG output:
    - "full": ~10-15 MB per map (detailed vector coastlines)
    - "simple": ~500 KB per map (simplified geometry)

    Rendering time impact:
    - "full": Baseline
    - "simple": 70-90% faster for complex coastlines

    Examples
    --------
    >>> ax = plt.axes(projection=ccrs.PlateCarree())
    >>> _add_cartopy_background(ax, detail_level="simple")
    """
    import warnings
    from cartopy.io import DownloadWarning

    # Suppress cartopy download warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DownloadWarning)

        if detail_level == "full":
            # Full detail: render all features with detailed polygons
            ax.add_feature(cfeature.LAND, zorder=0, edgecolor="black",
                          linewidth=0.2, facecolor="#f2f2f2")
            ax.add_feature(cfeature.OCEAN, zorder=0, facecolor="#d9edf7")
            ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

        elif detail_level == "simple":
            # Simplified: coastlines only with solid background colors
            # Use low-resolution coastlines for faster rendering
            ax.coastlines(resolution='110m', linewidth=0.5, color='black')
            # Set solid background colors (no complex polygons)
            ax.set_facecolor('#d9edf7')  # Ocean blue
            # Add land as patches (much faster than cfeature.LAND)
            ax.add_feature(cfeature.LAND.with_scale('110m'),
                          zorder=0, facecolor="#f2f2f2", edgecolor='none')

        else:
            raise ValueError(f"Unknown detail_level: {detail_level}. Use 'full' or 'simple'.")


def _point_sizes(counts: pd.Series, min_size: float, max_size: float,
                 equal_size: Optional[float] = None) -> pd.Series:
    """
    Marker sizes scaled linearly from min_size (smallest count) to max_size
    (largest). When every count is the same there is nothing to scale: all
    markers get equal_size (default min_size). Scaling then divided by zero,
    and the NaN sizes broke the legend when the figure was saved.
    """
    lo, hi = counts.min(), counts.max()
    if hi > lo:
        return min_size + (counts - lo) / (hi - lo) * (max_size - min_size)
    return pd.Series(min_size if equal_size is None else equal_size, index=counts.index, dtype=float)


# ============================================================================
# Species-Level Visualization Functions
# ============================================================================


# ============================================================================
# Species-Faceted Haplotype Visualization Functions
# ============================================================================


