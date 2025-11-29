"""
Plot Regeneration Module for BOLDGenotyper

Python-based plot regeneration from exported data with custom styling.
This module provides functions to regenerate all plots using matplotlib/seaborn/geopandas,
reading configuration from plot_config.yaml.

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import logging
import warnings

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Suppress expected warnings from shapely and numpy when handling missing/invalid coordinates
# These warnings are expected when data contains NaN values and don't affect plot generation
warnings.filterwarnings('ignore', category=RuntimeWarning, module='shapely')
warnings.filterwarnings('ignore', message='invalid value encountered')
warnings.filterwarnings('ignore', message='Tight layout not applied')
np.seterr(invalid='ignore')

# Optional imports with fallbacks
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

try:
    import geopandas as gpd
    from shapely.geometry import Point
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    CARTOPY_AVAILABLE = True
except ImportError:
    CARTOPY_AVAILABLE = False

logger = logging.getLogger(__name__)


def convert_r_color_to_matplotlib(color: str) -> str:
    """
    Convert R-style color names to matplotlib-compatible colors.

    Parameters
    ----------
    color : str
        Color name or hex code

    Returns
    -------
    str
        Matplotlib-compatible color
    """
    # Handle R's grayXX naming (gray0 to gray100)
    if color.startswith('gray') or color.startswith('grey'):
        try:
            # Extract number (e.g., "gray70" -> 70)
            num = int(color.replace('gray', '').replace('grey', ''))
            # Convert to 0-1 scale (e.g., 70 -> 0.7)
            return str(num / 100.0)
        except (ValueError, AttributeError):
            pass

    # Common R color name conversions
    r_to_mpl = {
        'lightgray': 'lightgray',
        'darkgray': 'darkgray',
        'gray': 'gray',
    }

    if color in r_to_mpl:
        return r_to_mpl[color]

    # Return as-is if already valid (hex codes, named colors)
    return color


def load_config(config_path: Path) -> Dict[str, Any]:
    """
    Load plot configuration from YAML file.

    Parameters
    ----------
    config_path : Path
        Path to plot_config.yaml

    Returns
    -------
    dict
        Configuration dictionary
    """
    if not YAML_AVAILABLE:
        raise ImportError("PyYAML is required for plot regeneration. Install with: pip install pyyaml")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    return config


def apply_genotype_filters(
    df: pd.DataFrame,
    filters: Dict[str, Any],
    genotype_col: str = 'genotype'
) -> pd.DataFrame:
    """
    Apply include/exclude filters to genotypes.

    Parameters
    ----------
    df : pd.DataFrame
        Data frame to filter
    filters : dict
        Filter configuration with 'include_genotypes' and 'exclude_genotypes' keys
    genotype_col : str
        Column name containing genotype labels

    Returns
    -------
    pd.DataFrame
        Filtered dataframe
    """
    result = df.copy()

    include = filters.get('include_genotypes', [])
    exclude = filters.get('exclude_genotypes', [])

    # Handle case where YAML parses single value as string instead of list
    if isinstance(include, str):
        include = [include]
    if isinstance(exclude, str):
        exclude = [exclude]

    # Convert to list if needed
    if include and not isinstance(include, list):
        include = list(include)
    if exclude and not isinstance(exclude, list):
        exclude = list(exclude)

    if include:
        result = result[result[genotype_col].isin(include)]

    if exclude:
        result = result[~result[genotype_col].isin(exclude)]

    return result


def get_palette(genotypes: List[str], color_config: Dict[str, str]) -> Dict[str, str]:
    """
    Build color palette for genotypes, using config colors and fallbacks.

    Parameters
    ----------
    genotypes : list
        List of genotype names
    color_config : dict
        Mapping of genotype names to hex colors

    Returns
    -------
    dict
        Complete color mapping for all genotypes
    """
    palette = {}

    # Use configured colors
    for genotype in genotypes:
        if genotype in color_config:
            palette[genotype] = color_config[genotype]

    # Generate fallback colors for missing genotypes
    missing = [g for g in genotypes if g not in palette]
    if missing:
        # Use tab20 for consistent, distinct colors
        n_missing = len(missing)
        if n_missing <= 20:
            colors = plt.cm.tab20(np.linspace(0, 1, 20))
        else:
            colors = plt.cm.tab20b(np.linspace(0, 1, n_missing))

        for i, genotype in enumerate(missing):
            palette[genotype] = '#%02x%02x%02x' % tuple(int(c * 255) for c in colors[i % len(colors)][:3])

    return palette


def regenerate_distribution_map(
    base_dir: Path,
    output_dir: Optional[Path] = None,
    config: Optional[Dict[str, Any]] = None
) -> List[Path]:
    """
    Regenerate distribution map from exported data.

    Parameters
    ----------
    base_dir : Path
        Base directory containing plot_config.yaml and data/ subdirectory
    output_dir : Path, optional
        Output directory for plots (default: base_dir/visualization)
    config : dict, optional
        Configuration dictionary (if not provided, will load from plot_config.yaml)

    Returns
    -------
    list of Path
        Paths to generated plot files
    """
    if not GEOPANDAS_AVAILABLE:
        raise ImportError("geopandas is required for map regeneration. Install with: conda install geopandas")

    if not CARTOPY_AVAILABLE:
        logger.warning("cartopy not available, using basic projection")

    # Load configuration
    if config is None:
        config_path = base_dir / "plot_config.yaml"
        config = load_config(config_path)

    colors_cfg = config.get('colors', {})
    filters = config.get('filters', {})
    map_cfg = config.get('map', {})
    general_cfg = config.get('general', {})

    # Load data
    data_dir = base_dir / "data"
    map_path = data_dir / "distribution_map.csv"

    if not map_path.exists():
        raise FileNotFoundError(f"distribution_map.csv not found in {data_dir}")

    df = pd.read_csv(map_path)

    if df.empty:
        raise ValueError("No rows in distribution_map.csv")

    # Determine genotype column
    if 'consensus_group_sp' in df.columns:
        df['genotype'] = df['consensus_group_sp'].fillna(df.get('consensus_group', ''))
    else:
        df['genotype'] = df.get('consensus_group', '')

    # Filter out rows with missing genotypes (NaN or empty strings)
    df = df[df['genotype'].notna() & (df['genotype'] != '')]

    # Filter out unknown geography if requested
    if filters.get('exclude_unknown_geography', False):
        if 'ocean_basin' in df.columns:
            unknown_mask = df['ocean_basin'].fillna('Unknown').str.lower().str.contains('unknown')
            df = df[~unknown_mask]
            logger.info(f"Excluded {unknown_mask.sum()} samples with unknown geography")

    # Store available genotypes before filtering for error messages
    available_genotypes = sorted(df['genotype'].unique())
    logger.info(f"Available genotypes in data: {available_genotypes}")
    logger.info(f"Filters - include: {filters.get('include_genotypes', [])}, exclude: {filters.get('exclude_genotypes', [])}")

    # Apply filters
    df = apply_genotype_filters(df, filters, 'genotype')

    if df.empty:
        include = filters.get('include_genotypes', [])
        exclude = filters.get('exclude_genotypes', [])
        error_msg = f"No data to plot after applying genotype filters.\n"
        error_msg += f"  Available genotypes in data: {available_genotypes}\n"
        error_msg += f"  Requested include: {include}\n"
        error_msg += f"  Requested exclude: {exclude}\n"
        error_msg += "  Tip: Check for exact spelling, spaces, and capitalization"
        raise ValueError(error_msg)

    # Build color palette
    genotypes = sorted(df['genotype'].unique())
    palette = get_palette(genotypes, colors_cfg)

    # Filter out rows with invalid coordinates before creating geometries
    # This prevents shapely warnings about NaN values
    valid_coords = df['lat'].notna() & df['lon'].notna() & \
                   (df['lat'] != '') & (df['lon'] != '')
    df = df[valid_coords].copy()

    if df.empty:
        raise ValueError("No valid coordinates found in data after filtering")

    # Create GeoDataFrame with valid coordinates only
    geometry = [Point(xy) for xy in zip(df['lon'], df['lat'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

    # Create figure
    projection = map_cfg.get('projection', 'robinson')
    center_lon = map_cfg.get('center_longitude', 0)

    if CARTOPY_AVAILABLE:
        # Map projection names to cartopy
        proj_map = {
            'robinson': ccrs.Robinson(central_longitude=center_lon),
            'mollweide': ccrs.Mollweide(central_longitude=center_lon),
            'mercator': ccrs.Mercator(),
            'platecarree': ccrs.PlateCarree(central_longitude=center_lon)
        }
        crs = proj_map.get(projection, ccrs.Robinson(central_longitude=center_lon))

        # Use map-specific dimensions if provided, otherwise fall back to general settings
        width = map_cfg.get('width_inches', general_cfg.get('width_inches', 10))
        height = map_cfg.get('height_inches', general_cfg.get('height_inches', 8))

        fig, ax = plt.subplots(figsize=(width, height), subplot_kw={'projection': crs})

        # Add features (convert R color names to matplotlib)
        land_color = convert_r_color_to_matplotlib(map_cfg.get('land_color', '#F5F5F5'))
        ocean_color = convert_r_color_to_matplotlib(map_cfg.get('ocean_color', '#E8F4F8'))
        border_color = convert_r_color_to_matplotlib(map_cfg.get('border_color', 'gray'))

        ax.add_feature(cfeature.LAND, facecolor=land_color)
        ax.add_feature(cfeature.OCEAN, facecolor=ocean_color)

        if map_cfg.get('show_country_borders', True):
            ax.add_feature(cfeature.BORDERS,
                          linewidth=map_cfg.get('border_width', 0.3),
                          edgecolor=border_color)

        ax.coastlines(linewidth=map_cfg.get('border_width', 0.3))

        # Add gridlines with optional labels for lat/lon coordinates
        if map_cfg.get('show_gridlines', True):
            gridline_alpha = map_cfg.get('gridline_alpha', 0.5)
            gridline_style = map_cfg.get('gridline_style', '--')
            show_labels = map_cfg.get('gridline_labels', True)

            gl = ax.gridlines(
                alpha=gridline_alpha,
                linestyle=gridline_style,
                draw_labels=show_labels
            )

            if show_labels:
                gl.top_labels = False
                gl.right_labels = False

        # Transform points to projection
        try:
            gdf_proj = gdf.to_crs(crs.proj4_init if hasattr(crs, 'proj4_init') else crs)
        except:
            gdf_proj = gdf  # Use original CRS if transformation fails

        # Plot points by genotype
        point_alpha = map_cfg.get('point_alpha', 0.7)
        point_stroke = map_cfg.get('point_stroke', 0.5)
        size_range = map_cfg.get('point_size_range', [2, 10])

        # Normalize sizes
        if 'n_at_location' in df.columns:
            sizes = df['n_at_location']
            size_min, size_max = size_range
            sizes_normalized = size_min + (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-6) * (size_max - size_min)
        else:
            sizes_normalized = np.full(len(df), np.mean(size_range))

        for genotype in genotypes:
            mask = gdf_proj['genotype'] == genotype
            if mask.sum() > 0:
                gdf_subset = gdf_proj[mask]
                sizes_subset = sizes_normalized[mask]

                ax.scatter(
                    gdf_subset.geometry.x,
                    gdf_subset.geometry.y,
                    s=sizes_subset**2,  # Square for area scaling
                    c=palette[genotype],
                    alpha=point_alpha,
                    edgecolors='black' if point_stroke > 0 else 'none',
                    linewidths=point_stroke,
                    label=genotype,
                    transform=crs,
                    zorder=10
                )
    else:
        # Fallback: simple matplotlib plot
        # Use map-specific dimensions if provided, otherwise fall back to general settings
        width = map_cfg.get('width_inches', general_cfg.get('width_inches', 10))
        height = map_cfg.get('height_inches', general_cfg.get('height_inches', 8))

        fig, ax = plt.subplots(figsize=(width, height))

        point_alpha = map_cfg.get('point_alpha', 0.7)
        size_range = map_cfg.get('point_size_range', [20, 200])

        if 'n_at_location' in df.columns:
            sizes = df['n_at_location']
            size_min, size_max = size_range
            sizes_normalized = size_min + (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-6) * (size_max - size_min)
        else:
            sizes_normalized = np.full(len(df), np.mean(size_range))

        for genotype in genotypes:
            mask = df['genotype'] == genotype
            if mask.sum() > 0:
                ax.scatter(
                    df.loc[mask, 'lon'],
                    df.loc[mask, 'lat'],
                    s=sizes_normalized[mask],
                    c=palette[genotype],
                    alpha=point_alpha,
                    edgecolors='black',
                    linewidths=0.5,
                    label=genotype
                )

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.grid(True, alpha=0.3)

    # Add legend
    if map_cfg.get('show_legend', True):
        legend_title = map_cfg.get('legend_title', 'Genotype')
        legend_pos = map_cfg.get('legend_position', 'bottom')

        # Configure legend based on position
        legend_kwargs = {
            'title': legend_title,
            'frameon': True,
            'fancybox': True,
            'shadow': True
        }

        if legend_pos == 'right':
            legend_kwargs['loc'] = 'center left'
            legend_kwargs['bbox_to_anchor'] = (1.05, 0.5)
        elif legend_pos == 'left':
            legend_kwargs['loc'] = 'center right'
            legend_kwargs['bbox_to_anchor'] = (-0.05, 0.5)
        elif legend_pos == 'top':
            legend_kwargs['loc'] = 'upper center'
            legend_kwargs['bbox_to_anchor'] = (0.5, 1.15)
            legend_kwargs['ncol'] = min(len(genotypes), 4)
        elif legend_pos == 'bottom':
            legend_kwargs['loc'] = 'upper center'
            legend_kwargs['bbox_to_anchor'] = (0.5, -0.05)
            legend_kwargs['ncol'] = min(len(genotypes), 4)
        else:
            legend_kwargs['loc'] = 'best'

        ax.legend(**legend_kwargs)

    plt.tight_layout()

    # Save outputs
    if output_dir is None:
        output_dir = base_dir / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    formats = general_cfg.get('output_format', ['pdf', 'png'])
    dpi = general_cfg.get('dpi', 300)

    for fmt in formats:
        outfile = output_dir / f"distribution_map_custom.{fmt}"
        plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
        output_files.append(outfile)
        logger.info(f"✓ Wrote {outfile}")

    plt.close()

    return output_files


def regenerate_bar_charts(
    base_dir: Path,
    output_dir: Optional[Path] = None,
    config: Optional[Dict[str, Any]] = None
) -> List[Path]:
    """
    Regenerate bar charts from exported data.

    Parameters
    ----------
    base_dir : Path
        Base directory containing plot_config.yaml and data/ subdirectory
    output_dir : Path, optional
        Output directory for plots (default: base_dir/visualization)
    config : dict, optional
        Configuration dictionary (if not provided, will load from plot_config.yaml)

    Returns
    -------
    list of Path
        Paths to generated plot files
    """
    # Load configuration
    if config is None:
        config_path = base_dir / "plot_config.yaml"
        config = load_config(config_path)

    colors_cfg = config.get('colors', {})
    filters = config.get('filters', {})
    bars_cfg = config.get('bars', {})
    general_cfg = config.get('general', {})

    # Load data
    data_dir = base_dir / "data"
    rel_path = data_dir / "distribution_bar_relative.csv"
    abs_path = data_dir / "distribution_bar_absolute.csv"

    if not rel_path.exists() or not abs_path.exists():
        raise FileNotFoundError(f"Bar chart data files not found in {data_dir}")

    df_rel = pd.read_csv(rel_path)
    df_abs = pd.read_csv(abs_path)

    # Determine genotype column - handle missing consensus_group_sp gracefully
    if 'consensus_group_sp' in df_rel.columns:
        df_rel['genotype'] = df_rel['consensus_group_sp']
    elif 'consensus_group' in df_rel.columns:
        df_rel['genotype'] = df_rel['consensus_group']
    else:
        raise ValueError("Neither 'consensus_group_sp' nor 'consensus_group' found in relative bar data")

    # Handle absolute bar data - create consensus_group_sp if missing
    if 'consensus_group_sp' in df_abs.columns:
        df_abs['genotype'] = df_abs['consensus_group_sp']
    elif 'consensus_group' in df_abs.columns:
        # Try to create consensus_group_sp from consensus_group + species
        if 'species' in df_abs.columns:
            # Remove "consensus_" prefix and add species prefix
            df_abs['genotype'] = df_abs['species'] + ' ' + df_abs['consensus_group'].str.replace('consensus_', '')
        else:
            df_abs['genotype'] = df_abs['consensus_group']
    else:
        raise ValueError("Neither 'consensus_group_sp' nor 'consensus_group' found in absolute bar data")

    # Filter out rows with missing genotypes
    df_rel = df_rel[df_rel['genotype'].notna() & (df_rel['genotype'] != '')]
    df_abs = df_abs[df_abs['genotype'].notna() & (df_abs['genotype'] != '')]

    # Filter out unknown geography if requested
    if filters.get('exclude_unknown_geography', False):
        if 'ocean_basin' in df_rel.columns:
            unknown_mask_rel = df_rel['ocean_basin'].fillna('Unknown').str.lower().str.contains('unknown')
            df_rel = df_rel[~unknown_mask_rel]
        if 'ocean_basin' in df_abs.columns:
            unknown_mask_abs = df_abs['ocean_basin'].fillna('Unknown').str.lower().str.contains('unknown')
            df_abs = df_abs[~unknown_mask_abs]

    # Apply filters
    df_rel = apply_genotype_filters(df_rel, filters, 'genotype')
    df_abs = apply_genotype_filters(df_abs, filters, 'genotype')

    if df_rel.empty or df_abs.empty:
        raise ValueError("No bar data to plot after applying filters")

    # Recalculate relative abundance percentages if requested
    # This ensures filtered genotypes sum to 100% within each ocean basin
    if bars_cfg.get('recalculate_relative', True):
        # Group by ocean basin and recalculate percentages
        if 'pct_of_basin' in df_rel.columns and 'ocean_basin' in df_rel.columns:
            # Calculate total for each basin after filtering
            basin_totals = df_rel.groupby('ocean_basin')['n_samples'].sum()

            # Recalculate percentage for each row
            def recalc_pct(row):
                basin_total = basin_totals.get(row['ocean_basin'], 1)
                if basin_total > 0:
                    return (row['n_samples'] / basin_total) * 100
                return 0

            df_rel['pct_of_basin'] = df_rel.apply(recalc_pct, axis=1)
            logger.info("Recalculated relative abundance percentages for filtered genotypes")

    # Build color palette
    genotypes = sorted(set(df_rel['genotype'].unique()) | set(df_abs['genotype'].unique()))
    palette = get_palette(genotypes, colors_cfg)

    # Plot parameters
    orientation = bars_cfg.get('orientation', 'vertical')
    bar_width = bars_cfg.get('bar_width', 0.8)
    axis_text_angle = bars_cfg.get('axis_text_angle', 45)
    axis_text_size = bars_cfg.get('axis_text_size', 10)
    bar_edgecolor = bars_cfg.get('bar_edgecolor', 'black')
    bar_edgewidth = bars_cfg.get('bar_edgewidth', 0.5)
    grid_behind = bars_cfg.get('grid_behind_bars', True)

    output_files = []

    # Helper function to create bar plot
    def create_bar_plot(df, value_col, value_label, title):
        # Use bar-specific dimensions if provided, otherwise fall back to general settings
        width = bars_cfg.get('width_inches', general_cfg.get('width_inches', 10))
        height = bars_cfg.get('height_inches', general_cfg.get('height_inches', 8))

        fig, ax = plt.subplots(figsize=(width, height))

        # Pivot data for stacking
        # Use sum() to aggregate multiple rows with same genotype in same basin
        pivot_df = df.pivot_table(
            index='ocean_basin',
            columns='genotype',
            values=value_col,
            fill_value=0,
            aggfunc='sum'
        )

        # Reorder columns by genotype order
        pivot_df = pivot_df[[g for g in genotypes if g in pivot_df.columns]]

        # Create stacked bar chart
        if orientation == 'horizontal':
            pivot_df.plot(
                kind='barh',
                stacked=True,
                ax=ax,
                color=[palette[g] for g in pivot_df.columns],
                width=bar_width,
                edgecolor=bar_edgecolor,
                linewidth=bar_edgewidth
            )
            ax.set_ylabel('Ocean Basin', fontsize=11)
            ax.set_xlabel(value_label, fontsize=11)
        else:
            pivot_df.plot(
                kind='bar',
                stacked=True,
                ax=ax,
                color=[palette[g] for g in pivot_df.columns],
                width=bar_width,
                edgecolor=bar_edgecolor,
                linewidth=bar_edgewidth
            )
            ax.set_xlabel('Ocean Basin', fontsize=11)
            ax.set_ylabel(value_label, fontsize=11)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=axis_text_angle, ha='right', fontsize=axis_text_size)

        ax.legend(
            title=config.get('map', {}).get('legend_title', 'Genotype'),
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            frameon=True
        )

        # Add grid and ensure it's behind bars if configured
        ax.grid(axis='y' if orientation == 'vertical' else 'x', alpha=0.3)
        if grid_behind:
            ax.set_axisbelow(True)

        plt.tight_layout()

        return fig

    # Generate relative abundance plot
    fig_rel = create_bar_plot(
        df_rel,
        'pct_of_basin',
        'Percentage of basin',
        'Relative Abundance by Ocean Basin'
    )

    # Generate absolute abundance plot
    fig_abs = create_bar_plot(
        df_abs,
        'n_samples',
        'Number of samples',
        'Absolute Abundance by Ocean Basin'
    )

    # Save outputs
    if output_dir is None:
        output_dir = base_dir / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)

    formats = general_cfg.get('output_format', ['pdf', 'png'])
    dpi = general_cfg.get('dpi', 300)

    for fmt in formats:
        # Relative
        outfile_rel = output_dir / f"distribution_bar_relative_custom.{fmt}"
        fig_rel.savefig(outfile_rel, dpi=dpi, bbox_inches='tight')
        output_files.append(outfile_rel)
        logger.info(f"✓ Wrote {outfile_rel}")

        # Absolute
        outfile_abs = output_dir / f"distribution_bar_absolute_custom.{fmt}"
        fig_abs.savefig(outfile_abs, dpi=dpi, bbox_inches='tight')
        output_files.append(outfile_abs)
        logger.info(f"✓ Wrote {outfile_abs}")

    plt.close('all')

    return output_files


def regenerate_identity_distribution(
    base_dir: Path,
    output_dir: Optional[Path] = None,
    config: Optional[Dict[str, Any]] = None
) -> List[Path]:
    """
    Regenerate identity distribution histogram from exported data.

    Parameters
    ----------
    base_dir : Path
        Base directory containing plot_config.yaml and data/ subdirectory
    output_dir : Path, optional
        Output directory for plots (default: base_dir/visualization)
    config : dict, optional
        Configuration dictionary (if not provided, will load from plot_config.yaml)

    Returns
    -------
    list of Path
        Paths to generated plot files
    """
    # Load configuration
    if config is None:
        config_path = base_dir / "plot_config.yaml"
        config = load_config(config_path)

    colors_cfg = config.get('colors', {})
    filters = config.get('filters', {})
    identity_cfg = config.get('identity', {})
    general_cfg = config.get('general', {})

    # Load data
    data_dir = base_dir / "data"
    id_path = data_dir / "identity_distribution.csv"

    if not id_path.exists():
        raise FileNotFoundError(f"identity_distribution.csv not found in {data_dir}")

    df = pd.read_csv(id_path)

    if df.empty:
        raise ValueError("No rows in identity_distribution.csv")

    # Convert identity to percentage if needed
    df['identity'] = pd.to_numeric(df['identity'], errors='coerce')
    if df['identity'].max() <= 1:
        df['identity'] = df['identity'] * 100

    # Determine genotype column - handle missing consensus_group_sp
    if 'consensus_group_sp' in df.columns:
        df['genotype'] = df['consensus_group_sp'].fillna(df.get('consensus_group', ''))
    elif 'consensus_group' in df.columns:
        # Try to create consensus_group_sp from consensus_group + species
        if 'species' in df.columns:
            # Remove "consensus_" prefix and add species prefix
            df['genotype'] = df['species'] + ' ' + df['consensus_group'].str.replace('consensus_', '')
        else:
            df['genotype'] = df['consensus_group']
    else:
        df['genotype'] = ''

    # Filter out rows with missing genotypes
    df = df[df['genotype'].notna() & (df['genotype'] != '')]

    # Apply filters
    df = apply_genotype_filters(df, filters, 'genotype')

    if df.empty:
        raise ValueError("No identity data to plot after applying filters")

    # Build color palette
    genotypes = sorted(df['genotype'].unique())
    palette = get_palette(genotypes, colors_cfg)

    # Create figure
    # Use identity-specific dimensions if provided, otherwise fall back to general settings
    width = identity_cfg.get('width_inches', general_cfg.get('width_inches', 10))
    height = identity_cfg.get('height_inches', general_cfg.get('height_inches', 8))

    fig, ax = plt.subplots(figsize=(width, height))

    # Plot parameters
    binwidth = identity_cfg.get('binwidth', 0.5)
    show_mean = identity_cfg.get('show_mean', True)
    show_median = identity_cfg.get('show_median', True)
    show_density = identity_cfg.get('show_density', True)
    density_alpha = identity_cfg.get('density_alpha', 0.3)
    stat_line_color = convert_r_color_to_matplotlib(identity_cfg.get('stat_line_color', 'red'))
    stat_line_type = identity_cfg.get('stat_line_type', 'dashed')
    x_limits = identity_cfg.get('x_limits', None)
    x_breaks = identity_cfg.get('x_breaks', None)

    # Determine bins
    if x_limits:
        bins = np.arange(x_limits[0], x_limits[1] + binwidth, binwidth)
    else:
        bins = np.arange(df['identity'].min(), df['identity'].max() + binwidth, binwidth)

    # Plot histograms by genotype
    for genotype in genotypes:
        subset = df[df['genotype'] == genotype]['identity']
        if len(subset) > 0:
            ax.hist(
                subset,
                bins=bins,
                alpha=0.8,
                color=palette[genotype],
                label=genotype,
                edgecolor='white',
                linewidth=0.5
            )

            # Optionally add density curve
            if show_density and len(subset) > 1:
                from scipy import stats
                density = stats.gaussian_kde(subset)
                x_dense = np.linspace(subset.min(), subset.max(), 200)
                y_dense = density(x_dense) * len(subset) * binwidth
                ax.plot(x_dense, y_dense, color=palette[genotype], alpha=density_alpha, linewidth=2)

    # Add mean/median lines
    if show_mean:
        mean_val = df['identity'].mean()
        ax.axvline(mean_val, color=stat_line_color, linestyle='dashed', linewidth=2, label=f'Mean: {mean_val:.1f}%')

    if show_median:
        median_val = df['identity'].median()
        ax.axvline(median_val, color=stat_line_color, linestyle='dotted', linewidth=2, label=f'Median: {median_val:.1f}%')

    # Styling
    ax.set_xlabel('Percent identity', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.grid(axis='y', alpha=0.3)

    if x_limits:
        ax.set_xlim(x_limits)

    if x_breaks:
        ax.set_xticks(x_breaks)

    ax.legend(
        title=config.get('map', {}).get('legend_title', 'Genotype'),
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=True
    )

    plt.tight_layout()

    # Save outputs
    if output_dir is None:
        output_dir = base_dir / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    formats = general_cfg.get('output_format', ['pdf', 'png'])
    dpi = general_cfg.get('dpi', 300)

    for fmt in formats:
        outfile = output_dir / f"identity_distribution_custom.{fmt}"
        plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
        output_files.append(outfile)
        logger.info(f"✓ Wrote {outfile}")

    plt.close()

    return output_files


def regenerate_phylogenetic_tree(
    base_dir: Path,
    output_dir: Optional[Path] = None,
    config: Optional[Dict[str, Any]] = None
) -> List[Path]:
    """
    Regenerate phylogenetic tree from exported Newick file.

    Parameters
    ----------
    base_dir : Path
        Base directory containing plot_config.yaml and data/ subdirectory
    output_dir : Path, optional
        Output directory for plots (default: base_dir/visualization)
    config : dict, optional
        Configuration dictionary (if not provided, will load from plot_config.yaml)

    Returns
    -------
    list of Path
        Paths to generated plot files
    """
    try:
        from Bio import Phylo
    except ImportError:
        logger.warning("BioPython not available, skipping tree regeneration")
        return []

    # Load configuration
    if config is None:
        config_path = base_dir / "plot_config.yaml"
        config = load_config(config_path)

    colors_cfg = config.get('colors', {})
    filters = config.get('filters', {})
    tree_cfg = config.get('tree', {})
    general_cfg = config.get('general', {})

    # Load tree data
    data_dir = base_dir / "data"
    tree_path = data_dir / "tree_data.nwk"

    if not tree_path.exists():
        logger.warning(f"tree_data.nwk not found in {data_dir}, skipping tree regeneration")
        return []

    # Load consensus_group to consensus_group_sp mapping
    # Try multiple data files to find the mapping
    label_map = {}
    for csv_file in ['distribution_bar_relative.csv', 'distribution_bar_absolute.csv',
                     'distribution_map.csv', 'identity_distribution.csv']:
        csv_path = data_dir / csv_file
        if csv_path.exists():
            try:
                df_mapping = pd.read_csv(csv_path)
                if 'consensus_group' in df_mapping.columns and 'consensus_group_sp' in df_mapping.columns:
                    # Create mapping, taking first occurrence for each consensus_group
                    mapping_df = df_mapping[['consensus_group', 'consensus_group_sp']].drop_duplicates()

                    # Build label map with fuzzy matching
                    # Tree tips may have short names like "consensus_c14" while mapping has "consensus_c14_n11"
                    for _, row in mapping_df.iterrows():
                        full_name = row['consensus_group']
                        label = row['consensus_group_sp']

                        # Add full name mapping
                        label_map[full_name] = label

                        # Also add short name mapping (e.g., "consensus_c14" from "consensus_c14_n11")
                        if '_n' in full_name:
                            short_name = full_name.split('_n')[0]
                            if short_name not in label_map:
                                label_map[short_name] = label

                    logger.info(f"Loaded {len(label_map)} tip label mappings from {csv_file}")
                    break
            except Exception as e:
                logger.debug(f"Could not load mapping from {csv_file}: {e}")
                continue

    # Read tree
    tree = Phylo.read(str(tree_path), "newick")

    # Get tree parameters from config (check for outgroup before modifying tree)
    outgroup = tree_cfg.get('outgroup', None)

    # Reroot tree with outgroup if specified
    if outgroup and label_map:
        # Find the outgroup tip - try exact match first, then fuzzy match
        outgroup_clade = None

        # Try to find matching tip
        for tip in tree.get_terminals():
            # Check if tip name matches outgroup (before or after relabeling)
            if tip.name == outgroup:
                outgroup_clade = tip
                break
            # Check if tip will match after relabeling
            elif tip.name in label_map and label_map[tip.name] == outgroup:
                outgroup_clade = tip
                break

        if outgroup_clade:
            tree.root_with_outgroup(outgroup_clade)
            logger.info(f"Rerooted tree with outgroup: {outgroup}")
        else:
            logger.warning(f"Outgroup '{outgroup}' not found in tree tips, skipping rerooting")

    # Apply label mapping to tree tips (consensus_group -> consensus_group_sp)
    relabeled_count = 0
    if label_map:
        for clade in tree.get_terminals():
            if clade.name in label_map:
                old_name = clade.name
                clade.name = label_map[clade.name]
                logger.debug(f"Relabeled tip: {old_name} -> {clade.name}")
                relabeled_count += 1
        logger.info(f"Applied consensus_group_sp labels to {relabeled_count}/{len(tree.get_terminals())} tree tips")

    # Get tree parameters from config
    layout = tree_cfg.get('layout', 'rectangular')
    show_bootstrap = tree_cfg.get('show_bootstrap', True)
    bootstrap_threshold = tree_cfg.get('bootstrap_threshold', 70)
    bootstrap_size = tree_cfg.get('bootstrap_size', 8)
    bootstrap_offset_x = tree_cfg.get('bootstrap_offset_x', 0.0)
    bootstrap_offset_y = tree_cfg.get('bootstrap_offset_y', 0.0)
    tip_label_size = tree_cfg.get('tip_label_size', 10)
    tip_label_offset = tree_cfg.get('tip_label_offset', 0.001)
    branch_width = tree_cfg.get('branch_width', 1.0)
    show_scale_bar = tree_cfg.get('show_scale_bar', True)
    highlight_groups = tree_cfg.get('highlight_groups', [])

    # Build color palette for tips (after relabeling to consensus_group_sp)
    terminals = tree.get_terminals()
    tip_names = [t.name for t in terminals]

    # Match tip names (now consensus_group_sp format) to configured colors
    palette = {}
    for tip_name in tip_names:
        # Check for exact match first (e.g., "Sphyrna lewini c18_n281")
        if tip_name in colors_cfg:
            palette[tip_name] = colors_cfg[tip_name]
            logger.debug(f"Color matched (exact): {tip_name} -> {colors_cfg[tip_name]}")
        else:
            # Try partial matching by cluster ID (e.g., "c18_n281" in "Sphyrna lewini c18_n281")
            for color_key, color_val in colors_cfg.items():
                # Extract cluster from tip name (e.g., "c18_n281" from "Sphyrna lewini c18_n281")
                parts = tip_name.split()
                if len(parts) >= 2 and parts[-1].startswith('c') and '_n' in parts[-1]:
                    tip_cluster = parts[-1]  # e.g., "c18_n281"
                    # Check if this cluster appears in the color key
                    if tip_cluster in color_key:
                        palette[tip_name] = color_val
                        logger.debug(f"Color matched (cluster): {tip_name} -> {color_val}")
                        break

    logger.info(f"Matched colors for {len(palette)}/{len(tip_names)} tree tips")

    # Compute x (depth) and y positions for each clade
    depths = tree.depths()
    if not depths:
        depths = tree.depths(unit_branch_lengths=True)

    # y positions: leaves are evenly spaced, internal nodes are midpoints
    y_pos: Dict[Any, float] = {
        clade: float(i) for i, clade in enumerate(terminals)
    }

    def calc_y(clade):
        if clade in y_pos:
            return y_pos[clade]
        if not clade.clades:
            y_pos[clade] = 0.0
            return 0.0
        child_ys = [calc_y(child) for child in clade.clades]
        y_val = (min(child_ys) + max(child_ys)) / 2.0
        y_pos[clade] = y_val
        return y_val

    calc_y(tree.root)

    # Calculate figure size
    # Use tree-specific dimensions if provided, otherwise fall back to general settings
    # If neither are specified, auto-calculate based on tree size
    n_tips = len(terminals)

    # Width: tree config > general config > default 10
    width = tree_cfg.get('width_inches', general_cfg.get('width_inches', 10))

    # Height: tree config > general config > auto-calculated based on number of tips
    if 'height_inches' in tree_cfg:
        height = tree_cfg['height_inches']
    elif 'height_inches' in general_cfg:
        height = general_cfg['height_inches']
    else:
        # Auto-calculate: 0.3 inches per tip, min 6, max 20
        height = max(6, min(20, n_tips * 0.3))

    logger.info(f"Tree figure size: {width}\" × {height}\" ({n_tips} tips)")

    # Create figure
    fig, ax = plt.subplots(figsize=(width, height))

    # Draw tree using Phylo
    # Map tip colors if available
    tip_label_colors = None
    if palette:
        tip_label_colors = {}
        for clade in terminals:
            if clade.name in palette:
                tip_label_colors[clade] = palette[clade.name]

    Phylo.draw(
        tree,
        do_show=False,
        axes=ax,
        label_colors=tip_label_colors,
        show_confidence=False  # We'll add bootstrap values manually for positioning control
    )

    # Extend x-axis to prevent tip labels from being cut off
    xlim = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    ax.set_xlim(xlim[0], xlim[1] + x_range * 0.15)

    # Add colored circular nodes at tips
    if palette:
        xs, ys, cs = [], [], []
        for clade in terminals:
            if clade.name in palette:
                x = depths.get(clade, 0.0)
                y = y_pos.get(clade, 0.0)
                xs.append(x)
                ys.append(y)
                cs.append(palette[clade.name])

        if xs:
            ax.scatter(
                xs, ys,
                s=80,
                c=cs,
                marker='o',
                edgecolors='black',
                linewidths=1.0,
                zorder=3,
                clip_on=False
            )

    # Add bootstrap values with configurable positioning
    if show_bootstrap:
        for clade in tree.get_nonterminals():
            if clade.confidence is not None:
                val = float(clade.confidence)
                if val >= bootstrap_threshold:
                    x = depths.get(clade, 0.0)
                    y = y_pos.get(clade, 0.0)

                    # Apply configurable offsets for fine-tuning label positions
                    label_x = x + bootstrap_offset_x
                    label_y = y + bootstrap_offset_y

                    ax.text(
                        label_x, label_y,
                        f'{int(val)}',
                        fontsize=bootstrap_size,
                        ha='center',
                        va='bottom',
                        color='red',
                        weight='bold'
                    )

    # Add scale bar if requested
    if show_scale_bar:
        # Get the x-axis limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Calculate a reasonable scale bar length (10-20% of x range)
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]

        # Round to a nice number
        scale_length = 10 ** (np.floor(np.log10(x_range * 0.1)))

        # Position scale bar at bottom left
        scale_x = xlim[0] + x_range * 0.05
        scale_y = ylim[0] - y_range * 0.05

        # Draw scale bar
        ax.plot([scale_x, scale_x + scale_length], [scale_y, scale_y],
                'k-', linewidth=2, solid_capstyle='butt')

        # Add label
        ax.text(scale_x + scale_length / 2, scale_y - y_range * 0.02,
                f'{scale_length:.3g}',
                ha='center', va='top', fontsize=8)

    # Remove axis labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tight_layout()

    # Save outputs
    if output_dir is None:
        output_dir = base_dir / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    formats = general_cfg.get('output_format', ['pdf', 'png'])
    dpi = general_cfg.get('dpi', 300)

    for fmt in formats:
        outfile = output_dir / f"phylogenetic_tree_custom.{fmt}"
        plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
        output_files.append(outfile)
        logger.info(f"✓ Wrote {outfile}")

    plt.close()

    return output_files


def regenerate_all_plots(base_dir: Path, output_dir: Optional[Path] = None) -> Dict[str, List[Path]]:
    """
    Regenerate all plots from exported data.

    Parameters
    ----------
    base_dir : Path
        Base directory containing plot_config.yaml and data/ subdirectory
    output_dir : Path, optional
        Output directory for plots (default: base_dir/visualization)

    Returns
    -------
    dict
        Dictionary mapping plot type to list of output files
    """
    logger.info("Regenerating all plots...")

    # Load config once
    config_path = base_dir / "plot_config.yaml"
    config = load_config(config_path)

    results = {}

    # Regenerate distribution map
    try:
        logger.info("  1. Regenerating distribution map...")
        results['map'] = regenerate_distribution_map(base_dir, output_dir, config)
    except Exception as e:
        logger.error(f"  Failed to regenerate map: {e}")
        results['map'] = []

    # Regenerate bar charts
    try:
        logger.info("  2. Regenerating bar charts...")
        results['bars'] = regenerate_bar_charts(base_dir, output_dir, config)
    except Exception as e:
        logger.error(f"  Failed to regenerate bars: {e}")
        results['bars'] = []

    # Regenerate identity distribution
    try:
        logger.info("  3. Regenerating identity distribution...")
        results['identity'] = regenerate_identity_distribution(base_dir, output_dir, config)
    except Exception as e:
        logger.error(f"  Failed to regenerate identity: {e}")
        results['identity'] = []

    # Regenerate phylogenetic tree
    try:
        logger.info("  4. Regenerating phylogenetic tree...")
        results['tree'] = regenerate_phylogenetic_tree(base_dir, output_dir, config)
    except Exception as e:
        logger.error(f"  Failed to regenerate tree: {e}")
        results['tree'] = []

    total_files = sum(len(files) for files in results.values())
    logger.info(f"✓ All plots regenerated successfully! ({total_files} files)")

    if output_dir is None:
        output_dir = base_dir / "visualization"
    logger.info(f"  Check {output_dir} for updated plots")

    return results
