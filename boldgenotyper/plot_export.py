"""
Plot Data Export Module for BOLDGenotyper

This module exports raw data and scripts to regenerate all plots with custom styling.
Users can modify plot aesthetics, colors, and formats for publication without rerunning
the entire analysis.

Key Features:
- Export plot data to CSV files
- Generate plot_config.yaml with all styling parameters
- Create R scripts for plot regeneration
- Provide examples and documentation

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Optional, Any
from pathlib import Path
import logging
import shutil

import pandas as pd

# Try to import YAML support
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

logger = logging.getLogger(__name__)


def export_plot_data(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    consensus_path: Optional[Path] = None,
    tree_path: Optional[Path] = None,
    diagnostics_path: Optional[Path] = None,
    color_map: Optional[Dict[str, str]] = None
) -> Dict[str, Path]:
    """
    Export all plot data to CSV files for custom regeneration.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with all analysis results
    output_dir : Path
        Base output directory
    organism : str
        Organism name
    consensus_path : Path, optional
        Path to consensus sequences FASTA
    tree_path : Path, optional
        Path to phylogenetic tree file
    diagnostics_path : Path, optional
        Path to assignment diagnostics CSV
    color_map : dict, optional
        Mapping of genotypes to colors

    Returns
    -------
    dict
        Dictionary mapping data type to output path
    """
    plot_data_dir = output_dir / "plots" / "data"
    plot_data_dir.mkdir(parents=True, exist_ok=True)

    exported_files = {}

    # 1. Distribution map data
    if 'lat' in df.columns and 'lon' in df.columns and 'consensus_group' in df.columns:
        logger.info("  Exporting distribution map data...")
        map_cols = ['processid', 'consensus_group', 'species', 'lat', 'lon']
        if 'ocean_basin' in df.columns:
            map_cols.append('ocean_basin')
        if 'consensus_group_sp' in df.columns:
            map_cols.append('consensus_group_sp')

        # Add n_at_location for point sizing
        map_df = df[map_cols].copy()
        map_df['location_key'] = map_df['lat'].round(4).astype(str) + '_' + map_df['lon'].round(4).astype(str)
        location_counts = map_df.groupby('location_key').size().rename('n_at_location')
        map_df = map_df.merge(location_counts, left_on='location_key', right_index=True)
        map_df = map_df.drop('location_key', axis=1)

        map_path = plot_data_dir / "distribution_map.csv"
        map_df.to_csv(map_path, index=False)
        exported_files['distribution_map'] = map_path

    # 2. Distribution bar data (relative)
    if 'ocean_basin' in df.columns and 'consensus_group' in df.columns:
        logger.info("  Exporting distribution bar data...")
        basin_counts = df.groupby(['ocean_basin', 'consensus_group', 'species']).size().reset_index(name='n_samples')

        # Calculate percentages within each basin
        basin_totals = basin_counts.groupby('ocean_basin')['n_samples'].transform('sum')
        basin_counts['pct_of_basin'] = (basin_counts['n_samples'] / basin_totals * 100).round(2)

        if 'consensus_group_sp' in df.columns:
            # Add species labels
            group_sp_map = df[['consensus_group', 'consensus_group_sp']].drop_duplicates()
            basin_counts = basin_counts.merge(group_sp_map, on='consensus_group', how='left')

        bar_rel_path = plot_data_dir / "distribution_bar_relative.csv"
        basin_counts.to_csv(bar_rel_path, index=False)
        exported_files['distribution_bar_relative'] = bar_rel_path

        # 3. Distribution bar data (absolute)
        bar_abs_path = plot_data_dir / "distribution_bar_absolute.csv"
        # Include consensus_group_sp if it exists
        abs_cols = ['ocean_basin', 'consensus_group', 'species', 'n_samples']
        if 'consensus_group_sp' in basin_counts.columns:
            abs_cols.append('consensus_group_sp')
        basin_counts[abs_cols].to_csv(bar_abs_path, index=False)
        exported_files['distribution_bar_absolute'] = bar_abs_path

    # 4. Identity distribution data
    if diagnostics_path and diagnostics_path.exists():
        logger.info("  Exporting identity distribution data...")
        diag_df = pd.read_csv(diagnostics_path)

        identity_cols = ['processid', 'consensus_group', 'identity']
        if 'is_tie' in diag_df.columns:
            identity_cols.append('is_tie')
        if 'is_low_confidence' in diag_df.columns:
            identity_cols.append('is_low_confidence')

        # Merge with species info and consensus_group_sp if available
        merge_cols = ['processid']
        if 'species' in df.columns:
            merge_cols.append('species')
            identity_cols.append('species')
        if 'consensus_group_sp' in df.columns:
            merge_cols.append('consensus_group_sp')
            identity_cols.append('consensus_group_sp')

        if len(merge_cols) > 1:
            # Need to get consensus_group_sp from main df by matching on processid + consensus_group
            merge_df = df[merge_cols + ['consensus_group']].drop_duplicates()
            diag_df = diag_df.merge(merge_df, on=['processid', 'consensus_group'], how='left', suffixes=('', '_y'))

        identity_path = plot_data_dir / "identity_distribution.csv"
        diag_df[identity_cols].to_csv(identity_path, index=False)
        exported_files['identity_distribution'] = identity_path

    # 5. Genotype colors
    if color_map:
        logger.info("  Exporting genotype color map...")
        color_df = pd.DataFrame([
            {'genotype': g, 'color_hex': c}
            for g, c in color_map.items()
        ])

        # Add species and sample counts
        if 'consensus_group_sp' in df.columns:
            group_info = df.groupby('consensus_group').agg({
                'species': lambda x: x.mode()[0] if len(x.mode()) > 0 else 'Unknown',
                'processid': 'count',
                'consensus_group_sp': lambda x: x.iloc[0] if len(x) > 0 else ''
            }).reset_index()
            group_info.columns = ['genotype', 'species', 'n_samples', 'genotype_label']

            # Merge with color map, using consensus_group_sp as genotype key in color map
            # First map consensus_group to consensus_group_sp
            genotype_sp_map = df[['consensus_group', 'consensus_group_sp']].drop_duplicates()
            color_df_expanded = []
            for _, row in genotype_sp_map.iterrows():
                cg = row['consensus_group']
                cg_sp = row['consensus_group_sp']
                if cg_sp in color_map:
                    color_hex = color_map[cg_sp]
                elif cg in color_map:
                    color_hex = color_map[cg]
                else:
                    color_hex = '#808080'  # default gray

                color_df_expanded.append({
                    'genotype': cg,
                    'genotype_label': cg_sp,
                    'color_hex': color_hex
                })

            color_df = pd.DataFrame(color_df_expanded)

            # Merge with group_info
            color_df = color_df.merge(
                group_info[['genotype', 'species', 'n_samples']],
                on='genotype',
                how='left'
            )

        color_path = plot_data_dir / "genotype_colors.csv"
        color_df.to_csv(color_path, index=False)
        exported_files['genotype_colors'] = color_path

    # 6. Tree data (copy if exists)
    if tree_path and tree_path.exists():
        logger.info("  Copying tree file...")
        tree_copy_path = plot_data_dir / "tree_data.nwk"
        shutil.copy(tree_path, tree_copy_path)
        exported_files['tree'] = tree_copy_path

    # 7. Consensus sequences (copy if exists)
    if consensus_path and consensus_path.exists():
        logger.info("  Copying consensus sequences...")
        consensus_copy_path = plot_data_dir / "consensus_sequences.fasta"
        shutil.copy(consensus_path, consensus_copy_path)
        exported_files['consensus'] = consensus_copy_path

    logger.info(f"  ✓ Exported {len(exported_files)} plot data files to {plot_data_dir}")

    return exported_files


def create_plot_config(
    output_dir: Path,
    organism: str,
    color_map: Optional[Dict[str, str]] = None
) -> Path:
    """
    Create plot_config.yaml with all plotting parameters.

    Parameters
    ----------
    output_dir : Path
        Base output directory
    organism : str
        Organism name
    color_map : dict, optional
        Mapping of genotypes to colors

    Returns
    -------
    Path
        Path to created config file
    """
    plots_dir = output_dir / "plots"
    config_path = plots_dir / "plot_config.yaml"

    # Default configuration
    config = {
        'general': {
            'output_format': ['pdf', 'png', 'svg'],  # Include SVG for vector graphics
            'dpi': 600,  # High resolution for publication
            'width_inches': 10,  # Default width (can be overridden per plot)
            'height_inches': 8   # Default height (can be overridden per plot)
        },
        'colors': color_map or {},
        'filters': {
            'include_genotypes': [],
            'exclude_genotypes': [],
            'exclude_unknown_geography': False  # Set to true to omit 'Unknown' ocean basins
        },
        'map': {
            # Map-specific dimensions (overrides general settings)
            'width_inches': 12,
            'height_inches': 8,
            # Projection settings
            'projection': 'robinson',
            'center_longitude': 0,
            'show_country_borders': True,
            'border_color': 'gray70',
            'border_width': 0.3,
            'ocean_color': '#E8F4F8',
            'land_color': '#F5F5F5',
            'point_alpha': 0.7,
            'point_size_range': [2, 10],
            'point_stroke': 0.5,
            # Legend settings
            'show_legend': True,
            'legend_position': 'bottom',  # Options: right, left, top, bottom
            'legend_title': 'Genotype',
            # Gridline settings for lat/lon coordinates
            'show_gridlines': True,
            'gridline_labels': True,  # Show lat/lon coordinate labels
            'gridline_alpha': 0.5,
            'gridline_style': '--'
        },
        'bars': {
            # Bar chart dimensions (overrides general settings)
            'width_inches': 10,
            'height_inches': 6,
            # Layout
            'orientation': 'vertical',
            'bar_width': 0.8,
            'show_percentages': True,
            'percentage_size': 3,
            'show_sample_counts': True,
            'facet_scales': 'free_y',
            'axis_text_angle': 45,
            'axis_text_size': 10,
            'color_palette_type': 'discrete',
            # Bar styling
            'bar_edgecolor': 'black',  # Outline color for bars
            'bar_edgewidth': 0.5,      # Outline width for bars
            'grid_behind_bars': True,  # Place grid behind bars (not on top)
            # Recalculate relative abundance for filtered genotypes
            'recalculate_relative': True  # Recalculate percentages to sum to 100% after filtering
        },
        'identity': {
            # Identity plot dimensions (overrides general settings)
            'width_inches': 8,
            'height_inches': 6,
            # Histogram settings
            'binwidth': 0.5,
            'show_mean': True,
            'show_median': True,
            'show_density': True,
            'density_alpha': 0.3,
            'stat_line_color': 'red',
            'stat_line_type': 'dashed',
            'x_limits': [50, 100],
            'x_breaks': [50, 60, 70, 80, 90, 100]
        },
        'tree': {
            'layout': 'rectangular',
            'show_bootstrap': True,
            'bootstrap_threshold': 70,
            'bootstrap_size': 3,
            'tip_label_size': 3,
            'tip_label_offset': 0.001,
            'branch_width': 0.5,
            'show_scale_bar': True,
            'highlight_groups': []
        }
    }

    # Write config file
    if YAML_AVAILABLE:
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    else:
        # Fallback: write manual YAML format
        with open(config_path, 'w') as f:
            f.write("# Plot Configuration for BOLDGenotyper\n")
            f.write("# Edit values to customize plots, then run regenerate_all.sh\n\n")

            def write_dict(d, indent=0):
                """Recursively write dictionary in YAML format."""
                for key, value in d.items():
                    if isinstance(value, dict):
                        f.write('  ' * indent + f"{key}:\n")
                        write_dict(value, indent + 1)
                    elif isinstance(value, list):
                        f.write('  ' * indent + f"{key}:\n")
                        for item in value:
                            if isinstance(item, str):
                                f.write('  ' * (indent + 1) + f"- {item}\n")
                            else:
                                f.write('  ' * (indent + 1) + f"- {item}\n")
                    elif isinstance(value, str):
                        f.write('  ' * indent + f"{key}: \"{value}\"\n")
                    elif isinstance(value, bool):
                        f.write('  ' * indent + f"{key}: {str(value).lower()}\n")
                    else:
                        f.write('  ' * indent + f"{key}: {value}\n")

            write_dict(config)

    logger.info(f"  ✓ Created plot configuration: {config_path}")

    return config_path


def create_python_scripts(output_dir: Path, organism: str) -> Dict[str, Path]:
    """
    Create Python scripts for regenerating plots.

    Parameters
    ----------
    output_dir : Path
        Base output directory
    organism : str
        Organism name

    Returns
    -------
    dict
        Dictionary mapping script name to path
    """
    scripts_dir = output_dir / "plots" / "scripts"
    scripts_dir.mkdir(parents=True, exist_ok=True)

    scripts = {}

    # README for regeneration
    readme_content = f"""# Plot Regeneration Guide for {organism}

## Quick Start

1. Ensure you have the boldgenotyper conda environment activated:
   ```bash
   conda activate boldgenotyper
   ```

2. Modify `plot_config.yaml` to customize colors, sizes, labels

3. Run all scripts:
   ```bash
   python regenerate_all.py
   ```

4. Or run individual scripts:
   ```bash
   python regenerate_map.py
   python regenerate_bars.py
   python regenerate_identity.py
   ```

## Files

- **regenerate_all.py**: Run all regeneration scripts
- **regenerate_map.py**: Regenerate distribution maps
- **regenerate_bars.py**: Regenerate bar charts
- **regenerate_identity.py**: Regenerate identity histograms

## Requirements

All required packages are included in the boldgenotyper conda environment:
- matplotlib
- seaborn
- pandas
- geopandas
- cartopy
- pyyaml

No additional installation needed!

## Customization

Edit `plot_config.yaml` to change:
- Color palettes
- Figure dimensions and DPI
- Map projections and styling
- Bar chart orientations
- Histogram bin widths
- Genotype filters (include/exclude specific consensus groups)

## Output

Regenerated plots are saved to the parent `visualization/` directory,
replacing the original plots.
"""

    readme_path = output_dir / "plots" / "README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    scripts['readme'] = readme_path

    # Master regeneration script
    regenerate_all = """#!/usr/bin/env python3
\"\"\"Regenerate all plots from exported data.\"\"\"

import sys
from pathlib import Path

# Add parent directory to path to import boldgenotyper
parent_dir = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(parent_dir))

from boldgenotyper.plot_regeneration import regenerate_all_plots

def main():
    # Base directory is parent of scripts/
    scripts_dir = Path(__file__).parent
    base_dir = scripts_dir.parent

    print("Regenerating all plots...")
    print()

    try:
        results = regenerate_all_plots(base_dir)

        total_files = sum(len(files) for files in results.values())
        print()
        print(f"✓ All plots regenerated successfully! ({total_files} files)")
        print(f"  Check {base_dir / 'visualization'} for updated plots")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
"""

    regenerate_all_path = scripts_dir / "regenerate_all.py"
    with open(regenerate_all_path, 'w') as f:
        f.write(regenerate_all)
    regenerate_all_path.chmod(0o755)  # Make executable
    scripts['regenerate_all'] = regenerate_all_path

    # Individual script templates
    map_script = """#!/usr/bin/env python3
\"\"\"Regenerate distribution map from exported data.\"\"\"

import sys
from pathlib import Path

# Add parent directory to path to import boldgenotyper
parent_dir = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(parent_dir))

from boldgenotyper.plot_regeneration import regenerate_distribution_map

def main():
    # Base directory is parent of scripts/
    scripts_dir = Path(__file__).parent
    base_dir = scripts_dir.parent

    print("Regenerating distribution map...")

    try:
        output_files = regenerate_distribution_map(base_dir)
        print()
        print(f"✓ Generated {len(output_files)} file(s)")
        for f in output_files:
            print(f"  - {f}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
"""

    bars_script = """#!/usr/bin/env python3
\"\"\"Regenerate bar charts from exported data.\"\"\"

import sys
from pathlib import Path

# Add parent directory to path to import boldgenotyper
parent_dir = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(parent_dir))

from boldgenotyper.plot_regeneration import regenerate_bar_charts

def main():
    # Base directory is parent of scripts/
    scripts_dir = Path(__file__).parent
    base_dir = scripts_dir.parent

    print("Regenerating bar charts...")

    try:
        output_files = regenerate_bar_charts(base_dir)
        print()
        print(f"✓ Generated {len(output_files)} file(s)")
        for f in output_files:
            print(f"  - {f}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
"""

    identity_script = """#!/usr/bin/env python3
\"\"\"Regenerate identity distribution histogram from exported data.\"\"\"

import sys
from pathlib import Path

# Add parent directory to path to import boldgenotyper
parent_dir = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(parent_dir))

from boldgenotyper.plot_regeneration import regenerate_identity_distribution

def main():
    # Base directory is parent of scripts/
    scripts_dir = Path(__file__).parent
    base_dir = scripts_dir.parent

    print("Regenerating identity distribution...")

    try:
        output_files = regenerate_identity_distribution(base_dir)
        print()
        print(f"✓ Generated {len(output_files)} file(s)")
        for f in output_files:
            print(f"  - {f}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
"""

    # Write Python scripts
    map_path = scripts_dir / "regenerate_map.py"
    with open(map_path, 'w') as f:
        f.write(map_script)
    map_path.chmod(0o755)
    scripts['regenerate_map'] = map_path

    bars_path = scripts_dir / "regenerate_bars.py"
    with open(bars_path, 'w') as f:
        f.write(bars_script)
    bars_path.chmod(0o755)
    scripts['regenerate_bars'] = bars_path

    identity_path = scripts_dir / "regenerate_identity.py"
    with open(identity_path, 'w') as f:
        f.write(identity_script)
    identity_path.chmod(0o755)
    scripts['regenerate_identity'] = identity_path

    logger.info(f"  ✓ Created Python regeneration scripts in {scripts_dir}")

    return scripts


def create_readme(output_dir: Path, organism: str, exported_files: Dict[str, Path]) -> Path:
    """
    Create README.md in plots directory with instructions.

    Parameters
    ----------
    output_dir : Path
        Base output directory
    organism : str
        Organism name
    exported_files : dict
        Dictionary of exported files

    Returns
    -------
    Path
        Path to README file
    """
    plots_dir = output_dir / "plots"
    readme_path = plots_dir / "README.md"

    file_list = "\n".join([f"- `{path.name}`" for path in exported_files.values()])

    content = f"""# Plot Data Export for {organism}

## Overview

This directory contains raw data and scripts to regenerate all plots with custom
styling for publication. You can modify colors, labels, sizes, and formats without
rerunning the entire BOLDGenotyper analysis.

## Directory Structure

```
plots/
├── README.md                   # This file
├── plot_config.yaml           # All plotting parameters
├── data/                      # Raw data for each plot
{file_list}
├── scripts/                   # R regeneration scripts
│   ├── regenerate_all.sh
│   ├── regenerate_map.R
│   ├── regenerate_bars.R
│   ├── regenerate_identity.R
│   └── requirements.txt
└── examples/                  # Example modifications (TODO)
```

## Quick Start

1. **Install R and required packages:**
   ```R
   install.packages(c("ggplot2", "dplyr", "yaml", "sf", "rnaturalearth"))
   ```

2. **Modify `plot_config.yaml`** to customize:
   - Colors for each genotype
   - Figure dimensions and DPI
   - Map projection and styling
   - Bar chart layout
   - Histogram binning

3. **Regenerate plots:**
   ```bash
   cd scripts/
   bash regenerate_all.sh
   ```

4. **Find updated plots** in `../visualization/`

## Customization Examples

### Change Color Palette

Edit `plot_config.yaml`:
```yaml
colors:
  c15_n386: "#E41A1C"  # Red
  c18_n255: "#377EB8"  # Blue
  c16_n86: "#4DAF4A"   # Green
```

### Adjust Figure Size

```yaml
general:
  width_inches: 12
  height_inches: 8
  dpi: 600  # High resolution for publication
```

### Modify Map Projection

```yaml
map:
  projection: "mollweide"  # or "mercator", "robinson"
  center_longitude: -180   # Center on Pacific
```

## Data Files

{file_list}

## Notes

- All data files are in CSV format for easy editing
- Tree data is in Newick format
- Consensus sequences are in FASTA format
- Original plots remain unchanged until you regenerate

## Support

For issues or questions, see the BOLDGenotyper documentation or file an issue on GitHub.
"""

    with open(readme_path, 'w') as f:
        f.write(content)

    logger.info(f"  ✓ Created README: {readme_path}")

    return readme_path


def export_plots_complete(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    consensus_path: Optional[Path] = None,
    tree_path: Optional[Path] = None,
    diagnostics_path: Optional[Path] = None,
    color_map: Optional[Dict[str, str]] = None
) -> Dict[str, Any]:
    """
    Complete plot export workflow: data, config, scripts, and documentation.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe
    output_dir : Path
        Base output directory
    organism : str
        Organism name
    consensus_path : Path, optional
        Path to consensus FASTA
    tree_path : Path, optional
        Path to tree file
    diagnostics_path : Path, optional
        Path to diagnostics CSV
    color_map : dict, optional
        Genotype to color mapping

    Returns
    -------
    dict
        Summary of exported files and paths
    """
    logger.info("Exporting plot data and regeneration kit...")

    results = {}

    # Export data files
    exported_files = export_plot_data(
        df, output_dir, organism,
        consensus_path=consensus_path,
        tree_path=tree_path,
        diagnostics_path=diagnostics_path,
        color_map=color_map
    )
    results['data_files'] = exported_files

    # Create config
    config_path = create_plot_config(output_dir, organism, color_map=color_map)
    results['config'] = config_path

    # Create Python scripts
    scripts = create_python_scripts(output_dir, organism)
    results['scripts'] = scripts

    # Create README
    readme_path = create_readme(output_dir, organism, exported_files)
    results['readme'] = readme_path

    logger.info(f"✓ Plot regeneration kit exported to {output_dir / 'plots'}")

    return results
