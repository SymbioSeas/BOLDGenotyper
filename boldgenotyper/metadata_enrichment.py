#!/usr/bin/env python3
"""
Metadata Enrichment Tools for BOLDGenotyper

This module allows users to:
1. Add custom metadata to analysis results
2. Update ocean basin assignments with new shapefiles
3. Add custom groupings for analysis
4. Regenerate visualizations with enriched data
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json
from datetime import datetime
from typing import Optional, Dict, List, Tuple

logger = logging.getLogger(__name__)


def load_annotated_csv(csv_path: Path) -> pd.DataFrame:
    """
    Load the annotated CSV from BOLDGenotyper output.

    Args:
        csv_path: Path to annotated CSV file

    Returns:
        DataFrame with annotated data
    """
    logger.info(f"Loading annotated CSV from {csv_path}")

    if not csv_path.exists():
        raise FileNotFoundError(f"Annotated CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    logger.info(f"Loaded {len(df)} samples with {len(df.columns)} columns")

    return df


def merge_metadata(
    base_df: pd.DataFrame,
    metadata_path: Path,
    join_column: str = 'processid',
    validate: bool = True
) -> Tuple[pd.DataFrame, Dict]:
    """
    Merge additional metadata into the base dataframe.

    Args:
        base_df: Base dataframe to enrich
        metadata_path: Path to CSV with additional metadata
        join_column: Column name to join on
        validate: Validate join operation

    Returns:
        Tuple of (enriched dataframe, report dict)
    """
    logger.info(f"Merging metadata from {metadata_path}")

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    # Load metadata
    metadata_df = pd.read_csv(metadata_path)
    logger.info(f"Loaded metadata: {len(metadata_df)} rows, {len(metadata_df.columns)} columns")

    # Check join column exists
    if join_column not in base_df.columns:
        raise ValueError(f"Join column '{join_column}' not found in base data")
    if join_column not in metadata_df.columns:
        raise ValueError(f"Join column '{join_column}' not found in metadata file")

    # Identify new columns
    new_columns = [col for col in metadata_df.columns if col != join_column and col not in base_df.columns]
    overlapping_columns = [col for col in metadata_df.columns if col in base_df.columns and col != join_column]

    logger.info(f"New columns to add: {new_columns}")
    if overlapping_columns:
        logger.warning(f"Overlapping columns (will be suffixed): {overlapping_columns}")

    # Perform merge
    merged_df = base_df.merge(
        metadata_df,
        on=join_column,
        how='left',
        suffixes=('', '_new'),
        validate='many_to_one' if validate else None
    )

    # Generate report
    n_matched = merged_df[new_columns[0]].notna().sum() if new_columns else 0
    n_unmatched = len(merged_df) - n_matched

    report = {
        'metadata_file': str(metadata_path),
        'join_column': join_column,
        'new_columns': new_columns,
        'overlapping_columns': overlapping_columns,
        'samples_in_base': len(base_df),
        'samples_in_metadata': len(metadata_df),
        'samples_matched': n_matched,
        'samples_unmatched': n_unmatched,
        'match_rate': n_matched / len(base_df) * 100 if len(base_df) > 0 else 0
    }

    logger.info(f"Merge complete: {n_matched}/{len(base_df)} samples matched ({report['match_rate']:.1f}%)")

    return merged_df, report


def update_geographic_regions(
    df: pd.DataFrame,
    shapefile_path: Optional[Path] = None,
    shapefile_field: str = 'name',
    geo_category: str = 'ocean_basin',
    force_recalculate: bool = False
) -> Tuple[pd.DataFrame, Dict]:
    """
    Update geographic region assignments using a shapefile.

    This function is flexible and works with any shapefile containing geographic
    regions (ocean basins, freshwater basins, ecoregions, countries, etc.).

    Args:
        df: DataFrame with lat/lon columns
        shapefile_path: Path to shapefile (if None, uses default ocean basins)
        shapefile_field: Name of shapefile attribute containing region names (default: 'name')
        geo_category: Name for geographic category (default: 'ocean_basin')
                     Examples: 'ecoregion', 'watershed', 'biome', 'country'
        force_recalculate: Recalculate all regions even if already assigned

    Returns:
        Tuple of (updated dataframe, report dict)
    """
    logger.info(f"Updating geographic region assignments (category: {geo_category})")

    # Check if geographic module is available
    try:
        from . import geographic
    except ImportError:
        logger.error("Geographic module not available - cannot update geographic regions")
        return df, {'error': 'Geographic module not available'}

    # Check for required columns
    if 'lat' not in df.columns or 'lon' not in df.columns:
        logger.warning("No lat/lon columns found - skipping region assignment")
        return df, {'error': 'No lat/lon columns found'}

    # Filter to samples with coordinates
    has_coords = df['lat'].notna() & df['lon'].notna()
    n_with_coords = has_coords.sum()

    if n_with_coords == 0:
        logger.warning("No samples with coordinates - skipping region assignment")
        return df, {'error': 'No samples with coordinates'}

    logger.info(f"Updating {geo_category} for {n_with_coords} samples with coordinates")

    # Store old assignments (check for both the custom category name and 'ocean_basin')
    old_regions = None
    old_column = None
    if geo_category in df.columns:
        old_regions = df[geo_category].copy()
        old_column = geo_category
    elif 'ocean_basin' in df.columns:
        old_regions = df['ocean_basin'].copy()
        old_column = 'ocean_basin'

    # Update regions
    try:
        if shapefile_path:
            logger.info(f"Using custom shapefile: {shapefile_path}")
            logger.info(f"  Field: {shapefile_field}")
            logger.info(f"  Category: {geo_category}")
            # Load and use custom shapefile with custom field
            df_updated = geographic.assign_regions_from_shapefile(
                df,
                shapefile_path=shapefile_path,
                shapefile_field=shapefile_field,
                output_column=geo_category
            )
        else:
            logger.info("Using default ocean basin shapefile")
            df_updated = geographic.assign_ocean_basins(df)
            # Rename to custom category if specified and not 'ocean_basin'
            if geo_category != 'ocean_basin' and 'ocean_basin' in df_updated.columns:
                df_updated = df_updated.rename(columns={'ocean_basin': geo_category})
    except Exception as e:
        logger.error(f"Failed to update geographic regions: {e}")
        return df, {'error': str(e)}

    # Generate report
    if old_regions is not None:
        # Compare old vs new assignments
        changed = (old_regions != df_updated[geo_category]) & old_regions.notna() & df_updated[geo_category].notna()
        n_changed = changed.sum()

        changes_df = df_updated[changed][['processid', 'lat', 'lon']].copy()
        changes_df[f'old_{geo_category}'] = old_regions[changed]
        changes_df[f'new_{geo_category}'] = df_updated.loc[changed, geo_category]

        report = {
            'geo_category': geo_category,
            'shapefile': str(shapefile_path) if shapefile_path else 'default',
            'shapefile_field': shapefile_field if shapefile_path else 'name',
            'samples_with_coords': n_with_coords,
            'samples_assigned': df_updated[geo_category].notna().sum(),
            'assignments_changed': n_changed,
            'change_rate': n_changed / n_with_coords * 100 if n_with_coords > 0 else 0,
            'changes': changes_df.to_dict('records') if n_changed > 0 else []
        }
    else:
        report = {
            'geo_category': geo_category,
            'shapefile': str(shapefile_path) if shapefile_path else 'default',
            'shapefile_field': shapefile_field if shapefile_path else 'name',
            'samples_with_coords': n_with_coords,
            'samples_assigned': df_updated[geo_category].notna().sum(),
            'assignments_changed': 0,
            'change_rate': 0,
            'changes': []
        }

    logger.info(f"{geo_category} update complete: {report['samples_assigned']} assigned, {report['assignments_changed']} changed")

    return df_updated, report


def add_custom_grouping(
    df: pd.DataFrame,
    grouping_column: str,
    validate: bool = True
) -> Tuple[pd.DataFrame, Dict]:
    """
    Add a custom grouping variable for analysis.

    Args:
        df: DataFrame with enriched data
        grouping_column: Name of column to use for grouping
        validate: Validate grouping variable

    Returns:
        Tuple of (dataframe with grouping flag, report dict)
    """
    logger.info(f"Adding custom grouping: {grouping_column}")

    if grouping_column not in df.columns:
        raise ValueError(f"Grouping column '{grouping_column}' not found in data")

    # Get grouping statistics
    n_total = len(df)
    n_with_group = df[grouping_column].notna().sum()
    n_groups = df[grouping_column].nunique()
    group_sizes = df[grouping_column].value_counts()

    report = {
        'grouping_column': grouping_column,
        'total_samples': n_total,
        'samples_with_group': n_with_group,
        'samples_without_group': n_total - n_with_group,
        'n_groups': n_groups,
        'group_sizes': group_sizes.to_dict(),
        'min_group_size': int(group_sizes.min()) if len(group_sizes) > 0 else 0,
        'max_group_size': int(group_sizes.max()) if len(group_sizes) > 0 else 0,
        'mean_group_size': float(group_sizes.mean()) if len(group_sizes) > 0 else 0
    }

    logger.info(f"Grouping defined: {n_groups} groups, {n_with_group}/{n_total} samples assigned")

    # Add flag indicating this is a custom grouping
    df['custom_grouping'] = grouping_column

    return df, report


def generate_enrichment_report(
    output_dir: Path,
    merge_reports: List[Dict],
    region_report: Optional[Dict],
    grouping_report: Optional[Dict]
) -> Path:
    """
    Generate a text report summarizing enrichment operations.

    Args:
        output_dir: Output directory for report
        merge_reports: List of merge operation reports
        region_report: Geographic region update report
        grouping_report: Custom grouping report

    Returns:
        Path to report file
    """
    report_path = output_dir / "enrichment_report.txt"

    with open(report_path, 'w') as f:
        f.write("Metadata Enrichment Report\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")

        # Metadata merging
        if merge_reports:
            f.write("METADATA MERGING\n")
            f.write("-" * 80 + "\n")
            for i, report in enumerate(merge_reports, 1):
                f.write(f"\nMerge Operation {i}:\n")
                f.write(f"  Metadata File: {report['metadata_file']}\n")
                f.write(f"  Join Column: {report['join_column']}\n")
                f.write(f"  New Columns Added: {len(report['new_columns'])}\n")
                if report['new_columns']:
                    for col in report['new_columns']:
                        f.write(f"    - {col}\n")
                f.write(f"\n  Match Statistics:\n")
                f.write(f"    Samples in base: {report['samples_in_base']}\n")
                f.write(f"    Samples in metadata: {report['samples_in_metadata']}\n")
                f.write(f"    Samples matched: {report['samples_matched']}\n")
                f.write(f"    Samples unmatched: {report['samples_unmatched']}\n")
                f.write(f"    Match rate: {report['match_rate']:.1f}%\n")

                if report.get('overlapping_columns'):
                    f.write(f"\n  ⚠ Warning: Overlapping columns (suffixed with '_new'):\n")
                    for col in report['overlapping_columns']:
                        f.write(f"    - {col}\n")
            f.write("\n")

        # Geographic region updates
        if region_report and 'error' not in region_report:
            geo_cat = region_report.get('geo_category', 'geographic_region').upper().replace('_', ' ')
            f.write(f"{geo_cat} UPDATES\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Category: {region_report.get('geo_category', 'unknown')}\n")
            f.write(f"  Shapefile: {region_report['shapefile']}\n")
            if region_report.get('shapefile_field'):
                f.write(f"  Shapefile field: {region_report['shapefile_field']}\n")
            f.write(f"  Samples with coordinates: {region_report['samples_with_coords']}\n")
            f.write(f"  Samples assigned: {region_report['samples_assigned']}\n")
            f.write(f"  Assignments changed: {region_report['assignments_changed']}\n")
            if region_report['assignments_changed'] > 0:
                f.write(f"  Change rate: {region_report['change_rate']:.1f}%\n")
                f.write(f"\n  Changed Assignments:\n")
                # Get the column names dynamically
                geo_category = region_report.get('geo_category', 'region')
                old_col = f'old_{geo_category}'
                new_col = f'new_{geo_category}'
                for change in region_report['changes'][:10]:  # Show first 10
                    old_val = change.get(old_col, 'Unknown')
                    new_val = change.get(new_col, 'Unknown')
                    f.write(f"    {change['processid']}: {old_val} → {new_val}\n")
                if len(region_report['changes']) > 10:
                    f.write(f"    ... and {len(region_report['changes']) - 10} more\n")
            f.write("\n")
        elif region_report and 'error' in region_report:
            f.write("GEOGRAPHIC REGION UPDATES\n")
            f.write("-" * 80 + "\n")
            f.write(f"  ⚠ Error: {region_report['error']}\n\n")

        # Custom grouping
        if grouping_report:
            f.write("CUSTOM GROUPING\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Grouping Column: {grouping_report['grouping_column']}\n")
            f.write(f"  Total Samples: {grouping_report['total_samples']}\n")
            f.write(f"  Samples with Group: {grouping_report['samples_with_group']}\n")
            f.write(f"  Number of Groups: {grouping_report['n_groups']}\n")
            f.write(f"  Group Size Range: {grouping_report['min_group_size']}-{grouping_report['max_group_size']}\n")
            f.write(f"  Mean Group Size: {grouping_report['mean_group_size']:.1f}\n")

            f.write(f"\n  Group Distribution:\n")
            for group, size in sorted(grouping_report['group_sizes'].items(),
                                    key=lambda x: x[1], reverse=True)[:10]:
                f.write(f"    {group}: {size} samples\n")
            if len(grouping_report['group_sizes']) > 10:
                f.write(f"    ... and {len(grouping_report['group_sizes']) - 10} more groups\n")
            f.write("\n")

        # Summary
        f.write("SUMMARY\n")
        f.write("-" * 80 + "\n")
        total_operations = len(merge_reports)
        if region_report and 'error' not in region_report:
            total_operations += 1
        if grouping_report:
            total_operations += 1
        f.write(f"Total enrichment operations: {total_operations}\n")

        total_new_columns = sum(len(r['new_columns']) for r in merge_reports)
        f.write(f"Total new columns added: {total_new_columns}\n")

        f.write("\n")
        f.write("Enriched data saved to: {organism}_enriched.csv\n")
        f.write("\n")

    logger.info(f"Enrichment report written to {report_path}")
    return report_path


def regenerate_visualizations(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    grouping_column: Optional[str] = None,
    geo_category: Optional[str] = None
) -> Dict[str, Path]:
    """
    Regenerate visualizations with enriched data.

    Args:
        df: Enriched dataframe
        output_dir: Output directory for visualizations
        organism: Organism name
        grouping_column: Optional custom grouping variable
        geo_category: Optional geographic category name (e.g., 'ocean_basin', 'ecoregion')

    Returns:
        Dictionary of generated plot paths
    """
    logger.info("Regenerating visualizations with enriched data")

    viz_dir = output_dir / "updated_visualizations"
    viz_dir.mkdir(parents=True, exist_ok=True)

    plots = {}

    # Try to import visualization module
    try:
        from . import visualization
    except ImportError:
        logger.warning("Visualization module not available - skipping plot regeneration")
        return plots

    # Distribution map (if geographic data available)
    if 'lat' in df.columns and 'lon' in df.columns:
        try:
            map_path = viz_dir / f"{organism}_distribution_map_enriched.pdf"
            visualization.create_distribution_map(
                df,
                output_path=map_path,
                title=f"{organism} Distribution (Enriched Data)"
            )
            plots['distribution_map'] = map_path
            logger.info(f"Generated distribution map: {map_path}")
        except Exception as e:
            logger.warning(f"Failed to generate distribution map: {e}")

    # Distribution by custom grouping (if specified)
    if grouping_column and grouping_column in df.columns:
        try:
            bar_path = viz_dir / f"{organism}_distribution_by_{grouping_column}.pdf"
            visualization.create_distribution_barplot(
                df,
                grouping_var=grouping_column,
                output_path=bar_path,
                title=f"{organism} Distribution by {grouping_column}"
            )
            plots[f'distribution_by_{grouping_column}'] = bar_path
            logger.info(f"Generated distribution barplot: {bar_path}")
        except Exception as e:
            logger.warning(f"Failed to generate distribution barplot: {e}")

    # Distribution by geographic region (detect dynamically or use specified category)
    # Auto-detect geographic category column if not specified
    if geo_category is None:
        # Try common geographic category names
        for potential_col in ['ocean_basin', 'ecoregion', 'watershed', 'biome', 'freshwater_basin', 'region']:
            if potential_col in df.columns:
                geo_category = potential_col
                break

    if geo_category and geo_category in df.columns:
        try:
            # Create human-readable label
            geo_label = geo_category.replace('_', ' ').title()
            safe_name = geo_category.replace('_', '-')

            region_path = viz_dir / f"{organism}_distribution_by_{safe_name}_enriched.pdf"
            visualization.create_distribution_barplot(
                df,
                grouping_var=geo_category,
                output_path=region_path,
                title=f"{organism} Distribution by {geo_label} (Enriched Data)"
            )
            plots[f'distribution_by_{safe_name}'] = region_path
            logger.info(f"Generated {geo_label} distribution plot: {region_path}")
        except Exception as e:
            logger.warning(f"Failed to generate {geo_category} distribution plot: {e}")

    logger.info(f"Generated {len(plots)} visualization(s)")

    return plots


def enrich_metadata(
    input_csv: Path,
    output_dir: Path,
    metadata_files: Optional[List[Path]] = None,
    join_column: str = 'processid',
    shapefile_path: Optional[Path] = None,
    shapefile_field: str = 'name',
    geo_category: str = 'ocean_basin',
    grouping_column: Optional[str] = None,
    recalculate_geography: bool = False,
    organism: Optional[str] = None
) -> Dict:
    """
    Main function to enrich metadata in BOLDGenotyper output.

    Args:
        input_csv: Path to input annotated CSV
        output_dir: Output directory for enriched data
        metadata_files: List of metadata CSV files to merge
        join_column: Column to join on
        shapefile_path: Optional path to custom shapefile
        shapefile_field: Shapefile attribute containing region names (default: 'name')
        geo_category: Geographic category name (default: 'ocean_basin')
        grouping_column: Optional custom grouping variable
        recalculate_geography: Whether to recalculate geographic summaries
        organism: Organism name (extracted from filename if not provided)

    Returns:
        Dictionary with paths to output files and summary statistics
    """
    logger.info("Starting metadata enrichment")

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load base data
    df = load_annotated_csv(input_csv)

    # Extract organism name if not provided
    if organism is None:
        organism = input_csv.stem.replace('_annotated', '')
        logger.info(f"Extracted organism name: {organism}")

    # Track all operations
    merge_reports = []
    region_report = None
    grouping_report = None

    # Merge metadata files
    if metadata_files:
        for metadata_file in metadata_files:
            try:
                df, report = merge_metadata(df, metadata_file, join_column)
                merge_reports.append(report)
            except Exception as e:
                logger.error(f"Failed to merge {metadata_file}: {e}")
                merge_reports.append({
                    'metadata_file': str(metadata_file),
                    'error': str(e)
                })

    # Update geographic regions
    if recalculate_geography or shapefile_path:
        try:
            df, region_report = update_geographic_regions(
                df,
                shapefile_path=shapefile_path,
                shapefile_field=shapefile_field,
                geo_category=geo_category,
                force_recalculate=recalculate_geography
            )

            # Save geographic updates if there were changes
            if region_report and region_report.get('changes'):
                changes_df = pd.DataFrame(region_report['changes'])
                changes_path = output_dir / "geographic_updates.csv"
                changes_df.to_csv(changes_path, index=False)
                logger.info(f"Saved geographic updates to {changes_path}")
        except Exception as e:
            logger.error(f"Failed to update geographic regions: {e}")
            region_report = {'error': str(e)}

    # Add custom grouping
    if grouping_column:
        try:
            df, grouping_report = add_custom_grouping(df, grouping_column)
        except Exception as e:
            logger.error(f"Failed to add custom grouping: {e}")
            grouping_report = {'error': str(e)}

    # Save enriched data
    enriched_path = output_dir / f"{organism}_enriched.csv"
    df.to_csv(enriched_path, index=False)
    logger.info(f"Saved enriched data to {enriched_path}")

    # Generate enrichment report
    report_path = generate_enrichment_report(
        output_dir,
        merge_reports,
        region_report,
        grouping_report
    )

    # Regenerate visualizations
    plots = regenerate_visualizations(
        df,
        output_dir,
        organism,
        grouping_column,
        geo_category
    )

    # Return summary
    result = {
        'enriched_csv': enriched_path,
        'report': report_path,
        'n_samples': len(df),
        'n_columns': len(df.columns),
        'merge_operations': len(merge_reports),
        'region_updated': region_report is not None and 'error' not in region_report,
        'grouping_added': grouping_report is not None and 'error' not in grouping_report,
        'plots_generated': len(plots),
        'plot_paths': plots
    }

    logger.info("Metadata enrichment complete")
    logger.info(f"  Enriched CSV: {enriched_path}")
    logger.info(f"  Report: {report_path}")
    logger.info(f"  Operations: {result['merge_operations']} merges, "
                f"region_updated={result['region_updated']}, "
                f"grouping_added={result['grouping_added']}")
    logger.info(f"  Visualizations: {result['plots_generated']} plot(s) generated")

    return result
