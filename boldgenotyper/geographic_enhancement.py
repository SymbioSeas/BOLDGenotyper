"""
Geographic Analysis Enhancement Module for BOLDGenotyper

This module enhances geographic analysis by providing better transparency about
data coverage, quality, and limitations. It helps users understand what geographic
interpretations are well-supported by their data.

Works with any geographic category (ocean basins, freshwater basins, ecoregions,
watersheds, etc.) via the ``geo_category`` parameter.

Key Features:
- Geographic coverage assessment per haplotype
- Region assignment confidence levels
- Missing data transparency reports
- Coordinate quality evaluation
- Geographic limitations documentation

Author: Steph Smith (symbioseas@outlook.com)
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class GeographicEnhancementError(Exception):
    """Base exception for geographic enhancement errors."""
    pass


def assess_geographic_coverage(
    df: pd.DataFrame,
    group_col: str = "haplotype_id",
    species_col: str = "species",
    geo_category: str = "ocean_basin"
) -> pd.DataFrame:
    """
    Assess geographic data coverage per haplotype.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with geographic data
    group_col : str
        Column name for haplotype groups
    species_col : str
        Column name for species
    geo_category : str
        Geographic category column name (default: 'ocean_basin')

    Returns
    -------
    pd.DataFrame
        Coverage statistics per haplotype
    """
    logger.info("  Assessing geographic coverage...")

    coverage_data = []

    for group in df[group_col].unique():
        if pd.isna(group):
            continue

        group_df = df[df[group_col] == group]
        n_total = len(group_df)

        # Get species
        species = group_df[species_col].mode()[0] if species_col in group_df.columns else "Unknown"

        # Count samples with different levels of geographic data
        has_lat_lon = 0
        has_region = 0
        has_country = 0

        if 'lat' in group_df.columns and 'lon' in group_df.columns:
            has_lat_lon = group_df[['lat', 'lon']].notna().all(axis=1).sum()

        if geo_category in group_df.columns:
            has_region = group_df[geo_category].notna().sum()
            # Exclude "Unknown" regions
            has_region = (
                (group_df[geo_category].notna()) &
                (group_df[geo_category] != 'Unknown') &
                (group_df[geo_category] != '')
            ).sum()

        if 'country' in group_df.columns:
            has_country = group_df['country'].notna().sum()

        # Calculate percentages
        pct_lat_lon = (has_lat_lon / n_total * 100) if n_total > 0 else 0
        pct_region = (has_region / n_total * 100) if n_total > 0 else 0
        pct_country = (has_country / n_total * 100) if n_total > 0 else 0

        # Assess representativeness
        if pct_region > 50:
            representativeness = 'excellent'
        elif pct_region >= 25:
            representativeness = 'good'
        elif pct_region >= 10 or pct_lat_lon >= 50:
            representativeness = 'moderate'
        elif pct_lat_lon >= 10:
            representativeness = 'poor'
        else:
            representativeness = 'very_poor'

        geo_label = geo_category.replace('_', ' ')
        coverage_data.append({
            'genotype': group,
            'species': species,
            'total_samples': n_total,
            'with_lat_lon': has_lat_lon,
            'pct_with_coords': pct_lat_lon,
            f'with_{geo_category}': has_region,
            'pct_with_region': pct_region,
            'with_country': has_country,
            'pct_with_country': pct_country,
            'representativeness': representativeness
        })

    return pd.DataFrame(coverage_data)


def calculate_basin_assignment_confidence(
    df: pd.DataFrame,
    goas_data: Optional[Any] = None,
    geo_category: str = "ocean_basin"
) -> pd.DataFrame:
    """
    Calculate confidence levels for geographic region assignments.

    Confidence is based on distance to region boundaries. Samples near
    boundaries have lower confidence.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with lat, lon, and geographic region columns
    goas_data : GeoDataFrame, optional
        Shapefile data for boundary calculations
    geo_category : str
        Geographic category column name (default: 'ocean_basin')

    Returns
    -------
    pd.DataFrame
        Enhanced dataframe with basin_confidence column
    """
    logger.info("  Calculating basin assignment confidence...")

    df = df.copy()

    # Initialize confidence column
    df['basin_confidence'] = 'none'
    df['distance_to_boundary_km'] = np.nan
    df['basin_notes'] = ''

    # For samples with coordinates and basin assignments
    has_coords = df['lat'].notna() & df['lon'].notna()
    has_basin = df[geo_category].notna() & (df[geo_category] != 'Unknown') if geo_category in df.columns else pd.Series([False] * len(df), index=df.index)

    valid_samples = has_coords & has_basin

    if not valid_samples.any():
        return df

    # If GOaS data available, calculate distances
    if goas_data is not None:
        try:
            import geopandas as gpd
            from shapely.geometry import Point

            # Convert samples to GeoDataFrame
            geometry = [Point(lon, lat) for lat, lon in
                       zip(df.loc[valid_samples, 'lat'],
                           df.loc[valid_samples, 'lon'])]

            samples_gdf = gpd.GeoDataFrame(
                df.loc[valid_samples],
                geometry=geometry,
                crs="EPSG:4326"
            )

            # Pre-compute boundary geometries per region (avoids redundant work)
            name_col = 'name' if 'name' in goas_data.columns else goas_data.columns[0]
            region_boundaries = {}
            for region_name in samples_gdf[geo_category].unique():
                basin_polygons = goas_data[goas_data[name_col] == region_name]
                if len(basin_polygons) > 0:
                    region_boundaries[region_name] = basin_polygons.unary_union.boundary

            # Vectorized distance calculation per region
            distances = pd.Series(np.nan, index=samples_gdf.index)
            confidences = pd.Series('unknown', index=samples_gdf.index)
            notes = pd.Series('', index=samples_gdf.index)

            for region_name, boundary in region_boundaries.items():
                mask = samples_gdf[geo_category] == region_name
                if not mask.any():
                    continue
                region_points = samples_gdf.loc[mask, 'geometry']
                dist_deg = region_points.distance(boundary)
                dist_km = dist_deg * 111  # Rough approximation

                distances.loc[mask] = dist_km
                confidences.loc[mask] = np.where(
                    dist_km > 50, 'high',
                    np.where(dist_km > 10, 'medium', 'low')
                )
                notes.loc[mask] = np.where(dist_km <= 10, 'near basin boundary', '')

            df.loc[valid_samples, 'distance_to_boundary_km'] = distances
            df.loc[valid_samples, 'basin_confidence'] = confidences
            df.loc[valid_samples, 'basin_notes'] = notes

        except Exception as e:
            logger.warning(f"Could not calculate basin distances: {e}")
            # Fallback: assign medium confidence to all
            df.loc[valid_samples, 'basin_confidence'] = 'medium'

    else:
        # Without GOaS data, assign medium confidence to all assigned basins
        df.loc[valid_samples, 'basin_confidence'] = 'medium'

    # Samples without coordinates
    no_coords = ~has_coords
    df.loc[no_coords, 'basin_confidence'] = 'none'
    df.loc[no_coords, 'basin_notes'] = 'no coordinates available'

    return df


def generate_missing_data_report(
    df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    output_path: Path,
    organism: str,
    geo_category: str = "ocean_basin"
) -> None:
    """
    Generate comprehensive missing data report.

    Parameters
    ----------
    df : pd.DataFrame
        Full annotated dataframe
    coverage_df : pd.DataFrame
        Coverage statistics per genotype
    output_path : Path
        Output path for report
    organism : str
        Organism name
    geo_category : str
        Geographic category name (default: 'ocean_basin')
    """
    logger.info("  Generating missing data report...")

    geo_label = geo_category.replace('_', ' ').title()
    region_col = f'with_{geo_category}'

    # Overall statistics
    n_total = len(df)
    n_with_coords = df[['lat', 'lon']].notna().all(axis=1).sum()
    n_with_region = ((df[geo_category].notna()) &
                    (df[geo_category] != 'Unknown') &
                    (df[geo_category] != '')).sum() if geo_category in df.columns else 0
    n_with_country = df['country'].notna().sum() if 'country' in df.columns else 0
    n_no_geo = ((~df[['lat', 'lon']].notna().all(axis=1)) &
                (~df[geo_category].notna()) &
                (~df['country'].notna())).sum() if geo_category in df.columns and 'country' in df.columns else n_total - n_with_coords

    pct_coords = (n_with_coords / n_total * 100) if n_total > 0 else 0
    pct_region = (n_with_region / n_total * 100) if n_total > 0 else 0
    pct_country = (n_with_country / n_total * 100) if n_total > 0 else 0
    pct_no_geo = (n_no_geo / n_total * 100) if n_total > 0 else 0

    # Generate report text
    report = f"""Geographic Data Availability Report
====================================

Dataset: {organism}
Total Samples: {n_total}
Geographic Category: {geo_label}

Overall Coverage:
-----------------
With Coordinates: {n_with_coords} ({pct_coords:.1f}%)
With {geo_label}: {n_with_region} ({pct_region:.1f}%)
With Country Only: {n_with_country} ({pct_country:.1f}%)
No Geographic Data: {n_no_geo} ({pct_no_geo:.1f}%)

By Genotype:
------------
"""

    # Add per-genotype statistics
    for _, row in coverage_df.iterrows():
        genotype = row['genotype']
        n_samples = row['total_samples']
        n_coords = row['with_lat_lon']
        n_region = row.get(region_col, row.get('with_ocean_basin', 0))
        n_country = row['with_country']
        representativeness = row['representativeness']

        pct_with_region = row.get('pct_with_region', row.get('pct_with_basin', 0))

        report += f"\n{genotype} (n={n_samples}):\n"
        if n_coords > 0:
            report += f"  + Coordinates: {n_coords} ({row['pct_with_coords']:.1f}%)\n"
        else:
            report += f"  - Coordinates: 0 (0.0%)\n"

        if n_region > 0:
            report += f"  + {geo_label}: {n_region} ({pct_with_region:.1f}%)\n"
        else:
            report += f"  - {geo_label}: 0 (0.0%)\n"

        if n_country > 0:
            report += f"  ~ Country Only: {n_country} ({row['pct_with_country']:.1f}%)\n"
        else:
            report += f"  - Country Only: 0 (0.0%)\n"

        report += f"  -> Representativeness: {representativeness.upper()}\n"

    # Add limitations section
    report += f"""
Limitations:
------------
"""

    if pct_region < 25:
        report += f"1. {geo_label} assignments available for only {pct_region:.1f}% of samples\n"

    if pct_coords < 50:
        report += f"2. Coordinate-based visualizations rely on {pct_coords:.1f}% of samples\n"

    if pct_country < 80 and pct_region < 25:
        report += "3. Country-level data insufficient for multi-region countries\n"

    if pct_region < 25 or pct_coords < 50:
        report += "4. Geographic interpretations should be considered preliminary\n"

    # Add recommendations
    report += """
Recommendations:
----------------
"""

    if pct_region < 25:
        report += f"1. Present {geo_label.lower()} data as supplemental (not primary findings)\n"

    if pct_coords >= 25:
        report += "2. Rely primarily on coordinate-based maps for geographic patterns\n"

    if pct_region < 50 or pct_coords < 50:
        report += "3. Acknowledge coverage limitations in methods/discussion\n"
        report += "4. Consider contacting depositors for missing coordinates\n"

    report += "5. Future sampling should prioritize GPS precision\n"

    # Data quality notes
    report += """
Data Quality Notes:
-------------------
"""

    if 'basin_confidence' in df.columns:
        confidence_counts = df['basin_confidence'].value_counts()
        report += "Region Assignment Confidence:\n"
        for conf, count in confidence_counts.items():
            pct = (count / n_total * 100) if n_total > 0 else 0
            report += f"  - {conf}: {count} ({pct:.1f}%)\n"

    # Save report
    with open(output_path, 'w') as f:
        f.write(report)

    logger.info(f"  ✓ Missing data report saved to {output_path}")


def print_geographic_warnings(
    coverage_df: pd.DataFrame,
    overall_pct_basin: float,
    overall_pct_coords: float,
    geo_category: str = "ocean_basin"
) -> None:
    """
    Print geographic coverage warnings to console.

    Parameters
    ----------
    coverage_df : pd.DataFrame
        Coverage statistics
    overall_pct_basin : float
        Overall percentage with geographic regions
    overall_pct_coords : float
        Overall percentage with coordinates
    geo_category : str
        Geographic category name (default: 'ocean_basin')
    """
    if overall_pct_basin < 25 or overall_pct_coords < 50:
        print("\n" + "="*70)
        print("⚠️  GEOGRAPHIC DATA WARNING")
        print("="*70)

        if overall_pct_basin < 25:
            print(f"\nOnly {overall_pct_basin:.1f}% of samples have {geo_category.replace('_', ' ')} assignments")

        if overall_pct_coords < 50:
            print(f"Only {overall_pct_coords:.1f}% have precise coordinates")

        print("\nThis may limit biogeographic conclusions.")
        print("Consider:")
        print("1. Focusing interpretation on samples with coordinates")
        print(f"2. Presenting {geo_category.replace('_', ' ')} analysis as preliminary")
        print("3. Acknowledging limitations in manuscript")

        # Identify problematic genotypes
        poor_genotypes = coverage_df[
            coverage_df['representativeness'].isin(['poor', 'very_poor'])
        ]

        if len(poor_genotypes) > 0:
            print(f"\n⚠️  {len(poor_genotypes)} genotype(s) with poor geographic coverage:")
            for _, row in poor_genotypes.head(5).iterrows():
                print(f"  - {row['genotype']}: {row['pct_with_region']:.1f}% with {geo_category.replace('_', ' ')}")

        print("\nSee: geographic_analysis/geographic_coverage.csv for per-haplotype details")
        print("="*70 + "\n")


def enhance_geographic_analysis(
    df: pd.DataFrame,
    output_dir: Path,
    organism: str,
    group_col: str = "haplotype_id",
    goas_data: Optional[Any] = None,
    geo_category: str = "ocean_basin"
) -> Dict[str, Any]:
    """
    Complete geographic enhancement workflow.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe
    output_dir : Path
        Output directory
    organism : str
        Organism name
    group_col : str
        Column name for haplotype groups
    goas_data : GeoDataFrame, optional
        GOaS data for basin confidence
    geo_category : str
        Geographic category name (default: 'ocean_basin')

    Returns
    -------
    dict
        Results and file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Enhancing geographic analysis...")

    results = {}

    # 1. Assess coverage
    coverage_df = assess_geographic_coverage(df, group_col=group_col, geo_category=geo_category)
    coverage_path = output_dir / "geographic_coverage.csv"
    coverage_df.to_csv(coverage_path, index=False)
    results['coverage'] = coverage_path
    logger.info(f"  ✓ Coverage assessment saved: {coverage_path}")

    # 2. Calculate basin confidence
    df_enhanced = calculate_basin_assignment_confidence(df, goas_data, geo_category=geo_category)
    results['enhanced_df'] = df_enhanced

    # Calculate overall statistics
    n_total = len(df)
    n_with_basin = ((df[geo_category].notna()) &
                    (df[geo_category] != 'Unknown') &
                    (df[geo_category] != '')).sum() if geo_category in df.columns else 0
    n_with_coords = df[['lat', 'lon']].notna().all(axis=1).sum()

    overall_pct_basin = (n_with_basin / n_total * 100) if n_total > 0 else 0
    overall_pct_coords = (n_with_coords / n_total * 100) if n_total > 0 else 0

    results['overall_stats'] = {
        'pct_region': overall_pct_basin,
        'pct_coords': overall_pct_coords
    }

    # 3. Print warnings
    print_geographic_warnings(coverage_df, overall_pct_basin, overall_pct_coords, geo_category=geo_category)

    # 4. Create README
    readme_path = output_dir / "README.md"
    readme_content = f"""# Geographic Analysis for {organism}

## Overview

Enhanced geographic analysis with coverage assessment and quality metrics.

## Files

- `geographic_coverage.csv` - Coverage statistics per haplotype
- Enhanced annotated CSV includes `basin_confidence` column

## Coverage Summary

- **Geographic Region**: {overall_pct_basin:.1f}% of samples with region assignment
- **Coordinates**: {overall_pct_coords:.1f}% of samples with lat/lon

## Confidence Levels

**Region Assignment Confidence:**
- `high`: >50 km from region boundary
- `medium`: 10-50 km from boundary
- `low`: <10 km from boundary (uncertain assignment)
- `none`: No coordinates available

## Representativeness Scores

**Per-Haplotype Geographic Coverage:**
- `excellent`: >50% with region assignment
- `good`: 25-50% with region assignment
- `moderate`: 10-25% with region OR >50% with coords
- `poor`: <10% with region AND <25% with coords
- `very_poor`: <10% with any geographic data

## Interpretation

Geographic patterns should be interpreted cautiously when coverage is
poor or moderate. High-confidence region assignments are preferred for
biogeographic conclusions.
"""

    with open(readme_path, 'w') as f:
        f.write(readme_content)

    results['readme'] = readme_path

    logger.info("  ✓ Geographic enhancement complete")

    return results
