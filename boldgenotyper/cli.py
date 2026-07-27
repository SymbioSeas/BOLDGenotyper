#!/usr/bin/env python3
"""
BOLDGenotyper Command-Line Interface

Unified pipeline for BOLD barcode data implementing haplotype-first workflow:
COI validation, haplotype discovery, assignment, geographic analysis,
phylogenetic reconstruction, and visualization.

Workflow:
1. Data Loading & Coordinate Quality
2. Pre-processing QC (Orientation normalization, ORF validation, Dynamic filtering)
3. Haplotype Discovery (Align, Extract core region, Identify ESVs, Flag suspects)
4. Haplotype Assignment (Match samples to haplotypes)
5. Taxonomy Assignment
6. Post-assignment QC (Contamination detection)
6.5. Species-Level Aggregation (Group by species, diversity metrics)
7. Geographic Enhancement
8. Phylogenetic Analysis (Optional)
9. Divergence Analysis (Haplotype-level, species-level, within-species)
9.5. Metadata Analysis (Optional)
10. Visualization (Haplotype, species-level, and species-faceted plots)
10.5. Plot Data Export (Optional)
11. Population Genetics Export (Optional)
12. Reports

Author: Steph Smith (symbioseas@outlook.com)
"""

import argparse
import sys
import os
import logging
import json
import shutil
from pathlib import Path
from typing import Optional
import pandas as pd
from Bio import SeqIO

# Local imports
from . import (
    utils, config, metadata, geographic, dereplication,
    haplotype_assignment, phylogenetics, visualization, reports,
    quality_control, plot_export,
    divergence_analysis, parameter_sweep, geographic_enhancement,
    popgen_export, species_analysis, msa_visualization,
    metadata_summary,
)
from . import __version__

logger = logging.getLogger(__name__)


def extract_organism_from_path(path: Path) -> str:
    """
    Extract organism name from TSV filename and sanitize for consistent file naming.

    This function extracts the organism name from common BOLD filename patterns
    and ensures the result uses underscores instead of spaces for consistent
    output file naming across all modules.

    Parameters
    ----------
    path : Path
        Path to input TSV file

    Returns
    -------
    str
        Sanitized organism name with underscores (no spaces)

    Examples
    --------
    >>> extract_organism_from_path(Path("Sphyrna_lewini.tsv"))
    'Sphyrna_lewini'
    >>> extract_organism_from_path(Path("Great White Shark.tsv"))
    'Great_White_Shark'
    """
    # Remove extension and path
    name = path.stem
    # Common patterns: "Genus_species.tsv", "Genus.tsv", "genus_species_data.tsv"
    # Take first part before underscore or use whole name if no underscore
    parts = name.split('_')
    if len(parts) >= 2 and parts[1] not in ['data', 'bold', 'samples', 'sequences']:
        # Likely "Genus_species" format
        result = '_'.join(parts[:2])
    else:
        result = parts[0].capitalize()

    # Always sanitize to ensure spaces are replaced with underscores
    return utils.sanitize_filename(result)


def remove_empty_directories(base_path: Path) -> None:
    """
    Remove empty directories recursively, starting from leaf directories.

    Parameters
    ----------
    base_path : Path
        Base directory to scan for empty subdirectories
    """
    # Walk bottom-up to remove leaf directories first
    for dirpath, dirnames, filenames in os.walk(base_path, topdown=False):
        dir_path = Path(dirpath)

        # Skip the base directory itself
        if dir_path == base_path:
            continue

        # Check if directory is empty (no files and no subdirectories)
        try:
            if not any(dir_path.iterdir()):
                dir_path.rmdir()
                logger.debug(f"Removed empty directory: {dir_path}")
        except OSError:
            # Directory not empty or can't be removed
            pass


def setup_directories(base_output: Path) -> dict:
    """Create organized output directory structure for haplotype-first workflow."""
    dirs = {
        'base': base_output,
        'intermediate': base_output / 'intermediate',
        'qc': base_output / 'intermediate' / 'quality_control',
        'haplotype_discovery': base_output / 'intermediate' / 'haplotype_discovery',
        'intermediate_phylo': base_output / 'intermediate' / 'phylogenetic',
        'intermediate_assignments': base_output / 'intermediate' / 'haplotype_assignments',
        'intermediate_geographic': base_output / 'intermediate' / 'geographic',
        'haplotypes': base_output / 'haplotypes',
        'assignments': base_output / 'haplotype_assignments',
        'taxonomy': base_output / 'taxonomy',
        'species_analysis': base_output / 'species_analysis',
        'quality_control': base_output / 'quality_control',
        'divergence_analysis': base_output / 'divergence_analysis',
        'geographic_analysis': base_output / 'geographic_analysis',
        'phylogenetic': base_output / 'phylogenetic',
        'visualization': base_output / 'visualization',
        'visualization_pdf': base_output / 'visualization' / 'pdf',
        'visualization_svg': base_output / 'visualization' / 'svg',
        'visualization_json': base_output / 'visualization' / 'json',
        'haplotype_plots': base_output / 'visualization' / 'haplotype_plots',
        'metadata_viz_pdf': base_output / 'visualization' / 'metadata' / 'pdf',
        'metadata_viz_svg': base_output / 'visualization' / 'metadata' / 'svg',
        'metadata_viz_png': base_output / 'visualization' / 'metadata' / 'png',
        'metadata_viz_json': base_output / 'visualization' / 'metadata' / 'json',
        'reports': base_output / 'reports',
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def get_viz_path(dirs: dict, filename: str, fmt: str, is_metadata: bool = False) -> Path:
    """
    Get the appropriate visualization output path based on format.

    Parameters
    ----------
    dirs : dict
        Directory dictionary from setup_directories
    filename : str
        Base filename without extension
    fmt : str
        Format (pdf, svg, png, etc.)
    is_metadata : bool
        If True, use metadata visualization subdirectory

    Returns
    -------
    Path
        Full path for the output file
    """
    if is_metadata:
        if fmt == 'pdf':
            return dirs['metadata_viz_pdf'] / f"{filename}.{fmt}"
        elif fmt == 'svg':
            return dirs['metadata_viz_svg'] / f"{filename}.{fmt}"
        elif fmt == 'png':
            return dirs['metadata_viz_png'] / f"{filename}.{fmt}"
        else:
            return dirs['metadata_viz_pdf'] / f"{filename}.{fmt}"
    else:
        if fmt == 'pdf':
            return dirs['visualization_pdf'] / f"{filename}.{fmt}"
        elif fmt == 'svg':
            return dirs['visualization_svg'] / f"{filename}.{fmt}"
        else:
            # PNG and other formats go to pdf directory
            return dirs['visualization_pdf'] / f"{filename}.{fmt}"


def run_pipeline(
    tsv_path: Path,
    organism: str,
    output_dir: Path,
    cfg: config.PipelineConfig,
    no_report: bool = False,
    skip_geo: bool = False,
    export_plot_data: bool = True,
    export_popgen_formats: Optional[list] = None,
    metadata_summary_enabled: bool = False,
    metadata_fields: Optional[list] = None,
    normalize_sex: bool = False,
    temporal_analysis: bool = False,
    shapefile_path: Optional[Path] = None,
    shapefile_field: str = 'name',
    geo_category: str = 'ocean_basin',
) -> bool:
    """
    Run the complete BOLDGenotyper pipeline.

    Parameters
    ----------
    tsv_path : Path
        Path to input BOLD TSV file
    organism : str
        Organism name
    output_dir : Path
        Output directory
    cfg : config.PipelineConfig
        Pipeline configuration
    no_report : bool, optional
        Skip HTML report generation (default: False)
    skip_geo : bool, optional
        Skip geographic analysis (default: False)
    export_plot_data : bool, optional
        Export raw plot data and Python regeneration scripts (default: True)
    export_popgen_formats : list, optional
        List of population genetics formats to export (default: None)
    shapefile_path : Path, optional
        Path to custom shapefile for geographic region assignment. Works with
        any polygon shapefile. If None, uses built-in GOaS ocean basins.
    shapefile_field : str, optional
        Shapefile attribute field containing region names (default: 'name')
    geo_category : str, optional
        Geographic category name for outputs. Defaults to 'ocean_basin' when
        using GOaS, or 'geographic_region' when using a custom shapefile.

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    logger.info("=" * 80)
    logger.info(f"BOLDGenotyper Pipeline - {organism}")
    logger.info("=" * 80)
    logger.info(f"Input TSV: {tsv_path}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("")

    # Setup directory structure
    dirs = setup_directories(output_dir)

    # Save pipeline parameters for reference and HTML report
    params = {
        'workflow': 'haplotype-first',
        'coi_validation': {
            'genetic_code': cfg.coi_validation.mitochondrial_code,
            'orf_min_coverage': cfg.coi_validation.orf_min_coverage,
            'orf_max_internal_stops': cfg.coi_validation.orf_max_internal_stops,
        },
        'qc': {
            'min_raw_length_abs': cfg.qc.min_raw_length_abs,
            'min_raw_length_frac_of_median': cfg.qc.min_raw_length_frac_of_median,
            'max_raw_N_fraction': cfg.qc.max_raw_N_fraction,
        },
        'core_region': {
            'min_coverage': cfg.core_region.core_min_coverage,
            'mask_gap_threshold': cfg.core_region.mask_gap_threshold,
            'min_length': cfg.core_region.core_min_length,
        },
        'haplotype': {
            'max_singleton_distance': cfg.haplotype.max_singleton_distance,
            'flag_suspect_haplotypes': cfg.haplotype.flag_suspect_haplotypes,
        },
        'assignment': {
            'min_identity': cfg.genotype_assignment.min_identity,
            'tie_margin': cfg.genotype_assignment.tie_margin,
        },
        'threads': cfg.n_threads,
        'build_tree': cfg.phylogenetic.build_tree,
        'geo_category': geo_category,
        'custom_shapefile': str(shapefile_path) if shapefile_path else None,
    }
    params_file = output_dir / f"{organism}_pipeline_parameters.json"
    with open(params_file, 'w') as f:
        json.dump(params, f, indent=2)
    logger.info(f"Saved pipeline parameters to {params_file}")

    # ========================================================================
    # PHASE 1: Data Loading and Quality Control
    # ========================================================================
    logger.info("PHASE 1: Data Loading and Quality Control")
    logger.info("-" * 80)

    try:
        # Parse BOLD TSV
        logger.info("1.1: Parsing BOLD TSV metadata...")
        df = metadata.parse_bold_tsv(tsv_path)
        logger.info(f"  ✓ Loaded {len(df)} samples with {len(df.columns)} columns")

        # Save parsed metadata
        parsed_tsv = dirs['intermediate'] / "01_parsed_metadata.tsv"
        df.to_csv(parsed_tsv, sep='\t', index=False)

        # Mark coordinate quality (don't filter - keep all samples for genotyping)
        logger.info("1.2: Marking coordinate quality...")
        df = metadata.mark_coordinate_quality(df)

        # Save metadata with quality markers
        marked_tsv = dirs['intermediate'] / "02_marked_metadata.tsv"
        df.to_csv(marked_tsv, sep='\t', index=False)

        # Assign geographic regions (only to samples with geographic quality)
        logger.info(f"1.3: Assigning geographic regions ({geo_category})...")
        geo_analysis_performed = False  # Track whether geographic analysis was successful

        if skip_geo:
            logger.info("  ⊘ Geographic analysis disabled (--no-geo flag)")
            df_with_basins = df.copy()
            df_with_basins[geo_category] = 'Unknown'
        elif shapefile_path:
            # Use custom shapefile
            logger.info(f"  → Using custom shapefile: {shapefile_path}")
            logger.info(f"  → Shapefile field: {shapefile_field}")
            logger.info(f"  → Geographic category: {geo_category}")

            try:
                # Only assign regions to samples with geographic quality coordinates
                df_geo_quality = df[df['is_geographic_quality']].copy()
                logger.info(f"  → {len(df_geo_quality)}/{len(df)} samples have geographic-quality coordinates")

                if len(df_geo_quality) > 0:
                    df_geo_assigned = geographic.assign_regions_from_shapefile(
                        df_geo_quality,
                        shapefile_path=shapefile_path,
                        shapefile_field=shapefile_field,
                        output_column=geo_category,
                        coord_col="coord"
                    )

                    # Merge back with full dataset, preserving lat/lon columns for visualization
                    df_with_basins = df.copy()
                    # First set all to Unknown
                    df_with_basins[geo_category] = 'Unknown'

                    # Initialize lat/lon columns if not present
                    if 'lat' not in df_with_basins.columns:
                        df_with_basins['lat'] = pd.NA
                    if 'lon' not in df_with_basins.columns:
                        df_with_basins['lon'] = pd.NA

                    # Extract columns to merge: geo_category and lat/lon (if they exist in df_geo_assigned)
                    cols_to_merge = ['processid', geo_category]
                    if 'lat' in df_geo_assigned.columns:
                        cols_to_merge.append('lat')
                    if 'lon' in df_geo_assigned.columns:
                        cols_to_merge.append('lon')

                    # Merge the geographic data back
                    merged_geo = df_with_basins[['processid']].merge(
                        df_geo_assigned[cols_to_merge], on='processid', how='left', suffixes=('', '_new')
                    )

                    # Update the columns
                    df_with_basins[geo_category] = merged_geo[geo_category].fillna('Unknown')
                    if 'lat_new' in merged_geo.columns:
                        df_with_basins['lat'] = merged_geo['lat_new']
                    elif 'lat' in merged_geo.columns:
                        df_with_basins['lat'] = merged_geo['lat']
                    if 'lon_new' in merged_geo.columns:
                        df_with_basins['lon'] = merged_geo['lon_new']
                    elif 'lon' in merged_geo.columns:
                        df_with_basins['lon'] = merged_geo['lon']

                    n_assigned = (df_with_basins[geo_category] != 'Unknown').sum()
                    logger.info(f"  ✓ Assigned {geo_category} to {n_assigned} samples")
                    logger.info(f"  → {len(df) - n_assigned} samples marked as 'Unknown' (missing coordinates)")
                    geo_analysis_performed = True
                else:
                    logger.warning("  ⚠ No samples with geographic-quality coordinates")
                    df_with_basins = df.copy()
                    df_with_basins[geo_category] = 'Unknown'

            except Exception as e:
                logger.error(f"  ✗ Failed to load custom shapefile: {e}")
                logger.warning("  Pipeline will continue without geographic analysis...")
                df_with_basins = df.copy()
                df_with_basins[geo_category] = 'Unknown'

        elif not cfg.geographic.goas_shapefile_path.exists():
            logger.warning(f"  ⚠ GOaS shapefile not found at: {cfg.geographic.goas_shapefile_path}")
            logger.warning("")
            logger.warning("  To enable geographic analysis, download the GOaS shapefile manually:")
            logger.warning("  1. Visit: https://www.marineregions.org/downloads.php")
            logger.warning("  2. Download: GOaS_v1_20211214.zip (registration form required)")
            logger.warning(f"  3. Extract to: {cfg.geographic.goas_shapefile_path.parent}")
            logger.warning("  4. Verify with: python -m boldgenotyper.goas_downloader")
            logger.warning("")
            logger.warning("  Pipeline will continue without geographic analysis...")
            df_with_basins = df.copy()
            df_with_basins[geo_category] = 'Unknown'
        else:
            # Use default GOaS shapefile
            logger.info("  → Using default GOaS ocean basin shapefile")
            try:
                goas_data = geographic.load_goas_data(cfg.geographic.goas_shapefile_path)

                # Only assign basins to samples with geographic quality coordinates
                df_geo_quality = df[df['is_geographic_quality']].copy()
                logger.info(f"  → {len(df_geo_quality)}/{len(df)} samples have geographic-quality coordinates")

                if len(df_geo_quality) > 0:
                    df_geo_assigned = geographic.assign_ocean_basins(
                        df_geo_quality, goas_data=goas_data, coord_col="coord"
                    )

                    # Merge back with full dataset, preserving lat/lon columns for visualization
                    df_with_basins = df.copy()
                    # First set all to Unknown
                    df_with_basins[geo_category] = 'Unknown'

                    # Initialize lat/lon columns if not present
                    if 'lat' not in df_with_basins.columns:
                        df_with_basins['lat'] = pd.NA
                    if 'lon' not in df_with_basins.columns:
                        df_with_basins['lon'] = pd.NA

                    # Extract columns to merge: ocean_basin and lat/lon (if they exist in df_geo_assigned)
                    cols_to_merge = ['processid', 'ocean_basin']
                    if 'lat' in df_geo_assigned.columns:
                        cols_to_merge.append('lat')
                    if 'lon' in df_geo_assigned.columns:
                        cols_to_merge.append('lon')

                    # Merge the geographic data back
                    merged_geo = df_with_basins[['processid']].merge(
                        df_geo_assigned[cols_to_merge], on='processid', how='left', suffixes=('', '_new')
                    )

                    # Update the columns (rename ocean_basin to geo_category if different)
                    if 'ocean_basin' in merged_geo.columns:
                        df_with_basins[geo_category] = merged_geo['ocean_basin'].fillna('Unknown')
                    if 'lat_new' in merged_geo.columns:
                        df_with_basins['lat'] = merged_geo['lat_new']
                    elif 'lat' in merged_geo.columns:
                        df_with_basins['lat'] = merged_geo['lat']
                    if 'lon_new' in merged_geo.columns:
                        df_with_basins['lon'] = merged_geo['lon_new']
                    elif 'lon' in merged_geo.columns:
                        df_with_basins['lon'] = merged_geo['lon']

                    n_assigned = (df_with_basins[geo_category] != 'Unknown').sum()
                    logger.info(f"  ✓ Assigned {geo_category} to {n_assigned} samples")
                    logger.info(f"  → {len(df) - n_assigned} samples marked as 'Unknown' (centroid/missing coordinates)")
                    geo_analysis_performed = True
                else:
                    logger.warning("  ⚠ No samples with geographic-quality coordinates")
                    df_with_basins = df.copy()
                    df_with_basins[geo_category] = 'Unknown'

            except Exception as e:
                logger.error(f"  ✗ Failed to load GOaS data: {e}")
                logger.warning("  Pipeline will continue without geographic analysis...")
                df_with_basins = df.copy()
                df_with_basins[geo_category] = 'Unknown'

        basins_tsv = dirs['intermediate_geographic'] / f"samples_with_{geo_category}.tsv"
        df_with_basins.to_csv(basins_tsv, sep='\t', index=False)

    except Exception as e:
        logger.error(f"Phase 1 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 2: Pre-processing Quality Control
    # ========================================================================
    logger.info("")
    logger.info("PHASE 2: Pre-processing Quality Control")
    logger.info("-" * 80)

    try:
        # Step 2.1: Generate initial FASTA from raw sequences
        logger.info("2.1: Generating initial FASTA from raw sequences...")
        logger.info(f"  Processing {len(df_with_basins)} samples")

        # Create dictionary: processid -> sequence
        sequences_dict = {}
        skipped_count = 0
        for _, row in df_with_basins.iterrows():
            processid = row['processid']
            sequence = row['nuc']

            if isinstance(sequence, str) and sequence.strip():
                sequences_dict[processid] = sequence.strip().upper()
            else:
                skipped_count += 1

        logger.info(f"  Loaded {len(sequences_dict)} sequences")
        if skipped_count > 0:
            logger.warning(f"  Skipped {skipped_count} samples with missing/empty sequences")

        # Step 2.2: Orientation normalization and ORF validation
        logger.info("2.2: Normalizing sequence orientation and validating ORF...")
        logger.info(f"  Using genetic code: {cfg.coi_validation.mitochondrial_code} (vertebrate mitochondrial)")
        logger.info(f"  Minimum ORF coverage: {cfg.coi_validation.orf_min_coverage:.0%}")
        logger.info(f"  Maximum internal stops: {cfg.coi_validation.orf_max_internal_stops}")

        corrected_sequences, orf_stats_df = quality_control.apply_orientation_normalization(
            sequences_dict=sequences_dict,
            genetic_code=cfg.coi_validation.mitochondrial_code,
            min_orf_coverage=cfg.coi_validation.orf_min_coverage,
            max_internal_stops=cfg.coi_validation.orf_max_internal_stops
        )

        # Save ORF validation results
        orf_stats_file = dirs['qc'] / f"{organism}_orf_validation.csv"
        orf_stats_df.to_csv(orf_stats_file, index=False)
        logger.info(f"  Saved ORF validation results: {orf_stats_file}")

        # Save oriented sequences
        oriented_fasta_records = [(pid, seq) for pid, seq in corrected_sequences.items()]
        oriented_fasta = dirs['qc'] / f"{organism}_oriented.fasta"
        utils.write_fasta(oriented_fasta_records, oriented_fasta)
        logger.info(f"  Saved oriented sequences: {oriented_fasta}")

        # Step 2.3: Dynamic QC filtering
        logger.info("2.3: Applying dynamic quality control filters...")
        logger.info(f"  Absolute minimum length: {cfg.qc.min_raw_length_abs} bp")
        logger.info(f"  Median-based minimum: {cfg.qc.min_raw_length_frac_of_median:.0%} of median")
        logger.info(f"  Maximum N fraction: {cfg.qc.max_raw_N_fraction:.1%}")
        logger.info(f"  Require valid ORF: {cfg.qc.require_valid_orf}")

        df_qc_passed, qc_stats = quality_control.apply_dynamic_qc_filters(
            df=df_with_basins,
            sequences_dict=corrected_sequences,
            orf_stats_df=orf_stats_df,
            min_raw_length_abs=cfg.qc.min_raw_length_abs,
            min_raw_length_frac_of_median=cfg.qc.min_raw_length_frac_of_median,
            max_raw_N_fraction=cfg.qc.max_raw_N_fraction,
            require_valid_orf=cfg.qc.require_valid_orf,
            processid_col='processid'
        )

        # Save QC-passed data
        qc_passed_tsv = dirs['qc'] / f"{organism}_qc_passed.tsv"
        df_qc_passed.to_csv(qc_passed_tsv, sep='\t', index=False)
        logger.info(f"  Saved QC-passed metadata: {qc_passed_tsv}")

        # Generate FASTA from QC-passed samples with corrected sequences
        logger.info("2.4: Generating FASTA from QC-passed samples...")
        qc_passed_fasta_records = []
        for _, row in df_qc_passed.iterrows():
            processid = row['processid']
            if processid in corrected_sequences:
                header = f"{organism}_{processid}.COI-5P"
                sequence = corrected_sequences[processid]
                qc_passed_fasta_records.append((header, sequence))

        qc_passed_fasta = dirs['qc'] / f"{organism}_qc_passed.fasta"
        utils.write_fasta(qc_passed_fasta_records, qc_passed_fasta)
        logger.info(f"  Created {len(qc_passed_fasta_records)} QC-passed FASTA records")
        logger.info(f"  Saved: {qc_passed_fasta}")

        # Log QC summary
        logger.info("QC Summary:")
        logger.info(f"  Total input sequences: {qc_stats['n_total']}")
        logger.info(f"  Passed all QC filters: {qc_stats['n_pass']} ({qc_stats['pass_rate']*100:.1f}%)")
        logger.info(f"  Failed QC filters: {qc_stats['n_fail']} ({(1-qc_stats['pass_rate'])*100:.1f}%)")

    except Exception as e:
        logger.error(f"Phase 2 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 3: Haplotype Discovery
    # ========================================================================
    logger.info("")
    logger.info("PHASE 3: Haplotype Discovery")
    logger.info("-" * 80)

    try:
        logger.info("3.1: Identifying haplotypes using ESV approach...")
        logger.info(f"  Core region minimum coverage: {cfg.core_region.core_min_coverage:.0%}")
        logger.info(f"  Gap masking threshold: {cfg.core_region.mask_gap_threshold:.0%}")
        logger.info(f"  Minimum core length: {cfg.core_region.core_min_length} bp")
        logger.info(f"  Singleton error filter threshold: >{cfg.haplotype.min_singleton_distance*100:.1f}% divergence")
        logger.info(f"  Singleton suspect flag threshold: >{cfg.haplotype.max_singleton_distance*100:.1f}% divergence")

        # Run haplotype discovery
        haplotype_records, haplotype_mapping, haplotype_stats = dereplication.identify_haplotypes(
            tsv_path=tsv_path,
            fasta_path=qc_passed_fasta,
            output_dir=dirs['haplotype_discovery'],
            min_core_coverage=cfg.core_region.core_min_coverage,
            gap_mask_threshold=cfg.core_region.mask_gap_threshold,
            min_core_length=cfg.core_region.core_min_length,
            min_singleton_distance=cfg.haplotype.min_singleton_distance,
            max_singleton_distance=cfg.haplotype.max_singleton_distance,
            orf_validation_df=orf_stats_df,
            cleanup_intermediates=False
        )

        logger.info(f"  Identified {len(haplotype_records)} high-confidence haplotypes (after error filtering)")

        # Copy haplotype outputs to main haplotypes directory
        # Move/copy haplotype FASTA
        source_haplotype_fasta = dirs['haplotype_discovery'] / f"{Path(tsv_path).stem}_haplotypes.fasta"
        haplotype_fasta = dirs['haplotypes'] / f"{organism}_haplotypes.fasta"
        if source_haplotype_fasta.exists():
            shutil.copy(source_haplotype_fasta, haplotype_fasta)
            logger.info(f"  Saved haplotypes: {haplotype_fasta}")

        # Move/copy mapping
        source_mapping = dirs['haplotype_discovery'] / f"{Path(tsv_path).stem}_haplotype_mapping.csv"
        haplotype_mapping_file = dirs['haplotypes'] / f"{organism}_haplotype_mapping.csv"
        if source_mapping.exists():
            shutil.copy(source_mapping, haplotype_mapping_file)
            logger.info(f"  Saved haplotype mapping: {haplotype_mapping_file}")

        # Move/copy stats
        source_stats = dirs['haplotype_discovery'] / f"{Path(tsv_path).stem}_haplotype_stats.csv"
        haplotype_stats_file = dirs['haplotypes'] / f"{organism}_haplotype_stats.csv"
        if source_stats.exists():
            shutil.copy(source_stats, haplotype_stats_file)
            logger.info(f"  Saved haplotype statistics: {haplotype_stats_file}")

        # Log haplotype summary
        n_suspect = haplotype_stats['is_suspect'].sum() if 'is_suspect' in haplotype_stats.columns else 0
        n_singletons = haplotype_stats['is_singleton'].sum() if 'is_singleton' in haplotype_stats.columns else 0

        logger.info("Haplotype Summary:")
        logger.info(f"  Total haplotypes: {len(haplotype_records)}")
        logger.info(f"  Singletons: {n_singletons}")
        logger.info(f"  Suspect haplotypes: {n_suspect}")

    except Exception as e:
        logger.error(f"Phase 3 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 4: Haplotype Assignment
    # ========================================================================
    logger.info("")
    logger.info("PHASE 4: Haplotype Assignment")
    logger.info("-" * 80)

    try:
        logger.info("4.1: Assigning samples to haplotypes...")
        logger.info(f"  Using {len(df_qc_passed)} QC-passed samples")
        logger.info("  Method: Direct mapping from ESV haplotype discovery (no tie detection)")

        annotated_tsv = dirs['intermediate_assignments'] / f"{organism}_with_haplotypes.tsv"
        diagnostics_csv = dirs['assignments'] / f"{organism}_diagnostics.csv"

        # For ESV approach, use the direct mapping from Phase 3 (haplotype_mapping)
        # This avoids false ties from identity-based re-assignment

        # The haplotype_mapping has full FASTA headers (organism_processid.marker)
        # Need to extract base processid to match metadata
        # Format: "Rhizoprionodon_ANGBF12411-15.COI-5P" → "ANGBF12411-15"
        haplotype_mapping_clean = haplotype_mapping.copy()
        haplotype_mapping_clean['processid_base'] = haplotype_mapping_clean['processid'].str.replace(
            f'^{organism}_', '', regex=True
        ).str.replace(r'\.COI-5P$', '', regex=True)

        # Merge haplotype assignments with QC-passed metadata
        df_with_haplotypes = df_qc_passed.merge(
            haplotype_mapping_clean[['processid_base', 'haplotype_id', 'haplotype_number', 'n_members', 'is_singleton']],
            left_on='processid',
            right_on='processid_base',
            how='left'
        ).drop(columns=['processid_base'])

        # Add assignment statistics columns for consistency with old approach
        df_with_haplotypes['assigned_haplotype'] = df_with_haplotypes['haplotype_id']
        df_with_haplotypes['is_tie'] = False  # No ties in ESV direct mapping
        df_with_haplotypes['is_low_confidence'] = False  # All assignments are exact matches

        # Count assignments
        n_assigned = df_with_haplotypes['assigned_haplotype'].notna().sum()
        n_total = len(df_with_haplotypes)

        # Save annotated data
        df_with_haplotypes.to_csv(annotated_tsv, sep="\t", index=False)
        logger.info(f"  ✓ Wrote annotated metadata: {annotated_tsv}")

        # Create diagnostics file (minimal for ESV direct mapping)
        # Status values: 'assigned' (successful), 'no_sequence' (no valid sequence), 'filtered' (failed QC)
        diagnostics_data = []
        for _, row in df_with_haplotypes.iterrows():
            has_assignment = pd.notna(row['assigned_haplotype'])
            # Determine status
            if has_assignment:
                status = 'assigned'
            else:
                # Sample didn't get a haplotype - likely had no sequence or failed core region extraction
                status = 'no_sequence'

            diagnostics_data.append({
                'processid': row['processid'],
                'assigned_haplotype': row['assigned_haplotype'],
                'status': status,
                'is_tie': False,
                'is_low_confidence': False,
                'assignment_method': 'ESV_direct_mapping',
                'note': 'Exact match from core region during haplotype discovery' if has_assignment else 'No valid core region extracted'
            })

        diagnostics_df = pd.DataFrame(diagnostics_data)
        diagnostics_df.to_csv(diagnostics_csv, index=False)
        logger.info(f"  ✓ Wrote diagnostics: {diagnostics_csv}")

        logger.info(f"  ✓ Assigned {n_assigned}/{n_total} samples to haplotypes")
        logger.info(f"  ✓ Assignment rate: {n_assigned/n_total*100:.1f}%")
        logger.info(f"  ✓ Ties: 0 (ESV direct mapping - no ambiguous assignments)")

    except Exception as e:
        logger.error(f"Phase 4 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 5: Taxonomy Assignment
    # ========================================================================
    logger.info("")
    logger.info("PHASE 5: Taxonomy Assignment to Haplotypes")
    logger.info("-" * 80)

    try:
        logger.info("5.1: Assigning taxonomy to haplotypes...")
        all_haplotype_ids = []
        if haplotype_fasta.exists():
            all_haplotype_ids = [rec.id for rec in SeqIO.parse(haplotype_fasta, "fasta")]

        assign_table, species_counts = utils.assign_consensus_taxonomy(
            df_with_haplotypes,
            group_col="haplotype_id",
            species_col="species",
            genus_col="genus",
            majority_threshold=cfg.taxonomy.majority_species_threshold,
            all_groups=all_haplotype_ids or None
        )

        # Build haplotype_sp labels (e.g., "Species name h1_n8")
        def _join_label(row):
            sp = row["assigned_sp"]
            hap_id = row["haplotype_id"]
            if not sp:
                return hap_id if isinstance(hap_id, str) else ""
            return f"{sp} {hap_id}"

        assign_table["haplotype_sp"] = assign_table.apply(_join_label, axis=1)

        # Save taxonomy files
        species_counts_out = dirs['taxonomy'] / f"{organism}_species_by_haplotype.csv"
        assign_table_out = dirs['taxonomy'] / f"{organism}_haplotype_taxonomy.csv"
        species_counts.to_csv(species_counts_out, index=False)
        assign_table.to_csv(assign_table_out, index=False)

        # Merge taxonomy back into main dataframe
        df_with_haplotypes = df_with_haplotypes.merge(
            assign_table[["haplotype_id", "assigned_sp", "haplotype_sp", "assignment_level", "assignment_notes", "majority_fraction"]],
            on="haplotype_id",
            how="left",
            validate="many_to_one"
        )

        # Handle cluster sequence taxonomy if available
        cluster_seq_path = dirs['taxonomy'] / f"{organism}_haplotype_taxonomy_seq.csv"
        if cluster_seq_path.exists():
            cluster_seq_df = pd.read_csv(cluster_seq_path)
            df_with_haplotypes = df_with_haplotypes.merge(
                cluster_seq_df, on="haplotype_id", how="left", validate="many_to_one"
            )

            # Apply final taxonomy decision
            def _final_label(row):
                final_sp, final_level, prov = utils.pick_final_group_taxon(
                    cluster_sp=row.get("cluster_seq_sp", ""),
                    cluster_level=row.get("cluster_seq_level", ""),
                    cluster_id=row.get("cluster_seq_best_identity", 0.0),
                    cluster_qcov=row.get("cluster_seq_qcov", 0.0),
                    majority_sp=row.get("assigned_sp", ""),
                    majority_level=row.get("assignment_level", ""),
                    majority_frac=row.get("majority_fraction", 0.0),
                    cfg_taxonomy=cfg.taxonomy,
                )
                hap_id = row["haplotype_id"] if isinstance(row.get("haplotype_id"), str) else ""
                label = f"{final_sp} {hap_id}".strip() if final_sp else hap_id
                return pd.Series({
                    "final_group_sp": final_sp,
                    "final_group_level": final_level,
                    "tax_provenance": prov,
                    "haplotype_sp": label
                })

            final_cols = df_with_haplotypes.apply(_final_label, axis=1)
            for c in final_cols.columns:
                df_with_haplotypes[c] = final_cols[c]

        logger.info(f"  ✓ Assigned taxonomy to {len(assign_table)} haplotypes")

        # Merge with geographic data (only if not already present in df_with_haplotypes)
        geo_cols_to_merge = []
        for col in ['lat', 'lon', geo_category]:
            if col in df_qc_passed.columns and col not in df_with_haplotypes.columns:
                geo_cols_to_merge.append(col)

        if geo_cols_to_merge:
            # Only merge if there are columns to add
            geo_keep = ['processid'] + geo_cols_to_merge
            df_final = df_with_haplotypes.merge(df_qc_passed[geo_keep], on='processid', how='left', validate='one_to_one')
        else:
            # df_with_haplotypes already has all geographic data
            df_final = df_with_haplotypes.copy()

        # Save final annotated file with proper CSV quoting to handle fields with commas
        annotated_csv = dirs['base'] / f"{organism}_annotated.csv"
        df_final.to_csv(annotated_csv, index=False, quoting=1)  # quoting=1 is csv.QUOTE_MINIMAL

        # Verify file was created
        if annotated_csv.exists():
            logger.info(f"  ✓ Saved annotated dataset: {annotated_csv}")
        else:
            logger.warning(f"  ⚠ Failed to save annotated dataset: {annotated_csv}")

    except Exception as e:
        logger.error(f"Phase 5 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 6: Post-assignment Quality Control
    # ========================================================================
    logger.info("")
    logger.info("PHASE 6: Post-assignment Quality Control")
    logger.info("-" * 80)

    try:
        logger.info("6.1: Adding contamination analysis columns to main output...")

        # Add contamination columns to the dataframe
        df_final_enhanced = quality_control.add_contamination_columns(
            df_final,
            group_col="haplotype_id",
            species_col="species",
            notes_col="notes" if "notes" in df_final.columns else None,
            majority_threshold=cfg.taxonomy.majority_species_threshold
        )

        # Save enhanced annotated file
        annotated_csv_enhanced = dirs['base'] / f"{organism}_annotated.csv"
        df_final_enhanced.to_csv(annotated_csv_enhanced, index=False, quoting=1)
        logger.info(f"  ✓ Enhanced annotated dataset with QC columns: {annotated_csv_enhanced}")

        # Update df_final to use the enhanced version
        df_final = df_final_enhanced

        logger.info("6.2: Generating quality control reports...")

        # Generate comprehensive QC report
        qc_results = quality_control.generate_quality_control_report(
            df_final,
            output_dir=dirs['quality_control'],
            organism=organism,
            group_col="haplotype_id",
            species_col="species"
        )

        logger.info(f"  ✓ Quality control reports saved to: {dirs['quality_control']}")

        # Print quality alerts to console
        quality_control.print_quality_alert(qc_results)

    except Exception as e:
        logger.warning(f"Quality control analysis failed (non-critical): {e}")
        logger.debug("QC error details:", exc_info=True)

    # ========================================================================
    # PHASE 6.5: Species-Level Aggregation
    # ========================================================================
    logger.info("")
    logger.info("PHASE 6.5: Species-Level Aggregation")
    logger.info("-" * 80)

    try:
        logger.info("6.5.1: Aggregating samples by assigned species...")

        # Check if we have taxonomy assignments
        if 'assigned_sp' not in df_final.columns or 'majority_fraction' not in df_final.columns:
            logger.warning("  ⊘ Skipping species-level analysis (taxonomy assignments not found)")
        else:
            # Create species_composition DataFrame from taxonomy data
            # Map from assign_table columns to expected species_composition columns
            species_composition = assign_table.copy()
            species_composition = species_composition.rename(columns={
                'assigned_sp': 'primary_species',
                'majority_fraction': 'primary_species_pct'
            })

            # Determine if haplotype is ambiguous
            # Ambiguous if: (1) assignment_level is not 'species', or (2) majority_fraction < 0.7
            # fillna(True) ensures NaN values are treated as ambiguous (safe default)
            species_composition['is_ambiguous'] = (
                (species_composition['assignment_level'] != 'species') |
                (species_composition['primary_species_pct'] < 0.7)
            ).fillna(True).astype(bool)

            # Determine if haplotype has multiple species (check species_counts for multi-species haplotypes)
            species_composition['is_multi_species'] = False  # Will be updated below

            # Count how many species each haplotype has (from species_counts)
            haplotype_species_count = species_counts.groupby('haplotype_id')['reported_species'].nunique()
            multi_species_haplotypes = haplotype_species_count[haplotype_species_count > 1].index
            species_composition.loc[
                species_composition['haplotype_id'].isin(multi_species_haplotypes),
                'is_multi_species'
            ] = True

            # Run species aggregation
            species_assignments, species_summary = species_analysis.aggregate_samples_by_species(
                annotated_metadata=df_final,
                species_composition=species_composition,
                min_confidence=0.7,
                output_dir=dirs['species_analysis']
            )

            logger.info("6.5.2: Calculating species-level diversity metrics...")
            species_diversity = species_analysis.calculate_species_diversity(
                species_assignments=species_assignments,
                output_dir=dirs['species_analysis']
            )

            logger.info("6.5.3: Generating species-level geographic summary...")
            species_geographic = species_analysis.generate_species_geographic_summary(
                species_assignments=species_assignments,
                output_dir=dirs['species_analysis']
            )

            logger.info("6.5.4: Generating species-faceted haplotype subsets...")
            species_subsets = species_analysis.generate_species_faceted_haplotype_subsets(
                species_assignments=species_assignments,
                output_dir=dirs['species_analysis'],
                min_haplotypes=2
            )
            logger.info(f"  → Created {len(species_subsets)} species-specific haplotype subsets")

            logger.info(f"  ✓ Species-level analysis complete: {dirs['species_analysis']}")

    except Exception as e:
        logger.warning(f"Species-level aggregation failed (non-critical): {e}")
        logger.debug("Species aggregation error details:", exc_info=True)

    # ========================================================================
    # PHASE 7: Geographic Analysis Enhancement
    # ========================================================================
    if geo_analysis_performed:
        logger.info("")
        logger.info("PHASE 7: Geographic Analysis Enhancement")
        logger.info("-" * 80)

        try:
            logger.info("7.1: Assessing geographic data coverage and quality...")

            # Run geographic enhancement
            geo_enhancement_results = geographic_enhancement.enhance_geographic_analysis(
                df=df_final,
                output_dir=dirs['geographic_analysis'],
                organism=organism,
                group_col="haplotype_id",
                goas_data=goas_data if 'goas_data' in locals() else None,
                geo_category=geo_category
            )

            # Update df_final with enhanced geographic data
            if 'enhanced_df' in geo_enhancement_results:
                df_final = geo_enhancement_results['enhanced_df']

                # Re-save the enhanced annotated CSV
                annotated_csv_enhanced = dirs['base'] / f"{organism}_annotated.csv"
                df_final.to_csv(annotated_csv_enhanced, index=False, quoting=1)
                logger.info(f"  ✓ Updated annotated dataset with geographic enhancements")

            logger.info(f"  ✓ Geographic enhancement complete: {dirs['geographic_analysis']}")

        except Exception as e:
            logger.warning(f"Geographic enhancement failed (non-critical): {e}")
            logger.debug("Geographic enhancement error details:", exc_info=True)

    # ========================================================================
    # PHASE 8: Phylogenetic Analysis (Optional)
    # ========================================================================
    logger.info("")
    logger.info("PHASE 8: Phylogenetic Analysis")
    logger.info("-" * 80)

    tree_path = None
    if cfg.phylogenetic.build_tree:
        try:
            logger.info("8.1: Building phylogenetic tree from haplotypes...")

            # Check for required tools
            if not utils.check_external_tool("mafft"):
                logger.warning("  ⚠ MAFFT not found, skipping tree building")
            elif not utils.check_external_tool("fasttree"):
                logger.warning("  ⚠ FastTree not found, skipping tree building")
            else:
                # Build phylogenetic tree using haplotypes
                # Save alignment files to intermediate, tree files to final phylogenetic directory
                intermediate_prefix = dirs['intermediate_phylo'] / organism
                output_prefix = dirs['phylogenetic'] / organism

                tree = phylogenetics.build_phylogeny(
                    consensus_fasta=str(haplotype_fasta),
                    output_prefix=str(intermediate_prefix),
                    threads=cfg.n_threads,
                    min_consensus_length=cfg.phylogenetic.min_consensus_length,
                    min_cluster_size=cfg.phylogenetic.min_cluster_size,
                    outgroup_fasta=str(cfg.phylogenetic.outgroup_fasta) if cfg.phylogenetic.outgroup_fasta else None,
                    trim_alignment=cfg.phylogenetic.trim_alignment,
                    trim_method=cfg.phylogenetic.trim_method
                )

                # Move tree files from intermediate to final directory
                intermediate_tree = f"{intermediate_prefix}_tree.nwk"
                tree_path = f"{output_prefix}_tree.nwk"
                if Path(intermediate_tree).exists():
                    Path(intermediate_tree).rename(tree_path)

                # Verify tree file was actually created
                if tree is not None and Path(tree_path).exists():
                    # Optional rerooting using in-tree outgroup label or taxon
                    if cfg.phylogenetic.outgroup_label:
                        phylogenetics.reroot_tree_by_label(
                            tree_file=tree_path,
                            outgroup_label=cfg.phylogenetic.outgroup_label,
                            output_file=tree_path
                        )
                    elif cfg.phylogenetic.outgroup_taxon:
                        phylogenetics.reroot_tree_by_taxon(
                            tree_file=tree_path,
                            taxonomy_csv=str(dirs['taxonomy'] / f"{organism}_haplotype_taxonomy.csv"),
                            taxon_query=cfg.phylogenetic.outgroup_taxon,
                            output_file=tree_path
                        )

                    logger.info(f"  ✓ Built phylogenetic tree: {tree_path}")

                    # Create relabeled versions with haplotype_sp labels
                    logger.info("8.2: Creating relabeled tree and alignment files...")
                    try:
                        alignment_path = f"{intermediate_prefix}_aligned.fasta"
                        taxonomy_csv_path = dirs['taxonomy'] / f"{organism}_haplotype_taxonomy.csv"

                        if Path(alignment_path).exists() and taxonomy_csv_path.exists():
                            relabeled_tree_path = f"{output_prefix}_tree_relabeled.nwk"
                            relabeled_alignment_path = f"{intermediate_prefix}_aligned_relabeled.fasta"

                            phylogenetics.relabel_tree_and_alignment(
                                tree_file=tree_path,
                                alignment_file=alignment_path,
                                taxonomy_csv=str(taxonomy_csv_path),
                                output_tree=relabeled_tree_path,
                                output_alignment=relabeled_alignment_path,
                                label_column="haplotype_sp",
                                id_column="haplotype_id",
                            )
                            logger.info(f"  ✓ Created relabeled tree: {relabeled_tree_path}")
                            logger.info(f"  ✓ Created relabeled alignment (intermediate): {relabeled_alignment_path}")
                        else:
                            logger.warning("  ⚠ Skipping relabeling: alignment or taxonomy file not found")
                    except Exception as e:
                        logger.warning(f"  ⚠ Relabeling failed (non-critical): {e}")

                    # Generate MSA visualization if enabled
                    if cfg.msa.enabled:
                        logger.info("8.3: Generating MSA visualization...")
                        try:
                            # Use the alignment file that exists (prefer trimmed if available)
                            alignment_file = None
                            if Path(f"{intermediate_prefix}_aligned_trimmed.fasta").exists():
                                alignment_file = Path(f"{intermediate_prefix}_aligned_trimmed.fasta")
                            elif Path(f"{intermediate_prefix}_aligned.fasta").exists():
                                alignment_file = Path(f"{intermediate_prefix}_aligned.fasta")

                            if alignment_file and alignment_file.exists():
                                msa_plots = msa_visualization.create_phylo_ordered_msa(
                                    alignment_path=alignment_file,
                                    tree_path=Path(tree_path),
                                    output_dir=dirs['phylogenetic'],
                                    organism=organism,
                                    chunk_size=cfg.msa.chunk_size,
                                    max_sequences=cfg.msa.max_sequences,
                                    color_scheme=cfg.msa.color_scheme,
                                    show_consensus=cfg.msa.show_consensus,
                                    show_logo=cfg.msa.show_logo,
                                    output_formats=cfg.msa.output_formats
                                )
                                if msa_plots:
                                    logger.info(f"  ✓ Generated {len(msa_plots)} MSA plot(s)")
                                else:
                                    logger.warning("  ⚠ No MSA plots generated")
                            else:
                                logger.warning("  ⚠ Alignment file not found for MSA visualization")
                        except ImportError:
                            logger.warning(
                                "  ⚠ MSA visualization skipped: pymsaviz not installed. "
                                "Install with: pip install pymsaviz>=0.4.0"
                            )
                        except Exception as e:
                            logger.warning(f"  ⚠ MSA visualization failed (non-critical): {e}")
                            logger.debug("Full traceback:", exc_info=True)
                    else:
                        logger.info("  ⊘ MSA visualization disabled in config")

                else:
                    logger.warning(f"  ⚠ Phylogenetic tree building completed but output file not found: {tree_path}")
                    tree_path = None

        except Exception as e:
            logger.warning(f"Phylogenetic analysis failed (non-critical): {e}")
    else:
        logger.info("  ⊘ Phylogenetic tree building disabled in config")

    # ========================================================================
    # PHASE 9: Divergence Analysis
    # ========================================================================
    logger.info("")
    logger.info("PHASE 9: Divergence Analysis")
    logger.info("-" * 80)

    try:
        logger.info("9.1: Calculating pairwise divergence between haplotypes...")

        # Check if we have the required files
        taxonomy_csv_path = dirs['taxonomy'] / f"{organism}_haplotype_taxonomy.csv"

        if haplotype_fasta.exists() and taxonomy_csv_path.exists():
            # Run divergence analysis
            divergence_results = divergence_analysis.generate_divergence_analysis(
                consensus_fasta=str(haplotype_fasta),
                taxonomy_csv=str(taxonomy_csv_path),
                output_dir=dirs['divergence_analysis'],
                organism=organism
            )

            logger.info(f"  ✓ Divergence analysis complete: {dirs['divergence_analysis']}")

            # Print summary to console
            divergence_analysis.print_divergence_summary(divergence_results)

            # Species-level divergence analysis (if species data available)
            logger.info("9.2: Calculating species-level divergence...")
            species_summary_csv = dirs['species_analysis'] / 'species_summary.csv'

            if species_summary_csv.exists() and 'divergence_matrix' in divergence_results:
                try:
                    # Load divergence matrix
                    div_matrix_path = divergence_results['divergence_matrix']
                    div_matrix = pd.read_csv(div_matrix_path, index_col=0)

                    # Run species-level divergence analysis
                    species_div_results = divergence_analysis.generate_species_divergence_analysis(
                        divergence_matrix=div_matrix,
                        species_summary_csv=species_summary_csv,
                        haplotype_taxonomy_csv=taxonomy_csv_path,
                        output_dir=dirs['species_analysis']
                    )

                    logger.info(f"  ✓ Species-level divergence analysis complete")

                    # Species-faceted divergence (within-species haplotype divergence)
                    logger.info("9.3: Calculating within-species divergence matrices...")
                    try:
                        species_faceted_div_results = species_analysis.generate_within_species_divergence_matrices(
                            divergence_matrix=div_matrix,
                            species_summary_csv=species_summary_csv,
                            haplotype_taxonomy_csv=taxonomy_csv_path,
                            output_dir=dirs['species_analysis'],
                            min_haplotypes=2
                        )
                        n_species_with_div = len(species_faceted_div_results)
                        logger.info(f"  ✓ Generated within-species divergence matrices for {n_species_with_div} species")
                    except Exception as e:
                        logger.warning(f"Within-species divergence analysis failed (non-critical): {e}")
                        logger.debug("Within-species divergence error details:", exc_info=True)

                except Exception as e:
                    logger.warning(f"Species-level divergence analysis failed (non-critical): {e}")
                    logger.debug("Species divergence error details:", exc_info=True)
            else:
                logger.info("  ⊘ Skipping species-level divergence (species data not available)")

        else:
            logger.warning("  ⊘ Skipping divergence analysis (haplotype or taxonomy files not found)")

    except Exception as e:
        logger.warning(f"Divergence analysis failed (non-critical): {e}")
        logger.debug("Divergence analysis error details:", exc_info=True)

    # ========================================================================
    # PHASE 9.5: Metadata Summarization
    # ========================================================================
    metadata_summary_results = None
    if metadata_summary_enabled:
        logger.info("")
        logger.info("PHASE 9.5: Metadata Summarization")
        logger.info("-" * 80)

        try:
            color_map = visualization.build_genotype_colors_dict(df_final, 'haplotype_sp')

            metadata_summary_results = metadata_summary.run_metadata_summary(
                annotated_df=df_final,
                output_dir=output_dir,
                organism=organism,
                fields=metadata_fields,
                normalize_sex=normalize_sex,
                temporal_analysis=temporal_analysis,
                color_map=color_map,
                haplotype_col='haplotype_sp'
            )

            logger.info(f"  ✓ Metadata summarization complete")
            logger.info(f"    - Fields tabulated: {metadata_summary_results.get('n_fields_tabulated', 0)}")
            logger.info(f"    - Visualizations created: {metadata_summary_results.get('n_visualizations', 0)}")
            if metadata_summary_results.get('temporal_analysis'):
                logger.info(f"    - Temporal analysis: enabled")
            if metadata_summary_results.get('missing_fields'):
                logger.info(f"    - Missing fields (skipped): {', '.join(metadata_summary_results.get('missing_fields', []))}")

        except Exception as e:
            logger.warning(f"Metadata summarization failed (non-critical): {e}")
            logger.debug("Metadata summarization error details:", exc_info=True)
    else:
        logger.info("")
        logger.info("PHASE 9.5: Metadata Summarization - SKIPPED (--no-metadata-summary was specified)")

    # ========================================================================
    # PHASE 10: Visualization
    # ========================================================================
    logger.info("")
    logger.info("PHASE 10: Visualization")
    logger.info("-" * 80)

    try:
        logger.info("10.1: Generating plots...")

        # Generate visualizations for each format
        for fmt in cfg.visualization.figure_format:
            # Geographic visualizations - only if geographic analysis was performed
            if geo_analysis_performed:
                # Distribution maps
                if 'lat' in df_final.columns and 'lon' in df_final.columns and 'haplotype_sp' in df_final.columns:
                    try:
                        # Determine background detail: use "full" for main plots
                        detail = "full" if cfg.visualization.map_background_detail == "full" else cfg.visualization.map_background_detail
                        if cfg.visualization.map_background_detail == "auto":
                            detail = "full"  # Main plot gets full detail

                        map_path, map_data = visualization.plot_distribution_map(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_distribution_map", fmt)),
                            genotype_column='haplotype_sp',
                            latitude_col='lat',
                            longitude_col='lon',
                            dpi=cfg.visualization.figure_dpi,
                            map_background_detail=detail
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if map_data:
                            json_path = dirs['visualization_json'] / f"{organism}_distribution_map_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(map_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.warning(f"Distribution map skipped: {e}")

                # Geographic region abundance bar plot (relative)
                if geo_category in df_final.columns and 'haplotype_sp' in df_final.columns:
                    try:
                        bar_path, bar_data = visualization.plot_ocean_basin_abundance(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_distribution_bar", fmt)),
                            genotype_column='haplotype_sp',
                            basin_column=geo_category
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if bar_data:
                            json_path = dirs['visualization_json'] / f"{organism}_distribution_bar_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(bar_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Geographic region bar plot skipped: {e}")

                # Geographic region abundance bar plot (total counts)
                if geo_category in df_final.columns and 'haplotype_sp' in df_final.columns:
                    try:
                        total_bar_path, total_bar_data = visualization.plot_ocean_basin_abundance_total(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_totaldistribution_bar", fmt)),
                            genotype_column='haplotype_sp',
                            basin_column=geo_category
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if total_bar_data:
                            json_path = dirs['visualization_json'] / f"{organism}_totaldistribution_bar_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(total_bar_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Geographic region total abundance bar plot skipped: {e}")

                # Species-level visualizations (if species data available)
                species_assignments_csv = dirs['species_analysis'] / 'species_assignments.csv'
                if species_assignments_csv.exists():
                    try:
                        # Load species assignments
                        species_df = pd.read_csv(species_assignments_csv)

                        # Species distribution map
                        if 'lat' in species_df.columns and 'lon' in species_df.columns and 'primary_species' in species_df.columns:
                            try:
                                # Use full detail for species-level main plots
                                detail = "full" if cfg.visualization.map_background_detail == "full" else cfg.visualization.map_background_detail
                                if cfg.visualization.map_background_detail == "auto":
                                    detail = "full"

                                species_map_path, species_map_data = visualization.plot_species_distribution_map(
                                    df=species_df,
                                    output_path=str(get_viz_path(dirs, f"{organism}_species_distribution_map", fmt)),
                                    species_column='primary_species',
                                    latitude_col='lat',
                                    longitude_col='lon',
                                    dpi=cfg.visualization.figure_dpi,
                                    map_background_detail=detail
                                )
                                if species_map_data:
                                    json_path = dirs['visualization_json'] / f"{organism}_species_distribution_map_data.json"
                                    with open(json_path, 'w') as f:
                                        json.dump(species_map_data, f, indent=2)
                                    logger.debug(f"Saved species map data to: {json_path}")
                            except Exception as e:
                                logger.debug(f"Species distribution map skipped: {e}")

                        # Species geographic region abundance (relative)
                        if geo_category in species_df.columns and 'primary_species' in species_df.columns:
                            try:
                                species_bar_path, species_bar_data = visualization.plot_species_ocean_basin_abundance(
                                    df=species_df,
                                    output_path=str(get_viz_path(dirs, f"{organism}_species_distribution_bar", fmt)),
                                    species_column='primary_species',
                                    basin_column=geo_category
                                )
                                if species_bar_data:
                                    json_path = dirs['visualization_json'] / f"{organism}_species_distribution_bar_data.json"
                                    with open(json_path, 'w') as f:
                                        json.dump(species_bar_data, f, indent=2)
                                    logger.debug(f"Saved species bar data to: {json_path}")
                            except Exception as e:
                                logger.debug(f"Species region bar plot skipped: {e}")

                        # Species geographic region abundance (total counts)
                        if geo_category in species_df.columns and 'primary_species' in species_df.columns:
                            try:
                                species_total_bar_path, species_total_bar_data = visualization.plot_species_ocean_basin_abundance_total(
                                    df=species_df,
                                    output_path=str(get_viz_path(dirs, f"{organism}_species_totaldistribution_bar", fmt)),
                                    species_column='primary_species',
                                    basin_column=geo_category
                                )
                                if species_total_bar_data:
                                    json_path = dirs['visualization_json'] / f"{organism}_species_totaldistribution_bar_data.json"
                                    with open(json_path, 'w') as f:
                                        json.dump(species_total_bar_data, f, indent=2)
                                    logger.debug(f"Saved species total bar data to: {json_path}")
                            except Exception as e:
                                logger.debug(f"Species geographic region total abundance bar plot skipped: {e}")

                        # Species-faceted haplotype visualizations
                        # (separate plots for each species showing haplotypes)
                        try:
                            logger.info("  → Generating species-faceted haplotype distribution maps...")
                            # Use simplified background for species-faceted plots (many maps)
                            detail = "simple" if cfg.visualization.map_background_detail == "auto" else cfg.visualization.map_background_detail

                            # For species-faceted plots, use format-specific subdirectory
                            species_facet_dir = dirs['visualization_pdf'] if fmt == 'pdf' else (
                                dirs['visualization_svg'] if fmt == 'svg' else dirs['visualization_pdf']
                            )
                            species_maps_results = visualization.plot_haplotypes_by_species_maps(
                                species_assignments=species_df,
                                output_dir=species_facet_dir,
                                species_column='primary_species',
                                haplotype_column='haplotype_sp',
                                latitude_col='lat',
                                longitude_col='lon',
                                min_haplotypes=2,
                                figure_format=fmt,
                                dpi=cfg.visualization.facet_dpi,
                                map_background_detail=detail
                            )
                            logger.info(f"  ✓ Created {len(species_maps_results)} species-faceted haplotype maps")
                        except Exception as e:
                            logger.debug(f"Species-faceted haplotype maps skipped: {e}")

                        try:
                            logger.info("  → Generating species-faceted haplotype basin charts...")
                            # For species-faceted plots, use format-specific subdirectory
                            species_bar_facet_dir = dirs['visualization_pdf'] if fmt == 'pdf' else (
                                dirs['visualization_svg'] if fmt == 'svg' else dirs['visualization_pdf']
                            )
                            species_bars_results = visualization.plot_haplotypes_by_species_basin_bars(
                                species_assignments=species_df,
                                output_dir=species_bar_facet_dir,
                                species_column='primary_species',
                                haplotype_column='haplotype_sp',
                                basin_column=geo_category,
                                min_haplotypes=2,
                                figure_format=fmt
                            )
                            logger.info(f"  ✓ Created {len(species_bars_results)} species-faceted haplotype bar charts")
                        except Exception as e:
                            logger.debug(f"Species-faceted haplotype bar charts skipped: {e}")

                    except Exception as e:
                        logger.debug(f"Species-level visualizations skipped: {e}")

                # Faceted distribution map by haplotype_sp
                if ('lat' in df_final.columns and 'lon' in df_final.columns and
                    'haplotype_sp' in df_final.columns and 'haplotype_id' in df_final.columns):
                    try:
                        # Use simplified background for faceted plots (better performance)
                        detail = "simple" if cfg.visualization.map_background_detail == "auto" else cfg.visualization.map_background_detail

                        # Only save individual facets if enabled (off by default for performance)
                        facet_formats = cfg.visualization.facet_formats if cfg.visualization.save_individual_facets else None

                        visualization.plot_distribution_map_faceted(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_distribution_map_faceted", fmt)),
                            genotype_column='haplotype_sp',
                            species_column='assigned_sp',
                            latitude_col='lat',
                            longitude_col='lon',
                            dpi=cfg.visualization.facet_dpi,
                            facet_by='genotype',
                            map_buffer_degrees=cfg.visualization.map_buffer_degrees,
                            show_unknown_annotation=cfg.visualization.show_unknown_geography_annotation,
                            show_scale_bar=cfg.visualization.show_scale_bar,
                            genotype_plots_dir=dirs['haplotype_plots'],
                            formats=facet_formats,
                            map_background_detail=detail
                        )

                        # Additional species-faceted map (explicit file name for clarity)
                        visualization.plot_distribution_map_faceted(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_distribution_map_species_faceted", fmt)),
                            genotype_column='haplotype_sp',
                            species_column='assigned_sp',
                            latitude_col='lat',
                            longitude_col='lon',
                            dpi=cfg.visualization.facet_dpi,
                            facet_by='species',
                            map_buffer_degrees=cfg.visualization.map_buffer_degrees,
                            show_unknown_annotation=cfg.visualization.show_unknown_geography_annotation,
                            show_scale_bar=cfg.visualization.show_scale_bar,
                            genotype_plots_dir=dirs['haplotype_plots'],
                            formats=facet_formats,
                            map_background_detail=detail
                        )
                    except Exception as e:
                        logger.warning(f"Faceted distribution map generation failed: {e}", exc_info=True)

                # Faceted geographic region bar plot by haplotype_sp
                if (geo_category in df_final.columns and 'haplotype_sp' in df_final.columns and
                    'haplotype_id' in df_final.columns):
                    try:
                        # Only save individual facets if enabled (off by default for performance)
                        facet_formats = cfg.visualization.facet_formats if cfg.visualization.save_individual_facets else None

                        faceted_bar_path, faceted_bar_data = visualization.plot_ocean_basin_abundance_faceted(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_distribution_bar_faceted", fmt)),
                            genotype_column='haplotype_id',
                            species_column='assigned_sp',
                            basin_column=geo_category,
                            dpi=cfg.visualization.facet_dpi,
                            facet_by='genotype',
                            genotype_plots_dir=dirs['haplotype_plots'],
                            formats=facet_formats
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if faceted_bar_data:
                            json_path = dirs['visualization_json'] / f"{organism}_distribution_bar_faceted_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(faceted_bar_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")

                        # Species-faceted version (explicit filename for clarity)
                        species_bar_path, species_bar_data = visualization.plot_ocean_basin_abundance_faceted(
                            df=df_final,
                            output_path=str(get_viz_path(dirs, f"{organism}_distribution_bar_species_faceted", fmt)),
                            genotype_column='haplotype_id',
                            species_column='assigned_sp',
                            basin_column=geo_category,
                            dpi=cfg.visualization.facet_dpi,
                            facet_by='species',
                            genotype_plots_dir=dirs['haplotype_plots'],
                            formats=facet_formats
                        )
                        if species_bar_data:
                            json_path = dirs['visualization_json'] / f"{organism}_distribution_bar_species_faceted_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(species_bar_data, f, indent=2)
                            logger.debug(f"Saved species-faceted plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Faceted geographic region bar plot skipped: {e}")
            else:
                logger.info("  ⊘ Skipping geographic visualizations (geographic analysis not performed)")

            # Identity distribution (always generated if diagnostics exist)
            if diagnostics_csv.exists():
                try:
                    # Use intermediate DPI for diagnostic plots
                    visualization.plot_identity_distribution(
                        diagnostics_csv=str(diagnostics_csv),
                        output_path=str(dirs['assignments'] / f"{organism}_identity_distribution.{fmt}"),
                        figsize=(10, 6),
                        dpi=cfg.visualization.intermediate_dpi
                    )
                except Exception as e:
                    logger.warning(f"Identity distribution plot failed: {e}")
                    logger.debug("Full traceback:", exc_info=True)

            # Phylogenetic tree
            if tree_path and Path(tree_path).exists():
                try:
                    # Use relabeled tree if it exists (so tips show haplotype_sp labels)
                    relabeled_tree_path = tree_path.replace("_tree.nwk", "_tree_relabeled.nwk")
                    if Path(relabeled_tree_path).exists():
                        tree_to_plot = relabeled_tree_path
                        logger.info(f"Using relabeled tree for visualization: {relabeled_tree_path}")
                    else:
                        tree_to_plot = tree_path

                    # Load color map if it exists
                    color_map_path = dirs['assignments'] / f"{organism}_haplotype_color_map.csv"
                    haplotype_colors = None
                    if color_map_path.exists():
                        color_df = pd.read_csv(color_map_path)
                        haplotype_colors = dict(zip(color_df['haplotype_sp'], color_df['color']))

                    visualization.plot_phylogenetic_tree(
                        tree_file=str(tree_to_plot),
                        output_path=str(dirs['phylogenetic'] / f"{organism}_tree.{fmt}"),
                        genotype_colors=haplotype_colors,
                        show_bootstrap=True,
                        bootstrap_threshold=cfg.visualization.show_bootstrap_threshold,
                        figsize=None,  # Auto-scale based on tree size
                        dpi=cfg.visualization.figure_dpi
                    )
                    logger.info(f"  ✓ Generated phylogenetic tree plot: {organism}_tree.{fmt}")
                except Exception as e:
                    logger.warning(f"Phylogenetic tree plot failed: {e}")
                    logger.debug("Full traceback:", exc_info=True)

        logger.info(f"  ✓ Generated visualization plots")

    except Exception as e:
        logger.warning(f"Visualization generation encountered errors (non-critical): {e}")

    # ========================================================================
    # PHASE 10.5: Plot Data Export (Optional)
    # ========================================================================
    if export_plot_data:
        logger.info("")
        logger.info("PHASE 10.5: Plot Data Export")
        logger.info("-" * 80)

        try:
            # Get color map if it exists
            color_map_path = dirs['assignments'] / f"{organism}_haplotype_color_map.csv"
            haplotype_colors = None
            if color_map_path.exists():
                color_df = pd.read_csv(color_map_path)
                haplotype_colors = dict(zip(color_df['haplotype_sp'], color_df['color']))

            # Export plot data and regeneration kit
            plot_export_results = plot_export.export_plots_complete(
                df=df_final,
                output_dir=output_dir,
                organism=organism,
                consensus_path=haplotype_fasta if haplotype_fasta.exists() else None,
                tree_path=Path(tree_path) if tree_path and Path(tree_path).exists() else None,
                diagnostics_path=diagnostics_csv if diagnostics_csv.exists() else None,
                color_map=haplotype_colors
            )

            logger.info(f"  ✓ Plot regeneration kit exported to {output_dir / 'plots'}")

        except Exception as e:
            logger.warning(f"Plot data export failed (non-critical): {e}")
            logger.debug("Plot export error details:", exc_info=True)

    # ========================================================================
    # PHASE 11: Population Genetics Export (Optional)
    # ========================================================================
    if export_popgen_formats:
        logger.info("")
        logger.info("PHASE 11: Population Genetics Export")
        logger.info("-" * 80)

        try:
            # Export to population genetics software formats
            popgen_results = popgen_export.export_population_genetics_formats(
                df=df_final,
                consensus_fasta_path=str(haplotype_fasta),
                output_dir=output_dir,
                organism=organism,
                formats=export_popgen_formats,
                pop_column=geo_category,
                haplotype_column='haplotype_id',
            )

            logger.info(f"  ✓ Population genetics formats exported to {output_dir / 'exports'}")

        except Exception as e:
            logger.warning(f"Population genetics export failed (non-critical): {e}")
            logger.debug("Popgen export error details:", exc_info=True)

    # ========================================================================
    # PHASE 12: Reports
    # ========================================================================
    logger.info("")
    logger.info("PHASE 12: Generating Reports")
    logger.info("-" * 80)

    try:
        # Generate assignment summary report
        summary_output = dirs['reports'] / f"{organism}_assignment_summary.csv"
        reports.generate_assignment_summary(
            annotated_csv=str(annotated_csv),
            diagnostics_csv=str(diagnostics_csv),
            output_csv=str(summary_output)
        )

        # Verify file was created
        if summary_output.exists():
            logger.info(f"  ✓ Generated assignment summary: {summary_output}")
        else:
            logger.warning(f"  ⚠ Failed to generate assignment summary: {summary_output}")

    except Exception as e:
        logger.warning(f"Report generation failed (non-critical): {e}")

    # ========================================================================
    # Cleanup: Remove Empty Directories
    # ========================================================================
    logger.info("")
    logger.info("Cleaning up empty directories...")
    try:
        remove_empty_directories(output_dir)
        logger.info("  ✓ Removed empty directories")
    except Exception as e:
        logger.debug(f"Directory cleanup encountered minor issues: {e}")

    # ========================================================================
    # HTML Summary Report
    # ========================================================================
    html_report_path = None
    if not no_report:
        logger.info("")
        logger.info("Generating HTML summary report...")
        try:
            html_report_path = reports.generate_html_report(
                organism=organism,
                output_dir=output_dir,
                version=__version__
            )
            if html_report_path:
                logger.info(f"  ✓ Generated HTML report: {html_report_path}")
            else:
                logger.warning("  ⚠ HTML report generation returned None")
        except ImportError as e:
            logger.warning(f"  ⚠ HTML report generation skipped: {e}")
            logger.warning("    Install jinja2 to enable HTML reports: pip install jinja2")
        except Exception as e:
            logger.warning(f"  ⚠ HTML report generation failed (non-critical): {e}")
            logger.debug("HTML report error details:", exc_info=True)
    else:
        logger.info("")
        logger.info("HTML report generation skipped (--no-report)")

    # ========================================================================
    # Cleanup: Remove Intermediate Directory (unless --keep-intermediates)
    # ========================================================================
    if not cfg.keep_intermediates:
        intermediate_dir = output_dir / 'intermediate'
        if intermediate_dir.exists():
            try:
                shutil.rmtree(intermediate_dir)
                logger.info("")
                logger.info("Removed intermediate directory")
            except Exception as e:
                logger.debug(f"Could not remove intermediate directory: {e}")

    # ========================================================================
    # Pipeline Complete
    # ========================================================================
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"✓ Pipeline completed successfully for {organism}")
    logger.info(f"  Output directory: {output_dir}")
    logger.info(f"  Key files:")
    logger.info(f"    - Annotated data: {annotated_csv}")
    logger.info(f"    - Haplotype sequences: {haplotype_fasta}")
    if tree_path:
        logger.info(f"    - Phylogenetic tree: {tree_path}")
    if html_report_path:
        logger.info(f"    - HTML summary report: {html_report_path}")
    logger.info("=" * 80)

    return True

def main_sweep(argv=None) -> int:
    """
    CLI entry point for the 'boldgenotyper-sweep' command.

    Sweep min_singleton_distance across a range of values to find the
    optimal singleton error-filtering threshold for the dataset.

    Example:
        boldgenotyper-sweep Sphyrna_lewini.tsv \
            --thresholds 0.01,0.015,0.02,0.03,0.05 \
            --output parameter_sweep/
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper-sweep",
        description="Sweep min_singleton_distance to optimise singleton error filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test default threshold range
  boldgenotyper-sweep data/Sphyrna_lewini.tsv

  # Test custom thresholds
  boldgenotyper-sweep data/samples.tsv --thresholds 0.005,0.01,0.015,0.02

  # Run in parallel with 8 threads
  boldgenotyper-sweep data/samples.tsv --threads 8

For more information: https://github.com/SymbioSeas/BOLDGenotyper
        """
    )

    parser.add_argument(
        'tsv',
        type=Path,
        help='Input BOLD TSV file'
    )

    parser.add_argument(
        '--thresholds',
        type=str,
        default='0.01,0.015,0.02,0.03,0.05',
        help='Comma-separated threshold values to test (default: 0.01,0.015,0.02,0.03,0.05)'
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=Path('parameter_sweep'),
        help='Output directory (default: parameter_sweep/)'
    )

    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of parallel threads (default: 4, sequential by default)'
    )

    parser.add_argument(
        '--keep-intermediates',
        action='store_true',
        help='Keep full output for each threshold (default: False)'
    )

    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging verbosity (default: INFO)'
    )

    args = parser.parse_args(argv)

    # Setup logging
    utils.setup_logging(log_level=args.log_level, log_file=None)

    # Parse thresholds
    try:
        thresholds = [float(t.strip()) for t in args.thresholds.split(',')]
    except ValueError:
        print(f"Error: Invalid threshold format: {args.thresholds}", file=sys.stderr)
        print("  Expected comma-separated floats, e.g.: 0.01,0.015,0.02", file=sys.stderr)
        return 1

    # Validate input
    if not args.tsv.exists():
        print(f"Error: Input TSV not found: {args.tsv}", file=sys.stderr)
        return 1

    # Run parameter sweep
    try:
        results = parameter_sweep.run_parameter_sweep(
            tsv_path=args.tsv,
            thresholds=thresholds,
            output_dir=args.output,
            organism=None,  # Will be inferred
            threads=args.threads,
            keep_intermediates=args.keep_intermediates
        )

        print("\n✓ Parameter sweep complete!")
        recommended = results.get('recommended_threshold', 'N/A')
        if isinstance(recommended, (int, float)):
            print(f"  Recommended threshold: {recommended:.3f}")
        else:
            print(f"  Recommended threshold: {recommended}")
        print(f"  Results saved to: {args.output}")
        print("\nKey files:")
        print(f"  - {args.output / 'recommendations.txt'}")
        print(f"  - {args.output / 'sweep_summary.csv'}")
        print(f"  - {args.output / 'elbow_plot.pdf'}")
        print()

        return 0

    except Exception as e:
        logger.error(f"Parameter sweep failed: {e}", exc_info=True)
        print(f"\nError: {e}", file=sys.stderr)
        return 1


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='BOLDGenotyper: haplotype-first COI genotyping with optional phylogeny, divergence, and mapping outputs.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (organism inferred from TSV filename)
  boldgenotyper data/Euprymna_scolopes.tsv

  # Specify organism name and output directory
  boldgenotyper data/Carcharhinus.tsv --organism Carcharhinus --output results/Carcharhinus_analysis

  # Tune assignment thresholds for diverse taxa
  boldgenotyper data/Carcharhinus.tsv --similarity-threshold 0.80 --tie-margin 0.005 --tie-threshold 0.97

  # Enable phylogeny and root on an in-tree outgroup label or taxon
  boldgenotyper data/Euprymna.tsv --build-tree --phylo-outgroup-label "Eusphyra blochii haplotype_h10_n21_n19"
  boldgenotyper data/Euprymna.tsv --build-tree --phylo-outgroup-taxon "Eusphyra blochii"

  # Include an external outgroup FASTA for rooting
  boldgenotyper data/Euprymna.tsv --build-tree --phylo-outgroup-fasta data/outgroup.fasta

  # Skip geographic plots or the HTML report
  boldgenotyper data/Euprymna.tsv --no-geo --no-report

  # Custom geographic analysis for freshwater organisms (HydroBASINS)
  boldgenotyper data/Salmonidae.tsv --build-tree \\
                --custom-shp shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev01-12_v1c.shp \\
                --shp-field HYBAS_ID --geo-category freshwater_basin

  # Custom geographic analysis for terrestrial organisms (Ecoregions2017)
  boldgenotyper data/Pieridae.tsv --build-tree \\
                --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp \\
                --shp-field ECO_NAME --geo-category ecoregion

Notes:
  - Phylogeny requires MAFFT and FastTree in PATH; trimAl is used if available.
  - Trees are built from haplotypes and relabeled with haplotype_sp by default.
  - Geographic analysis assigns samples to regions via point-in-polygon spatial joins.
  - Default: GOaS ocean basins (marine). Use --custom-shp for freshwater basins,
    terrestrial ecoregions, watersheds, biomes, or any polygon shapefile.
  - Plot regeneration kits are exported by default; use --no-export-plot-data to skip.
  - Reports go to the output folder (reports/, visualization/, phylogenetic/).
"""
    )

    # Required arguments
    parser.add_argument(
        'tsv',
        type=Path,
        help='Input BOLD TSV file with sequence and metadata'
    )

    # Optional arguments
    parser.add_argument(
        '-o', '--organism',
        type=str,
        default=None,
        help='Organism name (default: inferred from TSV filename)'
    )

    parser.add_argument(
        '--output', '--output-dir',
        type=Path,
        default=None,
        help='Output directory (default: {organism}_output in current directory)'
    )

    parser.add_argument(
        '--similarity-threshold',
        type=float,
        default=0.5,
        help='Minimum similarity for genotype assignment (default: 0.5)'
    )

    parser.add_argument(
        '--tie-margin',
        type=float,
        default=0.003,
        help='Maximum identity difference between top matches to call a tie. '
             'Samples with ambiguous assignments (identity difference < tie-margin) '
             'are flagged for manual review. '
             'Default: 0.003 (0.3%% difference)'
    )

    parser.add_argument(
        '--tie-threshold',
        type=float,
        default=0.95,
        help='Minimum identity required to consider tie detection. '
             'Ties are only called when best match identity >= tie-threshold. '
             'Prevents flagging low-quality matches as ties. '
             'Default: 0.95 (95%% identity)'
    )

    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of parallel threads (default: 4)'
    )

    parser.add_argument(
        '--build-tree',
        action='store_true',
        help='Build phylogenetic tree (requires MAFFT and FastTree)'
    )

    parser.add_argument(
        '--phylo-outgroup-fasta',
        type=Path,
        default=None,
        help='Optional FASTA of outgroup sequences to include for rooting the phylogeny'
    )

    parser.add_argument(
        '--phylo-outgroup-label',
        type=str,
        default=None,
        help='Optional tip label to reroot the tree on (when an outgroup sequence is already in the alignment)'
    )

    parser.add_argument(
        '--phylo-outgroup-taxon',
        type=str,
        default=None,
        help='Species/genus name to reroot by LCA of matching tips (uses haplotype taxonomy)'
    )

    parser.add_argument(
        '--no-export-plot-data',
        action='store_true',
        help='Skip exporting raw plot data and regeneration scripts'
    )

    parser.add_argument(
        '--export-format',
        nargs='+',
        default=['all'],
        choices=['arlequin', 'popart', 'dnasp', 'generic', 'all'],
        help='Population genetics formats to export (default: all). '
             'Options: arlequin (Arlequin .arp), popart (PopART/NEXUS), '
             'dnasp (DnaSP .fas), generic (CSV/FASTA), all. '
             'Populations in Arlequin and DnaSP are defined by --geo-category. '
             'Use --no-export-format to skip popgen export entirely. '
             'Example: --export-format arlequin popart'
    )
    parser.add_argument(
        '--no-export-format',
        action='store_true',
        help='Skip population genetics export entirely (overrides --export-format).'
    )

    parser.add_argument(
        '--no-report',
        action='store_true',
        help='Skip generating HTML summary report'
    )

    parser.add_argument(
        '--no-geo',
        action='store_true',
        help='Skip geographic region assignment and related visualizations. '
             'Use this if you only need genotyping and phylogeny without geographic analysis.'
    )

    parser.add_argument(
        '--custom-shp',
        type=Path,
        default=None,
        dest='shapefile_path',
        help='Path to custom shapefile (.shp) for geographic region assignment. '
             'Works with any polygon shapefile: ocean basins, freshwater basins, '
             'terrestrial ecoregions, watersheds, biomes, or user-defined regions. '
             'When provided, --geo-category defaults to "geographic_region" '
             '(override with --geo-category). Without this flag, the built-in '
             'GOaS marine ocean basin shapefile is used.'
    )

    parser.add_argument(
        '--shp-field',
        type=str,
        default='name',
        dest='shapefile_field',
        help='Name of shapefile attribute containing region labels (default: "name"). '
             'Examples: "ECO_NAME" for Ecoregions2017, "HYBAS_ID" for HydroBASINS. '
             'Use this to specify which field in your shapefile contains the region identifiers.'
    )

    parser.add_argument(
        '--geo-category',
        type=str,
        default=None,
        dest='geo_category',
        help='Name for geographic category in outputs. Defaults to "ocean_basin" '
             'when using GOaS (no --custom-shp), or "geographic_region" when using '
             'a custom shapefile. Examples: "ecoregion", "freshwater_basin", '
             '"watershed", "biome". This label is used in output files, plots, '
             'and column names.'
    )

    # Metadata Summarization Arguments (enabled by default).
    # Note: --no-metadata-analysis is retained as a deprecated alias of
    # --no-metadata-summary to keep older invocations working.
    parser.add_argument(
        '--no-metadata-summary',
        '--no-metadata-analysis',
        dest='no_metadata_summary',
        action='store_true',
        help='Disable metadata summarization (enabled by default). '
             'The legacy alias --no-metadata-analysis is also accepted.'
    )

    parser.add_argument(
        '--metadata-fields',
        nargs='+',
        default=None,
        help='Metadata fields to tabulate. Default: sex, life_stage, reproduction, '
             'country/ocean, country_iso, province/state, realm, biome, ecoregion, habitat, geoid'
    )

    parser.add_argument(
        '--no-normalize-sex',
        action='store_true',
        help='Do not normalize sex values. By default, sex values are normalized '
             '(M/Male/m -> Male, F/Female/f -> Female).'
    )

    parser.add_argument(
        '--no-temporal-analysis',
        action='store_true',
        help='Disable temporal analysis of collection dates (enabled by default)'
    )

    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging verbosity (default: INFO)'
    )

    parser.add_argument(
        '--singleton-distance',
        type=float,
        default=0.005,
        dest='singleton_distance',
        help='Minimum divergence threshold for retaining singleton haplotypes. '
             'Singleton ESVs diverging by less than this value from their nearest '
             'neighbour are removed as likely sequencing or PCR errors. '
             'Default: 0.005 (0.5%%). Use boldgenotyper-sweep to identify the '
             'optimal value for your dataset; the recommended threshold is reported '
             'in recommendations.txt and the elbow plot.'
    )

    parser.add_argument(
        '--keep-intermediates',
        action='store_true',
        help='Keep intermediate files (default: remove after pipeline completion)'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'BOLDGenotyper {__version__}'
    )

    args = parser.parse_args()

    # Validate input file
    if not args.tsv.exists():
        print(f"Error: Input TSV file not found: {args.tsv}", file=sys.stderr)
        return 1

    # Validate custom shapefile if provided
    if args.shapefile_path and not args.shapefile_path.exists():
        print(f"Error: Custom shapefile not found: {args.shapefile_path}", file=sys.stderr)
        return 1

    # Resolve geo_category default: 'geographic_region' with custom shp, 'ocean_basin' with GOaS
    if args.geo_category is None:
        args.geo_category = 'geographic_region' if args.shapefile_path else 'ocean_basin'

    # Determine organism name and sanitize for consistent file naming
    organism = args.organism if args.organism else extract_organism_from_path(args.tsv)
    # Always sanitize to ensure consistent filenames (replace spaces with underscores)
    organism = utils.sanitize_filename(organism)

    # Determine output directory
    output_dir = args.output if args.output else Path(f"{organism}_output")
    output_dir = output_dir.resolve()

    # Setup logging
    log_file = output_dir / f"{organism}_pipeline.log"
    output_dir.mkdir(parents=True, exist_ok=True)

    utils.setup_logging(log_level=args.log_level, log_file=str(log_file))

    # Load and configure pipeline
    cfg = config.get_default_config()
    cfg = cfg.update(
        genotype_assignment__min_identity=args.similarity_threshold,
        genotype_assignment__tie_margin=args.tie_margin,
        genotype_assignment__tie_threshold=args.tie_threshold,
        haplotype__min_singleton_distance=args.singleton_distance,
        n_threads=args.threads,
        output_dir=output_dir,
        log_level=args.log_level,
        phylogenetic__build_tree=args.build_tree,
        phylogenetic__outgroup_fasta=args.phylo_outgroup_fasta,
        phylogenetic__outgroup_label=args.phylo_outgroup_label,
        phylogenetic__outgroup_taxon=args.phylo_outgroup_taxon,
        keep_intermediates=args.keep_intermediates
    )

    # Print banner
    print("=" * 80)
    print("BOLDGenotyper Pipeline")
    print("=" * 80)
    print(f"Organism: {organism}")
    print(f"Input: {args.tsv}")
    print(f"Output: {output_dir}")
    print()
    print("Parameters:")
    print(f"  Similarity threshold: {args.similarity_threshold} ({args.similarity_threshold*100:.0f}% identity)")
    print(f"  Tie margin: {args.tie_margin} ({args.tie_margin*100:.1f}% difference)")
    print(f"  Tie threshold: {args.tie_threshold} ({args.tie_threshold*100:.0f}% identity)")
    print(f"  Singleton distance: {args.singleton_distance} ({args.singleton_distance*100:.1f}% — use boldgenotyper-sweep to optimise)")
    print(f"  Threads: {args.threads}")
    print(f"  Build tree: {args.build_tree}")
    print(f"  Phylo outgroup FASTA: {args.phylo_outgroup_fasta or 'None'}")
    print(f"  Phylo outgroup label: {args.phylo_outgroup_label or 'None'}")
    print(f"  Phylo outgroup taxon: {args.phylo_outgroup_taxon or 'None'}")
    if args.no_geo:
        print(f"  Geographic analysis: Disabled")
    elif args.shapefile_path:
        print(f"  Geographic analysis: Custom shapefile")
        print(f"    Shapefile: {args.shapefile_path}")
        print(f"    Field: {args.shapefile_field}")
        print(f"    Category: {args.geo_category}")
    else:
        print(f"  Geographic analysis: GOaS ocean basins (default)"
              f" — use --custom-shp for other regions")
    print("=" * 80)
    print()

    # Run pipeline
    try:
        success = run_pipeline(
            tsv_path=args.tsv,
            organism=organism,
            output_dir=output_dir,
            cfg=cfg,
            no_report=args.no_report,
            skip_geo=args.no_geo,
            export_plot_data=not args.no_export_plot_data,
            export_popgen_formats=(None if args.no_export_format else args.export_format),
            metadata_summary_enabled=not args.no_metadata_summary,
            metadata_fields=args.metadata_fields,
            normalize_sex=not args.no_normalize_sex,
            temporal_analysis=not args.no_temporal_analysis,
            shapefile_path=args.shapefile_path,
            shapefile_field=args.shapefile_field,
            geo_category=args.geo_category
        )

        return 0 if success else 1

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}", exc_info=True)
        print(f"\nError: Pipeline failed. Check log file: {log_file}", file=sys.stderr)
        return 1


def main_query(argv=None) -> int:
    """
    CLI entry point for the 'boldgenotyper-query' command.

    Query new COI sequences against previously identified haplotypes from
    a completed BOLD analysis. Enables reproducible haplotype assignment
    for new samples without re-running the full pipeline.

    Example:
        # Query single sequence
        boldgenotyper-query \
            --query new_sample.fasta \
            --haplotypes analysis/haplotypes/Organism_haplotypes.fasta \
            --output query_results/

        # Query with metadata enrichment
        boldgenotyper-query \
            --query new_samples.fasta \
            --haplotypes analysis/haplotypes/Organism_haplotypes.fasta \
            --analysis-dir analysis/ \
            --output query_results/
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper-query",
        description="Query new COI sequences against previously identified haplotypes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Query single sequence
  boldgenotyper-query --query new_sample.fasta \\
                       --haplotypes analysis/haplotypes/Organism_haplotypes.fasta

  # Query multiple sequences with metadata
  boldgenotyper-query --query samples.fasta \\
                       --haplotypes analysis/haplotypes/Organism_haplotypes.fasta \\
                       --analysis-dir analysis/ \\
                       --output results/

  # Adjust top matches and length filters
  boldgenotyper-query --query samples.fasta \\
                       --haplotypes haplotypes.fasta \\
                       --top-n 20 \\
                       --min-length 150 \\
                       --max-length 1500

Output Formats:
  All three formats are generated automatically:
  - query_results.csv: Machine-readable table
  - query_results.json: Structured data for programmatic access
  - query_results_detailed.txt: Human-readable report with alignments

Match Quality Classification:
  - perfect: 100% identity - exact haplotype match
  - high: ≥99.5% identity - likely same haplotype
  - good: ≥97% identity - same species, possibly different haplotype
  - moderate: ≥95% identity - same genus, divergent haplotype
  - low: <95% identity - different species or contamination

For more information:
  https://github.com/SymbioSeas/BOLDGenotyper
""")

    # Required arguments
    parser.add_argument(
        '--query',
        type=Path,
        required=True,
        help='Query FASTA file (single or multi-FASTA with COI sequences)'
    )

    parser.add_argument(
        '--haplotypes',
        type=Path,
        required=True,
        help='Haplotype consensus FASTA from previous analysis'
    )

    # Optional arguments
    parser.add_argument(
        '--analysis-dir',
        type=Path,
        default=None,
        help='Previous analysis directory (for metadata enrichment)'
    )

    parser.add_argument(
        '--output',
        '-o',
        type=Path,
        default=None,
        help='Output directory (default: haplotype_query_results/)'
    )

    parser.add_argument(
        '--top-n',
        type=int,
        default=10,
        help='Number of top matches to report per query (default: 10)'
    )

    parser.add_argument(
        '--min-length',
        type=int,
        default=100,
        help='Minimum query sequence length to process (default: 100)'
    )

    parser.add_argument(
        '--max-length',
        type=int,
        default=2000,
        help='Maximum query sequence length to process (default: 2000)'
    )

    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging verbosity (default: INFO)'
    )

    # Parse arguments
    args = parser.parse_args(argv)

    # Determine output directory
    output_dir = args.output if args.output else Path("haplotype_query_results")
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = output_dir / "query.log"
    utils.setup_logging(log_level=args.log_level, log_file=str(log_file))

    logger.info("=" * 80)
    logger.info("BOLDGenotyper Haplotype Query")
    logger.info("=" * 80)
    logger.info(f"Query file: {args.query}")
    logger.info(f"Haplotype file: {args.haplotypes}")
    if args.analysis_dir:
        logger.info(f"Analysis directory: {args.analysis_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("")

    try:
        # Import haplotype_query module
        from . import haplotype_query

        # Run query
        logger.info("Starting haplotype query analysis...")
        results, metadata = haplotype_query.query_against_haplotypes(
            query_fasta=args.query,
            haplotype_fasta=args.haplotypes,
            analysis_dir=args.analysis_dir,
            top_n=args.top_n,
            min_length=args.min_length,
            max_length=args.max_length
        )

        # Write results
        logger.info("Writing results...")
        haplotype_query.write_results(
            results=results,
            output_dir=output_dir,
            metadata=metadata,
            haplotype_file=args.haplotypes,
            analysis_dir=args.analysis_dir
        )

        logger.info("")
        logger.info("=" * 80)
        logger.info("Query Analysis Complete!")
        logger.info("=" * 80)
        logger.info(f"Results written to: {output_dir}")
        logger.info("")
        logger.info("Output files:")
        logger.info(f"  - {output_dir / 'query_results.csv'}")
        logger.info(f"  - {output_dir / 'query_results.json'}")
        logger.info(f"  - {output_dir / 'query_results_detailed.txt'}")
        logger.info("")

        return 0

    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        logger.error(f"Query analysis failed: {e}", exc_info=True)
        print(f"Error: Query analysis failed. Check log: {log_file}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
