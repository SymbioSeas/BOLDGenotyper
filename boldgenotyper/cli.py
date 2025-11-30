#!/usr/bin/env python3
"""
BOLDGenotyper Command-Line Interface

Unified pipeline for BOLD barcode data: genotype assignment, geographic analysis,
phylogenetic reconstruction, and visualization.

Author: Steph Smith (steph.smith@unc.edu)
"""

import argparse
import sys
import os
import logging
import json
from pathlib import Path
from typing import Optional
import pandas as pd

# Local imports
from . import (
    utils, config, metadata, geographic, dereplication,
    genotype_assignment, phylogenetics, visualization, reports,
    cluster_diagnostics, quality_control, plot_export, comparative_analysis,
    divergence_analysis, parameter_sweep, geographic_enhancement,
    metadata_enrichment, popgen_export,
)

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
    """Create organized output directory structure."""
    dirs = {
        'base': base_output,
        'intermediate': base_output / 'intermediate',
        'dereplication': base_output / 'intermediate' / 'dereplication',
        'intermediate_phylo': base_output / 'intermediate' / 'phylogenetic',
        'intermediate_assignments': base_output / 'intermediate' / 'genotype_assignments',
        'intermediate_geographic': base_output / 'intermediate' / 'geographic',
        'assignments': base_output / 'genotype_assignments',
        'taxonomy': base_output / 'taxonomy',
        'quality_control': base_output / 'quality_control',
        'divergence_analysis': base_output / 'divergence_analysis',
        'geographic_analysis': base_output / 'geographic_analysis',
        'phylogenetic': base_output / 'phylogenetic',
        'visualization': base_output / 'visualization',
        'genotype_plots': base_output / 'visualization' / 'genotype_plots',
        'reports': base_output / 'reports',
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def run_pipeline(
    tsv_path: Path,
    organism: str,
    output_dir: Path,
    cfg: config.PipelineConfig,
    no_report: bool = False,
    skip_geo: bool = False,
    export_plot_data: bool = False,
    export_popgen_formats: Optional[list] = None,
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
        Export raw plot data and R regeneration scripts (default: False)
    export_popgen_formats : list, optional
        List of population genetics formats to export (default: None)

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
    import json
    params = {
        'clustering_threshold': cfg.dereplication.clustering_threshold,
        'similarity_threshold': cfg.genotype_assignment.min_identity,
        'tie_margin': cfg.genotype_assignment.tie_margin,
        'tie_threshold': cfg.genotype_assignment.tie_threshold,
        'threads': cfg.n_threads,
        'build_tree': cfg.phylogenetic.build_tree
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

        # Assign ocean basins (only to samples with geographic quality)
        logger.info("1.3: Assigning ocean basins...")
        geo_analysis_performed = False  # Track whether geographic analysis was successful

        if skip_geo:
            logger.info("  ⊘ Geographic analysis disabled (--no-geo flag)")
            df_with_basins = df.copy()
            df_with_basins['ocean_basin'] = 'Unknown'
        elif not cfg.geographic.goas_shapefile_path.exists():
            logger.warning(f"  ⚠ GOaS shapefile not found at: {cfg.geographic.goas_shapefile_path}")
            logger.warning("")
            logger.warning("  To enable geographic analysis, download the GOaS shapefile:")
            logger.warning("  1. Run: python -m boldgenotyper.goas_downloader")
            logger.warning("  2. Or download manually from: https://www.marineregions.org/download_file.php?name=World_Seas_IHO_v3.zip")
            logger.warning(f"  3. Extract to: {cfg.geographic.goas_shapefile_path.parent}")
            logger.warning("  4. Ensure the .shp file is named: goas_v01.shp")
            logger.warning("")
            logger.warning("  Pipeline will continue without geographic analysis...")
            df_with_basins = df.copy()
            df_with_basins['ocean_basin'] = 'Unknown'
        else:
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
                    df_with_basins['ocean_basin'] = 'Unknown'

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

                    # Update the columns
                    df_with_basins['ocean_basin'] = merged_geo['ocean_basin'].fillna('Unknown')
                    if 'lat_new' in merged_geo.columns:
                        df_with_basins['lat'] = merged_geo['lat_new']
                    elif 'lat' in merged_geo.columns:
                        df_with_basins['lat'] = merged_geo['lat']
                    if 'lon_new' in merged_geo.columns:
                        df_with_basins['lon'] = merged_geo['lon_new']
                    elif 'lon' in merged_geo.columns:
                        df_with_basins['lon'] = merged_geo['lon']

                    n_assigned = (df_with_basins['ocean_basin'] != 'Unknown').sum()
                    logger.info(f"  ✓ Assigned ocean basins to {n_assigned} samples")
                    logger.info(f"  → {len(df) - n_assigned} samples marked as 'Unknown' (centroid/missing coordinates)")
                    geo_analysis_performed = True
                else:
                    logger.warning("  ⚠ No samples with geographic-quality coordinates")
                    df_with_basins = df.copy()
                    df_with_basins['ocean_basin'] = 'Unknown'

            except Exception as e:
                logger.error(f"  ✗ Failed to load GOaS data: {e}")
                logger.warning("  Pipeline will continue without geographic analysis...")
                df_with_basins = df.copy()
                df_with_basins['ocean_basin'] = 'Unknown'

        basins_tsv = dirs['intermediate_geographic'] / "samples_with_ocean_basins.tsv"
        df_with_basins.to_csv(basins_tsv, sep='\t', index=False)

    except Exception as e:
        logger.error(f"Phase 1 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 2: Sequence Dereplication and Consensus Generation
    # ========================================================================
    logger.info("")
    logger.info("PHASE 2: Sequence Dereplication and Consensus Generation")
    logger.info("-" * 80)

    try:
        # Generate FASTA from ALL samples (including those with centroid coordinates)
        logger.info("2.1: Generating FASTA from sequences...")
        logger.info(f"  Using all {len(df_with_basins)} samples (including those with centroid coordinates)")
        logger.info(f"  Minimum sequence length: {cfg.dereplication.min_sequence_length} bp")
        fasta_records = []
        skipped_count = 0
        skipped_reasons = {}
        for _, row in df_with_basins.iterrows():
            header = f"{organism}_{row['processid']}.COI-5P"
            sequence = row['nuc']

            # Skip if sequence is missing or not a string
            if not isinstance(sequence, str) or not sequence.strip():
                skipped_count += 1
                skipped_reasons['missing_or_empty'] = skipped_reasons.get('missing_or_empty', 0) + 1
                continue

            is_valid, reason = utils.validate_sequence(
                sequence,
                min_length=cfg.dereplication.min_sequence_length
            )
            if is_valid:
                fasta_records.append((header, sequence))
            else:
                skipped_count += 1
                skipped_reasons[reason] = skipped_reasons.get(reason, 0) + 1

        if skipped_count > 0:
            logger.warning(f"  ⚠ Skipped {skipped_count} samples with missing or invalid sequences")
            for reason, count in skipped_reasons.items():
                logger.debug(f"    - {reason}: {count}")

        fasta_path = dirs['intermediate'] / f"{organism}.fasta"
        utils.write_fasta(fasta_records, fasta_path)
        logger.info(f"  ✓ Created {len(fasta_records)} FASTA records")

        # Dereplicate sequences
        logger.info("2.2: Dereplicating sequences (clustering at 99% identity)...")
        consensus_records = dereplication.dereplicate_from_fasta(
            input_fasta=str(fasta_path),
            output_dir=str(dirs['dereplication']),
            threshold=cfg.dereplication.clustering_threshold,
            frequency_cutoff=cfg.dereplication.consensus_frequency_cutoff,
            min_post_trim_length=cfg.dereplication.min_post_trim_length,
            min_consensus_length_ratio=cfg.dereplication.min_consensus_length_ratio
        )
        logger.info(f"  ✓ Identified {len(consensus_records)} unique genotypes")

        # Consensus sequences are saved in intermediate/dereplication/
        consensus_path = dirs['dereplication'] / f"{organism}_consensus.fasta"

    except Exception as e:
        logger.error(f"Phase 2 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 3: Genotype Assignment
    # ========================================================================
    logger.info("")
    logger.info("PHASE 3: Genotype Assignment")
    logger.info("-" * 80)

    try:
        logger.info("3.1: Assigning samples to genotypes...")
        logger.info(f"  Using all {len(df_with_basins)} samples (including those with centroid coordinates)")
        annotated_tsv = dirs['intermediate_assignments'] / f"{organism}_with_genotypes.tsv"
        diagnostics_csv = dirs['assignments'] / f"{organism}_diagnostics.csv"

        stats = genotype_assignment.assign_genotypes(
            metadata_path=str(basins_tsv),  # Use full dataset with ocean basins
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(annotated_tsv),
            min_identity=cfg.genotype_assignment.min_identity,
            n_processes=cfg.n_threads,
            diagnostics_path=str(diagnostics_csv)
        )

        logger.info(f"  ✓ Assigned {stats['assigned']}/{stats['total']} samples to genotypes")
        logger.info(f"  ✓ Assignment rate: {stats['assigned']/stats['total']*100:.1f}%")

        # Load annotated data
        df_with_genotypes = pd.read_csv(annotated_tsv, sep="\t")

    except Exception as e:
        logger.error(f"Phase 3 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 4: Taxonomy Assignment
    # ========================================================================
    logger.info("")
    logger.info("PHASE 4: Taxonomy Assignment to Consensus Groups")
    logger.info("-" * 80)

    try:
        logger.info("4.1: Assigning taxonomy to consensus groups...")
        assign_table, species_counts = utils.assign_consensus_taxonomy(
            df_with_genotypes,
            group_col="consensus_group",
            species_col="species",
            genus_col="genus",
            majority_threshold=cfg.taxonomy.majority_species_threshold
        )

        # Build consensus_group_sp labels
        def _strip_prefix(s: str) -> str:
            return s.replace("consensus_", "") if isinstance(s, str) else s

        assign_table["consensus_group_short"] = assign_table["consensus_group"].map(_strip_prefix)

        def _join_label(row):
            sp = row["assigned_sp"]
            short = row["consensus_group_short"]
            if not sp:
                return short
            return f"{sp} {short}"

        assign_table["consensus_group_sp"] = assign_table.apply(_join_label, axis=1)

        # Save taxonomy files
        species_counts_out = dirs['taxonomy'] / f"{organism}_species_by_consensus.csv"
        assign_table_out = dirs['taxonomy'] / f"{organism}_consensus_taxonomy.csv"
        species_counts.to_csv(species_counts_out, index=False)
        assign_table.to_csv(assign_table_out, index=False)

        # Merge taxonomy back into main dataframe
        df_with_genotypes = df_with_genotypes.merge(
            assign_table[["consensus_group", "assigned_sp", "consensus_group_sp", "assignment_level", "assignment_notes", "majority_fraction"]],
            on="consensus_group",
            how="left",
            validate="many_to_one"
        )

        # Handle cluster sequence taxonomy if available
        cluster_seq_path = dirs['taxonomy'] / f"{organism}_consensus_taxonomy_seq.csv"
        if cluster_seq_path.exists():
            cluster_seq_df = pd.read_csv(cluster_seq_path)
            df_with_genotypes = df_with_genotypes.merge(
                cluster_seq_df, on="consensus_group", how="left", validate="many_to_one"
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
                short = row["consensus_group"].replace("consensus_", "") if isinstance(row.get("consensus_group"), str) else ""
                label = f"{final_sp} {short}".strip() if final_sp else short
                return pd.Series({
                    "final_group_sp": final_sp,
                    "final_group_level": final_level,
                    "tax_provenance": prov,
                    "consensus_group_sp": label
                })

            final_cols = df_with_genotypes.apply(_final_label, axis=1)
            for c in final_cols.columns:
                df_with_genotypes[c] = final_cols[c]

        logger.info(f"  ✓ Assigned taxonomy to {len(assign_table)} consensus groups")

        # Merge with geographic data (only if not already present in df_with_genotypes)
        # df_with_genotypes may already have lat/lon/ocean_basin from the TSV file
        geo_cols_to_merge = []
        for col in ['lat', 'lon', 'ocean_basin']:
            if col in df_with_basins.columns and col not in df_with_genotypes.columns:
                geo_cols_to_merge.append(col)

        if geo_cols_to_merge:
            # Only merge if there are columns to add
            geo_keep = ['processid'] + geo_cols_to_merge
            df_final = df_with_genotypes.merge(df_with_basins[geo_keep], on='processid', how='left', validate='one_to_one')
        else:
            # df_with_genotypes already has all geographic data
            df_final = df_with_genotypes.copy()

        # Save final annotated file with proper CSV quoting to handle fields with commas
        annotated_csv = dirs['base'] / f"{organism}_annotated.csv"
        df_final.to_csv(annotated_csv, index=False, quoting=1)  # quoting=1 is csv.QUOTE_MINIMAL

        # Verify file was created
        if annotated_csv.exists():
            logger.info(f"  ✓ Saved annotated dataset: {annotated_csv}")
        else:
            logger.warning(f"  ⚠ Failed to save annotated dataset: {annotated_csv}")

    except Exception as e:
        logger.error(f"Phase 4 failed: {e}", exc_info=True)
        return False

    # ========================================================================
    # PHASE 4.5: Quality Control and Contamination Detection
    # ========================================================================
    logger.info("")
    logger.info("PHASE 4.5: Quality Control and Contamination Detection")
    logger.info("-" * 80)

    try:
        logger.info("4.5.1: Adding contamination analysis columns to main output...")

        # Add contamination columns to the dataframe
        df_final_enhanced = quality_control.add_contamination_columns(
            df_final,
            group_col="consensus_group",
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

        logger.info("4.5.2: Generating quality control reports...")

        # Generate comprehensive QC report
        qc_results = quality_control.generate_quality_control_report(
            df_final,
            output_dir=dirs['quality_control'],
            organism=organism,
            group_col="consensus_group",
            species_col="species"
        )

        logger.info(f"  ✓ Quality control reports saved to: {dirs['quality_control']}")

        # Print quality alerts to console
        quality_control.print_quality_alert(qc_results)

    except Exception as e:
        logger.warning(f"Quality control analysis failed (non-critical): {e}")
        logger.debug("QC error details:", exc_info=True)

    # ========================================================================
    # PHASE 4.7: Geographic Analysis Enhancement
    # ========================================================================
    if geo_analysis_performed:
        logger.info("")
        logger.info("PHASE 4.7: Geographic Analysis Enhancement")
        logger.info("-" * 80)

        try:
            logger.info("4.7.1: Assessing geographic data coverage and quality...")

            # Run geographic enhancement
            geo_enhancement_results = geographic_enhancement.enhance_geographic_analysis(
                df=df_final,
                output_dir=dirs['geographic_analysis'],
                organism=organism,
                goas_data=goas_data if 'goas_data' in locals() else None
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
    # PHASE 5: Phylogenetic Analysis (Optional)
    # ========================================================================
    logger.info("")
    logger.info("PHASE 5: Phylogenetic Analysis")
    logger.info("-" * 80)

    tree_path = None
    if cfg.phylogenetic.build_tree:
        try:
            logger.info("5.1: Building phylogenetic tree...")

            # Check for required tools
            if not utils.check_external_tool("mafft"):
                logger.warning("  ⚠ MAFFT not found, skipping tree building")
            elif not utils.check_external_tool("fasttree"):
                logger.warning("  ⚠ FastTree not found, skipping tree building")
            else:
                # Build phylogenetic tree
                # Save alignment files to intermediate, tree files to final phylogenetic directory
                intermediate_prefix = dirs['intermediate_phylo'] / organism
                output_prefix = dirs['phylogenetic'] / organism

                tree = phylogenetics.build_phylogeny(
                    consensus_fasta=str(consensus_path),
                    output_prefix=str(intermediate_prefix),
                    threads=cfg.n_threads,
                    min_consensus_length=cfg.phylogenetic.min_consensus_length,
                    min_cluster_size=cfg.phylogenetic.min_cluster_size,
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
                    logger.info(f"  ✓ Built phylogenetic tree: {tree_path}")

                    # Create relabeled versions with consensus_group_sp labels
                    logger.info("5.2: Creating relabeled tree and alignment files...")
                    try:
                        alignment_path = f"{intermediate_prefix}_aligned.fasta"
                        taxonomy_csv_path = dirs['taxonomy'] / f"{organism}_consensus_taxonomy.csv"

                        if Path(alignment_path).exists() and taxonomy_csv_path.exists():
                            relabeled_tree_path = f"{output_prefix}_tree_relabeled.nwk"
                            relabeled_alignment_path = f"{intermediate_prefix}_aligned_relabeled.fasta"

                            phylogenetics.relabel_tree_and_alignment(
                                tree_file=tree_path,
                                alignment_file=alignment_path,
                                taxonomy_csv=str(taxonomy_csv_path),
                                output_tree=relabeled_tree_path,
                                output_alignment=relabeled_alignment_path,
                                label_column="consensus_group_sp",
                                id_column="consensus_group",
                            )
                            logger.info(f"  ✓ Created relabeled tree: {relabeled_tree_path}")
                            logger.info(f"  ✓ Created relabeled alignment (intermediate): {relabeled_alignment_path}")
                        else:
                            logger.warning("  ⚠ Skipping relabeling: alignment or taxonomy file not found")
                    except Exception as e:
                        logger.warning(f"  ⚠ Relabeling failed (non-critical): {e}")
                else:
                    logger.warning(f"  ⚠ Phylogenetic tree building completed but output file not found: {tree_path}")
                    tree_path = None

        except Exception as e:
            logger.warning(f"Phylogenetic analysis failed (non-critical): {e}")
    else:
        logger.info("  ⊘ Phylogenetic tree building disabled in config")

    # ========================================================================
    # PHASE 5.5: Divergence Analysis
    # ========================================================================
    logger.info("")
    logger.info("PHASE 5.5: Divergence Analysis")
    logger.info("-" * 80)

    try:
        logger.info("5.5.1: Calculating pairwise divergence between consensus groups...")

        # Check if we have the required files
        taxonomy_csv_path = dirs['taxonomy'] / f"{organism}_consensus_taxonomy.csv"

        if consensus_path.exists() and taxonomy_csv_path.exists():
            # Run divergence analysis
            divergence_results = divergence_analysis.generate_divergence_analysis(
                consensus_fasta=consensus_path,
                taxonomy_csv=taxonomy_csv_path,
                output_dir=dirs['divergence_analysis'],
                organism=organism
            )

            logger.info(f"  ✓ Divergence analysis complete: {dirs['divergence_analysis']}")

            # Print summary to console
            divergence_analysis.print_divergence_summary(divergence_results)

        else:
            logger.warning("  ⊘ Skipping divergence analysis (consensus or taxonomy files not found)")

    except Exception as e:
        logger.warning(f"Divergence analysis failed (non-critical): {e}")
        logger.debug("Divergence analysis error details:", exc_info=True)

    # ========================================================================
    # PHASE 6: Visualization
    # ========================================================================
    logger.info("")
    logger.info("PHASE 6: Visualization")
    logger.info("-" * 80)

    try:
        logger.info("6.1: Generating plots...")

        # Generate visualizations for each format
        for fmt in cfg.visualization.figure_format:
            # Geographic visualizations - only if geographic analysis was performed
            if geo_analysis_performed:
                # Distribution maps
                if 'lat' in df_final.columns and 'lon' in df_final.columns:
                    try:
                        map_path, map_data = visualization.plot_distribution_map(
                            df=df_final,
                            output_path=str(dirs['visualization'] / f"{organism}_distribution_map.{fmt}"),
                            genotype_column='consensus_group_sp',
                            latitude_col='lat',
                            longitude_col='lon'
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if map_data:
                            json_path = dirs['visualization'] / f"{organism}_distribution_map_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(map_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Distribution map skipped: {e}")

                # Ocean basin abundance bar plot (relative)
                if 'ocean_basin' in df_final.columns and 'consensus_group_sp' in df_final.columns:
                    try:
                        bar_path, bar_data = visualization.plot_ocean_basin_abundance(
                            df=df_final,
                            output_path=str(dirs['visualization'] / f"{organism}_distribution_bar.{fmt}"),
                            genotype_column='consensus_group_sp',
                            basin_column='ocean_basin'
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if bar_data:
                            json_path = dirs['visualization'] / f"{organism}_distribution_bar_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(bar_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Ocean basin bar plot skipped: {e}")

                # Ocean basin abundance bar plot (total counts)
                if 'ocean_basin' in df_final.columns and 'consensus_group_sp' in df_final.columns:
                    try:
                        total_bar_path, total_bar_data = visualization.plot_ocean_basin_abundance_total(
                            df=df_final,
                            output_path=str(dirs['visualization'] / f"{organism}_totaldistribution_bar.{fmt}"),
                            genotype_column='consensus_group_sp',
                            basin_column='ocean_basin'
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if total_bar_data:
                            json_path = dirs['visualization'] / f"{organism}_totaldistribution_bar_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(total_bar_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Ocean basin total abundance bar plot skipped: {e}")

                # Faceted distribution map by consensus_group_sp
                if ('lat' in df_final.columns and 'lon' in df_final.columns and
                    'consensus_group_sp' in df_final.columns and 'consensus_group' in df_final.columns):
                    try:
                        visualization.plot_distribution_map_faceted(
                            df=df_final,
                            output_path=str(dirs['visualization'] / f"{organism}_distribution_map_faceted.{fmt}"),
                            genotype_column='consensus_group',
                            species_column='consensus_group_sp',
                            latitude_col='lat',
                            longitude_col='lon',
                            facet_by=cfg.visualization.facet_by,
                            map_buffer_degrees=cfg.visualization.map_buffer_degrees,
                            show_unknown_annotation=cfg.visualization.show_unknown_geography_annotation,
                            show_scale_bar=cfg.visualization.show_scale_bar,
                            genotype_plots_dir=dirs['genotype_plots'],
                            formats=cfg.visualization.figure_format
                        )
                    except Exception as e:
                        logger.warning(f"Faceted distribution map generation failed: {e}", exc_info=True)

                # Faceted ocean basin bar plot by consensus_group_sp
                if ('ocean_basin' in df_final.columns and 'consensus_group_sp' in df_final.columns and
                    'consensus_group' in df_final.columns):
                    try:
                        faceted_bar_path, faceted_bar_data = visualization.plot_ocean_basin_abundance_faceted(
                            df=df_final,
                            output_path=str(dirs['visualization'] / f"{organism}_distribution_bar_faceted.{fmt}"),
                            genotype_column='consensus_group',
                            species_column='consensus_group_sp',
                            basin_column='ocean_basin',
                            facet_by=cfg.visualization.facet_by,
                            genotype_plots_dir=dirs['genotype_plots'],
                            formats=cfg.visualization.figure_format
                        )
                        # Save plot data as JSON for interactive plotting in HTML report
                        if faceted_bar_data:
                            json_path = dirs['visualization'] / f"{organism}_distribution_bar_faceted_data.json"
                            with open(json_path, 'w') as f:
                                json.dump(faceted_bar_data, f, indent=2)
                            logger.debug(f"Saved plot data to: {json_path}")
                    except Exception as e:
                        logger.debug(f"Faceted basin bar plot skipped: {e}")
            else:
                logger.info("  ⊘ Skipping geographic visualizations (geographic analysis not performed)")

            # Identity distribution (always generated if diagnostics exist)
            if diagnostics_csv.exists():
                visualization.plot_identity_distribution(
                    diagnostics_csv=str(diagnostics_csv),
                    output_path=str(dirs['assignments'] / f"{organism}_identity_distribution.{fmt}"),
                    figsize=(10, 6),
                    dpi=cfg.visualization.figure_dpi
                )

            # Phylogenetic tree
            if tree_path and Path(tree_path).exists():
                # Use relabeled tree if it exists (so tips show consensus_group_sp labels)
                relabeled_tree_path = tree_path.replace("_tree.nwk", "_tree_relabeled.nwk")
                if Path(relabeled_tree_path).exists():
                    tree_to_plot = relabeled_tree_path
                    logger.info(f"Using relabeled tree for visualization: {relabeled_tree_path}")
                else:
                    tree_to_plot = tree_path

                # Load color map if it exists
                color_map_path = dirs['assignments'] / f"{organism}_genotype_color_map.csv"
                genotype_colors = None
                if color_map_path.exists():
                    color_df = pd.read_csv(color_map_path)
                    genotype_colors = dict(zip(color_df['consensus_group_sp'], color_df['color']))

                visualization.plot_phylogenetic_tree(
                    tree_file=str(tree_to_plot),
                    output_path=str(dirs['phylogenetic'] / f"{organism}_tree.{fmt}"),
                    genotype_colors=genotype_colors,
                    show_bootstrap=True,
                    bootstrap_threshold=cfg.visualization.show_bootstrap_threshold,
                    figsize=None,  # Auto-scale based on tree size
                    dpi=cfg.visualization.figure_dpi
                )

        logger.info(f"  ✓ Generated visualization plots")

    except Exception as e:
        logger.warning(f"Visualization generation encountered errors (non-critical): {e}")

    # ========================================================================
    # PHASE 6.5: Plot Data Export (Optional)
    # ========================================================================
    if export_plot_data:
        logger.info("")
        logger.info("PHASE 6.5: Plot Data Export")
        logger.info("-" * 80)

        try:
            # Get color map if it exists
            color_map_path = dirs['assignments'] / f"{organism}_genotype_color_map.csv"
            genotype_colors = None
            if color_map_path.exists():
                color_df = pd.read_csv(color_map_path)
                genotype_colors = dict(zip(color_df['consensus_group_sp'], color_df['color']))

            # Export plot data and regeneration kit
            plot_export_results = plot_export.export_plots_complete(
                df=df_final,
                output_dir=output_dir,
                organism=organism,
                consensus_path=consensus_path if consensus_path.exists() else None,
                tree_path=Path(tree_path) if tree_path and Path(tree_path).exists() else None,
                diagnostics_path=diagnostics_csv if diagnostics_csv.exists() else None,
                color_map=genotype_colors
            )

            logger.info(f"  ✓ Plot regeneration kit exported to {output_dir / 'plots'}")

        except Exception as e:
            logger.warning(f"Plot data export failed (non-critical): {e}")
            logger.debug("Plot export error details:", exc_info=True)

    # ========================================================================
    # PHASE 6.6: Population Genetics Export (Optional)
    # ========================================================================
    if export_popgen_formats:
        logger.info("")
        logger.info("PHASE 6.6: Population Genetics Export")
        logger.info("-" * 80)

        try:
            # Export to population genetics software formats
            popgen_results = popgen_export.export_population_genetics_formats(
                df=df_final,
                consensus_fasta_path=consensus_path,
                output_dir=output_dir,
                organism=organism,
                formats=export_popgen_formats,
                group_by='consensus_group'
            )

            logger.info(f"  ✓ Population genetics formats exported to {output_dir / 'exports'}")

        except Exception as e:
            logger.warning(f"Population genetics export failed (non-critical): {e}")
            logger.debug("Popgen export error details:", exc_info=True)

    # ========================================================================
    # PHASE 7: Reports
    # ========================================================================
    logger.info("")
    logger.info("PHASE 7: Generating Reports")
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
                version="1.0.0"
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
    # Pipeline Complete
    # ========================================================================
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"✓ Pipeline completed successfully for {organism}")
    logger.info(f"  Output directory: {output_dir}")
    logger.info(f"  Key files:")
    logger.info(f"    - Annotated data: {annotated_csv}")
    logger.info(f"    - Consensus sequences: {consensus_path}")
    if tree_path:
        logger.info(f"    - Phylogenetic tree: {tree_path}")
    if html_report_path:
        logger.info(f"    - HTML summary report: {html_report_path}")
    logger.info("=" * 80)

    return True

def main_cluster_diagnostics(argv=None) -> int:
    """
    CLI entry point for the 'cluster-diagnostics' subcommand.

    Example:
        boldgenotyper cluster-diagnostics \
            --consensus data/Sphyrnidae_consensus.fasta \
            --diagnostics data/Sphyrnidae_assignment_diagnostics.tsv \
            --alignment data/Sphyrnidae_trimmed_alignment.fasta \
            --threshold 0.01 \
            --output-dir diagnostics/Sphyrnidae
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper cluster-diagnostics",
        description="Cluster-level diagnostics for dereplication and genotype assignment",
    )

    parser.add_argument(
        "--consensus",
        required=True,
        type=Path,
        help="Path to consensus FASTA (e.g., consensus_Sphyrnidae.fasta)",
    )
    parser.add_argument(
        "--diagnostics",
        required=True,
        type=Path,
        help="Path to genotype assignment diagnostics TSV",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to write diagnostics outputs",
    )
    parser.add_argument(
        "--alignment",
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional trimmed alignment FASTA used in dereplication. "
            "If provided, intra-cluster distance stats and a dendrogram "
            "will be generated."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help=(
            "Distance threshold used during dereplication clustering "
            "(must match your dereplication threshold). Default: 0.01"
        ),
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity (default: INFO)",
    )

    args = parser.parse_args(argv)

    # Resolve and create output directory
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging to a dedicated diagnostics log file
    log_file = output_dir / "cluster_diagnostics.log"
    utils.setup_logging(log_level=args.log_level, log_file=str(log_file))

    logger.info("=== Cluster diagnostics ===")
    logger.info(f"Consensus FASTA: {args.consensus}")
    logger.info(f"Diagnostics TSV: {args.diagnostics}")
    if args.alignment:
        logger.info(f"Alignment FASTA: {args.alignment}")
        logger.info(f"Dereplication threshold: {args.threshold}")
    logger.info(f"Output directory: {output_dir}")

    # Load inputs via cluster_diagnostics helpers
    consensus_df = cluster_diagnostics.parse_consensus_fasta(args.consensus)
    diag_df = cluster_diagnostics.load_diagnostics(args.diagnostics)

    # Per-cluster assignment statistics
    stats_df = cluster_diagnostics.compute_cluster_assignment_stats(
        consensus_df, diag_df
    )

    # Optional intra-cluster distances + dendrogram
    if args.alignment:
        (
            intra_df,
            dist_condensed,
            seq_ids,
            cluster_labels,
        ) = cluster_diagnostics.compute_intracluster_distance_stats(
            alignment_path=args.alignment,
            threshold=args.threshold,
        )

        # Merge intra-cluster stats into main table
        stats_df = stats_df.merge(
            intra_df,
            on="cluster_id",
            how="left",
            suffixes=("", "_intra"),
        )

        # Global dendrogram
        cluster_diagnostics.plot_global_dendrogram(
            dist_condensed=dist_condensed,
            seq_ids=seq_ids,
            output_dir=output_dir,
            title="Global sequence dendrogram (trimmed alignment)",
        )

    # Write summary table
    out_tsv = output_dir / "cluster_diagnostics.tsv"
    stats_df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Wrote cluster diagnostics to {out_tsv}")

    logger.info("Cluster diagnostics completed.")
    return 0

def main_compare(argv=None) -> int:
    """
    CLI entry point for the 'boldgenotyper-compare' command.

    Compare species-level and family-level analyses to detect contamination.

    Example:
        boldgenotyper-compare \
            --species-level Sphyrna_lewini_output/ \
            --family-level Sphyrnidae_output/ \
            --output comparative_analysis/ \
            --generate-reassignment-table
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper-compare",
        description="Compare species-level and family-level analyses for contamination detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two completed analyses
  boldgenotyper-compare --species-level Sphyrna_lewini_output/ \
                         --family-level Sphyrnidae_output/

  # Generate sample reassignment table
  boldgenotyper-compare --species-level species_analysis/ \
                         --family-level family_analysis/ \
                         --generate-reassignment-table

  # Specify custom output directory
  boldgenotyper-compare --species-level sp/ --family-level fam/ \
                         --output custom_comparison/

For more information: https://github.com/your-repo/boldgenotyper
        """
    )

    parser.add_argument(
        '--species-level',
        type=Path,
        required=True,
        help='Path to species-level analysis directory or annotated CSV'
    )

    parser.add_argument(
        '--family-level',
        type=Path,
        required=True,
        help='Path to family-level analysis directory or annotated CSV'
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=Path('comparative_analysis'),
        help='Output directory for comparison results (default: comparative_analysis/)'
    )

    parser.add_argument(
        '--generate-reassignment-table',
        action='store_true',
        help='Generate sample-level reassignment table (Table S4 equivalent)'
    )

    parser.add_argument(
        '--majority-threshold',
        type=float,
        default=0.7,
        help='Threshold for species-level assignment (default: 0.7)'
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

    # Validate inputs
    if not args.species_level.exists():
        print(f"Error: Species-level path not found: {args.species_level}", file=sys.stderr)
        return 1

    if not args.family_level.exists():
        print(f"Error: Family-level path not found: {args.family_level}", file=sys.stderr)
        return 1

    # Run comparison
    try:
        results = comparative_analysis.compare_analyses(
            species_path=args.species_level,
            family_path=args.family_level,
            output_dir=args.output,
            generate_reassignment_table=args.generate_reassignment_table,
            majority_threshold=args.majority_threshold
        )

        print("\n✓ Comparative analysis complete!")
        print(f"  Results saved to: {args.output}")
        print("\nGenerated files:")
        for key, path in results.items():
            if isinstance(path, Path):
                print(f"  - {path.name}")

        return 0

    except Exception as e:
        logger.error(f"Comparative analysis failed: {e}", exc_info=True)
        print(f"\nError: {e}", file=sys.stderr)
        return 1


def main_sweep(argv=None) -> int:
    """
    CLI entry point for the 'boldgenotyper-sweep' command.

    Test multiple parameter values to optimize clustering thresholds.

    Example:
        boldgenotyper-sweep Sphyrna_lewini.tsv \
            --thresholds 0.01,0.015,0.02,0.03,0.05 \
            --output parameter_sweep/
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper-sweep",
        description="Parameter sweep tool for threshold optimization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test default threshold range
  boldgenotyper-sweep data/Sphyrna_lewini.tsv

  # Test custom thresholds
  boldgenotyper-sweep data/samples.tsv --thresholds 0.005,0.01,0.015,0.02

  # Run in parallel with 8 threads
  boldgenotyper-sweep data/samples.tsv --threads 8

For more information: https://github.com/your-repo/boldgenotyper
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


def main_enrich(argv=None) -> int:
    """
    CLI entry point for the 'boldgenotyper-enrich' command.

    Add custom metadata or update geographic assignments in BOLDGenotyper output.

    Example:
        boldgenotyper-enrich Sphyrnidae_annotated.csv \
            --add-metadata my_sampling_data.csv \
            --join-on processid \
            --output enriched_analysis/
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper-enrich",
        description="Add custom metadata or update geographic assignments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add custom metadata
  boldgenotyper-enrich Sphyrnidae_annotated.csv \
                       --add-metadata sampling_data.csv \
                       --output enriched_analysis/

  # Update geographic regions with custom shapefile (e.g., freshwater basins)
  boldgenotyper-enrich Salmo_trutta_annotated.csv \
                       --custom-shp freshwater_basins.shp \
                       --shp-field basin_name \
                       --geo-category freshwater_basin \
                       --output enriched_analysis/

  # Add custom grouping variable
  boldgenotyper-enrich Sphyrnidae_annotated.csv \
                       --add-metadata expedition_data.csv \
                       --add-grouping sampling_expedition \
                       --output enriched_analysis/

  # Combine multiple operations
  boldgenotyper-enrich Sphyrnidae_annotated.csv \
                       --add-metadata data1.csv \
                       --add-metadata data2.csv \
                       --recalculate-geography \
                       --add-grouping expedition \
                       --output enriched_analysis/

For more information: https://github.com/your-repo/boldgenotyper
        """
    )

    parser.add_argument(
        'input_csv',
        type=Path,
        help='Input annotated CSV from BOLDGenotyper output'
    )

    parser.add_argument(
        '--add-metadata',
        type=Path,
        action='append',
        dest='metadata_files',
        help='CSV file(s) with additional metadata to merge (can specify multiple times)'
    )

    parser.add_argument(
        '--join-on',
        type=str,
        default='processid',
        help='Column name to join metadata on (default: processid)'
    )

    parser.add_argument(
        '--custom-shp',
        type=Path,
        dest='shapefile_path',
        help='Path to custom shapefile for geographic region assignment '
             '(ocean basins, freshwater basins, ecoregions, etc.)'
    )

    parser.add_argument(
        '--shp-field',
        type=str,
        default='name',
        dest='shapefile_field',
        help='Name of shapefile attribute containing region labels (default: name). '
             'Use this to specify which field in your shapefile contains the region names.'
    )

    parser.add_argument(
        '--geo-category',
        type=str,
        default='ocean_basin',
        dest='geo_category',
        help='Name for geographic category (default: ocean_basin). '
             'Examples: "ecoregion", "watershed", "freshwater_basin", "biome"'
    )

    parser.add_argument(
        '--add-grouping',
        type=str,
        dest='grouping_column',
        help='Column name to use for custom grouping (e.g., sampling_expedition)'
    )

    parser.add_argument(
        '--recalculate-geography',
        action='store_true',
        help='Recalculate all geographic summaries with updated data'
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=Path('enriched_analysis'),
        help='Output directory for enriched data (default: enriched_analysis/)'
    )

    parser.add_argument(
        '--organism',
        type=str,
        help='Organism name (extracted from filename if not provided)'
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

    # Validate input
    if not args.input_csv.exists():
        print(f"Error: Input CSV not found: {args.input_csv}", file=sys.stderr)
        return 1

    # Check that at least one enrichment operation is specified
    if not any([args.metadata_files, args.shapefile_path,
                args.grouping_column, args.recalculate_geography]):
        print("Error: No enrichment operations specified.", file=sys.stderr)
        print("  Specify at least one of:", file=sys.stderr)
        print("    --add-metadata", file=sys.stderr)
        print("    --custom-shp", file=sys.stderr)
        print("    --add-grouping", file=sys.stderr)
        print("    --recalculate-geography", file=sys.stderr)
        return 1

    # Validate metadata files
    if args.metadata_files:
        for metadata_file in args.metadata_files:
            if not metadata_file.exists():
                print(f"Error: Metadata file not found: {metadata_file}", file=sys.stderr)
                return 1

    # Validate shapefile
    if args.shapefile_path and not args.shapefile_path.exists():
        print(f"Error: Shapefile not found: {args.shapefile_path}", file=sys.stderr)
        return 1

    # Run enrichment
    try:
        results = metadata_enrichment.enrich_metadata(
            input_csv=args.input_csv,
            output_dir=args.output,
            metadata_files=args.metadata_files,
            join_column=args.join_on,
            shapefile_path=args.shapefile_path,
            shapefile_field=args.shapefile_field,
            geo_category=args.geo_category,
            grouping_column=args.grouping_column,
            recalculate_geography=args.recalculate_geography,
            organism=args.organism
        )

        print("\n✓ Metadata enrichment complete!")
        print(f"  Enriched CSV: {results['enriched_csv']}")
        print(f"  Report: {results['report']}")
        print(f"\nEnrichment summary:")
        print(f"  Samples: {results['n_samples']}")
        print(f"  Columns: {results['n_columns']}")
        print(f"  Merge operations: {results['merge_operations']}")
        print(f"  Geographic regions updated: {results['region_updated']}")
        print(f"  Grouping added: {results['grouping_added']}")
        print(f"  Visualizations generated: {results['plots_generated']}")

        if results['plot_paths']:
            print("\nGenerated plots:")
            for plot_name, plot_path in results['plot_paths'].items():
                print(f"  - {plot_path.name}")
        print()

        return 0

    except Exception as e:
        logger.error(f"Metadata enrichment failed: {e}", exc_info=True)
        print(f"\nError: {e}", file=sys.stderr)
        return 1


def main():
    """Main CLI entry point."""
    # Support subcommands such as:
    #   boldgenotyper cluster-diagnostics ...
    if len(sys.argv) > 1 and sys.argv[1] == "cluster-diagnostics":
        # Pass everything after 'cluster-diagnostics' into the subcommand parser
        return main_cluster_diagnostics(sys.argv[2:])

    parser = argparse.ArgumentParser(
        description='BOLDGenotyper: Automated genotyping pipeline for BOLD barcode data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (organism and output inferred from filename)
  boldgenotyper data/Euprymna_scolopes.tsv

  # Specify organism name
  boldgenotyper data/samples.tsv --organism Euprymna

  # Specify output directory
  boldgenotyper data/Carcharhinus.tsv --output results/Carcharhinus_analysis

  # Adjust similarity threshold for highly diverse taxa
  boldgenotyper data/Carcharhinus.tsv --similarity-threshold 0.80

  # Use stricter clustering for fine-scale genotyping
  boldgenotyper data/Population_study.tsv --clustering-threshold 0.005

  # Adjust tie detection parameters for ambiguous assignments
  boldgenotyper data/Complex_group.tsv --tie-margin 0.005 --tie-threshold 0.97

  # Enable phylogenetic tree building
  boldgenotyper data/Euprymna.tsv --build-tree

  # Skip HTML report generation
  boldgenotyper data/Euprymna.tsv --no-report

Notes:
  - An HTML summary report is generated by default in the reports/ directory
  - Use --no-report to skip HTML report generation
  - HTML reports require jinja2 (included in default installation)

For more information: https://github.com/your-repo/boldgenotyper
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
        '--clustering-threshold',
        type=float,
        default=0.03,
        help='Maximum genetic distance for clustering sequences into consensus groups. '
             'Lower values create more groups with tighter genetic similarity. '
             'Default: 0.03 (97%% identity)'
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
        '--export-plot-data',
        action='store_true',
        help='Export raw plot data and R regeneration scripts for custom publication figures'
    )

    parser.add_argument(
        '--export-format',
        nargs='+',
        choices=['arlequin', 'popart', 'dnasp', 'generic', 'all'],
        help='Export genotypes to population genetics software formats. '
             'Options: arlequin (Arlequin), popart (PopART/NEXUS), dnasp (DnaSP), '
             'generic (CSV/FASTA), all (export all formats). '
             'Example: --export-format arlequin popart'
    )

    parser.add_argument(
        '--no-report',
        action='store_true',
        help='Skip generating HTML summary report'
    )

    parser.add_argument(
        '--no-geo',
        action='store_true',
        help='Skip geographic analysis (ocean basin assignment and related visualizations). '
             'Use this if you only need genotyping and phylogeny without geographic distribution.'
    )

    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging verbosity (default: INFO)'
    )

    parser.add_argument(
        '--version',
        action='version',
        version='BOLDGenotyper 1.0.0'
    )

    args = parser.parse_args()

    # Validate input file
    if not args.tsv.exists():
        print(f"Error: Input TSV file not found: {args.tsv}", file=sys.stderr)
        return 1

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
        dereplication__clustering_threshold=args.clustering_threshold,
        n_threads=args.threads,
        output_dir=output_dir,
        log_level=args.log_level,
        phylogenetic__build_tree=args.build_tree,
        keep_intermediates=True
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
    print(f"  Clustering threshold: {args.clustering_threshold} ({(1-args.clustering_threshold)*100:.1f}% identity)")
    print(f"  Similarity threshold: {args.similarity_threshold} ({args.similarity_threshold*100:.0f}% identity)")
    print(f"  Tie margin: {args.tie_margin} ({args.tie_margin*100:.1f}% difference)")
    print(f"  Tie threshold: {args.tie_threshold} ({args.tie_threshold*100:.0f}% identity)")
    print(f"  Threads: {args.threads}")
    print(f"  Build tree: {args.build_tree}")
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
            export_plot_data=args.export_plot_data,
            export_popgen_formats=args.export_format
        )

        return 0 if success else 1

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}", exc_info=True)
        print(f"\nError: Pipeline failed. Check log file: {log_file}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
