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
from pathlib import Path
from typing import Optional
import pandas as pd

# Local imports
from . import (
    utils, config, metadata, geographic, dereplication,
    genotype_assignment, phylogenetics, visualization, reports
)

logger = logging.getLogger(__name__)


def extract_organism_from_path(path: Path) -> str:
    """Extract organism name from TSV filename."""
    # Remove extension and path
    name = path.stem
    # Common patterns: "Genus_species.tsv", "Genus.tsv", "genus_species_data.tsv"
    # Take first part before underscore or use whole name if no underscore
    parts = name.split('_')
    if len(parts) >= 2 and parts[1] not in ['data', 'bold', 'samples', 'sequences']:
        # Likely "Genus_species" format
        return '_'.join(parts[:2])
    return parts[0].capitalize()


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
        'phylogenetic': base_output / 'phylogenetic',
        'visualization': base_output / 'visualization',
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
) -> bool:
    """
    Run the complete BOLDGenotyper pipeline.

    Returns True if successful, False otherwise.
    """
    logger.info("=" * 80)
    logger.info(f"BOLDGenotyper Pipeline - {organism}")
    logger.info("=" * 80)
    logger.info(f"Input TSV: {tsv_path}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("")

    # Setup directory structure
    dirs = setup_directories(output_dir)

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

        # Filter coordinates
        logger.info("1.2: Filtering coordinate quality...")
        df_filtered = metadata.filter_by_coordinate_quality(df, cfg.geographic)
        logger.info(f"  ✓ Retained {len(df_filtered)}/{len(df)} samples after coordinate filtering")

        filtered_tsv = dirs['intermediate'] / "02_filtered_metadata.tsv"
        df_filtered.to_csv(filtered_tsv, sep='\t', index=False)

        # Assign ocean basins
        logger.info("1.3: Assigning ocean basins...")
        if cfg.geographic.goas_shapefile_path.exists():
            goas_data = geographic.load_goas_data(cfg.geographic.goas_shapefile_path)
            df_with_basins = geographic.assign_ocean_basins(
                df_filtered, goas_data=goas_data, coord_col="coord"
            )
            logger.info(f"  ✓ Assigned ocean basins to samples")
        else:
            logger.warning(f"  ⚠ GOaS shapefile not found, skipping basin assignment")
            df_with_basins = df_filtered.copy()
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
        # Generate FASTA
        logger.info("2.1: Generating FASTA from sequences...")
        logger.info(f"  Minimum sequence length: {cfg.dereplication.min_sequence_length} bp")
        fasta_records = []
        skipped_count = 0
        skipped_reasons = {}
        for _, row in df_filtered.iterrows():
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
        annotated_tsv = dirs['intermediate_assignments'] / f"{organism}_with_genotypes.tsv"
        diagnostics_csv = dirs['assignments'] / f"{organism}_diagnostics.csv"

        stats = genotype_assignment.assign_genotypes(
            metadata_path=str(filtered_tsv),
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

        # Merge with geographic data
        geo_keep = [c for c in ['processid', 'lat', 'lon', 'ocean_basin'] if c in df_with_basins.columns]
        df_final = df_with_genotypes.merge(df_with_basins[geo_keep], on='processid', how='left', validate='one_to_one')

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
                    threads=cfg.n_threads
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
                                output_alignment=relabeled_alignment_path
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
    # PHASE 6: Visualization
    # ========================================================================
    logger.info("")
    logger.info("PHASE 6: Visualization")
    logger.info("-" * 80)

    try:
        logger.info("6.1: Generating plots...")

        # Generate visualizations for each format
        for fmt in cfg.visualization.figure_format:
            # Distribution maps
            if 'lat' in df_final.columns and 'lon' in df_final.columns:
                try:
                    visualization.plot_distribution_map(
                        df=df_final,
                        output_path=str(dirs['visualization'] / f"{organism}_distribution_map.{fmt}"),
                        genotype_column='consensus_group_sp',
                        latitude_col='lat',
                        longitude_col='lon'
                    )
                except Exception as e:
                    logger.debug(f"Distribution map skipped: {e}")

            # Ocean basin abundance bar plot
            if 'ocean_basin' in df_final.columns and 'consensus_group_sp' in df_final.columns:
                try:
                    visualization.plot_ocean_basin_abundance(
                        df=df_final,
                        output_path=str(dirs['visualization'] / f"{organism}_distribution_bar.{fmt}"),
                        genotype_column='consensus_group_sp',
                        basin_column='ocean_basin'
                    )
                except Exception as e:
                    logger.debug(f"Ocean basin bar plot skipped: {e}")

            # Identity distribution
            if diagnostics_csv.exists():
                visualization.plot_identity_distribution(
                    diagnostics_csv=str(diagnostics_csv),
                    output_path=str(dirs['visualization'] / f"{organism}_identity_distribution.{fmt}"),
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
                color_map_path = dirs['visualization'] / f"{organism}_genotype_color_map.csv"
                genotype_colors = None
                if color_map_path.exists():
                    color_df = pd.read_csv(color_map_path)
                    genotype_colors = dict(zip(color_df['consensus_group_sp'], color_df['color']))

                visualization.plot_phylogenetic_tree(
                    tree_file=str(tree_to_plot),
                    output_path=str(dirs['visualization'] / f"{organism}_tree.{fmt}"),
                    genotype_colors=genotype_colors,
                    show_bootstrap=True,
                    bootstrap_threshold=cfg.visualization.show_bootstrap_threshold,
                    figsize=None,  # Auto-scale based on tree size
                    dpi=cfg.visualization.figure_dpi
                )

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
                        longitude_col='lon'
                    )
                except Exception as e:
                    logger.warning(f"Faceted distribution map generation failed: {e}", exc_info=True)

            # Faceted ocean basin bar plot by consensus_group_sp
            if ('ocean_basin' in df_final.columns and 'consensus_group_sp' in df_final.columns and
                'consensus_group' in df_final.columns):
                try:
                    visualization.plot_ocean_basin_abundance_faceted(
                        df=df_final,
                        output_path=str(dirs['visualization'] / f"{organism}_distribution_bar_faceted.{fmt}"),
                        genotype_column='consensus_group',
                        species_column='consensus_group_sp',
                        basin_column='ocean_basin'
                    )
                except Exception as e:
                    logger.debug(f"Faceted basin bar plot skipped: {e}")

        logger.info(f"  ✓ Generated visualization plots")

    except Exception as e:
        logger.warning(f"Visualization generation encountered errors (non-critical): {e}")

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
    logger.info("=" * 80)

    return True


def main():
    """Main CLI entry point."""
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

  # Enable phylogenetic tree building
  boldgenotyper data/Euprymna.tsv --build-tree

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
        default=0.90,
        help='Minimum similarity for genotype assignment (default: 0.90)'
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

    # Determine organism name
    organism = args.organism if args.organism else extract_organism_from_path(args.tsv)

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
    print(f"Similarity threshold: {args.similarity_threshold}")
    print(f"Threads: {args.threads}")
    print(f"Build tree: {args.build_tree}")
    print("=" * 80)
    print()

    # Run pipeline
    try:
        success = run_pipeline(
            tsv_path=args.tsv,
            organism=organism,
            output_dir=output_dir,
            cfg=cfg
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
