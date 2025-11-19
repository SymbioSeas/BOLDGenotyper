"""
Core Pipeline Orchestration for BOLDGenotyper

This module provides the main pipeline orchestration for automated COI sequence
genotyping and biogeographic analysis from BOLD data. It coordinates the entire
workflow from TSV input through visualization output.

The core pipeline integrates all analysis steps:
1. Data input and validation from BOLD TSV files
2. Sequence dereplication and consensus generation
3. Genotype assignment to samples
4. Geographic coordinate filtering and ocean basin assignment
5. Visualization of distribution patterns and relative abundances
6. Optional phylogenetic tree construction

Key Features:
- Configurable clustering thresholds (default: 99% sequence identity)
- Robust coordinate filtering to exclude country centroids
- Publication-ready figure generation (PNG and PDF formats)
- Comprehensive logging and progress tracking
- Graceful error handling with informative messages

The pipeline is designed to be flexible for both command-line usage and
programmatic integration into custom workflows.

Example Usage:
    >>> from boldgenotyper.core import run_pipeline
    >>> results = run_pipeline(
    ...     input_tsv="Sphyrna_lewini.tsv",
    ...     output_dir="results/",
    ...     cluster_threshold=0.01,
    ...     build_phylogeny=True
    ... )

Author: Steph Smith (steph.smith@unc.edu)
"""

from typing import Dict, Optional, Any
from pathlib import Path
import logging
import pandas as pd

# Local imports
from . import (
    utils, config, metadata, geographic, dereplication,
    genotype_assignment, phylogenetics, visualization, reports
)

# Configure logging
logger = logging.getLogger(__name__)


def run_pipeline(
    input_tsv: str,
    output_dir: str,
    cluster_threshold: float = 0.01,
    min_identity: float = 0.90,
    build_phylogeny: bool = False,
    outgroup_fasta: Optional[str] = None,
    threads: int = 1,
    organism_name: Optional[str] = None,
    skip_geographic: bool = False,
    config_obj: Optional[config.PipelineConfig] = None,
) -> Dict[str, Any]:
    """
    Run the complete BOLDGenotyper analysis pipeline.

    This function orchestrates the entire workflow from BOLD TSV input through
    visualization and report generation. It is designed to be called programmatically
    and returns structured results for further analysis.

    Pipeline Phases:
    1. Data loading and quality control
    2. Sequence dereplication and consensus generation
    3. Genotype assignment to samples
    4. Geographic analysis (ocean basin assignment)
    5. Phylogenetic analysis (optional)
    6. Visualization generation
    7. Report generation

    Parameters
    ----------
    input_tsv : str
        Path to BOLD TSV file containing sequence and metadata
    output_dir : str
        Directory for output files (will be created if it doesn't exist)
    cluster_threshold : float, optional
        Sequence distance threshold for clustering (default: 0.01 = 99% identity)
    min_identity : float, optional
        Minimum identity for genotype assignment (default: 0.90)
    build_phylogeny : bool, optional
        Whether to build phylogenetic tree (default: False)
    outgroup_fasta : str, optional
        Path to outgroup sequences for phylogeny rooting
    threads : int, optional
        Number of CPU threads to use (default: 1)
    organism_name : str, optional
        Organism name for file naming (default: extracted from input filename)
    skip_geographic : bool, optional
        Skip geographic analysis if GOaS data unavailable (default: False)
    config_obj : PipelineConfig, optional
        Custom configuration object (default: uses default configuration)

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'success': bool - Whether pipeline completed successfully
        - 'output_dir': Path - Output directory path
        - 'organism': str - Organism name
        - 'n_samples_input': int - Number of input samples
        - 'n_samples_filtered': int - Number after filtering
        - 'n_genotypes': int - Number of consensus genotypes
        - 'n_assigned': int - Number of successfully assigned samples
        - 'assignment_rate': float - Fraction of samples assigned
        - 'files': Dict[str, Path] - Paths to output files
        - 'geographic_analysis': bool - Whether geographic analysis was performed
        - 'phylogenetic_tree': Optional[Path] - Path to tree file if built
        - 'errors': List[str] - Any non-critical errors encountered

    Raises
    ------
    FileNotFoundError
        If input_tsv does not exist
    ValueError
        If parameters are invalid (e.g., threshold out of range)

    Examples
    --------
    >>> from boldgenotyper.core import run_pipeline
    >>>
    >>> # Basic usage
    >>> results = run_pipeline(
    ...     input_tsv="data/Euprymna_scolopes.tsv",
    ...     output_dir="results/euprymna",
    ...     threads=4
    ... )
    >>> print(f"Identified {results['n_genotypes']} genotypes")
    >>>
    >>> # Custom parameters with phylogeny
    >>> results = run_pipeline(
    ...     input_tsv="data/Sphyrna_lewini.tsv",
    ...     output_dir="results/sphyrna",
    ...     cluster_threshold=0.005,  # More stringent clustering
    ...     min_identity=0.95,
    ...     build_phylogeny=True,
    ...     threads=8
    ... )

    Notes
    -----
    - External tools required: MAFFT, trimAl (always), FastTree (if build_phylogeny=True)
    - Geographic analysis requires GOaS shapefile (can be downloaded via goas_downloader)
    - Progress and detailed logs are written to {output_dir}/{organism}_pipeline.log
    - Intermediate files are saved in {output_dir}/intermediate/ for debugging

    See Also
    --------
    boldgenotyper.cli.main : Command-line interface wrapper
    boldgenotyper.config.PipelineConfig : Configuration details
    """
    # Initialize results dictionary
    results = {
        'success': False,
        'output_dir': Path(output_dir),
        'errors': [],
        'files': {},
    }

    try:
        # =====================================================================
        # Setup and Validation
        # =====================================================================

        # Validate input file
        input_path = Path(input_tsv)
        if not input_path.exists():
            raise FileNotFoundError(f"Input TSV not found: {input_tsv}")

        # Validate parameters
        if not 0 < cluster_threshold < 1:
            raise ValueError(f"cluster_threshold must be between 0 and 1, got {cluster_threshold}")
        if not 0 < min_identity <= 1:
            raise ValueError(f"min_identity must be between 0 and 1, got {min_identity}")
        if threads < 1:
            raise ValueError(f"threads must be >= 1, got {threads}")

        # Determine organism name
        if organism_name is None:
            organism_name = _extract_organism_name(input_path)
        results['organism'] = organism_name

        # Create output directory structure
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        dirs = _setup_directories(output_path)

        # Setup logging
        log_file = output_path / f"{organism_name}_pipeline.log"
        utils.setup_logging(log_level='INFO', log_file=str(log_file))

        logger.info("=" * 80)
        logger.info(f"BOLDGenotyper Pipeline - {organism_name}")
        logger.info("=" * 80)
        logger.info(f"Input: {input_tsv}")
        logger.info(f"Output: {output_dir}")
        logger.info(f"Clustering threshold: {cluster_threshold}")
        logger.info(f"Assignment threshold: {min_identity}")
        logger.info(f"Threads: {threads}")
        logger.info("")

        # Get or create configuration
        if config_obj is None:
            cfg = config.get_default_config()
            cfg = cfg.update(
                dereplication__clustering_threshold=cluster_threshold,
                genotype_assignment__min_identity=min_identity,
                phylogenetic__build_tree=build_phylogeny,
                n_threads=threads,
                output_dir=output_path
            )
        else:
            cfg = config_obj

        # =====================================================================
        # Phase 1: Data Loading and Quality Control
        # =====================================================================

        logger.info("PHASE 1: Data Loading and Quality Control")
        logger.info("-" * 80)

        # Parse BOLD TSV
        logger.info("1.1: Parsing BOLD TSV metadata...")
        df = metadata.parse_bold_tsv(input_tsv)
        results['n_samples_input'] = len(df)
        logger.info(f"  ✓ Loaded {len(df)} samples")

        # Save parsed metadata
        parsed_tsv = dirs['intermediate'] / "01_parsed_metadata.tsv"
        df.to_csv(parsed_tsv, sep='\t', index=False)
        results['files']['parsed_metadata'] = parsed_tsv

        # Filter by coordinate quality
        logger.info("1.2: Filtering coordinate quality...")
        df_filtered = metadata.filter_by_coordinate_quality(df, cfg.geographic)
        results['n_samples_filtered'] = len(df_filtered)
        logger.info(f"  ✓ Retained {len(df_filtered)}/{len(df)} samples after filtering")

        filtered_tsv = dirs['intermediate'] / "02_filtered_metadata.tsv"
        df_filtered.to_csv(filtered_tsv, sep='\t', index=False)
        results['files']['filtered_metadata'] = filtered_tsv

        # Assign ocean basins (if not skipped)
        logger.info("1.3: Assigning ocean basins...")
        geographic_analysis_performed = False

        if skip_geographic:
            logger.info("  ⊘ Geographic analysis disabled")
            df_with_basins = df_filtered.copy()
            df_with_basins['ocean_basin'] = 'Unknown'
        elif not cfg.geographic.goas_shapefile_path.exists():
            logger.warning(f"  ⚠ GOaS shapefile not found, skipping geographic analysis")
            results['errors'].append("GOaS shapefile not found")
            df_with_basins = df_filtered.copy()
            df_with_basins['ocean_basin'] = 'Unknown'
        else:
            try:
                goas_data = geographic.load_goas_data(cfg.geographic.goas_shapefile_path)
                df_with_basins = geographic.assign_ocean_basins(
                    df_filtered, goas_data=goas_data, coord_col="coord"
                )
                logger.info(f"  ✓ Assigned ocean basins")
                geographic_analysis_performed = True
            except Exception as e:
                logger.warning(f"  ⚠ Geographic analysis failed: {e}")
                results['errors'].append(f"Geographic analysis failed: {e}")
                df_with_basins = df_filtered.copy()
                df_with_basins['ocean_basin'] = 'Unknown'

        results['geographic_analysis'] = geographic_analysis_performed
        basins_tsv = dirs['geographic'] / "samples_with_ocean_basins.tsv"
        df_with_basins.to_csv(basins_tsv, sep='\t', index=False)
        results['files']['samples_with_basins'] = basins_tsv

        # =====================================================================
        # Phase 2: Sequence Dereplication and Consensus Generation
        # =====================================================================

        logger.info("")
        logger.info("PHASE 2: Sequence Dereplication and Consensus Generation")
        logger.info("-" * 80)

        # Generate FASTA from sequences
        logger.info("2.1: Generating FASTA from sequences...")
        logger.info(f"  Minimum sequence length: {cfg.dereplication.min_sequence_length} bp")
        fasta_records = []
        skipped = 0

        for _, row in df_filtered.iterrows():
            header = f"{organism_name}_{row['processid']}.COI-5P"
            sequence = row['nuc']

            if not isinstance(sequence, str) or not sequence.strip():
                skipped += 1
                continue

            is_valid, reason = utils.validate_sequence(
                sequence,
                min_length=cfg.dereplication.min_sequence_length
            )
            if is_valid:
                fasta_records.append((header, sequence))
            else:
                skipped += 1

        if skipped > 0:
            logger.warning(f"  ⚠ Skipped {skipped} invalid sequences")
            results['errors'].append(f"Skipped {skipped} invalid sequences")

        fasta_path = dirs['intermediate'] / f"{organism_name}.fasta"
        utils.write_fasta(fasta_records, fasta_path)
        results['files']['sequences_fasta'] = fasta_path
        logger.info(f"  ✓ Created {len(fasta_records)} FASTA records")

        # Dereplicate sequences
        logger.info("2.2: Dereplicating sequences...")
        consensus_records = dereplication.dereplicate_from_fasta(
            input_fasta=str(fasta_path),
            output_dir=str(dirs['dereplication']),
            threshold=cfg.dereplication.clustering_threshold,
            frequency_cutoff=cfg.dereplication.consensus_frequency_cutoff,
            min_post_trim_length=cfg.dereplication.min_post_trim_length,
            min_consensus_length_ratio=cfg.dereplication.min_consensus_length_ratio
        )
        results['n_genotypes'] = len(consensus_records)
        logger.info(f"  ✓ Identified {len(consensus_records)} consensus genotypes")

        # Save consensus sequences
        consensus_path = dirs['consensus'] / f"{organism_name}_consensus.fasta"
        utils.write_fasta(
            [(rec.id, str(rec.seq)) for rec in consensus_records.values()],
            consensus_path
        )
        results['files']['consensus_sequences'] = consensus_path

        # =====================================================================
        # Phase 3: Genotype Assignment
        # =====================================================================

        logger.info("")
        logger.info("PHASE 3: Genotype Assignment")
        logger.info("-" * 80)

        logger.info("3.1: Assigning samples to genotypes...")
        annotated_tsv = dirs['assignments'] / f"{organism_name}_with_genotypes.tsv"
        diagnostics_csv = dirs['assignments'] / f"{organism_name}_diagnostics.csv"

        stats = genotype_assignment.assign_genotypes(
            metadata_path=str(filtered_tsv),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(annotated_tsv),
            min_identity=cfg.genotype_assignment.min_identity,
            n_processes=cfg.n_threads,
            diagnostics_path=str(diagnostics_csv),
            tie_margin=cfg.genotype_assignment.tie_margin,
            tie_min_identity=cfg.genotype_assignment.tie_threshold
        )

        results['n_assigned'] = stats['assigned']
        results['assignment_rate'] = stats['assigned'] / stats['total'] if stats['total'] > 0 else 0.0
        logger.info(f"  ✓ Assigned {stats['assigned']}/{stats['total']} samples")
        logger.info(f"  ✓ Assignment rate: {results['assignment_rate']*100:.1f}%")

        results['files']['genotype_assignments'] = annotated_tsv
        results['files']['assignment_diagnostics'] = diagnostics_csv

        # Load annotated data
        df_with_genotypes = pd.read_csv(annotated_tsv, sep="\t")

        # =====================================================================
        # Phase 4: Taxonomy Assignment
        # =====================================================================

        logger.info("")
        logger.info("PHASE 4: Taxonomy Assignment")
        logger.info("-" * 80)

        logger.info("4.1: Assigning taxonomy to consensus groups...")
        assign_table, species_counts = utils.assign_consensus_taxonomy(
            df_with_genotypes,
            group_col="consensus_group",
            species_col="species",
            genus_col="genus",
            majority_threshold=cfg.taxonomy.majority_species_threshold
        )

        # Build consensus_group_sp labels
        assign_table["consensus_group_short"] = assign_table["consensus_group"].str.replace(
            "consensus_", "", regex=False
        )
        assign_table["consensus_group_sp"] = (
            assign_table["assigned_sp"].fillna("") + " " +
            assign_table["consensus_group_short"]
        ).str.strip()

        # Save taxonomy files
        species_counts_out = dirs['taxonomy'] / f"{organism_name}_species_by_consensus.csv"
        assign_table_out = dirs['taxonomy'] / f"{organism_name}_consensus_taxonomy.csv"
        species_counts.to_csv(species_counts_out, index=False)
        assign_table.to_csv(assign_table_out, index=False)

        results['files']['species_counts'] = species_counts_out
        results['files']['consensus_taxonomy'] = assign_table_out

        # Merge taxonomy back into dataframe
        df_with_genotypes = df_with_genotypes.merge(
            assign_table[["consensus_group", "assigned_sp", "consensus_group_sp",
                         "assignment_level", "assignment_notes", "majority_fraction"]],
            on="consensus_group",
            how="left"
        )

        logger.info(f"  ✓ Assigned taxonomy to {len(assign_table)} consensus groups")

        # Merge with geographic data
        geo_cols = [c for c in ['processid', 'lat', 'lon', 'ocean_basin']
                   if c in df_with_basins.columns]
        df_final = df_with_genotypes.merge(
            df_with_basins[geo_cols],
            on='processid',
            how='left'
        )

        # Save final annotated file
        annotated_csv = dirs['base'] / f"{organism_name}_annotated.csv"
        df_final.to_csv(annotated_csv, index=False)
        results['files']['annotated_data'] = annotated_csv
        logger.info(f"  ✓ Saved annotated dataset: {annotated_csv}")

        # =====================================================================
        # Phase 5: Phylogenetic Analysis (Optional)
        # =====================================================================

        logger.info("")
        logger.info("PHASE 5: Phylogenetic Analysis")
        logger.info("-" * 80)

        tree_path = None
        if cfg.phylogenetic.build_tree:
            # Check for required tools
            if not utils.check_external_tool("mafft"):
                logger.warning("  ⚠ MAFFT not found, skipping phylogeny")
                results['errors'].append("MAFFT not found for phylogeny")
            elif not utils.check_external_tool("fasttree"):
                logger.warning("  ⚠ FastTree not found, skipping phylogeny")
                results['errors'].append("FastTree not found for phylogeny")
            else:
                try:
                    logger.info("5.1: Building phylogenetic tree...")
                    output_prefix = dirs['phylogenetic'] / organism_name
                    tree = phylogenetics.build_phylogeny(
                        consensus_fasta=str(consensus_path),
                        output_prefix=str(output_prefix),
                        threads=cfg.n_threads
                    )

                    tree_path = Path(f"{output_prefix}_tree.nwk")
                    if tree_path.exists():
                        logger.info(f"  ✓ Built phylogenetic tree")
                        results['files']['phylogenetic_tree'] = tree_path
                        results['phylogenetic_tree'] = tree_path

                        # Relabel tree
                        logger.info("5.2: Relabeling tree with taxonomy...")
                        alignment_path = Path(f"{output_prefix}_aligned.fasta")
                        if alignment_path.exists():
                            relabeled_tree = Path(f"{output_prefix}_tree_relabeled.nwk")
                            relabeled_aln = Path(f"{output_prefix}_aligned_relabeled.fasta")

                            phylogenetics.relabel_tree_and_alignment(
                                tree_file=str(tree_path),
                                alignment_file=str(alignment_path),
                                taxonomy_csv=str(assign_table_out),
                                output_tree=str(relabeled_tree),
                                output_alignment=str(relabeled_aln)
                            )
                            logger.info(f"  ✓ Created relabeled tree")
                            results['files']['relabeled_tree'] = relabeled_tree
                    else:
                        logger.warning("  ⚠ Tree file not created")
                        results['errors'].append("Phylogenetic tree file not created")

                except Exception as e:
                    logger.warning(f"  ⚠ Phylogenetic analysis failed: {e}")
                    results['errors'].append(f"Phylogenetic analysis failed: {e}")
        else:
            logger.info("  ⊘ Phylogenetic tree building disabled")

        results['phylogenetic_tree'] = tree_path

        # =====================================================================
        # Phase 6: Visualization
        # =====================================================================

        logger.info("")
        logger.info("PHASE 6: Visualization")
        logger.info("-" * 80)

        try:
            logger.info("6.1: Generating visualizations...")

            # Create consistent genotype color map based on abundance
            # This ensures all visualizations use the same colors for each genotype
            from boldgenotyper.visualization import _get_genotype_order_by_abundance, get_genotype_colors

            if 'consensus_group' in df_final.columns and 'consensus_group_sp' in df_final.columns:
                # Order genotypes by abundance (most abundant gets first color)
                genotype_order = _get_genotype_order_by_abundance(df_final, 'consensus_group')
                colors = get_genotype_colors(len(genotype_order))

                # Create color map DataFrame
                color_map_data = []
                for i, geno in enumerate(genotype_order):
                    # Get corresponding consensus_group_sp label
                    sp_label = df_final[df_final['consensus_group'] == geno]['consensus_group_sp'].iloc[0]
                    color_map_data.append({
                        'consensus_group': geno,
                        'consensus_group_sp': sp_label,
                        'color': colors[i],
                        'rank': i + 1  # 1-indexed rank by abundance
                    })

                # Save color map to multiple directories for easy access
                color_map_df = pd.DataFrame(color_map_data)
                for dir_key in ['visualization', 'phylogenetic', 'geographic']:
                    if dir_key in dirs:
                        color_map_path = dirs[dir_key] / f"{organism_name}_genotype_color_map.csv"
                        color_map_df.to_csv(color_map_path, index=False)

                logger.debug(f"Created genotype color map with {len(genotype_order)} genotypes ordered by abundance")

            for fmt in cfg.visualization.figure_format:
                # Geographic visualizations (if analysis was performed)
                if geographic_analysis_performed:
                    # Distribution map
                    if 'lat' in df_final.columns and 'lon' in df_final.columns:
                        try:
                            map_path = dirs['geographic'] / f"{organism_name}_distribution_map.{fmt}"
                            visualization.plot_distribution_map(
                                df=df_final,
                                output_path=str(map_path),
                                genotype_column='consensus_group_sp',
                                latitude_col='lat',
                                longitude_col='lon'
                            )
                            results['files'][f'distribution_map_{fmt}'] = map_path
                        except Exception as e:
                            logger.debug(f"Distribution map skipped: {e}")

                    # Ocean basin abundance
                    if 'ocean_basin' in df_final.columns:
                        try:
                            bar_path = dirs['geographic'] / f"{organism_name}_distribution_bar.{fmt}"
                            visualization.plot_ocean_basin_abundance(
                                df=df_final,
                                output_path=str(bar_path),
                                genotype_column='consensus_group_sp',
                                basin_column='ocean_basin'
                            )
                            results['files'][f'basin_abundance_{fmt}'] = bar_path
                        except Exception as e:
                            logger.debug(f"Basin abundance skipped: {e}")

                # Identity distribution (always generated if diagnostics exist)
                if diagnostics_csv.exists():
                    try:
                        identity_path = dirs['assignments'] / f"{organism_name}_identity_distribution.{fmt}"
                        visualization.plot_identity_distribution(
                            diagnostics_csv=str(diagnostics_csv),
                            output_path=str(identity_path),
                            dpi=cfg.visualization.figure_dpi
                        )
                        results['files'][f'identity_distribution_{fmt}'] = identity_path
                    except Exception as e:
                        logger.debug(f"Identity distribution skipped: {e}")

                # Phylogenetic tree visualization
                if tree_path and tree_path.exists():
                    try:
                        # Use relabeled tree if it exists (so tips show consensus_group_sp labels)
                        relabeled_tree_path = Path(str(tree_path).replace("_tree.nwk", "_tree_relabeled.nwk"))
                        if relabeled_tree_path.exists():
                            tree_to_plot = relabeled_tree_path
                            logger.debug(f"Using relabeled tree for visualization: {relabeled_tree_path}")
                        else:
                            tree_to_plot = tree_path
                            logger.debug(f"Using original tree for visualization: {tree_path}")

                        # Load genotype color map if it exists
                        genotype_colors = None
                        color_map_path = dirs['phylogenetic'] / f"{organism_name}_genotype_color_map.csv"
                        if color_map_path.exists():
                            try:
                                import pandas as pd
                                color_df = pd.read_csv(color_map_path)
                                if 'consensus_group_sp' in color_df.columns and 'color' in color_df.columns:
                                    genotype_colors = dict(zip(color_df['consensus_group_sp'], color_df['color']))
                                    logger.debug(f"Loaded {len(genotype_colors)} genotype colors")
                            except Exception as e:
                                logger.debug(f"Could not load color map: {e}")

                        tree_plot_path = dirs['phylogenetic'] / f"{organism_name}_tree.{fmt}"
                        visualization.plot_phylogenetic_tree(
                            tree_file=str(tree_to_plot),
                            output_path=str(tree_plot_path),
                            genotype_colors=genotype_colors,
                            show_bootstrap=True,
                            dpi=cfg.visualization.figure_dpi
                        )
                        results['files'][f'tree_plot_{fmt}'] = tree_plot_path
                    except Exception as e:
                        logger.debug(f"Tree visualization skipped: {e}")

            logger.info(f"  ✓ Generated visualizations")

        except Exception as e:
            logger.warning(f"  ⚠ Visualization generation had errors: {e}")
            results['errors'].append(f"Visualization errors: {e}")

        # =====================================================================
        # Phase 7: Reports
        # =====================================================================

        logger.info("")
        logger.info("PHASE 7: Generating Reports")
        logger.info("-" * 80)

        try:
            summary_output = dirs['reports'] / f"{organism_name}_assignment_summary.csv"
            reports.generate_assignment_summary(
                annotated_csv=str(annotated_csv),
                diagnostics_csv=str(diagnostics_csv),
                output_csv=str(summary_output)
            )

            if summary_output.exists():
                logger.info(f"  ✓ Generated assignment summary")
                results['files']['assignment_summary'] = summary_output
            else:
                logger.warning(f"  ⚠ Failed to generate summary report")
                results['errors'].append("Summary report not created")

        except Exception as e:
            logger.warning(f"  ⚠ Report generation failed: {e}")
            results['errors'].append(f"Report generation failed: {e}")

        # =====================================================================
        # Pipeline Complete
        # =====================================================================

        results['success'] = True

        logger.info("")
        logger.info("=" * 80)
        logger.info(f"✓ Pipeline completed successfully for {organism_name}")
        logger.info(f"  Samples: {results['n_samples_filtered']} filtered / {results['n_samples_input']} total")
        logger.info(f"  Genotypes: {results['n_genotypes']}")
        logger.info(f"  Assignment rate: {results['assignment_rate']*100:.1f}%")
        logger.info(f"  Output: {output_dir}")
        logger.info("=" * 80)

        return results

    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}", exc_info=True)
        results['success'] = False
        results['errors'].append(str(e))
        raise


def _extract_organism_name(path: Path) -> str:
    """
    Extract organism name from TSV filename.

    Parameters
    ----------
    path : Path
        Path to TSV file

    Returns
    -------
    str
        Extracted organism name

    Examples
    --------
    >>> _extract_organism_name(Path("Genus_species.tsv"))
    'Genus_species'
    >>> _extract_organism_name(Path("data/Genus.tsv"))
    'Genus'
    """
    name = path.stem
    parts = name.split('_')

    # Check for common patterns
    if len(parts) >= 2 and parts[1] not in ['data', 'bold', 'samples', 'sequences']:
        return '_'.join(parts[:2])

    return parts[0].capitalize()


def _setup_directories(base_output: Path) -> Dict[str, Path]:
    """
    Create organized output directory structure.

    Parameters
    ----------
    base_output : Path
        Base output directory

    Returns
    -------
    Dict[str, Path]
        Dictionary mapping directory names to paths
    """
    dirs = {
        'base': base_output,
        'intermediate': base_output / 'intermediate',
        'dereplication': base_output / 'intermediate' / 'dereplication',
        'consensus': base_output / 'consensus_sequences',
        'assignments': base_output / 'genotype_assignments',
        'taxonomy': base_output / 'taxonomy',
        'geographic': base_output / 'geographic',
        'phylogenetic': base_output / 'phylogenetic',
        'reports': base_output / 'reports',
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs
