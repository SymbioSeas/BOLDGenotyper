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
6. Record labels (majority labels, label conflicts, depositor-note keyword hits)
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
import importlib.util
import sys
import os
import logging
import shlex
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO

# Local imports
from . import (
    utils, config, metadata, geographic, dereplication,
    phylogenetics, reports,
    quality_control, plot_export,
    divergence_analysis, parameter_sweep, geographic_enhancement,
    popgen_export, species_analysis, msa_visualization,
    metadata_summary, shapefile_registry, provenance, record_status, taxon_labels, layout,
    run_readme, haplotype_ids, taxon_evidence, query_outputs, query_overlays, query_report,
    run_figures, label_variants, error_filter_diagnostics,
)
from . import report as run_report
from . import taxon_curation as taxon_curation_mod
from . import __version__

logger = logging.getLogger(__name__)

# Flags removed in 2.0.0, per command, with what to do instead. A removed flag
# fails with one line instead of being accepted and ignored.
_MIGRATION = "see docs/guides/MIGRATION_v1_to_v2.md"
REMOVED_FLAGS = {
    "boldgenotyper": {
        "--similarity-threshold": "records are mapped directly from haplotype discovery",
        "--tie-margin": "records are mapped directly from haplotype discovery",
        "--tie-threshold": "records are mapped directly from haplotype discovery",
        "--no-export-plot-data": "the plot kit is written only with --plot-kit",
        "--no-export-format": "exports are written only with --export-format",
        "--no-metadata-analysis": "use --no-metadata-summary",
    },
    "boldgenotyper-query": {
        "--tie-margin": "use --max-mismatch-delta (integer differences)",
    },
}


def removed_flag_error(argv, command):
    """One-line error for the first flag in argv that was removed in 2.0.0, else None."""
    removed = REMOVED_FLAGS.get(command, {})
    for arg in argv:
        flag = str(arg).split("=", 1)[0]
        if flag in removed:
            return f"{flag} was removed in 2.0.0; {removed[flag]} ({_MIGRATION})"
    return None


# Shared help footer listing every installed BOLDGenotyper command, so running
# --help on any one command surfaces the others. Referenced by each parser's
# epilog; keep this the single source of truth for the command list.
AVAILABLE_COMMANDS_EPILOG = """
Available BOLDGenotyper commands:
  boldgenotyper         Run the full haplotype-first genotyping pipeline
  boldgenotyper-sweep   Optimise the singleton-distance threshold for a dataset
  boldgenotyper-query   Assign new COI sequences to haplotypes from a prior run
  Run '<command> --help' for the options and examples of any command.
"""


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
    """
    Create the working folders. Step outputs are written under
    intermediate/staging/ ('base') and placed at their v2 paths at the end of
    the run by layout.plan_moves/apply_moves (P2); 'output' is the run folder.
    """
    stage = base_output / 'intermediate' / 'staging'
    dirs = {
        'output': base_output,
        'base': stage,
        'intermediate': base_output / 'intermediate',
        'qc': base_output / 'intermediate' / 'quality_control',
        'haplotype_discovery': base_output / 'intermediate' / 'haplotype_discovery',
        'intermediate_phylo': base_output / 'intermediate' / 'phylogenetic',
        'intermediate_assignments': base_output / 'intermediate' / 'haplotype_assignments',
        'intermediate_geographic': base_output / 'intermediate' / 'geographic',
        'haplotypes': stage / 'haplotypes',
        'assignments': stage / 'haplotype_assignments',
        'taxonomy': stage / 'taxonomy',
        'species_analysis': stage / 'species_analysis',
        'quality_control': stage / 'quality_control',
        'divergence_analysis': stage / 'divergence_analysis',
        'geographic_analysis': stage / 'geographic_analysis',
        'phylogenetic': stage / 'phylogenetic',
        'reports': stage / 'reports',
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def validate_shapefile_field(shapefile_path: Path, shapefile_field: Optional[str]) -> None:
    """
    Validate that shapefile_field exists in the given shapefile.

    Prints available fields and exits with error if the field is not found.
    When shapefile_field is None, returns immediately (auto-detection will occur later).

    Parameters
    ----------
    shapefile_path : Path
        Path to the .shp file
    shapefile_field : str or None
        The --shp-field value. None means auto-detect (always valid).

    Raises
    ------
    SystemExit
        If the field is specified but not found in the shapefile.
    """
    if shapefile_field is None:
        return  # Auto-detection will handle this

    try:
        import geopandas as gpd
        gdf = gpd.read_file(shapefile_path)
        available = [c for c in gdf.columns if c != gdf.geometry.name]
        if shapefile_field not in available:
            available_str = ', '.join(available) if available else '(none — shapefile may be invalid)'
            suggestion = (
                f"\nRe-run with one of those fields, e.g.:\n"
                f"  --shp-field {available[0]}"
            ) if available else ""
            print(
                f"\nError: --shp-field '{shapefile_field}' not found in shapefile.\n"
                f"  Shapefile: {shapefile_path}\n"
                f"  Available fields: {available_str}"
                f"{suggestion}",
                file=sys.stderr
            )
            sys.exit(1)
    except ImportError:
        pass  # geopandas not available; geographic.py will raise a clear error later
    except SystemExit:
        raise  # Don't catch our own sys.exit()
    except Exception as e:
        print(
            f"\nWarning: Could not pre-validate shapefile field: {e}\n"
            f"Pipeline will attempt geographic analysis and report any errors.",
            file=sys.stderr
        )


def _print_geo_failure_to_stderr(shapefile_path, error) -> None:
    """Print a visible geographic analysis failure message to stderr."""
    print(
        f"\n{'='*60}\n"
        f"WARNING: Geographic analysis was skipped.\n"
        f"  Shapefile: {shapefile_path}\n"
        f"  Reason: {error}\n"
        f"\nNo distribution maps or geographic plots will be generated.\n"
        f"To fix: check available shapefile fields with geopandas:\n"
        f"  python -c \"import geopandas as gpd; "
        f"gdf = gpd.read_file('{shapefile_path}'); print(gdf.columns.tolist())\"\n"
        f"Then re-run with: --shp-field <correct_field_name>\n"
        f"{'='*60}\n",
        file=sys.stderr
    )


def mappable_coordinates(frame: pd.DataFrame, mappable_ids) -> pd.DataFrame:
    """
    Copy of frame for the maps and the plot kit: lat/lon are blanked for records
    whose coordinates are not of geographic quality (centroid, zero, out of range;
    P0-4). Those records keep lat/lon in records.csv but are not mapped.
    """
    out = frame.copy()
    if {'processid', 'lat', 'lon'} <= set(out.columns):
        out.loc[~out['processid'].isin(mappable_ids), ['lat', 'lon']] = np.nan
    return out


def region_column_clash(input_columns, geo_category: str) -> Optional[str]:
    """
    Warning text when the region column (--geo-category) has the name of an
    input field, which it then replaces; None otherwise.
    """
    if geo_category not in set(input_columns):
        return None
    return (f"--geo-category '{geo_category}' is also the name of an input field: the region "
            f"column replaces that field in records.csv and the metadata summary. Use another "
            f"name (e.g. '{geo_category}_region') to keep both.")


def _build_tree(cfg, organism, dirs, haplotype_fasta, haplotype_stats, prov, step_warnings):
    """
    Build and root the haplotype tree (pipeline step 5.3; taxon evidence needs it
    before any step reads a taxon). Returns (tree path or None, intermediate
    prefix, output prefix). Failures are recorded as warnings, never raised.
    """
    intermediate_prefix = dirs['intermediate_phylo'] / organism
    output_prefix = dirs['phylogenetic'] / organism
    if not cfg.phylogenetic.build_tree:
        logger.info("  Phylogenetic tree building disabled in config")
        return None, intermediate_prefix, output_prefix
    try:
        logger.info("5.3: Building phylogenetic tree from haplotypes...")
        if not utils.check_external_tool("mafft"):
            step_warnings.record("5 tree", "MAFFT not found; tree not built")
            return None, intermediate_prefix, output_prefix
        if not utils.check_external_tool("fasttree"):
            step_warnings.record("5 tree", "FastTree not found; tree not built")
            return None, intermediate_prefix, output_prefix

        # Choose tips explicitly (P0-7): every haplotype with at least
        # --tree-min-members members and the minimum length, except
        # suspect haplotypes. Exclusions are listed with their reasons.
        records = list(SeqIO.parse(haplotype_fasta, "fasta"))
        suspect = {}
        if 'is_suspect' in haplotype_stats.columns:
            flagged = haplotype_stats[haplotype_stats['is_suspect'].fillna(False).astype(bool)]
            reasons = flagged['suspect_reasons'] if 'suspect_reasons' in flagged.columns else None
            suspect = {
                hid: (reasons.loc[i] if reasons is not None and isinstance(reasons.loc[i], str)
                      and reasons.loc[i] else 'flagged')
                for i, hid in flagged['haplotype_id'].items()
            }
        kept, excluded = phylogenetics.select_tree_haplotypes(
            records, cfg.phylogenetic.min_cluster_size,
            cfg.phylogenetic.min_consensus_length, suspect)
        excluded_path = dirs['phylogenetic'] / f"{organism}_excluded_tips.txt"
        with open(excluded_path, 'w') as f:
            f.write("haplotype_id\treason\n")
            for hid, reason in excluded:
                f.write(f"{hid}\t{reason}\n")
        tree_input = dirs['intermediate_phylo'] / f"{organism}_tree_input.fasta"
        SeqIO.write(kept, tree_input, "fasta")
        logger.info(f"  Tree tips: {len(kept)} of {len(records)} haplotypes "
                    f"({len(excluded)} excluded; see {excluded_path.name}); "
                    f"--tree-min-members {cfg.phylogenetic.min_cluster_size}, "
                    f"minimum length {cfg.phylogenetic.min_consensus_length} bp")
        prov.applied.update(
            tree_min_members=cfg.phylogenetic.min_cluster_size,
            tree_min_length_bp=cfg.phylogenetic.min_consensus_length,
            n_tree_tips=len(kept),
            n_tree_excluded=len(excluded),
        )

        # Tips are already selected, so build_phylogeny's own filter keeps everything.
        tree = phylogenetics.build_phylogeny(
            consensus_fasta=str(tree_input),
            output_prefix=str(intermediate_prefix),
            threads=cfg.n_threads,
            min_consensus_length=0,
            min_cluster_size=1,
            outgroup_fasta=str(cfg.phylogenetic.outgroup_fasta) if cfg.phylogenetic.outgroup_fasta else None,
            trim_alignment=cfg.phylogenetic.trim_alignment,
            trim_method=cfg.phylogenetic.trim_method
        )
        intermediate_tree = f"{intermediate_prefix}_tree.nwk"
        tree_path = f"{output_prefix}_tree.nwk"
        if Path(intermediate_tree).exists():
            Path(intermediate_tree).rename(tree_path)
        if tree is None or not Path(tree_path).exists():
            step_warnings.record("5 tree", f"tree building finished but {tree_path} was not written")
            return None, intermediate_prefix, output_prefix

        # Root on the outgroup when one is given, else at the midpoint;
        # support values stay on their branches (phylogenetics.root_on_tips)
        phy = cfg.phylogenetic
        rooting, rooted = {"method": "none"}, None
        if phy.outgroup_label:
            rooting = {"method": "outgroup_label", "outgroup": phy.outgroup_label}
            rooted = phylogenetics.reroot_tree_by_label(
                tree_file=tree_path, outgroup_label=phy.outgroup_label, output_file=tree_path)
        elif phy.outgroup_taxon:
            rooting = {"method": "outgroup_taxon", "outgroup": phy.outgroup_taxon}
            rooted = phylogenetics.reroot_tree_by_taxon(
                tree_file=tree_path,
                taxonomy_csv=str(dirs['taxonomy'] / f"{organism}_haplotype_taxonomy.csv"),
                taxon_query=phy.outgroup_taxon, output_file=tree_path)
        elif phy.outgroup_fasta:
            rooting = {"method": "outgroup_fasta", "outgroup": str(phy.outgroup_fasta)}
            rooted = phylogenetics.reroot_tree_by_tips(
                tree_path, [r.id for r in SeqIO.parse(str(phy.outgroup_fasta), "fasta")],
                output_file=tree_path, description="outgroup FASTA")
        if rooting["method"] != "none" and rooted is None:
            step_warnings.record(
                "5 tree reroot",
                f"rooting on {rooting['method']} '{rooting['outgroup']}' failed"
                + ("; midpoint-rooted instead" if phy.midpoint_root else "; tree left unrooted"))
            rooting = {"method": "none"}
        if rooting["method"] == "none" and phy.midpoint_root:
            if phylogenetics.midpoint_root_tree_file(tree_path):
                rooting = {"method": "midpoint"}
        prov.applied["tree_rooting"] = rooting
        logger.info(f"  Built phylogenetic tree: {tree_path}")
        return tree_path, intermediate_prefix, output_prefix
    except Exception as e:
        step_warnings.record("5 tree", e)
        return None, intermediate_prefix, output_prefix


def _add_input_row_arguments(parser: argparse.ArgumentParser) -> None:
    """--marker and --allow-malformed-rows, shared by boldgenotyper and boldgenotyper-sweep."""
    parser.add_argument(
        '--marker',
        default='COI-5P',
        metavar='MARKER',
        help='Marker to analyze (BOLD marker_code; default: COI-5P). Rows of other markers '
             'are removed before duplicate processids, and counted.'
    )
    parser.add_argument(
        '--allow-malformed-rows',
        action='store_true',
        help='Drop data rows whose number of fields differs from the header\'s, and count '
             'them, instead of stopping (default: stop and name the first such line).'
    )
    parser.add_argument(
        '--label-map',
        type=Path,
        default=None,
        metavar='FILE',
        help='CSV/TSV of confirmed label replacements (from, to, reason; optional source), '
             'applied before labels are pooled and recorded in run_provenance.json. '
             'label_variants.csv lists candidates; nothing is merged without this file.'
    )


def _write_error_filter_diagnostics(records_df, df_input, details, output_dir, summary_csv,
                                    prov, cfg) -> None:
    """sources.csv, the error-filter summary values and their warnings (pre-tag A4, A4b)."""
    efd = error_filter_diagnostics
    groups = details.haplotype_only_groups if details is not None else None
    groups = groups if groups is not None else pd.DataFrame(columns=efd.GROUP_COLUMNS)
    by_record = efd.record_groups(df_input).merge(records_df[["processid", "status"]], on="processid")
    efd.sources_table(groups, by_record).to_csv(output_dir / "sources.csv", index=False)
    statuses = records_df["status"]
    share = efd.filtered_share(int((statuses == "assigned").sum()),
                               int((statuses == "assigned_corrected").sum()),
                               int((statuses == "excluded_error_singleton").sum()))
    n_absent = int((records_df["label_in_retained"] == False).sum())  # noqa: E712 (object column)
    group_values = efd.summary_values(groups, details.n_singletons_exempt if details is not None else 0)
    if group_values["n_haplotype_only_groups"]:
        logger.warning(f"  {group_values['n_haplotype_only_groups']} source-and-label group(s) deposit one "
                       f"sequence per haplotype (sources.csv); {group_values['n_singletons_exempt']} of "
                       f"their singletons were kept as haplotypes, not treated as errors.")
    if share > efd.FILTERED_SHARE_WARNING:
        logger.warning(f"  {share * 100:.1f}% of the records the error filter saw were treated as errors "
                       f"(corrected or left unassigned); if studies deposit one sequence per haplotype, "
                       f"true rare haplotypes may have been treated as errors (see sources.csv).")
    if n_absent:
        logger.warning(f"  {n_absent} error-filtered records carry a species label found in no "
                       f"retained haplotype; the filter may have removed species (records.csv "
                       f"label_in_retained).")
    prov.data_notes["error_filter"] = {
        "pct_error_filtered": round(share * 100, 2),
        "n_singleton_label_absent": n_absent,
        **group_values,
        "filtered_share_warning": efd.FILTERED_SHARE_WARNING,
        "haplotype_only_min_records": cfg.haplotype.haplotype_only_min_records,
        "haplotype_only_min_distinct": cfg.haplotype.haplotype_only_min_distinct,
    }
    updated = reports.update_summary(summary_csv, {"n_singleton_label_absent": n_absent, **group_values})
    if updated:
        prov.applied["assignment_summary"] = updated
    logger.info(f"  sources.csv: {len(groups)} source-and-label groups, "
                f"{group_values['n_haplotype_only_groups']} depositing one sequence per haplotype")


def run_pipeline(
    tsv_path: Path,
    organism: str,
    output_dir: Path,
    cfg: config.PipelineConfig,
    no_report: bool = False,
    skip_geo: bool = False,
    export_plot_data: bool = False,
    export_popgen_formats: Optional[list] = None,
    metadata_summary_enabled: bool = False,
    metadata_fields: Optional[list] = None,
    normalize_sex: bool = True,
    temporal_analysis: bool = False,
    shapefile_path: Optional[Path] = None,
    shapefile_field: Optional[str] = None,
    geo_category: str = 'ocean_basin',
    argv: Optional[list] = None,
    download_date: Optional[str] = None,
    taxon_curation: Optional[Path] = None,
    marker: Optional[str] = 'COI-5P',
    allow_malformed_rows: bool = False,
    label_map: Optional[Path] = None,
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
        Write plots/: plot data and regeneration scripts (default: False; --plot-kit)
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
    argv : list, optional
        Command line recorded in run_provenance.json (default: sys.argv)
    download_date : str, optional
        BOLD download date (YYYY-MM-DD) recorded in run_provenance.json
    taxon_curation : Path, optional
        Curation file of accepted haplotype taxa (see taxon_curation), applied
        after the label rule
    marker : str, optional
        Target marker (BOLD marker_code; default 'COI-5P'). Rows of other
        markers are removed before duplicate processids and counted.
    allow_malformed_rows : bool, optional
        Drop and count rows whose field count differs from the header's
        instead of stopping (default: False)
    label_map : Path, optional
        Confirmed label map (label_variants.load_label_map): deposited
        species labels replaced before labels are pooled; recorded by md5

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

    # A curation file that cannot be used stops the run before it writes anything
    label_map_df = None
    if label_map:
        try:
            label_map_df = label_variants.load_label_map(label_map)
        except label_variants.LabelMapError as e:
            logger.error(f"Label map: {e}")
            return False

    curation = None
    if taxon_curation:
        try:
            curation = taxon_curation_mod.load_curation(taxon_curation)
        except taxon_curation_mod.CurationError as e:
            logger.error(f"Taxon curation file: {e}")
            return False

    # Look for an earlier run before this run creates any folders
    previous_run = provenance.find_previous_run(output_dir)

    # Setup directory structure
    dirs = setup_directories(output_dir)

    # Run provenance: the fully resolved configuration, environment and input
    # checksum (P0-1). Its StepWarnings collects step and figure failures (P0-6).
    prov = provenance.RunProvenance(cfg, argv if argv is not None else sys.argv, tsv_path,
                                    output_dir, download_date=download_date)
    prov.applied.update(
        geographic_analysis=not skip_geo,
        geo_category=geo_category,
        custom_shapefile=str(shapefile_path) if shapefile_path else None,
        shapefile_field=shapefile_field,
        run_options={
            'html_report': not no_report,
            'export_plot_data': export_plot_data,
            'export_popgen_formats': export_popgen_formats,
            'metadata_summary': metadata_summary_enabled,
            'metadata_fields': metadata_fields,
            'normalize_sex': normalize_sex,
            'temporal_analysis': temporal_analysis,
        },
    )
    step_warnings = prov.warnings.attach()
    prov_path = output_dir / "run_provenance.json"
    if previous_run:
        step_warnings.record(
            "0 output folder",
            f"Output folder already holds a previous run ({', '.join(previous_run)}). "
            "Files from that run that this run does not overwrite remain and may be stale; "
            "use a new --output folder for a clean run.")
    prov.write(prov_path)
    logger.info(f"Saved run provenance to {prov_path}")

    # ========================================================================
    # PHASE 1: Data Loading and Quality Control
    # ========================================================================
    logger.info("PHASE 1: Data Loading and Quality Control")
    logger.info("-" * 80)

    try:
        # Parse BOLD TSV
        logger.info("1.1: Parsing BOLD TSV metadata...")
        df = metadata.parse_bold_tsv(tsv_path, marker=marker,
                                     allow_malformed_rows=allow_malformed_rows)
        input_rows = {"marker": marker}
        input_rows.update({k: df.attrs[k] for k in (
            "n_other_marker_rows", "n_malformed_rows", "malformed_lines", "n_duplicate_rows")})
        prov.applied["input_rows"] = input_rows
        n_cleaned = 0
        for col in ("species", "genus"):
            if col in df.columns:
                given = df[col]
                cleaned = given.where(given.isna(), given.map(label_variants.clean_label))
                n_cleaned += int((cleaned.fillna("") != given.fillna("")).sum())
                df[col] = cleaned
        prov.data_notes["n_labels_cleaned"] = n_cleaned
        if label_map_df is not None:
            df, usage = label_variants.apply_label_map(df, label_map_df)
            prov.applied["label_map"] = {
                "path": str(label_map),
                "md5": provenance.file_md5(label_map),
                "n_rows": len(label_map_df),
                "n_records_changed": int(usage["n_records"].sum()),
                "unused": usage.loc[usage["n_records"] == 0, "from"].tolist(),
            }
            logger.info(f"  Label map: {int(usage['n_records'].sum())} record label(s) replaced "
                        f"({len(label_map_df)} map rows)")
        logger.info(f"  Loaded {len(df)} samples with {len(df.columns)} columns")
        clash = None if skip_geo else region_column_clash(df.columns, geo_category)
        if clash:
            logger.warning(f"  {clash}")

        # Save parsed metadata
        parsed_tsv = dirs['intermediate'] / "01_parsed_metadata.tsv"
        df.to_csv(parsed_tsv, sep='\t', index=False)

        # Parse coordinates for every record and mark their quality. Records are
        # never removed here; only region assignment and maps use the flag (P0-4).
        logger.info("1.2: Parsing coordinates and marking coordinate quality...")
        df = metadata.mark_coordinate_quality(metadata.parse_coord_column(df))
        mappable_ids = set(df.loc[df['is_geographic_quality'], 'processid'])
        prov.data_notes["coord_status"] = {
            str(k): int(v) for k, v in df['coord_status'].value_counts().items()
        }
        n_centroid = int((df['coord_status'] == 'centroid').sum())
        if n_centroid:
            logger.info(f"  {n_centroid} records have country- or province-centroid coordinates; "
                        f"they keep lat/lon but are not assigned a region ({geo_category}) or mapped")

        # Save metadata with quality markers
        marked_tsv = dirs['intermediate'] / "02_marked_metadata.tsv"
        df.to_csv(marked_tsv, sep='\t', index=False)

        # Assign geographic regions (only to samples with geographic quality).
        # A record without a region has NaN (empty in CSV), never 'Unknown'.
        logger.info(f"1.3: Assigning geographic regions ({geo_category})...")
        geo_analysis_performed = False  # Track whether geographic analysis was successful
        df_with_basins = df.copy()
        df_with_basins[geo_category] = np.nan

        def _attach_regions(assigned_df, region_col):
            """Copy region assignments onto df_with_basins by processid."""
            regions = assigned_df.set_index('processid')[region_col].replace('Unknown', np.nan)
            df_with_basins[geo_category] = df_with_basins['processid'].map(regions)
            n_quality = int(df['is_geographic_quality'].sum())
            n_assigned = int(df_with_basins[geo_category].notna().sum())
            logger.info(f"  Assigned {geo_category} to {n_assigned} samples")
            logger.info(f"  {n_quality - n_assigned} geographic-quality samples fall outside all regions; "
                        f"{len(df) - n_quality} samples lack geographic-quality coordinates "
                        f"(missing, unparseable, centroid, zero or out of range)")

        if skip_geo:
            logger.info("  Geographic analysis disabled (--no-geo flag)")
        elif shapefile_path:
            # Use custom shapefile
            logger.info(f"  Using custom shapefile: {shapefile_path}")
            logger.info(f"  Shapefile field: {shapefile_field}")
            logger.info(f"  Geographic category: {geo_category}")

            try:
                df_geo_quality = df[df['is_geographic_quality']].copy()
                logger.info(f"  {len(df_geo_quality)}/{len(df)} samples have geographic-quality coordinates")

                if len(df_geo_quality) > 0:
                    df_geo_assigned = geographic.assign_regions_from_shapefile(
                        df_geo_quality,
                        shapefile_path=shapefile_path,
                        shapefile_field=shapefile_field,
                        output_column=geo_category,
                        coord_col="coord"
                    )
                    _attach_regions(df_geo_assigned, geo_category)
                    geo_analysis_performed = True
                else:
                    logger.warning("  No samples with geographic-quality coordinates")

            except Exception as e:
                step_warnings.record("1.3 region assignment",
                                     f"Custom shapefile failed: {e}; continuing without geographic analysis")
                _print_geo_failure_to_stderr(shapefile_path, e)
                df_with_basins[geo_category] = np.nan

        elif not cfg.geographic.goas_shapefile_path.exists():
            logger.warning(f"  GOaS shapefile not found at: {cfg.geographic.goas_shapefile_path}")
            logger.warning("")
            logger.warning("  To enable geographic analysis, download the GOaS shapefile manually:")
            logger.warning("  1. Visit: https://www.marineregions.org/downloads.php")
            logger.warning("  2. Download: GOaS_v1_20211214.zip (registration form required)")
            logger.warning(f"  3. Extract to: {cfg.geographic.goas_shapefile_path.parent}")
            logger.warning("  4. Verify with: python -m boldgenotyper.goas_downloader")
            logger.warning("")
            logger.warning("  Pipeline will continue without geographic analysis...")
            step_warnings.record("1.3 region assignment",
                                 f"GOaS shapefile not found at {cfg.geographic.goas_shapefile_path}; "
                                 "geographic analysis skipped")
        else:
            # Use default GOaS shapefile
            logger.info("  Using default GOaS ocean basin shapefile")
            try:
                goas_data = geographic.load_goas_data(cfg.geographic.goas_shapefile_path)

                df_geo_quality = df[df['is_geographic_quality']].copy()
                logger.info(f"  {len(df_geo_quality)}/{len(df)} samples have geographic-quality coordinates")

                if len(df_geo_quality) > 0:
                    df_geo_assigned = geographic.assign_ocean_basins(
                        df_geo_quality, goas_data=goas_data, coord_col="coord"
                    )
                    _attach_regions(df_geo_assigned, 'ocean_basin')
                    geo_analysis_performed = True
                else:
                    logger.warning("  No samples with geographic-quality coordinates")

            except Exception as e:
                step_warnings.record("1.3 region assignment",
                                     f"GOaS data failed: {e}; continuing without geographic analysis")
                df_with_basins[geo_category] = np.nan

        basins_tsv = dirs['intermediate_geographic'] / f"samples_with_{geo_category}.tsv"
        df_with_basins.to_csv(basins_tsv, sep='\t', index=False)

    except metadata.MalformedTSVError as e:
        logger.error(f"Phase 1 stopped: {e}")
        step_warnings.detach()
        prov.finish(success=False)
        prov.write(prov_path)
        return False
    except Exception as e:
        logger.error(f"Phase 1 failed: {e}", exc_info=True)
        step_warnings.detach()
        prov.finish(success=False)
        prov.write(prov_path)
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
        _code_names = {2: "vertebrate mitochondrial", 5: "invertebrate mitochondrial"}
        _code_label = _code_names.get(cfg.coi_validation.mitochondrial_code, "NCBI genetic code")
        logger.info(f"  Using genetic code: {cfg.coi_validation.mitochondrial_code} ({_code_label})")
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

        df_qc_all, qc_stats = quality_control.apply_dynamic_qc_filters(
            df=df_with_basins,
            sequences_dict=corrected_sequences,
            orf_stats_df=orf_stats_df,
            min_raw_length_abs=cfg.qc.min_raw_length_abs,
            min_raw_length_frac_of_median=cfg.qc.min_raw_length_frac_of_median,
            max_raw_N_fraction=cfg.qc.max_raw_N_fraction,
            require_valid_orf=cfg.qc.require_valid_orf,
            processid_col='processid',
            return_all=True,
        )
        df_qc_all['qc_pass'] = df_qc_all['qc_pass'].fillna(False).astype(bool)
        df_qc_passed = df_qc_all[df_qc_all['qc_pass']].copy()

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
                header = processid  # bare processid, as in every other table (P0-10)
                sequence = corrected_sequences[processid]
                qc_passed_fasta_records.append((header, sequence))

        qc_passed_fasta = dirs['qc'] / f"{organism}_qc_passed.fasta"
        utils.write_fasta(qc_passed_fasta_records, qc_passed_fasta)
        # Kept in the run folder for boldgenotyper-query --search-records (Q10)
        utils.write_fasta(qc_passed_fasta_records, output_dir / "records_qc_passed.fasta")
        logger.info(f"  Created {len(qc_passed_fasta_records)} QC-passed FASTA records")
        logger.info(f"  Saved: {qc_passed_fasta}")

        # Log QC summary
        logger.info("QC Summary:")
        logger.info(f"  Total input sequences: {qc_stats['n_total']}")
        logger.info(f"  Passed all QC filters: {qc_stats['n_pass']} ({qc_stats['pass_rate']*100:.1f}%)")
        logger.info(f"  Failed QC filters: {qc_stats['n_fail']} ({(1-qc_stats['pass_rate'])*100:.1f}%)")

    except Exception as e:
        logger.error(f"Phase 2 failed: {e}", exc_info=True)
        step_warnings.detach()
        prov.finish(success=False)
        prov.write(prov_path)
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
        haplotype_records, haplotype_mapping, haplotype_stats, details = dereplication.identify_haplotypes(
            tsv_path=tsv_path,
            fasta_path=qc_passed_fasta,
            output_dir=dirs['haplotype_discovery'],
            min_core_coverage=cfg.core_region.core_min_coverage,
            gap_mask_threshold=cfg.core_region.mask_gap_threshold,
            min_core_length=cfg.core_region.core_min_length,
            min_singleton_distance=cfg.haplotype.min_singleton_distance,
            max_singleton_distance=cfg.haplotype.max_singleton_distance,
            orf_validation_df=orf_stats_df,
            cleanup_intermediates=False,
            min_sequence_coverage=cfg.core_region.min_sequence_coverage,
            assign_corrected=cfg.haplotype.assign_corrected,
            return_details=True,
            record_groups=error_filter_diagnostics.record_groups(df_qc_passed),
            haplotype_only_min_records=cfg.haplotype.haplotype_only_min_records,
            haplotype_only_min_distinct=cfg.haplotype.haplotype_only_min_distinct,
        )

        logger.info(f"  Identified {len(haplotype_records)} high-confidence haplotypes (after error filtering)")
        prov.applied.update(
            core_region_length_bp=details.core_region_length_bp,
            core_reference_length_bp=details.core_reference_length_bp,
            min_sequence_length_bp=details.min_sequence_length_bp,
            n_excluded_core_coverage=len(details.excluded_core_coverage),
            n_error_singleton_records=len(details.error_singletons),
            n_excluded_error_tied=details.n_excluded_error_tied,
            n_excluded_error_far=details.n_excluded_error_far,
            n_haplotypes=len(haplotype_records),
        )

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
        step_warnings.detach()
        prov.finish(success=False)
        prov.write(prov_path)
        return False

    # ========================================================================
    # PHASE 4: Haplotype Assignment
    # ========================================================================
    logger.info("")
    logger.info("PHASE 4: Haplotype Assignment")
    logger.info("-" * 80)

    try:
        logger.info("4.1: Assigning a status to every input record...")
        logger.info(f"  Using {len(df_qc_passed)} QC-passed samples")
        logger.info("  Method: direct mapping from ESV discovery; error-filtered singletons "
                    + ("assigned to their nearest retained haplotype"
                       if cfg.haplotype.assign_corrected else "left unassigned (--no-assign-corrected)"))

        annotated_tsv = dirs['intermediate_assignments'] / f"{organism}_with_haplotypes.tsv"
        diagnostics_csv = dirs['assignments'] / f"{organism}_diagnostics.csv"

        # One status per unique input record (P0-2)
        status_df = record_status.build_record_status(
            unique_ids=list(df_with_basins['processid']),
            sequence_ids=set(sequences_dict),
            qc=df_qc_all[['processid', 'qc_pass', 'qc_fail_reason']],
            exact_mapping=haplotype_mapping.loc[haplotype_mapping['assignment'] == 'exact',
                                                ['processid', 'haplotype_id']],
            excluded_core=details.excluded_core_coverage,
            error_singletons=details.error_singletons,
            singleton_distance=cfg.haplotype.min_singleton_distance,
        )

        # Attach statuses and haplotype membership to the QC-passed records
        df_with_haplotypes = df_qc_passed.merge(
            status_df, on='processid', how='left', validate='one_to_one'
        ).merge(
            haplotype_mapping[['processid', 'seq_md5', 'n_members', 'is_singleton', 'assignment']],
            on='processid', how='left', validate='one_to_one'
        )
        df_with_haplotypes['assigned_haplotype'] = df_with_haplotypes['haplotype_id']

        # Save annotated data
        df_with_haplotypes.to_csv(annotated_tsv, sep="\t", index=False)
        logger.info(f"  Wrote annotated metadata: {annotated_tsv}")

        # Diagnostics: every unique record, with its status
        diagnostics_df = status_df.rename(columns={'haplotype_id': 'assigned_haplotype'})
        diagnostics_df['assignment_method'] = diagnostics_df['status'].map({
            'assigned': 'ESV_direct_mapping',
            'assigned_corrected': 'nearest_retained_haplotype',
        }).fillna('')
        diagnostics_df.to_csv(diagnostics_csv, index=False)
        logger.info(f"  Wrote diagnostics: {diagnostics_csv}")

        status_counts = status_df['status'].value_counts().reindex(record_status.STATUSES, fill_value=0)
        logger.info("  Record statuses:")
        for status, n in status_counts.items():
            logger.info(f"    {status}: {n}")

    except Exception as e:
        logger.error(f"Phase 4 failed: {e}", exc_info=True)
        step_warnings.detach()
        prov.finish(success=False)
        prov.write(prov_path)
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
        assign_table["assignment_source"] = "labels"
        if curation is not None:
            hap_meta = pd.DataFrame({
                "haplotype_id": [r.id for r in haplotype_records],
                "seq_md5": [haplotype_ids.seq_md5(str(r.seq)) for r in haplotype_records]})
            applied_rows, not_applied_rows = taxon_curation_mod.match_curation(
                curation, hap_meta, {r.id: str(r.seq) for r in haplotype_records})
            assign_table = taxon_curation_mod.apply_curation(assign_table, applied_rows)
            prov.applied["taxon_curation"] = {
                "path": str(taxon_curation),
                "md5": provenance.file_md5(taxon_curation),
                "n_rows": len(curation),
                "n_applied": len(applied_rows),
                "n_not_applied": len(not_applied_rows),
                "not_applied": not_applied_rows[["row", "taxon", "reason_not_applied"]]
                .to_dict("records"),
            }
            logger.info(f"  Curated taxa applied to {len(applied_rows)} haplotype(s)")
            for r in not_applied_rows.to_dict("records"):
                logger.warning(f"  Curation row {r['row']} ({r['taxon']}) not applied: "
                               f"{r['reason_not_applied']}")

        # Build haplotype_sp labels (e.g., "Species name H0001")
        def _join_label(row):
            sp = row["assigned_sp"]
            hap_id = row["haplotype_id"]
            if not sp:
                return hap_id if isinstance(hap_id, str) else ""
            return f"{sp} {hap_id}"

        def _write_taxonomy():
            assign_table["haplotype_sp"] = assign_table.apply(_join_label, axis=1)
            species_counts.to_csv(dirs['taxonomy'] / f"{organism}_species_by_haplotype.csv", index=False)
            assign_table.to_csv(dirs['taxonomy'] / f"{organism}_haplotype_taxonomy.csv", index=False)

        _write_taxonomy()                      # read by outgroup-taxon rooting

        # 5.3: the tree (taxon evidence needs it before any step reads a taxon)
        tree_path, intermediate_prefix, output_prefix = _build_tree(
            cfg, organism, dirs, haplotype_fasta, haplotype_stats, prov, step_warnings)

        # 5.4: taxon evidence, then evidence taxa (design addendum of 2026-10-05).
        # Evidence uses the label rule's and curated taxa only (one pass).
        evidence, evidence_not_applied = None, None
        try:
            logger.info("5.4: Taxon evidence...")
            evidence_records = layout.build_records_table(
                df_with_basins, status_df, df_qc_all, assign_table, geo_category)
            evidence = taxon_evidence.compute_evidence(
                assign_table, evidence_records, {r.id: str(r.seq) for r in haplotype_records},
                tree_path)
            assign_table, evidence_not_applied = taxon_evidence.apply_evidence(
                assign_table, evidence, apply=cfg.taxonomy.apply_evidence)
            n_evidence = int((assign_table["assignment_source"] == taxon_evidence.EVIDENCE_SOURCE).sum())
            prov.applied["taxon_evidence"] = {
                "apply": bool(cfg.taxonomy.apply_evidence),
                "tree": tree_path is not None,
                "n_applied": n_evidence,
                "n_suggested_not_applied": len(evidence_not_applied),
                "suggestions_file": (taxon_evidence.SUGGESTIONS_FILE
                                     if len(evidence_not_applied) else None),
            }
            logger.info(f"  Evidence taxa: {n_evidence} haplotype(s) set by evidence; "
                        f"{len(evidence_not_applied)} suggestion(s) not applied")
        except Exception as e:
            step_warnings.record("5 taxon evidence", e)
            logger.debug("Taxon evidence error details:", exc_info=True)
        _write_taxonomy()

        # Merge taxonomy back into main dataframe
        df_with_haplotypes = df_with_haplotypes.merge(
            assign_table[["haplotype_id", "assigned_sp", "haplotype_sp", "assignment_level",
                          "assignment_notes", "assignment_source", "majority_fraction"]],
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

        logger.info(f"  Assigned taxonomy to {len(assign_table)} haplotypes")

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
            logger.info(f"  Saved annotated dataset: {annotated_csv}")
        else:
            step_warnings.record("5 annotated table", f"Failed to save annotated dataset: {annotated_csv}")

    except Exception as e:
        logger.error(f"Phase 5 failed: {e}", exc_info=True)
        step_warnings.detach()
        prov.finish(success=False)
        prov.write(prov_path)
        return False

    # ========================================================================
    # PHASE 6: Post-assignment Quality Control
    # ========================================================================
    logger.info("")
    logger.info("PHASE 6: Post-assignment Quality Control")
    logger.info("-" * 80)

    try:
        logger.info("6.1: Adding record labels and label conflicts to the annotated table...")

        # Record label, haplotype majority label, label conflicts, depositor-note keyword hits
        df_final_enhanced = quality_control.add_label_columns(
            df_final,
            group_col="haplotype_id",
            species_col="species",
            genus_col="genus",
            notes_col="notes" if "notes" in df_final.columns else None,
            curated_taxa=taxon_curation_mod.taxa_not_from_labels(assign_table),
        )

        # Save enhanced annotated file
        annotated_csv_enhanced = dirs['base'] / f"{organism}_annotated.csv"
        df_final_enhanced.to_csv(annotated_csv_enhanced, index=False, quoting=1)
        logger.info(f"  Annotated dataset with label columns: {annotated_csv_enhanced}")

        # Update df_final to use the enhanced version
        df_final = df_final_enhanced

        logger.info("6.2: Writing the label-conflict report...")
        qc_results = quality_control.generate_quality_control_report(
            df_final,
            output_dir=dirs['quality_control'],
            organism=organism,
            group_col="haplotype_id",
            species_col="species",
            haplotype_sequences={r.id: str(r.seq) for r in haplotype_records},
        )

        logger.info(f"  Label report saved to: {dirs['quality_control']}")

        # Print quality alerts to console
        quality_control.print_quality_alert(qc_results)

    except Exception as e:
        step_warnings.record("6 post-assignment QC", e)
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
            step_warnings.record("6.5 species aggregation", "skipped: taxonomy assignments not found")
        else:
            # Create species_composition DataFrame from taxonomy data
            # Map from assign_table columns to expected species_composition columns
            species_composition = assign_table.copy()
            species_composition = species_composition.rename(columns={
                'assigned_sp': 'primary_species',
                'majority_fraction': 'primary_species_pct'
            })

            # A haplotype counts toward its species when it is assigned at species level,
            # by labels, evidence or curation (a curated or evidence taxon can have a low
            # majority_fraction)
            species_composition['is_ambiguous'] = (
                species_composition['assignment_level'] != 'species').astype(bool)

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
                output_dir=dirs['species_analysis']
            )

            if cfg.extended:
                logger.info("6.5.2: Calculating haplotype diversity per species label...")
                species_analysis.calculate_species_diversity(
                    species_assignments=species_assignments,
                    output_dir=dirs['species_analysis']
                )

            logger.info("6.5.3: Generating species-level geographic summary...")
            species_geographic = species_analysis.generate_species_geographic_summary(
                species_assignments=species_assignments,
                output_dir=dirs['species_analysis']
            )

            # Per-species record files are written from records.csv at the end (--extended)

            logger.info(f"  Species-level analysis complete: {dirs['species_analysis']}")

    except Exception as e:
        step_warnings.record("6.5 species aggregation", e)
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
                logger.info(f"  Updated annotated dataset with geographic enhancements")

            logger.info(f"  Geographic enhancement complete: {dirs['geographic_analysis']}")

        except Exception as e:
            step_warnings.record("7 geographic enhancement", e)
            logger.debug("Geographic enhancement error details:", exc_info=True)

    # ========================================================================
    # PHASE 8: Phylogenetic Analysis (Optional)
    # ========================================================================
    logger.info("")
    logger.info("PHASE 8: Phylogenetic Analysis")
    logger.info("-" * 80)

    # The tree is built in step 5.3 (_build_tree); this phase writes its other outputs
    if tree_path:
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
                logger.info(f"  Created relabeled tree: {relabeled_tree_path}")
                logger.info(f"  Created relabeled alignment (intermediate): {relabeled_alignment_path}")
            else:
                step_warnings.record("8.2 tree relabeling", "alignment or taxonomy file not found")
        except Exception as e:
            step_warnings.record("8.2 tree relabeling", e)

        # MSA chunk plots: --extended only (P2)
        if cfg.msa.enabled and cfg.extended:
            logger.info("8.3: Generating MSA visualization...")
            try:
                # Use the alignment file that exists (prefer trimmed if available)
                alignment_file = None
                if Path(f"{intermediate_prefix}_aligned_trimmed.fasta").exists():
                    alignment_file = Path(f"{intermediate_prefix}_aligned_trimmed.fasta")
                elif Path(f"{intermediate_prefix}_aligned.fasta").exists():
                    alignment_file = Path(f"{intermediate_prefix}_aligned.fasta")

                if importlib.util.find_spec("pymsaviz") is None:
                    # msa_visualization swallows its own ImportError;
                    # check first so the optional extra is not a warning
                    raise ImportError("pymsaviz")
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
                        logger.info(f"  Generated {len(msa_plots)} MSA plot(s)")
                    else:
                        step_warnings.record("8.3 MSA visualization", "no MSA plots generated")
                else:
                    step_warnings.record("8.3 MSA visualization", "alignment file not found")
            except ImportError:
                # Optional dependency (extras_require["msa"]): not a step failure
                logger.info("  MSA visualization skipped: pymsaviz is not installed "
                            "(optional; pip install pymsaviz>=0.4.0)")
            except Exception as e:
                step_warnings.record("8.3 MSA visualization", e)
                logger.debug("Full traceback:", exc_info=True)
        elif cfg.msa.enabled:
            logger.info("  MSA chunk plots are written with --extended")
        else:
            logger.info("  MSA visualization disabled in config")
    else:
        logger.info("  No tree: relabeled tree and alignment plots skipped")

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
                organism=organism,
                extended=cfg.extended,
                records=df_final,
            )

            logger.info(f"  Divergence analysis complete: {dirs['divergence_analysis']}")

            # Print summary to console
            divergence_analysis.print_divergence_summary(divergence_results)

            # Species-level divergence analysis (if species data available)
            logger.info("9.2: Calculating species-level divergence...")
            species_summary_csv = dirs['species_analysis'] / 'species_summary.csv'

            if not cfg.extended:
                logger.info("  Species-level divergence summary written only with --extended")
            elif species_summary_csv.exists() and 'divergence_matrix' in divergence_results:
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

                    logger.info(f"  Species-level divergence analysis complete")

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
                        n_species_with_div = len(species_faceted_div_results['matrices'])
                        logger.info(f"  Generated within-species divergence matrices for {n_species_with_div} species")
                    except Exception as e:
                        step_warnings.record("9.3 within-species divergence", e)
                        logger.debug("Within-species divergence error details:", exc_info=True)

                except Exception as e:
                    step_warnings.record("9.2 species divergence", e)
                    logger.debug("Species divergence error details:", exc_info=True)
            else:
                logger.info("  Skipping species-level divergence (species data not available)")

        else:
            step_warnings.record("9 divergence", "skipped: haplotype or taxonomy file not found")

    except Exception as e:
        step_warnings.record("9 divergence", e)
        logger.debug("Divergence analysis error details:", exc_info=True)

    # ========================================================================
    # PHASE 9.5: Metadata Summarization
    # ========================================================================
    # Coverage is always written; --no-metadata-summary only drops the
    # cross-tabs that --extended adds.
    metadata_summary_results = None
    logger.info("")
    logger.info("PHASE 9.5: Metadata Summarization")
    logger.info("-" * 80)

    try:
        metadata_summary_results = metadata_summary.run_metadata_summary(
            annotated_df=df_final,
            output_dir=dirs['base'],
            organism=organism,
            fields=metadata_fields,
            normalize_sex=normalize_sex,
            temporal_analysis=temporal_analysis,
            color_map=None,
            haplotype_col='haplotype_sp',
            extended=cfg.extended,
            cross_tabs=metadata_summary_enabled
        )

        logger.info(f"  Metadata summarization complete")
        logger.info(f"    - Fields tabulated: {metadata_summary_results.get('n_fields_tabulated', 0)}")
        logger.info(f"    - Visualizations created: {metadata_summary_results.get('n_visualizations', 0)}")
        if metadata_summary_results.get('temporal_analysis'):
            logger.info(f"    - Temporal analysis: enabled")
        if metadata_summary_results.get('missing_fields'):
            logger.info(f"    - Missing fields (skipped): {', '.join(metadata_summary_results.get('missing_fields', []))}")

    except Exception as e:
        step_warnings.record("9.5 metadata summary", e)
        logger.debug("Metadata summarization error details:", exc_info=True)

    # ========================================================================
    # PHASE 10: Figures
    # ========================================================================
    # Drawn from the records.csv table and the rooted tree, by taxon, straight
    # into the run folder (run_figures); each figure folder gets a manifest.json.
    logger.info("")
    logger.info("PHASE 10: Figures")
    logger.info("-" * 80)

    records_df = None
    try:
        confidence_cols = [c for c in ('processid', 'basin_confidence', 'distance_to_boundary_km')
                           if c in df_final.columns]
        records_df = layout.build_records_table(
            df_with_basins, status_df, df_qc_all, assign_table, geo_category,
            region_confidence=df_final[confidence_cols] if len(confidence_cols) > 1 else None,
        )
        manifests, figure_warnings = run_figures.write_run_figures(
            records_df, output_dir,
            tree_path=tree_path,
            region_col=geo_category if geo_analysis_performed else None,
            geography=geo_analysis_performed,
            metadata_fields=metadata_fields,
            extended=cfg.extended,
            cross_tabs=metadata_summary_enabled,
            temporal=temporal_analysis,
            normalize_sex=normalize_sex,
            settings=cfg.visualization,
        )
        for folder in manifests:
            logger.info(f"  Figures: {folder}/ (manifest.json)")
        for w in figure_warnings:
            step_warnings.record("10 figures", w)
    except Exception as e:
        step_warnings.record("10 figures", e)
        logger.debug("Figure error details:", exc_info=True)

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
            # No FASTA or tree copies: the kit reads ../haplotypes.fasta and
            # ../phylogeny/tree.nwk (P2)
            plot_export_results = plot_export.export_plots_complete(
                df=mappable_coordinates(df_final, mappable_ids),
                output_dir=output_dir,
                organism=organism,
                color_map=haplotype_colors
            )

            logger.info(f"  Plot regeneration kit exported to {output_dir / 'plots'}")

        except Exception as e:
            step_warnings.record("10.5 plot data export", e)
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

            logger.info(f"  Population genetics formats exported to {output_dir / 'exports'}")

        except Exception as e:
            step_warnings.record("11 popgen export", e)
            logger.debug("Popgen export error details:", exc_info=True)

    # ========================================================================
    # PHASE 12: Reports
    # ========================================================================
    logger.info("")
    logger.info("PHASE 12: Generating Reports")
    logger.info("-" * 80)

    try:
        # Assignment summary and species counts (P0-3, P0-9)
        summary_output = output_dir / "summary.csv"
        species_counts = taxon_labels.species_count_summary(
            input_labels=df_with_basins['species'] if 'species' in df_with_basins.columns else [],
            qc_passed_labels=df_qc_passed['species'] if 'species' in df_qc_passed.columns else [],
            haplotype_taxonomy=assign_table,
            majority_threshold=cfg.taxonomy.majority_species_threshold,
        )
        prov.applied["assignment_summary"] = reports.generate_assignment_summary(
            annotated_csv=str(annotated_csv),
            diagnostics_csv=str(diagnostics_csv),
            output_csv=str(summary_output),
            n_input=provenance.count_data_rows(tsv_path),
            n_haplotypes=len(haplotype_records),
            species_counts=species_counts,
            n_excluded_error_tied=details.n_excluded_error_tied,
            n_excluded_error_far=details.n_excluded_error_far,
            input_counts=prov.applied.get("input_rows"),
        )

        # Verify file was created
        if summary_output.exists():
            logger.info(f"  Generated assignment summary: {summary_output}")
        else:
            step_warnings.record("12 assignment summary", f"{summary_output} was not written")

    except Exception as e:
        step_warnings.record("12 assignment summary", e)

    # ========================================================================
    # PHASE 13: Output layout (P2)
    # ========================================================================
    logger.info("")
    logger.info("PHASE 13: Output layout")
    logger.info("-" * 80)

    records_csv = output_dir / "records.csv"
    moves = {}
    try:
        if records_df is None:                      # Phase 10 builds it; rebuilt if that failed
            confidence_cols = [c for c in ('processid', 'basin_confidence', 'distance_to_boundary_km')
                               if c in df_final.columns]
            records_df = layout.build_records_table(
                df_with_basins, status_df, df_qc_all, assign_table, geo_category,
                region_confidence=df_final[confidence_cols] if len(confidence_cols) > 1 else None,
            )
        records_df.to_csv(records_csv, index=False)
        logger.info(f"  records.csv: {len(records_df)} records")

        tree_tips = None
        if tree_path and Path(tree_path).exists():
            tree_tips = {tip.name for tip in Phylo.read(str(tree_path), "newick").get_terminals()}
        haplotypes_df = layout.build_haplotypes_table(
            haplotype_stats, assign_table, records_df, geo_category, tree_tips=tree_tips,
            evidence=evidence)                 # evidence: step 5.4
        haplotypes_df.to_csv(output_dir / "haplotypes.csv", index=False)
        logger.info(f"  haplotypes.csv: {len(haplotypes_df)} haplotypes")
        if evidence_not_applied is not None and len(evidence_not_applied):
            taxon_evidence.suggestions_table(
                evidence_not_applied, haplotypes_df, {r.id: str(r.seq) for r in haplotype_records}
            ).to_csv(output_dir / taxon_evidence.SUGGESTIONS_FILE, index=False)
            logger.info(f"  {taxon_evidence.SUGGESTIONS_FILE}: {len(evidence_not_applied)} "
                        f"suggestion(s) not applied")

        try:
            variants = label_variants.detect_label_variants(records_df, label_map=label_map_df)
            variants.to_csv(output_dir / "label_variants.csv", index=False)
            logger.info(f"  label_variants.csv: {len(variants)} possible spelling variant(s)")
        except Exception as e:
            step_warnings.record("13 label variants", e)

        try:
            _write_error_filter_diagnostics(records_df, df_with_basins, details, output_dir,
                                            summary_output, prov, cfg)
        except Exception as e:
            step_warnings.record("13 error-filter diagnostics", e)

        if geo_analysis_performed:
            (output_dir / "geography").mkdir(exist_ok=True)
            layout.build_haplotype_by_region(records_df, geo_category).to_csv(
                output_dir / "geography" / "haplotype_by_region.csv", index=False)

        moves = layout.plan_moves(dirs['base'], organism, extended=cfg.extended)
    except Exception as e:
        step_warnings.record("13 output layout", e)
        logger.debug("Output layout error details:", exc_info=True)

    # Place staged files at their v2 paths; per-species records with --extended
    try:
        placed = layout.apply_moves(moves, output_dir)
        logger.info(f"  Placed {len(placed)} files in the run folder")
        if cfg.extended and records_csv.exists():
            layout.write_per_species_records(
                pd.read_csv(records_csv, low_memory=False),
                output_dir / "extended" / "species" / "per_species")
    except Exception as e:
        step_warnings.record("13 output layout", e)

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
    try:
        remove_empty_directories(output_dir)
    except Exception as e:
        logger.debug(f"Directory cleanup encountered minor issues: {e}")

    # HTML report from the final run folder; written before the README, which lists it
    html_report_path = None
    if not no_report:
        logger.info("")
        logger.info("Generating HTML summary report...")
        try:
            html_report_path = run_report.generate_html_report(output_dir, organism, prov.to_dict(),
                                                               __version__)
            logger.info(f"  Generated HTML report: {html_report_path}")
        except Exception as e:
            step_warnings.record("13 HTML report", e)
            logger.debug("HTML report error details:", exc_info=True)
    else:
        logger.info("")
        logger.info("HTML report generation skipped (--no-report)")

    # Run README: every file and column (P2)
    try:
        _, undocumented = run_readme.write_run_readme(output_dir, organism, region_col=geo_category)
        if undocumented:
            logger.debug(f"README has no definition for: {', '.join(undocumented)}")
    except Exception as e:
        step_warnings.record("13 run README", e)

    # ========================================================================
    # Pipeline Complete
    # ========================================================================
    step_warnings.detach()
    prov.finish()
    prov.write(prov_path)
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"  Output directory: {output_dir}")
    logger.info(f"  Run provenance: {prov_path}")
    logger.info(f"  Records: {records_csv}")
    logger.info(f"  Haplotype sequences: {output_dir / 'haplotypes.fasta'}")
    if tree_path:
        logger.info(f"  Phylogenetic tree: {output_dir / 'phylogeny' / 'tree.nwk'}")
    if html_report_path:
        logger.info(f"  HTML summary report: {html_report_path}")
    logger.info("=" * 80)
    step_warnings.log_summary(logger, organism)

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
        epilog=AVAILABLE_COMMANDS_EPILOG + """
Examples:
  # Test default threshold range
  boldgenotyper-sweep data/Sphyrna_lewini.tsv

  # Test custom thresholds
  boldgenotyper-sweep data/samples.tsv --thresholds 0.005,0.01,0.015,0.02

  # Run in parallel with 8 threads
  boldgenotyper-sweep data/samples.tsv --threads 8

  # Invertebrates: use the genetic code of the final run
  boldgenotyper-sweep data/samples.tsv --genetic-code 5

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
        '--genetic-code',
        type=int,
        default=2,
        dest='genetic_code',
        help='NCBI genetic code table for ORF validation (default: 2 = vertebrate mitochondrial). '
             'Use the code of the final run (5 for most invertebrates).'
    )

    _add_input_row_arguments(parser)

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
            keep_intermediates=args.keep_intermediates,
            genetic_code=args.genetic_code,
            marker=args.marker,
            allow_malformed_rows=args.allow_malformed_rows,
            label_map=args.label_map,
        )

        print("\n✓ Parameter sweep complete!")
        recommended = results.get('recommended_threshold', 'N/A')
        if isinstance(recommended, (int, float)):
            print(f"  Recommended threshold: {recommended:.3f}")
        elif recommended is None:
            print(f"  No elbow: sweep at least {parameter_sweep.MIN_THRESHOLDS_FOR_ELBOW} "
                  "thresholds to locate one")
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


def _iso_date(value: str) -> str:
    """argparse type: a YYYY-MM-DD date, returned as text."""
    try:
        return datetime.strptime(value, "%Y-%m-%d").date().isoformat()
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected YYYY-MM-DD, got {value!r}")


def _build_main_parser() -> argparse.ArgumentParser:
    """Build the main boldgenotyper argument parser."""
    parser = argparse.ArgumentParser(
        description='BOLDGenotyper: haplotype-first COI genotyping with optional phylogeny, divergence, and mapping outputs.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=AVAILABLE_COMMANDS_EPILOG + """
Examples:
  # Basic usage (organism inferred from TSV filename)
  boldgenotyper data/Euprymna_scolopes.tsv

  # Specify organism name and output directory
  boldgenotyper data/Carcharhinus.tsv --organism Carcharhinus --output results/Carcharhinus_analysis

  # Use the sweep-recommended singleton threshold and record the download date
  boldgenotyper data/Carcharhinus.tsv --singleton-distance 0.015 --download-date 2026-09-30

  # Enable phylogeny and root on an in-tree outgroup label or taxon
  boldgenotyper data/Euprymna.tsv --build-tree --phylo-outgroup-label "Eusphyra blochii H0010"
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
  - Trees are built from haplotypes (phylogeny/tree.nwk; tree_labeled.nwk has species labels).
  - Geographic analysis assigns samples to regions via point-in-polygon spatial joins.
  - Default: GOaS ocean basins (marine). Use --custom-shp for freshwater basins,
    terrestrial ecoregions, watersheds, biomes, or any polygon shapefile.
  - The run folder's README.md describes every file and column. --extended adds
    extended/, --export-format adds exports/, --plot-kit adds plots/.
"""
    )

    # Required arguments
    parser.add_argument(
        'tsv',
        type=Path,
        nargs='?',
        default=None,
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
        '--tree-min-members',
        type=int,
        default=None,
        help='Minimum members for a haplotype to be a tree tip (default: 1, every '
             'non-suspect haplotype). Suspect haplotypes are always excluded. '
             'Excluded haplotypes are listed with their reason in excluded_tips.txt.'
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
        '--plot-kit',
        action='store_true',
        help='Write plots/: plot data, a plot configuration and scripts that '
             'regenerate the map, bar charts and tree'
    )

    parser.add_argument(
        '--export-format',
        nargs='+',
        default=None,
        choices=['arlequin', 'popart', 'dnasp', 'generic', 'all'],
        help='Write exports/ in these population genetics formats (default: none). '
             'Options: arlequin (Arlequin .arp), popart (PopART/NEXUS), '
             'dnasp (DnaSP .fas), generic (CSV), all. '
             'Populations in Arlequin and DnaSP are defined by --geo-category. '
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
        default=None,
        dest='shapefile_field',
        help='Name of shapefile attribute containing region labels. '
             'If not specified, the first available non-geometry field is auto-detected '
             'and a warning is logged. '
             'Examples: "ECO_NAME" for Ecoregions2017, "HYBAS_ID" for HydroBASINS, '
             '"name" for GOaS. '
             'Use this to specify which field contains the region identifiers.'
    )

    parser.add_argument(
        '--geo-category',
        type=str,
        default=None,
        dest='geo_category',
        help='Name for geographic category in outputs. Defaults to "ocean_basin" '
             'when using GOaS (no --custom-shp), or "geographic_region" when using '
             'a custom shapefile. Examples: "freshwater_basin", "biogeographic_realm", '
             '"watershed". This label is used in output files, plots, and column '
             'names; avoid names of BOLD fields (e.g. "realm", "biome", "ecoregion"), '
             'which the region column would replace.'
    )

    parser.add_argument(
        '--shp',
        type=str,
        default=None,
        dest='shp_name',
        help='Shorthand name of a shapefile registered in your shapefiles.yaml '
             '(see --init-shapefile-config). Resolves its path, field, and '
             'geo-category. Explicit --shp-field/--geo-category override the '
             'registered values. Cannot be combined with --custom-shp.'
    )

    parser.add_argument(
        '--shapefile-config',
        type=Path,
        default=None,
        dest='shapefile_config',
        help='Path to a shapefile registry file. Overrides the '
             'BOLDGENOTYPER_SHAPEFILE_CONFIG env var and the default '
             '(~/.config/boldgenotyper/shapefiles.yaml).'
    )

    parser.add_argument(
        '--init-shapefile-config',
        action='store_true',
        dest='init_shapefile_config',
        help='Write a commented shapefiles.yaml template to the target config '
             'path (--shapefile-config if given, else the default location) and exit.'
    )

    # Metadata Summarization Arguments (enabled by default).
    parser.add_argument(
        '--no-metadata-summary',
        dest='no_metadata_summary',
        action='store_true',
        help='With --extended, skip the metadata cross-tabs (metadata coverage is always '
             'written).'
    )

    parser.add_argument(
        '--metadata-fields',
        nargs='+',
        default=None,
        help='Metadata fields to tabulate. Default: sex, life_stage, reproduction, habitat, '
             'country/ocean, province/state. Other BOLD fields can be named, e.g. realm, '
             'biome, ecoregion (BOLD derives these from the coordinates using terrestrial '
             'ecoregions), geoid, country_iso'
    )

    parser.add_argument(
        '--no-normalize-sex',
        action='store_true',
        help='Keep sex values as given. By default sex codes are normalized '
             '(M/m/male -> Male, F/f/female -> Female, U -> missing); blank values are always missing.'
    )

    parser.add_argument(
        '--no-temporal-analysis',
        action='store_true',
        help='With --extended, skip the temporal and first-sampled summaries of '
             'collection dates.'
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
             'Default: 0.005 (0.5%%). Use boldgenotyper-sweep to choose a value for '
             'your dataset (recommendations.txt reports the elbow and how stable '
             'assignments are between thresholds), and always pass it explicitly.'
    )

    parser.add_argument(
        '--no-assign-corrected',
        action='store_true',
        help='Leave records of error-filtered singletons unassigned (status '
             'excluded_error_singleton) instead of assigning them to their nearest '
             'retained haplotype (status assigned_corrected).'
    )

    parser.add_argument(
        '--keep-intermediates',
        action='store_true',
        help='Keep intermediate files (default: remove after pipeline completion)'
    )

    parser.add_argument(
        '--genetic-code',
        type=int,
        default=2,
        dest='genetic_code',
        help='NCBI genetic code table for ORF validation (default: 2 = vertebrate mitochondrial). '
             'Use 5 for most invertebrates (invertebrate mitochondrial). '
             'See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for full list.'
    )

    parser.add_argument(
        '--extended',
        action='store_true',
        help='Also write extended/: metadata cross-tabs (normalized vocabularies), '
             'temporal and first-sampled summaries, the label-based barcoding gap, '
             'species tables and diversity, species-level divergence, alignment (MSA) '
             'plots, and region bar charts and faceted maps. Several depend on sampling '
             'and labels; see the run README.'
    )

    parser.add_argument(
        '--download-date',
        type=_iso_date,
        default=None,
        help='Date the BOLD TSV was downloaded (YYYY-MM-DD); recorded in run_provenance.json.'
    )

    parser.add_argument(
        '--taxon-curation',
        type=Path,
        default=None,
        metavar='FILE',
        help='CSV/TSV of accepted haplotype taxa (seq_md5 or sequence, taxon, reason; '
             'optional curator, date), applied after the label rule. Each row is recorded '
             'in run_provenance.json; see the run README.'
    )

    _add_input_row_arguments(parser)

    parser.add_argument(
        '--no-evidence-taxa',
        action='store_true',
        help='Report taxon evidence without applying it: a genus-level haplotype keeps the '
             'genus even when two evidence lines name its leading species label (see the '
             'run README, Taxonomy).'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'BOLDGenotyper {__version__}'
    )

    return parser


def main():
    """Main CLI entry point."""
    msg = removed_flag_error(sys.argv[1:], "boldgenotyper")
    if msg:
        print(f"boldgenotyper: error: {msg}", file=sys.stderr)
        sys.exit(2)
    parser = _build_main_parser()
    args = parser.parse_args()

    # Handle --init-shapefile-config: write template and exit (no TSV needed).
    if args.init_shapefile_config:
        target = args.shapefile_config or shapefile_registry.default_config_path()
        try:
            written = shapefile_registry.write_template(target)
        except FileExistsError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        print(f"Wrote shapefile config template to: {written}")
        print("Edit it to register your shapefiles, then use: --shp <name>")
        return 0

    # A TSV is required for every run except --init-shapefile-config.
    if args.tsv is None:
        print("Error: the following argument is required: tsv (input BOLD TSV file)",
              file=sys.stderr)
        return 1

    # Resolve a registered shapefile shorthand (--shp), if given.
    goas_path_override = None
    try:
        selection = shapefile_registry.resolve_selection(
            shp_name=args.shp_name,
            custom_shp=args.shapefile_path,
            explicit_field=args.shapefile_field,
            explicit_geo_category=args.geo_category,
            config_path=args.shapefile_config,
        )
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if selection is not None:
        args.geo_category = selection.geo_category
        if selection.is_goas:
            # Route through the existing GOaS ocean-basin logic; override the
            # GOaS shapefile path so the registered copy is used.
            args.shapefile_path = None
            goas_path_override = selection.path
        else:
            args.shapefile_path = selection.path
            args.shapefile_field = selection.field

    # Validate input file
    if not args.tsv.exists():
        print(f"Error: Input TSV file not found: {args.tsv}", file=sys.stderr)
        return 1

    # Validate custom shapefile if provided
    if args.shapefile_path and not args.shapefile_path.exists():
        print(f"Error: Custom shapefile not found: {args.shapefile_path}", file=sys.stderr)
        return 1

    # Pre-validate shapefile field before running pipeline
    if args.shapefile_path:
        validate_shapefile_field(args.shapefile_path, args.shapefile_field)

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
    logger.info(f"BOLDGenotyper {__version__} run started "
                f"{datetime.now().astimezone().isoformat(timespec='seconds')}")
    logger.info(f"Command: {shlex.join(sys.argv)}")

    # Load and configure pipeline
    cfg = config.get_default_config()
    cfg = cfg.update(
        haplotype__min_singleton_distance=args.singleton_distance,
        haplotype__assign_corrected=not args.no_assign_corrected,
        extended=args.extended,
        n_threads=args.threads,
        output_dir=output_dir,
        log_level=args.log_level,
        phylogenetic__build_tree=args.build_tree,
        taxonomy__apply_evidence=not args.no_evidence_taxa,
        **({'phylogenetic__min_cluster_size': args.tree_min_members}
           if args.tree_min_members is not None else {}),
        phylogenetic__outgroup_fasta=args.phylo_outgroup_fasta,
        phylogenetic__outgroup_label=args.phylo_outgroup_label,
        phylogenetic__outgroup_taxon=args.phylo_outgroup_taxon,
        keep_intermediates=args.keep_intermediates,
        coi_validation__mitochondrial_code=args.genetic_code
    )

    if goas_path_override is not None:
        cfg = cfg.update(geographic__goas_shapefile_path=goas_path_override)

    # Print banner
    print("=" * 80)
    print("BOLDGenotyper Pipeline")
    print("=" * 80)
    print(f"Organism: {organism}")
    print(f"Input: {args.tsv}")
    print(f"Output: {output_dir}")
    print()
    print("Parameters:")
    print(f"  Singleton distance: {args.singleton_distance} ({args.singleton_distance*100:.1f}% — use boldgenotyper-sweep to optimise)")
    print(f"  Threads: {args.threads}")
    print(f"  Genetic code: {args.genetic_code}")
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
            export_plot_data=args.plot_kit,
            export_popgen_formats=args.export_format,
            metadata_summary_enabled=not args.no_metadata_summary,
            metadata_fields=args.metadata_fields,
            normalize_sex=not args.no_normalize_sex,
            temporal_analysis=not args.no_temporal_analysis,
            shapefile_path=args.shapefile_path,
            shapefile_field=args.shapefile_field,
            geo_category=args.geo_category,
            argv=sys.argv,
            download_date=args.download_date,
            taxon_curation=args.taxon_curation,
            marker=args.marker,
            allow_malformed_rows=args.allow_malformed_rows,
            label_map=args.label_map,
        )

        return 0 if success else 1

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}", exc_info=True)
        print(f"\nError: Pipeline failed. Check log file: {log_file}", file=sys.stderr)
        return 1


def main_labels(argv=None) -> int:
    """
    CLI entry point for 'boldgenotyper-labels'.

    boldgenotyper-labels suggest RUN_DIR [--backbone gbif|worms] [--output FILE]
    checks a run's species labels against GBIF or WoRMS and writes a suggested
    label map (spelling fixes and synonyms) to review before passing it to
    --label-map. Needs network access; nothing in the run changes.
    """
    from . import label_backbone

    parser = argparse.ArgumentParser(
        prog='boldgenotyper-labels',
        description='Suggest a label map by checking a run\'s species labels against a '
                    'taxonomic backbone. Review the file, then pass it to --label-map.',
    )
    sub = parser.add_subparsers(dest='command', required=True)
    suggest = sub.add_parser('suggest', help='Write a suggested label map for a run folder')
    suggest.add_argument('run_dir', type=Path, help='A v2 run folder (with records.csv)')
    suggest.add_argument('--backbone', choices=sorted(label_backbone.BACKBONES), default='gbif',
                         help='Taxonomic backbone (default: gbif; worms for marine taxa)')
    suggest.add_argument('--output', type=Path, default=None,
                         help='Output CSV (default: RUN_DIR/label_map_suggested.csv)')
    args = parser.parse_args(argv)

    output = args.output or args.run_dir / "label_map_suggested.csv"
    try:
        suggestions = label_backbone.suggest_label_map(args.run_dir, backbone=args.backbone)
    except (label_backbone.BackboneError, FileNotFoundError) as e:
        print(f"boldgenotyper-labels: error: {e}", file=sys.stderr)
        return 1
    label_backbone.write_suggestions(suggestions, output)
    print(f"{len(suggestions)} suggestion(s) written to {output}. Review them; the confirmed "
          f"rows are a label map for --label-map.")
    return 0


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
            --haplotypes analysis/ \
            --output query_results/

        # Query with metadata enrichment
        boldgenotyper-query \
            --query new_samples.fasta \
            --haplotypes analysis/ \
            --analysis-dir analysis/ \
            --output query_results/
    """
    parser = argparse.ArgumentParser(
        prog="boldgenotyper-query",
        description="Query new COI sequences against previously identified haplotypes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=AVAILABLE_COMMANDS_EPILOG + """
Examples:
  # Query single sequence
  boldgenotyper-query --query new_sample.fasta \\
                       --haplotypes analysis/

  # Query multiple sequences with metadata
  boldgenotyper-query --query samples.fasta \\
                       --haplotypes analysis/ \\
                       --analysis-dir analysis/ \\
                       --output results/

  # Adjust top matches and length filters
  boldgenotyper-query --query samples.fasta \\
                       --haplotypes haplotypes.fasta \\
                       --top-n 20 \\
                       --min-length 150 \\
                       --max-length 1500

Output Files:
  - query_calls.csv: one row per query (call, equivalence set, species call and
    its basis, pooled member labels, flags): the result to report
  - query_candidates.csv: the top --top-n haplotypes per query with their numbers
  - query_alignments.txt: each query's call and its best alignment
  - query_provenance.json: reference run, parameters, curation file, outputs
  - query_results.json: everything as one JSON document (--json only)

How queries are called:
  Each query is aligned over its whole overlap with every haplotype. A match is
  eligible when it covers --min-query-coverage of the shorter sequence and
  reaches --assign-identity. The equivalence set is every eligible haplotype
  within --max-mismatch-delta differences (mismatches + gap columns) of the
  best. species_call is the leading species of the set's pooled member labels
  when species_support >= --min-species-support. No identity value is
  interpreted as "same species" or "same genus"; the numbers are reported.

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
        help='A BOLDGenotyper run folder (v2: haplotypes.fasta at its root; v1: '
             'haplotypes/*_haplotypes.fasta) or a haplotype FASTA file. Species calls need a v2 run folder.'
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
        '--assign-identity',
        type=float,
        default=0.97,
        help='Minimum identity (fraction, 0-1) to call a match assignable '
             '(default: 0.97).'
    )

    parser.add_argument(
        '--min-query-coverage',
        type=float,
        default=0.95,
        help='Minimum fraction of the shorter sequence (query or haplotype) that the '
             'alignment must cover for a match to be eligible (default: 0.95).'
    )

    parser.add_argument(
        '--max-mismatch-delta',
        type=int,
        default=0,
        help='Haplotypes within this many differences (mismatches + gap columns) of '
             'the best eligible match form the equivalence set (default: 0).'
    )

    parser.add_argument(
        '--min-species-support',
        type=float,
        default=0.70,
        help='Minimum share of pooled member labels for a species call (default: 0.70).'
    )

    parser.add_argument(
        '--taxon-curation',
        type=Path,
        default=None,
        metavar='FILE',
        help='CSV/TSV of accepted haplotype taxa (seq_md5 or sequence, taxon, reason; '
             'optional curator, date) applied to the reference run\'s haplotype taxa for '
             'this query. Needs a v2 run folder as --haplotypes. Curated taxa already in '
             'the run folder are used without it.'
    )

    parser.add_argument(
        '--query-metadata',
        type=Path,
        default=None,
        metavar='FILE',
        help='CSV/TSV with query_id and any of lat/lon or coord, collection_date, '
             'country/ocean, region. Adds the queries to the figures; a query\'s region is '
             'its region value, else its coordinates placed with the reference run\'s shapefile.'
    )

    parser.add_argument(
        '--json',
        action='store_true',
        help='Also write query_results.json (calls, candidates, reference run and '
             'parameters as one JSON document).'
    )

    parser.add_argument(
        '--search-records',
        action='store_true',
        help='Also search every QC-passed reference record (records_qc_passed.fasta in a '
             'v2 run folder) and report exact-match records with their labels, '
             'haplotype and status, including records not assigned to any haplotype.'
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
    msg = removed_flag_error(argv if argv is not None else sys.argv[1:], "boldgenotyper-query")
    if msg:
        print(f"boldgenotyper-query: error: {msg}", file=sys.stderr)
        return 2
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
        # Import query modules
        from . import haplotype_query, query_assignment

        # Resolve the reference: --haplotypes may be a FASTA file or a run directory
        try:
            reference_fasta, analysis_dir = query_assignment.resolve_reference(
                args.haplotypes, args.analysis_dir)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

        reference_run = query_assignment.load_reference_run(analysis_dir)
        if reference_run is None:
            logger.warning("Reference run not identified (no run_provenance.json); "
                           "outputs record the haplotype file only")
        if not args.query.exists():
            raise FileNotFoundError(f"Query file not found: {args.query}")

        # Haplotype taxa of the reference run, with query-time curation
        haplotype_taxa = query_assignment.load_haplotype_taxa(analysis_dir)
        curation_parameters, curation_provenance = {}, None
        if args.taxon_curation:
            if haplotype_taxa is None:
                print("Error: --taxon-curation needs a v2 run folder as --haplotypes "
                      "(haplotypes.csv and run_provenance.json)", file=sys.stderr)
                return 1
            try:
                curation = taxon_curation_mod.load_curation(args.taxon_curation)
            except taxon_curation_mod.CurationError as e:
                print(f"Error: taxon curation file: {e}", file=sys.stderr)
                return 1
            ref_seqs = {r.id: str(r.seq) for r in SeqIO.parse(reference_fasta, "fasta")}
            hap_meta = pd.DataFrame({"haplotype_id": list(ref_seqs),
                                     "seq_md5": [haplotype_ids.seq_md5(s) for s in ref_seqs.values()]})
            applied_rows, not_applied_rows = taxon_curation_mod.match_curation(
                curation, hap_meta, ref_seqs)
            for r in applied_rows.to_dict("records"):
                haplotype_taxa[r["haplotype_id"]] = query_assignment.HaplotypeTaxon(
                    r["taxon"], taxon_curation_mod.taxon_level(r["taxon"]), "curated")
            logger.info(f"Curated taxa applied to {len(applied_rows)} haplotype(s)")
            for r in not_applied_rows.to_dict("records"):
                logger.warning(f"Curation row {r['row']} ({r['taxon']}) not applied: "
                               f"{r['reason_not_applied']}")
            curation_parameters = {"taxon_curation": str(args.taxon_curation)}
            curation_provenance = {
                "path": str(args.taxon_curation),
                "md5": provenance.file_md5(args.taxon_curation),
                "n_applied": len(applied_rows),
                "n_not_applied": len(not_applied_rows),
                "not_applied": not_applied_rows[["row", "taxon", "reason_not_applied"]]
                .to_dict("records"),
            }

        query_records = haplotype_query.load_query_sequences(
            args.query, min_length=args.min_length, max_length=args.max_length)
        if not query_records:
            raise ValueError("No valid query sequences after length filtering")
        query_meta = None
        if args.query_metadata:
            try:
                query_meta = query_overlays.load_query_metadata(
                    args.query_metadata, [q.id for q in query_records])
            except query_overlays.QueryMetadataError as e:
                print(f"Error: query metadata: {e}", file=sys.stderr)
                return 1

        # Run query
        logger.info("Starting haplotype query analysis...")
        results, metadata, verdicts = haplotype_query.query_against_haplotypes(
            query_fasta=args.query,
            haplotype_fasta=reference_fasta,
            analysis_dir=analysis_dir,
            top_n=args.top_n,
            min_length=args.min_length,
            max_length=args.max_length,
            min_identity=args.assign_identity,
            min_query_coverage=args.min_query_coverage,
            max_mismatch_delta=args.max_mismatch_delta,
            min_species_support=args.min_species_support,
            query_records=query_records,
            haplotype_taxa=haplotype_taxa,
        )

        # Exact-match reference records (Q10)
        record_matches = None
        if args.search_records:
            records_fasta = Path(analysis_dir) / "records_qc_passed.fasta" if analysis_dir else None
            if records_fasta is None or not records_fasta.exists():
                logger.warning("--search-records needs a v2 run folder with records_qc_passed.fasta; skipped")
            else:
                record_table, _ = query_assignment.load_record_table(analysis_dir)
                matches_df = haplotype_query.search_records(query_records, records_fasta, record_table)
                matches_path = output_dir / "query_record_matches.csv"
                matches_df.to_csv(matches_path, index=False)
                record_matches = {q.id: int((matches_df["query_id"] == q.id).sum()) for q in query_records}
                logger.info(f"Exact-match records written to {matches_path}")
        parameters = {
            "assign_identity": args.assign_identity,
            "min_query_coverage": args.min_query_coverage,
            "max_mismatch_delta": args.max_mismatch_delta,
            "min_species_support": args.min_species_support,
            "top_n": args.top_n,
            "search_records": args.search_records,
            **curation_parameters,
        }

        # Write results (query_outputs); provenance last, so it lists every file
        logger.info("Writing results...")
        n_reference = (verdicts[0].n_reference_haplotypes if verdicts
                       else sum(1 for _ in SeqIO.parse(str(reference_fasta), "fasta")))
        qprov = query_outputs.build_query_provenance(
            argv if argv is not None else sys.argv, args.query, len(query_records),
            reference_fasta, n_reference, reference_run, parameters, curation_provenance)
        calls = query_outputs.calls_table(verdicts, reference_run, record_matches)
        query_outputs.write_calls(calls, output_dir)
        candidates = query_outputs.candidates_table(
            results, verdicts, haplotypes=metadata, haplotype_taxa=haplotype_taxa,
            min_identity=args.assign_identity, min_query_coverage=args.min_query_coverage)
        query_outputs.write_candidates(candidates, output_dir)
        query_outputs.write_alignments(results, verdicts, output_dir, reference_run=reference_run,
                                       haplotypes=metadata, haplotype_file=reference_fasta)
        if args.json:
            query_outputs.write_json(results, verdicts, output_dir, parameters=parameters,
                                     reference_run=reference_run, haplotypes=metadata)

        # Place the queries in the reference run (v2 run folders)
        reference = query_overlays.load_reference(analysis_dir)
        ref_curation = (reference.provenance.get("applied", {}).get("taxon_curation")
                        if reference is not None else None)
        if ref_curation:
            qprov["reference_taxon_curation"] = {
                k: ref_curation.get(k) for k in ("path", "md5", "n_applied", "n_not_applied")}
        ref_evidence = (reference.provenance.get("applied", {}).get("taxon_evidence")
                        if reference is not None else None)
        if ref_evidence:
            qprov["reference_taxon_evidence"] = {k: ref_evidence.get(k) for k in ("apply", "n_applied")}
        if query_meta is not None:
            unmatched = query_meta.attrs.get("unmatched_ids", [])
            n_rows = query_meta.attrs.get("n_rows", len(query_meta))
            query_meta = query_overlays.assign_query_regions(query_meta, reference)
            qprov["query_metadata"] = {
                "path": str(args.query_metadata),
                "md5": provenance.file_md5(args.query_metadata),
                "n_rows": n_rows,
                "n_matched": n_rows - len(unmatched),
                "unmatched_ids": unmatched,
            }
        figures = {}
        if reference is not None:
            figures, figure_warnings = query_overlays.write_query_figures(
                calls, query_meta, reference, output_dir / "query_figures")
            qprov["figures"] = {name: [str(p.relative_to(output_dir)) for p in paths]
                                for name, paths in figures.items()}
            qprov["warnings"] += figure_warnings
        best_alignments = {}
        for r in results:                       # results are in rank order per query
            best_alignments.setdefault(r.query_id, haplotype_query.format_alignment_display(r, width=100))
        query_report.generate_query_report(output_dir, calls, candidates, qprov, figures,
                                           best_alignments)
        query_outputs.write_query_provenance(qprov, output_dir)
        query_outputs.write_query_readme(output_dir, qprov)
        query_outputs.write_query_provenance(qprov, output_dir)      # now lists README.md too

        logger.info("")
        logger.info("=" * 80)
        logger.info("Query Analysis Complete!")
        logger.info("=" * 80)
        logger.info(f"Results written to: {output_dir}")
        logger.info(f"  Calls to report: {output_dir / 'query_calls.csv'}")
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
