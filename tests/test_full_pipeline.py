#!/usr/bin/env python3
"""
Complete BOLDGenotyper Pipeline with Organized Outputs

Runs phylogenetic analysis and visualization generation on pre-processed BOLD data,
organizing outputs into clean subdirectories.

This script expects that data has already been processed through steps 1-4:
  1. TSV parsing → FASTA + annotated CSV
  2. Dereplication → consensus sequences
  3. Genotype assignment → diagnostics
  4. Taxonomy assignment → species annotations

Input Requirements:
  - {organism}_consensus.fasta
  - {organism}_annotated.csv
  - {organism}_diagnostics.csv

Usage:
    python tests/test_full_pipeline.py -i input_dir -o output_dir -n organism_name [options]

Example:
    python tests/test_full_pipeline.py \\
        -i tests/output \\
        -o results/Sphyrnidae_organized \\
        -n Sphyrnidae_test \\
        --threads 8

Output Structure:
    output_dir/
    ├── sequences/
    │   ├── {organism}_consensus.fasta
    │   ├── {organism}_consensus_aligned.fasta
    │   └── {organism}_consensus_trimmed.fasta
    ├── trees/
    │   ├── {organism}_tree.nwk
    │   └── {organism}_tree_relabeled.nwk
    ├── reports/
    │   ├── {organism}_annotated.csv
    │   ├── {organism}_diagnostics.csv
    │   ├── {organism}_taxonomy_conflicts.csv
    │   ├── {organism}_assignment_summary.csv
    │   ├── {organism}_consensus_characterization.csv
    │   └── {organism}_sequence_quality.csv
    ├── visualizations/
    │   ├── maps/
    │   │   ├── {organism}_distribution_map.{png,pdf}
    │   │   └── {organism}_distribution_map_faceted.{png,pdf}
    │   └── plots/
    │       ├── {organism}_tree.{png,pdf}
    │       ├── {organism}_ocean_basin_abundance.{png,pdf}
    │       ├── {organism}_ocean_basin_abundance_faceted.{png,pdf}
    │       ├── {organism}_identity_distribution.{png,pdf}
    │       ├── {organism}_identity_by_status.{png,pdf}
    │       └── {organism}_assignment_status.{png,pdf}
    └── tables/
        ├── {organism}_genotype_color_map.csv
        ├── {organism}_genotype_by_basin.csv
        ├── {organism}_species_by_consensus.csv
        └── {organism}_basin_summary.csv
"""

import argparse
import sys
import shutil
from pathlib import Path
import subprocess
import logging
import re
import pandas as pd
from Bio import Phylo

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from boldgenotyper import phylogenetics, visualization, reports

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("full_pipeline")


def check_dependencies():
    """Check for required external tools."""
    required = ['mafft', 'trimal']
    optional = ['iqtree2', 'fasttree']

    missing = []
    for tool in required:
        if not shutil.which(tool):
            missing.append(tool)

    if missing:
        logger.error(f"Missing required tools: {', '.join(missing)}")
        logger.error("Install with: conda install -c bioconda mafft trimal")
        return False

    has_tree_builder = shutil.which('iqtree2') or shutil.which('fasttree')
    if not has_tree_builder:
        logger.warning("No tree builder found (iqtree2 or fasttree); skipping phylogenetics")

    return True


def organize_outputs(input_dir: Path, output_dir: Path, organism: str):
    """
    Copy and organize input files into output directory structure.

    Returns paths dict with all organized file locations.
    """
    # Create directory structure
    dirs = {
        'sequences': output_dir / 'sequences',
        'trees': output_dir / 'trees',
        'reports': output_dir / 'reports',
        'viz_maps': output_dir / 'visualizations' / 'maps',
        'viz_plots': output_dir / 'visualizations' / 'plots',
        'tables': output_dir / 'tables',
    }

    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # Required input files
    consensus_fasta = input_dir / f"{organism}_consensus.fasta"
    annotated_csv = input_dir / f"{organism}_annotated.csv"
    diagnostics_csv = input_dir / f"{organism}_diagnostics.csv"

    # Check inputs exist
    for f in [consensus_fasta, annotated_csv, diagnostics_csv]:
        if not f.exists():
            raise FileNotFoundError(f"Required input file not found: {f}")

    # Copy input files to organized locations
    shutil.copy(consensus_fasta, dirs['sequences'] / consensus_fasta.name)
    shutil.copy(annotated_csv, dirs['reports'] / annotated_csv.name)
    shutil.copy(diagnostics_csv, dirs['reports'] / diagnostics_csv.name)

    logger.info(f"Organized input files into {output_dir}")

    return dirs


def run_phylogenetics(dirs, organism: str, threads: int):
    """Run phylogenetic analysis pipeline."""
    logger.info("=" * 60)
    logger.info("PHYLOGENETIC ANALYSIS")
    logger.info("=" * 60)

    consensus_fa = dirs['sequences'] / f"{organism}_consensus.fasta"
    aligned_fa = dirs['sequences'] / f"{organism}_consensus_aligned.fasta"
    trimmed_fa = dirs['sequences'] / f"{organism}_consensus_trimmed.fasta"
    tree_nwk = dirs['trees'] / f"{organism}_tree.nwk"

    # MAFFT alignment
    logger.info("Running MAFFT alignment...")
    with subprocess.Popen(
        ["mafft", "--auto", str(consensus_fa)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True
    ) as p:
        aligned_fa.write_text(p.stdout.read())
    logger.info(f"✓ Alignment: {aligned_fa.name}")

    # trimAl
    logger.info("Running trimAl...")
    subprocess.run(
        ["trimal", "-automated1", "-in", str(aligned_fa), "-out", str(trimmed_fa)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    logger.info(f"✓ Trimmed: {trimmed_fa.name}")

    # Tree building
    if shutil.which("iqtree2"):
        logger.info("Building tree with IQ-TREE2...")
        subprocess.run(
            ["iqtree2", "-s", str(trimmed_fa), "-m", "GTR+G", "-B", "1000",
             "-T", str(threads), "--quiet", "-redo"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        treefile = trimmed_fa.with_suffix(".treefile")
        tree_nwk.write_text(treefile.read_text())
    elif shutil.which("fasttree"):
        logger.info("Building tree with FastTree...")
        with subprocess.Popen(
            ["fasttree", "-nt", "-gtr", "-gamma", str(trimmed_fa)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True
        ) as p:
            tree_nwk.write_text(p.stdout.read())
    else:
        logger.warning("No tree builder available; skipping")
        return None

    logger.info(f"✓ Tree: {tree_nwk.name}")
    return tree_nwk


def relabel_tree(tree_nwk: Path, dirs, organism: str, df: pd.DataFrame):
    """Relabel tree tips with safe pretty names."""
    tree_relabeled = dirs['trees'] / f"{organism}_tree_relabeled.nwk"

    # Build relabeling maps
    raw_to_pretty = {}
    if {'consensus_group', 'consensus_group_sp'}.issubset(df.columns):
        tmp = df[['consensus_group', 'consensus_group_sp']].dropna().drop_duplicates()
        raw_to_pretty = dict(zip(tmp['consensus_group'], tmp['consensus_group_sp']))

    raw_to_safe = {k: re.sub(r"[^A-Za-z0-9._-]+", "_", v) for k, v in raw_to_pretty.items()}

    try:
        tree_obj = Phylo.read(str(tree_nwk), "newick")
        for clade in tree_obj.get_terminals():
            if clade.name in raw_to_safe:
                clade.name = raw_to_safe[clade.name]
        Phylo.write(tree_obj, str(tree_relabeled), "newick")
        logger.info(f"✓ Relabeled tree: {tree_relabeled.name}")
        return tree_relabeled, raw_to_safe
    except Exception as e:
        logger.warning(f"Tree relabeling failed: {e}")
        return tree_nwk, raw_to_safe


def generate_reports(dirs, organism: str):
    """Generate Phase 1 and Phase 2 reports."""
    logger.info("=" * 60)
    logger.info("GENERATING QC REPORTS")
    logger.info("=" * 60)

    annotated_csv = dirs['reports'] / f"{organism}_annotated.csv"
    diagnostics_csv = dirs['reports'] / f"{organism}_diagnostics.csv"

    # Phase 1 reports
    conflicts_csv = dirs['reports'] / f"{organism}_taxonomy_conflicts.csv"
    reports.generate_taxonomy_conflicts_report(
        str(annotated_csv), str(diagnostics_csv), str(conflicts_csv))
    logger.info(f"✓ {conflicts_csv.name}")

    summary_csv = dirs['reports'] / f"{organism}_assignment_summary.csv"
    reports.generate_assignment_summary(
        str(annotated_csv), str(diagnostics_csv), str(summary_csv))
    logger.info(f"✓ {summary_csv.name}")

    consensus_char = dirs['reports'] / f"{organism}_consensus_characterization.csv"
    reports.generate_consensus_characterization(
        str(annotated_csv), str(diagnostics_csv), str(consensus_char))
    logger.info(f"✓ {consensus_char.name}")

    seqq_csv = dirs['reports'] / f"{organism}_sequence_quality.csv"
    reports.generate_sequence_quality_metrics(str(annotated_csv), str(seqq_csv))
    logger.info(f"✓ {seqq_csv.name}")


def generate_visualizations(dirs, organism: str, tree_relabeled, raw_to_safe: dict):
    """Generate all visualizations."""
    logger.info("=" * 60)
    logger.info("GENERATING VISUALIZATIONS")
    logger.info("=" * 60)

    annotated_csv = dirs['reports'] / f"{organism}_annotated.csv"
    diagnostics_csv = dirs['reports'] / f"{organism}_diagnostics.csv"
    df = pd.read_csv(annotated_csv)

    # Genotype colors
    if 'consensus_group_sp' in df.columns:
        genotypes = sorted(df['consensus_group_sp'].dropna().unique())
        colors = visualization.get_genotype_colors(len(genotypes))
        color_map = dict(zip(genotypes, colors))

        colormap_csv = dirs['tables'] / f"{organism}_genotype_color_map.csv"
        pd.DataFrame({'genotype': list(color_map.keys()), 'color': list(color_map.values())}).to_csv(colormap_csv, index=False)
        logger.info(f"✓ {colormap_csv.name}")

    # Distribution maps
    if {'lat', 'lon', 'consensus_group_sp'}.issubset(df.columns):
        logger.info("Generating distribution maps...")
        for fmt in ['png', 'pdf']:
            visualization.plot_distribution_map(
                df, str(dirs['viz_maps'] / f"{organism}_distribution_map.{fmt}"),
                'consensus_group_sp', 'lat', 'lon', (10, 6), 300)
        logger.info("✓ Distribution maps")

        if 'assigned_sp' in df.columns:
            for fmt in ['png', 'pdf']:
                visualization.plot_distribution_map_faceted(
                    df, str(dirs['viz_maps'] / f"{organism}_distribution_map_faceted.{fmt}"),
                    'consensus_group_sp', 'assigned_sp', 'lat', 'lon', 10, 5, 300)
            logger.info("✓ Faceted distribution maps")

    # Basin plots
    if {'ocean_basin', 'consensus_group_sp'}.issubset(df.columns):
        logger.info("Generating basin plots...")
        for fmt in ['png', 'pdf']:
            visualization.plot_ocean_basin_abundance(
                df, str(dirs['viz_plots'] / f"{organism}_ocean_basin_abundance.{fmt}"),
                'consensus_group_sp', 'ocean_basin', (10, 6), 300)
        logger.info("✓ Basin plots")

        if 'assigned_sp' in df.columns:
            for fmt in ['png', 'pdf']:
                visualization.plot_ocean_basin_abundance_faceted(
                    df, str(dirs['viz_plots'] / f"{organism}_ocean_basin_abundance_faceted.{fmt}"),
                    'consensus_group_sp', 'assigned_sp', 'ocean_basin', 9, 5, 300)
            logger.info("✓ Faceted basin plots")

    # Tree visualization
    if tree_relabeled and tree_relabeled.exists():
        logger.info("Generating tree plots...")
        safe_to_pretty = {v: k for k, v in raw_to_safe.items()}
        for fmt in ['png', 'pdf']:
            visualization.plot_phylogenetic_tree(
                tree_file=str(tree_relabeled),
                output_path=str(dirs['viz_plots'] / f"{organism}_tree.{fmt}"),
                genotype_colors=color_map if 'color_map' in locals() else None,
                show_bootstrap=True,
                bootstrap_threshold=70,
                figsize=None,  # Auto-scale based on number of tips
                dpi=300,
                label_map=safe_to_pretty
            )
        logger.info("✓ Tree plots")

    # Phase 2 plots
    logger.info("Generating Phase 2 plots...")
    for fmt in ['png', 'pdf']:
        visualization.plot_identity_distribution(
            str(diagnostics_csv), str(dirs['viz_plots'] / f"{organism}_identity_distribution.{fmt}"), (10, 6), 300)
        visualization.plot_identity_by_status(
            str(diagnostics_csv), str(dirs['viz_plots'] / f"{organism}_identity_by_status.{fmt}"), (12, 6), 300)
        visualization.plot_assignment_status(
            str(diagnostics_csv), str(dirs['viz_plots'] / f"{organism}_assignment_status.{fmt}"), (10, 6), 300)
    logger.info("✓ Phase 2 plots")

    # Summary tables
    if 'ocean_basin' in df.columns and 'consensus_group_sp' in df.columns:
        ct = pd.crosstab(df['consensus_group_sp'].fillna("Unassigned"), df['ocean_basin'].fillna("Unknown"))
        ct.to_csv(dirs['tables'] / f"{organism}_genotype_by_basin.csv")
        logger.info(f"✓ Genotype by basin table")


def main():
    parser = argparse.ArgumentParser(
        description="Run complete BOLDGenotyper pipeline with organized outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('-i', '--input', required=True, help='Input directory with consensus/annotated files')
    parser.add_argument('-o', '--output', required=True, help='Output directory (will be created)')
    parser.add_argument('-n', '--name', required=True, help='Organism name (e.g., Sphyrnidae_test)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for phylogenetics (default: 4)')
    parser.add_argument('--skip-phylo', action='store_true', help='Skip phylogenetic analysis')

    args = parser.parse_args()

    input_dir = Path(args.input).resolve()
    output_dir = Path(args.output).resolve()

    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        sys.exit(1)

    logger.info(f"Input: {input_dir}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Organism: {args.name}")

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    # Organize outputs
    try:
        dirs = organize_outputs(input_dir, output_dir, args.name)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    # Run pipeline
    try:
        # Phylogenetics
        tree_nwk = None
        tree_relabeled = None
        raw_to_safe = {}
        if not args.skip_phylo:
            tree_nwk = run_phylogenetics(dirs, args.name, args.threads)
            if tree_nwk:
                df = pd.read_csv(dirs['reports'] / f"{args.name}_annotated.csv")
                tree_relabeled, raw_to_safe = relabel_tree(tree_nwk, dirs, args.name, df)

        # Reports
        generate_reports(dirs, args.name)

        # Visualizations
        generate_visualizations(dirs, args.name, tree_relabeled, raw_to_safe)

        logger.info("\n" + "=" * 60)
        logger.info("PIPELINE COMPLETE!")
        logger.info("=" * 60)
        logger.info(f"Output directory: {output_dir}")

    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
