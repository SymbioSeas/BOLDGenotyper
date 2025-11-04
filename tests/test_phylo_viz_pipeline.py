#!/usr/bin/env python3
"""
Test pipeline for phylogenetics + visualization on prior genotyping outputs.

Inputs expected from previous pipeline:
  - tests/output/{organism}_consensus.fasta
  - tests/output/{organism}_annotated.csv   (processid, lat, lon, ocean_basin, consensus_group, consensus_group_sp, ...)

Outputs produced here:
  - tests/output/{organism}_consensus_aligned.fasta
  - tests/output/{organism}_consensus_trimmed.fasta
  - tests/output/{organism}_tree.nwk
  - tests/output/{organism}_tree.relabeled.nwk           # NEW: tips relabeled to safe pretty names
  - tests/output/{organism}_tree.png / .pdf
  - tests/output/{organism}_distribution_map.png / .pdf
  - tests/output/{organism}_ocean_basin_abundance.png / .pdf
  - tests/output/{organism}_genotype_color_map.csv
  - tests/output/{organism}_genotype_by_basin.csv
"""

import os
import sys
import shutil
import logging
from pathlib import Path
import subprocess
import re
import pandas as pd
from Bio import Phylo  # used for relabeling/reading/writing Newick

# local imports
from boldgenotyper import phylogenetics, visualization, reports

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("test_phylo_viz")

# ---------- Config ----------
ORGANISM = os.environ.get("BG_ORGANISM", "Euprymna")
OUT_DIR  = Path("tests/Euprymna_output").resolve()
ANNOT    = OUT_DIR / f"{ORGANISM}_annotated.csv"
CONSFA   = OUT_DIR / f"{ORGANISM}_consensus.fasta"

ALIGN    = OUT_DIR / f"{ORGANISM}_consensus_aligned.fasta"
TRIM     = OUT_DIR / f"{ORGANISM}_consensus_trimmed.fasta"
TREE     = OUT_DIR / f"{ORGANISM}_tree.nwk"

DIST_PNG = OUT_DIR / f"{ORGANISM}_distribution_map.png"
DIST_PDF = OUT_DIR / f"{ORGANISM}_distribution_map.pdf"
DIST_FACET_PNG = OUT_DIR / f"{ORGANISM}_distribution_map_faceted.png"
DIST_FACET_PDF = OUT_DIR / f"{ORGANISM}_distribution_map_faceted.pdf"
BARS_PNG = OUT_DIR / f"{ORGANISM}_ocean_basin_abundance.png"
BARS_PDF = OUT_DIR / f"{ORGANISM}_ocean_basin_abundance.pdf"
BARS_FACET_PNG = OUT_DIR / f"{ORGANISM}_ocean_basin_abundance_faceted.png"
BARS_FACET_PDF = OUT_DIR / f"{ORGANISM}_ocean_basin_abundance_faceted.pdf"
TREE_PNG = OUT_DIR / f"{ORGANISM}_tree.png"
TREE_PDF = OUT_DIR / f"{ORGANISM}_tree.pdf"

COLORMAP_CSV = OUT_DIR / f"{ORGANISM}_genotype_color_map.csv"
CT_CSV       = OUT_DIR / f"{ORGANISM}_genotype_by_basin.csv"

# Phase 1 Reports
DIAG_CSV     = OUT_DIR / f"{ORGANISM}_diagnostics.csv"
CONFLICTS_CSV = OUT_DIR / f"{ORGANISM}_taxonomy_conflicts.csv"
SUMMARY_CSV   = OUT_DIR / f"{ORGANISM}_assignment_summary.csv"
CONSENSUS_CHAR_CSV = OUT_DIR / f"{ORGANISM}_consensus_characterization.csv"

# Phase 2 Reports and Visualizations
IDENTITY_DIST_PNG = OUT_DIR / f"{ORGANISM}_identity_distribution.png"
IDENTITY_DIST_PDF = OUT_DIR / f"{ORGANISM}_identity_distribution.pdf"
IDENTITY_STATUS_PNG = OUT_DIR / f"{ORGANISM}_identity_by_status.png"
IDENTITY_STATUS_PDF = OUT_DIR / f"{ORGANISM}_identity_by_status.pdf"
STATUS_BREAKDOWN_PNG = OUT_DIR / f"{ORGANISM}_assignment_status.png"
STATUS_BREAKDOWN_PDF = OUT_DIR / f"{ORGANISM}_assignment_status.pdf"
SEQQ_CSV = OUT_DIR / f"{ORGANISM}_sequence_quality.csv"

# columns in annotated CSV (adjust if yours differ)
LAT_COL   = "lat"
LON_COL   = "lon"
BASIN_COL = "ocean_basin"
GENO_COL  = "consensus_group_sp"   # we color/legend by the species-augmented label

def check_inputs():
    if not CONSFA.exists():
        raise FileNotFoundError(
            f"Missing consensus FASTA: {CONSFA}"
            f"Hint: either set BG_ORGANISM to match the output folder "
            f"(currently ORGANISM='{ORGANISM}', OUT_DIR='{OUT_DIR}') "
            f"or rerun the main pipeline to generate {CONSFA.name} in that folder."
        )
    if not ANNOT.exists():
        raise FileNotFoundError(f"Missing annotated CSV: {ANNOT}")

def which(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def run(cmd, cwd=None):
    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)

def _safe_label(s: str) -> str:
    """Replace any character not in [A-Za-z0-9._-] with underscore for Newick."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)

def main() -> bool:
    check_inputs()
    df = pd.read_csv(ANNOT)
    if GENO_COL not in df.columns:
        raise ValueError(f"Column '{GENO_COL}' not found in {ANNOT}")

    # final display names (pretty) = consensus_group_sp
    genotypes_pretty = sorted(df[GENO_COL].dropna().unique().tolist())
    colors = visualization.get_genotype_colors(len(genotypes_pretty))
    color_by_pretty = dict(zip(genotypes_pretty, colors))

    # --------- PHYLOGENETICS PIPE ----------
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Align consensus sequences
    try:
        if hasattr(phylogenetics, "align_mafft"):
            ALIGN.write_text("")  # ensure path exists if module writes in-place
            phylogenetics.align_mafft(
                input_fasta=str(CONSFA),
                output_fasta=str(ALIGN),
                mafft_options=["--auto"]
            )
        else:
            if not which("mafft"):
                raise RuntimeError("MAFFT not found and phylogenetics.align_mafft() not available.")
            # capture stdout to file
            with subprocess.Popen(["mafft", "--auto", str(CONSFA)],
                                  stdout=subprocess.PIPE, text=True) as p:
                ALIGN.write_text(p.stdout.read())
        logger.info("✓ Alignment written: %s", ALIGN)
    except Exception as e:
        logger.exception("Alignment failed: %s", e)
        return False

    # 2) Trim alignment
    try:
        if hasattr(phylogenetics, "trim_trimal"):
            phylogenetics.trim_trimal(
                input_fasta=str(ALIGN),
                output_fasta=str(TRIM),
                trimal_options=["-automated1"]
            )
        else:
            if not which("trimal"):
                raise RuntimeError("trimAl not found and phylogenetics.trim_trimal() not available.")
            run(["trimal", "-automated1", "-in", str(ALIGN), "-out", str(TRIM)])
        logger.info("✓ Trimmed alignment: %s", TRIM)
    except Exception as e:
        logger.exception("Trimming failed: %s", e)
        return False

    # 3) Build tree (prefer IQ-TREE; fallback to FastTree)
    try:
        built = False
        if hasattr(phylogenetics, "build_tree_iqtree") and which("iqtree2"):
            phylogenetics.build_tree_iqtree(
                input_alignment=str(TRIM),
                output_tree=str(TREE),
                threads=int(os.environ.get("BG_THREADS", "4")),
                model="GTR+G"
            )
            built = True
        elif hasattr(phylogenetics, "build_tree_fasttree") and which("fasttree"):
            phylogenetics.build_tree_fasttree(
                input_alignment=str(TRIM),
                output_tree=str(TREE)
            )
            built = True
        else:
            if which("iqtree2"):
                run(["iqtree2", "-s", str(TRIM), "-m", "GTR+G", "-T", "AUTO", "-B", "1000", "-nt", "AUTO"])
                treefile = TRIM.with_suffix(".treefile")
                if not treefile.exists():
                    raise FileNotFoundError(f"Expected {treefile} not found")
                TREE.write_text(treefile.read_text())
                built = True
            elif which("fasttree"):
                with subprocess.Popen(["fasttree", "-nt", "-gtr", "-gamma", str(TRIM)],
                                      stdout=subprocess.PIPE, text=True) as p:
                    TREE.write_text(p.stdout.read())
                built = True

        if not built:
            raise RuntimeError("No tree builder available (iqtree2 or fasttree).")

        logger.info("✓ Tree written: %s", TREE)
    except Exception as e:
        logger.exception("Tree building failed: %s", e)
        return False

    # --------- VISUALIZATION PIPE ----------
    # 3.5) Relabel Newick tips to SAFE(pretty) names on disk, and build SAFE->PRETTY map for plotting
    TREE_REL = OUT_DIR / f"{ORGANISM}_tree.relabeled.nwk"

    # raw consensus -> pretty label
    raw_to_pretty = {}
    if {"consensus_group", "consensus_group_sp"}.issubset(df.columns):
        tmp = df[["consensus_group", "consensus_group_sp"]].dropna().drop_duplicates()
        raw_to_pretty = dict(zip(tmp["consensus_group"], tmp["consensus_group_sp"]))

    # raw -> SAFE(pretty)
    raw_to_safe = {k: _safe_label(v) for k, v in raw_to_pretty.items()}
    # SAFE(pretty) -> PRETTY (for figure labels)
    safe_to_pretty = {raw_to_safe[k]: v for k, v in raw_to_pretty.items()}

    try:
        tree_obj = Phylo.read(str(TREE), "newick")
        for clade in tree_obj.get_terminals():
            # Newick tips are raw consensus IDs at this point
            if clade.name in raw_to_safe:
                clade.name = raw_to_safe[clade.name]
        Phylo.write(tree_obj, str(TREE_REL), "newick")  # <-- FIXED format ("newick")
        logger.info("✓ Relabeled Newick written: %s", TREE_REL)
    except Exception as e:
        logger.exception("Relabeling Newick failed: %s", e)
        TREE_REL = TREE  # fall back

    # Build genotype→color map keyed by PRETTY names for legend/labels
    geno_to_color = {g: c for g, c in zip(genotypes_pretty, colors)}
    pd.DataFrame({"genotype": list(geno_to_color.keys()),
                  "color": list(geno_to_color.values())}).to_csv(COLORMAP_CSV, index=False)

    # 4) Distribution map (lat/lon points colored by genotype)
    try:
        visualization.plot_distribution_map(
            df=df,
            output_path=str(DIST_PNG),
            genotype_column=GENO_COL,
            latitude_col=LAT_COL if LAT_COL in df.columns else "lat",
            longitude_col=LON_COL if LON_COL in df.columns else "lon",
            figsize=(10, 6),
            dpi=300,
        )
        visualization.plot_distribution_map(
            df=df,
            output_path=str(DIST_PDF),
            genotype_column=GENO_COL,
            latitude_col=LAT_COL if LAT_COL in df.columns else "lat",
            longitude_col=LON_COL if LON_COL in df.columns else "lon",
            figsize=(10, 6),
            dpi=300,
        )
        logger.info("✓ Distribution map: %s / %s", DIST_PNG, DIST_PDF)
    except Exception as e:
        logger.exception("Distribution map failed: %s", e)

    # 4b) Distribution map - FACETED by species
    try:
        if "assigned_sp" in df.columns:
            visualization.plot_distribution_map_faceted(
                df=df,
                output_path=str(DIST_FACET_PNG),
                genotype_column=GENO_COL,
                species_column="assigned_sp",
                latitude_col=LAT_COL if LAT_COL in df.columns else "lat",
                longitude_col=LON_COL if LON_COL in df.columns else "lon",
                width=10,
                height_per_species=5,
                dpi=300,
            )
            visualization.plot_distribution_map_faceted(
                df=df,
                output_path=str(DIST_FACET_PDF),
                genotype_column=GENO_COL,
                species_column="assigned_sp",
                latitude_col=LAT_COL if LAT_COL in df.columns else "lat",
                longitude_col=LON_COL if LON_COL in df.columns else "lon",
                width=10,
                height_per_species=5,
                dpi=300,
            )
            logger.info("✓ Faceted distribution map: %s / %s", DIST_FACET_PNG, DIST_FACET_PDF)
        else:
            logger.warning("Column 'assigned_sp' not found; skipping faceted distribution map")
    except Exception as e:
        logger.exception("Faceted distribution map failed: %s", e)

    # 5) Ocean basin stacked bars
    try:
        visualization.plot_ocean_basin_abundance(
            df=df,
            output_path=str(BARS_PNG),
            genotype_column=GENO_COL,
            basin_column=BASIN_COL if BASIN_COL in df.columns else "ocean_basin",
            figsize=(10, 6),
            dpi=300,
        )
        visualization.plot_ocean_basin_abundance(
            df=df,
            output_path=str(BARS_PDF),
            genotype_column=GENO_COL,
            basin_column=BASIN_COL if BASIN_COL in df.columns else "ocean_basin",
            figsize=(10, 6),
            dpi=300,
        )
        logger.info("✓ Basin abundance plots: %s / %s", BARS_PNG, BARS_PDF)
    except Exception as e:
        logger.exception("Basin abundance plot failed: %s", e)

    # 5b) Ocean basin stacked bars - FACETED by species
    try:
        if "assigned_sp" in df.columns:
            visualization.plot_ocean_basin_abundance_faceted(
                df=df,
                output_path=str(BARS_FACET_PNG),
                genotype_column=GENO_COL,
                species_column="assigned_sp",
                basin_column=BASIN_COL if BASIN_COL in df.columns else "ocean_basin",
                width=9,
                height_per_species=5,
                dpi=300,
            )
            visualization.plot_ocean_basin_abundance_faceted(
                df=df,
                output_path=str(BARS_FACET_PDF),
                genotype_column=GENO_COL,
                species_column="assigned_sp",
                basin_column=BASIN_COL if BASIN_COL in df.columns else "ocean_basin",
                width=9,
                height_per_species=5,
                dpi=300,
            )
            logger.info("✓ Faceted basin abundance plots: %s / %s", BARS_FACET_PNG, BARS_FACET_PDF)
        else:
            logger.warning("Column 'assigned_sp' not found; skipping faceted basin plots")
    except Exception as e:
        logger.exception("Faceted basin abundance plot failed: %s", e)

    # 6) Tree visualization (labels/colors by PRETTY; map SAFE->PRETTY for figure)
    try:
        visualization.plot_phylogenetic_tree(
            tree_file=str(TREE_REL),
            output_path=str(TREE_PNG),
            genotype_colors=geno_to_color,    # keys are PRETTY names
            show_bootstrap=True,
            bootstrap_threshold=70,
            figsize=None,  # Auto-scale based on number of tips
            dpi=300,
            label_map=safe_to_pretty,         # SAFE -> PRETTY
        )
        visualization.plot_phylogenetic_tree(
            tree_file=str(TREE_REL),
            output_path=str(TREE_PDF),
            genotype_colors=geno_to_color,
            show_bootstrap=True,
            bootstrap_threshold=70,
            figsize=None,  # Auto-scale based on number of tips
            dpi=300,
            label_map=safe_to_pretty,
        )
        logger.info("✓ Tree figure: %s / %s", TREE_PNG, TREE_PDF)
    except Exception as e:
        logger.exception("Tree visualization failed: %s", e)

    # 7) QC: genotype × basin crosstab
    try:
        if BASIN_COL in df.columns and GENO_COL in df.columns:
            ct = pd.crosstab(
                df[GENO_COL].fillna("Unassigned"),
                df[BASIN_COL].fillna("Unknown")
            )
            ct.to_csv(CT_CSV)
            logger.info("✓ Genotype × Basin crosstab: %s", CT_CSV)
    except Exception as e:
        logger.exception("Crosstab failed: %s", e)

    # --------- PHASE 1 REPORTS ----------
    # 8) Taxonomy conflicts report
    try:
        if DIAG_CSV.exists():
            reports.generate_taxonomy_conflicts_report(
                annotated_csv=str(ANNOT),
                diagnostics_csv=str(DIAG_CSV),
                output_csv=str(CONFLICTS_CSV)
            )
            logger.info("✓ Taxonomy conflicts report: %s", CONFLICTS_CSV)
        else:
            logger.warning(f"Diagnostics file not found: {DIAG_CSV}; skipping conflict report")
    except Exception as e:
        logger.exception("Taxonomy conflicts report failed: %s", e)

    # 9) Assignment summary statistics
    try:
        if DIAG_CSV.exists():
            reports.generate_assignment_summary(
                annotated_csv=str(ANNOT),
                diagnostics_csv=str(DIAG_CSV),
                output_csv=str(SUMMARY_CSV)
            )
            logger.info("✓ Assignment summary: %s", SUMMARY_CSV)
        else:
            logger.warning(f"Diagnostics file not found: {DIAG_CSV}; skipping summary")
    except Exception as e:
        logger.exception("Assignment summary failed: %s", e)

    # 10) Consensus group characterization
    try:
        if DIAG_CSV.exists():
            reports.generate_consensus_characterization(
                annotated_csv=str(ANNOT),
                diagnostics_csv=str(DIAG_CSV),
                output_csv=str(CONSENSUS_CHAR_CSV)
            )
            logger.info("✓ Consensus characterization: %s", CONSENSUS_CHAR_CSV)
        else:
            logger.warning(f"Diagnostics file not found: {DIAG_CSV}; skipping characterization")
    except Exception as e:
        logger.exception("Consensus characterization failed: %s", e)

    # --------- PHASE 2 VISUALIZATIONS & REPORTS ----------
    # 11) Identity distribution plot
    try:
        if DIAG_CSV.exists():
            visualization.plot_identity_distribution(
                diagnostics_csv=str(DIAG_CSV),
                output_path=str(IDENTITY_DIST_PNG),
                figsize=(10, 6),
                dpi=300
            )
            visualization.plot_identity_distribution(
                diagnostics_csv=str(DIAG_CSV),
                output_path=str(IDENTITY_DIST_PDF),
                figsize=(10, 6),
                dpi=300
            )
            logger.info("✓ Identity distribution plots: %s / %s", IDENTITY_DIST_PNG, IDENTITY_DIST_PDF)
        else:
            logger.warning(f"Diagnostics file not found: {DIAG_CSV}; skipping identity distribution")
    except Exception as e:
        logger.exception("Identity distribution plot failed: %s", e)

    # 12) Identity by status plot
    try:
        if DIAG_CSV.exists():
            visualization.plot_identity_by_status(
                diagnostics_csv=str(DIAG_CSV),
                output_path=str(IDENTITY_STATUS_PNG),
                figsize=(12, 6),
                dpi=300
            )
            visualization.plot_identity_by_status(
                diagnostics_csv=str(DIAG_CSV),
                output_path=str(IDENTITY_STATUS_PDF),
                figsize=(12, 6),
                dpi=300
            )
            logger.info("✓ Identity by status plots: %s / %s", IDENTITY_STATUS_PNG, IDENTITY_STATUS_PDF)
        else:
            logger.warning(f"Diagnostics file not found: {DIAG_CSV}; skipping identity by status")
    except Exception as e:
        logger.exception("Identity by status plot failed: %s", e)

    # 13) Assignment status breakdown
    try:
        if DIAG_CSV.exists():
            visualization.plot_assignment_status(
                diagnostics_csv=str(DIAG_CSV),
                output_path=str(STATUS_BREAKDOWN_PNG),
                figsize=(10, 6),
                dpi=300
            )
            visualization.plot_assignment_status(
                diagnostics_csv=str(DIAG_CSV),
                output_path=str(STATUS_BREAKDOWN_PDF),
                figsize=(10, 6),
                dpi=300
            )
            logger.info("✓ Assignment status plots: %s / %s", STATUS_BREAKDOWN_PNG, STATUS_BREAKDOWN_PDF)
        else:
            logger.warning(f"Diagnostics file not found: {DIAG_CSV}; skipping status breakdown")
    except Exception as e:
        logger.exception("Assignment status plot failed: %s", e)

    # 14) Sequence quality metrics
    try:
        reports.generate_sequence_quality_metrics(
            annotated_csv=str(ANNOT),
            output_csv=str(SEQQ_CSV)
        )
        logger.info("✓ Sequence quality metrics: %s", SEQQ_CSV)
    except Exception as e:
        logger.exception("Sequence quality metrics failed: %s", e)

    logger.info("All done.")
    return True


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)