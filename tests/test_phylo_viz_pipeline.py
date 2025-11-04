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
from boldgenotyper import phylogenetics, visualization

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("test_phylo_viz")

# ---------- Config ----------
ORGANISM = os.environ.get("BG_ORGANISM", "Sphyrnidae_test")
OUT_DIR  = Path("tests/output").resolve()
ANNOT    = OUT_DIR / f"{ORGANISM}_annotated.csv"
CONSFA   = OUT_DIR / f"{ORGANISM}_consensus.fasta"

ALIGN    = OUT_DIR / f"{ORGANISM}_consensus_aligned.fasta"
TRIM     = OUT_DIR / f"{ORGANISM}_consensus_trimmed.fasta"
TREE     = OUT_DIR / f"{ORGANISM}_tree.nwk"

DIST_PNG = OUT_DIR / f"{ORGANISM}_distribution_map.png"
DIST_PDF = OUT_DIR / f"{ORGANISM}_distribution_map.pdf"
BARS_PNG = OUT_DIR / f"{ORGANISM}_ocean_basin_abundance.png"
BARS_PDF = OUT_DIR / f"{ORGANISM}_ocean_basin_abundance.pdf"
TREE_PNG = OUT_DIR / f"{ORGANISM}_tree.png"
TREE_PDF = OUT_DIR / f"{ORGANISM}_tree.pdf"

COLORMAP_CSV = OUT_DIR / f"{ORGANISM}_genotype_color_map.csv"
CT_CSV       = OUT_DIR / f"{ORGANISM}_genotype_by_basin.csv"

# columns in annotated CSV (adjust if yours differ)
LAT_COL   = "lat"
LON_COL   = "lon"
BASIN_COL = "ocean_basin"
GENO_COL  = "consensus_group_sp"   # we color/legend by the species-augmented label

def check_inputs():
    if not CONSFA.exists():
        raise FileNotFoundError(f"Missing consensus FASTA: {CONSFA}")
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

    # 6) Tree visualization (labels/colors by PRETTY; map SAFE->PRETTY for figure)
    try:
        visualization.plot_phylogenetic_tree(
            tree_file=str(TREE_REL),
            output_path=str(TREE_PNG),
            genotype_colors=geno_to_color,    # keys are PRETTY names
            show_bootstrap=True,
            bootstrap_threshold=70,
            figsize=(8, 10),
            dpi=300,
            label_map=safe_to_pretty,         # SAFE -> PRETTY
        )
        visualization.plot_phylogenetic_tree(
            tree_file=str(TREE_REL),
            output_path=str(TREE_PDF),
            genotype_colors=geno_to_color,
            show_bootstrap=True,
            bootstrap_threshold=70,
            figsize=(8, 10),
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

    logger.info("All done.")
    return True


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)