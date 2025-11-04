#!/usr/bin/env python3
"""
Test pipeline for phylogenetics + visualization on prior genotyping outputs.

Inputs expected from previous pipeline:
  - tests/output/{organism}_consensus.fasta
  - tests/output/{organism}_annotated.csv   (processid, lat, lon, ocean_basin, consensus_group, ...)

Outputs produced here:
  - tests/output/{organism}_consensus_aligned.fasta
  - tests/output/{organism}_consensus_trimmed.fasta
  - tests/output/{organism}_tree.nwk
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
import pandas as pd

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
LAT_COL  = "lat"
LON_COL  = "lon"
BASIN_COL = "ocean_basin"
GENO_COL = "consensus_group_sp"

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

def main() -> bool:
    check_inputs()
    df = pd.read_csv(ANNOT)
    if GENO_COL not in df.columns:
        raise ValueError(f"Column '{GENO_COL}' not found in {ANNOT}")
    genotypes = sorted([g for g in df[GENO_COL].dropna().unique().tolist()])

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
            # fallback shell call
            if not which("mafft"):
                raise RuntimeError("MAFFT not found and phylogenetics.align_mafft() not available.")
            run(["mafft", "--auto", str(CONSFA)], cwd=None)
            # MAFFT writes to stdout; capture to file
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
                # IQ-TREE writes *.treefile next to alignment
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
    # Build genotype→color map
    colors = visualization.get_genotype_colors(len(genotypes))
    geno_to_color = {g: c for g, c in zip(genotypes, colors)}
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
        # also a PDF
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

    # Map consensus_group -> consensus_group_sp so tree tips get species-augmented labels
    label_map = {}
    if {"consensus_group", "consensus_group_sp"}.issubset(df.columns):
        tmp = df[["consensus_group", "consensus_group_sp"]].dropna().drop_duplicates()
        label_map = dict(zip(tmp["consensus_group"], tmp["consensus_group_sp"]))

    # 6) Tree visualization (colors by genotype)
    try:
        # For consensus trees, tip labels are consensus IDs (genotypes).
        visualization.plot_phylogenetic_tree(
            tree_file=str(TREE),
            output_path=str(TREE_PNG),
            genotype_colors=geno_to_color,
            show_bootstrap=True,
            bootstrap_threshold=70,
            figsize=(8, 10),
            dpi=300,
            label_map=label_map,
        )
        visualization.plot_phylogenetic_tree(
            tree_file=str(TREE),
            output_path=str(TREE_PDF),
            genotype_colors=geno_to_color,
            show_bootstrap=True,
            bootstrap_threshold=70,
            figsize=(8, 10),
            dpi=300,
            label_map=label_map,
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