#!/bin/bash
# Run full dataset analysis for all 5 case study organisms with appropriate shapefiles
# Usage: bash scripts/run_full_datasets.sh

set -e  # Exit on error
set -u  # Exit on undefined variable

# Configuration
THREADS=8
OUTPUT_DIR="output"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "BOLDGenotyper Full Dataset Analysis"
echo "=========================================="
echo ""
echo "Running 5 datasets with appropriate geographic shapefiles:"
echo "  1. Sphyrnidae (marine - GOaS)"
echo "  2. Carcharhiniformes (marine - GOaS)"
echo "  3. Panulirus (marine - GOaS)"
echo "  4. Salmonidae (freshwater - HydroBASINS)"
echo "  5. Pieridae (terrestrial - Ecoregions2017)"
echo ""
echo "Threads: ${THREADS}"
echo "Output: ${OUTPUT_DIR}/"
echo ""

# Function to run dataset and report status
run_dataset() {
    local name=$1
    shift
    echo "=========================================="
    echo "Starting: ${name}"
    echo "=========================================="
    start_time=$(date +%s)

    if "$@"; then
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        echo "✓ ${name} completed in ${elapsed}s"
    else
        echo "✗ ${name} failed"
        exit 1
    fi
    echo ""
}

# 1. Sphyrnidae (Marine - Hammerhead Sharks)
run_dataset "Sphyrnidae" \
    boldgenotyper data/Sphyrnidae.tsv \
    --output "${OUTPUT_DIR}/Sphyrnidae_full" \
    --organism Sphyrnidae \
    --build-tree \
    --threads "${THREADS}"

# 2. Carcharhiniformes (Marine - Requiem Sharks)
run_dataset "Carcharhiniformes" \
    boldgenotyper data/Carcharhiniformes.tsv \
    --output "${OUTPUT_DIR}/Carcharhiniformes_full" \
    --organism Carcharhiniformes \
    --build-tree \
    --threads "${THREADS}"

# 3. Panulirus (Marine - Spiny Lobsters)
run_dataset "Panulirus" \
    boldgenotyper data/Panulirus.tsv \
    --output "${OUTPUT_DIR}/Panulirus_full" \
    --organism Panulirus \
    --build-tree \
    --threads "${THREADS}"

# 4. Salmonidae (Freshwater - Salmon and Trout)
run_dataset "Salmonidae" \
    boldgenotyper data/Salmonidae.tsv \
    --output "${OUTPUT_DIR}/Salmonidae_full" \
    --organism Salmonidae \
    --build-tree \
    --custom-shp hybas_pour_lev07_v1_shp/hybas_pour_lev07_v1.shp \
    --shp-field HYBAS_ID \
    --geo-category freshwater_basin \
    --threads "${THREADS}"

# 5. Pieridae (Terrestrial - Butterflies)
run_dataset "Pieridae" \
    boldgenotyper data/Pieridae.tsv \
    --output "${OUTPUT_DIR}/Pieridae_full" \
    --organism Pieridae \
    --build-tree \
    --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp \
    --shp-field ECO_NAME \
    --geo-category ecoregion \
    --threads "${THREADS}"

echo "=========================================="
echo "All datasets completed successfully!"
echo "=========================================="
echo ""
echo "Output directories:"
echo "  - ${OUTPUT_DIR}/Sphyrnidae_full/"
echo "  - ${OUTPUT_DIR}/Carcharhiniformes_full/"
echo "  - ${OUTPUT_DIR}/Panulirus_full/"
echo "  - ${OUTPUT_DIR}/Salmonidae_full/"
echo "  - ${OUTPUT_DIR}/Pieridae_full/"
echo ""
echo "Key outputs per dataset:"
echo "  - haplotypes/*_haplotype_stats.csv"
echo "  - species_analysis/species_diversity.csv"
echo "  - species_analysis/species_divergence_summary.csv"
echo "  - phylogeny/*_tree.nwk"
echo "  - geographic_analysis/*_map.html"
echo "  - *_report.html"
