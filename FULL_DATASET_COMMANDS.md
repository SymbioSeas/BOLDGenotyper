# Full Dataset Analysis Commands for Manuscript

This document provides the exact commands to run all 5 case study datasets with appropriate geographic shapefiles for manuscript presentation.

## Overview of Datasets and Shapefiles

| Dataset | Habitat | Shapefile | Field Name | Geographic Category |
|---------|---------|-----------|------------|---------------------|
| Sphyrnidae | Marine | GOaS (default) | name | ocean_basin |
| Carcharhiniformes | Marine | GOaS (default) | name | ocean_basin |
| Panulirus | Marine | GOaS (default) | name | ocean_basin |
| Salmonidae | Freshwater | HydroBASINS level 7 | HYBAS_ID | freshwater_basin |
| Pieridae | Terrestrial | Ecoregions2017 | ECO_NAME | ecoregion |

## Commands

### 1. Sphyrnidae (Marine - Hammerhead Sharks)

```bash
boldgenotyper data/Sphyrnidae.tsv \
  --output output/Sphyrnidae_full \
  --organism Sphyrnidae \
  --build-tree \
  --threads 8
```

**Notes:**
- Uses default GOaS shapefile for ocean basin assignment
- Output: `output/Sphyrnidae_full/`
- ~1.4M input file, ~750 sequences based on benchmarking

### 2. Carcharhiniformes (Marine - Requiem Sharks)

```bash
boldgenotyper data/Carcharhiniformes.tsv \
  --output output/Carcharhiniformes_full \
  --organism Carcharhiniformes \
  --build-tree \
  --threads 8
```

**Notes:**
- Uses default GOaS shapefile for ocean basin assignment
- Output: `output/Carcharhiniformes_full/`
- ~14M input file, largest dataset (9244 sequences)

### 3. Panulirus (Marine - Spiny Lobsters)

```bash
boldgenotyper data/Panulirus.tsv \
  --output output/Panulirus_full \
  --organism Panulirus \
  --build-tree \
  --threads 8
```

**Notes:**
- Uses default GOaS shapefile for ocean basin assignment
- Output: `output/Panulirus_full/`
- ~2.3M input file

### 4. Salmonidae (Freshwater - Salmon and Trout)

```bash
boldgenotyper data/Salmonidae.tsv \
  --output output/Salmonidae_full \
  --organism Salmonidae \
  --build-tree \
  --shp hybas_pour_lev07_v1_shp/hybas_pour_lev07_v1.shp \
  --shp-field HYBAS_ID \
  --geo-category freshwater_basin \
  --threads 8
```

**Notes:**
- Uses HydroBASINS level 7 pour points for freshwater basin assignment
- HYBAS_ID provides unique basin identifiers
- Output: `output/Salmonidae_full/`
- ~9.5M input file, second largest dataset

### 5. Pieridae (Terrestrial - Butterflies)

```bash
boldgenotyper data/Pieridae.tsv \
  --output output/Pieridae_full \
  --organism Pieridae \
  --build-tree \
  --shp Ecoregions2017/Ecoregions2017.shp \
  --shp-field ECO_NAME \
  --geo-category ecoregion \
  --threads 8
```

**Notes:**
- Uses Ecoregions2017 for terrestrial ecoregion assignment
- ECO_NAME provides human-readable ecoregion names (e.g., "Tropical and subtropical moist broadleaf forests")
- Output: `output/Pieridae_full/`
- ~31M input file, largest file size (most metadata)

## Running All Datasets

To run all datasets sequentially:

```bash
# Create output directory
mkdir -p output

# Marine datasets (GOaS default)
boldgenotyper data/Sphyrnidae.tsv --output output/Sphyrnidae_full --organism Sphyrnidae --build-tree --threads 8
boldgenotyper data/Carcharhiniformes.tsv --output output/Carcharhiniformes_full --organism Carcharhiniformes --build-tree --threads 8
boldgenotyper data/Panulirus.tsv --output output/Panulirus_full --organism Panulirus --build-tree --threads 8

# Freshwater dataset (HydroBASINS)
boldgenotyper data/Salmonidae.tsv --output output/Salmonidae_full --organism Salmonidae --build-tree --shp hybas_pour_lev07_v1_shp/hybas_pour_lev07_v1.shp --shp-field HYBAS_ID --geo-category freshwater_basin --threads 8

# Terrestrial dataset (Ecoregions2017)
boldgenotyper data/Pieridae.tsv --output output/Pieridae_full --organism Pieridae --build-tree --shp Ecoregions2017/Ecoregions2017.shp --shp-field ECO_NAME --geo-category ecoregion --threads 8
```

## Expected Outputs

Each dataset will produce the following key outputs in its respective output directory:

### Core Analysis Files
- `haplotypes/{organism}_haplotype_stats.csv` - Per-haplotype statistics
- `haplotypes/{organism}_sequences.fasta` - Consensus sequences
- `species_analysis/species_diversity.csv` - Species-level diversity metrics
- `species_analysis/species_divergence_summary.csv` - Divergence metrics and barcode gaps
- `species_analysis/within_species_divergence_summary.csv` - Within-species variation
- `genotypes/{organism}_genotype_assignments.csv` - Sample-to-haplotype assignments
- `phylogeny/{organism}_tree.nwk` - Phylogenetic tree (from --build-tree)

### Geographic Analysis Files
For marine datasets (GOaS):
- `geographic_analysis/ocean_basin_diversity.csv`
- `geographic_analysis/ocean_basin_divergence.csv`
- `geographic_analysis/ocean_basin_tree.nwk`
- `geographic_analysis/ocean_basin_map.html`

For Salmonidae (HydroBASINS):
- `geographic_analysis/freshwater_basin_diversity.csv`
- `geographic_analysis/freshwater_basin_divergence.csv`
- `geographic_analysis/freshwater_basin_tree.nwk`
- `geographic_analysis/freshwater_basin_map.html`

For Pieridae (Ecoregions2017):
- `geographic_analysis/ecoregion_diversity.csv`
- `geographic_analysis/ecoregion_divergence.csv`
- `geographic_analysis/ecoregion_tree.nwk`
- `geographic_analysis/ecoregion_map.html`

### Report Files
- `{organism}_report.html` - Interactive HTML report with all visualizations

## Parameter Explanations

- `--output`: Output directory for all results
- `--organism`: Organism name for file labeling (extracted from filename if omitted)
- `--build-tree`: Build phylogenetic tree using MAFFT + FastTree (required for manuscript)
- `--shp`: Path to custom shapefile (overrides default GOaS)
- `--shp-field`: Shapefile attribute field containing region names
- `--geo-category`: Label for geographic category in outputs
- `--threads`: Number of parallel threads (adjust based on your system)

## Computational Requirements

- **Runtime**: Varies by dataset size (5-60 minutes per dataset on 8 cores)
- **Memory**: ~4-8GB RAM depending on dataset
- **Disk space**: ~100-500MB per dataset output
- **Dependencies**: MAFFT and FastTree must be installed for --build-tree

## Verifying Shapefile Setup

To confirm shapefiles are correctly configured:

```bash
# Check GOaS (should show 'name' field)
ogrinfo -so GOaS_v1_20211214/goas_v01.shp goas_v01 | grep "name:"

# Check HydroBASINS (should show 'HYBAS_ID' field)
ogrinfo -so hybas_pour_lev07_v1_shp/hybas_pour_lev07_v1.shp hybas_pour_lev07_v1 | grep "HYBAS_ID"

# Check Ecoregions (should show 'ECO_NAME' field)
ogrinfo -so Ecoregions2017/Ecoregions2017.shp Ecoregions2017 | grep "ECO_NAME"
```
