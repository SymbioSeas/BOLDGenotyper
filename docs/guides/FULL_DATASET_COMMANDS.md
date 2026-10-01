# Full Dataset Analysis Commands for Manuscript

This document provides the exact commands to run all 5 case study datasets with appropriate geographic shapefiles for manuscript presentation.

Every command passes `--singleton-distance` explicitly. The values below are the
elbows of the sweeps in `benchmarking/parameter_sweep/<dataset>/recommendations.txt`
(run with BOLDGenotyper 1.0.0). Error-singleton handling changed in 2.0.0
(records of removed singletons are now assigned to their nearest haplotype when it
is within the threshold), so re-run `boldgenotyper-sweep` with 2.0.0 and check the
stability table before a final run (see `PARAMETER_SWEEP_GUIDE.md`). Pass
`--download-date` with the date of the BOLD download so it is recorded in
`run_provenance.json`, and use a new `--output` folder for each run.

The benchmark runs used the default genetic code (2, vertebrate mitochondrial) for all
five datasets. Panulirus (crustaceans) and Pieridae (insects) are invertebrates, for which
`--genetic-code 5` (invertebrate mitochondrial) is the matching code; decide which code the
final runs use and record it in the methods.

## Overview of Datasets and Shapefiles

| Dataset | Habitat | Shapefile | Field Name | Geographic Category | `--singleton-distance` |
|---------|---------|-----------|------------|---------------------|------------------------|
| Sphyrnidae | Marine | GOaS (default) | name | ocean_basin | 0.015 |
| Carcharhiniformes | Marine | GOaS (default) | name | ocean_basin | 0.015 |
| Panulirus | Marine | GOaS (default) | name | ocean_basin | 0.015 |
| Salmonidae | Freshwater | HydroBASINS level 7 | HYBAS_ID | freshwater_basin | 0.015 |
| Pieridae | Terrestrial | Ecoregions2017 | ECO_NAME | ecoregion | 0.085 |

## Commands

### 1. Sphyrnidae (Marine - Hammerhead Sharks)

```bash
boldgenotyper data/Sphyrnidae.tsv \
  --output output/Sphyrnidae_full \
  --organism Sphyrnidae \
  --singleton-distance 0.015 \
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
  --singleton-distance 0.015 \
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
  --singleton-distance 0.015 \
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
  --singleton-distance 0.015 \
  --build-tree \
  --custom-shp hybas_pour_lev07_v1_shp/hybas_pour_lev07_v1.shp \
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
  --singleton-distance 0.085 \
  --build-tree \
  --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp \
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
boldgenotyper data/Sphyrnidae.tsv --output output/Sphyrnidae_full --organism Sphyrnidae --singleton-distance 0.015 --build-tree --threads 8
boldgenotyper data/Carcharhiniformes.tsv --output output/Carcharhiniformes_full --organism Carcharhiniformes --singleton-distance 0.015 --build-tree --threads 8
boldgenotyper data/Panulirus.tsv --output output/Panulirus_full --organism Panulirus --singleton-distance 0.015 --build-tree --threads 8

# Freshwater dataset (HydroBASINS)
boldgenotyper data/Salmonidae.tsv --output output/Salmonidae_full --organism Salmonidae --singleton-distance 0.015 --build-tree --custom-shp shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev01-12_v1c.shp --shp-field HYBAS_ID --geo-category freshwater_basin --threads 8

# Terrestrial dataset (Ecoregions2017)
boldgenotyper data/Pieridae.tsv --output output/Pieridae_full --organism Pieridae --singleton-distance 0.085 --build-tree --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp --shp-field ECO_NAME --geo-category ecoregion --threads 8
```

## Expected Outputs

Each run folder holds a generated `README.md` that lists every file and defines every
column. The main files:

- `summary.csv` - record counts at each step, assignment percentages, species counts
- `records.csv` - one row per input record: status, haplotype, label, coordinates, region
- `haplotypes.csv` - one row per haplotype: members, labels, assigned taxon, nearest haplotype, suspect flags, tree membership
- `haplotypes.fasta` - haplotype sequences (`>H0001 seq_md5=... n_members=...`)
- `label_conflicts.csv` - label conflicts with their evidence
- `divergence/pairwise_p_distance.csv` - p-distances between haplotypes
- `phylogeny/tree.nwk`, `tree_labeled.nwk`, `excluded_tips.txt` - tree (from `--build-tree`)
- `geography/haplotype_by_region.csv`, `distribution_map.pdf` - records per region (`ocean_basin`, `freshwater_basin` or `ecoregion`) and the map
- `run_provenance.json` - configuration, tool versions, input md5, run status
- `report.html` - interactive report (open it from the run folder)

With `--extended`: `extended/species/species_diversity.csv` and species tables,
`extended/divergence/` (label-based barcoding gap, species-level divergence,
per-species matrices, alignment plots), `extended/metadata/`, and extra maps and bars.

## Parameter Explanations

- `--output`: Output directory for all results
- `--organism`: Organism name for file labeling (extracted from filename if omitted)
- `--build-tree`: Build phylogenetic tree using MAFFT + FastTree (required for manuscript)
- `--singleton-distance`: error-filter threshold, chosen with `boldgenotyper-sweep`
- `--custom-shp`: Path to custom shapefile (overrides default GOaS); `--shp` takes the name of a shapefile registered in `shapefiles.yaml`
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
ogrinfo -so shapefiles/GOaS_v1_20211214/goas_v01.shp goas_v01 | grep "name:"

# Check HydroBASINS (should show 'HYBAS_ID' field)
ogrinfo -so shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev01-12_v1c.shp hybas_na_lev01-12_v1c | grep "HYBAS_ID"

# Check Ecoregions (should show 'ECO_NAME' field)
ogrinfo -so shapefiles/Ecoregions2017/Ecoregions2017.shp Ecoregions2017 | grep "ECO_NAME"
```
