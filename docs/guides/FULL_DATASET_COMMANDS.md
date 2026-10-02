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

Panulirus (crustaceans) and Pieridae (insects) are invertebrates: their commands pass
`--genetic-code 5` (invertebrate mitochondrial), to `boldgenotyper-sweep` as well, so the
sweep checks the same records as the final run. The other datasets use the default code 2
(vertebrate mitochondrial). The benchmark runs used code 2 for all five datasets.

## Overview of Datasets and Shapefiles

| Dataset | Habitat | Shapefile | Field Name | Geographic Category | `--singleton-distance` | `--genetic-code` |
|---------|---------|-----------|------------|---------------------|------------------------|------------------|
| Sphyrnidae | Marine | GOaS (default) | name | ocean_basin | 0.015 | 2 (default) |
| Carcharhiniformes | Marine | GOaS (default) | name | ocean_basin | 0.015 | 2 (default) |
| Panulirus | Marine | GOaS (default) | name | ocean_basin | 0.015 | 5 |
| Salmonidae | Freshwater | BasinATLAS v10 level 3 | PFAF_ID | freshwater_basin | 0.015 | 2 (default) |
| Pieridae | Terrestrial | Ecoregions2017 | REALM | biogeographic_realm | 0.085 | 5 |

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
  --genetic-code 5 \
  --build-tree \
  --threads 8
```

**Notes:**
- Uses default GOaS shapefile for ocean basin assignment
- Genetic code 5 (invertebrate mitochondrial)
- Output: `output/Panulirus_full/`
- ~2.3M input file

### 4. Salmonidae (Freshwater - Salmon and Trout)

```bash
boldgenotyper data/Salmonidae.tsv \
  --output output/Salmonidae_full \
  --organism Salmonidae \
  --singleton-distance 0.015 \
  --build-tree \
  --custom-shp shapefiles/BasinATLAS_Data_v10_shp/BasinATLAS_v10_shp/BasinATLAS_v10_lev03.shp \
  --shp-field PFAF_ID \
  --geo-category freshwater_basin \
  --threads 8
```

**Notes:**
- Uses BasinATLAS (HydroSHEDS/HydroATLAS, https://www.hydrosheds.org/hydroatlas) level 3
  for freshwater basin assignment: 292 major drainage basins, global
- PFAF_ID (Pfafstetter code) identifies each basin
- Output: `output/Salmonidae_full/`
- ~9.5M input file, second largest dataset

### 5. Pieridae (Terrestrial - Butterflies)

```bash
boldgenotyper data/Pieridae.tsv \
  --output output/Pieridae_full \
  --organism Pieridae \
  --singleton-distance 0.085 \
  --genetic-code 5 \
  --build-tree \
  --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp \
  --shp-field REALM \
  --geo-category biogeographic_realm \
  --threads 8
```

**Notes:**
- Uses Ecoregions2017 for terrestrial ecoregion assignment
- Genetic code 5 (invertebrate mitochondrial)
- REALM groups the ecoregions into 8 biogeographic realms (e.g., "Palearctic", "Neotropic");
  ECO_NAME would give the 847 individual ecoregions
- The region column is `biogeographic_realm` (`--geo-category`), a name that is not a BOLD
  field, so BOLD's own `realm` and `ecoregion` fields stay in `records.csv` (the 1.0.0
  runs used `ecoregion`, which replaced BOLD's field)
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
boldgenotyper data/Panulirus.tsv --output output/Panulirus_full --organism Panulirus --singleton-distance 0.015 --genetic-code 5 --build-tree --threads 8

# Freshwater dataset (BasinATLAS level 3)
boldgenotyper data/Salmonidae.tsv --output output/Salmonidae_full --organism Salmonidae --singleton-distance 0.015 --build-tree --custom-shp shapefiles/BasinATLAS_Data_v10_shp/BasinATLAS_v10_shp/BasinATLAS_v10_lev03.shp --shp-field PFAF_ID --geo-category freshwater_basin --threads 8

# Terrestrial dataset (Ecoregions2017)
boldgenotyper data/Pieridae.tsv --output output/Pieridae_full --organism Pieridae --singleton-distance 0.085 --genetic-code 5 --build-tree --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp --shp-field REALM --geo-category biogeographic_realm --threads 8
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
- `--genetic-code`: NCBI genetic code for reading-frame QC (2 vertebrate mitochondrial, the default; 5 invertebrate mitochondrial); pass the same code to `boldgenotyper-sweep`
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

# Check BasinATLAS level 3 (should show 'PFAF_ID' field)
ogrinfo -so shapefiles/BasinATLAS_Data_v10_shp/BasinATLAS_v10_shp/BasinATLAS_v10_lev03.shp BasinATLAS_v10_lev03 | grep "PFAF_ID"

# Check Ecoregions (should show 'REALM' field)
ogrinfo -so shapefiles/Ecoregions2017/Ecoregions2017.shp Ecoregions2017 | grep "REALM"
```
