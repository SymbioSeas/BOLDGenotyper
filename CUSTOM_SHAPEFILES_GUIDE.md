# Custom Shapefiles Guide for BOLDGenotyper

**Author**: Steph Smith
**Last Updated**: 2025-11-25
**Version**: 1.0.0

---

## Table of Contents

- [Overview](#overview)
- [Why Use Custom Shapefiles?](#why-use-custom-shapefiles)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Detailed Workflow](#detailed-workflow)
- [Shapefile Sources](#shapefile-sources)
- [Use Case Examples](#use-case-examples)
- [Troubleshooting](#troubleshooting)
- [Citation Guidelines](#citation-guidelines)

---

## Overview

By default, BOLDGenotyper uses the GOaS (Global Oceans and Seas) shapefile for geographic analysis of marine organisms. However, many research projects focus on freshwater or terrestrial organisms that require different geographic frameworks.

The `boldgenotyper-enrich --custom-shp` command allows you to apply **any shapefile** to your genotyping results, enabling:

- Freshwater basin assignments (e.g., HydroBASINS, river networks)
- Terrestrial ecoregion assignments (e.g., WWF ecoregions, biomes)
- Custom study area boundaries
- Political boundaries (countries, states, provinces)
- Conservation areas (protected areas, national parks)
- Any other polygon-based geographic classification

This feature transforms BOLDGenotyper from a marine-only tool into a universal genotyping pipeline for **any organism type**.

---

## Why Use Custom Shapefiles?

### Scientific Applications

1. **Freshwater organisms**: Assign samples to watersheds, river basins, or hydrological units
   - Example: Salmon (*Salmonidae*) distributed across Pacific and Atlantic drainage basins
   - Example: Freshwater mussels with basin-specific populations

2. **Terrestrial organisms**: Assign samples to biomes, ecoregions, or biogeographic zones
   - Example: Forest species distributed across temperate and tropical zones
   - Example: Desert-adapted species with range-restricted populations

3. **Conservation studies**: Map genotypes to protected areas or management zones
   - Example: Population structure within and between national parks
   - Example: Connectivity between marine protected areas

4. **Political boundaries**: Assign samples to jurisdictions for management purposes
   - Example: Cross-border populations requiring international cooperation
   - Example: State-level genetic diversity assessments

### Advantages Over Manual Assignment

- **Automated**: No manual coordinate lookups required
- **Reproducible**: Exact shapefile and field documented
- **Scalable**: Works for datasets of any size
- **Visual**: Automatic regeneration of distribution maps and bar charts
- **Integrated**: Results available in interactive HTML reports

---

## Requirements

### Software Requirements

- BOLDGenotyper installed (see main README.md)
- GeoPandas and dependencies (installed with `pip install boldgenotyper[geo]`)
- Optional: QGIS, ArcGIS, or ogrinfo for shapefile inspection

### Shapefile Requirements

Your custom shapefile must include:

1. **Complete shapefile components**:
   - `.shp` (geometries)
   - `.shx` (index)
   - `.dbf` (attribute table)
   - `.prj` (projection, recommended)

2. **Polygon geometries**: Points or lines not currently supported

3. **Attribute field with region names**: The field containing your classification labels (e.g., basin names, ecoregion names)

4. **Coordinate Reference System**: Any CRS supported by GeoPandas (will be automatically reprojected to WGS84 if needed)

---

## Quick Start

### Step 1: Run Main Pipeline

First, run the standard BOLDGenotyper pipeline. You can skip geographic analysis initially:

```bash
boldgenotyper data/Salmonidae.tsv \
  --no-geo \
  --output Salmonidae_output/
```

This generates `Salmonidae_output/Salmonidae_annotated.csv` with genotype assignments but no geographic regions.

### Step 2: Apply Custom Shapefile

Apply your custom shapefile to assign samples to geographic regions:

```bash
boldgenotyper-enrich Salmonidae_output/Salmonidae_annotated.csv \
  --custom-shp data/shapefiles/hydrobasins_level4.shp \
  --shp-field HYBAS_ID \
  --geo-category freshwater_basin \
  --output Salmonidae_enriched/
```

### Step 3: Review Results

```bash
# View enriched dataset
head Salmonidae_enriched/Salmonidae_enriched.csv

# Open updated visualizations
open Salmonidae_enriched/visualization/freshwater_basin_distribution_map.pdf

# Review geographic summary
cat Salmonidae_enriched/geographic_summary.csv
```

---

## Detailed Workflow

### 1. Selecting a Shapefile

Choose a shapefile appropriate for your organism and research question:

| Organism Type | Recommended Shapefile | Source |
|--------------|----------------------|--------|
| Freshwater fish | HydroBASINS (level 3-5) | HydroSHEDS |
| Salmon/anadromous | Freshwater Ecoregions of the World | WWF |
| Terrestrial species | Terrestrial Ecoregions | WWF |
| Migratory birds | Bird Conservation Regions | BirdLife International |
| Forest species | Forest biomes | FAO |
| Custom study | Your own shapefile | - |

See [Shapefile Sources](#shapefile-sources) section for download links.

### 2. Inspecting Your Shapefile

Before using a shapefile, identify the field containing region names:

#### Option A: Using QGIS (Recommended for beginners)

1. Open QGIS
2. Layer → Add Vector Layer → Select your shapefile
3. Right-click layer → Open Attribute Table
4. Identify the column with region labels (e.g., "NAME", "ECO_NAME", "HYBAS_ID")

#### Option B: Using ogrinfo (Command-line)

```bash
# List all fields
ogrinfo -al -so path/to/shapefile.shp

# Example output:
# Layer name: hydrobasins
# Geometry: Polygon
# Feature Count: 5847
# ...
# HYBAS_ID: Integer (10.0)
# PFAF_ID: Integer64 (15.0)
# MAIN_BAS: Integer64 (15.0)
```

Use the field name exactly as shown (case-sensitive).

#### Option C: Using GeoPandas (Python)

```python
import geopandas as gpd

# Load shapefile
gdf = gpd.read_file("path/to/shapefile.shp")

# Show column names
print(gdf.columns)

# Preview attribute table
print(gdf.head())

# Check unique values in a field
print(gdf['HYBAS_ID'].unique())
```

### 3. Understanding Parameters

#### `--custom-shp` (required)

Path to your shapefile (.shp file). All companion files (.shx, .dbf, .prj) must be in the same directory.

```bash
--custom-shp data/shapefiles/hydrobasins.shp
```

#### `--shp-field` (required or default: "name")

The attribute field containing region labels. **Case-sensitive**.

Common field names:
- HydroBASINS: `HYBAS_ID`, `PFAF_ID`, `MAIN_BAS`
- WWF Ecoregions: `ECO_NAME`, `BIOME_NAME`
- Political boundaries: `NAME`, `NAME_EN`, `ISO_A3`
- Custom shapefiles: Check your attribute table

```bash
--shp-field ECO_NAME
```

#### `--geo-category` (default: "ocean_basin")

The name for your geographic classification. This appears in:
- Column headers in output CSV
- Axis labels in visualizations
- HTML report sections

Choose a descriptive name:
- For freshwater: `freshwater_basin`, `watershed`, `drainage`
- For terrestrial: `ecoregion`, `biome`, `habitat_type`
- For custom: `study_area`, `management_zone`, `population`

```bash
--geo-category freshwater_basin
```

### 4. Running the Enrichment

#### Basic Example: Freshwater Basins

```bash
boldgenotyper-enrich Salmonidae_annotated.csv \
  --custom-shp data/hydrobasins_lev04.shp \
  --shp-field MAIN_BAS \
  --geo-category freshwater_basin \
  --output Salmonidae_enriched/
```

#### Advanced Example: Multiple Operations

Combine custom shapefile with metadata merging:

```bash
boldgenotyper-enrich Salmonidae_annotated.csv \
  --custom-shp data/hydrobasins.shp \
  --shp-field HYBAS_ID \
  --geo-category watershed \
  --add-metadata field_collections.csv \
  --add-metadata lab_measurements.csv \
  --add-grouping sampling_year \
  --output Salmonidae_enriched/
```

### 5. Interpreting Results

#### Enriched Dataset

`{organism}_enriched.csv` contains all original columns plus:

- New column: `{geo_category}` (e.g., `freshwater_basin`, `ecoregion`)
- Values: Region labels from your shapefile
- Unassigned samples: `"Unknown"` (coordinates outside shapefile polygons)

#### Geographic Summary

`geographic_summary.csv` provides counts by region:

```csv
freshwater_basin,n_samples,n_genotypes,dominant_genotype
Pacific_Northwest,234,12,c1_n156
Atlantic_Northeast,189,8,c3_n98
Great_Lakes,67,5,c7_n34
```

#### Updated Visualizations

New visualizations are generated in `visualization/`:

- `{geo_category}_distribution_map.png/pdf`: Global map colored by custom regions
- `{geo_category}_distribution_bar.png/pdf`: Genotype composition by region
- `{geo_category}_summary_table.html`: Interactive summary

---

## Shapefile Sources

### Freshwater

**HydroBASINS (HydroSHEDS)**
- Description: Global watershed boundaries at 12 hierarchical levels
- Download: https://www.hydrosheds.org/products/hydrobasins
- Citation: Lehner & Grill (2013) Global river hydrography and network routing. *Hydrolog. Earth Syst. Sci.* 17:1093-1123.
- Best for: Freshwater fish, amphibians, aquatic invertebrates
- Recommended levels: Level 4-6 for regional studies, Level 7-8 for local studies

**Freshwater Ecoregions of the World (FEOW)**
- Description: 426 freshwater ecoregions based on fish distributions
- Download: https://www.feow.org/
- Citation: Abell et al. (2008) Freshwater ecoregions of the world. *BioScience* 58:403-414.
- Best for: Biogeographic studies, conservation planning

### Terrestrial

**WWF Terrestrial Ecoregions**
- Description: 867 terrestrial ecoregions across 14 biomes
- Download: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
- Citation: Olson et al. (2001) Terrestrial ecoregions of the world. *BioScience* 51:933-938.
- Best for: Terrestrial species, habitat associations

**RESOLVE Ecoregions (Updated)**
- Description: Updated ecoregions with improved boundaries
- Download: https://ecoregions.appspot.com/
- Citation: Dinerstein et al. (2017) An ecoregion-based approach to protecting half the terrestrial realm. *BioScience* 67:534-545.
- Best for: Modern conservation assessments

### Political Boundaries

**Natural Earth**
- Description: Public domain political boundaries at multiple scales
- Download: https://www.naturalearthdata.com/
- Scales: 1:10m (high detail), 1:50m (medium), 1:110m (low)
- Best for: Countries, states/provinces, administrative regions

**GADM**
- Description: Global administrative boundaries (5 levels)
- Download: https://gadm.org/
- Citation: GADM database version 4.1 (2022)
- Best for: Fine-scale administrative units

### Conservation

**World Database on Protected Areas (WDPA)**
- Description: Global protected area boundaries
- Download: https://www.protectedplanet.net/
- Citation: UNEP-WCMC & IUCN (2023) Protected Planet
- Best for: Conservation genetics, protected area assessments

---

## Use Case Examples

### Example 1: Freshwater Fish (Salmonidae)

**Research Question**: Do salmon genotypes correspond to Pacific vs. Atlantic drainage basins?

**Workflow**:

```bash
# 1. Run main pipeline (skip geographic analysis for now)
boldgenotyper data/Salmonidae.tsv \
  --clustering-threshold 0.02 \
  --no-geo \
  --output Salmonidae_output/

# 2. Download HydroBASINS level 3 (major drainage basins)
# https://www.hydrosheds.org/products/hydrobasins

# 3. Apply custom shapefile
boldgenotyper-enrich Salmonidae_output/Salmonidae_annotated.csv \
  --custom-shp data/hydrobasins_lev03_v1c.shp \
  --shp-field MAIN_BAS \
  --geo-category drainage_basin \
  --output Salmonidae_enriched/

# 4. Review results
cat Salmonidae_enriched/geographic_summary.csv
open Salmonidae_enriched/visualization/drainage_basin_distribution_map.pdf
```

**Interpretation**:
- Genotypes clustering by major drainage systems (Pacific, Atlantic, Arctic)
- Evidence for basin-specific populations
- Potential barriers to gene flow at continental divides

### Example 2: Terrestrial Species (Forest Songbirds)

**Research Question**: Are warbler genotypes associated with specific forest biomes?

**Workflow**:

```bash
# 1. Run main pipeline
boldgenotyper data/Parulidae.tsv \
  --clustering-threshold 0.03 \
  --no-geo \
  --output Parulidae_output/

# 2. Download WWF Terrestrial Ecoregions
# https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world

# 3. Apply ecoregions
boldgenotyper-enrich Parulidae_output/Parulidae_annotated.csv \
  --custom-shp data/wwf_terr_ecos/wwf_terr_ecos.shp \
  --shp-field ECO_NAME \
  --geo-category ecoregion \
  --output Parulidae_enriched/

# 4. Analyze by biome instead of ecoregion
# Re-run with BIOME field for broader classification
boldgenotyper-enrich Parulidae_output/Parulidae_annotated.csv \
  --custom-shp data/wwf_terr_ecos/wwf_terr_ecos.shp \
  --shp-field BIOME_NAME \
  --geo-category biome \
  --output Parulidae_biomes/
```

**Interpretation**:
- Genotype turnover between temperate and tropical forests
- Habitat specialists vs. generalists
- Conservation implications for habitat-specific genotypes

### Example 3: Conservation Area Assessment

**Research Question**: Is genetic diversity maintained within protected area networks?

**Workflow**:

```bash
# 1. Run main pipeline
boldgenotyper data/Endangered_species.tsv \
  --clustering-threshold 0.015 \
  --output Endangered_output/

# 2. Download WDPA protected areas
# https://www.protectedplanet.net/

# 3. Apply protected area shapefile
boldgenotyper-enrich Endangered_output/Endangered_annotated.csv \
  --custom-shp data/WDPA_polygons.shp \
  --shp-field NAME \
  --geo-category protected_area \
  --output Endangered_PA/

# 4. Check genetic diversity by protected area
cat Endangered_PA/geographic_summary.csv
```

**Interpretation**:
- Number of genotypes per protected area
- Protected areas with unique genotypes (conservation priorities)
- Connectivity between protected areas

### Example 4: Custom Study Area

**Research Question**: Population structure within our multi-site sampling transect

**Workflow**:

```bash
# 1. Create custom shapefile in QGIS:
#    - Draw polygons for each sampling site or region
#    - Add attribute field "site_name" with descriptive labels
#    - Export as shapefile

# 2. Run main pipeline
boldgenotyper data/Field_samples.tsv \
  --no-geo \
  --output Field_output/

# 3. Apply custom site boundaries
boldgenotyper-enrich Field_output/Field_annotated.csv \
  --custom-shp data/study_sites.shp \
  --shp-field site_name \
  --geo-category sampling_site \
  --add-metadata field_data.csv \
  --output Field_enriched/
```

**Interpretation**:
- Genotype frequencies at each sampling site
- Geographic patterns across your transect
- Integration with field metadata (habitat type, environmental variables)

---

## Troubleshooting

### Problem: "Shapefile not found" error

**Solution**:
```bash
# Verify file exists and path is correct
ls -la path/to/shapefile.shp

# Ensure all components present (.shp, .shx, .dbf, .prj)
ls path/to/shapefile.*

# Use absolute path if relative path fails
boldgenotyper-enrich data.csv \
  --custom-shp /full/path/to/shapefile.shp
```

### Problem: "Field 'NAME' not found in shapefile"

**Solution**:
```bash
# List all available fields
ogrinfo -al -so shapefile.shp | grep ":"

# Or use GeoPandas
python3 -c "import geopandas as gpd; print(gpd.read_file('shapefile.shp').columns)"

# Use exact field name (case-sensitive)
--shp-field HYBAS_ID  # Correct
--shp-field hybas_id  # Wrong (case mismatch)
```

### Problem: All samples assigned to "Unknown"

**Causes and solutions**:

1. **Coordinate mismatch**: Shapefile and sample coordinates use different projections
   ```bash
   # Check shapefile CRS
   ogrinfo -al -so shapefile.shp | grep "GEOGCS"

   # BOLDGenotyper expects WGS84 (EPSG:4326) for input coordinates
   # Shapefiles in other CRS are automatically reprojected
   ```

2. **Coordinates outside shapefile extent**:
   ```bash
   # Check sample coordinate range
   python3 -c "import pandas as pd; df = pd.read_csv('data.csv'); print(df[['lat', 'lon']].describe())"

   # Compare to shapefile extent
   ogrinfo -al -so shapefile.shp | grep "Extent"

   # Verify samples actually fall within shapefile polygons
   ```

3. **Invalid geometries in shapefile**:
   ```bash
   # Fix invalid geometries in QGIS:
   # Vector → Geometry Tools → Fix Geometries

   # Or use ogr2ogr
   ogr2ogr -f "ESRI Shapefile" fixed.shp original.shp -makevalid
   ```

### Problem: Visualization not generated

**Solution**:
```bash
# Check log file for errors
cat enriched_analysis/enrichment_report.txt

# Verify GeoPandas and Cartopy installed
python3 -c "import geopandas, cartopy; print('OK')"

# Install if missing
pip install geopandas cartopy

# Or conda
conda install -c conda-forge geopandas cartopy
```

### Problem: Shapefile too large / slow processing

**Solutions**:

1. **Simplify geometries**:
   ```bash
   # Use ogr2ogr to simplify
   ogr2ogr -f "ESRI Shapefile" simplified.shp original.shp -simplify 0.01
   ```

2. **Clip to study region**:
   ```bash
   # In QGIS: Vector → Geoprocessing Tools → Clip
   # Clip shapefile to bounding box of your samples
   ```

3. **Use lower resolution**:
   ```bash
   # For HydroBASINS: Use level 3-4 instead of 8-10
   # For political boundaries: Use Natural Earth 1:50m instead of 1:10m
   ```

---

## Citation Guidelines

### Citing BOLDGenotyper

When using custom shapefiles, cite BOLDGenotyper and mention the enrichment feature:

```
We used BOLDGenotyper v0.1.0 (Smith, 2025) for COI sequence genotyping
and assigned samples to geographic regions using the boldgenotyper-enrich
command with custom shapefiles.
```

### Citing Shapefiles

Always cite the shapefile source in your methods:

**HydroBASINS example**:
```
Samples were assigned to freshwater drainage basins using HydroBASINS
level 4 polygons (Lehner & Grill, 2013).

Reference:
Lehner, B., & Grill, G. (2013). Global river hydrography and network
routing: baseline data and new approaches to study the world's large
river systems. Hydrological Processes, 27(15), 2171-2186.
```

**WWF Ecoregions example**:
```
Genotypes were mapped to terrestrial ecoregions using the WWF Terrestrial
Ecoregions of the World dataset (Olson et al., 2001).

Reference:
Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell,
G. V. N., Underwood, E. C., ... & Kassem, K. R. (2001). Terrestrial ecoregions
of the world: A new map of life on Earth. BioScience, 51(11), 933-938.
```

### Example Methods Section

```
DNA barcoding data were downloaded from BOLD (Ratnasingham & Hebert, 2007)
and processed with BOLDGenotyper v0.1.0 (Smith, 2025). Sequences were
clustered at 97% identity (clustering threshold = 0.03) to define consensus
genotypes. Samples were assigned to freshwater ecoregions using the Freshwater
Ecoregions of the World shapefile (Abell et al., 2008) via the
boldgenotyper-enrich command. Geographic visualizations and summary statistics
were generated automatically.

References:
- Abell, R., et al. (2008). Freshwater ecoregions of the world: A new map
  of biogeographic units for freshwater biodiversity conservation. BioScience,
  58(5), 403-414.
- Ratnasingham, S., & Hebert, P. D. N. (2007). BOLD: The Barcode of Life Data
  System. Molecular Ecology Notes, 7(3), 355-364.
- Smith, S. (2025). BOLDGenotyper: Automated COI sequence genotyping and
  biogeographic analysis. https://github.com/SymbioSeas/BOLDGenotyper
```

---

## Advanced Tips

### Combining Multiple Shapefiles

Assign samples to multiple geographic frameworks:

```bash
# First enrichment: Ecoregions
boldgenotyper-enrich original_annotated.csv \
  --custom-shp ecoregions.shp \
  --shp-field ECO_NAME \
  --geo-category ecoregion \
  --output step1/

# Second enrichment: Protected areas (using output from first)
boldgenotyper-enrich step1/{organism}_enriched.csv \
  --custom-shp protected_areas.shp \
  --shp-field PA_NAME \
  --geo-category protected_area \
  --output step2/
```

Result: Dataset with both `ecoregion` and `protected_area` columns.

### Creating Custom Shapefiles

Use QGIS to create study-specific shapefiles:

1. **New shapefile**: Layer → Create Layer → New Shapefile Layer
2. **Draw polygons**: Digitize your study regions
3. **Add attributes**: Create field for region names
4. **Populate data**: Enter descriptive labels for each polygon
5. **Save**: Export as ESRI Shapefile
6. **Use in BOLDGenotyper**: Apply with `--custom-shp`

### Nested Geographic Hierarchies

For hierarchical analyses (e.g., basin → sub-basin → stream):

```bash
# Coarse scale (major basins)
boldgenotyper-enrich data.csv \
  --custom-shp hydrobasins_lev03.shp \
  --shp-field MAIN_BAS \
  --geo-category major_basin \
  --output level3/

# Fine scale (sub-basins)
boldgenotyper-enrich data.csv \
  --custom-shp hydrobasins_lev06.shp \
  --shp-field HYBAS_ID \
  --geo-category sub_basin \
  --output level6/

# Compare hierarchical levels
# Analyze genetic structure at multiple spatial scales
```

---

## Support and Feedback

For questions or issues with custom shapefiles:

1. **Check the log**: Review `enriched_analysis/enrichment_report.txt`
2. **Validate shapefile**: Use QGIS or ogrinfo to inspect
3. **Search issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
4. **Ask for help**: Open a new issue with:
   - Command used
   - Error message
   - Shapefile source and field name
   - Sample coordinate range

---

## Future Enhancements

Planned features for custom shapefile support:

- Automatic shapefile download for common sources (HydroBASINS, WWF ecoregions)
- Support for raster-based classifications
- Multi-polygon assignment (samples overlapping multiple regions)
- Hierarchical geographic analysis (nested regions)
- Integration with Environmental DNA (eDNA) sampling metadata

---

**Created**: 2025-11-25
**Documentation version**: 1.0.0
**BOLDGenotyper version**: 0.1.0+

For general BOLDGenotyper documentation, see [README.md](README.md).
For population genetics export formats, see [POPGEN_EXPORT_GUIDE.md](POPGEN_EXPORT_GUIDE.md).
