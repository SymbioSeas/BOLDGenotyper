# Custom Shapefile Support - Current Status and Plan

## Problem Summary

The `--shp`, `--shp-field`, and `--geo-category` arguments **do NOT exist in the main `boldgenotyper` command**. These parameters are only available in the separate `boldgenotyper-enrich` post-processing command, which:

- ✓ Can add custom geographic region columns to the output CSV
- ✗ Does NOT run the full geographic analysis pipeline
- ✗ Does NOT create geographic maps, trees, diversity/divergence summaries by region

## Current Architecture

### Main `boldgenotyper` Command
- **Supports**: GOaS (Global Ocean and Seas) shapefile ONLY (hardcoded)
- **Arguments**: `--no-geo` to skip geographic analysis entirely
- **Output**: Full pipeline with geographic analysis (maps, trees, summaries by ocean basin)
- **Location**: `boldgenotyper/cli.py` lines 305-400, `main()` function

### `boldgenotyper-enrich` Command
- **Supports**: Custom shapefiles via `--custom-shp`, `--shp-field`, `--geo-category`
- **Purpose**: Post-processing enrichment to add metadata or update geographic assignments
- **Output**: Enriched CSV with custom geographic column added
- **Limitation**: Does NOT regenerate full geographic analysis outputs
- **Location**: `boldgenotyper/cli.py` lines 2009-2200, `main_enrich()` function

## Options for Full Dataset Analysis

### Option 1: Add Custom Shapefile Support to Main Command (RECOMMENDED)

**What needs to be done:**

1. **Add CLI arguments to main parser** (`boldgenotyper/cli.py` ~line 2313)
   ```python
   parser.add_argument('--custom-shp', type=Path, dest='shapefile_path', ...)
   parser.add_argument('--shp-field', type=str, default='name', ...)
   parser.add_argument('--geo-category', type=str, default='ocean_basin', ...)
   ```

2. **Modify geographic analysis section** (`boldgenotyper/cli.py` lines 305-400)
   - Check if custom shapefile is provided
   - If yes: use `geographic.assign_regions_from_shapefile()` instead of `assign_ocean_basins()`
   - Pass custom `geo_category` name through to all downstream functions

3. **Update geographic enhancement** (`boldgenotyper/geographic_enhancement.py`)
   - Accept `geo_category` parameter (currently hardcoded to 'ocean_basin')
   - Use dynamic category name in:
     - Diversity calculations
     - Divergence calculations
     - Map generation
     - Tree generation
     - Output file naming

4. **Test with all three habitat types**
   - Marine (GOaS) - existing functionality
   - Freshwater (HydroBASINS) - new
   - Terrestrial (Ecoregions2017) - new

**Estimated complexity:** Moderate (3-4 hours)
**Benefits:** Clean, maintainable, manuscript-ready solution

### Option 2: Use GOaS for All Datasets (TEMPORARY WORKAROUND)

Run all 5 datasets with default GOaS shapefile. Accept that:
- Marine datasets (Sphyrnidae, Carcharhiniformes, Panulirus) work correctly
- Salmonidae (freshwater) will show ocean basin assignments (incorrect but harmless)
- Pieridae (terrestrial) will mostly show "Unknown" for non-coastal samples

**Commands:**
```bash
# All datasets use default GOaS
boldgenotyper data/Sphyrnidae.tsv --output output/Sphyrnidae_full --organism Sphyrnidae --build-tree --threads 8
boldgenotyper data/Carcharhiniformes.tsv --output output/Carcharhiniformes_full --organism Carcharhiniformes --build-tree --threads 8
boldgenotyper data/Panulirus.tsv --output output/Panulirus_full --organism Panulirus --build-tree --threads 8
boldgenotyper data/Salmonidae.tsv --output output/Salmonidae_full --organism Salmonidae --build-tree --threads 8
boldgenotyper data/Pieridae.tsv --output output/Pieridae_full --organism Pieridae --build-tree --threads 8
```

**Manuscript note:** "Geographic analysis used ocean basin assignments for all datasets as proof-of-concept. Future versions will support freshwater basin and terrestrial ecoregion shapefiles."

### Option 3: Skip Geographic Analysis for Non-Marine (NOT RECOMMENDED)

Run marine datasets with GOaS, skip geo for others:

```bash
# Marine (with geo)
boldgenotyper data/Sphyrnidae.tsv --output output/Sphyrnidae_full --organism Sphyrnidae --build-tree --threads 8
boldgenotyper data/Carcharhiniformes.tsv --output output/Carcharhiniformes_full --organism Carcharhiniformes --build-tree --threads 8
boldgenotyper data/Panulirus.tsv --output output/Panulirus_full --organism Panulirus --build-tree --threads 8

# Non-marine (no geo)
boldgenotyper data/Salmonidae.tsv --output output/Salmonidae_full --organism Salmonidae --build-tree --no-geo --threads 8
boldgenotyper data/Pieridae.tsv --output output/Pieridae_full --organism Pieridae --build-tree --no-geo --threads 8
```

**Drawback:** Loses valuable geographic visualization/analysis for 40% of case studies

## Recommendation

**For manuscript quality:** Implement **Option 1** (add custom shapefile support to main command)

**For immediate progress:** Use **Option 2** (run all with GOaS) to generate results while implementing Option 1

## Implementation Priority

If implementing Option 1, tackle in this order:

1. **High Priority (required for basic functionality):**
   - CLI argument parsing
   - Geographic assignment with custom shapefile
   - Pass `geo_category` through pipeline

2. **Medium Priority (for complete geographic analysis):**
   - Update geographic enhancement to use dynamic category names
   - Update map generation for custom categories
   - Update diversity/divergence summaries

3. **Low Priority (polish):**
   - Update plot labels and legends
   - Documentation updates
   - Examples in help text

## Files That Need Modification

1. `boldgenotyper/cli.py` (main command)
   - Add arguments: ~line 2313
   - Modify geographic analysis: lines 305-400
   - Pass parameters to runner: ~line 450

2. `boldgenotyper/geographic_enhancement.py`
   - Accept `geo_category` parameter
   - Use dynamic column names
   - Update file naming

3. `boldgenotyper/reports.py` (if geographic plots are embedded)
   - Update plot titles/labels if hardcoded

4. Documentation:
   - Update README with custom shapefile examples
   - Update CLI help text
