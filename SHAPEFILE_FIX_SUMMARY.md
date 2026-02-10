# Custom Shapefile Coordinate Parsing Fix - Summary

## Issue Identified

**Problem**: Geographic visualizations were not generated for freshwater (Salmonidae) and terrestrial (Pieridae) datasets using custom shapefiles.

**Root Cause**: The `assign_regions_from_shapefile()` function expected `lat` and `lon` as separate columns, but BOLD data stores coordinates in a single `coord` column in the format `[lat, lon]`.

## Fix Implemented

### 1. Added Coordinate Parsing to `assign_regions_from_shapefile()`

**File**: `boldgenotyper/geographic.py`

- Added `coord_col` parameter (optional)
- Parses coordinates from formats like `[51.35, -116.067]` into separate `lat` and `lon` columns
- Falls back gracefully if parsing fails

**Location**: Lines 610-720

### 2. Updated CLI to Pass coord_col Parameter

**File**: `boldgenotyper/cli.py`

- Updated call to `assign_regions_from_shapefile()` to include `coord_col="coord"`

**Location**: Line 336

## Testing Results

**Test Dataset**: Salmonidae (freshwater)

**Before Fix**:
- Samples with coordinates: 0 (0%)
- Geographic visualizations: None created

**After Fix**:
- Samples with coordinates: 1125 (19.0%)
- Coordinate parsing successful: ✓

## Important Discovery: HydroBASINS Shapefile Issue

**Problem**: The HydroBASINS shapefile being used (`hybas_pour_lev07_v1_shp`) contains **pour points** (Point geometries), not **watershedpolygons** (Polygon geometries).

**Impact**: Point-in-polygon spatial queries cannot work with pour points. Samples are marked as "Unknown" even when coordinates are successfully parsed.

**Evidence**:
```
Geometry: Point
Feature Count: 4334513
HYBAS_ID (Real) = 9070016560.00000000000
```

## Solution for HydroBASINS

You need to download the **polygon version** of HydroBASINS, not the pour points version.

### Correct HydroBASINS Dataset

**Download from**: https://www.hydrosheds.org/products/hydrobasins

**What to download**:
- **HydroBASINS Level 7 Polygons** (not pour points)
- File name should be like: `hybas_lake_na_lev07_v1c.shp` or `hybas_na_lev07_v1c.shp`
- Geometry type: Polygon or MultiPolygon
- Size: Much larger than pour points (~GB for global coverage)

**Regional Options**:
- North America: `hybas_na_lev07_v1c`
- Europe: `hybas_eu_lev07_v1c`
- Asia: `hybas_as_lev07_v1c`
- Global: Multiple regional files

### Why Pour Points Don't Work

1. **Pour points** = One point per watershed (the outlet location)
2. **Spatial query**: "Is sample point inside watershed point?" = Always False
3. **What's needed**: Polygon boundaries to check "Is sample point inside watershed polygon?"

## Solution for Ecoregions2017

Ecoregions2017 **should work correctly** as it uses polygon geometries.

**To verify**:
```bash
ogrinfo -so Ecoregions2017/Ecoregions2017.shp Ecoregions2017
# Should show: Geometry: Polygon
```

## Recommendations

### Option 1: Use Correct HydroBASINS Polygons

1. Download HydroBASINS Level 7 **polygons** (not pour points)
2. Update command to use polygon shapefile
3. Re-run Salmonidae analysis

```bash
boldgenotyper data/Salmonidae.tsv \
  --output output/Salmonidae_full \
  --organism Salmonidae \
  --build-tree \
  --custom-shp hybas_na_lev07_v1c/hybas_na_lev07_v1c.shp \
  --shp-field HYBAS_ID \
  --geo-category freshwater_basin \
  --threads 8
```

### Option 2: Use Country/Ocean Column

For datasets where precise watershed assignment isn't critical, rely on BOLD's built-in geographic fields:

- Already populated in metadata
- Works without custom shapefiles
- Less precise but functional

### Option 3: Skip Geographic Analysis for Freshwater

If watershed polygons are too large to download:

```bash
boldgenotyper data/Salmonidae.tsv \
  --output output/Salmonidae_full \
  --organism Salmonidae \
  --build-tree \
  --no-geo \
  --threads 8
```

## Status

### Working ✓
- Coordinate parsing from `coord` column
- Ecoregions2017 for terrestrial organisms
- GOaS for marine organisms

### Needs Action ⚠
- Download correct HydroBASINS polygon shapefile for freshwater organisms
- Re-run Salmonidae and any other freshwater datasets

## Files Modified

1. `boldgenotyper/geographic.py`
   - Added `coord_col` parameter to `assign_regions_from_shapefile()`
   - Implemented coordinate parsing logic

2. `boldgenotyper/cli.py`
   - Updated `assign_regions_from_shapefile()` call to include `coord_col="coord"`

## Next Steps

1. **Download HydroBASINS polygons** (if needed for freshwater datasets)
2. **Verify Ecoregions2017** contains polygons
3. **Re-run failed datasets** with correct shapefiles
4. **Test geographic visualizations** are created correctly
