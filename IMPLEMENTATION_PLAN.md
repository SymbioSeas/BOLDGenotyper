# Implementation Plan: Publication Readiness Fixes

## Overview

This plan addresses all issues identified in the publication readiness review, with specific file locations, line numbers, and exact changes needed.

## Issue 1: Hardcoded 'ocean_basin' in Visualization Generation

### Problem
Bar chart visualizations are only generated when `'ocean_basin'` column exists, causing custom shapefile datasets (Salmonidae, Pieridae) to skip these visualizations.

### Root Cause
**File**: `boldgenotyper/cli.py`
**Lines**: 1281, 1299, 1349, 1366

**Current code**:
```python
if 'ocean_basin' in df_final.columns and 'haplotype_sp' in df_final.columns:
```

**Problem**: This check fails for `freshwater_basin`, `ecoregion`, or any custom `geo_category`.

### Solution

Replace hardcoded `'ocean_basin'` checks with dynamic `geo_category` parameter.

#### Fix 1.1: Line 1281 - Relative abundance bar chart
**Current**:
```python
if 'ocean_basin' in df_final.columns and 'haplotype_sp' in df_final.columns:
    try:
        bar_path, bar_data = visualization.plot_ocean_basin_abundance(
            df=df_final,
            output_path=str(get_viz_path(dirs, f"{organism}_distribution_bar", fmt)),
            genotype_column='haplotype_sp',
            basin_column='ocean_basin'
        )
```

**Fixed**:
```python
if geo_category in df_final.columns and 'haplotype_sp' in df_final.columns:
    try:
        bar_path, bar_data = visualization.plot_ocean_basin_abundance(
            df=df_final,
            output_path=str(get_viz_path(dirs, f"{organism}_distribution_bar", fmt)),
            genotype_column='haplotype_sp',
            basin_column=geo_category
        )
```

#### Fix 1.2: Line 1299 - Total abundance bar chart
**Current**:
```python
if 'ocean_basin' in df_final.columns and 'haplotype_sp' in df_final.columns:
    try:
        total_bar_path, total_bar_data = visualization.plot_ocean_basin_abundance_total(
            df=df_final,
            output_path=str(get_viz_path(dirs, f"{organism}_totaldistribution_bar", fmt)),
            genotype_column='haplotype_sp',
            basin_column='ocean_basin'
        )
```

**Fixed**:
```python
if geo_category in df_final.columns and 'haplotype_sp' in df_final.columns:
    try:
        total_bar_path, total_bar_data = visualization.plot_ocean_basin_abundance_total(
            df=df_final,
            output_path=str(get_viz_path(dirs, f"{organism}_totaldistribution_bar", fmt)),
            genotype_column='haplotype_sp',
            basin_column=geo_category
        )
```

#### Fix 1.3: Line 1349 - Species relative abundance
**Current**:
```python
if 'ocean_basin' in species_df.columns and 'primary_species' in species_df.columns:
    try:
        species_bar_path, species_bar_data = visualization.plot_ocean_basin_abundance(
            df=species_df,
            output_path=str(get_viz_path(dirs, f"{organism}_species_distribution_bar", fmt)),
            genotype_column='primary_species',
            basin_column='ocean_basin'
        )
```

**Fixed**:
```python
if geo_category in species_df.columns and 'primary_species' in species_df.columns:
    try:
        species_bar_path, species_bar_data = visualization.plot_ocean_basin_abundance(
            df=species_df,
            output_path=str(get_viz_path(dirs, f"{organism}_species_distribution_bar", fmt)),
            genotype_column='primary_species',
            basin_column=geo_category
        )
```

#### Fix 1.4: Line 1366 - Species total abundance
**Current**:
```python
if 'ocean_basin' in species_df.columns and 'primary_species' in species_df.columns:
    try:
        species_total_bar_path, species_total_bar_data = visualization.plot_ocean_basin_abundance_total(
            df=species_df,
            output_path=str(get_viz_path(dirs, f"{organism}_species_totaldistribution_bar", fmt)),
            genotype_column='primary_species',
            basin_column='ocean_basin'
        )
```

**Fixed**:
```python
if geo_category in species_df.columns and 'primary_species' in species_df.columns:
    try:
        species_total_bar_path, species_total_bar_data = visualization.plot_ocean_basin_abundance_total(
            df=species_df,
            output_path=str(get_viz_path(dirs, f"{organism}_species_totaldistribution_bar", fmt)),
            genotype_column='primary_species',
            basin_column=geo_category
        )
```

### Additional Changes Needed

#### Fix 1.5: Pass geo_category to visualization phase
The `geo_category` parameter must be available in the visualization phase.

**Check location**: Search for where visualization phase is executed in `run_pipeline()`

**Ensure**: `geo_category` variable is accessible in visualization code block

---

## Issue 2: Faceted Visualizations with Custom Shapefiles

### Problem
Faceted bar charts (lines 1487, 1506) also have hardcoded column references.

### Solution

#### Fix 2.1: Line ~1480-1500 - Faceted visualizations
Need to check these sections and replace:
- Column name checks for `ocean_basin`
- Parameter passing for `basin_column='ocean_basin'`

**Action**: Read lines 1480-1520 and apply same fixes as above

---

## Issue 3: HTML Report File Path Documentation

### Problem
File path documentation in HTML report shows incorrect paths.

### Location
**File**: `boldgenotyper/reports.py`
**Line**: 3585

### Current Code
```python
html += f'<p style="color: #666; font-size: 0.85em; margin-bottom: 10px;">📁 Files: <code>{file_dir}/{image_name}</code>, <code>{pdf_name}</code>'
```

### Fixed Code
```python
# Determine actual file locations based on visualization structure
pdf_path = f'visualization/pdf/{pdf_name}'
svg_name = Path(image_name).with_suffix('.svg').name
svg_path = f'visualization/svg/{svg_name}'
html += f'<p style="color: #666; font-size: 0.85em; margin-bottom: 10px;">📁 Files: <code>{pdf_path}</code>, <code>{svg_path}</code>'
```

---

## Issue 4: Visualization Axis Labels and Titles

### Problem
Visualization functions likely have hardcoded "Ocean Basin" labels.

### Files to Check
- `boldgenotyper/visualization.py` - Look for axis labels, titles
- `boldgenotyper/plots.py` - Check plot generation functions

### Search Command
```bash
grep -n "Ocean Basin\|ocean basin" boldgenotyper/visualization.py boldgenotyper/plots.py
```

### Solution
Update axis labels to use dynamic category names:
- "Ocean Basin" → "Geographic Region" or use geo_category value
- "ocean_basin" → geo_category parameter

---

## Issue 5: Species CSV Column Names

### Problem
Species assignments CSV may have hardcoded `ocean_basin` column.

### Files to Check
- Where species_assignments.csv is created
- Ensure it uses dynamic geo_category column name

### Search Command
```bash
grep -rn "species_assignments\.csv" boldgenotyper/*.py | grep -i write
```

---

## Implementation Order

### Phase 1: Core Visualization Fixes (Highest Priority)
1. ✅ Fix cli.py lines 1281, 1299, 1349, 1366 (ocean_basin → geo_category)
2. ✅ Check and fix faceted visualization code (~lines 1480-1520)
3. ✅ Verify geo_category is passed to visualization phase
4. ✅ Test with Salmonidae dataset

### Phase 2: Visualization Function Updates
1. ✅ Update visualization.py axis labels/titles
2. ✅ Update plots.py function parameters
3. ✅ Ensure all plot functions accept basin_column parameter
4. ✅ Test with all datasets

### Phase 3: HTML Report Fixes
1. ✅ Fix file path documentation in reports.py
2. ✅ Verify visualization embedding works
3. ✅ Test HTML reports in browser

### Phase 4: Final Testing
1. ✅ Re-run Salmonidae with fixes
2. ✅ Re-run Pieridae with fixes
3. ✅ Verify all 5 datasets produce complete outputs
4. ✅ Check all HTML reports render correctly

---

## Testing Checklist

### After Each Fix
- [ ] Code compiles without errors
- [ ] Test dataset runs without crashes
- [ ] Visualizations are created
- [ ] HTML report loads without errors

### After All Fixes
- [ ] Sphyrnidae: All 7 visualizations present
- [ ] Carcharhiniformes: All 7 visualizations present
- [ ] Panulirus: Expected visualizations present
- [ ] Salmonidae: All 7 visualizations present (including bars)
- [ ] Pieridae: All 7 visualizations present (including bars)

### Expected Visualizations Per Dataset
All datasets should have:
1. Distribution map
2. Distribution map (faceted)
3. Distribution map (species-faceted)
4. Distribution bar (relative abundance)
5. Total distribution bar (absolute counts)
6. Distribution bar (faceted)
7. Distribution bar (species-faceted)

**Note**: Some may be skipped if insufficient data, but the code should attempt to create them.

---

## Verification Commands

### Check visualization files
```bash
for dataset in Sphyrnidae Carcharhiniformes Panulirus Salmonidae Pieridae; do
    echo "$dataset:"
    ls output/${dataset}_full/visualization/pdf/*.pdf | wc -l
    ls output/${dataset}_full/visualization/svg/*.svg | wc -l
done
```

### Check for missing bar charts
```bash
for dataset in Salmonidae Pieridae; do
    echo "Checking $dataset for bar charts:"
    ls output/${dataset}_full/visualization/pdf/*bar*.pdf 2>/dev/null || echo "  MISSING: Bar charts"
done
```

### Test HTML reports
```bash
for dataset in Sphyrnidae Carcharhiniformes Panulirus Salmonidae Pieridae; do
    open output/${dataset}_full/${dataset}_summary_report.html
    sleep 2
done
```

---

## Success Criteria

✅ **Complete** when:
1. All 5 datasets generate complete visualization sets
2. Custom shapefile datasets have bar charts
3. HTML reports display all visualizations
4. No hardcoded 'ocean_basin' references in visualization code
5. File paths in HTML are accurate
6. All plots have appropriate axis labels

---

## Documentation Updates Needed

After implementation:
1. Update README.md with complete feature list
2. Add visualization gallery to documentation
3. Document geo_category parameter usage
4. Add troubleshooting section for custom shapefiles

---

## Notes

- Keep backward compatibility - ocean_basin should still work
- geo_category parameter should have default value 'ocean_basin'
- Error messages should be informative about missing columns
- Log which geographic category is being used
