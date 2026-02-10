# BOLDGenotyper Publication Readiness Review

## Date: February 5, 2026

## Executive Summary

Comprehensive review of all 5 case study datasets (Sphyrnidae, Carcharhiniformes, Panulirus, Salmonidae, Pieridae) reveals several issues that need addressing before publication.

## Dataset Output Review

### All Datasets Successfully Generated

✓ Sphyrnidae (marine - GOaS) - 750 samples
✓ Carcharhiniformes (marine - GOaS) - 9,244 samples
✓ Panulirus (marine - GOaS) - 640 samples
✓ Salmonidae (freshwater - custom shapefile) - 5,546 samples
✓ Pieridae (terrestrial - custom shapefile) - 3,484 samples

### Core Outputs Present

All datasets successfully produced:
- ✓ Annotated CSV with genotype assignments
- ✓ Haplotype statistics and FASTA files
- ✓ Species-level analysis (diversity, divergence)
- ✓ Phylogenetic trees
- ✓ Quality control reports
- ✓ HTML summary reports
- ✓ Metadata analysis
- ✓ Taxonomy assignment tables

## Critical Issues Identified

### 1. **HTML Report Visualization Paths** (HIGH PRIORITY)

**Problem**: Visualization search patterns in `reports.py` expect PNG files, but visualizations are saved as PDF/SVG in subdirectories.

**Evidence**:
```python
# Pattern looks for:
'pattern': f'{organism}_distribution_bar.png'

# But files actually exist at:
visualization/pdf/{organism}_distribution_bar.pdf
visualization/svg/{organism}_distribution_bar.svg
```

**Impact**:
- Visualizations may not appear in HTML reports
- Users see "No visualization images found" message

**Location**: `boldgenotyper/reports.py` lines 3453-3508

**Status**: Code DOES attempt fallback to PDF/SVG, needs verification

---

### 2. **Missing Visualizations for Custom Shapefiles** (HIGH PRIORITY)

**Problem**: Datasets using custom shapefiles (Salmonidae, Pieridae) are missing distribution bar charts.

**Evidence**:
- Sphyrnidae (GOaS): 7 PDF visualizations + 7 SVG
- Salmonidae (custom): 3 PDF visualizations + 3 SVG (missing 4 bar charts)
- Pieridae (custom): 3 PDF visualizations + 3 SVG (missing 4 bar charts)

**Missing for custom shapefiles**:
- `{organism}_distribution_bar.pdf` - Relative abundance by region
- `{organism}_totaldistribution_bar.pdf` - Total abundance by region
- `{organism}_distribution_bar_faceted.pdf` - Faceted by species/haplotype
- `{organism}_distribution_bar_species_faceted.pdf` - Species-faceted stacked

**Only present for custom shapefiles**:
- `{organism}_distribution_map.pdf` - Geographic map
- `{organism}_distribution_map_faceted.pdf` - Faceted map
- `{organism}_distribution_map_species_faceted.pdf` - Species-faceted map

**Root Cause**: Visualization code likely has hardcoded references to `ocean_basin` column

**Location**: Likely in visualization generation code (needs investigation)

---

### 3. **HydroBASINS Pour Points Issue** (MEDIUM PRIORITY)

**Problem**: Salmonidae using HydroBASINS pour points (Point geometries) instead of watershed polygons.

**Evidence**:
- 1125 samples parsed with coordinates (19%)
- But only "Unknown" basin assignments (spatial queries fail on points)

**Impact**:
- No meaningful watershed assignments for Salmonidae
- Geographic analysis incomplete

**Solution Required**:
- Download HydroBASINS polygon shapefile
- Re-run Salmonidae analysis
- OR document limitation in manuscript

**Status**: Documented in SHAPEFILE_FIX_SUMMARY.md

---

### 4. **Visualization File Path Documentation** (LOW PRIORITY)

**Problem**: HTML report shows incorrect file paths in documentation text.

**Example**:
```html
📁 Files: <code>visualization/{image_name}</code>
```

**Should be**:
```html
📁 Files: <code>visualization/pdf/{image_name}</code>, <code>visualization/svg/{image_name}</code>
```

**Location**: `boldgenotyper/reports.py` line 3585

---

### 5. **Metadata Visualization Directory Structure** (LOW PRIORITY)

**Observation**: Metadata visualizations are in `visualization/metadata/pdf/` subdirectory.

**Current structure**:
```
visualization/
├── pdf/                    # Main visualizations (maps, bar charts)
├── svg/                    # Main visualizations (SVG format)
├── json/                   # Interactive plot data
└── metadata/
    └── pdf/                # Metadata visualizations
```

**Question**: Should metadata visualizations also have SVG versions?

**Status**: Works but inconsistent

---

## Comprehensive Action Plan

### Phase 1: Critical Fixes (Before Publication)

#### Task 1.1: Fix Visualization Generation for Custom Shapefiles
- [ ] Investigate visualization code for hardcoded `ocean_basin` references
- [ ] Update to use dynamic `geo_category` parameter
- [ ] Ensure bar charts are generated for all geographic categories
- [ ] Re-run Salmonidae and Pieridae to generate missing visualizations

**Files to check**:
- `boldgenotyper/visualization.py` or equivalent
- `boldgenotyper/plots.py` or equivalent
- Any code that generates bar charts

#### Task 1.2: Verify HTML Report Visualization Paths
- [ ] Test that HTML reports correctly embed all available visualizations
- [ ] Verify fallback from PNG → PDF → SVG works correctly
- [ ] Fix file path documentation in report footer
- [ ] Ensure interactive plot JSON data is loaded correctly

**Files to modify**:
- `boldgenotyper/reports.py` (lines 3580-3590)

#### Task 1.3: Verify All Dataset HTML Reports
- [ ] Open each HTML report in browser
- [ ] Check that all expected visualizations appear
- [ ] Verify interactive plots work (if applicable)
- [ ] Check for any broken images or missing sections

**Datasets to verify**:
- Sphyrnidae_summary_report.html
- Carcharhiniformes_summary_report.html
- Panulirus_summary_report.html
- Salmonidae_summary_report.html
- Pieridae_summary_report.html

---

### Phase 2: Data Quality Improvements

#### Task 2.1: Resolve HydroBASINS Issue (Optional)
- [ ] Download HydroBASINS Level 7 polygon shapefile
- [ ] Re-run Salmonidae with correct shapefile
- [ ] Verify watershed assignments are populated
- [ ] OR document limitation in manuscript methods

#### Task 2.2: Verify Ecoregions2017 for Pieridae
- [ ] Check if ecoregion assignments worked
- [ ] Verify map visualizations show regional patterns
- [ ] Confirm polygon shapefile was used correctly

---

### Phase 3: Code Quality & Documentation

#### Task 3.1: Code Review for Hardcoded References
- [ ] Search codebase for hardcoded `ocean_basin` strings
- [ ] Replace with dynamic `geo_category` parameter
- [ ] Ensure consistent parameter passing throughout pipeline

**Search command**:
```bash
grep -r "ocean_basin" --include="*.py" boldgenotyper/ | grep -v "default.*ocean_basin"
```

#### Task 3.2: Update Documentation
- [ ] README.md - Verify custom shapefile examples
- [ ] Docstrings - Check all geo-related functions
- [ ] Error messages - Ensure they're helpful and accurate

#### Task 3.3: Test Suite
- [ ] Create unit tests for coordinate parsing
- [ ] Test visualization generation with different geo_categories
- [ ] Test HTML report generation with various visualization sets

---

### Phase 4: Publication Preparation

#### Task 4.1: Generate Final Outputs
- [ ] Re-run any datasets with fixes applied
- [ ] Verify all visualizations are publication-quality
- [ ] Export high-resolution PDFs for manuscript figures
- [ ] Create supplementary file bundles

#### Task 4.2: Quality Assurance Checklist
- [ ] All 5 datasets have complete outputs
- [ ] All HTML reports render correctly
- [ ] All visualizations are scientifically accurate
- [ ] No broken links or missing files
- [ ] Consistent styling across all outputs

---

## Investigation Priorities

### Immediate Investigation Required

1. **Where are bar chart visualizations generated?**
   - Need to find the code that creates distribution_bar plots
   - Check if it has hardcoded `ocean_basin` column references
   - Determine why custom shapefile datasets skip bar charts

2. **Why do custom shapefiles only generate maps?**
   - Maps work (suggests coordinate parsing successful)
   - Bar charts missing (suggests column name issue)
   - Need to trace visualization generation logic

3. **Are all expected files being created?**
   - Check if JSON files exist for all visualizations
   - Verify PDF and SVG versions are both created
   - Ensure file naming is consistent

---

## File Structure Expectations

### Expected visualization directory structure:
```
visualization/
├── pdf/
│   ├── {organism}_distribution_bar.pdf
│   ├── {organism}_totaldistribution_bar.pdf
│   ├── {organism}_distribution_bar_faceted.pdf
│   ├── {organism}_distribution_bar_species_faceted.pdf
│   ├── {organism}_distribution_map.pdf
│   ├── {organism}_distribution_map_faceted.pdf
│   └── {organism}_distribution_map_species_faceted.pdf
├── svg/
│   └── (same files as pdf/)
├── json/
│   ├── {organism}_distribution_bar_data.json
│   ├── {organism}_totaldistribution_bar_data.json
│   ├── {organism}_distribution_bar_faceted_data.json
│   └── {organism}_distribution_map_data.json
└── metadata/
    └── pdf/
        └── (metadata visualizations)
```

### Currently missing from custom shapefile datasets:
- All `*_bar.pdf` and `*_bar.svg` files
- Corresponding JSON data files for bar charts

---

## Next Steps

1. **Investigate visualization generation code** to find hardcoded references
2. **Create detailed issue list** with file locations and line numbers
3. **Implement fixes** systematically with testing
4. **Re-run affected datasets** to verify fixes
5. **Final QA review** of all outputs

---

## Notes

- Coordinate parsing fix (Task #1 from earlier) is complete and working
- Custom shapefile support is functional for geographic assignment
- Main issue is visualization generation not adapting to custom categories
- HTML report embedding appears to work (uses base64 encoding)
