# Carcharhinus Analysis Issues - Fixes Applied

## Issues Identified

1. ✅ **Species-level analysis completely failed**
   - Error: `bad operand type for unary ~: 'float'`
   - Cause: `is_ambiguous` column contained NaN floats, `~` operator can't be applied to floats

2. ✅ **Identity distribution plot failed**
   - Error: `'identity'` column missing from diagnostics
   - Cause: ESV-based workflow doesn't generate identity scores (uses exact matching)

3. ✅ **Assignment summary report generation failed**
   - Error: `'status'` column missing from diagnostics
   - Cause: New diagnostics format uses flags (`is_tie`, `is_low_confidence`) instead of single `status` column

4. ✅ **HTML report diagnostics path incorrect**
   - Path: Used `genotype_assignments` but actual path is `haplotype_assignments`

5. ⏳ **HTML report missing sections** (PARTIALLY FIXED - need re-run to verify)
   - Analysis summary data
   - Visualizations section
   - QC metrics
   - Species composition tables

6. ❌ **Species-faceted visualizations missing** (NOT YET IMPLEMENTED)
   - Need plots faceted by assigned species (not haplotype)
   - Example: `Carcharhinus_distribution_bar_species_faceted.pdf`

---

## Fixes Applied

### 1. Fixed Species Analysis Boolean Error

**File**: `boldgenotyper/species_analysis.py:81-88`

**Before**:
```python
confident_assignments = species_assignments[
    (species_assignments['primary_species_pct'] >= min_confidence) &
    (~species_assignments['is_ambiguous'])  # ERROR: can't negate float NaN
].copy()
```

**After**:
```python
# Ensure is_ambiguous is boolean (may contain NaN floats)
is_not_ambiguous = species_assignments['is_ambiguous'].fillna(False) == False

confident_assignments = species_assignments[
    (species_assignments['primary_species_pct'] >= min_confidence) &
    is_not_ambiguous
].copy()
```

**Impact**: Species-level analysis will now run successfully and generate:
- `species_analysis/species_assignments.csv`
- `species_analysis/species_summary.csv`
- `species_analysis/species_diversity.csv`
- `species_analysis/species_geographic_summary.csv`
- Species-level divergence matrices
- Species-specific haplotype subsets

---

### 2. Fixed Identity Distribution Plot

**File**: `boldgenotyper/visualization.py:1601-1607`

**Added check**:
```python
# Check if identity column exists (not present in ESV-based workflows)
if 'identity' not in df.columns:
    logger.warning(
        f"Identity column not found in diagnostics. "
        f"Skipping identity distribution plot (ESV-based workflow uses exact matching)."
    )
    return
```

**Impact**: Plot will skip gracefully instead of crashing when identity column is missing.

---

### 3. Fixed Assignment Summary Report

**File**: `boldgenotyper/reports.py:177-226`

**Changes**:
- Added detection for new vs old diagnostics format
- New format: derives status from flags (`is_tie`, `is_low_confidence`)
- Handles missing `identity` column gracefully

**Impact**: Assignment summary CSV will be generated successfully with correct metrics:
- Total samples
- Assignment status breakdown (assigned/tie/low_confidence)
- Haplotype counts
- Species counts
- Taxonomy conflict statistics

---

### 4. Fixed HTML Report Diagnostics Path

**File**: `boldgenotyper/reports.py:3912-3916`

**Before**:
```python
diagnostics_file = output_dir / 'genotype_assignments' / f"{organism}_diagnostics.csv"
```

**After**:
```python
# Note: diagnostics may be in 'haplotype_assignments' or 'genotype_assignments' (legacy)
diagnostics_file = output_dir / 'haplotype_assignments' / f"{organism}_diagnostics.csv"
if not diagnostics_file.exists():
    # Try legacy path
    diagnostics_file = output_dir / 'genotype_assignments' / f"{organism}_diagnostics.csv"
```

**Impact**: HTML report will correctly load diagnostics data.

---

## Testing Required

To verify fixes, re-run the analysis:

```bash
boldgenotyper data/Carcharhinus.tsv --output data/Carcharhinus_3
```

### Expected Results

#### 1. Species Analysis Directory Created
```
data/Carcharhinus_3/species_analysis/
├── species_assignments.csv
├── species_summary.csv
├── species_diversity.csv
├── species_geographic_summary.csv
└── species_faceted_haplotypes/
    ├── Carcharhinus_sp1_haplotypes.csv
    ├── Carcharhinus_sp2_haplotypes.csv
    └── ...
```

#### 2. Reports Directory Created
```
data/Carcharhinus_3/reports/
└── Carcharhinus_assignment_summary.csv
```

#### 3. HTML Report Sections Populated
- ✅ Analysis Summary (with actual metrics)
- ⏳ Visualizations section (needs verification)
- ⏳ QC metrics (needs verification)
- ⏳ Species composition tables (needs verification)

#### 4. Log Messages
Should see:
```
[INFO] PHASE 6.5: Species-Level Aggregation
[INFO] 6.5.1: Aggregating samples by assigned species...
[INFO]   Samples with confident species: X/Y (Z%)
[INFO]   Number of species: N
[INFO]   ✓ Saved species assignments: ...
[INFO]   ✓ Saved species summary: ...
[INFO] 6.5.2: Calculating species-level diversity metrics...
[INFO] 6.5.3: Generating species-level geographic summary...
[INFO] 6.5.4: Generating species-faceted haplotype subsets...
[INFO]   ✓ Species-level analysis complete
```

No errors for:
- Identity distribution plot (should skip gracefully)
- Report generation (should succeed)
- Species-level divergence (should run)

---

## Remaining Work

### Priority 1: Species-Faceted Visualizations (NOT YET IMPLEMENTED)

Need to create plots where:
- **Each facet = one assigned species**
- **Within each facet**: all haplotypes of that species, each colored uniquely

Example output files needed:
- `Carcharhinus_distribution_map_species_faceted.pdf`
- `Carcharhinus_distribution_bar_species_faceted.pdf`

**Implementation Approach**:

1. Add new visualization functions:
```python
def plot_species_faceted_distribution_maps(
    species_assignments: pd.DataFrame,
    output_path: str,
    species_column: str = "primary_species",
    haplotype_column: str = "haplotype_id",
    ...
) -> None:
    """
    Create faceted maps where each facet is one species showing all its haplotypes.

    Similar to plot_distribution_map_faceted but faceting by species instead of
    genotype.
    """
    # Group by species
    # For each species, plot all haplotypes with unique colors
    # Similar structure to existing faceted plots
```

2. Call from CLI in Phase 10 (Visualization):
```python
# Species-faceted visualizations (if species data available)
species_assignments_csv = dirs['species_analysis'] / 'species_assignments.csv'
if species_assignments_csv.exists():
    species_df = pd.read_csv(species_assignments_csv)

    # Species-faceted distribution map
    visualization.plot_species_faceted_distribution_maps(
        species_assignments=species_df,
        output_path=str(dirs['visualization'] / f"{organism}_distribution_map_species_faceted.{fmt}"),
        species_column='primary_species',
        haplotype_column='haplotype_id',
        latitude_col='lat',
        longitude_col='lon',
        ...
    )

    # Species-faceted basin chart
    visualization.plot_species_faceted_basin_charts(...)
```

### Priority 2: HTML Report Section Population

The HTML report builder needs to be verified to ensure it's actually populating all sections. This requires:

1. **Visualizations Section**:
   - Embed generated plots (SVG/PNG)
   - Add interactive maps if JSON data exists
   - Show phylogenetic tree

2. **QC Metrics Section**:
   - Load QC results from `quality_control/` directory
   - Display purity distribution
   - Display contamination heatmap
   - Show ORF validation stats

3. **Species Composition Section**:
   - Load from `taxonomy/species_by_haplotype.csv`
   - Display species composition per haplotype
   - Show mixed-species groups

**Files to check/update**:
- `boldgenotyper/reports.py` - `generate_html_report()` function
- `boldgenotyper/reports.py` - `HTMLReportBuilder` class

---

## Known Limitations

1. **Identity Distribution Plot**: Not applicable for ESV-based workflows (exact matching, no identity scores)

2. **Diagnostics Format**: Two formats exist:
   - Old: `status`, `identity`, `runner_up_identity` columns
   - New: `is_tie`, `is_low_confidence` flags only

3. **Assignment Summary**: When using new diagnostics format:
   - `n_assigned` = total - ties - low_confidence
   - `n_below_threshold` always 0 (not applicable in ESV workflow)
   - All identity statistics = NaN

---

## Files Modified

1. `boldgenotyper/species_analysis.py` - Fixed boolean negation error
2. `boldgenotyper/visualization.py` - Fixed identity distribution plot
3. `boldgenotyper/reports.py` - Fixed assignment summary and diagnostics path

---

## Next Steps

1. **Re-run analysis** to verify fixes work
2. **Implement species-faceted visualizations** (new feature)
3. **Verify HTML report sections** are populated correctly
4. **Test with multiple datasets** to ensure robustness

---

**Author**: Claude Code
**Date**: December 10, 2024
**Related Issue**: Carcharhinus analysis failures
