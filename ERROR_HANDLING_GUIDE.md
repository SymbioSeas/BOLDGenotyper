# BOLDGenotyper Error Handling Guide

**Version**: 1.0
**Date**: 2025-11-18
**Purpose**: Comprehensive guide to error handling, validation, and troubleshooting

---

## Table of Contents

1. [Overview](#overview)
2. [Error Categories](#error-categories)
3. [Input Data Validation](#input-data-validation)
4. [Common Errors and Solutions](#common-errors-and-solutions)
5. [Error Messages Reference](#error-messages-reference)
6. [Recovery Mechanisms](#recovery-mechanisms)
7. [Developer Guidelines](#developer-guidelines)

---

## Overview

BOLDGenotyper implements a multi-layered error handling strategy designed to:

1. **Fail Fast**: Detect and report critical errors early
2. **Recover Gracefully**: Continue processing when possible
3. **Inform Users**: Provide actionable error messages
4. **Log Everything**: Maintain detailed logs for debugging

### Error Handling Philosophy

- **Critical Errors**: Stop pipeline immediately (e.g., missing input file, invalid sequences)
- **Recoverable Errors**: Log warning and continue with fallback (e.g., missing GOaS data)
- **Silent Handling**: Skip problematic data points with DEBUG logging (e.g., invalid coordinates)

---

## Error Categories

### 1. Critical Errors (Pipeline Stops)

**Characteristics**:
- Raise exceptions
- Pipeline cannot continue
- User must fix issue and re-run

**Examples**:
- Missing input TSV file
- Missing required columns (`processid`, `nuc`)
- All sequences invalid
- MAFFT/trimAl failures
- Invalid parameter values

### 2. Recoverable Errors (Pipeline Continues)

**Characteristics**:
- Log WARNING message
- Use fallback value or skip feature
- Recorded in `results['errors']` list
- Pipeline continues to completion

**Examples**:
- Missing GOaS shapefile → Set ocean_basin='Unknown'
- MAFFT not found for phylogeny → Skip tree building
- Visualization failures → Skip affected plots
- Report generation failures → Continue without report

### 3. Data Quality Issues (Logged, Data Skipped)

**Characteristics**:
- Log WARNING or DEBUG message
- Skip affected data points
- Pipeline continues with remaining data
- Statistics logged at end

**Examples**:
- Duplicate `processid` → Keep first occurrence
- Invalid coordinates → Set to null
- Sequences <100bp → Skip sequence
- Invalid nucleotides → Skip sequence
- Excessive N content (>50%) → Skip sequence

---

## Input Data Validation

### BOLD TSV File Validation

#### Required Columns

**Checked in**: `metadata.parse_bold_tsv()`
**Error Type**: Critical (raises `ValueError`)

```python
REQUIRED_COLUMNS = ['processid', 'nuc']
```

**Error Message**:
```
BOLD TSV is missing required columns: ['processid']
Available columns: ['sampleid', 'species', 'genus', 'lat', 'lon', ...]
```

**Solution**: Ensure TSV has headers and includes `processid` and `nuc` columns

---

#### Duplicate Sample IDs

**Checked in**: `metadata.parse_bold_tsv()`
**Error Type**: Recoverable (WARNING, auto-deduplicated)

**Behavior**:
- Keeps **first occurrence** of each duplicate `processid`
- Logs count and examples
- Removes subsequent duplicates

**Log Message**:
```
WARNING: Found 15 duplicate processids.
Examples: ['SAMPLE001', 'SAMPLE042', 'SAMPLE133', 'SAMPLE201', 'SAMPLE305'].
Keeping only first occurrence of each duplicate.
```

**Statistics Logged**:
```
Data quality summary:
  total_rows: 1500
  unique_processids: 1485  # 15 duplicates removed
  sequences_present: 1485
```

**Solution**: Review upstream data for duplicate sample IDs. First occurrence is retained.

---

#### Empty TSV File

**Checked in**: `metadata.parse_bold_tsv()`
**Error Type**: Critical (raises `pd.errors.EmptyDataError`)

**Error Message**:
```
BOLD TSV file is empty: /path/to/file.tsv
```

**Solution**: Ensure TSV file contains data rows, not just headers

---

#### Encoding Issues

**Checked in**: `metadata.parse_bold_tsv()`
**Error Type**: Recoverable (WARNING, fallback to latin-1)

**Behavior**:
- Attempts UTF-8 first
- Falls back to latin-1 if UTF-8 fails
- Logs warning on fallback

**Log Message**:
```
WARNING: UTF-8 encoding failed, trying latin-1
```

**Solution**: Save TSV in UTF-8 encoding if possible. Pipeline handles latin-1 automatically.

---

### Sequence Validation

#### Minimum Length Check

**Checked in**: `utils.validate_sequence()` and `dereplication.dereplicate_from_fasta()`
**Error Type**: Data Quality (WARNING, sequence skipped)

**Threshold**: 100 bp (default)

**Log Message**:
```
WARNING: Skipping SAMPLE001: Sequence too short (87 < 100)
```

**Statistics**:
```
Skipped 12 sequences: 8 too short, 4 invalid characters
```

**Solution**:
- Exclude short sequences from input TSV
- Sequences <100bp are typically too short for reliable COI analysis

---

#### Invalid Nucleotide Characters

**Checked in**: `utils.validate_sequence()` and `dereplication.dereplicate_from_fasta()`
**Error Type**: Data Quality (WARNING, sequence skipped)

**Valid Characters**: `ACGTN-` (IUPAC ambiguity codes accepted in dereplication)
**Extended in dereplication**: `ACGTURYKMSWBDHVN-`

**Log Message**:
```
WARNING: Skipping SAMPLE042: invalid characters {'X', '?', '#'}
```

**Common Invalid Characters**:
- `X`: Unknown amino acid (protein sequence)
- `?`: Missing data marker
- `#`: Comment character
- Numbers: Accession numbers in sequence
- Spaces: Formatting issues

**Solution**:
- Check that sequences are DNA, not protein
- Remove formatting characters
- Replace invalid ambiguity codes with 'N'

---

#### Excessive N Content

**Checked in**: `utils.validate_sequence()`
**Error Type**: Data Quality (WARNING, sequence skipped)

**Threshold**: 50% N bases

**Log Message**:
```
WARNING: Skipping SAMPLE133: Excessive ambiguous bases (65.3% N)
```

**Solution**:
- Exclude low-quality sequences from input
- Sequences with >50% N cannot be reliably aligned or clustered

---

#### Excessive Gaps

**Checked in**: `utils.validate_sequence()` (for aligned sequences)
**Error Type**: Data Quality (WARNING, sequence skipped)

**Threshold**: 80% gap characters (-)

**Log Message**:
```
WARNING: Skipping SAMPLE201: Excessive gaps (85.2%)
```

**Solution**:
- Typically only applies to pre-aligned sequences
- BOLDGenotyper performs alignment internally - input should be unaligned
- If providing aligned sequences, ensure gaps are from biological indels, not padding

---

### Geographic Coordinate Validation

#### Missing Coordinates

**Checked in**: `metadata.extract_coordinates()`, `geographic.create_points_geodataframe()`
**Error Type**: Data Quality (DEBUG/INFO, null geometry created)

**Behavior**:
- Missing lat/lon → null geometry
- Excluded from spatial joins
- Fallback to `country/ocean` column if available

**Log Message** (DEBUG):
```
DEBUG: Coordinate string is empty or null for row 42
```

**Statistics Logged** (INFO):
```
Geographic analysis results:
  Samples with coordinates: 1200/1485 (80.8%)
  Missing coordinates: 285
  Fallback to country/ocean: 150
  Final unknown: 135
```

**Solution**:
- Acceptable to have samples without coordinates
- Pipeline uses `country/ocean` column as fallback
- Remaining samples assigned `ocean_basin='Unknown'`

---

#### Out-of-Range Coordinates

**Checked in**: `metadata.extract_coordinates()`, `geographic.create_points_geodataframe()`
**Error Type**: Data Quality (WARNING/DEBUG, null geometry created)

**Valid Ranges**:
- Latitude: -90 to +90
- Longitude: -180 to +180

**Log Message** (geographic.py, WARNING):
```
WARNING: Latitude 127.5 out of range [-90, 90] for row 142, creating null geometry
WARNING: Longitude -195.3 out of range [-180, 180] for row 156, creating null geometry
```

**Common Causes**:
- Latitude/longitude swapped (lon, lat instead of lat, lon)
- Incorrect decimal placement (1275 instead of 12.75)
- Wrong sign (127.5 instead of -127.5)
- Coordinate system confusion (UTM vs lat/lon)

**Solution**:
- Check coordinate order in source data
- Verify decimal point placement
- Ensure coordinates are in WGS84 decimal degrees format

---

#### Coordinate Order Ambiguity

**Checked in**: `metadata.infer_coord_order()`
**Error Type**: Data Quality (WARNING if conflicting evidence)

**Behavior**:
- Examines first 100 coordinate pairs
- Detects if values are in lat/lon or lon/lat order
- Handles conflicting evidence by defaulting to most common pattern

**Log Message**:
```
WARNING: Conflicting coordinate order evidence detected:
  Evidence for lat,lon order: 45 samples
  Evidence for lon,lat order: 38 samples
  Defaulting to lat,lon based on majority
```

**Solution**:
- Review source data for coordinate order consistency
- Manually verify a few samples on a map
- Consider standardizing input data to single format

---

#### Centroid Coordinates

**Checked in**: `metadata.filter_by_coordinate_quality()`
**Error Type**: Data Quality (INFO, filtered out)

**Behavior**:
- Detects centroids by checking `coord_source` column for keyword 'centroid'
- Removes centroid coordinates (county/country centroids are imprecise)
- Logs count of removed samples

**Log Message**:
```
INFO: Coordinate filtering: 1450/1485 samples retained (97.6%), 35 excluded
INFO: Excluded coordinates: 35 centroids, 0 [0,0] points
```

**Solution**:
- Acceptable behavior - centroids are imprecise
- If you need to include centroid data, can skip this filtering step via config

---

#### [0, 0] Coordinates

**Checked in**: `metadata.filter_by_coordinate_quality()`
**Error Type**: Data Quality (INFO, filtered out)

**Behavior**:
- Removes coordinates at exactly [0, 0] (null island)
- Typically indicates missing data encoded as zeros

**Log Message**:
```
INFO: Excluded coordinates: 0 centroids, 12 [0,0] points
```

**Solution**:
- Acceptable behavior - [0, 0] is typically a data error
- Fix upstream data to use null/empty for missing coordinates

---

#### Terrestrial Coordinates

**Checked in**: `geographic.validate_marine_coordinates()` (if called)
**Error Type**: Data Quality (WARNING if >20% terrestrial)

**Behavior**:
- Uses Natural Earth land polygons to identify terrestrial samples
- Logs warning if >20% of samples appear terrestrial
- **Default**: No automatic filtering (keeps terrestrial samples)

**Log Message**:
```
WARNING: 245/1200 samples (20.4%) fall outside GOaS basins.
This may indicate terrestrial coordinates or incomplete basin coverage.
```

**Common Causes**:
- Freshwater collection sites
- Near-shore samples slightly on land due to resolution
- Data entry errors (coordinate swaps, decimal errors)

**Solution**:
- Review samples with `ocean_basin='Unknown'`
- Check if coordinates are truly terrestrial or just near-shore
- Consider using `classify_and_snap_coordinates()` to snap to nearest coast

---

### External Tool Validation

#### MAFFT Not Found

**Checked in**: `dereplication.run_mafft_alignment()`
**Error Type**: Critical (raises `AlignmentError`)

**Error Message**:
```
AlignmentError: MAFFT not found in PATH. Please install MAFFT:
  - macOS: brew install mafft
  - Ubuntu/Debian: sudo apt-get install mafft
  - conda: conda install -c bioconda mafft
```

**Solution**: Install MAFFT following one of the provided instructions

---

#### trimAl Not Found

**Checked in**: `dereplication.run_trimal_trimming()`
**Error Type**: Critical (raises `TrimmingError`)

**Error Message**:
```
TrimmingError: trimAl not found in PATH. Please install trimAl:
  - macOS: brew install trimal
  - Ubuntu/Debian: sudo apt-get install trimal
  - conda: conda install -c bioconda trimal
```

**Solution**: Install trimAl following one of the provided instructions

---

#### FastTree Not Found (for phylogeny)

**Checked in**: `core.run_pipeline()` (before phylogenetics phase)
**Error Type**: Recoverable (WARNING, phylogeny skipped)

**Behavior**:
- Checks for FastTree before attempting tree building
- Skips phylogeny if not found
- Logs warning and continues pipeline

**Log Message**:
```
WARNING:   ⚠ FastTree not found, skipping phylogeny
```

**Solution**:
- Install FastTree: `conda install -c bioconda fasttree`
- Or run pipeline with `--build-tree` disabled

---

#### GOaS Shapefile Not Found

**Checked in**: `geographic.load_goas_data()`, `core.run_pipeline()`
**Error Type**: Recoverable (WARNING in core, FileNotFoundError if called directly)

**Behavior in Pipeline**:
- Checks for shapefile before geographic analysis
- Sets all `ocean_basin='Unknown'` if missing
- Logs warning with installation instructions
- Continues pipeline

**Error Message** (if `load_goas_data()` called directly):
```
FileNotFoundError: GOaS shapefile not found: /path/to/goas_v01.shp

BOLDGenotyper requires the GOaS (Global Oceans and Seas) shapefile
for ocean basin assignment.

Option 1 - Automated download (recommended):
    python -m boldgenotyper.goas_downloader

Option 2 - Manual download:
    1. Download from: https://www.marineregions.org/download_file.php?name=World_Seas_IHO_v3.zip
    2. Extract to: /path/to/data/
    3. Rename the .shp file to: goas_v01.shp

Option 3 - Skip geographic analysis:
    Run BOLDGenotyper with the --no-geo flag
```

**Solution**: Follow one of the three options provided

---

### Parameter Validation

#### Invalid Clustering Threshold

**Checked in**: `core.run_pipeline()`, `dereplication.cluster_sequences()`
**Error Type**: Critical (raises `ValueError` or `ClusteringError`)

**Valid Range**: 0 < threshold < 1

**Error Message**:
```
ValueError: cluster_threshold must be between 0 and 1, got 1.5
```

**Common Values**:
- 0.01 = 99% identity (default, standard for COI species delimitation)
- 0.005 = 99.5% identity (more stringent, for cryptic species)
- 0.02 = 98% identity (more permissive, for divergent taxa)

**Solution**: Use a value between 0 and 1 representing maximum distance

---

#### Invalid Identity Threshold

**Checked in**: `core.run_pipeline()`
**Error Type**: Critical (raises `ValueError`)

**Valid Range**: 0 < min_identity ≤ 1

**Error Message**:
```
ValueError: min_identity must be between 0 and 1, got 0.0
```

**Common Values**:
- 0.90 = 90% identity (default)
- 0.95 = 95% identity (more stringent)
- 0.85 = 85% identity (more permissive)

**Solution**: Use a value between 0 and 1 representing minimum identity for assignment

---

#### Invalid Thread Count

**Checked in**: `core.run_pipeline()`
**Error Type**: Critical (raises `ValueError`)

**Valid Range**: threads ≥ 1

**Error Message**:
```
ValueError: threads must be >= 1, got 0
```

**Recommendation**: Use number of CPU cores minus 1 for optimal performance

**Solution**: Provide a positive integer

---

## Common Errors and Solutions

### Error: "No sequences found after filtering"

**Cause**: All sequences failed validation (length, characters, N content)

**Check**:
1. Are sequences DNA (not protein)?
2. Are sequences at least 100bp?
3. Do sequences have <50% N content?

**Solution**:
- Review sequence quality in source data
- Check for protein sequences (wrong molecule type)
- Verify sequences are COI barcode region (~650bp)

---

### Error: "MAFFT alignment failed"

**Causes**:
1. Sequences too divergent
2. Sequences contain invalid characters
3. MAFFT crashed (memory, disk space)

**Log Check**:
```
ERROR: MAFFT alignment failed:
<stderr output from MAFFT>
```

**Solutions**:
- Check MAFFT stderr output for specific error
- Verify sequences are DNA and properly formatted
- Ensure sufficient disk space and memory
- Try with fewer sequences to isolate issue

---

### Error: "Clustering failed"

**Causes**:
1. Empty distance matrix
2. All sequences identical or all maximally different
3. Memory issues with large datasets

**Solutions**:
- Check that alignment succeeded
- Verify sequences have some variation
- Reduce dataset size or increase RAM

---

### Warning: "Geographic analysis failed"

**Non-critical** - pipeline continues

**Common Causes**:
1. GOaS shapefile missing
2. Coordinate column missing
3. GeoPandas not installed
4. CRS/projection issues

**Result**: All samples assigned `ocean_basin='Unknown'`

**Solutions**:
- Install GOaS data: `python -m boldgenotyper.goas_downloader`
- Install GeoPandas: `pip install geopandas`
- Use `--no-geo` flag to skip geographic analysis

---

### Warning: "Phylogenetic analysis failed"

**Non-critical** - pipeline continues

**Common Causes**:
1. FastTree not installed
2. Too few sequences (need at least 3-4)
3. Alignment failed
4. Tree building timeout

**Result**: No tree file generated, pipeline continues

**Solutions**:
- Install FastTree: `conda install -c bioconda fasttree`
- Ensure at least 3-4 consensus sequences
- Check that alignment files exist in intermediate/

---

### Warning: "Visualization generation had errors"

**Non-critical** - pipeline continues

**Common Causes**:
1. Matplotlib/Cartopy import issues
2. No data to plot (empty dataframes)
3. Missing coordinate columns
4. Font rendering issues

**Result**: Some or all plots not generated, pipeline continues

**Solutions**:
- Install visualization dependencies: `pip install -e ".[all]"`
- Check that data exists for plotting
- Review logs for specific plot failures

---

## Error Messages Reference

### Error Format

BOLDGenotyper uses consistent error message formatting:

```
<LEVEL>: <COMPONENT>: <MESSAGE>
  <CONTEXT/DETAILS>
  <SOLUTION/ACTION>
```

**Example**:
```
ERROR: boldgenotyper.metadata: BOLD TSV is missing required columns: ['processid']
  Available columns: ['sampleid', 'species', 'genus', 'lat', 'lon', 'nuc', 'coord']

  Solution: Ensure your TSV file has a 'processid' column identifying each sample.
```

### Logging Levels

| Level | Usage | Visibility | Action Required |
|-------|-------|------------|-----------------|
| DEBUG | Detailed diagnostic info | Log file only | None (informational) |
| INFO | Progress updates, statistics | Console + log | None (informational) |
| WARNING | Recoverable issues, skipped data | Console + log | Review if unexpected |
| ERROR | Critical failures | Console + log | Fix and re-run |

### Progress Symbols

```
✓  - Success
⚠  - Warning (non-critical)
⊘  - Skipped/Disabled
✗  - Error (critical)
```

**Example Log Output**:
```
PHASE 1: Data Loading and Quality Control
  ✓ Loaded 1485 samples
  ✓ Retained 1450/1485 samples after filtering
  ⚠ GOaS shapefile not found, skipping geographic analysis

PHASE 2: Sequence Dereplication
  ✓ Created 1450 FASTA records
  ✓ Identified 127 consensus genotypes

PHASE 3: Genotype Assignment
  ✓ Assigned 1423/1450 samples (98.1%)
```

---

## Recovery Mechanisms

### Graceful Degradation

BOLDGenotyper implements graceful degradation for optional features:

1. **Geographic Analysis**:
   - Missing GOaS → Set all basins to 'Unknown'
   - Missing coordinates → Use country/ocean fallback
   - All fallback → Assign 'Unknown'

2. **Phylogenetic Analysis**:
   - Missing FastTree → Skip tree building
   - Tree build fails → Continue without tree
   - Relabeling fails → Use original tree

3. **Visualization**:
   - Missing Cartopy → Skip map plots
   - Individual plot fails → Skip that plot, continue others
   - All plots fail → Continue to reports

### Data Filtering vs. Failing

**Filtered (Continues)**:
- Duplicate processids (keeps first)
- Invalid coordinates (sets to null)
- Short sequences (skips)
- Invalid characters (skips)
- Excessive N/gaps (skips)
- Centroids (removes)

**Fails (Stops)**:
- Missing input file
- Missing required columns
- Empty input file
- All sequences invalid
- MAFFT/trimAl failures
- Invalid parameters

### Intermediate Files for Recovery

If pipeline fails mid-run, intermediate files are preserved in `output/intermediate/`:

```
output/
├── intermediate/
│   ├── 01_parsed_metadata.tsv
│   ├── 02_filtered_metadata.tsv
│   ├── organism.fasta
│   └── dereplication/
│       ├── aligned.fasta
│       ├── trimmed.fasta
│       └── ...
```

**Manual Recovery**:
1. Check which phase failed (last logged phase)
2. Examine intermediate files for that phase
3. Fix input data issue
4. Re-run pipeline (intermediate files will be overwritten)

**Future Enhancement**: Checkpoint-based resumption from last successful phase

---

## Developer Guidelines

### Adding New Validation

When adding new validation checks, follow this pattern:

```python
def validate_new_field(value, context=""):
    """
    Validate new field.

    Parameters
    ----------
    value : type
        Value to validate
    context : str
        Context for error messages (e.g., row number, sample ID)

    Returns
    -------
    is_valid : bool
        Whether value passes validation
    message : str
        Error message if invalid, "Valid" if valid
    """
    # Check 1
    if not condition1:
        return False, f"Reason for failure in {context}"

    # Check 2
    if not condition2:
        return False, f"Another reason in {context}"

    return True, "Valid"
```

**Usage**:
```python
is_valid, message = validate_new_field(value, context=f"sample {processid}")
if not is_valid:
    logger.warning(f"Skipping {processid}: {message}")
    continue
```

### Custom Exception Guidelines

When to create custom exceptions:

1. **Module-specific failures**: `AlignmentError`, `GOaSDataError`
2. **User-fixable issues**: Provide installation/configuration instructions
3. **Chained failures**: Use `raise NewError(...) from original_error`

**Pattern**:
```python
class ModuleError(Exception):
    """Base exception for module errors."""
    pass

class SpecificError(ModuleError):
    """Specific error type."""
    pass

try:
    # operation
except OriginalError as e:
    raise SpecificError(f"User-friendly message:\n{instructions}") from e
```

### Error Logging Best Practices

1. **Use appropriate levels**:
   - DEBUG: Diagnostic details (coordinate parsing details)
   - INFO: Progress, statistics (samples loaded, genotypes identified)
   - WARNING: Recoverable issues (missing optional data, skipped samples)
   - ERROR: Critical failures (missing files, invalid parameters)

2. **Include context**:
   ```python
   logger.warning(f"Invalid coordinate for sample {processid}: {coord}")
   # Not: logger.warning("Invalid coordinate")
   ```

3. **Log statistics**:
   ```python
   logger.info(f"Filtered {n_filtered}/{n_total} samples ({pct:.1f}%)")
   ```

4. **Use exc_info for exceptions**:
   ```python
   except Exception as e:
       logger.error(f"Operation failed: {e}", exc_info=True)
   ```

---

## Validation Checklist

Before submitting data to BOLDGenotyper, verify:

### TSV File
- [ ] File exists and is readable
- [ ] Contains `processid` column
- [ ] Contains `nuc` column with DNA sequences
- [ ] No duplicate `processid` values (or first occurrence is preferred)
- [ ] UTF-8 or latin-1 encoding

### Sequences
- [ ] DNA sequences (ACGTN), not protein
- [ ] Length ≥ 100bp (ideally 500-700bp for COI)
- [ ] <50% N content
- [ ] No formatting characters (#, ?, X, spaces)

### Coordinates (optional)
- [ ] Decimal degrees format (not DMS)
- [ ] WGS84 coordinate system
- [ ] Latitude: -90 to +90
- [ ] Longitude: -180 to +180
- [ ] Consistent order (lat,lon or lon,lat)
- [ ] No centroid coordinates (or acceptable if removed)

### Environment
- [ ] MAFFT installed and in PATH
- [ ] trimAl installed and in PATH
- [ ] FastTree installed (if building phylogeny)
- [ ] GeoPandas installed (if using geographic analysis)
- [ ] GOaS shapefile downloaded (if using geographic analysis)

---

## Support and Debugging

### Log Files

All pipeline runs create a log file: `<output_dir>/<organism>_pipeline.log`

**Contents**:
- All INFO, WARNING, ERROR messages
- Timestamps for each operation
- Full exception tracebacks
- Diagnostic statistics

**For Support Requests**, include:
1. Command run
2. Relevant section of log file
3. Input file format (first few rows)
4. Python and package versions

### Verbose Logging

For detailed debugging, run with DEBUG logging:

```python
from boldgenotyper import utils
utils.setup_logging(log_level='DEBUG')
```

Or via CLI (if supported):
```bash
boldgenotyper input.tsv --log-level DEBUG
```

### Common Debug Checks

1. **Check input file**:
   ```bash
   head -n 5 input.tsv
   wc -l input.tsv
   ```

2. **Check tool availability**:
   ```bash
   which mafft
   which trimal
   which fasttree
   ```

3. **Check Python environment**:
   ```bash
   python -c "import geopandas; print(geopandas.__version__)"
   python -c "import cartopy; print(cartopy.__version__)"
   ```

4. **Check intermediate files**:
   ```bash
   ls -lh output/intermediate/
   head -n 10 output/intermediate/01_parsed_metadata.tsv
   ```

---

## Summary

BOLDGenotyper's error handling is designed to be:

✅ **Robust**: Handles edge cases and recovers from non-critical failures
✅ **Informative**: Provides clear error messages with solutions
✅ **Transparent**: Logs all decisions and skipped data
✅ **Forgiving**: Continues processing when possible

**Key Principles**:
1. Validate early (fail fast on critical errors)
2. Log everything (transparent decision-making)
3. Recover gracefully (continue when possible)
4. Guide users (actionable error messages)

---

**Document Version**: 1.0
**Last Updated**: 2025-11-18
**Maintainer**: BOLDGenotyper Development Team
