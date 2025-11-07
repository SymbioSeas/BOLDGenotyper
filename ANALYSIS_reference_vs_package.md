# BOLDGenotyper: Reference Scripts vs Primary Package Analysis

**Date**: 2025-11-06
**Purpose**: Compare legacy reference scripts with primary package implementation to identify gaps, essential features, and recommendations for broad applicability.

---

## Current State

### Primary Package (`boldgenotyper/`)
- **Status**: Empty skeleton
- **Contents**: Single empty `__init__.py` file
- **Functionality**: None implemented

### Reference Scripts (`reference_scripts/`)
- **Status**: Fully functional but not modularized
- **Contents**: 5 standalone scripts implementing complete pipeline
- **Functionality**: All analysis capabilities present

---

## Detailed Script Analysis

### 1. `BOLD_tsv_to_fasta.py` (68 lines)

#### Current Implementation
**Purpose**: Convert BOLD TSV to FASTA format

**Calculations**:
- None (pure data transformation)

**Key Features**:
- Required columns: `record_id`, `species`, `nuc`
- Header format: `>{species}_{record_id}`
- Optional sequence wrapping (default: 80 chars)

**Assumptions**:
- TSV is well-formed with tab delimiters
- Species names can have spaces (replaced with underscores)
- Empty sequences are skipped

**Parameters**:
- `tsv`: Input file path (required)
- `--wrap`: Sequence wrapping width (default: 80, 0=no wrap)

**Issues for Broad Applicability**:
1. ❌ **Hard-coded column names**: Requires exactly `record_id`, `species`, `nuc`
   - BOLD downloads may have `processid`, `nucleotides`, `bin_uri`, etc.
   - Should support column mapping or auto-detection

2. ❌ **No organism name extraction**: Doesn't parse filename for organism
   - Critical for pipeline automation

3. ❌ **No metadata preservation**: Discards all other columns
   - Users may want lat/lon, country, etc. in downstream analysis

4. ❌ **No sequence validation**: Doesn't check for valid DNA characters
   - Could include invalid bases that break alignment

---

### 2. `msa_to_consensus.py` (158 lines)

#### Current Implementation
**Purpose**: Align sequences, trim, cluster, and generate consensus

**Calculations**:

1. **Distance Calculation** (lines 32-40):
```python
def compute_distance(seq1: str, seq2: str) -> float:
    matches = sites = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a in "ACGT" and b in "ACGT":
            sites += 1
            if a == b:
                matches += 1
    return 1.0 if sites == 0 else 1.0 - (matches / sites)
```
- **Formula**: `distance = 1 - (matches / informative_sites)`
- **Ignores**: Gaps (`-`), ambiguous bases (`N`), other IUPAC codes
- **Edge case**: Returns 1.0 if no informative sites

2. **Clustering** (lines 85-87):
```python
Z = linkage(dist_vec, method="average")
labels = fcluster(Z, t=threshold, criterion="distance")
```
- **Method**: Hierarchical clustering with average linkage (UPGMA)
- **Threshold**: Distance cutoff (default: 0.01)

3. **Consensus Generation** (lines 99-106):
```python
for i in range(L):
    col = [r.seq[i] for r in members]
    legit = [b.upper() for b in col if b.upper() in "ACGT"]
    if not legit:
        cons.append("N")
    else:
        base, cnt = Counter(legit).most_common(1)[0]
        cons.append(base if cnt/len(legit) >= freq_cutoff else "N")
```
- **Method**: Majority rule with frequency cutoff
- **Default freq_cutoff**: 0.7 (70%)
- **Ambiguous positions**: Called as `N`

**Parameters**:
- `indir`: Directory with FASTA files (required)
- `--threshold`: Distance cutoff (default: 0.01)
- `--freq-cutoff`: Minimum frequency for consensus base (default: 0.7)
- `--wrap`: Sequence wrapping (default: 80)

**Assumptions**:
1. MAFFT and trimAl available in PATH
2. All `.fa`, `.fasta`, `.fna` files in directory should be processed
3. Sufficient memory for distance matrix (O(n²) space)
4. Average linkage clustering is appropriate

**Issues for Broad Applicability**:

1. ✅ **Good**: Distance calculation ignores gaps/Ns (appropriate for COI)

2. ❌ **No IUPAC ambiguity support**: Only recognizes ACGT
   - Some sequences may have R, Y, W, S, etc.
   - Should optionally handle IUPAC codes

3. ❌ **Fixed clustering method**: Hard-coded to average linkage
   - Single linkage might be better for some datasets
   - Complete linkage might reduce over-clustering

4. ⚠️ **Memory scaling**: O(n²) distance matrix
   - Will fail on very large datasets (>10,000 sequences)
   - Consider sampling or progressive alignment for large datasets

5. ❌ **No minimum cluster size**: Small clusters (n=1, 2) may not be reliable
   - Could add `--min-cluster-size` parameter

6. ❌ **No gap handling in consensus**: Consensus ignores gaps entirely
   - If all sequences have gaps at a position, should it be gap or N?

7. ❌ **Batch processing only**: Processes entire directory
   - Should support single-file mode

8. ✅ **Good**: Creates organized subdirectories for outputs

---

### 3. `consensus_group_to_metadata.py` (304 lines)

#### Current Implementation
**Purpose**: Assign samples to consensus genotypes using sequence similarity

**Calculations**:

1. **ProcessID Extraction** (lines 91-102):
```python
PID_REGEX = re.compile(r"_(?P<pid>[^.\s_]+)(?:[.\s]|$)")
def extract_pid_from_header(header):
    m = PID_REGEX.search(header)
    if m:
        return m.group("pid")
    return None
```
- **Pattern**: Extracts last underscore-delimited field before dot/whitespace
- **Example**: `Sphyrna_lewini_ANGBF11456-15.COI-5P` → `ANGBF11456-15`

2. **Levenshtein Distance** (lines 105-119):
```python
def levenshtein_distance(a, b):
    if len(a) > len(b):
        a, b = b, a
    previous = list(range(len(a) + 1))
    for j, ch_b in enumerate(b, start=1):
        current = [j]
        for i, ch_a in enumerate(a, start=1):
            ins = previous[i] + 1
            dele = current[i - 1] + 1
            sub = previous[i - 1] + (ch_a != ch_b)
            current.append(min(ins, dele, sub))
        previous = current
    return previous[-1]
```
- **Method**: Dynamic programming implementation
- **Complexity**: O(m×n) time, O(min(m,n)) space
- **Operations**: Insertion, deletion, substitution (all cost 1)

3. **Identity Calculation** (lines 122-127):
```python
def identity_from_distance(d, len_a, len_b):
    denom = max(len_a, len_b)
    if denom == 0:
        return 1.0
    return 1.0 - (d / denom)
```
- **Formula**: `identity = 1 - (edit_distance / max_length)`
- **Normalization**: By longer sequence length

4. **Best Match Assignment** (lines 130-157):
```python
for gid, cseq in consensus_list:
    if use_edlib:
        res = edlib.align(seq, cseq, mode="NW", task="distance")
        d = res["editDistance"]
    else:
        d = levenshtein_distance(seq, cseq)
    iden = identity_from_distance(d, len(seq), len(cseq))

    if iden > best_iden:
        second_gid, second_iden = best_gid, best_iden
        best_gid, best_iden = gid, iden
    elif iden > second_iden:
        second_gid, second_iden = gid, iden

if best_iden < min_identity:
    return (None, best_iden, second_gid, second_iden)
```
- **Alignment**: Global (Needleman-Wunsch)
- **Keeps**: Best and runner-up matches
- **Threshold**: Minimum identity required (default: 0.90)

**Parameters**:
- `--metadata`: Input TSV with processid column (required)
- `--fasta`: Raw FASTA with processid in headers (required)
- `--consensus`: Consensus FASTA (required)
- `--out`: Output TSV path (required)
- `--min-identity`: Minimum identity for assignment (default: 0.90)
- `--threads`: Parallel processes (default: 1)
- `--diag`: Optional diagnostics CSV path

**Assumptions**:
1. ProcessID in FASTA header follows specific pattern
2. ProcessID in metadata matches exactly
3. Consensus sequences are representative of clusters
4. Global alignment appropriate (sequences full-length)

**Issues for Broad Applicability**:

1. ❌ **Rigid processID extraction**: Single regex pattern
   - BOLD formats vary: `>processid|gene`, `>bin|processid`, etc.
   - Should support multiple patterns or custom regex

2. ❌ **Global alignment only**: Needleman-Wunsch assumes full-length
   - Partial sequences or fragments need local alignment (Smith-Waterman)
   - COI has conserved regions - semi-global might be better

3. ⚠️ **Identity normalization**: Uses max(len1, len2)
   - Alternative: alignment length, min length, average length
   - Different normalizations affect threshold interpretation

4. ❌ **No gap penalties**: Edit distance treats gaps same as mismatches
   - COI may have real indels (rare but possible)
   - Biological alignment should weight gaps differently

5. ❌ **Hard assignment only**: Assigns to single best group
   - Ambiguous samples (similar to multiple groups) not flagged
   - Could add "ambiguous" category if best and 2nd too close

6. ✅ **Good**: Parallel processing support

7. ✅ **Good**: Diagnostics output for QC

8. ❌ **No unaligned sequence handling**: Assumes all same gene region
   - If mixing markers or gene regions, should validate

---

### 4. `plot_shark_genotypes_ocean_basins_complex.py` (763 lines)

#### Current Implementation
**Purpose**: Visualize genotype distributions with geographic and ocean basin context

**Key Calculations**:

1. **Coordinate Parsing** (lines 47-56):
```python
def parse_coord_raw(val: str) -> Tuple[float, float]:
    if pd.isna(val):
        return float("nan"), float("nan")
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(val))
    if len(nums) >= 2:
        try:
            return float(nums[0]), float(nums[1])
        except ValueError:
            return float("nan"), float("nan")
    return float("nan"), float("nan")
```
- **Method**: Regex extraction of first two numbers
- **Handles**: Various formats (comma-separated, space-separated, etc.)

2. **Coordinate Order Inference** (lines 65-79):
```python
def infer_coord_order(a: pd.Series, b: pd.Series) -> Literal["lonlat", "latlon"]:
    valid = ~(a.isna() | b.isna())
    if valid.sum() == 0:
        return "lonlat"
    aa = a[valid].abs()
    bb = b[valid].abs()
    cond_latlon = (aa <= 90).mean() > 0.7 and (bb <= 180).mean() > 0.7
    cond_lonlat = (aa <= 180).mean() > 0.7 and (bb <= 90).mean() > 0.7
    # ... more logic
```
- **Method**: Statistical inference based on value ranges
- **Latitude**: -90 to 90
- **Longitude**: -180 to 180

3. **Land Detection** (lines 104-114):
```python
def is_on_land(row):
    pt = Point(float(row["lon"]), float(row["lat"]))
    try:
        return land_union.contains(pt)
    except Exception:
        return False
```
- **Method**: Point-in-polygon test using Natural Earth land geometries
- **Resolution**: 110m (coarse but fast)

4. **Ocean Basin Assignment** (lines 211-225):
```python
def get_basin(lon, lat):
    if np.isnan(lon) or np.isnan(lat):
        return np.nan
    p = Point(float(lon), float(lat))
    for prepped, geom, name in basin_recs:
        try:
            if prepped.contains(p) or geom.contains(p):
                return name
        except Exception:
            continue
    return np.nan
```
- **Method**: Point-in-polygon with prepared geometries for speed
- **Source**: External shapefile (user-provided)

5. **Territory Polygon Generation** (lines 289-305):
```python
def territory_from_points(points_lonlat: np.ndarray,
                          method: Literal["buffer", "convex"] = "buffer",
                          buffer_deg: float = 0.5):
    pts = [Point(float(lon), float(lat)) for lon, lat in points_lonlat]
    if not pts:
        return None
    if method == "convex":
        if len(pts) < 3:
            return MultiPoint(pts).buffer(max(1e-3, buffer_deg * 0.5))
        return MultiPoint(pts).convex_hull
    # buffer method
    return unary_union([p.buffer(buffer_deg) for p in pts])
```
- **Buffer method**: Union of circular buffers (smooth blobs)
- **Convex hull**: Minimum convex polygon containing all points
- **Default buffer**: 0.5 degrees (≈55 km at equator)

6. **Bubble Size Calculation** (line 539):
```python
sub["marker_size"] = size_min + size_scale * np.sqrt(sub["n_at_coord"].astype(float))
```
- **Formula**: `size = size_min + size_scale × √count`
- **Rationale**: Square root for perceptual linearity

**Parameters**: 29 command-line arguments including:
- `--csv`: Input data (required)
- `--outdir`: Output directory (required)
- `--plot-mode`: "simple" or "basins"
- `--basin-shp`: Shapefile path for ocean basins
- `--alpha`, `--size-min`, `--size-scale`: Visual aesthetics
- `--land-policy`: "mark", "drop", "snap"
- `--draw-territories`: Enable territory polygons
- `--territory-*`: 9 parameters for territory customization
- `--legend-*`: 5 parameters for legend customization
- Color, coordinate, and styling options

**Assumptions**:
1. Coordinates in `coord` column or separate `lat`/`lon`
2. `consensus_group` column exists
3. `COI_common` column for species names
4. Cartopy and shapely available
5. Sufficient memory for territory polygons

**Issues for Broad Applicability**:

1. ❌ **Hard-coded column names**: Requires `coord`, `consensus_group`, `COI_common`
   - Should support column mapping configuration
   - Different BOLD downloads have different column names

2. ⚠️ **Territory calculation in geographic coordinates**: Buffers in degrees
   - Distortion at high latitudes (1° lon ≠ 1° lat distance)
   - Should use proper projected coordinate systems for accuracy
   - Could project to equal-area CRS, buffer, then reproject

3. ❌ **No coordinate validation**: Doesn't check for valid ranges
   - Lat outside [-90, 90] or lon outside [-180, 180] should error
   - Invalid coords can break projections

4. ❌ **Ocean basin shapefile required**: Not bundled
   - Users must provide their own shapefile
   - Should include default basins or download automatically

5. ❌ **Shark-specific naming**: Script name and some variables reference sharks
   - Should be renamed for general use

6. ⚠️ **Complexity**: 763 lines, many interdependent functions
   - Difficult to maintain and extend
   - Should be refactored into classes/modules

7. ✅ **Good**: Extensive customization options

8. ✅ **Good**: Handles multiple location types (precise, centroid, snapped)

9. ✅ **Good**: Creates both simple and basin visualizations

10. ❌ **No temporal analysis**: If data has dates, could show change over time
    - Useful for tracking populations, invasive species, etc.

---

### 5. `genotype_basin_separation.Rmd` (210 lines)

#### Current Implementation
**Purpose**: Statistical analysis of genotype × ocean basin associations

**Key Calculations**:

1. **Bias-Corrected Cramér's V** (lines 26-37):
```r
cramers_v_bc <- function(chi2, n, r, k){
  phi2 <- chi2 / n
  r_corr <- r - 1
  k_corr <- k - 1
  if(r_corr == 0 || k_corr == 0) return(NA_real_)
  phi2corr <- max(0, phi2 - (r_corr * k_corr) / (n - 1))
  rcorr <- r_corr - ((r_corr^2 - 1) / (n - 1))
  kcorr <- k_corr - ((k_corr^2 - 1) / (n - 1))
  denom <- min(rcorr, kcorr)
  if(denom <= 0) return(NA_real_)
  sqrt(phi2corr / denom)
}
```
- **Method**: Bergsma (2013) bias correction for small samples
- **Advantages**: Better than standard Cramér's V for small n
- **Range**: 0 (no association) to 1 (perfect association)

2. **Chi-Square Test** (line 104):
```r
ct <- suppressWarnings(chisq.test(mat, correct = FALSE))
```
- **Method**: Pearson's chi-square test of independence
- **H₀**: Genotype and ocean basin are independent
- **No continuity correction**: Yates correction disabled

3. **Adjusted Standardized Residuals** (lines 114-116):
```r
row_prob <- rowSums(mat) / n
col_prob <- colSums(mat) / n
std_resid <- (mat - expected) / sqrt(expected * (1 - row_prob) %o% (1 - col_prob))
```
- **Method**: Haberman's adjusted residuals
- **Interpretation**: |residual| > 2 indicates significant over/under-representation
- **Advantage**: Accounts for marginal probabilities

**Outputs**:
- Contingency table (counts)
- Row-normalized proportions
- Chi-square test results
- Cramér's V effect size
- Standardized residuals
- Stacked bar charts (proportions)
- Residual heatmap

**Assumptions**:
1. Sufficient sample size for chi-square (expected frequencies ≥ 5)
2. Samples are independent
3. Genotype and basin are categorical variables
4. R packages installed: readr, dplyr, tidyr, ggplot2, ggthemes, scales, svglite

**Issues for Broad Applicability**:

1. ❌ **R-specific**: Requires R and multiple packages
   - Many Python users may not have R configured
   - Should implement Python equivalents using scipy, statsmodels

2. ❌ **Hard-coded basin order**: Lines 40-47 specify exact basin names
   - Won't work for other basins or organisms
   - Should auto-detect from data

3. ❌ **Parametrized file paths**: Lines 6-8 use absolute paths
   - Not portable
   - Should use relative paths or command-line arguments

4. ⚠️ **No sample size checks**: Doesn't warn if chi-square assumptions violated
   - Small expected frequencies → invalid p-values
   - Should check and recommend Fisher's exact test if needed

5. ❌ **No permutation tests**: Parametric assumptions may not hold
   - Consider permutation-based p-values for robustness

6. ⚠️ **Binary only**: Tests genotype × basin association
   - Could also test: genotype × depth, genotype × temperature, etc.
   - Should generalize to arbitrary categorical variables

7. ❌ **No temporal analysis**: Doesn't incorporate collection dates
   - Could test if genotype distributions changed over time

8. ✅ **Good**: Bias-corrected effect size

9. ✅ **Good**: Post-hoc residual analysis

10. ✅ **Good**: Publication-ready visualizations

---

## Major Gaps in Current Implementation

### 1. **No Primary Package Implementation**
- All functionality in standalone scripts
- No importable modules
- No unified API
- No pip/conda installable package

### 2. **Hard-coded Assumptions**
- Column names not configurable
- ProcessID regex pattern fixed
- File naming conventions rigid
- Many parameters not exposed

### 3. **Limited Input Format Support**
- Only works with specific BOLD TSV format
- No support for GenBank, FASTA with metadata, etc.
- No column mapping or auto-detection

### 4. **No Sequence Quality Control**
- No length filtering (too short/long)
- No stop codon detection
- No frame-shift detection
- No primer trimming
- No quality score filtering (if available)

### 5. **Limited Statistical Analysis**
- Only chi-square test implemented
- No phylogenetic analysis
- No population genetics metrics (FST, AMOVA, etc.)
- No diversity indices (π, h, etc.)

### 6. **No Validation or Testing**
- No unit tests
- No integration tests
- No example datasets
- No benchmarks

### 7. **No Documentation Beyond README**
- No API documentation
- No tutorials
- No troubleshooting guides
- No parameter guidelines

---

## Essential Modules for Primary Package

### Core Modules (Critical)

#### 1. `boldgenotyper.io` - Input/Output
**Functions**:
```python
def read_bold_tsv(path, column_map=None, validate=True)
    """Read BOLD TSV with flexible column mapping"""

def write_fasta(sequences, path, wrap=80, header_format=">{id}")
    """Write sequences to FASTA with customizable headers"""

def read_fasta(path, parse_header=True, header_parser=None)
    """Read FASTA with optional header parsing"""

def read_metadata(path, format='tsv', required_columns=None)
    """Read metadata from TSV/CSV with validation"""
```

**Key Features**:
- Flexible column mapping (user can specify their column names)
- Auto-detection of common column name variants
- Validation of required data
- Support for multiple input formats

---

#### 2. `boldgenotyper.sequence` - Sequence Operations
**Functions**:
```python
def validate_sequence(seq, alphabet='DNA', allow_ambiguous=True)
    """Check if sequence contains valid characters"""

def compute_pairwise_distance(seq1, seq2, method='identity',
                               ignore_gaps=True, ignore_ambiguous=True)
    """Calculate distance between two sequences"""

def align_sequences(sequences, method='mafft', options=None)
    """Align sequences using external tool"""

def trim_alignment(alignment, method='trimal', options=None)
    """Trim alignment using external tool"""

def quality_filter(sequences, min_length=None, max_length=None,
                   check_stopcodons=False, check_frameshifts=False)
    """Filter sequences by quality criteria"""
```

**Key Features**:
- Multiple distance metrics (identity, Hamming, K2P, JC69)
- Flexible gap/ambiguity handling
- Tool wrappers with error handling
- Quality control functions

---

#### 3. `boldgenotyper.clustering` - Genotype Identification
**Functions**:
```python
def cluster_sequences(alignment, threshold=0.01, method='average',
                     distance_metric='identity')
    """Cluster sequences by similarity"""

def generate_consensus(cluster_members, freq_cutoff=0.7,
                      handle_gaps='ignore', ambiguity_codes=True)
    """Generate consensus sequence from cluster members"""

def assign_to_clusters(sequences, consensus_seqs, min_identity=0.9,
                      alignment_mode='global', return_diagnostics=True)
    """Assign sequences to existing consensus groups"""
```

**Key Features**:
- Multiple clustering algorithms (average, single, complete linkage)
- Configurable distance metrics
- Flexible consensus calling
- Diagnostic outputs

---

#### 4. `boldgenotyper.geography` - Geographic Analysis
**Functions**:
```python
def parse_coordinates(coord_string, order='auto')
    """Parse coordinate strings with flexible formats"""

def validate_coordinates(lat, lon, strict=True)
    """Check if coordinates are valid"""

def filter_coordinates(df, policy='precise', land_union=None)
    """Filter samples by coordinate precision"""

def assign_ocean_basins(df, basin_shapefile=None, use_default=True)
    """Assign samples to ocean basins"""

def calculate_distances(coords1, coords2, method='haversine')
    """Calculate geographic distances"""
```

**Key Features**:
- Flexible coordinate parsing
- Multiple filtering policies
- Bundled default basin shapefile
- Distance calculations

---

#### 5. `boldgenotyper.visualization` - Plotting
**Functions**:
```python
def plot_distribution_map(df, color_by='genotype', size_by='count',
                         projection='platecarree', add_basins=False)
    """Create geographic distribution map"""

def plot_genotype_abundance(df, group_by='basin', plot_type='stacked_bar')
    """Plot genotype abundance by group"""

def plot_alignment(sequences, highlight_polymorphic=True)
    """Visualize sequence alignment"""

def plot_distance_matrix(distance_matrix, labels=None, method='heatmap')
    """Visualize distance/similarity matrix"""
```

**Key Features**:
- Multiple projection systems
- Customizable aesthetics
- Interactive plots (optional)
- Publication-ready defaults

---

#### 6. `boldgenotyper.statistics` - Statistical Analysis
**Functions**:
```python
def chi_square_test(contingency_table, correct=False)
    """Chi-square test of independence"""

def cramers_v(chi2, n, r, c, bias_corrected=True)
    """Effect size for categorical association"""

def standardized_residuals(contingency_table, method='haberman')
    """Post-hoc residual analysis"""

def diversity_indices(genotype_counts, measures=['shannon', 'simpson'])
    """Calculate diversity indices"""

def fst_analysis(genotypes_by_population)
    """Calculate FST between populations"""
```

**Key Features**:
- Comprehensive test suite
- Multiple effect size measures
- Population genetics metrics
- Permutation tests for robustness

---

### Utility Modules (Important)

#### 7. `boldgenotyper.utils` - Helper Functions
```python
def extract_organism_from_filename(filename, pattern=None)
    """Parse organism name from filename"""

def check_external_tools(required=['mafft', 'trimal'])
    """Verify external tools are available"""

def setup_logger(name, level='INFO', log_file=None)
    """Configure logging"""

def validate_config(config, schema=None)
    """Validate configuration parameters"""
```

---

#### 8. `boldgenotyper.config` - Configuration Management
```python
class Config:
    """Configuration object for pipeline parameters"""

def load_config(path)
    """Load configuration from YAML/JSON file"""

def get_default_config()
    """Return default pipeline configuration"""
```

---

#### 9. `boldgenotyper.pipeline` - Workflow Integration
```python
class Pipeline:
    """Complete analysis pipeline"""

    def run_full_pipeline(input_tsv, output_dir, config=None)
        """Execute all analysis steps"""

    def run_step(step_name, inputs, params=None)
        """Run individual pipeline step"""
```

---

## Critical Parameters for Broad Applicability

### 1. **Sequence Quality Control**
```python
# Essential parameters missing from current implementation
MIN_SEQUENCE_LENGTH = 400       # Minimum COI fragment length
MAX_SEQUENCE_LENGTH = 1000      # Maximum expected length
CHECK_STOP_CODONS = True        # Detect premature stops
CHECK_FRAME_SHIFTS = False      # Detect indels (rare in COI)
MAX_N_CONTENT = 0.05           # Maximum 5% ambiguous bases
TRIM_PRIMERS = False           # Remove primer sequences
```

**Rationale**: Poor quality sequences can create false genotypes

---

### 2. **Flexible Clustering**
```python
# Additional parameters needed
CLUSTERING_METHOD = 'average'    # average|single|complete|ward
DISTANCE_METRIC = 'identity'     # identity|k2p|jc69|hamming
MIN_CLUSTER_SIZE = 3            # Minimum samples per genotype
MERGE_RARE_GENOTYPES = False    # Combine rare genotypes
RARE_THRESHOLD = 0.01           # < 1% of samples = rare
```

**Rationale**: Different organisms/markers need different thresholds

---

### 3. **Assignment Flexibility**
```python
# More nuanced assignment
MIN_IDENTITY = 0.90             # Minimum for assignment
AMBIGUITY_THRESHOLD = 0.05      # If best and 2nd within 5% = ambiguous
ALIGNMENT_MODE = 'global'       # global|local|semiglobal
ALLOW_MULTIPLE_ASSIGNMENT = False  # Assign to multiple groups if close
```

**Rationale**: Some samples may be ambiguous or chimeric

---

### 4. **Geographic Precision**
```python
# Coordinate filtering
COORDINATE_PRECISION = 'precise'  # precise|country|region|any
MIN_COORDINATE_UNCERTAINTY = 10   # Maximum uncertainty in km
SNAP_TO_COAST_DISTANCE = 50       # Maximum distance for snapping (km)
REQUIRE_DEPTH_DATA = False        # Require bathymetry
```

**Rationale**: Different studies need different precision levels

---

### 5. **Metadata Flexibility**
```python
# Column mapping (user-configurable)
COLUMN_MAP = {
    'processid': ['processid', 'sample_id', 'id', 'record_id'],
    'sequence': ['nucleotides', 'nuc', 'seq', 'dna'],
    'lat': ['lat', 'latitude', 'dec_lat'],
    'lon': ['lon', 'long', 'longitude', 'dec_lon'],
    'species': ['species_name', 'species', 'COI_species'],
    'common_name': ['COI_common', 'common_name', 'vernacular'],
}
```

**Rationale**: BOLD downloads have varying column names

---

## Recommended Implementation Priority

### Phase 1: Core Functionality (Weeks 1-2)
1. ✅ `boldgenotyper.io` - Read/write with flexible formats
2. ✅ `boldgenotyper.sequence` - Distance calculation, alignment wrappers
3. ✅ `boldgenotyper.clustering` - Clustering and consensus
4. ✅ `boldgenotyper.utils` - Helper functions

**Goal**: Reproduce reference script functionality in modular form

---

### Phase 2: Enhanced Analysis (Weeks 3-4)
5. ✅ `boldgenotyper.geography` - Coordinate handling, basin assignment
6. ✅ `boldgenotyper.statistics` - Statistical tests
7. ✅ `boldgenotyper.quality` - New QC module
8. ✅ Unit tests for all modules

**Goal**: Add missing QC and statistical features

---

### Phase 3: Visualization & Integration (Weeks 5-6)
9. ✅ `boldgenotyper.visualization` - Refactored plotting
10. ✅ `boldgenotyper.pipeline` - Complete workflow
11. ✅ `boldgenotyper.config` - Configuration system
12. ✅ Command-line interface
13. ✅ Integration tests

**Goal**: User-friendly interface and complete pipeline

---

### Phase 4: Polish & Documentation (Weeks 7-8)
14. ✅ Comprehensive documentation
15. ✅ Tutorial notebooks
16. ✅ Example datasets
17. ✅ Performance optimization
18. ✅ Package for PyPI/conda

**Goal**: Production-ready, publishable package

---

## Calculation Differences Summary

| Feature | Reference Scripts | Should Implement |
|---------|------------------|------------------|
| **Distance metric** | Identity only | Identity, K2P, JC69 |
| **Gap handling** | Ignore | Ignore, penalize, or treat as 5th base |
| **IUPAC codes** | Only ACGT | Support full IUPAC alphabet |
| **Clustering method** | Average linkage | Average, single, complete, Ward |
| **Consensus calling** | Majority rule | Majority, plurality, threshold |
| **Alignment mode** | Global (NW) | Global, local (SW), semi-global |
| **Identity normalization** | Max length | Max, min, alignment, query length |
| **Coordinate projection** | Geographic (degrees) | Geographic, projected (meters) |
| **Statistical tests** | Chi-square only | Chi-square, Fisher, permutation |
| **Effect sizes** | Cramér's V | Cramér's V, odds ratios, φ |

---

## Critical Missing Features

### 1. **Sequence Quality Control** (HIGH PRIORITY)
- Length filtering
- Stop codon detection
- Ambiguous base content
- Frame-shift detection
- Primer trimming

**Why**: Bad sequences create false genotypes → invalid conclusions

---

### 2. **Phylogenetic Analysis** (MEDIUM PRIORITY)
- Tree building (ML, NJ, UPGMA)
- Bootstrap support
- Tree visualization
- Root finding
- Clade assignment

**Why**: Important for understanding evolutionary relationships

---

### 3. **Population Genetics** (MEDIUM PRIORITY)
- FST calculations
- AMOVA
- Tajima's D
- Haplotype diversity
- Gene flow estimates

**Why**: Standard analyses for population structure studies

---

### 4. **Temporal Analysis** (LOW-MEDIUM PRIORITY)
- Change over time
- Seasonal patterns
- Invasion dynamics
- Demographic history

**Why**: Many studies have temporal component

---

### 5. **Batch Processing** (HIGH PRIORITY)
- Multiple organisms at once
- Comparative analysis
- Meta-barcoding support
- HPC cluster submission

**Why**: Real-world use often involves many species

---

### 6. **Output Standardization** (MEDIUM PRIORITY)
- Darwin Core format
- GenBank submission format
- GBIF format
- Standard metadata schemas

**Why**: Interoperability with other databases

---

### 7. **Web Interface** (LOW PRIORITY - Future)
- Upload data
- Configure parameters
- Run pipeline
- Download results
- Share analyses

**Why**: Accessibility for non-programmers

---

## Recommendations

### Immediate Actions (Next 2 Weeks)

1. **Create core package structure**:
```
boldgenotyper/
├── __init__.py
├── io.py           # I/O functions
├── sequence.py     # Sequence operations
├── clustering.py   # Clustering and consensus
├── geography.py    # Geographic functions
├── statistics.py   # Statistical tests
├── visualization.py # Plotting
├── quality.py      # Quality control (NEW)
├── utils.py        # Helper functions
├── config.py       # Configuration
└── pipeline.py     # Workflow integration
```

2. **Implement flexible I/O**:
   - Column name mapping
   - Multiple input formats
   - Validation with helpful errors

3. **Add quality control**:
   - Length filtering
   - Stop codon detection
   - Ambiguous base filtering

4. **Refactor clustering**:
   - Multiple distance metrics
   - Multiple clustering methods
   - Minimum cluster size parameter

5. **Write comprehensive tests**:
   - Unit tests for each module
   - Integration tests for pipeline
   - Example datasets for testing

### Short-term (Month 1-2)

6. **Enhance statistical analysis**:
   - Implement in Python (remove R dependency)
   - Add permutation tests
   - Add diversity indices

7. **Improve geographic analysis**:
   - Bundle default ocean basin shapefile
   - Support custom regions
   - Proper projected coordinate systems

8. **Create command-line interface**:
   - Simple one-command pipeline
   - Individual commands for each step
   - Configuration file support

### Medium-term (Month 3-6)

9. **Add phylogenetic analysis**:
   - Tree building wrappers
   - Tree visualization
   - Clade-based analysis

10. **Add population genetics**:
    - FST and AMOVA
    - Diversity indices
    - Gene flow estimation

11. **Documentation**:
    - API documentation (Sphinx)
    - Tutorial notebooks
    - Best practices guide

12. **Package distribution**:
    - PyPI package
    - Conda package
    - Docker container

### Long-term (Month 6+)

13. **Web interface** (Streamlit or Django)
14. **Database integration** (local caching)
15. **Cloud deployment** options
16. **GPU acceleration** for large datasets

---

## Conclusion

The reference scripts provide a **solid foundation** but have significant **limitations for broad applicability**:

### Strengths
✅ Core algorithms are sound
✅ Good default parameters for COI
✅ Comprehensive visualization
✅ Statistical rigor (bias-corrected effect sizes)

### Critical Gaps
❌ No modular package structure
❌ Hard-coded assumptions
❌ No sequence quality control
❌ Limited statistical analysis
❌ No phylogenetic capabilities
❌ Rigid input format requirements

### Priority Recommendations
1. **Refactor into proper package** with modular design
2. **Add quality control** to prevent bad genotypes
3. **Flexible I/O** with column mapping
4. **Multiple distance/clustering methods**
5. **Comprehensive testing** with example data
6. **Remove R dependency** (Python implementation of statistics)

This will create a **robust, flexible, well-tested package** suitable for any barcoding study across diverse organisms.
