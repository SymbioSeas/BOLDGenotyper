# BOLDGenotyper

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://zenodo.org/)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/SymbioSeas/BOLDGenotyper/releases)

**Automated COI sequence genotyping and biogeographic analysis from BOLD database data**

BOLDGenotyper is a bioinformatics pipeline that enables researchers to identify and analyze COI (Cytochrome Oxidase I) genotypes from the BOLD (Barcode of Life Database) for any taxonomic group. The pipeline performs sequence dereplication, genotype assignment, geographic filtering, ocean basin classification, and visualization of biogeographic patterns.

This package enables reproducible analysis of mitochondrial COI genotypes and their geographic distributions, as demonstrated in the companion manuscript analyzing *Sphyrna lewini* (scalloped hammerhead shark) genotypes separated by ocean basin.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Analysis Pipeline](#data-analysis-pipeline)
- [Processing Assumptions](#processing-assumptions)
- [Usage Guide](#usage-guide)
- [Output Files](#output-files)
- [Advanced Options](#advanced-options)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Version History](#version-history)

---

## Features

- **Automated genotyping**: Identifies unique COI genotypes through sequence clustering and consensus generation
- **Geographic analysis**: Filters coordinates and assigns samples to ocean basins (marine organisms only in v1.0)
- **Visualization**: Creates publication-ready maps showing genotype distributions with territory polygons
- **Statistical analysis**: Chi-square tests and standardized residuals for genotype × ocean basin associations
- **Reproducibility**: Standardized workflow for any organism in BOLD database
- **Flexibility**: Customizable parameters for clustering, filtering, and visualization
- **Modular design**: Skip geographic analysis with `--no-geo` flag for non-marine organisms or when GOaS is unavailable

---

## Installation

### Prerequisites

BOLDGenotyper requires:
- **Python**: ≥3.8
- **Bioinformatics tools**: MAFFT v7, trimAl
- **Conda** (recommended for dependency management)

### Step 1: Clone the Repository

```bash
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper
```

### Step 2: Create Conda Environment

```bash
# Create environment from the provided YAML file
conda env create -f boldgenotyper_env.yml

# Activate the environment
conda activate depredation
```

The environment includes all required dependencies:
- **Python packages**: biopython, pandas, numpy, scipy, matplotlib, cartopy, edlib
- **External tools**: MAFFT (v7.526), trimAl (v1.5.0)

### Step 3: Verify Installation

Check that required tools are available:

```bash
# Check MAFFT
mafft --version

# Check trimAl
trimal --version

# Check Python packages
python -c "import Bio, pandas, numpy, scipy, matplotlib, cartopy, edlib; print('All packages installed successfully')"
```

---

## GOaS Shapefile Setup (Required for Geographic Analysis)

### What is GOaS?

GOaS (Global Oceans and Seas) is a shapefile dataset that defines standardized ocean basin boundaries. BOLDGenotyper uses this dataset to:
- Assign samples to ocean basins (e.g., North Atlantic, South Pacific)
- Create geographic distribution maps
- Analyze genotype-by-basin patterns

**Important Notes:**
- The GOaS shapefile (~100-200 MB) is **not included** in this repository due to size constraints
- Geographic analysis is **currently designed for marine organisms only**
- If you're only interested in genotyping and phylogeny (not geography), you can skip this setup using the `--no-geo` flag

### Downloading and Installing GOaS

Choose one of the following options:

#### Option 1: Automated Download (Recommended)

```bash
# Navigate to the BOLDGenotyper directory
cd BOLDGenotyper

# Activate your conda environment
conda activate depredation

# Run the GOaS downloader
python -m boldgenotyper.goas_downloader

# This will:
# 1. Download World_Seas_IHO_v3.zip from Marine Regions
# 2. Extract it to boldgenotyper/GOaS_v1_20211214/
# 3. Verify all required files are present
```

#### Option 2: Manual Download

If the automated download fails or you prefer manual setup:

1. **Download the shapefile:**
   - Visit: https://www.marineregions.org/download_file.php?name=World_Seas_IHO_v3.zip
   - Or alternative: https://github.com/iobis/mregions-static/raw/master/shapefiles/World_Seas_IHO_v3.zip

2. **Extract the ZIP file:**
   ```bash
   unzip World_Seas_IHO_v3.zip
   ```

3. **Move files to the correct location:**
   ```bash
   # Create the directory if it doesn't exist
   mkdir -p boldgenotyper/GOaS_v1_20211214

   # Move all files (must include .shp, .shx, .dbf, .prj, .cpg)
   mv World_Seas_IHO_v3.* boldgenotyper/GOaS_v1_20211214/

   # Rename the main shapefile
   cd boldgenotyper/GOaS_v1_20211214/
   ln -s World_Seas_IHO_v3.shp goas_v01.shp
   ln -s World_Seas_IHO_v3.shx goas_v01.shx
   ln -s World_Seas_IHO_v3.dbf goas_v01.dbf
   ln -s World_Seas_IHO_v3.prj goas_v01.prj
   ln -s World_Seas_IHO_v3.cpg goas_v01.cpg
   ```

4. **Verify the installation:**
   ```bash
   # Check that the files exist
   ls -lh boldgenotyper/GOaS_v1_20211214/

   # You should see:
   # goas_v01.shp (or World_Seas_IHO_v3.shp)
   # goas_v01.shx (or World_Seas_IHO_v3.shx)
   # goas_v01.dbf (or World_Seas_IHO_v3.dbf)
   # goas_v01.prj (or World_Seas_IHO_v3.prj)
   # goas_v01.cpg (or World_Seas_IHO_v3.cpg)
   ```

#### Option 3: Skip Geographic Analysis

If you don't need geographic distribution analysis, you can run BOLDGenotyper without GOaS:

```bash
# Add the --no-geo flag to skip geographic analysis
boldgenotyper data/Euprymna_scolopes.tsv --no-geo

# This will:
# - Perform sequence clustering and genotyping
# - Generate phylogenetic trees (if --build-tree is specified)
# - Skip ocean basin assignment
# - Skip geographic visualizations
```

### Geographic Limitations

**Current Version (v1.0):**
- BOLDGenotyper is currently designed for **marine organisms only**
- Ocean basin assignment uses marine-specific GOaS boundaries
- Terrestrial and freshwater organisms are not yet supported for geographic analysis

**Future Development:**
- We plan to add modules for terrestrial environments
- Additional GOaS packages will support freshwater and terrestrial distributions
- Geographic analysis will be extended to non-marine taxa

**For Non-Marine Organisms:**
- Use the `--no-geo` flag to skip geographic modules
- You can still perform genotyping and phylogenetic analysis
- Geographic visualizations will be skipped

### Verifying GOaS Setup

To verify that GOaS is correctly installed, run a quick test:

```bash
python -c "from boldgenotyper import config, geographic; \
cfg = config.get_default_config(); \
print(f'GOaS path: {cfg.geographic.goas_shapefile_path}'); \
print(f'GOaS exists: {cfg.geographic.goas_shapefile_path.exists()}'); \
if cfg.geographic.goas_shapefile_path.exists(): \
    goas = geographic.load_goas_data(cfg.geographic.goas_shapefile_path); \
    print(f'Loaded {len(goas)} ocean basins')"
```

Expected output:
```
GOaS path: /path/to/BOLDGenotyper/boldgenotyper/GOaS_v1_20211214/goas_v01.shp
GOaS exists: True
Loaded [number] ocean basins
```

---

## Quick Start

**Prerequisites:**
- Complete the [Installation](#installation) steps
- (Optional) Set up [GOaS shapefile](#goas-shapefile-setup-required-for-geographic-analysis) for geographic analysis
  - Skip this if using `--no-geo` flag

### 1. Download Data from BOLD

1. Visit the [BOLD database](http://www.boldsystems.org/)
2. Search for your organism of interest
3. Download the TSV file with sequence data
4. **Important**: Name the file following the pattern: `Genus_species_commonname.tsv`
   - Example: `Sphyrna_lewini_scallopedhammerhead.tsv`

### 2. Run the Pipeline

**With Geographic Analysis (requires GOaS setup):**

The pipeline consists of sequential steps that can be run individually:

```bash
# Navigate to reference_scripts directory
cd reference_scripts

# Step 1: Convert TSV to FASTA
python BOLD_tsv_to_fasta.py ../data/Sphyrna_lewini_scallopedhammerhead.tsv

# Step 2: Generate consensus sequences (clustering and dereplication)
python msa_to_consensus.py ../data/

# Step 3: Assign genotypes to metadata
python consensus_group_to_metadata.py \
  --metadata ../data/Sphyrna_lewini_scallopedhammerhead.tsv \
  --fasta ../data/Sphyrna_lewini_scallopedhammerhead.fasta \
  --consensus ../data/consensus/Sphyrna_lewini_scallopedhammerhead_consensus.fasta \
  --out ../data/Sphyrna_lewini_scallopedhammerhead_with_consensus.tsv

# Step 4: Generate visualizations
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv ../data/Sphyrna_lewini_scallopedhammerhead_with_consensus.tsv \
  --outdir ../results/
```

**Without Geographic Analysis (no GOaS required):**

```bash
# Navigate to reference_scripts directory
cd reference_scripts

# Steps 1-3 are the same as above...

# Step 4: Skip geographic visualizations
# (Ocean basin plots will not be generated)
```

**Note:** When running without GOaS, the pipeline will automatically skip geographic analysis and continue with genotyping and phylogeny.

---

## Data Analysis Pipeline

The BOLDGenotyper pipeline consists of the following stages:

### 1. Data Parsing (`BOLD_tsv_to_fasta.py`)

**Purpose**: Convert BOLD TSV data to FASTA format for sequence analysis

**Process**:
- Extracts organism name from filename using pattern: `Genus_species_commonname`
- Parses TSV to extract sequences and metadata
- Validates required columns: `record_id`, `species`, `nuc`
- Generates FASTA file with headers: `>species_record_id`

**Output**: `{organism}.fasta`

**Example**:
```bash
python BOLD_tsv_to_fasta.py Sphyrna_lewini_scallopedhammerhead.tsv
# Creates: Sphyrna_lewini_scallopedhammerhead.fasta
```

---

### 2. Sequence Dereplication and Consensus Generation (`msa_to_consensus.py`)

**Purpose**: Identify unique genotypes through sequence clustering and generate representative consensus sequences

**Process**:
1. **Alignment**: Sequences aligned with MAFFT (`--auto` mode for optimal algorithm selection)
2. **Trimming**: Alignments trimmed with trimAl (`-automated1` for automatic threshold selection)
3. **Distance Calculation**: Pairwise distance matrix computed as `1 - % identity`, ignoring gaps (`-`) and ambiguous bases (`N`)
4. **Hierarchical Clustering**: Average linkage clustering with distance threshold (default: 0.01 = 99% identity)
5. **Consensus Generation**: Majority-rule consensus per cluster (most common base at each position)

**Output**:
- `alignments/{organism}.aln.fasta` - Raw alignment
- `trimmed/{organism}.trimmed.fasta` - Trimmed alignment
- `consensus/{organism}_consensus.fasta` - Consensus sequences with headers: `>consensus_c{cluster_id}_n{sample_count}`

**Example**:
```bash
python msa_to_consensus.py data/ --threshold 0.01 --freq-cutoff 0.7
# Creates consensus sequences in data/consensus/
```

---

### 3. Genotype Assignment (`consensus_group_to_metadata.py`)

**Purpose**: Map each sample (processid) to its best-matching consensus genotype

**Process**:
1. **ProcessID Extraction**: Extracts processid from FASTA headers (e.g., `>Sphyrna_lewini_ANGBF11456-15.COI-5P` → `ANGBF11456-15`)
2. **Sequence Alignment**: Computes global edit distance between each sample and all consensus sequences
   - Uses `edlib` library if available (fast C++ implementation)
   - Falls back to pure Python Levenshtein distance if edlib not installed
3. **Identity Calculation**: `identity = 1 - (edit_distance / max_length)`
4. **Assignment**: Sample assigned to consensus group if identity ≥ 90% (configurable)
5. **Diagnostics**: Optional CSV output with identity scores and runner-up assignments

**Output**:
- `{organism}_with_consensus.tsv` - Original TSV with added `consensus_group` column
- `{organism}_diagnostics.csv` (optional) - Identity scores and assignment details

**Example**:
```bash
python consensus_group_to_metadata.py \
  --metadata input.tsv \
  --fasta input.fasta \
  --consensus input_consensus.fasta \
  --out output_with_consensus.tsv \
  --min-identity 0.90 \
  --threads 8 \
  --diag diagnostics.csv
```

---

### 4. Geographic Filtering and Ocean Basin Assignment

**Geographic Filtering**:
The visualization pipeline automatically filters coordinates to ensure geographic precision:

- **Removes samples with**:
  - Missing lat/lon coordinates (NaN, None, empty)
  - Country-level only coordinates
  - Coordinates marked as "country centroid"

- **Rationale**: Avoids ambiguity in ocean basin assignments for countries bordering multiple ocean basins

**Land Detection and Handling**:
The pipeline detects points on land and provides three handling policies:

1. **mark** (default): Labels land-based coordinates as "country_centroid"
2. **drop**: Removes all land-based coordinates
3. **snap**: Snaps land-based points to nearest coastline

**Ocean Basin Assignment**:
Samples are classified into ocean basins using coordinate boundaries and shapefile polygons:

- North Atlantic Ocean
- South Atlantic Ocean
- North Pacific Ocean
- South Pacific Ocean
- Indian Ocean
- South China and Eastern Archipelagic Seas
- Southern Ocean
- Arctic Ocean

---

### 5. Visualization (`plot_shark_genotypes_ocean_basins_complex.py`)

**Purpose**: Generate publication-ready maps and plots showing genotype distributions

**Features**:

1. **Global Distribution Maps**:
   - Sample coordinates plotted on world map
   - Points colored by genotype (consensus group)
   - Point size scaled by sample count at location
   - Different markers for location types (precise, country centroid, snapped)

2. **Territory Polygons** (optional):
   - Buffer-based or convex hull territories around genotype locations
   - Clipping options: coastal ring, full ocean, or none
   - Polygon simplification to reduce file size
   - Overlap detection with hatching

3. **Ocean Basin Mode**:
   - Basin polygons overlaid on maps
   - Basin-specific legends and labels
   - Basin assignment for statistical analysis

4. **Multiple Output Formats**:
   - SVG (vector graphics, publication-ready)
   - PDF (high-resolution)
   - Separate "strict" versions with precise coordinates only

**Output Files**:
- `{species}_genotypes.svg` - Vector map
- `{species}_genotypes.pdf` - PDF map
- `genotype_color_map.csv` - Color assignments for genotypes
- `genotype_counts_by_species.csv` - Summary statistics
- `genotype_basin_summary_long.csv` - Genotype × basin contingency table

**Example**:
```bash
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input_with_consensus.tsv \
  --outdir results/ \
  --plot-mode basins \
  --basin-shp ocean_basins.shp \
  --draw-territories \
  --territory-buffer-deg 2.0
```

---

### 6. Statistical Analysis (R Markdown)

**Purpose**: Chi-square test of independence between genotype and ocean basin

**Process**:
1. **Contingency Table**: Builds genotype × ocean basin count matrix
2. **Chi-square Test**: Tests for non-random association
3. **Effect Size**: Calculates bias-corrected Cramér's V
4. **Standardized Residuals**: Identifies over/under-represented cells (Haberman's adjusted residuals)
5. **Visualization**: Stacked bar charts and residual heatmaps

**Output**:
- `chi2_summary.txt` - Statistical test results
- `contingency_table.csv` - Raw counts
- `proportions_by_basin.csv` - Row-normalized proportions
- `stacked_proportion_bars.pdf/svg` - Proportion visualizations
- `standardized_residuals_heatmap.pdf/svg` - Post-hoc analysis

**Example**:
```bash
Rscript -e "rmarkdown::render('reference_scripts/genotype_basin_separation.Rmd', \
  params = list(csv = 'data/plotting_data.csv', outdir = 'results/stats/'))"
```

---

## Processing Assumptions

Understanding these assumptions is critical for interpreting results:

### 1. Sequence Similarity and Distance Calculation

**Assumption**: Gaps (`-`) and ambiguous bases (`N`) are ignored in distance calculations

**Rationale**:
- Gaps may result from sequencing errors or alignment artifacts rather than true biological variation
- Ambiguous bases (`N`) lack information and shouldn't penalize similarity
- Distance calculated only over informative sites (A, C, G, T)

**Formula**:
```
distance = 1 - (matches / informative_sites)
where informative_sites = positions where both sequences have A, C, G, or T
```

---

### 2. Clustering Threshold

**Default**: 0.01 distance (99% identity)

**Rationale**:
- COI is highly conserved within species but varies between populations
- 1% divergence threshold captures intraspecific variation while grouping similar haplotypes
- Conservative approach that balances resolution with biological variation

**Customization**: Adjust with `--threshold` parameter
```bash
# More stringent (fewer, larger clusters)
python msa_to_consensus.py data/ --threshold 0.005  # 99.5% identity

# More relaxed (more, smaller clusters)
python msa_to_consensus.py data/ --threshold 0.02   # 98% identity
```

---

### 3. Consensus Generation

**Default**: Majority rule with 70% frequency cutoff

**Assumption**: Most common base at each position represents the true genotype

**Process**:
- At each alignment position, count base frequencies (ignoring gaps/Ns)
- If most common base ≥ 70% frequency → use that base
- Otherwise → assign `N` (ambiguous)

**Customization**:
```bash
# Require 80% agreement for consensus
python msa_to_consensus.py data/ --freq-cutoff 0.8

# Accept simple majority (50%)
python msa_to_consensus.py data/ --freq-cutoff 0.5
```

---

### 4. Genotype Assignment Threshold

**Default**: 90% identity minimum

**Assumption**: Samples below 90% identity to all consensus groups are unassigned

**Rationale**:
- Consensus sequences may not exactly match any individual read
- 90% threshold allows for sequencing errors and intra-cluster variation
- Stricter than clustering (99%) to ensure confident assignments

**Customization**:
```bash
python consensus_group_to_metadata.py \
  --metadata input.tsv \
  --fasta input.fasta \
  --consensus consensus.fasta \
  --out output.tsv \
  --min-identity 0.85  # More permissive
```

---

### 5. Geographic Coordinate Filtering

**Assumption**: Country-level coordinates introduce ambiguity for countries bordering multiple ocean basins

**Examples of Excluded Coordinates**:
- Missing lat/lon (NaN, empty)
- Metadata containing "country centroid"
- Metadata containing "country" without "precise"

**Rationale**:
- A "centroid" coordinate for Mexico could fall in Pacific or Atlantic
- Only precise GPS coordinates ensure accurate ocean basin assignment
- Conservative approach prioritizes accuracy over sample size

**Alternative**: Use `--land-policy snap` to snap land points to nearest coast
```bash
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input.tsv \
  --outdir results/ \
  --land-policy snap  # Keep land points by snapping to coast
```

---

### 6. Ocean Basin Boundaries

**Assumption**: Standard oceanographic definitions with custom handling for marginal seas

**Basin Definitions**:
- **Atlantic/Pacific divide**: 20°W longitude (approximate)
- **Hemispheric divide**: Equator (0° latitude)
- **Marginal seas**:
  - Mediterranean Sea → Atlantic Ocean
  - Caribbean Sea → North Atlantic Ocean
  - South China Sea → Separate category or Pacific (shapefile-dependent)

**Note**: Basin boundaries defined by shapefile polygons (if provided) or coordinate rules

---

### 7. Visualization Defaults

**Territory polygons**:
- **Method**: Buffer (smooth blobs) vs. convex hull
- **Buffer radius**: 0.5 degrees (≈55 km at equator)
- **Clipping**: Coastal ring by default (territories only in ocean near coast)

**Point sizing**:
- Bubble size = `size_min + size_scale × √(count)`
- Square root transformation for perceptual scaling

---

## Usage Guide

### Complete Pipeline Example

Here's a full workflow from BOLD download to publication figures:

```bash
# 1. Prepare your data directory
mkdir -p data results

# 2. Convert TSV to FASTA
python reference_scripts/BOLD_tsv_to_fasta.py \
  data/Sphyrna_lewini_scallopedhammerhead.tsv \
  --wrap 80

# 3. Generate consensus genotypes
python reference_scripts/msa_to_consensus.py \
  data/ \
  --threshold 0.01 \
  --freq-cutoff 0.7 \
  --wrap 80

# 4. Assign genotypes to samples
python reference_scripts/consensus_group_to_metadata.py \
  --metadata data/Sphyrna_lewini_scallopedhammerhead.tsv \
  --fasta data/Sphyrna_lewini_scallopedhammerhead.fasta \
  --consensus data/consensus/Sphyrna_lewini_scallopedhammerhead_consensus.fasta \
  --out data/Sphyrna_lewini_scallopedhammerhead_with_consensus.tsv \
  --min-identity 0.90 \
  --threads 8 \
  --diag results/diagnostics.csv

# 5a. Simple visualization (coastlines only)
python reference_scripts/plot_shark_genotypes_ocean_basins_complex.py \
  --csv data/Sphyrna_lewini_scallopedhammerhead_with_consensus.tsv \
  --outdir results/simple/ \
  --plot-mode simple \
  --alpha 0.6 \
  --size-min 12 \
  --size-scale 6 \
  --land-policy drop

# 5b. Basin visualization with territories
python reference_scripts/plot_shark_genotypes_ocean_basins_complex.py \
  --csv data/Sphyrna_lewini_scallopedhammerhead_with_consensus.tsv \
  --outdir results/basins/ \
  --plot-mode basins \
  --basin-shp ocean_basins/ne_10m_geography_marine_polys.shp \
  --basin-name-field name \
  --label-basins \
  --draw-territories \
  --territory-method buffer \
  --territory-buffer-deg 2.0 \
  --territory-alpha-fill 0.25

# 6. Statistical analysis (R Markdown)
Rscript -e "rmarkdown::render('reference_scripts/genotype_basin_separation.Rmd', \
  params = list(
    csv = 'results/basins/genotype_basin_summary_long.csv',
    outdir = 'results/statistics/',
    title_prefix = 'Sphyrna lewini'
  ))"
```

---

### Input File Requirements

#### TSV File Format

Your BOLD download must be a tab-delimited file with **required columns**:

| Column | Description | Example |
|--------|-------------|---------|
| `processid` | BOLD process ID | `ANGBF11456-15` |
| `species` | Species name | `Sphyrna lewini` |
| `nuc` or `nucleotides` | DNA sequence | `ACGTACGT...` |
| `coord` | Coordinates | `25.5, -80.2` or `25.5 -80.2` |

**Optional but recommended columns**:
- `lat`, `lon` - Separate latitude/longitude
- `country` - Country of collection
- `COI_common` - Common name
- Any metadata you want preserved

#### File Naming Convention

**Critical**: Input TSV must follow this pattern:
```
Genus_species_commonname.tsv
```

Examples:
- ✅ `Sphyrna_lewini_scallopedhammerhead.tsv`
- ✅ `Carcharodon_carcharias_greatwhiteshark.tsv`
- ❌ `shark_data.tsv` (missing organism info)
- ❌ `Sphyrna-lewini.tsv` (hyphens instead of underscores)

The pipeline extracts the organism name from the filename for all output files.

---

## Output Files

### Sequence Analysis Outputs

**Location**: `data/` (or your specified directory)

| File | Description |
|------|-------------|
| `{organism}.fasta` | Raw sequences in FASTA format |
| `alignments/{organism}.aln.fasta` | MAFFT alignment |
| `trimmed/{organism}.trimmed.fasta` | Trimmed alignment |
| `consensus/{organism}_consensus.fasta` | Consensus genotypes |
| `{organism}_with_consensus.tsv` | Metadata with genotype assignments |

---

### Visualization Outputs

**Location**: `results/simple/` or `results/basins/`

| File | Description |
|------|-------------|
| `{species}_genotypes.svg` | Vector map (publication-ready) |
| `{species}_genotypes.pdf` | PDF map (high-resolution) |
| `genotype_color_map.csv` | Color assignments for each genotype |
| `genotype_counts_by_species.csv` | Sample counts per genotype |
| `genotype_counts_by_species_coord.csv` | Counts by location and genotype |
| `species_totals_by_location_type.csv` | Summary by location type |
| `genotype_basin_summary_long.csv` | Genotype × basin long format |
| `genotype_basin_summary_pivot.csv` | Genotype × basin pivot table |

**Strict subdirectory** (`strict_precise_only/`):
- Same files but using only precise coordinates (no country centroids)

---

### Statistical Outputs

**Location**: `results/statistics/` (if R Markdown run)

| File | Description |
|------|-------------|
| `chi2_summary.txt` | Chi-square test results |
| `contingency_table.csv` | Genotype × basin counts |
| `proportions_by_basin.csv` | Row-normalized proportions |
| `stacked_proportion_bars.pdf` | Proportion visualization |
| `standardized_residuals_heatmap.pdf` | Post-hoc residual plot |

---

## Advanced Options

### Skipping Geographic Analysis

If you don't need geographic distribution analysis or don't have the GOaS shapefile installed, use the `--no-geo` flag:

```bash
# Run pipeline without geographic analysis
boldgenotyper data/Euprymna_scolopes.tsv --no-geo

# Combine with other options
boldgenotyper data/Euprymna_scolopes.tsv --no-geo --build-tree --threads 8
```

**What happens when using `--no-geo`:**
- ✓ Sequence clustering and genotyping proceed normally
- ✓ Phylogenetic tree building works (if `--build-tree` specified)
- ✓ Identity distribution plots are generated
- ✗ Ocean basin assignment is skipped (all samples marked as "Unknown")
- ✗ Geographic distribution maps are not generated
- ✗ Basin-specific visualizations are skipped

**Use cases for `--no-geo`:**
- Working with non-marine organisms (terrestrial, freshwater)
- GOaS shapefile not available or failed to download
- Only interested in genotype identification and phylogeny
- Samples lack precise geographic coordinates

### Parallel Processing

Speed up genotype assignment with multiple threads:

```bash
python consensus_group_to_metadata.py \
  --metadata input.tsv \
  --fasta input.fasta \
  --consensus consensus.fasta \
  --out output.tsv \
  --threads 16  # Use 16 CPU cores
```

---

### Custom Color Palettes

Control genotype colors in visualizations:

```bash
# Use matplotlib colormap
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input.tsv \
  --outdir results/ \
  --palette Set2  # or tab20, Dark2, etc.

# Provide explicit colors
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input.tsv \
  --outdir results/ \
  --colors "#FF0000,#00FF00,#0000FF,#FFFF00"  # Red, Green, Blue, Yellow
```

---

### Territory Polygon Customization

Fine-tune genotype territory visualization:

```bash
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input.tsv \
  --outdir results/ \
  --draw-territories \
  --territory-method convex \           # Use convex hull instead of buffer
  --territory-buffer-deg 3.0 \          # 3° buffer radius (if method=buffer)
  --territory-alpha-fill 0.15 \         # More transparent territories
  --territory-clip ocean \              # Clip to all ocean (not just coastal)
  --territory-simplify-deg 0.1 \        # Simplify polygons to reduce file size
  --territory-overlap-hatch             # Show overlapping territories with hatching
```

---

### Coordinate Order Detection

The pipeline auto-detects lat/lon vs lon/lat order, but you can force it:

```bash
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input.tsv \
  --outdir results/ \
  --coord-order latlon  # Force latitude, longitude order
```

---

### Legend Customization

Control legend appearance:

```bash
python plot_shark_genotypes_ocean_basins_complex.py \
  --csv input.tsv \
  --outdir results/ \
  --legend-counts "1,10,50,100" \       # Bubble legend example counts
  --color-legend-ncol 4 \               # Genotype legend in 4 columns
  --legend-fontsize 10 \                # Legend text size
  --legend-title-fontsize 12            # Legend title size
```

---

## Troubleshooting

### Common Issues

**1. "GOaS shapefile not found"**
```bash
# Solution 1: Run automated download
python -m boldgenotyper.goas_downloader

# Solution 2: Use --no-geo flag to skip geographic analysis
boldgenotyper data/input.tsv --no-geo

# Solution 3: Verify GOaS path
python -c "from boldgenotyper import config; \
cfg = config.get_default_config(); \
print(f'Expected GOaS location: {cfg.geographic.goas_shapefile_path}'); \
print(f'File exists: {cfg.geographic.goas_shapefile_path.exists()}')"
```

**2. "Pipeline continues but no geographic plots"**
- This is expected behavior when GOaS file is not found
- Pipeline will skip geographic analysis and continue with genotyping
- Check the log for warnings about missing GOaS file
- Install GOaS shapefile or use `--no-geo` flag explicitly

**3. "MAFFT not found in PATH"**
```bash
# Check if MAFFT is installed
conda list mafft

# If missing, install it
conda install -c bioconda mafft
```

**4. "Column 'processid' not found"**
- Verify your TSV has required columns: `processid`, `species`, `nuc`
- Check that column headers are exact (case-sensitive)

**5. "No consensus sequences generated"**
- Check that FASTA file has sequences
- Verify MAFFT and trimAl are working: `mafft --version`, `trimal --version`
- Try relaxing threshold: `--threshold 0.02`

**6. "Unassigned samples"**
- Lower the minimum identity: `--min-identity 0.85`
- Check diagnostics file to see identity scores
- Verify consensus sequences are representative

**7. "Empty territory polygons"**
- Increase buffer radius: `--territory-buffer-deg 2.0`
- Change clip mode: `--territory-clip ocean` or `--territory-clip none`
- Use convex hull method: `--territory-method convex`

---

## Citation

If you use BOLDGenotyper in your research, please cite both the software and the publication:

### Software Citation

```bibtex
@software{boldgenotyper2025,
  author = {Smith, Steph},
  title = {BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis},
  year = {2025},
  version = {0.1.0},
  url = {https://github.com/SymbioSeas/BOLDGenotyper},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

**Note**: Once published, a DOI will be assigned via Zenodo. Please use the DOI for citation when available.

### Primary Publication

```bibtex
@article{smith2025ocean,
  author = {Smith, S. and Black, C.},
  title = {Ocean basin-scale genetic partitioning in Sphyrna lewini revealed through COI sequence analysis},
  journal = {[Journal TBD]},
  year = {2025},
  note = {In preparation}
}
```

### Key Dependencies to Cite

- **BOLD Database**: Ratnasingham, S. & Hebert, P.D.N. (2007). BOLD: The Barcode of Life Data System. *Molecular Ecology Notes*, 7(3), 355-364. doi:[10.1111/j.1471-8286.2007.01678.x](https://doi.org/10.1111/j.1471-8286.2007.01678.x)

- **MAFFT**: Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772-780. doi:[10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

- **trimAl**: Capella-Gutiérrez, S., Silla-Martínez, J.M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*, 25(15), 1972-1973. doi:[10.1093/bioinformatics/btp348](https://doi.org/10.1093/bioinformatics/btp348)

- **COI Barcoding**: Hebert, P.D.N., Cywinska, A., Ball, S.L., & deWaard, J.R. (2003). Biological identifications through DNA barcodes. *Proceedings of the Royal Society B*, 270(1512), 313-321. doi:[10.1098/rspb.2002.2218](https://doi.org/10.1098/rspb.2002.2218)

---

## Contributing

We welcome contributions from the community! BOLDGenotyper is open-source and benefits from user feedback, bug reports, and code contributions.

**Quick Start for Contributors**:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes with tests and documentation
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

For detailed information on development setup, coding standards, and testing guidelines, see [CONTRIBUTING.md](CONTRIBUTING.md) (coming soon).

**Reporting Issues**: Please use [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues) to report bugs or request features.

---

## Support

**Issues**: Report bugs and request features at [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues)

**Contact**:
- Steph Smith ([@SymbioSeas](https://github.com/SymbioSeas))
- Email: steph.smith@unc.edu

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**MIT License Summary**: You are free to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of this software, provided that the copyright notice and permission notice are included in all copies or substantial portions of the software.

---

## Acknowledgments

- **BOLD Systems** for providing open access to COI sequence data
- **BioConda** community for bioinformatics tool packaging
- All contributors to MAFFT, trimAl, Biopython, and Cartopy projects

---

## Version History

### Current Release: v0.1.0 (2025-11-17)

**Core Features**:
- Complete pipeline from BOLD TSV to annotated results
- Sequence dereplication and consensus generation
- Genotype assignment with CIGAR-based identity
- Geographic analysis with GOaS integration
- Phylogenetic tree building (optional)
- Publication-ready visualizations

**For detailed version history and changelog**, see [CHANGELOG.md](CHANGELOG.md).

**Planned Features** (v0.2.0 and beyond):
- Terrestrial and freshwater organism support
- Enhanced phylogenetic analysis
- Performance optimizations for large datasets
- Docker container support

---

**Built with ❤️ for open science and marine conservation**

**Need Help?** Check out our [documentation](https://github.com/SymbioSeas/BOLDGenotyper), file an [issue](https://github.com/SymbioSeas/BOLDGenotyper/issues), or contact steph.smith@unc.edu
