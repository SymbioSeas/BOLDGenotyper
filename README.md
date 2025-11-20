# BOLDGenotyper

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://zenodo.org/)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/SymbioSeas/BOLDGenotyper/releases)

**Automated COI sequence genotyping and biogeographic analysis from BOLD database data**

BOLDGenotyper is a comprehensive bioinformatics pipeline that enables researchers to identify and analyze COI (Cytochrome Oxidase I) genotypes from the Barcode of Life Database (BOLD) for any taxonomic group. The pipeline performs sequence dereplication, consensus generation, genotype assignment, phylogenetic analysis, geographic filtering, ocean basin classification, and publication-ready visualization of biogeographic patterns.

This package enables reproducible analysis of mitochondrial COI genotypes and their geographic distributions, as demonstrated in the companion manuscript analyzing *Sphyrna lewini* (scalloped hammerhead shark) genotypes separated by ocean basin.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Obtaining Input Data from BOLD](#obtaining-input-data-from-bold)
- [Pipeline Overview](#pipeline-overview)
- [Usage Guide](#usage-guide)
- [Parameter Reference](#parameter-reference)
- [Output Files](#output-files)
- [Biological Context for Threshold Selection](#biological-context-for-threshold-selection)
- [Assumptions and Limitations](#assumptions-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Features

- **Unified Command-Line Interface**: Single command runs the complete pipeline from TSV input to publication-ready outputs
- **Automated Genotyping**: Identifies unique COI genotypes through hierarchical clustering and consensus sequence generation
- **Robust Quality Control**: Multi-stage sequence filtering, coordinate quality assessment, and centroid exclusion
- **Intelligent Genotype Assignment**: Edit distance-based matching with tie detection and low-confidence flagging
- **Geographic Analysis**: Ocean basin assignment using GOaS (Global Oceans and Seas) reference data for marine organisms
- **Phylogenetic Tree Building**: Optional phylogenetic reconstruction with FastTree (GTR+Gamma model)
- **Interactive HTML Reports**: Comprehensive summary reports with interactive visualizations using Plotly.js
- **Publication-Ready Visualizations**:
  - Identity distribution histograms
  - Phylogenetic trees with tip labels
  - Geographic distribution maps (global and faceted by species)
  - Relative and total abundance bar charts (by ocean basin)
  - All outputs in PNG and PDF formats
- **Reproducibility**: Complete parameter logging and standardized workflow
- **Flexibility**: Customizable thresholds for different taxonomic groups and research questions
- **Modular Design**: Skip geographic analysis with `--no-geo` flag for non-marine organisms or when GOaS is unavailable

---

## Installation

### Prerequisites

**Required Software**:
- **Python**: ≥3.8
- **MAFFT**: v7+ (multiple sequence alignment)
- **trimAl**: v1.4+ (alignment trimming)
- **FastTree**: v2+ (optional, for phylogenetic trees)

**Recommended**: Conda/Mamba for dependency management

### Option 1: Install from Source (Recommended for Development)

```bash
# Clone the repository
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper

# Create and activate conda environment (includes all dependencies)
conda env create -f environment.yml
conda activate boldgenotyper

# Install package in editable mode
pip install -e .

# Verify installation
boldgenotyper --version

# Alternative: If you have the old boldgenotyper_env.yml file
# conda env create -f boldgenotyper_env.yml
# conda activate depredation  # Note: old env name
```

### Option 2: Install with pip (Coming Soon)

```bash
# Once published to PyPI
pip install boldgenotyper

# Install with geographic analysis support
pip install boldgenotyper[geo]

# Install with all optional dependencies
pip install boldgenotyper[all]
```

### Verifying External Dependencies

```bash
# Check MAFFT
mafft --version

# Check trimAl
trimal --version

# Check FastTree (optional, for phylogenetic analysis)
fasttree 2>&1 | head -1

# Verify Python packages
python -c "import Bio, pandas, scipy, matplotlib, seaborn; print('Core packages OK')"

# Verify geographic packages (if installed)
python -c "import geopandas, cartopy, shapely; print('Geographic packages OK')"
```

### Installing External Tools with Conda

If you don't have MAFFT, trimAl, or FastTree installed:

```bash
conda activate boldgenotyper
conda install -c bioconda mafft trimal fasttree
```

---

## GOaS Reference Data Setup (Required for Geographic Analysis)

### What is GOaS?

GOaS (Global Oceans and Seas) is a shapefile dataset from Marine Regions that defines standardized ocean basin boundaries. BOLDGenotyper uses this dataset to:
- Assign samples to ocean basins (e.g., North Atlantic, South Pacific, Indian Ocean)
- Create geographic distribution maps with basin boundaries
- Analyze genotype-by-basin patterns for biogeographic studies

**Important Notes**:
- The GOaS shapefile (~150-200 MB) is **not included** in this repository due to size constraints
- Geographic analysis is **designed for marine organisms only** in v1.0
- If you're only interested in genotyping and phylogeny, skip this setup using the `--no-geo` flag

### Automated Download (Recommended)

```bash
# Navigate to the BOLDGenotyper directory
cd BOLDGenotyper

# Activate your conda environment
conda activate boldgenotyper

# Run the automated downloader
python -m boldgenotyper.goas_downloader

# This will:
# 1. Download World_Seas_IHO_v3.zip from Marine Regions
# 2. Extract to boldgenotyper/GOaS_v1_20211214/
# 3. Verify all required files (.shp, .shx, .dbf, .prj, .cpg)
# 4. Create a citation file
```

**Expected output**:
```
[INFO] Downloading from https://www.marineregions.org/...
[INFO] Downloaded: 100.0% (156.3 MB)
[INFO] Extracting World_Seas_IHO_v3.zip...
[INFO] ✓ All GOaS files present
[INFO] ===============================================
[INFO] GOaS setup complete!
[INFO] Data location: boldgenotyper/GOaS_v1_20211214
[INFO] ===============================================
```

### Manual Download (If Automated Fails)

1. **Download the shapefile**:
   - Visit: https://www.marineregions.org/download_file.php?name=World_Seas_IHO_v3.zip
   - Alternative: https://github.com/iobis/mregions-static/raw/master/shapefiles/World_Seas_IHO_v3.zip

2. **Extract and setup**:
   ```bash
   # Extract the ZIP file
   unzip World_Seas_IHO_v3.zip

   # Create the directory
   mkdir -p boldgenotyper/GOaS_v1_20211214

   # Move all files (must include .shp, .shx, .dbf, .prj, .cpg)
   mv World_Seas_IHO_v3.* boldgenotyper/GOaS_v1_20211214/
   ```

3. **Verify installation**:
   ```bash
   python -c "from boldgenotyper import config, geographic; \
   cfg = config.get_default_config(); \
   print(f'GOaS path: {cfg.geographic.goas_shapefile_path}'); \
   print(f'GOaS exists: {cfg.geographic.goas_shapefile_path.exists()}')"
   ```

### Skipping Geographic Analysis

If you don't need geographic distribution analysis or are working with non-marine organisms:

```bash
# Use the --no-geo flag to skip geographic analysis
boldgenotyper data/my_species.tsv --no-geo

# This will:
# ✓ Perform sequence clustering and genotyping
# ✓ Generate phylogenetic trees (if --build-tree specified)
# ✓ Create identity distribution plots
# ✗ Skip ocean basin assignment (all samples marked as "Unknown")
# ✗ Skip geographic distribution maps
# ✗ Skip basin-specific visualizations
```

**Use cases for `--no-geo`**:
- Non-marine organisms (terrestrial, freshwater)
- GOaS shapefile not available or failed to download
- Only interested in genotype identification and phylogeny
- Samples lack precise geographic coordinates

---

## Quick Start

### Step 1: Obtain Data from BOLD (Detailed instructions in next section)

Visit [boldsystems.org](http://www.boldsystems.org/) and download sequence data for your organism of interest as a TSV file.

### Step 2: Run the Pipeline

**Basic usage** (organism name inferred from filename):
```bash
boldgenotyper data/Sphyrna_lewini.tsv
```

**With phylogenetic tree**:
```bash
boldgenotyper data/Sphyrna_lewini.tsv --build-tree
```

**Specify custom output directory**:
```bash
boldgenotyper data/Sphyrna_lewini.tsv --output results/Sphyrna_analysis
```

**Without geographic analysis**:
```bash
boldgenotyper data/Euprymna_scolopes.tsv --no-geo
```

**Adjust parameters for highly diverse taxa**:
```bash
boldgenotyper data/Carcharhinus.tsv \
  --clustering-threshold 0.05 \
  --similarity-threshold 0.80 \
  --threads 8 \
  --build-tree
```

### Step 3: Review Results

Results are organized in the output directory:
```
{organism}_output/
├── {organism}_annotated.csv           # Full annotated dataset
├── {organism}_summary_report.html     # Interactive HTML report
├── {organism}_pipeline.log            # Complete log file
├── {organism}_pipeline_parameters.json # Parameters used
├── genotype_assignments/              # Assignment results
├── taxonomy/                          # Taxonomic summaries
├── phylogenetic/                      # Tree files (if --build-tree)
├── visualization/                     # PNG/PDF figures
└── reports/                           # CSV summaries
```

**Open the HTML report** to explore your results interactively:
```bash
open {organism}_output/{organism}_summary_report.html
```

---

## Obtaining Input Data from BOLD

### Step-by-Step Guide

1. **Visit the BOLD Database**:
   - Go to http://www.boldsystems.org/
   - Click on "Public Data Portal" or use direct search

2. **Search for Your Organism**:
   - Use the search box to find your species (e.g., "Sphyrna lewini")
   - You can search by:
     - Species name
     - Genus name
     - Family name
     - BIN (Barcode Index Number)
     - Process ID

3. **Select and Download Data**:
   - Click on your organism in the search results
   - Navigate to "Sequences" or "Public Records"
   - Click "Download" or "Export"
   - **Format**: Select "TSV" or "Tab-Separated Values"
   - **Options**: Ensure these columns are included:
     - `processid` (required)
     - `nucleotides` or `nuc` (required)
     - `species` or `species_name` (required)
     - `lat` and `lon` or `coord` (recommended for geographic analysis)
     - `country`, `province_state`, `region` (recommended)
     - Any other metadata you want to preserve

4. **Rename Your File** (Important):
   ```bash
   # BOLDGenotyper expects this naming convention
   # Format: Genus_species.tsv
   # Example:
   mv bold_data.tsv Sphyrna_lewini.tsv
   ```

### Required Columns

Your BOLD TSV file **must** contain these columns:

| Column | Alternative Names | Description |
|--------|------------------|-------------|
| `processid` | `process_id` | BOLD process ID (e.g., "ANGBF11456-15") |
| `nucleotides` | `nuc`, `sequence` | DNA sequence |
| `species` | `species_name` | Species name |

### Recommended Columns for Full Analysis

| Column | Purpose |
|--------|---------|
| `lat`, `lon` | Geographic coordinates (required for ocean basin assignment) |
| `coord` | Alternative format for coordinates (e.g., "25.5, -80.2") |
| `country` | Country of collection |
| `province_state` | State/province |
| `region` | Geographic region |
| `coord_accuracy` | Coordinate precision indicator |
| `bin_uri` | Barcode Index Number |
| `marker_code` | Genetic marker (usually "COI-5P") |

### Example BOLD Query

For reproducible downloads, you can use BOLD's API or direct URLs:

```bash
# Example: Download all Sphyrna lewini COI sequences
# URL format: http://www.boldsystems.org/index.php/API_Public/combined?taxon=Species_name&format=tsv

wget -O Sphyrna_lewini.tsv \
  "http://www.boldsystems.org/index.php/API_Public/combined?taxon=Sphyrna%20lewini&format=tsv"
```

### Data Quality Tips

- **Sequence Quality**: BOLD data quality varies. BOLDGenotyper includes quality filters, but starting with high-quality sequences improves results.
- **Geographic Precision**: For accurate ocean basin assignment, ensure samples have precise coordinates (not country centroids).
- **Taxonomic Consistency**: Check that species names are consistently formatted in the BOLD data.
- **Marker**: Ensure all sequences are from the same genetic marker (typically COI-5P).

---

## Pipeline Overview

BOLDGenotyper runs a comprehensive 7-phase pipeline:

### Phase 1: Data Loading and Quality Control
- Parses BOLD TSV file and validates required columns
- Filters samples based on coordinate quality (excludes country centroids)
- Assigns samples to ocean basins using GOaS shapefile (if available)
- Logs summary statistics (total samples, samples with coordinates, basin assignments)

### Phase 2: Sequence Dereplication
- Filters sequences by length (default: ≥400 bp) and N content (default: ≤10%)
- Aligns sequences with MAFFT (auto algorithm selection)
- Trims alignment with trimAl (automated1 method)
- Calculates pairwise genetic distances (ignoring gaps and ambiguous bases)
- Performs hierarchical clustering (average linkage, default threshold: 0.03)
- Generates consensus sequences for each cluster
- Filters short consensus sequences (default: ≥75% of median length)

**Output**: Consensus FASTA file with representative genotypes

### Phase 3: Genotype Assignment
- Computes edit distance between each sample and all consensus sequences
- Assigns samples to best-matching genotype above similarity threshold (default: 50%)
- Flags ambiguous assignments (tie margin default: 0.003)
- Flags low-confidence assignments (below tie threshold: 0.95)
- Generates diagnostic CSV with identity scores and runner-up matches

**Output**: Annotated dataset with genotype assignments and diagnostics

### Phase 4: Taxonomic Analysis
- Aggregates species names within each consensus group
- Determines majority species for each genotype
- Identifies taxonomy conflicts (genotypes spanning multiple species)
- Generates species composition tables

**Output**: Consensus taxonomy and species-by-genotype tables

### Phase 5: Phylogenetic Analysis (Optional, `--build-tree`)
- Aligns consensus sequences with MAFFT
- Constructs maximum-likelihood tree with FastTree (GTR+Gamma model)
- Generates tree visualizations with tip labels
- Outputs Newick format for further analysis

**Output**: Tree files (.nwk), visualizations (.png, .pdf)

### Phase 6: Visualization
- **Identity Distribution**: Histogram of sequence identity scores for assigned samples
- **Geographic Distribution Maps**:
  - Global map with genotypes color-coded
  - Faceted maps (one per species or genotype)
  - Sample points sized by abundance
- **Ocean Basin Bar Charts**:
  - Relative abundance (stacked, normalized by basin)
  - Total abundance (stacked counts)
  - Faceted versions (one per species)
- All visualizations in PNG (high-res) and PDF (vector) formats
- Interactive visualizations in HTML report

**Output**: Multiple PNG/PDF files, JSON data files for interactive plots

### Phase 7: Report Generation
- Aggregates all results into interactive HTML report
- Summary statistics (total samples, genotypes, assignment rate, etc.)
- Detailed tables (assignment status, identity scores, taxonomy, geography)
- Methods section with parameters suitable for publication
- Interactive visualizations (filterable by genotype, basin, threshold)
- Export functionality (PNG, SVG, CSV)

**Output**: Comprehensive HTML report, assignment summary CSVs

---

## Usage Guide

### Basic Command Structure

```bash
boldgenotyper <input_tsv> [options]
```

### Complete Example

```bash
# Full pipeline with all options
boldgenotyper data/Sphyrna_lewini.tsv \
  --organism "Sphyrna_lewini" \
  --output results/Sphyrna_analysis \
  --clustering-threshold 0.03 \
  --similarity-threshold 0.50 \
  --tie-margin 0.003 \
  --tie-threshold 0.95 \
  --threads 8 \
  --build-tree \
  --log-level INFO
```

### Common Use Cases

**1. Standard Analysis (Marine Organism)**:
```bash
boldgenotyper data/Sphyrna_lewini.tsv --build-tree
```

**2. Non-Marine Organism (Skip Geographic Analysis)**:
```bash
boldgenotyper data/Anopheles_gambiae.tsv --no-geo --build-tree
```

**3. Highly Diverse Taxon (Relaxed Clustering)**:
```bash
# Use looser clustering for species complexes
boldgenotyper data/Carcharhinus_complex.tsv \
  --clustering-threshold 0.05 \
  --similarity-threshold 0.80 \
  --threads 16
```

**4. Fine-Scale Population Study (Strict Clustering)**:
```bash
# Use stringent clustering for within-species genotyping
boldgenotyper data/Population_samples.tsv \
  --clustering-threshold 0.01 \
  --similarity-threshold 0.95 \
  --build-tree
```

**5. Skip HTML Report (Faster, CI/CD)**:
```bash
boldgenotyper data/my_species.tsv --no-report --threads 16
```

**6. Custom Output Location**:
```bash
boldgenotyper data/my_species.tsv \
  --output /mnt/storage/results/my_analysis_2024
```

### Input File Requirements

**Filename Convention**: Your input TSV should follow this pattern for automatic organism detection:
```
Genus_species.tsv
```

Examples:
- ✅ `Sphyrna_lewini.tsv`
- ✅ `Carcharodon_carcharias.tsv`
- ✅ `Euprymna_scolopes.tsv`
- ❌ `shark_data.tsv` (no organism info)
- ❌ `Sphyrna-lewini.tsv` (hyphens instead of underscores)

**Override organism name** with `--organism` flag if needed:
```bash
boldgenotyper data/bold_download_2024.tsv --organism Sphyrna_lewini
```

---

## Parameter Reference

### All Available Parameters

```bash
boldgenotyper --help
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tsv` | Path | Required | Input BOLD TSV file with sequences and metadata |
| `--organism` | String | Auto | Organism name (inferred from filename if not specified) |
| `--output` | Path | `./{organism}_output` | Output directory |
| `--clustering-threshold` | Float | 0.03 | Maximum genetic distance for clustering (0-1) |
| `--similarity-threshold` | Float | 0.50 | Minimum identity for genotype assignment (0-1) |
| `--tie-margin` | Float | 0.003 | Maximum identity difference to call a tie (0-1) |
| `--tie-threshold` | Float | 0.95 | Minimum identity to consider tie detection (0-1) |
| `--threads` | Integer | 4 | Number of parallel CPU threads |
| `--build-tree` | Flag | False | Build phylogenetic tree with FastTree |
| `--no-report` | Flag | False | Skip HTML report generation |
| `--no-geo` | Flag | False | Skip geographic analysis |
| `--log-level` | String | INFO | Logging verbosity (DEBUG, INFO, WARNING, ERROR) |
| `--version` | Flag | - | Show version and exit |

### Parameter Details

#### `--clustering-threshold` (default: 0.03)

**What it controls**: Maximum genetic distance for grouping sequences into consensus genotypes.

**Formula**: `distance = 1 - (sequence_identity)`
- `0.03` = sequences with ≥97% identity cluster together
- `0.01` = sequences with ≥99% identity cluster together
- `0.05` = sequences with ≥95% identity cluster together

**When to adjust**:
- **Lower (0.01-0.02)**: Fine-scale population studies, within-species genotyping
- **Higher (0.04-0.10)**: Species complexes, divergent lineages, genus-level analyses

#### `--similarity-threshold` (default: 0.50)

**What it controls**: Minimum sequence identity required for a sample to be assigned to a genotype.

**Interpretation**:
- Samples below this threshold are marked as "unassigned"
- Should be lower than clustering threshold to account for within-cluster variation
- Default of 50% is permissive to capture most samples

**When to adjust**:
- **Higher (0.85-0.95)**: High-confidence assignments only, willing to lose ambiguous samples
- **Lower (0.40-0.60)**: Retain more samples, accept some misassignments

#### `--tie-margin` (default: 0.003)

**What it controls**: Maximum identity difference between best and second-best matches to flag as ambiguous.

**Example**:
- Best match identity: 0.980
- Second-best identity: 0.978
- Difference: 0.002 < 0.003 → flagged as "tie"

**When to adjust**:
- **Lower (0.001)**: Only flag very close ties, accept more definitive assignments
- **Higher (0.01)**: Flag more ambiguous cases for manual review

#### `--tie-threshold` (default: 0.95)

**What it controls**: Minimum best-match identity required to even consider tie detection.

**Purpose**: Prevents flagging low-quality matches as ties (if both matches are poor, it's not a meaningful tie).

**Example**:
- Best match: 0.85, Second: 0.84 → Not flagged as tie (both below 0.95)
- Best match: 0.97, Second: 0.96 → Flagged as tie (both above 0.95)

#### `--threads`

**What it controls**: Number of CPU cores used for parallel processing.

**Parallelized steps**:
- Genotype assignment (edit distance calculations)
- MAFFT alignment
- FastTree construction

**Recommendation**: Set to number of available cores (check with `nproc` or `sysctl -n hw.ncpu`)

---

## Output Files

### Directory Structure

```
{organism}_output/
├── {organism}_annotated.csv              # Main output: all samples with genotype assignments
├── {organism}_summary_report.html        # Interactive HTML report
├── {organism}_pipeline.log               # Complete log with timestamps
├── {organism}_pipeline_parameters.json   # Parameters used for reproducibility
│
├── genotype_assignments/
│   ├── {organism}_diagnostics.csv        # Identity scores, ties, low-confidence flags
│   ├── {organism}_identity_distribution.png  # Histogram of identity scores
│   └── {organism}_identity_distribution.pdf
│
├── taxonomy/
│   ├── {organism}_consensus_taxonomy.csv       # Species assignments for each genotype
│   └── {organism}_species_by_consensus.csv     # Species composition tables
│
├── phylogenetic/                         # (only if --build-tree)
│   ├── {organism}_tree.nwk               # Newick format tree
│   ├── {organism}_tree_relabeled.nwk     # Tree with readable tip labels
│   ├── {organism}_tree.png               # Tree visualization
│   └── {organism}_tree.pdf
│
├── visualization/
│   ├── {organism}_identity_distribution.*       # Identity histogram
│   ├── {organism}_distribution_map.*            # Global distribution map
│   ├── {organism}_distribution_map_faceted.*    # Faceted by species
│   ├── {organism}_distribution_bar.*            # Relative abundance by basin
│   ├── {organism}_distribution_bar_faceted.*    # Faceted bar charts
│   ├── {organism}_totaldistribution_bar.*       # Total counts by basin
│   ├── *_data.json                              # Data for interactive HTML plots
│   └── {organism}_tree.* (if --build-tree)
│
├── reports/
│   └── {organism}_assignment_summary.csv        # Summary statistics
│
└── intermediate/                         # Intermediate files (for debugging)
    ├── dereplication/                    # Consensus generation outputs
    ├── genotype_assignments/             # Assignment intermediates
    ├── phylogenetic/                     # Tree intermediates
    └── geographic/                       # Basin assignment intermediates
```

### Key Output Files

#### 1. `{organism}_annotated.csv`

**Primary output** containing all samples with their genotype assignments, geographic annotations, and metadata.

**Key columns**:
- `processid`: BOLD process ID
- `consensus_group`: Assigned genotype (e.g., "Consensus_1")
- `identity`: Sequence identity to assigned genotype (0-1)
- `assignment_status`: "assigned", "low_confidence", "tie", "below_threshold", "no_sequence"
- `ocean_basin`: Assigned ocean basin (e.g., "North Atlantic Ocean")
- `species`: Original BOLD species name
- Original BOLD metadata columns preserved

#### 2. `{organism}_summary_report.html`

**Interactive HTML report** with:
- Summary statistics dashboard
- Pipeline parameters used
- Interactive visualizations (Plotly.js)
  - Filter by genotype, ocean basin, minimum sample count
  - Toggle between static and interactive views
  - Export plots as PNG/SVG
  - Download filtered data as CSV
- Methods section with publication-ready text
- Complete tables (assignment status, taxonomy, geography)

**To view**: Open in any modern web browser

#### 3. `{organism}_diagnostics.csv`

**Detailed assignment diagnostics** for troubleshooting and quality control.

**Columns**:
- `processid`: Sample ID
- `assigned_consensus`: Best-matching genotype
- `best_identity`: Identity to best match
- `runner_up_consensus`: Second-best match
- `runner_up_identity`: Identity to second-best match
- `is_tie`: Boolean indicating ambiguous assignment
- `is_low_confidence`: Boolean indicating low-confidence assignment

#### 4. Visualization Files

All visualization outputs are provided in multiple formats:
- **PNG**: High-resolution (300 DPI) for presentations, web
- **PDF**: Vector graphics for publication, posters
- **JSON**: Data files for interactive HTML plots

**Types of visualizations**:
- **Identity Distribution** (`*_identity_distribution.*`): Histogram showing distribution of sequence identity scores for assigned samples
- **Geographic Maps** (`*_distribution_map.*`): World map with genotypes color-coded and sample points
- **Faceted Maps** (`*_distribution_map_faceted.*`): Separate map for each species or major genotype
- **Relative Abundance** (`*_distribution_bar.*`): Stacked bar chart showing genotype proportions within each ocean basin
- **Total Abundance** (`*_totaldistribution_bar.*`): Stacked bar chart showing raw sample counts by basin
- **Faceted Bar Charts** (`*_distribution_bar_faceted.*`): Separate bar chart for each species
- **Phylogenetic Tree** (`*_tree.*`): Maximum-likelihood tree with tip labels (if `--build-tree`)

#### 5. `{organism}_pipeline.log`

**Complete pipeline log** with timestamps, including:
- Input validation and data loading statistics
- Dereplication metrics (clusters, consensus sequences)
- Assignment summary (success rate, ties, low-confidence)
- Geographic analysis results (basin assignments)
- Warnings and errors encountered
- Runtime for each phase

---

## Biological Context for Threshold Selection

### Understanding Clustering Threshold

The **clustering threshold** determines how genetically similar sequences must be to group into the same consensus genotype.

**COI as a Molecular Marker**:
- COI is typically conserved within species (>97% identity)
- Can vary between populations (1-5% divergence)
- Varies significantly between species (>5-10% divergence)

**Choosing Your Threshold**:

| Research Question | Recommended Threshold | Rationale |
|------------------|----------------------|-----------|
| **Within-species population structure** | 0.005-0.01 (99-99.5%) | Capture fine-scale haplotype variation |
| **Species delimitation** | 0.03-0.05 (95-97%) | Standard COI barcoding gap |
| **Genus-level diversity** | 0.05-0.10 (90-95%) | Group closely related species |
| **Species complex** | 0.05-0.10 (90-95%) | Handle cryptic species |

**Example**: For *Sphyrna lewini* population study across ocean basins, `--clustering-threshold 0.03` captures intraspecific genotypes while avoiding oversplitting due to sequencing errors.

### Understanding Similarity Threshold

The **similarity threshold** determines the minimum identity required for a sample to be assigned to a genotype.

**Why Different from Clustering?**:
- Consensus sequences may not exactly match any individual sequence
- Allows for within-cluster variation
- Should be lower than clustering threshold

**Recommended Values**:
- **Conservative (0.85-0.95)**: High-confidence assignments, willing to exclude ambiguous samples
- **Moderate (0.70-0.85)**: Balance between coverage and accuracy
- **Permissive (0.50-0.70)**: Maximize sample retention, accept more ambiguity

**Default (0.50)**: Retains most samples while still requiring ≥50% identity to the genotype representative.

### Understanding Tie Detection

**Tie margin** and **tie threshold** work together to identify ambiguous assignments.

**When are ties important?**:
- Detecting potential cryptic species
- Identifying samples at genotype boundaries
- Quality control for borderline cases

**Example Scenario**:
- Sample matches Genotype_A at 97.0% identity
- Sample matches Genotype_B at 96.8% identity
- Difference: 0.2% < tie-margin (0.3%) → Flagged for manual review
- Both above tie-threshold (95%) → Meaningful comparison

**Adjusting for Your Data**:
- **Closely related genotypes**: Lower tie-margin (0.001-0.002) to flag only very close calls
- **Distinct genotypes**: Higher tie-margin (0.005-0.01) to catch more ambiguous cases

### Guidelines by Taxonomic Group

| Taxonomic Group | Clustering | Similarity | Notes |
|----------------|------------|-----------|-------|
| **Elasmobranchs (sharks, rays)** | 0.02-0.03 | 0.80-0.90 | Slow molecular evolution, lower divergence |
| **Teleost fish** | 0.03-0.05 | 0.80-0.90 | Standard barcoding thresholds |
| **Marine invertebrates** | 0.03-0.05 | 0.70-0.85 | Higher COI variation |
| **Cephalopods (squid, octopus)** | 0.03-0.05 | 0.80-0.90 | Similar to fish |
| **Crustaceans** | 0.05-0.10 | 0.70-0.85 | Often more divergent |
| **Marine mammals** | 0.01-0.03 | 0.85-0.95 | Low COI variation |

**Literature Support**:
- Ward et al. (2005): 97-98% threshold for fish species delimitation
- Hebert et al. (2003): 3% divergence as standard barcoding gap
- Costa et al. (2012): 2% threshold for elasmobranch species identification

---

## Assumptions and Limitations

Understanding these assumptions is critical for interpreting results and planning your analysis.

### Assumptions

#### 1. Sequence Quality
- **Assumption**: Input sequences are of reasonable quality
- **Pipeline behavior**: Filters sequences <400bp and >10% N content
- **Implication**: Low-quality sequences are excluded; BOLD quality varies by contributor
- **Recommendation**: Start with "public records" or sequences with BIN assignments

#### 2. Gap and Ambiguity Handling
- **Assumption**: Gaps (`-`) and ambiguous bases (`N`) are ignored in distance calculations
- **Rationale**: Gaps may result from alignment artifacts; Ns lack information
- **Calculation**: `distance = 1 - (matches / informative_sites)`
  where informative sites = positions with A, C, G, or T in both sequences
- **Implication**: Fragmented sequences with many gaps aren't penalized artificially

#### 3. Consensus Generation
- **Assumption**: Majority-rule consensus (70% frequency cutoff) represents true genotype
- **Process**: At each alignment position, if most common base ≥70% → use that base; otherwise → N
- **Implication**: Intra-cluster variation collapsed; minority variants not represented
- **Adjustment**: Change with `consensus_frequency_cutoff` in config (not exposed in CLI)

#### 4. Geographic Coordinate Quality
- **Assumption**: Country-level coordinates introduce ambiguity for multi-basin countries
- **Pipeline behavior**: Excludes:
  - Missing lat/lon (NaN, empty)
  - Samples with "country centroid" in metadata
  - Country-level coordinates without "precise" indicator
- **Rationale**: A centroid for Mexico could fall in Pacific or Atlantic
- **Implication**: Conservative approach prioritizes accuracy over sample size
- **Override**: Not currently possible via CLI (feature request welcome)

#### 5. Ocean Basin Boundaries
- **Assumption**: GOaS polygon definitions are accurate and appropriate
- **Implementation**:
  - Shapefile defines basin polygons
  - Coordinate-to-basin assignment via spatial join
  - Marginal seas (Mediterranean, Caribbean) included in parent basins
- **Limitation**: Disputed waters, transitional zones not explicitly handled
- **Resolution**: Samples near boundaries assigned to single basin based on polygon overlap

#### 6. Phylogenetic Model
- **Assumption**: GTR+Gamma model appropriate for COI sequences (if `--build-tree`)
- **Justification**: GTR allows different substitution rates; Gamma handles rate heterogeneity
- **Limitation**: FastTree is approximate ML (not exhaustive search like RAxML/IQ-TREE)
- **Use case**: Phylogenetic tree is for visualization and preliminary analysis, not deep phylogenetics

### Limitations

#### 1. Marine Organisms Only (v1.0)
- **Current**: Geographic analysis designed for marine taxa using GOaS
- **Limitation**: No terrestrial or freshwater basin/region definitions
- **Workaround**: Use `--no-geo` flag for non-marine organisms
- **Planned**: Future versions will support terrestrial/freshwater environments

#### 2. Coordinate Precision Requirements
- **Requirement**: Geographic analysis requires precise GPS coordinates
- **Problem**: Many BOLD records have country-level or imprecise coordinates
- **Impact**: Samples with low-quality coordinates excluded from basin assignment
- **Quantification**: Check HTML report for "samples excluded by coordinate filter"

#### 3. Single-Locus Limitation
- **Current**: Pipeline designed for single-locus data (COI)
- **Limitation**: Cannot integrate nuclear markers, multiple mitochondrial genes
- **Rationale**: BOLD primarily stores COI data; most users analyzing COI
- **Extension**: Could be adapted for other markers with appropriate thresholds

#### 4. No Haplotype Network Construction
- **Current**: Phylogenetic tree only (if `--build-tree`)
- **Missing**: Haplotype networks (TCS, median-joining) for population genetics
- **Workaround**: Export genotypes and use PopART, Network, or POPART
- **Planned**: May add in future versions

#### 5. No Population Genetics Statistics
- **Missing**:
  - Fst between populations/basins
  - AMOVA (Analysis of Molecular Variance)
  - Nucleotide diversity (π), haplotype diversity (Hd)
  - Tajima's D, Fu's Fs
- **Rationale**: Focused on genotype identification and biogeography, not population genetics
- **Workaround**: Export data and use dedicated tools (Arlequin, DnaSP, MEGA)

#### 6. No Outgroup Specification (Yet)
- **Current**: Trees unrooted or midpoint-rooted
- **Missing**: Ability to specify outgroup sequences for rooting
- **Impact**: Tree topology may be misleading without proper root
- **Workaround**: Re-root manually in FigTree, iTOL, or R (ape package)
- **Planned**: `--outgroup` parameter in future version

#### 7. Performance on Very Large Datasets
- **Tested**: Up to ~10,000 sequences
- **Potential bottlenecks**:
  - MAFFT alignment for >20,000 sequences
  - Pairwise distance matrix calculation (O(n²))
  - Clustering (O(n² log n))
- **Mitigation**: Use `--threads` for parallelization
- **Future**: Implement subsampling or approximate methods for >50,000 sequences

#### 8. Visualization Scalability
- **Challenge**: Maps with >100 genotypes become cluttered
- **Current**: Color palette recycles if >20 genotypes
- **Impact**: Difficult to distinguish genotypes visually
- **Workaround**: Filter HTML report interactively to focus on subset
- **Future**: Implement genotype grouping/collapsing options

### Known Issues and Feature Requests

See [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues) for:
- Bug reports
- Feature requests
- Planned enhancements
- Community contributions

**Top feature requests** (as of v0.1.0):
1. Outgroup specification for tree rooting
2. Terrestrial/freshwater geographic modules
3. Haplotype network construction
4. Population genetics statistics integration
5. Multi-locus support
6. Docker containerization
7. Galaxy tool integration
8. Web interface

---

## Troubleshooting

### Installation Issues

**Problem**: `command not found: boldgenotyper`

```bash
# Solution 1: Ensure package is installed
pip install -e .

# Solution 2: Check that conda environment is activated
conda activate boldgenotyper

# Solution 3: Run as Python module
python -m boldgenotyper data/my_species.tsv
```

**Problem**: `ImportError: No module named 'Bio'`

```bash
# Install core dependencies
pip install biopython pandas scipy matplotlib seaborn

# Or reinstall from environment file
conda env update -f environment.yml
```

**Problem**: `MAFFT not found in PATH`

```bash
# Check if MAFFT is installed
conda list mafft

# Install if missing
conda install -c bioconda mafft

# Verify installation
mafft --version
```

### Runtime Errors

**Problem**: `GOaS shapefile not found`

```bash
# Solution 1: Run automated downloader
python -m boldgenotyper.goas_downloader

# Solution 2: Skip geographic analysis
boldgenotyper data/my_species.tsv --no-geo

# Solution 3: Verify GOaS path
python -c "from boldgenotyper import config; \
cfg = config.get_default_config(); \
print(cfg.geographic.goas_shapefile_path)"
```

**Problem**: `ValueError: Column 'processid' not found in TSV`

```bash
# Check your TSV has required columns
head -1 data/my_species.tsv

# Verify column names (case-sensitive)
# Required: processid, nucleotides (or nuc), species

# If columns have different names, rename them:
# Option 1: Edit TSV manually
# Option 2: Use sed/awk to rename headers
```

**Problem**: `No consensus sequences generated`

```bash
# Check that input sequences exist
grep -c ">" data/my_species.fasta  # (if FASTA exists)

# Check log file for details
tail -100 {organism}_output/{organism}_pipeline.log

# Possible causes:
# 1. All sequences too short (<400bp) or too many Ns (>10%)
# 2. MAFFT/trimAl failed
# 3. All sequences filtered after trimming

# Try more permissive settings (edit config or future CLI flags)
```

**Problem**: `Most samples unassigned`

```bash
# Check diagnostics file
less {organism}_output/genotype_assignments/{organism}_diagnostics.csv

# Look at identity scores - are they all low?
cut -d',' -f3 {organism}_output/genotype_assignments/{organism}_diagnostics.csv | sort -n

# If identities are 0.4-0.6, consensus sequences may not represent data well
# Solutions:
# 1. Adjust clustering threshold (e.g., --clustering-threshold 0.05)
# 2. Lower similarity threshold (e.g., --similarity-threshold 0.40)
# 3. Check if sequences are actually from the same marker/region
```

**Problem**: `FastTree failed` (if using `--build-tree`)

```bash
# Verify FastTree is installed
fasttree 2>&1 | head -1

# Install if missing
conda install -c bioconda fasttree

# Check that consensus sequences exist
ls {organism}_output/intermediate/dereplication/

# Try running pipeline without tree first
boldgenotyper data/my_species.tsv  # omit --build-tree
```

### Data Quality Issues

**Problem**: Too many genotypes (oversplitting)

```bash
# Increase clustering threshold to group more sequences
boldgenotyper data/my_species.tsv --clustering-threshold 0.05

# This creates fewer, larger genotype clusters
```

**Problem**: Too few genotypes (undersplitting)

```bash
# Decrease clustering threshold for finer resolution
boldgenotyper data/my_species.tsv --clustering-threshold 0.01

# This creates more, smaller genotype clusters
```

**Problem**: Many samples flagged as ties

```bash
# Option 1: Accept ties as real biological signal (genotypes not clearly distinct)
# Option 2: Increase tie-margin to reduce tie flagging
boldgenotyper data/my_species.tsv --tie-margin 0.001

# Option 3: Manual review - inspect diagnostics file
```

**Problem**: No samples in certain ocean basins

```bash
# Check if GOaS assignment worked
grep -c "Unknown" {organism}_output/{organism}_annotated.csv

# If many "Unknown", check:
# 1. Coordinate quality (might be filtered out)
# 2. Sample actually in ocean (could be land-based)
# 3. GOaS shapefile loaded correctly
```

### Performance Issues

**Problem**: Pipeline very slow

```bash
# Use more threads
boldgenotyper data/large_dataset.tsv --threads 16

# Check which step is slow (see log file)
tail -f {organism}_output/{organism}_pipeline.log

# Bottlenecks:
# - MAFFT (for >10,000 sequences): Use --threads
# - Genotype assignment: Use --threads
# - Tree building: FastTree is already fast
```

**Problem**: Out of memory

```bash
# Reduce thread count (paradoxically helps for memory-limited systems)
boldgenotyper data/my_species.tsv --threads 2

# Or subsample your data before running
head -10000 data/my_species.tsv > data/my_species_subset.tsv
```

### Visualization Issues

**Problem**: Blank or empty maps

```bash
# Check if samples have coordinates
cut -d',' -f<lat_column>,<lon_column> {organism}_output/{organism}_annotated.csv

# Check if samples assigned to ocean basins
cut -d',' -f<ocean_basin_column> {organism}_output/{organism}_annotated.csv | sort | uniq -c

# If all "Unknown", geographic analysis was skipped or failed
```

**Problem**: HTML report plots not interactive

```bash
# Possible causes:
# 1. JavaScript disabled in browser
# 2. JSON data files missing (check visualization/ directory)
# 3. Plotly.js CDN blocked (requires internet connection)

# Check browser console for errors (F12)
```

### Getting Help

1. **Check the log file**: `{organism}_output/{organism}_pipeline.log`
2. **Enable debug logging**: `boldgenotyper data/my_species.tsv --log-level DEBUG`
3. **Search existing issues**: [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues)
4. **Ask a question**: Open a new issue with:
   - Command you ran
   - Error message (full traceback)
   - Relevant portions of log file
   - BOLDGenotyper version (`boldgenotyper --version`)
   - Operating system and Python version

---

## Example Dataset: *Sphyrna lewini*

The repository includes a complete example analysis of scalloped hammerhead shark (*Sphyrna lewini*) in `data/Sphyrna_lewini/`.

### Reproducing the Example

```bash
# Navigate to the repository root
cd BOLDGenotyper

# Run the analysis (results already included)
boldgenotyper data/Sphyrna_lewini_input.tsv --build-tree --output data/Sphyrna_lewini

# Open the HTML report
open data/Sphyrna_lewini/Sphyrna\ lewini_summary_report.html
```

### What's Included

```
data/Sphyrna_lewini/
├── Sphyrna lewini_summary_report.html    # Interactive report
├── Sphyrna lewini_annotated.csv          # Full annotated dataset (900+ samples)
├── genotype_assignments/                 # Assignment diagnostics
├── taxonomy/                             # Taxonomic summaries
├── phylogenetic/                         # Tree files and visualizations
└── visualization/                        # Publication-ready figures
```

### Key Findings from Example

- **Input**: 937 *Sphyrna lewini* COI sequences from BOLD
- **Genotypes**: 26 consensus genotypes identified (clustering threshold: 0.03)
- **Assignment Rate**: ~85% of samples assigned to genotypes
- **Geographic Distribution**: Samples from North Atlantic, South Atlantic, North Pacific, South Pacific, and Indian Oceans
- **Biogeographic Pattern**: Distinct genotypes show basin-specific distributions, supporting ocean basin-scale population structure

### Using This Example for Learning

1. **Explore the HTML report** to understand all output types
2. **Compare parameter choices** in `Sphyrna lewini_pipeline_parameters.json`
3. **Examine visualizations** to see publication-quality figure examples
4. **Check diagnostics CSV** to understand assignment quality metrics
5. **Review log file** to see pipeline progress and summary statistics

---

## Citation

If you use BOLDGenotyper in your research, please cite:

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

### Primary Publication (In Preparation)

```bibtex
@article{smith2025ocean,
  author = {Smith, S. and Black, C.},
  title = {Ocean basin-scale genetic partitioning in Sphyrna lewini revealed through COI sequence analysis},
  journal = {[Journal TBD]},
  year = {2025},
  note = {In preparation}
}
```

### Key Dependencies

Please also cite these foundational tools:

- **BOLD Database**: Ratnasingham, S. & Hebert, P.D.N. (2007). BOLD: The Barcode of Life Data System. *Molecular Ecology Notes*, 7(3), 355-364. [doi:10.1111/j.1471-8286.2007.01678.x](https://doi.org/10.1111/j.1471-8286.2007.01678.x)

- **MAFFT**: Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772-780. [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

- **trimAl**: Capella-Gutiérrez, S., Silla-Martínez, J.M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*, 25(15), 1972-1973. [doi:10.1093/bioinformatics/btp348](https://doi.org/10.1093/bioinformatics/btp348)

- **FastTree**: Price, M.N., Dehal, P.S., & Arkin, A.P. (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. *PLoS ONE*, 5(3), e9490. [doi:10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)

- **COI Barcoding**: Hebert, P.D.N., Cywinska, A., Ball, S.L., & deWaard, J.R. (2003). Biological identifications through DNA barcodes. *Proceedings of the Royal Society B*, 270(1512), 313-321. [doi:10.1098/rspb.2002.2218](https://doi.org/10.1098/rspb.2002.2218)

- **GOaS Dataset**: Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at https://www.marineregions.org/

---

## Contributing

We welcome contributions from the community! BOLDGenotyper is open-source and benefits from user feedback, bug reports, and code contributions.

### How to Contribute

1. **Report bugs**: Use [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues)
2. **Request features**: Open an issue with the "enhancement" label
3. **Improve documentation**: Submit pull requests for README, docstrings, or examples
4. **Add test cases**: Help us improve test coverage
5. **Contribute code**:
   - Fork the repository
   - Create a feature branch (`git checkout -b feature/amazing-feature`)
   - Make your changes with tests and documentation
   - Commit your changes (`git commit -m 'Add amazing feature'`)
   - Push to the branch (`git push origin feature/amazing-feature`)
   - Open a Pull Request

### Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR-USERNAME/BOLDGenotyper.git
cd BOLDGenotyper

# Create development environment
conda env create -f environment.yml
conda activate boldgenotyper

# Install in editable mode with development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/

# Check code style
black boldgenotyper/
flake8 boldgenotyper/
```

### Code Standards

- **Python Style**: Follow PEP 8 (enforced by black and flake8)
- **Docstrings**: NumPy style for all functions and classes
- **Type Hints**: Use type annotations for function signatures
- **Tests**: Add tests for new features (pytest)
- **Documentation**: Update README and docstrings

### Areas Needing Contribution

- **Terrestrial/freshwater modules** for non-marine organisms
- **Haplotype network construction** using TCS or median-joining methods
- **Population genetics statistics** (Fst, AMOVA, diversity indices)
- **Performance optimization** for large datasets (>50,000 sequences)
- **Additional output formats** (BEAST XML, Arlequin, STRUCTURE)
- **Web interface** using Flask or Streamlit
- **Docker container** for reproducibility
- **Galaxy tool wrapper** for integration

---

## Support

### Getting Help

- **Documentation**: You're reading it! Check [Usage Guide](#usage-guide) and [Troubleshooting](#troubleshooting)
- **GitHub Issues**: [Report bugs or ask questions](https://github.com/SymbioSeas/BOLDGenotyper/issues)
- **Email**: Steph Smith (steph.smith@unc.edu)

### Frequently Asked Questions

**Q: Can I use BOLDGenotyper for non-COI markers?**
A: Yes, but you may need to adjust thresholds. The default clustering threshold (0.03) is COI-specific. For more variable markers, increase the threshold; for more conserved markers, decrease it.

**Q: How do I analyze terrestrial or freshwater organisms?**
A: Use the `--no-geo` flag to skip ocean basin assignment. Geographic modules for terrestrial environments are planned for v0.2.0.

**Q: Can I combine data from multiple BOLD downloads?**
A: Yes, concatenate TSV files (keeping headers from first file) before running the pipeline.

**Q: How long does the pipeline take?**
A: Depends on dataset size:
- ~1,000 sequences: 5-15 minutes
- ~5,000 sequences: 30-60 minutes
- ~10,000 sequences: 1-3 hours
Use `--threads` to speed up parallelizable steps.

**Q: Can I customize visualizations (colors, labels, etc.)?**
A: Currently limited CLI options. Advanced customization requires editing config or directly calling visualization module. Interactive HTML report allows filtering and export.

**Q: What if my species is not in BOLD?**
A: You can use any FASTA file with similar format. See documentation for custom input preparation (future feature).

**Q: How do I export data for other software (Arlequin, POPART, etc.)?**
A: Use the annotated CSV and consensus FASTA files. Conversion scripts for popular formats are planned.

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

**MIT License Summary**: You are free to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of this software, provided that the copyright notice and permission notice are included in all copies or substantial portions of the software.

---

## Acknowledgments

- **BOLD Systems** for providing open access to COI sequence data and maintaining a critical resource for DNA barcoding research
- **BioConda** community for bioinformatics tool packaging and distribution
- **Marine Regions** for providing the GOaS dataset and maintaining marine geographic reference data
- All contributors to MAFFT, trimAl, FastTree, BioPython, Cartopy, GeoPandas, Plotly, and related open-source projects
- Dr. Chelsea Black for collaboration on the *Sphyrna lewini* case study
- The University of North Carolina at Chapel Hill for supporting this work

---

## Version History

### v0.1.0 (2025-01-20) - Initial Release

**Core Features**:
- Complete unified pipeline from BOLD TSV to annotated results
- Automated sequence dereplication and consensus generation
- Edit distance-based genotype assignment with diagnostics
- Geographic analysis with GOaS ocean basin integration
- Optional phylogenetic tree building with FastTree
- Comprehensive interactive HTML reports
- Publication-ready visualizations (PNG/PDF)
- Multi-threaded processing support

**Modules**:
- `boldgenotyper.cli`: Command-line interface
- `boldgenotyper.metadata`: TSV parsing and coordinate filtering
- `boldgenotyper.dereplication`: Sequence clustering and consensus generation
- `boldgenotyper.genotype_assignment`: Sample-to-genotype matching
- `boldgenotyper.geographic`: Ocean basin assignment
- `boldgenotyper.phylogenetics`: Tree construction and visualization
- `boldgenotyper.visualization`: Distribution maps and abundance charts
- `boldgenotyper.reports`: HTML report generation
- `boldgenotyper.config`: Configuration management

**Known Limitations**:
- Marine organisms only (GOaS-based)
- No terrestrial/freshwater support
- No outgroup specification for tree rooting
- No haplotype network construction
- No population genetics statistics

---

**Planned Features** (v0.2.0 and beyond):
1. Terrestrial and freshwater organism support with appropriate reference datasets
2. Outgroup specification for phylogenetic tree rooting (`--outgroup` flag)
3. Haplotype network construction (TCS, median-joining)
4. Population genetics statistics (Fst, AMOVA, diversity indices)
5. Multi-locus support for concatenated/combined markers
6. Performance optimizations for datasets >50,000 sequences
7. Docker container for reproducibility and portability
8. Galaxy tool integration for workflow platforms
9. Web interface for browser-based analysis
10. Additional export formats (BEAST, Arlequin, STRUCTURE)

---

**Built with ❤️ for open science and marine conservation**

**Created by**: Steph Smith ([@SymbioSeas](https://github.com/SymbioSeas))

**Questions or feedback?**
- GitHub: https://github.com/SymbioSeas/BOLDGenotyper
- Email: steph.smith@unc.edu
- Issues: https://github.com/SymbioSeas/BOLDGenotyper/issues
