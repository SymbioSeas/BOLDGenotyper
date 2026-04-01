# BOLDGenotyper

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://zenodo.org/)
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/SymbioSeas/BOLDGenotyper/releases)

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
- [Customizing Plots for Publication](#customizing-plots-for-publication)
- [Biological Context for Threshold Selection](#biological-context-for-threshold-selection)
- [Assumptions and Limitations](#assumptions-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Features

- **Unified Command-Line Interface**: Single command runs the complete pipeline from TSV input to publication-ready outputs
- **Automated Haplotype Discovery**: Identifies unique COI haplotypes using Exact Sequence Variants (ESVs) from a core alignment region
- **Robust Quality Control**: Multi-stage sequence filtering, coordinate quality assessment, and centroid exclusion
- **Intelligent Genotype Assignment**: Edit distance-based matching with tie detection and low-confidence flagging
- **Geographic Analysis**: Region assignment using GOaS ocean basins (marine), HydroBASINS (freshwater), Ecoregions2017 (terrestrial), or any custom shapefile
- **Metadata Analysis**: Automated analysis of haplotype–metadata associations with chi-square testing, coverage assessment, and temporal emergence timelines
- **Phylogenetic Tree Building**: Optional phylogenetic reconstruction with FastTree (GTR+Gamma model)
- **Interactive HTML Reports**: Comprehensive summary reports with interactive visualizations using Plotly.js
- **Publication-Ready Visualizations**:
  - Identity distribution histograms
  - Phylogenetic trees with tip labels
  - Geographic distribution maps (global and faceted by species)
  - Relative and total abundance bar charts by geographic region
  - All outputs in PNG and PDF formats
- **Plot Customization**: YAML-configured plot regeneration — modify colors, dimensions, projections, and filters without re-running the analysis
- **Population Genetics Export**: Optional export to Arlequin, PopART/NEXUS, DnaSP, and generic formats with geographic populations defined automatically
- **Reproducibility**: Complete parameter logging and standardized workflow
- **Flexibility**: Customizable thresholds for different taxonomic groups and research questions
- **Modular Design**: Skip geographic analysis with `--no-geo` flag for non-marine organisms or when GOaS is unavailable

---

## Active Modules

### Layer 1: Core Infrastructure
- **cli.py**: Parses CLI args, calls other modules in sequence
- **config.py**: Dataclass-based config with nested sections for each module
- **utils.py**: Logging setup, external tool checking (MAFFT), file I/O, string sanitation

### Layer 2: QC & Preprocessing
- **metadata.py**: Parses BOLD TSV files, validates required columns, coordinate quality marking
- **quality_control.py**: COI orientation normalization, ORF validation, dynamic length/N-fraction filtering

### Layer 3: Haplotype Discovery & Assignment
- **dereplication.py**: MAFFT alignment, core region extraction, ESV identification, suspect flagging
- **haplotype_assignment.py**: Edit-distance matching of samples to haplotypes, confidence scoring

### Layer 4: Analysis
- **divergence_analysis.py**: Pairwise divergence matrices, barcoding gap, within/between-species divergence
- **species_analysis.py**: Species-level aggregation, diversity metrics, geographic summaries
- **scripts/compare_analyses.py**: Standalone utility — compares species vs. family-level runs, contamination crosswalk
- **parameter_sweep.py**: Tests multiple singleton thresholds, tracks group membership stability

### Layer 5: Geographic & Phylogenetic
- **geographic.py**: GOaS shapefile loading, point-in-polygon ocean basin assignment
- **geographic_enhancement.py**: Coverage assessment, basin confidence metrics
- **goas_downloader.py**: Verifies Global Oceans and Seas reference shapefile installation
- **phylogenetics.py**: MAFFT alignment, trimAl, FastTree tree building, outgroup rooting
- **msa_visualization.py**: Phylogeny-ordered MSA display with nucleotide coloring

### Layer 6: Visualization & Reporting
- **visualization.py**: Distribution maps, ocean basin bar charts, haplotype pie charts, faceted plots
- **metadata_analysis.py**: Categorical metadata field analysis and plotting
- **reports.py**: Assignment summary CSV, HTML report with embedded images
- **plot_export.py**: Exports plot data CSVs + config YAML + Python regeneration scripts
- **plot_regeneration.py**: Regenerates plots from exported data using customized config

### Layer 7: Export
- **popgen_export.py**: Exports to Arlequin, PopART/NEXUS, DnaSP, generic CSV/FASTA
- **haplotype_query.py**: Queries user sequences against haplotype database

## Installation

### Quick Install (Recommended)

BOLDGenotyper uses conda to manage all dependencies, including external tools (MAFFT, trimAl, FastTree) and Python packages.

**For Mac/Linux**:
```bash
# Clone the repository
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper

# Run installation script (installs everything automatically)
bash install.sh
```

**For Windows (requires WSL2)**:

BOLDGenotyper depends on bioconda packages (MAFFT, trimAl, FastTree) which are not available on native Windows. You must use [WSL2 (Windows Subsystem for Linux)](https://learn.microsoft.com/en-us/windows/wsl/install).

```powershell
# Step 1: Install WSL2 (if not already installed)
# Open PowerShell as Administrator and run:
wsl --install

# Step 2: After restarting, open your WSL terminal and install Miniconda:
# https://docs.conda.io/en/latest/miniconda.html#linux-installers

# Step 3: Clone and install inside WSL
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper
bash install.sh
```

Alternatively, from a Windows Command Prompt (if WSL2 is already set up with conda):
```
install.bat
```

**Manual Installation** (if scripts don't work):
```bash
# Clone the repository
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper

# Create and activate conda environment (includes ALL dependencies)
conda env create -f environment.yml
conda activate boldgenotyper

# Install package
pip install -e .

# Verify installation
boldgenotyper --version
```

### What Gets Installed

The conda environment includes:

**External Tools**:
- MAFFT (multiple sequence alignment)
- trimAl (alignment trimming)
- FastTree (phylogenetic tree building)

**Python Packages**:
- Core: biopython, pandas, numpy, scipy, matplotlib, seaborn
- Geographic: geopandas, cartopy, shapely (for maps and polygon analysis)
- Utilities: pyyaml, jinja2

**No additional installation needed!** All dependencies are managed by conda.

### Verifying Installation

```bash
# Activate environment
conda activate boldgenotyper

# Check boldgenotyper
boldgenotyper --version

# Check external tools
mafft --version
trimal --version
fasttree 2>&1 | head -1
```

### System Requirements

- **Operating Systems**: macOS, Linux, Windows (via WSL2 — see above)
- **Python**: 3.9 - 3.13
- **RAM**: 4GB minimum, 8GB+ recommended for large datasets
- **Disk**: ~2GB for conda environment + space for your data

---

## GOaS Reference Data Setup (Required for Geographic Analysis)

### What is GOaS?

GOaS (Global Oceans and Seas) is a shapefile dataset from Marine Regions that defines standardized ocean basin boundaries. BOLDGenotyper uses this dataset to:
- Assign samples to ocean basins (e.g., North Atlantic, South Pacific, Indian Ocean)
- Create geographic distribution maps with basin boundaries
- Analyze genotype-by-basin patterns for biogeographic studies

**Important Notes**:
- The GOaS shapefile (~150-200 MB) is **not included** in this repository due to size constraints
- GOaS is the **default for marine organisms**; use `--custom-shp` for freshwater/terrestrial (see below)
- If you're only interested in genotyping and phylogeny, skip this setup using the `--no-geo` flag

### Manual Download (Required)

> **Note**: Marine Regions requires a registration form before download, so the GOaS shapefile cannot be downloaded automatically. You must obtain it manually.

1. **Download the shapefile**:
   - Visit: https://www.marineregions.org/downloads.php
   - Search for **GOaS_v1_20211214** and fill in the required registration form
   - Download `GOaS_v1_20211214.zip`

2. **Extract and setup**:
   ```bash
   # Create the directory
   mkdir -p shapefiles/GOaS_v1_20211214

   # Extract directly into the directory
   unzip GOaS_v1_20211214.zip -d shapefiles/GOaS_v1_20211214/
   ```

   The directory must contain: `World_Seas_IHO_v3.shp`, `.shx`, `.dbf`, `.prj`, `.cpg`

3. **Verify installation**:
   ```bash
   python -m boldgenotyper.goas_downloader
   ```

   Expected output: `[INFO] All GOaS files present`

### Custom Shapefiles for Non-Marine Organisms

BOLDGenotyper supports custom shapefiles for **freshwater and terrestrial organisms**, enabling the same biogeographic analysis capabilities across all environments.

**Supported Geographic Categories**:
- **Marine**: GOaS ocean basins (default — no flag needed)
- **Freshwater**: BasinATLAS or HydroBASINS watershed polygons
- **Terrestrial**: Ecoregions2017 biogeographic regions
- **Custom**: Any polygon shapefile

**Example: Freshwater organisms (Salmonidae)**
```bash
boldgenotyper data/Salmonidae.tsv \
  --build-tree \
  --custom-shp shapefiles/BasinATLAS/BasinATLAS_lev03.shp \
  --shp-field PFAF_ID \
  --geo-category drainage_basin
```

**Example: Terrestrial organisms (Pieridae butterflies)**
```bash
boldgenotyper data/Pieridae.tsv \
  --build-tree \
  --custom-shp shapefiles/Ecoregions2017/Ecoregions2017.shp \
  --shp-field REALM \
  --geo-category biogeographic_realm
```

**Parameters**:
- `--custom-shp`: Path to shapefile (.shp file)
- `--shp-field`: Attribute field containing region labels (see field guidance below)
- `--geo-category`: Label for the region type in outputs (e.g., `drainage_basin`, `biogeographic_realm`)

#### Choosing `--shp-field`: What BOLDGenotyper handles vs. what you decide

BOLDGenotyper performs the spatial assignment automatically given a shapefile and field name. For the three common use cases, recommended fields and the trade-offs between them are:

| Shapefile | Recommended `--shp-field` | Regions (N) | Alternative fields | Notes |
|-----------|--------------------------|-------------|-------------------|-------|
| GOaS (marine default) | *(handled internally)* | 10 ocean basins | — | No flag needed; used by default |
| BasinATLAS / HydroBASINS | `PFAF_ID` at level 03–04 | ~180–500 | `HYBAS_ID` (numeric IDs) | Choose level based on taxon range: level 03 for circumpolar taxa, level 04–05 for regional taxa |
| Ecoregions2017 | `REALM` | 8 biogeographic realms | `BIOME_NAME` (14 biome types); `ECO_NAME` (847 individual ecoregions) | `REALM` is best for family/order-level or cosmopolitan taxa; `BIOME_NAME` for habitat-type framing; `ECO_NAME` only for highly localized taxa |

**For your own shapefile**: Inspect the attribute table to identify region name fields:
```python
import geopandas as gpd
gdf = gpd.read_file("your_shapefile.shp")
print(gdf.columns.tolist())  # See available fields
print(gdf.head())             # Preview values
```

The choice of field (and aggregation level) has biological consequences — it determines how samples are grouped for the bar charts, pop-gen exports, and Arlequin populations. Refer to your shapefile's documentation for field definitions.

**Where to get shapefiles**:
- BasinATLAS (global freshwater basins): https://www.hydrosheds.org/products/basinsatlas
- HydroBASINS (by continent): https://www.hydrosheds.org/products/hydrobasins
- Ecoregions2017: https://ecoregions.appspot.com/
- Any compatible polygon shapefile with `.shp`, `.shx`, `.dbf`, and `.prj` components

### Skipping Geographic Analysis

If you don't need geographic analysis:

```bash
# Use the --no-geo flag to skip geographic analysis entirely
boldgenotyper data/my_species.tsv --no-geo
```

**Use cases for `--no-geo`**:
- Only interested in genotype identification and phylogeny
- Samples lack geographic coordinates
- No suitable shapefile available

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

### Step 4: Advanced Workflows (Optional)

**Optimise singleton filtering with parameter sweep, then apply the result**:
```bash
# Step 1: Find the optimal singleton distance threshold
boldgenotyper-sweep data/Sphyrna_lewini.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.03,0.05 \
  --output parameter_sweep/Sphyrna_lewini/

# Step 2: Read the recommended threshold from the output
cat parameter_sweep/Sphyrna_lewini/recommendations.txt

# Step 3: Re-run the full pipeline with the recommended threshold
boldgenotyper data/Sphyrna_lewini.tsv \
  --singleton-distance 0.015 \
  --build-tree \
  --output results/Sphyrna_lewini/
```

**Compare species vs. family analysis** for contamination detection:
```bash
# First run both analyses
boldgenotyper data/Sphyrna_lewini.tsv --output Sphyrna_lewini_output/
boldgenotyper data/Sphyrnidae.tsv --output Sphyrnidae_output/

# Then compare
python scripts/compare_analyses.py \
  --species-level Sphyrna_lewini_output/ \
  --family-level Sphyrnidae_output/ \
  --output comparative_analysis/
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

BOLDGenotyper runs a comprehensive 8-phase pipeline:

### Phase 1: Data Loading and Quality Control
- Parses BOLD TSV file and validates required columns
- **Marks coordinate quality** (identifies centroid, missing, or invalid coordinates)
  - Samples with centroid/missing coordinates are **retained for genotyping**
  - Quality markers used only for geographic analysis (not for filtering)
- Assigns ocean basins using GOaS shapefile (only to samples with high-quality coordinates)
  - Samples with centroid/missing coordinates marked as "Unknown" ocean basin
- Logs summary statistics (total samples, geographic quality samples, basin assignments)

### Phase 2: Haplotype Discovery
- **Uses ALL samples** (including those with centroid coordinates)
- Filters sequences by length (default: ≥200 bp) and N content (default: ≤10%)
- Aligns sequences with MAFFT (auto algorithm selection)
- Extracts core alignment region (trims ragged ends)
- Identifies Exact Sequence Variants (ESVs) as unique haplotypes
- Flags error singletons (singletons within 1 substitution of a larger haplotype)
- Flags suspect haplotypes based on pairwise distance from majority cluster

**Output**: Haplotype FASTA file with representative ESV sequences

### Phase 3: Haplotype Assignment
- **Uses ALL samples** (including those with centroid coordinates)
- Computes edit distance between each sample and all haplotype sequences
- Assigns samples to best-matching haplotype above similarity threshold (default: 50%)
- Flags ambiguous assignments (tie margin default: 0.003)
- Flags low-confidence assignments (below tie threshold: 0.95)
- Generates diagnostic CSV with identity scores and runner-up matches

**Output**: Annotated dataset with haplotype assignments and diagnostics (all samples included)

### Phase 4: Taxonomic Analysis
- Aggregates species names within each haplotype group
- Determines majority species for each haplotype
- Identifies taxonomy conflicts (haplotypes spanning multiple species)
- Generates species composition tables

**Output**: Taxonomy summaries and species-by-haplotype tables

### Phase 5: Phylogenetic Analysis (Optional, `--build-tree`)
- Aligns haplotype sequences with MAFFT
- Constructs maximum-likelihood tree with FastTree (GTR+Gamma model)
- Generates tree visualizations with tip labels
- Outputs Newick format for further analysis in tree editors (e.g., [TreeViewer](https://treeviewer.org/), FigTree, iTOL)

**Output**: Tree files (.nwk), relabeled tree with readable tip labels (_relabeled.nwk), visualizations (.png, .pdf)

### Phase 6: Metadata Analysis
- Analyzes associations between haplotype assignments and specimen metadata fields
- **Coverage analysis**: Assesses availability and completeness of each metadata field
- **Categorical analysis**: Examines haplotype distributions across metadata categories (sex, life stage, country, habitat, realm, biome, ecoregion)
- **Statistical testing**: Chi-square tests of haplotype–metadata associations with significance reporting
- **Temporal analysis** (if collection dates available):
  - Parses collection dates from various BOLD formats
  - Calculates haplotype emergence timelines
  - Analyzes temporal distribution patterns
- **Visualizations**: Coverage bar charts, categorical stacked bars, heatmaps, emergence timelines, species-faceted and country-faceted temporal plots
- Enabled by default; disable with `--no-metadata-analysis`

**Output**: Coverage JSON, field analysis CSVs, association test results, temporal summaries, multiple PNG/PDF plots

### Phase 7: Visualization
- **Identity Distribution**: Histogram of sequence identity scores for assigned samples
- **Geographic Distribution Maps**:
  - Global map with genotypes color-coded
  - Faceted maps (one per species or genotype)
  - Sample points sized by abundance
- **Region Abundance Bar Charts**:
  - Relative abundance (stacked, normalized by region)
  - Total abundance (stacked counts)
  - Faceted versions (one per species)
- All visualizations in PNG (high-res) and PDF (vector) formats
- Interactive visualizations in HTML report

**Output**: Multiple PNG/PDF files, JSON data files for interactive plots

### Plot Data Export (automatic)
- Exports raw data for each plot type as CSV files
- Creates `plot_config.yaml` with all styling parameters (colors, DPI, projections, dimensions)
- Generates Python regeneration scripts for each plot type
- Enables users to customize plots for publication without re-running the pipeline
- Disable with `--no-export-plot-data`

**Output**: `plots/` directory with data CSVs, YAML config, and Python scripts (see [Customizing Plots](#customizing-plots-for-publication))

### Phase 8: Report Generation
- Aggregates all results into interactive HTML report
- Summary statistics (total samples, genotypes, assignment rate, etc.)
- Detailed tables (assignment status, identity scores, taxonomy, geography)
- Methods section with parameters suitable for publication
- Interactive visualizations (filterable by genotype, basin, threshold)
- Export functionality (PNG, SVG, CSV)

**Output**: Comprehensive HTML report, assignment summary CSVs

### Advanced Workflows

Beyond the main pipeline, BOLDGenotyper provides three additional commands for specialized analyses:

#### Parameter Sweep (`boldgenotyper-sweep`)

Systematically tests multiple `--singleton-distance` thresholds to identify the optimal setting for your dataset. The sweep analyzes:
- Number of haplotypes vs. singleton threshold (elbow-point detection)
- Assignment coverage vs. threshold
- Sample haplotype stability across thresholds (group membership tracking)
- Identity score distributions

**When to use**: Before running your final analysis. The recommended workflow is: run the sweep → read `recommendations.txt` → pass the elbow-point threshold to `boldgenotyper --singleton-distance`.

**Output**: `recommendations.txt` (threshold + rationale), `elbow_plot.pdf`, `threshold_stability.pdf`, `sweep_summary.csv`, `group_membership_tracking.csv`

See `parameter_sweep/README.md` for interpretation guide.

#### Comparative Analysis (`scripts/compare_analyses.py`)

Compares species-level and family-level analyses to detect potential contamination, mislabeling, or cryptic species. Identifies samples that cluster with unexpected species at the family level. This is a standalone utility script, not a package subcommand.

**When to use**: Quality control for published datasets, contamination detection, multi-species studies

**Usage**:
```bash
python scripts/compare_analyses.py \
  --species-level Sphyrna_lewini_output/ \
  --family-level Sphyrnidae_output/ \
  --generate-reassignment-table
```

**Output**: Comparison metrics, genotype crosswalk table, sample reassignment recommendations, publication-ready methods text

See [`docs/guides/COMPARATIVE_ANALYSIS_GUIDE.md`](docs/guides/COMPARATIVE_ANALYSIS_GUIDE.md) for detailed interpretation.

---

## Usage Guide

### Basic Command Structure

```bash
boldgenotyper <input_tsv> [options]
```

### Complete Example

```bash
# Full pipeline with all options
# Run boldgenotyper-sweep first to determine --singleton-distance for your dataset
boldgenotyper data/Sphyrna_lewini.tsv \
  --organism "Sphyrna_lewini" \
  --output results/Sphyrna_analysis \
  --similarity-threshold 0.50 \
  --tie-margin 0.003 \
  --tie-threshold 0.95 \
  --singleton-distance 0.015 \
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

**3. Highly Diverse Taxon (Lower Assignment Threshold)**:
```bash
# Use lower similarity threshold for species complexes with higher divergence
boldgenotyper data/Carcharhinus_complex.tsv \
  --similarity-threshold 0.80 \
  --threads 16
```

**4. Fine-Scale Population Study (Strict Assignment)**:
```bash
# Use stringent similarity threshold for within-species genotyping
boldgenotyper data/Population_samples.tsv \
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

### Querying New Sequences

After completing a BOLD analysis, you can query new COI sequences against the identified haplotypes without re-running the full pipeline. This is useful for:
- Assigning new specimens to existing haplotypes
- Validating reference genomes against field samples
- Cross-study comparison and longitudinal monitoring

**Basic Query (Single or Multi-FASTA)**:
```bash
# Query new sequence(s) against identified haplotypes
boldgenotyper-query \
  --query new_samples.fasta \
  --haplotypes results/Sphyrna_analysis/haplotypes/Sphyrna_haplotypes.fasta \
  --output query_results/
```

**Query with Metadata Enrichment**:
```bash
# Include species composition and haplotype metadata
boldgenotyper-query \
  --query new_sample.fasta \
  --haplotypes results/Sphyrna_analysis/haplotypes/Sphyrna_haplotypes.fasta \
  --analysis-dir results/Sphyrna_analysis/ \
  --output query_results/
```

**Customize Top Matches and Length Filters**:
```bash
# Report top 20 matches, filter by sequence length
boldgenotyper-query \
  --query samples.fasta \
  --haplotypes haplotypes.fasta \
  --top-n 20 \
  --min-length 150 \
  --max-length 1500 \
  --output results/
```

**Query Output Files**:
All three formats are generated automatically:
- `query_results.csv`: Machine-readable table with identity metrics
- `query_results.json`: Structured data for programmatic access
- `query_results_detailed.txt`: Human-readable report with alignments

**Match Quality Classification**:
- **perfect** (100%): Exact haplotype match
- **high** (≥99.5%): Likely same haplotype, minor sequencing variation
- **good** (≥97%): Same species, possibly different haplotype
- **moderate** (≥95%): Same genus, divergent haplotype
- **low** (<95%): Different species or contamination

**Important Notes**:
- Query sequences must be **pre-extracted COI regions** (not multi-gene genomes)
- Multi-FASTA files are supported for batch processing
- No identity threshold - all matches are reported and ranked by quality
- Local alignment handles length differences between query and consensus

### Input File Requirements

**Filename Convention**: Your input TSV should follow this pattern for automatic organism detection:
```
Genus_species.tsv
```

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
| `--similarity-threshold` | Float | 0.50 | Minimum identity for haplotype assignment (0-1) |
| `--tie-margin` | Float | 0.003 | Maximum identity difference to call a tie (0-1) |
| `--tie-threshold` | Float | 0.95 | Minimum identity to consider tie detection (0-1) |
| `--singleton-distance` | Float | 0.005 | Min divergence to retain a singleton haplotype (0-1); use `boldgenotyper-sweep` to optimise |
| `--threads` | Integer | 4 | Number of parallel CPU threads |
| `--build-tree` | Flag | False | Build phylogenetic tree with FastTree |
| `--no-report` | Flag | False | Skip HTML report generation |
| `--no-geo` | Flag | False | Skip geographic analysis |
| `--custom-shp` | Path | None | Custom shapefile for geographic region assignment |
| `--shp-field` | String | name | Shapefile attribute field containing region labels |
| `--geo-category` | String | Auto | Category name for geographic regions (auto-detected) |
| `--no-metadata-analysis` | Flag | False | Skip metadata analysis |
| `--metadata-fields` | List | All | Metadata fields to analyze (space-separated) |
| `--no-normalize-sex` | Flag | False | Do not normalize sex values |
| `--no-temporal-analysis` | Flag | False | Skip temporal analysis of collection dates |
| `--no-export-plot-data` | Flag | False | Skip exporting plot data and regeneration scripts |
| `--log-level` | String | INFO | Logging verbosity (DEBUG, INFO, WARNING, ERROR) |
| `--version` | Flag | - | Show version and exit |

### Parameter Details

#### `--similarity-threshold` (default: 0.50)

**What it controls**: Minimum sequence identity required for a sample to be assigned to a haplotype.

**Interpretation**:
- Samples below this threshold are marked as "unassigned"
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

#### `--singleton-distance` (default: 0.005)

**What it controls**: The minimum divergence a singleton haplotype must have from its nearest neighbour to be retained. Singletons below this threshold are removed as likely sequencing or PCR errors.

**Background**: Because ESV discovery treats every unique sequence as a distinct haplotype, error-derived sequences that differ by 1–2 bp from a true haplotype appear as spurious singletons. The default (0.5%) removes most such errors while retaining biologically real rare variants.

**The recommended workflow is to run `boldgenotyper-sweep` first and use the reported elbow-point value here**:
```bash
# Step 1 — find the optimal threshold for your dataset
boldgenotyper-sweep data/MySpecies.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.03,0.05 \
  --output sweep_results/

# Step 2 — read the recommendation
cat sweep_results/recommendations.txt   # e.g. "PRIMARY: 0.015"

# Step 3 — apply it to the full pipeline run
boldgenotyper data/MySpecies.tsv \
  --singleton-distance 0.015 \
  --build-tree
```

**When to adjust**:
- **Lower (0.001–0.005)**: Retain more rare variants; risk including sequencing errors
- **Higher (0.02–0.10)**: Aggressively remove singletons; risk losing biologically real rare haplotypes
- Values >0.05 are rarely needed and should be supported by the sweep elbow plot

**Taxon-specific guidance**: Optimal values tend to be consistent within major taxonomic groups (marine vertebrates ≈ 0.015; insect families may require up to 0.085). Always validate with `boldgenotyper-sweep` for new datasets.

#### `--threads`

**What it controls**: Number of CPU cores used for parallel processing.

**Parallelized steps**:
- Genotype assignment (edit distance calculations)
- MAFFT alignment
- FastTree construction

**Recommendation**: Set to number of available cores (check with `nproc` or `sysctl -n hw.ncpu`)

### Advanced Command Parameters

#### `boldgenotyper-sweep` Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tsv` | Path | Required | Input BOLD TSV file |
| `--thresholds` | String | 0.01,0.015,0.02,0.03,0.05 | Comma-separated threshold values to test |
| `--output` | Path | parameter_sweep/ | Output directory |
| `--threads` | Integer | 4 | Number of parallel threads (sequential by default) |
| `--keep-intermediates` | Flag | False | Keep full output for each threshold tested |
| `--log-level` | String | INFO | Logging verbosity |

**Example**:
```bash
boldgenotyper-sweep data/samples.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.025,0.03 \
  --threads 8 \
  --output sweep_results/
```

#### `scripts/compare_analyses.py` Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--species-level` | Path | Required | Path to species-level analysis directory or annotated CSV |
| `--family-level` | Path | Required | Path to family-level analysis directory or annotated CSV |
| `--output` | Path | comparative_analysis/ | Output directory |
| `--generate-reassignment-table` | Flag | False | Generate detailed sample-level reassignment table |
| `--majority-threshold` | Float | 0.7 | Threshold for majority species assignment |
| `--log-level` | String | INFO | Logging verbosity |

**Example**:
```bash
python scripts/compare_analyses.py \
  --species-level Sphyrna_lewini_output/ \
  --family-level Sphyrnidae_output/ \
  --generate-reassignment-table \
  --output comparison_results/
```

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
│   ├── {organism}_tree_relabeled.nwk     # Tree with readable tip labels (for tree editors like TreeViewer)
│   ├── {organism}_tree.png               # Tree visualization
│   └── {organism}_tree.pdf
│   #	Open _relabeled.nwk files in tree editors such as TreeViewer (https://treeviewer.org/),
│   #   FigTree, or iTOL for re-rooting, customization, and publication-quality figures
│
├── metadata_analysis/                           # Metadata analysis outputs
│   ├── {organism}_metadata_coverage.json        # Field availability statistics
│   ├── {organism}_metadata_summary.csv          # Summary of all fields
│   ├── {organism}_{field}_analysis.csv          # Per-field haplotype distributions
│   ├── {organism}_association_tests.csv         # Chi-square test results
│   └── {organism}_temporal_*.csv                # Temporal analysis (if dates available)
│
├── visualization/
│   ├── {organism}_identity_distribution.*       # Identity histogram
│   ├── {organism}_distribution_map.*            # Global distribution map
│   ├── {organism}_distribution_map_faceted.*    # Faceted by species
│   ├── {organism}_distribution_bar.*            # Relative abundance by region
│   ├── {organism}_distribution_bar_faceted.*    # Faceted bar charts
│   ├── {organism}_totaldistribution_bar.*       # Total counts by region
│   ├── metadata/                                # Metadata analysis plots
│   │   ├── {organism}_metadata_coverage.*       # Coverage bar chart
│   │   ├── {organism}_{field}_by_haplotype.*    # Categorical stacked bars
│   │   ├── {organism}_metadata_heatmap.*        # Haplotype × field heatmap
│   │   └── {organism}_haplotype_emergence.*     # Temporal emergence timeline
│   ├── *_data.json                              # Data for interactive HTML plots
│   └── {organism}_tree.* (if --build-tree)
│
├── plots/                                       # Plot customization kit
│   ├── plot_config.yaml                         # Editable styling parameters
│   ├── README.md                                # Regeneration instructions
│   ├── data/                                    # Raw plot data (CSV)
│   │   ├── distribution_map.csv
│   │   ├── distribution_bar_relative.csv
│   │   ├── distribution_bar_absolute.csv
│   │   ├── identity_distribution.csv
│   │   └── genotype_colors.csv
│   └── scripts/                                 # Python regeneration scripts
│       ├── regenerate_all.py
│       ├── regenerate_map.py
│       ├── regenerate_bars.py
│       └── regenerate_identity.py
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
- `ocean_basin`: Assigned ocean basin (e.g., "North Atlantic Ocean" or "Unknown")
- `species`: Original BOLD species name
- **Coordinate quality markers** (new):
  - `has_centroid_coords`: True if coordinates are country/region centroids
  - `has_missing_coords`: True if latitude or longitude is missing
  - `has_zero_coords`: True if coordinates are [0, 0]
  - `has_invalid_coords`: True if coordinates are out of valid range
  - `is_geographic_quality`: True if coordinates suitable for ocean basin assignment
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

**Newick Tree Files** (`.nwk`): Phylogenetic trees in standard Newick format can be found in `phylogenetic/{organism}_tree_relabeled.nwk`. These files can be opened in any tree editor of your choice for re-rooting, annotation, and customization:
- **[TreeViewer](https://treeviewer.org/)** (Bianchini & Sánchez-Baracaldo, 2024) - Modern, flexible, modular tree visualization software
- **FigTree** - Popular tree viewer with extensive annotation options
- **iTOL** (Interactive Tree of Life) - Web-based tree display and annotation
- **R (ape/ggtree)** - For programmatic tree manipulation and publication figures

#### 5. `{organism}_pipeline.log`

**Complete pipeline log** with timestamps, including:
- Input validation and data loading statistics
- Dereplication metrics (haplotypes identified, singletons filtered)
- Assignment summary (success rate, ties, low-confidence)
- Geographic analysis results (basin assignments)
- Warnings and errors encountered
- Runtime for each phase

---

## Customizing Plots for Publication

BOLDGenotyper exports all plot data and Python scripts for easy customization without re-running the entire analysis.

**📚 Additional Resources:**
- **Auto-generated guide**: After running the pipeline, see `{organism}_output/plots/README.md` for analysis-specific instructions
- **Interactive tutorial**: See `notebooks/06_plot_customization_tutorial.ipynb` for detailed examples, common errors, and troubleshooting

### Plot Regeneration Directory

After running the pipeline, you'll find:

```
{organism}_output/plots/
├── plot_config.yaml          # Customization settings (colors, sizes, formats)
├── README.md                 # Quick start guide
├── data/                     # Exported plot data (CSV files)
│   ├── distribution_map.csv
│   ├── distribution_bar_relative.csv
│   ├── distribution_bar_absolute.csv
│   ├── identity_distribution.csv
│   └── genotype_colors.csv
└── scripts/                  # Python regeneration scripts
    ├── regenerate_all.py     # Run all plots at once
    ├── regenerate_map.py     # Regenerate distribution map
    ├── regenerate_bars.py    # Regenerate bar charts
    └── regenerate_identity.py # Regenerate identity histogram
```

### Quick Start: Customize Your Plots

**1. Edit the configuration file:**

```bash
cd {organism}_output/plots
nano plot_config.yaml  # or use any text editor
```

**2. Modify settings:**

```yaml
# Example customizations
general:
  output_format: ['pdf', 'png']
  dpi: 600  # High resolution for publication
  width_inches: 12
  height_inches: 8

colors:
  Consensus_1_S._lewini: "#E41A1C"    # Custom hex color
  Consensus_2_S._lewini: "#377EB8"
  Consensus_3_S._lewini: "#4DAF4A"

filters:
  include_genotypes: []              # Leave empty to include all
  exclude_genotypes: ["Consensus_10"] # Hide specific genotypes

map:
  projection: "mollweide"             # robinson, mollweide, mercator
  center_longitude: -180              # Center on Pacific Ocean
  point_alpha: 0.8
  point_size_range: [3, 12]
  legend_position: "right"

bars:
  orientation: "horizontal"           # vertical or horizontal
  bar_width: 0.7

identity:
  binwidth: 1.0
  show_mean: true
  show_median: true
  x_limits: [95, 100]                # Focus on high-identity range
```

**3. Regenerate plots:**

```bash
# Regenerate all plots with new settings
python scripts/regenerate_all.py

# Or regenerate individual plots
python scripts/regenerate_map.py
python scripts/regenerate_bars.py
python scripts/regenerate_identity.py
```

**4. Find updated plots:**

Regenerated plots are saved to `{organism}_output/visualization/` with `_custom` suffix:
- `distribution_map_custom.pdf`
- `distribution_bar_relative_custom.pdf`
- `identity_distribution_custom.pdf`

### Customization Options

**Colors**: Specify hex colors for each genotype
**Figure size**: Control width, height, DPI for different journals
**Map projection**: Robinson, Mollweide, Mercator, PlateCarree
**Filters**: Include or exclude specific genotypes from plots
**Bar orientation**: Vertical or horizontal stacking
**Histogram binning**: Adjust bin width and statistical overlays

### Requirements

All required packages are in the boldgenotyper environment - no additional installation needed!
- matplotlib, seaborn (plotting)
- geopandas, cartopy (maps)
- pandas (data handling)
- pyyaml (configuration)

### Advanced Customization and Troubleshooting

For more detailed guidance including:
- Common YAML syntax errors and solutions
- 5+ example configurations (publication-quality, Pacific-focused, filtered genotypes, etc.)
- Complete parameter reference for all settings
- Three methods to regenerate plots (CLI, Python module, individual functions)
- Troubleshooting guide

See the **interactive tutorial**: `notebooks/06_plot_customization_tutorial.ipynb`

---

### Advanced Command Outputs

#### Parameter Sweep Output (`parameter_sweep/`)

```
parameter_sweep/
├── sweep_summary.csv              # Metrics for all tested thresholds
├── threshold_stability.pdf        # Multi-panel stability visualization
├── elbow_plot.pdf                # Optimal threshold detection plot
├── group_membership_tracking.csv  # Sample clustering stability across thresholds
├── recommendations.txt            # Detailed recommendations and rationale
├── README.md                      # Interpretation guide
└── runs/                          # Full outputs for each threshold (if --keep-intermediates)
    ├── threshold_0_010/
    ├── threshold_0_015/
    └── ...
```

**Key files**:
- `sweep_summary.csv`: Number of groups, assignment rates, mean identity for each threshold
- `elbow_plot.pdf`: Visualizes elbow point for optimal threshold selection
- `group_membership_tracking.csv`: Tracks which samples cluster together across thresholds (uses Jaccard similarity, not arbitrary group names)
- `recommendations.txt`: Data-driven threshold recommendation with biological interpretation

See `parameter_sweep/README.md` for detailed interpretation guide.

#### Comparative Analysis Output (`comparative_analysis/`)

```
comparative_analysis/
├── comparison_summary.csv         # High-level comparison metrics
├── genotype_crosswalk.csv        # Mapping of species groups to family groups
├── sample_reassignments.csv      # Sample-level reassignment recommendations
├── methods_text.md               # Publication-ready methods section
└── README.md                     # Interpretation guide
```

**Key files**:
- `comparison_summary.csv`: Total samples, consensus groups, species detected, identity metrics, contamination signals
- `genotype_crosswalk.csv`: Shows how species-level genotypes map to family-level genotypes
- `sample_reassignments.csv` (if `--generate-reassignment-table`): Detailed sample-level analysis identifying potential contamination
- `methods_text.md`: Ready-to-paste methods text for publications

**Interpretation**: Samples that cluster with unexpected species at the family level may represent contamination, mislabeling, or introgression.


#### Population Genetics Exports (`exports/`)

When using `--export-format`, genotypes are exported to formats compatible with population genetics software:

```
exports/
├── README.md                      # Format descriptions and usage
├── arlequin/
│   ├── {organism}.arp            # Arlequin project file
│   └── populations.txt           # Population definitions
├── popart/
│   ├── {organism}.nexus          # NEXUS alignment with traits
│   └── populations.csv           # Population mapping
├── dnasp/
│   └── {organism}.fas            # FASTA with population labels
└── generic/
    ├── alignment.fasta           # Consensus sequences
    ├── genotype_membership.csv   # Sample-to-genotype mapping
    └── haplotypes.csv            # Haplotype summary table
```

See [`docs/guides/POPGEN_EXPORT_GUIDE.md`](docs/guides/POPGEN_EXPORT_GUIDE.md) for detailed format specifications and usage examples.

---

## Biological Context for Threshold Selection

### Understanding Haplotype Discovery (ESV Approach)

BOLDGenotyper uses the **Exact Sequence Variant (ESV)** approach (Porter & Hajibabaei 2020) for haplotype discovery. Each unique sequence from the core alignment region is treated as a distinct haplotype — no clustering threshold is required.

**COI as a Molecular Marker**:
- COI is typically conserved within species (>97% identity)
- Can vary between populations (1-5% divergence)
- Varies significantly between species (>5-10% divergence)

**Error Filtering and `--singleton-distance`**: Singleton ESVs that diverge from their nearest neighbour by less than `--singleton-distance` (default: 0.5%) are removed as probable sequencing or PCR errors. This threshold is the single most impactful parameter for controlling haplotype count. **Always run `boldgenotyper-sweep` first to identify the elbow-point threshold for your dataset**, then pass that value to `boldgenotyper --singleton-distance`. The sweep plots haplotype count vs. threshold and identifies the inflection point at which further increasing the threshold stops meaningfully reducing haplotype count — this is the recommended operating point. Optimal thresholds vary by taxon and sequence quality (empirically: marine vertebrates ≈ 0.015; some insect families ≈ 0.085).

### Understanding Similarity Threshold

The **similarity threshold** determines the minimum identity required for a sample to be assigned to a haplotype.

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

| Taxonomic Group | Similarity Threshold | Notes |
|----------------|---------------------|-------|
| **Elasmobranchs (sharks, rays)** | 0.80-0.90 | Slow molecular evolution, lower divergence |
| **Teleost fish** | 0.80-0.90 | Standard barcoding thresholds |
| **Marine invertebrates** | 0.70-0.85 | Higher COI variation |
| **Cephalopods (squid, octopus)** | 0.80-0.90 | Similar to fish |
| **Crustaceans** | 0.70-0.85 | Often more divergent |
| **Marine mammals** | 0.85-0.95 | Low COI variation |

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
- **Pipeline behavior**: **Marks** (not excludes) samples with:
  - Missing lat/lon (NaN, empty)
  - Centroid coordinates (coord_source contains "centroid")
  - Invalid coordinates ([0, 0] or out of range)
- **Key Distinction**:
  - Samples with centroid/missing coordinates **ARE INCLUDED** in clustering and genotyping
  - These samples **ARE EXCLUDED** from ocean basin assignment (marked as "Unknown")
  - Users can see which samples have coordinate quality issues via marker columns
- **Rationale**: A centroid for Mexico could fall in Pacific or Atlantic, but the sequence is still valid for genotyping
- **Implication**: Maximizes samples for genotype identification while maintaining geographic accuracy
- **Quality Markers**: Added columns include `has_centroid_coords`, `has_missing_coords`, `is_geographic_quality`

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

#### 1. Custom Geographic Regions Supported
- **Default**: Geographic analysis uses GOaS (Global Oceans and Seas) shapefile for marine taxa
- **Custom shapefiles**: Use `--custom-shp` to apply any polygon shapefile (freshwater basins, ecoregions, biomes, etc.)
- **Examples**: HydroBASINS for freshwater organisms, WWF ecoregions for terrestrial species, custom study area boundaries
- **No geographic analysis needed**: Use `--no-geo` flag to skip geographic analysis entirely

#### 2. Coordinate Precision Requirements
- **Requirement**: Geographic analysis requires precise GPS coordinates
- **Problem**: Many BOLD records have country-level or imprecise coordinates
- **Impact**: Samples with low-quality coordinates are:
  - **Included in genotyping and clustering** (sequences are still valid)
  - **Marked as "Unknown" for ocean basin** (geographic location uncertain)
  - Flagged with quality markers (`has_centroid_coords`, `has_missing_coords`)
- **Benefit**: Maximizes samples for genotype identification while maintaining geographic accuracy
- **Quantification**: Check HTML report for samples with "Unknown" ocean basin and coordinate quality markers

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
- **Workaround**: Use `--export-format all` to export directly to Arlequin, DnaSP, and PopART formats

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

**Top feature requests** (as of v1.0.0):
1. Haplotype network construction (TCS, median-joining)
2. Haplotype network construction
3. Population genetics statistics integration
4. Multi-locus support
5. Docker containerization
6. Galaxy tool integration
7. Web interface

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
# Solution 1: Download manually from Marine Regions (requires registration form)
# Visit: https://www.marineregions.org/downloads.php
# Download: GOaS_v1_20211214.zip
# Extract to: shapefiles/GOaS_v1_20211214/
# Then verify:
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

**Problem**: `No haplotypes identified`

```bash
# Check log file for details
tail -100 {organism}_output/{organism}_pipeline.log

# Possible causes:
# 1. All sequences too short (<200bp) or too many Ns (>10%)
# 2. MAFFT alignment failed
# 3. No sequences passed quality control

# Check sequence quality
grep -c ">" data/my_species.fasta  # (if FASTA exists)
```

**Problem**: `Most samples unassigned`

```bash
# Check diagnostics file
less {organism}_output/genotype_assignments/{organism}_diagnostics.csv

# Look at identity scores - are they all low?
cut -d',' -f3 {organism}_output/genotype_assignments/{organism}_diagnostics.csv | sort -n

# If identities are 0.4-0.6, haplotypes may not represent data well
# Solutions:
# 1. Lower similarity threshold (e.g., --similarity-threshold 0.40)
# 2. Check if sequences are actually from the same marker/region
# 3. Run parameter sweep to assess haplotype structure
```

**Problem**: `FastTree failed` (if using `--build-tree`)

```bash
# Verify FastTree is installed
fasttree 2>&1 | head -1

# Install if missing
conda install -c bioconda fasttree

# Check that haplotype sequences exist
ls {organism}_output/haplotypes/

# Try running pipeline without tree first
boldgenotyper data/my_species.tsv  # omit --build-tree
```

### Data Quality Issues

**Problem**: Too many haplotypes (many singletons)

```bash
# Check the log for singleton filtering statistics
tail -50 {organism}_output/{organism}_pipeline.log

# Run parameter sweep to find the optimal threshold
boldgenotyper-sweep data/my_species.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.03,0.05 \
  --output parameter_sweep/

# Read the recommended threshold, then re-run with it applied
cat parameter_sweep/recommendations.txt
boldgenotyper data/my_species.tsv --singleton-distance 0.015
```

**Problem**: Too few haplotypes (expected more diversity)

```bash
# Inspect haplotype sequences and suspect flags
ls {organism}_output/haplotypes/

# Check if suspect haplotypes were filtered; review log
grep -i "suspect\|flagged" {organism}_output/{organism}_pipeline.log
```

**Problem**: Many samples flagged as ties

```bash
# Option 1: Accept ties as real biological signal (genotypes not clearly distinct)
# Option 2: Increase tie-margin to reduce tie flagging
boldgenotyper data/my_species.tsv --tie-margin 0.001

# Option 3: Manual review - inspect diagnostics file
```

**Problem**: Many samples marked as "Unknown" ocean basin

```bash
# Check how many samples have "Unknown" basin
grep -c "Unknown" {organism}_output/{organism}_annotated.csv

# If many "Unknown", check coordinate quality markers in annotated CSV:
# - has_centroid_coords: True if coordinates are centroids (common in BOLD)
# - has_missing_coords: True if lat/lon is missing
# - has_zero_coords: True if coordinates are [0, 0]
# - is_geographic_quality: False means one of the above issues

# Note: Samples with "Unknown" basin are STILL INCLUDED in genotyping!
# They're only excluded from geographic distribution analysis.

# Other potential causes:
# 1. Sample coordinates in ocean (correctly assigned)
# 2. Coordinates on land (no basin to assign)
# 3. GOaS shapefile not loaded correctly
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

### Advanced Command Troubleshooting

**Problem**: `boldgenotyper-sweep` - How do I interpret the elbow plot?

```bash
# The elbow plot shows metrics vs. threshold
# Look for the "elbow point" where curves level off
# This indicates where adding more groups provides diminishing returns

# Check the recommendations file:
cat parameter_sweep/recommendations.txt

# Examine stability plot:
open parameter_sweep/threshold_stability.pdf
```

**Problem**: `boldgenotyper-sweep` - Samples have high `n_changes` in group membership tracking

```bash
# High n_changes indicates unstable clustering
# These samples may be:
# 1. Intermediates between genotypes
# 2. Low-quality sequences
# 3. Potential contamination

# Review samples with n_changes > 2:
# Sort by n_changes descending
head -20 parameter_sweep/group_membership_tracking.csv

# Consider excluding or manually reviewing these samples
```

**Problem**: `compare_analyses.py` - No contamination detected but I expected some

```bash
# Check that you're comparing the right datasets:
# 1. Species-level should be a subset of family-level
# 2. Both analyses should use same or similar parameters

# Review the comparison summary:
cat comparative_analysis/comparison_summary.csv

# Check if samples are genuinely clean or majority-threshold is too permissive
```

**Problem**: Custom shapefile not working

```bash
# Verify shapefile has all required components (.shp, .shx, .dbf, .prj)
ls -la path/to/shapefile.*

# Check shapefile field names:
# Use QGIS, ogrinfo, or GeoPandas to inspect
ogrinfo -al -so path/to/shapefile.shp | grep -A 20 "Layer name"

# Verify the field name matches your --shp-field argument
# Field names are case-sensitive

# Check coordinate reference system (CRS):
# Shapefiles should use WGS84 (EPSG:4326) or will be reprojected automatically

# Verify column names match exactly (case-sensitive)
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
- **Haplotypes**: 26 ESV haplotypes identified
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
  version = {1.0.0},
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

- **TreeViewer** (optional, for tree customization): Bianchini, G. & Sánchez-Baracaldo, P. (2024). TreeViewer: Flexible, modular software to visualise and manipulate phylogenetic trees. *Ecology and Evolution*, 14, e10873. [doi:10.1002/ece3.10873](https://doi.org/10.1002/ece3.10873)

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
- **Email**: Steph Smith (symbioseas@outlook.com)

### Frequently Asked Questions

**Q: Can I use BOLDGenotyper for non-COI markers?**
A: Yes. The ESV-based haplotype discovery requires no clustering threshold and works for any aligned marker. You may need to adjust the similarity threshold (`--similarity-threshold`) for assignment based on the expected divergence range of your marker.

**Q: How do I analyze terrestrial or freshwater organisms?**
A: Use `--custom-shp` to supply any polygon shapefile for geographic analysis (e.g., HydroBASINS for freshwater, WWF Ecoregions for terrestrial). Use `--no-geo` if you don't need geographic analysis at all.

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

### v1.0.0 (2025) - Publication Release

**New in v1.0.0**:
- Outgroup rooting: `--phylo-outgroup-fasta`, `--phylo-outgroup-label`, `--phylo-outgroup-taxon`
- Custom shapefile support: `--custom-shp` for freshwater basins, ecoregions, or any polygon shapefile
- Metadata analysis module: haplotype-metadata associations, chi-square tests, temporal analysis
- Population genetics export: Arlequin, PopART/NEXUS, DnaSP, generic CSV/FASTA (`--export-format`)
- Plot customization: YAML-configured regeneration from exported data
- Haplotype query tool: `boldgenotyper-query` for assigning new sequences to existing haplotypes
- Parameter sweep tool: `boldgenotyper-sweep` for optimizing singleton filtering thresholds

**Core Features**:
- Complete unified pipeline from BOLD TSV to annotated results
- Automated haplotype discovery using Exact Sequence Variants (ESVs)
- Edit distance-based haplotype assignment with tie detection and confidence scoring
- Geographic analysis with any polygon shapefile (GOaS ocean basins, HydroBASINS, Ecoregions2017, etc.)
- Phylogenetic tree building with FastTree (GTR+Gamma), outgroup rooting, MSA visualization
- Comprehensive interactive HTML reports with Plotly.js
- Publication-ready visualizations (PDF/SVG)

**Known Limitations**:
- No haplotype network construction (TCS, median-joining)
- No population genetics statistics (Fst, AMOVA, diversity indices)
- Single-locus only (COI)

---

**Planned Features**:
1. Haplotype network construction (TCS, median-joining)
2. Population genetics statistics (Fst, AMOVA, diversity indices)
3. Multi-locus support for concatenated/combined markers
4. Docker container for reproducibility and portability

---

**Built with ❤️ for open science and marine conservation**

**Created by**: Steph Smith ([@SymbioSeas](https://github.com/SymbioSeas))

**Questions or feedback?**
- GitHub: https://github.com/SymbioSeas/BOLDGenotyper
- Email: symbioseas@outlook.com
- Issues: https://github.com/SymbioSeas/BOLDGenotyper/issues
