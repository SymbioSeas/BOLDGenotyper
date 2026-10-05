# BOLDGenotyper

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://zenodo.org/)
[![Version](https://img.shields.io/badge/version-2.0.0-green.svg)](https://github.com/SymbioSeas/BOLDGenotyper/releases)

**Automated COI sequence genotyping and biogeographic analysis from BOLD database data**

BOLDGenotyper is a bioinformatics pipeline that identifies COI (cytochrome c oxidase subunit I) haplotypes in Barcode of Life Data System (BOLD) records for any taxonomic group. It performs quality control, haplotype discovery as exact sequence variants of a core alignment region, record-level assignment with an explicit status for every record, taxonomic labeling with label-conflict evidence, optional phylogenetic analysis, region assignment from geographic-quality coordinates, and figures of the haplotypes' distribution.

**Upgrading from 1.0.0?** Output files, columns and some flags changed in 2.0.0. See the [migration guide](docs/guides/MIGRATION_v1_to_v2.md).

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

- **One command** from a BOLD TSV to a documented run folder: `records.csv` (one row per record), `haplotypes.csv` (one row per haplotype), `haplotypes.fasta`, a summary, and a generated `README.md` that defines every file and column
- **Haplotype discovery**: exact sequence variants (ESVs) of a core alignment region, keyed by a stable sequence hash (`seq_md5`) and numbered `H0001`, `H0002`, ...
- **A status for every record**: assigned, assigned after error correction, excluded as an error singleton, excluded for core coverage, failed QC, or no sequence
- **Quality control**: orientation and reading-frame checks with a configurable genetic code, length and ambiguous-base filters, and an error filter for singleton haplotypes (`--singleton-distance`)
- **Labels with evidence**: each haplotype's member labels, majority label and fraction, and every label conflict with the nearest haplotype carrying that label
- **Geographic analysis**: region assignment from geographic-quality coordinates using GOaS ocean basins (marine), HydroBASINS (freshwater), Ecoregions2017 (terrestrial), or any polygon shapefile; centroid coordinates are kept but not assigned a region
- **Phylogeny**: optional FastTree tree (GTR+Gamma) of every non-suspect haplotype, with excluded tips listed
- **Reproducibility**: `run_provenance.json` records the full configuration, tool versions, seeds and input checksum; the log ends with a warnings summary
- **HTML report** built from the run folder: key findings, what needs attention (haplotypes left at genus with their evidence, label conflicts, curated taxa), figures with species and genus menus, the methods with a ready-to-edit statement, and every output file
- **Optional extras** with `--extended` (metadata cross-tabulations, collection-date summaries, a label-based barcoding gap, species tables, alignment plots, region and metadata figures by species), `--export-format` (Arlequin, PopART, DnaSP, generic) and `--plot-kit` (plot data and regeneration scripts)
- **Query tool** (`boldgenotyper-query`) to place new sequences on a run's haplotypes, and a **parameter sweep** (`boldgenotyper-sweep`) for `--singleton-distance`
- **Skip geography** with `--no-geo` when no shapefile is available

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
- **dereplication.py**: MAFFT alignment, core region extraction, ESV identification, error-singleton filtering and correction, suspect flagging
- **haplotype_ids.py**: `seq_md5` keys and `H0001` display IDs
- **record_status.py**: one status per input record
- **taxon_labels.py**: label levels, pooled member labels, species counts

### Layer 4: Analysis
- **divergence_analysis.py**: Pairwise p-distance matrices, label-based barcoding gap, species-level divergence
- **species_analysis.py**: Species-level aggregation, diversity metrics, geographic summaries
- **parameter_sweep.py**: Tests multiple singleton thresholds, tracks group membership stability

### Layer 5: Geographic & Phylogenetic
- **geographic.py**: GOaS shapefile loading, point-in-polygon ocean basin assignment
- **geographic_enhancement.py**: Coverage counts, region-assignment confidence
- **goas_downloader.py**: Verifies Global Oceans and Seas reference shapefile installation
- **phylogenetics.py**: MAFFT alignment, trimAl, FastTree tree building, outgroup rooting
- **msa_visualization.py**: Phylogeny-ordered MSA display with nucleotide coloring

### Layer 6: Visualization & Reporting
- **tree_figures.py**, **run_figures.py**, **metadata_figures.py**: the run's figures (tree, maps, regions, metadata views), an overview plus one page per species
- **figure_pages.py**: color capping, per-species pages and figure manifests; **visualization.py**: shared palette and map helpers
- **metadata_summary.py**: Metadata coverage, normalized cross-tabulations and collection-date summaries
- **reports.py**: Assignment summary (`summary.csv`); **report/**: the HTML report, one module per section
- **layout.py**: `records.csv`, `haplotypes.csv` and the run-folder file map
- **run_readme.py**: the run folder's `README.md`
- **provenance.py**: `run_provenance.json` and end-of-run warnings
- **plot_export.py**: Exports plot data CSVs + config YAML + Python regeneration scripts
- **plot_regeneration.py**: Regenerates plots from exported data using customized config

### Layer 7: Export
- **popgen_export.py**: Exports to Arlequin, PopART/NEXUS, DnaSP, generic CSV
- **haplotype_query.py**, **query_assignment.py**: Place query sequences on a run's haplotypes

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
- Alignment plots: pymsaviz (multiple sequence alignment figures with `--extended`)
- Utilities: pyyaml, jinja2

**No additional installation needed!** All dependencies are managed by conda.

With pip instead of conda, install the alignment-plot extra with `pip install -e ".[msa]"`
(MAFFT, trimAl and FastTree must then be installed separately).

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

**Which records get a region**: coordinates are parsed for every record and given a
`coord_status` in `records.csv` (`ok`, `missing`, `unparseable`, `centroid`, `zero`,
`out_of_range`). Only records with `coord_status` ok are assigned a region and mapped.
Country- and province-centroid coordinates keep their `lat`/`lon` but are not assigned a
region or mapped, because a centroid says nothing about where in the country a specimen
was collected. A record without a region has an empty value (v1 wrote `Unknown`).

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

#### Registering shapefiles (recommended)

Rather than passing `--custom-shp`, `--shp-field`, and `--geo-category` on every run, register each shapefile once under a shorthand name and select it with a single flag. This also works from any working directory, so you can run BOLDGenotyper across different projects without editing paths.

Create a starter config:
```bash
boldgenotyper --init-shapefile-config
```
This writes a commented template to `~/.config/boldgenotyper/shapefiles.yaml` (`%APPDATA%\boldgenotyper\shapefiles.yaml` on Windows). Edit it to register your shapefiles:
```yaml
# Optional base dir for resolving the relative paths below.
# If omitted, relative paths resolve against this file's own directory.
shapefiles_dir: ~/boldgenotyper-data/shapefiles

shapefiles:
  GOaS_v1_20211214:
    path: GOaS_v1_20211214/goas_v01.shp
    field: name
    geo_category: ocean_basin
    type: goas
  Ecoregions2017:
    path: Ecoregions2017/Ecoregions2017.shp
    field: REALM
    geo_category: biogeographic_realm
  BasinATLAS_v10:
    path: BasinATLAS/BasinATLAS_lev03.shp
    field: PFAF_ID
    geo_category: drainage_basin
```
Then select a registered shapefile by name:
```bash
boldgenotyper data/Salmonidae.tsv --shp BasinATLAS_v10 --build-tree
```

**How the config is found** (first match wins):
1. `--shapefile-config <path>`
2. `BOLDGENOTYPER_SHAPEFILE_CONFIG` environment variable
3. `~/.config/boldgenotyper/shapefiles.yaml` (the default)

**Path resolution**: absolute paths and `~`/environment variables are used as given; relative paths resolve against `shapefiles_dir` if set, otherwise against the config file's own directory. Keep the config and the (large) shapefile data wherever you like.

**Notes**:
- Explicit `--shp-field` or `--geo-category` override the registered values, so you can register defaults and tweak per run: `--shp Ecoregions2017 --geo-category terrestrial_realm`.
- `--shp` and `--custom-shp` cannot be combined.
- GOaS is pre-registered so `--shp GOaS_v1_20211214` works uniformly; GOaS is still used automatically when no `--shp` is given.

#### Passing flags directly (alternative)

If you prefer not to keep a config, pass the three flags on each run.

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
- `--geo-category`: Label for the region type in outputs (e.g., `drainage_basin`, `biogeographic_realm`).
  Avoid names of BOLD fields such as `realm`, `biome` or `ecoregion`: the region column would
  replace that field (the log warns when it does)

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

**Set the error-filter threshold explicitly and add the extras**:
```bash
boldgenotyper data/Carcharhinus.tsv \
  --singleton-distance 0.015 \
  --threads 8 \
  --build-tree \
  --extended
```

Always pass `--singleton-distance` (chosen with `boldgenotyper-sweep`, Step 4) so the
value is explicit; it is recorded in `run_provenance.json` either way.

### Step 3: Review Results

Results are organized in the output directory:
```
{organism}_output/
├── README.md                  # every file and column in this folder, defined
├── run_provenance.json        # configuration, tool versions, input checksum, run status
├── {organism}_pipeline.log    # log; ends with a warnings summary
├── summary.csv                # record counts, assignment percentages, species counts
├── records.csv                # one row per input record, with its status and haplotype
├── haplotypes.csv             # one row per haplotype
├── haplotypes.fasta
├── label_conflicts.csv
├── metadata_coverage.csv
├── divergence/                # aligned haplotypes and pairwise p-distances
├── phylogeny/                 # tree files and figure (--build-tree)
├── geography/                 # records per region and the distribution map
├── report.html                # summary report: key findings, figures, methods, every file
└── extended/, exports/, plots/  (with --extended, --export-format, --plot-kit)
```

Start with `README.md` in the run folder, then **open the HTML report** from the run
folder (it links to the figures there):
```bash
open {organism}_output/report.html
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

**Compare a species-level and a family-level run** (labels and haplotypes):
```bash
# Run both analyses with the same --singleton-distance
boldgenotyper data/Sphyrna_lewini.tsv --singleton-distance 0.015 --output Sphyrna_lewini_output/
boldgenotyper data/Sphyrnidae.tsv --singleton-distance 0.015 --output Sphyrnidae_output/
```
Then match records between the two `records.csv` files by `processid`
(see [Comparing Species- and Family-Level Runs](#comparing-species--and-family-level-runs)).

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
| `coord` | Coordinates as BOLD portal exports write them (`[lat, lon]`); `lat`/`lon` columns also work |
| `country/ocean` | Country or ocean of collection (country counts) |
| `province/state`, `region`, `site` | Locality |
| `coord_accuracy`, `coord_source` | Coordinate precision and source |
| `genus` | Used as the record's label when the species is missing |
| `notes` | Depositor notes (searched for keywords such as "juvenile" or "cf.") |
| `collection_date_start` | Collection-date summaries (`--extended`) |
| `sex`, `life_stage`, `reproduction`, `habitat` | Metadata coverage and cross-tabulations (default fields, with `country/ocean` and `province/state`) |
| `realm`, `biome`, `ecoregion`, `geoid` | Summarized only when named in `--metadata-fields`: BOLD derives `realm`, `biome` and `ecoregion` from the coordinates using terrestrial ecoregions, and `geoid` is a numeric location code |
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

The run folder's `README.md` defines every output; this section describes the steps.

### 1. Input and coordinates
- Parses the BOLD TSV once; duplicate process IDs keep their first record
- Parses coordinates for every record and gives each a `coord_status`
  (`ok`, `missing`, `unparseable`, `centroid`, `zero`, `out_of_range`)
- Assigns regions (GOaS ocean basins by default, or a custom shapefile) only to records
  with `coord_status` ok; centroid coordinates keep `lat`/`lon` but get no region

### 2. Quality control
- Orients sequences and checks the open reading frame under the chosen genetic code
  (`--genetic-code`; default 2, vertebrate mitochondrial)
- Removes sequences that are too short (absolute and median-relative minimum length),
  have too many Ns, or lack a valid reading frame (`failed_qc`, with the reason)

### 3. Haplotype discovery
- Aligns QC-passed sequences with MAFFT and extracts the core region (columns covered by
  most sequences); records whose core is too short are excluded (`excluded_core_coverage`)
- Each distinct core sequence is a haplotype (exact sequence variant), keyed by `seq_md5`
  and numbered `H0001`, `H0002`, ... by member count
- Error filter: singleton haplotypes within `--singleton-distance` of another haplotype are
  treated as likely sequencing errors and removed; their records are assigned to the
  nearest retained haplotype when it is unique and within the same distance
  (`assigned_corrected`), otherwise left unassigned (`excluded_error_singleton`)
- Haplotypes that differ only at gaps or Ns are merged; distant singletons, distance
  outliers and reading-frame failures are flagged as suspect

### 4. Record statuses and labels
- Every unique record gets one status: `assigned`, `assigned_corrected`,
  `excluded_error_singleton`, `excluded_core_coverage`, `failed_qc` or `no_sequence`
- A record's label is its species, or its genus when the species is missing
- Each haplotype gets its member labels with counts, a majority label and fraction
  (members without a species label left out), and an assigned taxon: the majority species
  when its fraction is above 0.70, otherwise the genus
- Evidence for each haplotype's taxon (nearest species-level haplotype, BIN, tree sister,
  label sources) and, for haplotypes left at genus, a suggested taxon that is never applied
- `--taxon-curation FILE` applies taxa you have checked; see
  [Reviewing and Curating Haplotype Taxa](#reviewing-and-curating-haplotype-taxa)
- Label conflicts are listed with their evidence in `label_conflicts.csv`

### 5. Phylogeny (`--build-tree`)
- Every non-suspect haplotype with at least `--tree-min-members` members (default 1) is a
  tip; excluded haplotypes are listed in `phylogeny/excluded_tips.txt`
- MAFFT alignment (trimmed with trimAl when installed) and a FastTree tree (GTR+Gamma)
- Rooted on the outgroup when one is given (`--phylo-outgroup-fasta`, `--phylo-outgroup-label`
  or `--phylo-outgroup-taxon`), otherwise at the midpoint of the longest path between two
  tips; support values stay on their branches when the tree is rooted. Set
  `phylogenetic.midpoint_root: false` in a config file to keep FastTree's unrooted tree
- `tree.nwk` (haplotype IDs) and `tree_labeled.nwk` (haplotype labels) for tree editors
  such as [TreeViewer](https://treeviewer.org/), FigTree or iTOL

### 6. Divergence and geography
- Uncorrected p-distances between aligned haplotypes (`divergence/`)
- Records per haplotype and region, and a distribution map of records with
  geographic-quality coordinates (`geography/`)
- Metadata field coverage (`metadata_coverage.csv`), after normalizing values

### 7. Report and run folder
- `summary.csv`, `records.csv`, `haplotypes.csv`, the run `README.md`, and `report.html`
  (interactive; links to the figures in the run folder)
- `run_provenance.json` with the configuration, tool versions and run status; the log
  ends with a warnings summary

### Optional outputs
- `--extended` writes `extended/`: metadata cross-tabulations and collection-date
  summaries (these show sampling, not trends), a barcoding gap computed from record
  labels, species tables and diversity, species-level divergence, alignment plots, and
  region and metadata figures (an overview plus one page per species)
- `--export-format` writes `exports/` (Arlequin, PopART/NEXUS, DnaSP, generic CSV)
- `--plot-kit` writes `plots/`: plot data, a YAML configuration and scripts that
  regenerate the map, bar charts and tree (see [Customizing Plots](#customizing-plots-for-publication))

### Advanced Workflows

Beyond the main pipeline:

#### Parameter Sweep (`boldgenotyper-sweep`)

Runs the pipeline at several `--singleton-distance` thresholds and reports:
- Number of haplotypes and singletons per threshold, with the elbow of the curve
- Assignment coverage per threshold
- For each pair of adjacent thresholds, the share of records whose haplotype changes
  (compared by sequence), from `group_membership_tracking.csv`

**When to use**: Before running your final analysis: run the sweep, read
`recommendations.txt`, choose a threshold, and pass it to `boldgenotyper --singleton-distance`.

**Output**: `recommendations.txt`, `elbow_plot.pdf`, `threshold_stability.pdf`, `sweep_summary.csv`, `group_membership_tracking.csv`

See `parameter_sweep/README.md` for interpretation guide.

#### Comparing Species- and Family-Level Runs

A family-level run already reports label conflicts within its haplotypes
(`label_conflicts.csv`). To see which records of one species fall in family-level
haplotypes dominated by another species' label, run the species and its family with the
same settings and match the two `records.csv` files by `processid` (haplotype IDs and
core regions differ between runs on different datasets, so match records, not IDs). No
separate command is needed: [`docs/guides/COMPARATIVE_ANALYSIS_GUIDE.md`](docs/guides/COMPARATIVE_ANALYSIS_GUIDE.md)
gives the pandas code and how to read the evidence, and Tutorial 3 runs it.

**When to use**: Checking labels in multi-species datasets

---

## Usage Guide

### Basic Command Structure

```bash
boldgenotyper <input_tsv> [options]
```

### Complete Example

```bash
# Run boldgenotyper-sweep first to choose --singleton-distance for your dataset
boldgenotyper data/Sphyrna_lewini.tsv \
  --organism Sphyrna_lewini \
  --output results/Sphyrna_analysis \
  --singleton-distance 0.015 \
  --download-date 2026-09-30 \
  --threads 8 \
  --build-tree \
  --extended \
  --log-level INFO
```

### Common Use Cases

**1. Standard Analysis (Marine Organism)**:
```bash
boldgenotyper data/Sphyrna_lewini.tsv --singleton-distance 0.015 --build-tree
```

**2. Invertebrate Organism**:
```bash
# Use --genetic-code 5 (invertebrate mitochondrial) for molluscs, crustaceans,
# insects, worms, echinoderms, etc.
boldgenotyper data/Crassostrea_virginica.tsv --genetic-code 5 --build-tree
```

**3. Non-Marine Invertebrate (Skip Geographic Analysis)**:
```bash
boldgenotyper data/Anopheles_gambiae.tsv --genetic-code 5 --no-geo --build-tree
```

**4. Keep Error-Filtered Singletons Unassigned**:
```bash
# Records of singleton haplotypes removed as likely errors stay unassigned
# (status excluded_error_singleton) instead of joining their nearest haplotype
boldgenotyper data/Population_samples.tsv --no-assign-corrected --build-tree
```

**5. Exports and Plot Kit**:
```bash
boldgenotyper data/my_species.tsv --export-format arlequin popart --plot-kit
```

**6. Skip HTML Report (Faster, CI/CD)**:
```bash
boldgenotyper data/my_species.tsv --no-report --threads 16
```

**7. Custom Output Location**:
```bash
boldgenotyper data/my_species.tsv --output /mnt/storage/results/my_analysis_2026
```

Use a new `--output` folder for each run. A run that writes into a folder that already
holds a run continues, but warns (in the log, the end-of-run summary and
`run_provenance.json`), because files from the earlier run could be mistaken for new ones.

### Querying New Sequences

After a run, `boldgenotyper-query` places new COI sequences on that run's haplotypes
without re-running the pipeline: assigning new specimens, checking field samples, or
comparing studies. The [query guide](docs/guides/QUERY_GUIDE.md) covers the outputs,
reporting a query, figures, curated taxa and repeat queries in detail.

**Query input**: one consensus sequence per specimen in a single (multi-)FASTA, not raw
forward/reverse Sanger reads. Gaps (`-`, `.`) and lowercase are cleaned automatically.

**Query against a run folder (recommended)**:
```bash
# --haplotypes takes the run folder: it finds haplotypes.fasta and uses the run's
# record labels for species support
boldgenotyper-query \
  --query new_samples.fasta \
  --haplotypes results/Sphyrna_analysis \
  --output query_results/
```

**Also search every QC-passed record** (including records not assigned to a haplotype):
```bash
boldgenotyper-query --query new_samples.fasta --haplotypes results/Sphyrna_analysis \
  --search-records --output query_results/
```

**With taxa curated for this query** (a run folder's own curated taxa are used without it):
```bash
boldgenotyper-query --query new_samples.fasta --haplotypes results/Sphyrna_analysis \
  --taxon-curation curation.csv --output query_results/
```

**Against a FASTA only** (haplotype matching without species support):
```bash
boldgenotyper-query --query samples.fasta --haplotypes haplotypes.fasta --top-n 20 \
  --min-length 150 --max-length 1500 --output results/
```

A v1 run folder (no `run_provenance.json`) also works as `--haplotypes`, for haplotype
matching only.

**How queries are matched**: each query is aligned over its whole overlap with each
haplotype (global alignment with free end gaps), so differences at the query's ends are
counted. `differences` = mismatches + gap columns. A match is eligible when it covers at
least `--min-query-coverage` (0.95) of the shorter of query and haplotype and reaches
`--assign-identity` (0.97). The equivalence set is every eligible haplotype within
`--max-mismatch-delta` (0) differences of the best one.

**With sample metadata** (adds the queries' locations, regions and dates to the figures):
```bash
boldgenotyper-query --query new_samples.fasta --haplotypes results/Sphyrna_analysis \
  --query-metadata samples.csv --output query_results/
```
`samples.csv` (CSV, or TSV for `.tsv`/`.txt`) has a `query_id` column matching the FASTA
IDs and any of `lat`/`lon` or `coord` (`[lat, lon]`), `collection_date`, `country/ocean` and
`region`. A query's region is its `region` value; otherwise its coordinates are placed with
the reference run's own shapefile and field, when that shapefile is still at the path the
reference run recorded.

**Query output files**:
- `query_calls.csv`: one row per query; **the result to report**
- `query_candidates.csv`: the top `--top-n` haplotypes per query with the numbers behind the
  call, and each haplotype's taxon, `assignment_source` and `n_members`
- `query_alignments.txt`: each query's call and its best alignment
- `query_report.html`: the calls, one evidence card per query, the figures and the methods
- `query_figures/`: the queries placed in the reference run (v2 run folders): the reference
  tree pruned to each called genus with the queries beside their haplotypes (when the
  reference built a tree), a map per called taxon, reference records per region with the
  queries marked (`query_regions.csv`), and sampling years with the query dates; each as a
  multi-page PDF plus one SVG per page
- `query_provenance.json`: the reference run, parameters, curation and metadata files (with
  md5), figures, warnings and the files written
- `README.md`: which file to report, a methods paragraph with the run's values, and every
  file and column
- `query_record_matches.csv` (`--search-records`): records that contain the query
  unchanged, with their label, haplotype and status
- `query_results.json` (`--json`): everything as one JSON document

The CSVs have no comment lines, so spreadsheets open them directly.

**Reading `query_calls.csv`**:
- `call_level`: `haplotype` (one haplotype in the equivalence set), `haplotype_set`
  (several, listed in `equivalence_set`; the flag `set_not_identical` says they differ over
  the query's region) or `none` (no eligible match)
- `pooled_member_labels`: the labels of every record of the set's haplotypes, with counts
- `species_call` and `species_support`: the leading species and its share of the counted
  labels (genus-only labels of the leading genus are compatible); `species_call` is empty
  unless the support reaches `--min-species-support` (0.70)
- `call_basis`: `pooled labels` (the call above), or `haplotype taxa` when a haplotype in
  the set has a curated taxon: the call is then the set's shared haplotype taxon, and no
  call (flag `curated_taxa_disagree`) when the set's taxa differ or one is not at species
  level
- `flags`: e.g. `label_conflict_in_set`, `indel_in_query`, `no_reference_labels`
- `reference_run`: the md5 of the reference run's `run_provenance.json`

**Repeat queries**: keep one reference run per reference dataset and query the full set of
samples each time; a query run takes seconds. Every output records the reference run's md5,
so results from different dates can be compared.

**Important Notes**:
- Query sequences must be COI regions comparable to the reference (not whole genomes)
- Multi-FASTA files are supported for batch processing
- Every match is reported and ranked; eligibility decides the call

### Input File Requirements

**Filename Convention**: Your input TSV should follow this pattern for automatic organism detection:
```
Genus_species.tsv
```

**Override organism name** with `--organism` flag if needed:
```bash
boldgenotyper data/bold_download_2026.tsv --organism Sphyrna_lewini
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
| `--output` | Path | `./{organism}_output` | Output directory (use a new folder per run) |
| `--singleton-distance` | Float | 0.005 | Error filter: singleton haplotypes within this p-distance of another haplotype are removed as likely errors; choose with `boldgenotyper-sweep` and pass explicitly |
| `--no-assign-corrected` | Flag | False | Leave records of error-filtered singletons unassigned instead of assigning them to their nearest haplotype |
| `--genetic-code` | Integer | 2 | NCBI genetic code for ORF validation (2 = vertebrate mitochondrial; 5 = invertebrate mitochondrial) |
| `--threads` | Integer | 4 | Number of parallel CPU threads |
| `--build-tree` | Flag | False | Build a phylogenetic tree with FastTree |
| `--tree-min-members` | Integer | 1 | Minimum members for a haplotype to be a tree tip (suspect haplotypes are always excluded) |
| `--phylo-outgroup-fasta`, `--phylo-outgroup-label`, `--phylo-outgroup-taxon` | | None | Root the tree on an outgroup |
| `--extended` | Flag | False | Write `extended/` (cross-tabs, collection-date summaries, label-based gap, species tables, MSA plots, extra maps and bars) |
| `--export-format` | List | None | Write `exports/` in these formats: `arlequin`, `popart`, `dnasp`, `generic`, `all` |
| `--plot-kit` | Flag | False | Write `plots/` (plot data, config and regeneration scripts) |
| `--no-report` | Flag | False | Skip HTML report generation |
| `--no-geo` | Flag | False | Skip geographic analysis |
| `--custom-shp` | Path | None | Custom shapefile for geographic region assignment |
| `--shp` | String | None | Name of a registered shapefile (see `--init-shapefile-config`) |
| `--shp-field` | String | Auto | Shapefile attribute field containing region labels |
| `--geo-category` | String | `ocean_basin` | Name of the region column in outputs; avoid BOLD field names (`realm`, `biome`, `ecoregion`), which it would replace |
| `--metadata-fields` | List | sex, life_stage, reproduction, habitat, country/ocean, province/state | Metadata fields for coverage and cross-tabs |
| `--no-normalize-sex` | Flag | False | Keep sex values as given (by default M/F → Male/Female, U → missing) |
| `--no-metadata-summary` | Flag | False | With `--extended`, skip the metadata cross-tabs (coverage is always written) |
| `--no-temporal-analysis` | Flag | False | With `--extended`, skip the collection-date summaries |
| `--download-date` | Date | None | Date the BOLD TSV was downloaded (YYYY-MM-DD); recorded in `run_provenance.json` |
| `--taxon-curation` | Path | None | CSV/TSV of checked haplotype taxa applied after the label rule (see [Reviewing and Curating Haplotype Taxa](#reviewing-and-curating-haplotype-taxa)) |
| `--keep-intermediates` | Flag | False | Keep `intermediate/` |
| `--log-level` | String | INFO | Logging verbosity (DEBUG, INFO, WARNING, ERROR) |
| `--version` | Flag | - | Show version and exit |

### Parameter Details

#### `--singleton-distance` (default: 0.005)

**What it controls**: singleton haplotypes (one member) within this p-distance of another
haplotype are treated as likely sequencing or PCR errors and removed. Their records are
assigned to the nearest retained haplotype when it is unique and within the same distance
(`assigned_corrected`); otherwise they stay unassigned (`excluded_error_singleton`).

**Background**: Because ESV discovery treats every unique sequence as a distinct haplotype,
error-derived sequences that differ by 1–2 bp from a true haplotype appear as spurious
singletons. The default (0.5%) removes most such errors while retaining biologically real
rare variants.

**Run `boldgenotyper-sweep` first, choose a value, and pass it explicitly**:
```bash
# Step 1: test thresholds for your dataset
boldgenotyper-sweep data/MySpecies.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.03,0.05 \
  --output sweep_results/

# Step 2: read the haplotype-count table (elbow marked) and the stability table
cat sweep_results/recommendations.txt

# Step 3: apply the chosen value
boldgenotyper data/MySpecies.tsv --singleton-distance 0.015 --build-tree
```

**When to adjust**:
- **Lower (0.001–0.005)**: Retain more rare variants; risk including sequencing errors
- **Higher (0.02–0.10)**: Remove more singletons; risk losing biologically real rare haplotypes
- Prefer values where few records change haplotype between neighboring thresholds
  (the stability table in `recommendations.txt`)

#### `--genetic-code` (default: 2)

**What it controls**: The NCBI genetic code table used for ORF validation during quality control. Sequences are translated in all six reading frames and the result is used to confirm they are valid COI amplicons and to select the correct orientation.

**When to change**: The default (code 2, vertebrate mitochondrial) is correct for fish, sharks, mammals, birds, and reptiles. **Most invertebrates require code 5** (invertebrate mitochondrial). With the wrong code, codons such as AGA/AGG translate as stops, so many (sometimes all) sequences can fail ORF validation, and those that pass were checked against the wrong code.

| Code | Name | Example taxa |
|------|------|-------------|
| 2 | Vertebrate mitochondrial (default) | Fish, sharks, mammals, birds, reptiles |
| 5 | Invertebrate mitochondrial | Molluscs, crustaceans, insects, worms, echinoderms |
| 4 | Mold/Protozoan/Coelenterate mitochondrial | Cnidaria (jellyfish, corals), some fungi |

Full table: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

```bash
# Invertebrates (molluscs, crustaceans, insects, etc.)
boldgenotyper data/Crassostrea_virginica.tsv --genetic-code 5

# Cnidarians
boldgenotyper data/Acropora_millepora.tsv --genetic-code 4

# Vertebrates (default — flag not required)
boldgenotyper data/Sphyrna_lewini.tsv
```

#### `--tree-min-members` (default: 1)

Every non-suspect haplotype with at least this many members (and the minimum length)
is a tree tip. Haplotypes left out are listed with the reason in `phylogeny/excluded_tips.txt`,
and `haplotypes.csv` has `in_tree` for each haplotype.

#### `--threads`

**What it controls**: Number of CPU cores used by MAFFT and FastTree.

**Recommendation**: Set to number of available cores (check with `nproc` or `sysctl -n hw.ncpu`)

### Advanced Command Parameters

#### `boldgenotyper-query` Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--query` | Path | Required | Query (multi-)FASTA |
| `--haplotypes` | Path | Required | A run folder (v2 or v1) or a haplotype FASTA |
| `--output` | Path | `haplotype_query_results/` | Output directory |
| `--top-n` | Integer | 10 | Matches reported per query |
| `--assign-identity` | Float | 0.97 | Minimum identity for an eligible match |
| `--min-query-coverage` | Float | 0.95 | Minimum coverage of the shorter of query and haplotype |
| `--max-mismatch-delta` | Integer | 0 | Differences from the best match allowed in the equivalence set |
| `--min-species-support` | Float | 0.70 | Minimum share of counted labels for a species call |
| `--search-records` | Flag | False | Also search QC-passed records (`records_qc_passed.fasta`) |
| `--taxon-curation` | Path | None | Curation file applied to the reference run's haplotype taxa for this query (needs a v2 run folder) |
| `--query-metadata` | Path | None | CSV/TSV of query locations, dates and regions for the figures |
| `--json` | Flag | False | Also write `query_results.json` |
| `--min-length`, `--max-length` | Integer | 100, 2000 | Query length limits |

#### `boldgenotyper-sweep` Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tsv` | Path | Required | Input BOLD TSV file |
| `--thresholds` | String | 0.01,0.015,0.02,0.03,0.05 | Comma-separated threshold values to test (at least 4 to locate an elbow) |
| `--output` | Path | parameter_sweep/ | Output directory |
| `--threads` | Integer | 4 | Number of parallel threads (sequential by default) |
| `--keep-intermediates` | Flag | False | Keep full output for each threshold tested |
| `--genetic-code` | Integer | 2 | NCBI genetic code for reading-frame QC; use the code of the final run (5 for most invertebrates) |
| `--log-level` | String | INFO | Logging verbosity |

**Example**:
```bash
boldgenotyper-sweep data/samples.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.025,0.03 \
  --threads 8 \
  --output sweep_results/
```

---

## Output Files

Every run folder has a generated `README.md` that lists each file and defines every
column. This section is an overview. Upgrading from 1.0.0: see the
[migration guide](docs/guides/MIGRATION_v1_to_v2.md).

### Directory Structure

```
{organism}_output/
├── README.md                     # every file and column, defined
├── run_provenance.json           # configuration, tool versions, seeds, input md5, run status
├── {organism}_pipeline.log       # log; ends with a warnings summary
├── summary.csv                   # one row: record counts, assignment percentages, species counts
├── records.csv                   # one row per unique input record
├── records_qc_passed.fasta       # QC-passed sequences (boldgenotyper-query --search-records)
├── haplotypes.fasta              # >H0001 seq_md5=... n_members=...
├── haplotypes.csv                # one row per haplotype
├── label_conflicts.csv           # records whose label conflicts with their haplotype's majority
├── metadata_coverage.csv         # how many records have each metadata field
├── divergence/
│   ├── haplotypes_aligned.fasta
│   └── pairwise_p_distance.csv
├── phylogeny/                    # with --build-tree
│   ├── tree.nwk                  # tips named by haplotype ID
│   ├── tree_labeled.nwk          # tips named "Genus species H0001"
│   ├── excluded_tips.txt         # haplotypes left out, with the reason
│   ├── tree.pdf, tree.svg        # every tip
│   ├── tree_overview.pdf, .svg   # taxa that form a clade collapsed (above 40 haplotypes)
│   ├── tree_by_genus.pdf         # one page per genus; by_genus/ has one SVG per page
│   └── manifest.json             # the figures in this folder and their pages
├── geography/                    # when records have geographic-quality coordinates
│   ├── haplotype_by_region.csv
│   ├── distribution_map.pdf, .svg            # overview
│   ├── distribution_map_by_species.pdf       # one page per taxon; by_species/ has the SVGs
│   └── manifest.json
├── report.html                   # open from this folder (figures are linked)
├── extended/                     # with --extended
│   ├── metadata/                 # cross-tabulations, first-sampled and temporal summaries; figures/
│   ├── divergence/               # label-based gap, species-level divergence, heatmap, msa/ plots
│   ├── species/                  # species summary and diversity, per-species record files
│   ├── geography/                # records per region (relative, total), by species
│   └── labels/                   # depositor keyword hits, label composition figures
├── exports/                      # with --export-format (arlequin/, popart/, dnasp/, generic/)
├── plots/                        # with --plot-kit
└── intermediate/                 # with --keep-intermediates
```

### Figures

Figures scale to large runs. Each has an **overview** and, for most, **one page per
species** (or genus-level group): a multi-page PDF for publication, and one SVG per page
in `by_species/<figure>/` for the report and web use. Every figure folder has a
`manifest.json` listing each figure's overview and pages.

- Overviews show each haplotype for up to 40 haplotypes and each taxon above
  (`visualization.max_overview_haplotypes`).
- A panel colors at most 12 items (`max_colored_items`); the rest are gray as
  "Other (n)". Taxon overviews color 15 taxa (`max_overview_species`). Within a taxon a
  haplotype keeps its color on every page.
- The tree: `tree.pdf/.svg` keeps every tip; `tree_overview` collapses each taxon whose
  haplotypes form a clade and lists those that do not; `tree_by_genus.pdf` has the clade
  at each genus's most recent common ancestor.
- Captions give the denominators (for example, how many records have geographic-quality
  coordinates or a region).
- "By region" metadata views use the run's region column (`--geo-category`: ocean basin,
  or the category of a custom shapefile), or `country/ocean` when the run has fewer than
  two regions; the caption says which.
- Default runs draw the tree and the distribution map; `--extended` adds region
  abundance and the metadata views. No PNG files are written except by `--plot-kit`.

### Key Output Files

#### `records.csv`

One row per unique input record, in input order, including records that failed QC or
have no sequence. Columns include `status` (see below), `status_note`, `haplotype_id`,
`haplotype_label`, `assigned_sp`, `label`, `label_conflict`, `nearest_haplotype` and
`distance_to_nearest` (for error-filtered singletons), `qc_pass`, `qc_fail_reason`,
`depositor_keyword_hits`, `lat`, `lon`, `coord_status`, the region column
(e.g. `ocean_basin`), and key BOLD fields. Join any other input column on `processid`.

| `status` | Meaning |
|---|---|
| `assigned` | core sequence identical to a retained haplotype |
| `assigned_corrected` | error-filtered singleton assigned to its nearest retained haplotype, within `--singleton-distance` |
| `excluded_error_singleton` | error-filtered singleton left unassigned (tied, too far, or `--no-assign-corrected`) |
| `excluded_core_coverage` | core region too short |
| `failed_qc` | failed length, N or reading-frame checks (see `qc_fail_reason`) |
| `no_sequence` | no sequence in the input |

#### `haplotypes.csv`

One row per haplotype: `haplotype_id`, `haplotype_label`, `seq_md5` (stable across runs
on the same data), `length`, `n_members` (`n_exact` + `n_corrected`), the assigned taxon,
member `labels` with counts, `majority_label` and `majority_fraction`, `n_label_conflicts`,
the nearest other haplotype and its p-distance, suspect flags, `in_tree`, the members
with coordinates, a region and a country, `assignment_source` (`labels` or `curated`), and
the taxon evidence columns (see
[Reviewing and Curating Haplotype Taxa](#reviewing-and-curating-haplotype-taxa)).

#### `summary.csv`

Record counts at each step with explicit denominators (percentages are of QC-passed
records), and species counts under named definitions (`n_species_labels_input`,
`n_species_labels_qc_passed`, `n_species_haplotype_assigned`, `n_species_confident`, ...).

#### `label_conflicts.csv`

Records whose species-level label differs from their haplotype's majority label (or whose
haplotype has tied top labels), with the evidence: members, labeled members, the majority
label and fraction, and the nearest other haplotype carrying the record's label. For a
haplotype with a curated taxon, labels are compared with that taxon (`haplotype_taxon`,
`assignment_source`). A conflict is evidence to look at, not a verdict on which label is
right.

#### `run_provenance.json`

The full resolved configuration (including `--singleton-distance`), the values the run
applied (core length, cutoffs, tree tips), tool versions, seeds, the input path, md5 and
row count, `--download-date`, warnings and the run status. `boldgenotyper-query` reports
its md5 as `reference_run`.

#### `report.html`

Built from the final run folder after every file is in place, in nine sections:
Overview (counts and **key findings**: a fixed list of counts, each shown when it applies,
none interpreting biology), Records and haplotypes (statuses, error filter), Taxonomy and
labels (**needs attention**: haplotypes not assigned at species level with their evidence
and any suggested taxon, label conflicts per haplotype and label, curated haplotypes,
curation rows not applied), Phylogeny (tree overview with a genus menu), Geography (map
with a species menu, region bars, coordinate quality), Sampling and metadata (coverage,
collection years, timeline, field views with field and view menus), Divergence
(`--extended`), Methods (the values the run used and a methods statement to edit), and
Output files (every file, described as in the run README). Tables show ten rows and say how
many more are in which file. Figures are linked, not embedded: open the report from the run
folder.

#### Trees

`phylogeny/tree_labeled.nwk` opens in any tree editor for re-rooting and annotation:
**[TreeViewer](https://treeviewer.org/)** (Bianchini & Sánchez-Baracaldo, 2024),
**FigTree**, **iTOL**, or **R (ape/ggtree)**.

#### `{organism}_pipeline.log`

Written fresh for each run: input and QC counts, core region and error-filter results,
record statuses, tree tips, and every step that failed or was skipped, summarized at the
end with the run status.

---

## Customizing Plots for Publication

With `--plot-kit`, BOLDGenotyper writes the plot data and Python scripts that regenerate
the distribution map, region bar charts and tree with your own styling, without re-running
the analysis.

**Additional resources:**
- `{organism}_output/plots/README.md` (written with the kit)
- `notebooks/06_plot_customization_tutorial.ipynb`

### Plot Regeneration Directory

```
{organism}_output/plots/
├── plot_config.yaml          # settings (colors, sizes, formats, projection, filters)
├── README.md
├── data/                     # plot data (CSV)
│   ├── distribution_map.csv
│   ├── distribution_bar_relative.csv
│   └── distribution_bar_absolute.csv
└── scripts/
    ├── regenerate_all.py     # map, bars and tree
    ├── regenerate_map.py
    └── regenerate_bars.py
```

The scripts read the haplotype sequences and tree from the run folder
(`../haplotypes.fasta`, `../phylogeny/tree.nwk`); the kit has no copies of them. The data
files keep v1 column names (`consensus_group` is the haplotype ID, `consensus_group_sp` its
label) so existing scripts work.

### Quick Start: Customize Your Plots

**1. Edit the configuration file:**

```bash
cd {organism}_output/plots
nano plot_config.yaml  # or use any text editor
```

**2. Modify settings:**

```yaml
general:
  output_format: ['pdf', 'png']
  dpi: 600  # High resolution for publication
  width_inches: 12
  height_inches: 8

colors:
  "Sphyrna lewini H0001": "#E41A1C"   # keys are haplotype labels
  "Sphyrna lewini H0002": "#377EB8"

filters:
  include_genotypes: []                 # Leave empty to include all
  exclude_genotypes: ["Sphyrna lewini H0010"]

map:
  projection: "mollweide"             # robinson, mollweide, mercator
  center_longitude: -180              # Center on Pacific Ocean
  point_alpha: 0.8
  point_size_range: [3, 12]
  legend_position: "right"

bars:
  orientation: "horizontal"           # vertical or horizontal
  bar_width: 0.7
```

**3. Regenerate plots:**

```bash
python scripts/regenerate_all.py      # or regenerate_map.py, regenerate_bars.py
```

**4. Find updated plots** in `{organism}_output/plots/visualization/`, with a `_custom` suffix.

### Customization Options

**Colors**: hex colors per haplotype label
**Figure size**: width, height, DPI for different journals
**Map projection**: Robinson, Mollweide, Mercator, PlateCarree
**Filters**: include or exclude haplotypes
**Bar orientation**: vertical or horizontal stacking

### Requirements

All required packages are in the boldgenotyper environment:
- matplotlib, seaborn (plotting)
- geopandas, cartopy (maps)
- pandas (data handling)
- pyyaml (configuration)

### Advanced Customization and Troubleshooting

See the **interactive tutorial**: `notebooks/06_plot_customization_tutorial.ipynb`
(YAML syntax errors, example configurations, all settings, ways to regenerate plots).

---

### Advanced Command Outputs

#### Parameter Sweep Output (`parameter_sweep/`)

```
parameter_sweep/
├── sweep_summary.csv              # metrics for each tested threshold
├── threshold_stability.pdf        # metrics against threshold
├── elbow_plot.pdf                 # haplotype count against threshold, elbow marked
├── group_membership_tracking.csv  # each record's haplotype (ID and seq_md5) per threshold
├── recommendations.txt            # haplotype-count table and stability table
├── README.md                      # interpretation guide
└── runs/                          # record and haplotype tables per threshold (full run folders with --keep-intermediates)
```

**Key files**:
- `sweep_summary.csv`: haplotypes, singletons, suspect haplotypes, assignment percentage,
  label-conflict share and mean majority fraction per threshold
- `group_membership_tracking.csv`: `t_<threshold>` (haplotype ID) and `md5_<threshold>`
  (its `seq_md5`) per record, plus how often its co-members change
- `recommendations.txt`: the elbow of the haplotype-count curve and, for each pair of
  adjacent thresholds, the share of records whose haplotype changes

#### Population Genetics Exports (`exports/`)

With `--export-format`:

```
exports/
├── README.md                      # format descriptions and usage
├── arlequin/
│   ├── {organism}.arp             # Arlequin project file
│   └── populations.csv            # populations (regions)
├── popart/
│   ├── {organism}.nexus           # NEXUS alignment with traits
│   └── popart_traits.csv
├── dnasp/
│   └── {organism}.fas             # FASTA with population labels
└── generic/
    ├── haplotype_membership.csv   # record-to-haplotype table
    └── haplotype_summary.csv      # one row per haplotype
```

See [`docs/guides/POPGEN_EXPORT_GUIDE.md`](docs/guides/POPGEN_EXPORT_GUIDE.md) for format details.

---

## Biological Context for Threshold Selection

### Understanding Haplotype Discovery (ESV Approach)

BOLDGenotyper uses the **Exact Sequence Variant (ESV)** approach (Porter & Hajibabaei 2020) for haplotype discovery. Each distinct sequence of the core alignment region is a haplotype; no clustering threshold is involved, and records are mapped to haplotypes directly (there is no identity-based re-assignment).

**Error Filtering and `--singleton-distance`**: Singleton ESVs within `--singleton-distance` (default: 0.5%) of another haplotype are removed as probable sequencing or PCR errors, and their records are assigned to the nearest retained haplotype when it is unique and within the same distance. This threshold is the parameter that most affects the number of haplotypes. **Run `boldgenotyper-sweep` first**: it reports the haplotype count at each threshold with the elbow of the curve marked, and for each pair of neighboring thresholds the share of records whose haplotype changes. Choose a value where the count has leveled off and few records change, and pass it to `boldgenotyper --singleton-distance`. Suitable values vary with taxon and sequence quality (empirically: marine vertebrates ≈ 0.015; some insect families ≈ 0.085).

### Labels and the Haplotype Taxon

BOLD species labels come from depositors and are sometimes wrong, missing, or recorded
only to genus. BOLDGenotyper reports labels rather than judging them:

- A record's label is its species, or its genus when the species is missing.
- Each haplotype lists all member labels with counts. The majority label's share
  (`majority_fraction`) leaves out members without a species label, so a haplotype with
  ten "*Genus species*" records and five genus-only records has a fraction of 1.0.
- The haplotype's taxon (`assigned_sp`) is the majority species when its fraction is above
  0.70; otherwise the genus.
- Members whose species-level label differs from the majority are label conflicts, listed
  with the nearest haplotype that carries their label, so you can judge them.

Divergences are reported as uncorrected p-distances without species-delimitation
thresholds; `--extended` adds a barcoding gap computed from record labels.

### Reviewing and Curating Haplotype Taxa

When a haplotype's labels split, the label rule leaves it at genus. `haplotypes.csv`
reports evidence that does not depend on the haplotype's own majority, for every haplotype:

- **Distance**: `nearest_species_haplotype`, `nearest_species_taxon` and
  `nearest_species_distance` (the nearest other haplotype assigned at species level).
  `taxon_max_within_distance` is the largest distance between the other haplotypes of that
  taxon, and `within_taxon_spread` says whether this haplotype is as close to the taxon as
  its haplotypes are to one another.
- **BIN**: `member_bins` (the members' BOLD BINs with counts) and `bin_taxon` (for the
  most common BIN, the most common species-level label among the run's QC-passed records
  in that BIN, with its fraction). Closely related species can share a BIN; a fraction
  below 0.70 is not counted as evidence.
- **Tree**: `tree_sister_taxon` and `tree_sister_support`. Walking from the haplotype's tip
  toward the root, the first clade that holds another species-level haplotype; its taxon
  when those haplotypes share one. A clade with support below 0.70 is not counted as
  evidence. The tree is midpoint-rooted unless you set an outgroup, so the sister clade of a
  haplotype next to the root depends on where the root falls.
- **Label sources**: `label_sources` counts the records, BOLD projects (processid prefix;
  GenBank accession prefix for records mined from GenBank) and institutions behind each
  label. Many records from one project may be one identification repeated. The counts are
  approximate: one study can span several projects.

For a haplotype neither at species level nor curated, `suggested_sp` is a taxon that at
least two of the distance, BIN and tree lines name and none contradicts;
`suggestion_basis` names the lines that agree, or why there is no suggestion. Suggestions
are never applied. A label split can have several causes (misidentified specimens, shared
haplotypes between species, a lineage not yet described), and weighing them is your call.

To apply a taxon you have checked, list it in a curation file (CSV, or TSV for
`.tsv`/`.txt`) and pass it with `--taxon-curation`:

```csv
seq_md5,taxon,reason,curator,date
0f3c...,Genus species,"within Genus species spread; sister in tree (0.88); conflicting labels from one project",AB,2026-10-04
```

- `seq_md5` or `sequence` names the haplotype (one is required); `taxon` and `reason` are
  required; `curator` and `date` are optional.
- A row applies to the haplotype with that `seq_md5`, or whose core sequence equals,
  contains or is contained in `sequence`. `seq_md5` is the md5 of the core sequence, so it
  matches only runs with the same core region; `sequence` carries a call to another
  dataset. Rows matching no haplotype or several, and rows that disagree for one
  haplotype, are not applied; the log and `run_provenance.json` list each with the reason.
- A curated haplotype gets `assigned_sp` = `taxon`, `assignment_source` = `curated` and
  `assignment_notes` starting with `curated: <reason>`; its label conflicts are judged
  against the curated taxon. `run_provenance.json` records the file's path and md5.
- `boldgenotyper-query` uses a run folder's curated taxa, and `--taxon-curation` applies a
  file at query time; see `call_basis` under [Querying New Sequences](#querying-new-sequences).

---

## Assumptions and Limitations

Understanding these assumptions is critical for interpreting results and planning your analysis.

### Assumptions

#### 1. Sequence Quality
- **Assumption**: Input sequences are of reasonable quality
- **Pipeline behavior**: Removes sequences shorter than 200 bp or shorter than 70% of the median length, with more than 5% Ns, or without a valid open reading frame (status `failed_qc`, with the reason); thresholds are in `run_provenance.json`
- **Implication**: BOLD quality varies by contributor; check `qc_fail_reason` in `records.csv`
- **Recommendation**: Start with public records or sequences with BIN assignments

#### 2. Gap and Ambiguity Handling
- **Assumption**: Gaps (`-`) and ambiguous bases (`N`) are ignored in distance calculations
- **Rationale**: Gaps may result from alignment artifacts; Ns lack information
- **Calculation**: `distance = 1 - (matches / informative_sites)`
  where informative sites = positions with A, C, G, or T in both sequences
- **Implication**: Fragmented sequences with many gaps aren't penalized artificially

#### 3. Haplotypes Are Exact Core Sequences
- **Assumption**: The core alignment region (columns covered by most sequences) is comparable across records
- **Process**: Records whose core is too short are excluded (`excluded_core_coverage`); each distinct core sequence is a haplotype; haplotypes that differ only at gaps or Ns are merged
- **Implication**: Variation outside the core region is not used; the core length is in `run_provenance.json`

#### 4. Geographic Coordinate Quality
- **Assumption**: Country-level coordinates do not locate a specimen
- **Pipeline behavior**: Every record gets a `coord_status`: `ok`, `missing`, `unparseable`, `centroid`, `zero` or `out_of_range`
- **Key Distinction**:
  - All records with sequences are used for haplotypes, whatever their coordinates
  - Only records with `coord_status` ok are assigned a region and mapped
  - Country- and province-centroid coordinates keep `lat`/`lon` but are not assigned a region or mapped; a record without a region has an empty region value
- **Rationale**: A centroid for Mexico could fall in the Pacific or the Atlantic, but the sequence is still valid for haplotypes

#### 5. Ocean Basin Boundaries
- **Assumption**: GOaS polygon definitions are accurate and appropriate
- **Implementation**:
  - Shapefile defines basin polygons
  - Coordinate-to-basin assignment via spatial join
  - Marginal seas (Mediterranean, Caribbean) included in parent basins
- **Limitation**: Disputed waters, transitional zones not explicitly handled; `region_confidence` and `distance_to_region_boundary_km` in `records.csv` show records near a boundary

#### 6. Phylogenetic Model
- **Assumption**: GTR+Gamma model appropriate for COI sequences (if `--build-tree`)
- **Justification**: GTR allows different substitution rates; Gamma handles rate heterogeneity
- **Limitation**: FastTree is approximate ML (not exhaustive search like RAxML/IQ-TREE)
- **Use case**: Phylogenetic tree is for visualization and preliminary analysis, not deep phylogenetics

### Limitations

#### 1. Custom Geographic Regions Supported
- **Default**: Geographic analysis uses GOaS (Global Oceans and Seas) shapefile for marine taxa
- **Custom shapefiles**: Use `--custom-shp` (or a registered `--shp`) to apply any polygon shapefile (freshwater basins, ecoregions, biomes, etc.)
- **No geographic analysis needed**: Use `--no-geo` flag to skip geographic analysis entirely

#### 2. Coordinate Precision Requirements
- **Requirement**: Region assignment and maps need precise coordinates
- **Problem**: Many BOLD records have country-level or missing coordinates
- **Impact**: Those records are used for haplotypes but have no region; `haplotypes.csv` reports, per haplotype, the members with coordinates, with a region and with a country
- **Quantification**: `coord_status` counts are in `run_provenance.json` and the log

#### 3. Single-Locus Limitation
- **Current**: Pipeline designed for single-locus data (COI)
- **Limitation**: Cannot integrate nuclear markers, multiple mitochondrial genes
- **Rationale**: BOLD primarily stores COI data; most users analyzing COI
- **Extension**: Could be adapted for other markers with appropriate thresholds

#### 4. No Haplotype Network Construction
- **Current**: Phylogenetic tree only (if `--build-tree`)
- **Missing**: Haplotype networks (TCS, median-joining) for population genetics
- **Workaround**: Export with `--export-format popart` and use PopART or Network

#### 5. Few Population Genetics Statistics
- **Included**: Haplotype diversity per species label (`--extended`; label-based and sampling-dependent)
- **Missing**: Fst, AMOVA, nucleotide diversity, Tajima's D, Fu's Fs
- **Workaround**: Use `--export-format all` to export to Arlequin, DnaSP, and PopART formats

#### 6. Tree Rooting
- **Current**: Trees are midpoint-rooted unless an outgroup is given (`--phylo-outgroup-fasta`, `--phylo-outgroup-label` or `--phylo-outgroup-taxon`). A midpoint root assumes roughly clock-like rates; it is a drawing convention, not an inferred root
- **Workaround**: Give an outgroup, or re-root in TreeViewer, FigTree, iTOL, or R (ape package)

#### 7. Performance on Very Large Datasets
- **Tested**: Up to ~10,000 sequences
- **Potential bottlenecks**:
  - MAFFT alignment for >20,000 sequences
  - Pairwise distances between haplotypes (O(n²))
- **Mitigation**: Use `--threads` for parallelization

#### 8. Visualization Scalability
- **Challenge**: Maps with >100 haplotypes become cluttered
- **Current**: Overviews show taxa above 40 haplotypes, color at most 12 items, and each taxon has its own page (see Figures)

### Known Issues and Feature Requests

See [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues) for:
- Bug reports
- Feature requests
- Planned enhancements
- Community contributions

**Top feature requests** (as of v2.0.0):
1. Haplotype network construction (TCS, median-joining)
2. Population genetics statistics integration
3. Multi-locus support
4. Docker containerization
5. Galaxy tool integration
6. Web interface

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

**Problem**: Many or all sequences fail ORF validation

This is the most common issue when running BOLDGenotyper on invertebrates. The pipeline defaults to genetic code 2 (vertebrate mitochondrial), which mistranslates invertebrate COI (e.g. AGA/AGG read as stops) and can make many or all sequences fail ORF validation.

```bash
# Count the QC failure reasons
python -c "import pandas as pd; r = pd.read_csv('{organism}_output/records.csv'); print(r.qc_fail_reason.value_counts())"

# Fix: specify the correct genetic code for your taxon
boldgenotyper data/my_species.tsv --genetic-code 5   # most invertebrates
boldgenotyper data/my_species.tsv --genetic-code 4   # cnidarians

# Common codes:
#   2 = vertebrate mitochondrial (default): fish, sharks, mammals, birds, reptiles
#   5 = invertebrate mitochondrial: molluscs, crustaceans, insects, worms, echinoderms
#   4 = mold/protozoan/coelenterate mitochondrial: cnidarians, some fungi
# Full list: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
```

**Problem**: `No haplotypes identified`

```bash
# Check log file for details
tail -100 {organism}_output/{organism}_pipeline.log

# Possible causes:
# 1. All sequences too short (<200 bp or <70% of the median) or too many Ns (>5%)
# 2. MAFFT alignment failed
# 3. No sequences passed quality control (check --genetic-code for invertebrates)
```

**Problem**: Many records not assigned to a haplotype

```bash
# Status counts and denominators
cat {organism}_output/summary.csv

# Which statuses, and why
python -c "import pandas as pd; r = pd.read_csv('{organism}_output/records.csv'); print(r.status.value_counts()); print(r.status_note.value_counts().head(20))"

# excluded_core_coverage: sequences that do not cover the core region (mixed markers or
#   short fragments); check that all records are from the same marker region
# failed_qc: see qc_fail_reason (often --genetic-code for invertebrates)
# excluded_error_singleton: error-filtered singletons that were tied between haplotypes
#   or had no haplotype within --singleton-distance (see status_note)
```

**Problem**: `FastTree failed` (if using `--build-tree`)

```bash
# Verify FastTree is installed
fasttree 2>&1 | head -1

# Install if missing
conda install -c bioconda fasttree

# Check that haplotype sequences exist
head {organism}_output/haplotypes.fasta

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
# Suspect flags and reasons per haplotype
python -c "import pandas as pd; h = pd.read_csv('{organism}_output/haplotypes.csv'); print(h[h.is_suspect][['haplotype_id','suspect_reasons']])"

# A high --singleton-distance merges more rare variants into their neighbors;
# compare thresholds with boldgenotyper-sweep
```

**Problem**: Many records without a region

```bash
# Coordinate status counts (also in run_provenance.json under data_notes)
python -c "import pandas as pd; r = pd.read_csv('{organism}_output/records.csv'); print(r.coord_status.value_counts())"

# Only records with coord_status ok are assigned a region:
# - centroid: country or province centroid (common in BOLD); lat/lon kept, no region
# - missing / unparseable / zero / out_of_range: no usable coordinates
# These records are still used for haplotypes.

# ok records without a region fall outside every polygon (e.g. on land for GOaS),
# or the shapefile did not load (see the warnings at the end of the log)
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
# - figures: one map page per taxon; with --extended, region and metadata pages per taxon
#   and alignment plots
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
# Records with geographic-quality coordinates and a region
python -c "import pandas as pd; r = pd.read_csv('{organism}_output/records.csv'); print((r.coord_status == 'ok').sum(), r.ocean_basin.notna().sum())"

# If none have a region, geographic analysis was skipped (--no-geo, GOaS not installed)
# or failed: see the warnings summary at the end of the log
```

**Problem**: HTML report shows no figures, or the menus do nothing

```bash
# Figures are linked from the run folder: open report.html from that folder,
# and keep the folder's subfolders (phylogeny/, geography/, extended/) next to it.
# The species and genus menus need JavaScript enabled; no internet connection is needed.
```

### Advanced Command Troubleshooting

**Problem**: `boldgenotyper-sweep` - How do I interpret the elbow plot?

```bash
# The elbow plot shows the haplotype count against threshold
# Look for the "elbow point" where the curve levels off, then check the
# stability table: the share of records whose haplotype changes between
# neighboring thresholds

# Check the recommendations file:
cat parameter_sweep/recommendations.txt

# Examine stability plot:
open parameter_sweep/threshold_stability.pdf
```

**Problem**: `boldgenotyper-sweep` - Samples have high `n_changes` in group membership tracking

```bash
# High n_changes: the record's co-members change across thresholds
# These records may be:
# 1. Error-filtered singletons corrected to different haplotypes
# 2. Low-quality sequences
# 3. Records between closely related haplotypes

# Review samples with n_changes > 2:
# Sort by n_changes descending
head -20 parameter_sweep/group_membership_tracking.csv

# Consider excluding or manually reviewing these samples
```

**Problem**: Comparing a species-level and a family-level run finds no shared records, or no label differences I expected

```bash
# 1. Both downloads should hold the same records: the species-level records
#    should be a subset of the family-level download (same download date)
# 2. Match records by processid, not haplotype IDs (IDs differ between runs)
# 3. A record can pass QC in one run and fail core coverage in the other;
#    compare status and status_note in both records.csv files
```
See [`docs/guides/COMPARATIVE_ANALYSIS_GUIDE.md`](docs/guides/COMPARATIVE_ANALYSIS_GUIDE.md#troubleshooting).

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

The repository includes an example analysis of the scalloped hammerhead shark
(*Sphyrna lewini*) in `data/Sphyrna_lewini/`. Its outputs were produced with
BOLDGenotyper 1.0.0 and use the v1 layout (see the
[migration guide](docs/guides/MIGRATION_v1_to_v2.md)).

### Reproducing the Example with 2.0.0

```bash
cd BOLDGenotyper

# Write to a new folder (the included v1 outputs stay as they are)
boldgenotyper data/Sphyrna_lewini_input.tsv --organism Sphyrna_lewini \
  --singleton-distance 0.015 --build-tree --output Sphyrna_lewini_v2

# Start with the run README, then open the report from the run folder
less Sphyrna_lewini_v2/README.md
open Sphyrna_lewini_v2/report.html
```

### Key Findings from the v1.0.0 Example

- **Input**: 937 *Sphyrna lewini* COI sequences from BOLD
- **Haplotypes**: 26 ESV haplotypes identified
- **Geographic Distribution**: Samples from North Atlantic, South Atlantic, North Pacific, South Pacific, and Indian Oceans
- **Biogeographic Pattern**: Distinct haplotypes show basin-specific distributions

Re-running with 2.0.0 gives different haplotype IDs and status counts (record
statuses, error-singleton correction and coordinate parsing changed); compare
haplotypes by sequence.

---

## Citation

If you use BOLDGenotyper in your research, please cite:

### Software Citation

```bibtex
@software{boldgenotyper2026,
  author = {Smith, Steph},
  title = {BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis},
  year = {2026},
  version = {2.0.0},
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

# Build the API documentation (docs/build/html/index.html)
sphinx-build -b html docs/source docs/build/html

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
A: Possibly. The ESV-based haplotype discovery requires no clustering threshold and works for any aligned marker, but quality control checks a protein-coding reading frame (`--genetic-code`), so non-coding markers would fail QC. Choose `--singleton-distance` with `boldgenotyper-sweep` for your marker.

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
A: Yes: run with `--plot-kit`, edit `plots/plot_config.yaml`, and regenerate (see [Customizing Plots](#customizing-plots-for-publication)). The HTML report also allows filtering and export.

**Q: What if my species is not in BOLD?**
A: You can use any FASTA file with similar format. See documentation for custom input preparation (future feature).

**Q: How do I export data for other software (Arlequin, POPART, etc.)?**
A: Use `--export-format arlequin popart dnasp generic` (or `all`); see [`docs/guides/POPGEN_EXPORT_GUIDE.md`](docs/guides/POPGEN_EXPORT_GUIDE.md). `records.csv` and `haplotypes.fasta` serve any other tool.

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

### v2.0.0

A status for every record, error-singleton correction, stable haplotype keys
(`seq_md5`) and count-free IDs (`H0001`), coordinates parsed for every record,
label conflicts with evidence, `run_provenance.json`, a documented run folder
(`records.csv`, `haplotypes.csv`, run `README.md`) with optional `--extended`,
`--export-format` and `--plot-kit` outputs, and a reworked query tool. Outputs and
flags changed: see [CHANGELOG.md](CHANGELOG.md) and the
[migration guide](docs/guides/MIGRATION_v1_to_v2.md).

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
