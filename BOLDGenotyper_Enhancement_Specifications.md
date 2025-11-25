# BOLDGenotyper Enhancement Specifications
## Feature Additions Based on *Sphyrna lewini* Analysis

**Document Version:** 1.0  
**Date:** November 25, 2024  
**Purpose:** Comprehensive specifications for enhancing BOLDGenotyper based on real-world comparative taxonomic analysis

---

## Executive Summary

This document specifies enhancements to BOLDGenotyper identified through comparative analysis of *Sphyrna lewini* (species-level, 617 samples) versus Sphyrnidae (family-level, 995 samples). The primary finding—that hierarchical family-level analysis enables detection of ~10% database contamination invisible to single-species queries—should be made easily accessible to all users.

**Key Insight:** The methodological innovation (hierarchical taxonomic comparison for quality control) required manually running BOLDGenotyper twice and comparing outputs. This should be a one-command operation.

**Implementation Priority:** Features are categorized as HIGH (essential), MEDIUM (important), or LOW (nice-to-have) based on impact on reproducibility, usability, and scientific rigor.

---

## Table of Contents

1. [HIGH PRIORITY: Essential for Reproducibility](#high-priority-essential-for-reproducibility)
   - [Comparative Analysis Module](#1-comparative-analysis-module)
   - [Plot Regeneration Kit](#2-plot-regeneration-kit)
   - [Parameter Sweep Tool](#3-parameter-sweep-tool)
2. [MEDIUM PRIORITY: Enhanced Quality Control](#medium-priority-enhanced-quality-control)
   - [Enhanced Contamination Detection](#4-enhanced-contamination-detection-report)
   - [Divergence Analysis Module](#5-divergence-analysis-module)
   - [Geographic Analysis Enhancement](#6-geographic-analysis-enhancement)
3. [MEDIUM-LOW PRIORITY: Export & Integration](#medium-low-priority-export--integration)
   - [Population Genetics Export Formats](#7-population-genetics-export-formats)
   - [Metadata Enrichment Tools](#8-metadata-enrichment-tools)
4. [LOW PRIORITY: Nice-to-Have](#low-priority-nice-to-have)
   - [Interactive Dashboard](#9-interactive-dashboard-enhancement)
   - [Citation Helper](#10-citation-helper)
5. [Quick Wins: Minimal Changes, High Impact](#quick-wins-minimal-changes-high-impact)
6. [Implementation Guidance](#implementation-guidance)

---

# HIGH PRIORITY: Essential for Reproducibility

## 1. Comparative Analysis Module

### Overview
Enable users to easily compare species-level and family/genus-level analyses to detect database contamination through hierarchical taxonomic context.

### Command-Line Interface

**New Tool:** `boldgenotyper-compare`

```bash
boldgenotyper-compare \
  --species-level Sphyrna_lewini.tsv \
  --family-level Sphyrnidae.tsv \
  --output comparative_analysis/ \
  --generate-reassignment-table \
  --threshold 0.015
```

**Alternative: Integrated Flag**
```bash
boldgenotyper Sphyrnidae.tsv \
  --compare-to Sphyrna_lewini_output/ \
  --output comparison_with_species_level/
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--species-level` | Path | Required | TSV file or previous output directory for species-level analysis |
| `--family-level` | Path | Required | TSV file or previous output directory for family-level analysis |
| `--output` | Path | `./comparative_analysis/` | Output directory for comparison results |
| `--generate-reassignment-table` | Flag | False | Create sample-level reassignment table (Table S4 equivalent) |
| `--threshold` | Float | 0.015 | Clustering threshold (must match if comparing existing outputs) |
| `--majority-threshold` | Float | 0.70 | Threshold for species-level assignment |

### Output Directory Structure

```
comparative_analysis/
├── README.md                              # Overview of comparison
├── comparison_summary.csv                 # High-level metrics comparison
├── genotype_crosswalk.csv                # How species groups map to family groups
├── contamination_analysis.csv            # Mixed-species groups identified
├── sample_reassignments.csv              # Detailed reassignment table (if --generate-reassignment-table)
├── contamination_heatmap.pdf             # Visual: groups × species
├── identity_comparison.pdf               # Distributions: species vs family
├── venn_diagram.pdf                      # Sample overlap visualization
└── reports/
    ├── contamination_report.html         # Interactive contamination explorer
    └── methods_text.md                   # Ready-to-paste methods section
```

### Output File Specifications

#### A. `comparison_summary.csv`

**Purpose:** High-level comparison of analysis outcomes

**Columns:**
```csv
metric,species_level,family_level,difference,pct_change
total_samples,617,995,+378,+61.3
consensus_groups,10,17,+7,+70.0
species_detected,1,8,+7,+700.0
mean_identity,98.7,95.2,-3.5,-3.5
groups_with_100pct_majority,10,10,0,0.0
mixed_species_groups,0,7,+7,N/A
potential_misidentifications,0,66,+66,N/A
samples_in_mixed_groups,0,641,+641,N/A
contamination_rate_in_major_groups,0.0,10.3,+10.3,N/A
```

**Key Metrics to Include:**
- Total samples analyzed
- Number of consensus groups formed
- Number of species detected
- Mean assignment identity
- Number of groups with 100% majority
- Number of mixed-species groups
- Potential misidentifications detected
- Contamination rate calculation

#### B. `genotype_crosswalk.csv`

**Purpose:** Track how species-level groups split/merge in family-level analysis

**Columns:**
```csv
species_group,species_group_size,family_group,family_group_size,n_samples_mapped,pct_of_species_group,pct_of_family_group,mapping_type
c1_n211,211,c15_n386,386,195,92.4,50.5,majority
c1_n211,211,c18_n255,255,14,6.6,5.5,minor
c1_n211,211,c16_n86,86,2,0.9,2.3,contamination
c3_n146,146,c18_n255,255,140,95.9,54.9,majority
c3_n146,146,c15_n386,386,6,4.1,1.6,minor
```

**Mapping Types:**
- `majority`: >50% of species-level group maps to this family-level group
- `minor`: <50% but >5% maps here
- `contamination`: <5% maps here (likely misidentification)
- `singleton`: Only 1-2 samples
- `split`: Species group split evenly across multiple family groups

**Use Case:** Understand whether species-level groups are genuinely homogeneous or artificially lumping distinct taxa.

#### C. `contamination_analysis.csv`

**Purpose:** Detailed breakdown of mixed-species consensus groups

**Columns:**
```csv
family_group,total_samples,n_species,majority_species,majority_count,majority_pct,contaminating_species,contamination_count,contamination_pct,contamination_breakdown
c15_n386,386,4,Sphyrna lewini,344,89.1,"S. zygaena, S. tiburo, S. mokarran",42,10.9,"S. zygaena: 30 (7.8%), S. tiburo: 6 (1.6%), S. mokarran: 4 (1.0%)"
c18_n255,255,3,Sphyrna lewini,231,90.6,"S. zygaena, S. mokarran",24,9.4,"S. zygaena: 18 (7.1%), S. mokarran: 3 (1.2%)"
c16_n86,86,2,Sphyrna zygaena,81,94.2,Sphyrna lewini,5,5.8,"S. lewini: 5 (5.8%)"
c8_n98,98,2,Sphyrna mokarran,97,99.0,Sphyrna lewini,1,1.0,"S. lewini: 1 (1.0%)"
c1_n75,75,3,Sphyrna tiburo,73,97.3,"S. mokarran, S. lewini",2,2.7,"S. mokarran: 1 (1.3%), S. lewini: 1 (1.3%)"
c6_n12,12,2,Eusphyra blochii,10,83.3,Sphyrna zygaena,2,16.7,"S. zygaena: 2 (16.7%)"
c3_n3,3,2,Sphyrna tudes,2,66.7,Sphyrna zygaena,1,33.3,"S. zygaena: 1 (33.3%)"
```

**Summary Statistics to Include:**
```csv
summary_metric,value
total_mixed_groups,7
total_pure_groups,10
pct_groups_contaminated,41.2
total_samples_in_mixed_groups,815
total_contaminating_samples,77
overall_contamination_rate,9.4
most_common_contaminant,Sphyrna zygaena
contaminant_sample_count,57
```

#### D. `sample_reassignments.csv` (Table S4 equivalent)

**Purpose:** Sample-level documentation of potential misidentifications with supporting evidence

**Columns:**
```csv
processid,database_species_label,family_consensus_group,family_group_majority_species,family_group_majority_pct,assignment_identity,runner_up_group,runner_up_identity,is_tie,depositor_notes,reassignment_confidence,suggested_species
ANGBF11456-15,Sphyrna lewini,c16_n86,Sphyrna zygaena,94.2,98.5,c15_n386,96.2,FALSE,"",high,Sphyrna zygaena
GBSHC1234-18,Sphyrna zygaena,c15_n386,Sphyrna lewini,89.1,97.8,c16_n86,95.1,FALSE,"potential cryptic sp",medium,Sphyrna lewini
ANGBF22789-16,Sphyrna lewini,c15_n386,Sphyrna lewini,89.1,99.2,c18_n255,97.5,FALSE,"",high,Sphyrna lewini - CORRECT
```

**Reassignment Confidence Levels:**
- `high`: >95% identity, clear majority species (>90%), not flagged as tie
- `medium`: 90-95% identity, moderate majority (70-90%), or has depositor uncertainty notes
- `low`: <90% identity, weak majority (<70%), or flagged as tie
- `ambiguous`: Multiple equally good assignments, requires manual review

**Notes Column Extraction:**
Automatically search depositor notes for keywords indicating uncertainty:
- "cryptic", "uncertain", "juvenile", "damaged", "cf.", "aff.", "sp.", "complex"

#### E. `contamination_heatmap.pdf`

**Visual Specifications:**
- **Type:** Heatmap
- **Rows:** Family-level consensus groups (sorted by size)
- **Columns:** Species detected
- **Cell color:** Intensity based on sample count (log scale recommended)
- **Annotations:** Cell values showing counts
- **Highlight:** Cells where species ≠ group majority (border or asterisk)
- **Color scheme:** Sequential (e.g., Blues) for within-species; diverging for mixed

**Purpose:** Instantly visualize which groups are contaminated and by which species.

#### F. `identity_comparison.pdf`

**Visual Specifications:**
- **Type:** Overlapping histograms or violin plots
- **X-axis:** Assignment identity (%)
- **Y-axis:** Frequency or density
- **Colors:** Species-level (blue), Family-level (red)
- **Statistics:** Mean, median, SD annotated
- **Interpretation guide:** Text box explaining why family-level may have lower identity

**Purpose:** Show that high identity alone doesn't guarantee data quality.

#### G. `venn_diagram.pdf`

**Visual Specifications:**
- **Type:** Venn diagram or UpSet plot
- **Sets:** 
  - Species-level assigned samples
  - Family-level assigned samples
  - Samples with database label matching assignment
  - Samples flagged for reassignment
- **Annotations:** Sample counts in each intersection

**Purpose:** Visualize overlap and discrepancies between analyses.

#### H. `reports/methods_text.md`

**Content:** Auto-generated methods section describing the comparative analysis

**Template:**
```markdown
### Comparative Taxonomic Analysis

To assess database quality and detect potential misidentifications, we performed 
hierarchical taxonomic analysis by comparing species-level and family-level 
clustering approaches. Species-level analysis queried BOLD for all samples 
labeled "[SPECIES]" (n=[N] samples), while family-level analysis retrieved all 
[FAMILY] samples (n=[N] samples representing [N] species). Both analyses used 
identical BOLDGenotyper parameters (clustering threshold=[THRESHOLD], assignment 
threshold=[THRESHOLD]).

Species-level analysis generated [N] consensus genotypes, all assigned to [SPECIES] 
with 100% majority fractions, superficially suggesting homogeneous, high-quality 
data. However, family-level comparative analysis revealed [N] of [N] consensus 
groups (XX%) contained multiple species, indicating taxonomic heterogeneity 
invisible to single-species queries. The two dominant [SPECIES] genotypes 
([GROUP1] and [GROUP2]) exhibited [XX]% and [XX]% majority fractions respectively, 
revealing ~XX% contamination by morphologically similar congeners (primarily 
[CONTAMINANT SPECIES]).

Samples occurring in mixed-species consensus groups were evaluated for potential 
misidentification by comparing database-reported species labels against the 
majority species of their genetic cluster. Samples in [SPECIES]-majority groups 
with non-[SPECIES] labels were classified as likely misidentifications. Complete 
sample-level reassignments with supporting evidence are documented in Table SX.
```

### Implementation Considerations

**Performance:**
- Comparative analysis should reuse existing outputs when available (check for `{organism}_annotated.csv`)
- Don't re-run clustering if parameters match
- Only regenerate necessary comparison tables

**Error Handling:**
- Verify both analyses used identical parameters
- Warn if clustering thresholds differ
- Handle missing consensus groups gracefully
- Validate that species-level is subset of family-level (or explain if not)

**User Experience:**
- Progress bar for comparison steps
- Summary printed to console highlighting key findings
- Automatic opening of HTML report (optional flag)

**Testing:**
- Unit tests with small mock datasets
- Integration test with actual Sphyrnidae data
- Verify numbers match published S. lewini analysis

---

## 2. Plot Regeneration Kit

### Overview
Provide users with raw data and scripts to regenerate all plots with custom styling for publication.

### Implementation Strategy

Add `--export-plot-data` flag to main pipeline:
```bash
boldgenotyper Sphyrna_lewini.tsv --build-tree --export-plot-data
```

Creates `plots/` directory within output:

```
{organism}_output/
└── plots/
    ├── README.md                        # Instructions for regeneration
    ├── plot_config.yaml                 # All plotting parameters
    ├── data/                            # Raw data for each plot
    │   ├── distribution_map.csv
    │   ├── distribution_bar_relative.csv
    │   ├── distribution_bar_absolute.csv
    │   ├── identity_distribution.csv
    │   ├── tree_data.nwk
    │   ├── tree_tip_labels.csv
    │   └── genotype_colors.csv
    ├── scripts/                         # Regeneration scripts
    │   ├── regenerate_all.sh           # Convenience script
    │   ├── regenerate_map.R            # R script for maps
    │   ├── regenerate_bars.R           # R script for bar charts
    │   ├── regenerate_identity.R       # R script for histograms
    │   ├── regenerate_tree.R           # R script for tree
    │   └── requirements.txt            # R package dependencies
    └── examples/                        # Example modifications
        ├── custom_colors_example.yaml
        ├── publication_style_example.yaml
        └── multi_panel_figure_example.R
```

### File Specifications

#### `README.md`

**Content:**
```markdown
# Plot Regeneration Guide

## Quick Start

1. Modify `plot_config.yaml` to customize colors, sizes, labels
2. Run `bash scripts/regenerate_all.sh` to recreate all plots
3. Find updated PDFs in parent `visualization/` directory

## Individual Plot Regeneration

```bash
# Regenerate only the distribution map
Rscript scripts/regenerate_map.R

# Regenerate with custom config
Rscript scripts/regenerate_bars.R --config custom_style.yaml
```

## Requirements

Install required R packages:
```R
install.packages(c("ggplot2", "sf", "rnaturalearth", "ape", "ggtree"))
```

## Customization Examples

See `examples/` for common modifications:
- Changing color palettes
- Adjusting map projections
- Combining multiple plots
- Publication-ready formatting
```

#### `plot_config.yaml`

**Structure:**
```yaml
# General settings
output_format: ["pdf", "png", "svg"]
dpi: 300
width_inches: 10
height_inches: 8

# Color palette for genotypes
colors:
  c15_n386: "#E41A1C"
  c18_n255: "#377EB8"
  c14_n10: "#4DAF4A"
  c16_n86: "#984EA3"
  c8_n98: "#FF7F00"
  c1_n75: "#FFFF33"
  c4_n22: "#A65628"

# Distribution map settings
map:
  projection: "robinson"  # robinson, mercator, mollweide
  center_longitude: 0
  show_country_borders: true
  border_color: "gray70"
  border_width: 0.3
  ocean_color: "#E8F4F8"
  land_color: "#F5F5F5"
  point_alpha: 0.7
  point_size_range: [2, 10]  # min, max
  point_stroke: 0.5
  show_legend: true
  legend_position: "right"
  legend_title: "Genotype"

# Bar chart settings
bars:
  orientation: "vertical"  # vertical or horizontal
  bar_width: 0.8
  show_percentages: true
  percentage_size: 3
  show_sample_counts: true
  facet_scales: "free_y"  # free_y, free_x, free
  axis_text_angle: 45
  axis_text_size: 10
  color_palette_type: "discrete"  # discrete or gradient

# Identity distribution settings
identity:
  binwidth: 0.5  # percent
  show_mean: true
  show_median: true
  show_density: true
  density_alpha: 0.3
  stat_line_color: "red"
  stat_line_type: "dashed"
  x_limits: [50, 100]
  x_breaks: [50, 60, 70, 80, 90, 100]

# Tree settings
tree:
  layout: "rectangular"  # rectangular, circular, fan
  show_bootstrap: true
  bootstrap_threshold: 70  # Only show values above this
  bootstrap_size: 3
  tip_label_size: 3
  tip_label_offset: 0.001
  branch_width: 0.5
  show_scale_bar: true
  highlight_groups: ["c15_n386", "c18_n255", "c14_n10"]  # Optional
  highlight_color: "#FFE5E5"
```

#### `data/distribution_map.csv`

**Columns:**
```csv
processid,genotype,species,lat,lon,ocean_basin,n_at_location,is_centroid
ANGBF11456-15,c15_n386,Sphyrna lewini,-25.3456,35.2341,Indian Ocean,1,FALSE
GBSHC1234-18,c15_n386,Sphyrna lewini,-25.3456,35.2341,Indian Ocean,1,FALSE
ANGBF22789-16,c18_n255,Sphyrna lewini,22.4123,114.1234,South China Sea,1,FALSE
```

**Note:** `n_at_location` enables stacking/sizing when multiple samples at same coordinates

#### `data/distribution_bar_relative.csv`

**Columns:**
```csv
ocean_basin,genotype,species,n_samples,pct_of_basin
Indian Ocean,c15_n386,Sphyrna lewini,20,74.1
Indian Ocean,c18_n255,Sphyrna lewini,1,3.7
South China Sea,c18_n255,Sphyrna lewini,37,92.5
South China Sea,c15_n386,Sphyrna lewini,2,5.0
```

#### `data/distribution_bar_absolute.csv`

**Columns:**
```csv
ocean_basin,genotype,species,n_samples
Indian Ocean,c15_n386,Sphyrna lewini,20
Indian Ocean,c18_n255,Sphyrna lewini,1
South China Sea,c18_n255,Sphyrna lewini,37
South China Sea,c15_n386,Sphyrna lewini,2
```

#### `data/identity_distribution.csv`

**Columns:**
```csv
processid,genotype,species,identity_pct,is_tie,is_low_confidence
ANGBF11456-15,c15_n386,Sphyrna lewini,99.2,FALSE,FALSE
GBSHC1234-18,c15_n386,Sphyrna lewini,97.8,FALSE,FALSE
ANGBF22789-16,c18_n255,Sphyrna lewini,98.5,FALSE,FALSE
```

#### `data/genotype_colors.csv`

**Columns:**
```csv
genotype,color_hex,species,n_samples,is_primary_for_species
c15_n386,#E41A1C,Sphyrna lewini,386,TRUE
c18_n255,#377EB8,Sphyrna lewini,255,TRUE
c14_n10,#4DAF4A,Sphyrna lewini,10,FALSE
c16_n86,#984EA3,Sphyrna zygaena,86,TRUE
```

#### `scripts/regenerate_map.R`

**Template Structure:**
```r
#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(yaml)
library(dplyr)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) > 0, args[1], "../plot_config.yaml")

# Load configuration
config <- yaml::read_yaml(config_file)

# Load data
map_data <- read.csv("../data/distribution_map.csv", stringsAsFactors = FALSE)
colors <- read.csv("../data/genotype_colors.csv", stringsAsFactors = FALSE)

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create map
p <- ggplot() +
  geom_sf(data = world, 
          fill = config$map$land_color,
          color = config$map$border_color,
          size = config$map$border_width) +
  geom_point(data = map_data,
             aes(x = lon, y = lat, 
                 color = genotype, 
                 size = n_at_location),
             alpha = config$map$point_alpha) +
  scale_color_manual(values = setNames(colors$color_hex, colors$genotype)) +
  scale_size_continuous(range = config$map$point_size_range) +
  coord_sf(crs = st_crs(config$map$projection)) +
  theme_minimal() +
  theme(
    legend.position = config$map$legend_position,
    panel.background = element_rect(fill = config$map$ocean_color)
  ) +
  labs(color = config$map$legend_title, size = "Samples")

# Save outputs
ggsave("../../visualization/distribution_map.pdf", p, 
       width = config$width_inches, height = config$height_inches)
ggsave("../../visualization/distribution_map.png", p, 
       width = config$width_inches, height = config$height_inches, dpi = config$dpi)

message("Map regenerated successfully!")
```

### User Benefits

1. **Easy Customization:** Change colors by editing YAML, not hunting through code
2. **Journal Requirements:** Adjust DPI, dimensions, formats per journal specs
3. **Combined Figures:** Users can create multi-panel figures by modifying scripts
4. **Reproducibility:** Scripts document exactly how plots were made
5. **Learning:** Users can see how plots are constructed and learn R/ggplot2

### Implementation Notes

- Scripts should be well-commented and beginner-friendly
- Include error handling (file not found, missing packages)
- Validate config values (e.g., DPI must be positive integer)
- Provide fallback to defaults if config incomplete
- Test on Windows, Mac, Linux

---

## 3. Parameter Sweep Tool

### Overview
Help users optimize clustering thresholds by testing multiple values and visualizing stability.

### Command-Line Interface

```bash
boldgenotyper-sweep Sphyrna_lewini.tsv \
  --thresholds 0.01,0.015,0.02,0.03,0.05 \
  --parameter clustering-threshold \
  --output parameter_sweep/ \
  --threads 8
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--thresholds` | List[Float] | 0.01,0.015,0.02,0.03,0.05 | Comma-separated threshold values |
| `--parameter` | String | clustering-threshold | Which parameter to sweep (future: similarity-threshold, etc.) |
| `--output` | Path | ./parameter_sweep/ | Output directory |
| `--threads` | Integer | 4 | Parallel threads |
| `--keep-intermediates` | Flag | False | Keep full output for each threshold |

### Output Structure

```
parameter_sweep/
├── README.md
├── sweep_summary.csv                    # Comparison across thresholds
├── threshold_stability.pdf              # Line plots of key metrics
├── group_membership_tracking.csv        # How samples move between thresholds
├── elbow_plot.pdf                       # Identify optimal threshold
├── recommendations.txt                  # Automated recommendations
└── runs/                                # Individual run outputs
    ├── threshold_0.010/
    ├── threshold_0.015/
    ├── threshold_0.020/
    ├── threshold_0.030/
    └── threshold_0.050/
```

### Output File Specifications

#### `sweep_summary.csv`

**Columns:**
```csv
threshold,n_groups,n_singletons,mean_group_size,median_group_size,max_group_size,mean_identity,median_identity,pct_assigned,pct_tie,pct_low_confidence,n_species_detected,n_mixed_groups,avg_group_purity
0.010,25,8,24.7,12,145,99.2,99.5,91.5,2.3,1.2,1,0,100.0
0.015,10,4,61.7,35,211,98.7,99.1,94.3,1.8,0.8,1,0,100.0
0.020,8,2,77.1,58,246,97.8,98.3,95.1,1.5,0.5,1,0,100.0
0.030,5,1,123.4,92,312,96.5,97.1,96.2,1.2,0.3,1,0,100.0
0.050,3,0,205.7,187,398,93.2,94.8,97.8,0.8,0.2,1,0,100.0
```

**Metrics Explained:**
- `n_groups`: Number of consensus groups formed
- `n_singletons`: Groups with only 1 sample
- `mean/median_group_size`: Central tendency of group sizes
- `mean/median_identity`: Assignment quality
- `pct_assigned`: % of samples successfully assigned
- `pct_tie/low_confidence`: Quality flags
- `n_mixed_groups`: Groups with multiple species (if family-level)
- `avg_group_purity`: Average majority fraction across groups

#### `threshold_stability.pdf`

**Visual Specifications:**
- **Type:** Multi-panel line plot
- **Panel 1:** Number of groups (Y) vs threshold (X)
  - Show "elbow" where grouping stabilizes
  - Mark recommended range
- **Panel 2:** Mean identity (Y) vs threshold (X)
  - Show where identity drops below acceptable
- **Panel 3:** % assigned (Y) vs threshold (X)
  - Show trade-off between stringency and coverage
- **Panel 4:** Number of singletons (Y) vs threshold (X)
  - Too many singletons = over-splitting

#### `group_membership_tracking.csv`

**Purpose:** Track which samples change groups as threshold varies

**Columns:**
```csv
processid,species,t_0.010,t_0.015,t_0.020,t_0.030,t_0.050,n_changes,stability_score
ANGBF11456,Sphyrna lewini,c1,c1,c1,c1_merged,c1_merged,1,high
GBSHC1234,Sphyrna lewini,c5,c5,c1_merged,c1_merged,c1_merged,2,medium
ANGBF22789,Sphyrna lewini,c12,c8,c8,c8,c1_merged,3,low
```

**Stability Score:**
- `high`: Same group across all thresholds (or only merges)
- `medium`: 1-2 changes
- `low`: 3+ changes (unstable assignment)

#### `elbow_plot.pdf`

**Visual Specifications:**
- **Type:** Scatter + line plot
- **X-axis:** Clustering threshold
- **Y-axis:** Number of groups (primary)
- **Secondary Y-axis:** Mean identity (optional)
- **Overlay:** Rate of change (derivative)
- **Annotation:** "Elbow" point marked with vertical line
- **Recommendation zone:** Shaded region showing optimal range

#### `recommendations.txt`

**Template:**
```
Parameter Sweep Analysis Results
=================================

Dataset: Sphyrna_lewini.tsv
Total Samples: 617
Threshold Range Tested: 0.010 - 0.050

SUMMARY:
--------
Tested 5 threshold values across 0.010 to 0.050 range.

KEY FINDINGS:
-------------
1. Grouping stability achieved at threshold >= 0.015
   - At 0.010: 25 groups (over-splitting, many singletons)
   - At 0.015: 10 groups (stable)
   - At 0.020: 8 groups (minor changes)
   - At 0.030+: <5 groups (under-splitting, low identity)

2. Assignment quality optimal between 0.015 - 0.020
   - Mean identity: 98.7% - 97.8%
   - Assignment rate: 94.3% - 95.1%
   - Low false positive rate (<1% low confidence)

3. Group membership highly stable (0.015 - 0.020)
   - 95% of samples remain in same group
   - Changes primarily mergers, not reassignments
   - No unstable samples (>3 changes)

RECOMMENDATIONS:
----------------
PRIMARY: threshold = 0.015
  Rationale: 
  - Stable grouping structure (elbow point)
  - High identity (98.7% mean)
  - Minimal singletons (4/10 groups)
  - Appropriate for species-level analysis

ACCEPTABLE RANGE: 0.015 - 0.020
  Rationale:
  - Minor differences in grouping
  - All maintain >97% identity
  - Similar assignment rates (94-95%)

AVOID:
  - < 0.015: Over-splitting (many small groups)
  - > 0.030: Under-splitting (low identity, lumping divergent lineages)

SPECIES-SPECIFIC NOTES:
------------------------
For elasmobranchs (sharks/rays), COI literature suggests:
  - Intraspecific variation: 0-2% (distance 0-0.02)
  - Interspecific divergence: >3-4% (distance >0.03-0.04)
  
Your optimal threshold (0.015) aligns with species-level analysis.
If analyzing family-level or genus-level, consider 0.03-0.05.

NEXT STEPS:
-----------
1. Run full analysis with recommended threshold:
   boldgenotyper Sphyrna_lewini.tsv --clustering-threshold 0.015

2. If results seem over/under-split, adjust within acceptable range

3. Document threshold choice in methods with justification
```

### Implementation Considerations

**Parallelization:**
- Run threshold tests in parallel
- Each threshold is independent
- Collect results at end

**Efficiency:**
- Reuse alignments across thresholds (MAFFT is slow)
- Only re-run clustering and downstream steps
- Cache consensus sequences if possible

**User Guidance:**
- Print real-time progress for each threshold
- Estimate time remaining
- Auto-open recommendations.txt when done

---

# MEDIUM PRIORITY: Enhanced Quality Control

## 4. Enhanced Contamination Detection Report

### Overview
Automatically identify and flag contamination patterns that currently require manual inspection.

### Implementation

Add to existing pipeline output (no new command needed):

```
{organism}_output/
├── quality_control/                     # NEW DIRECTORY
│   ├── README.md
│   ├── mixed_species_summary.csv       # Groups with multiple species
│   ├── contamination_heatmap.pdf       # Visual overview
│   ├── depositor_flags_summary.csv     # Samples with uncertainty notes
│   ├── potential_misidentifications.csv # Flagged samples
│   └── purity_distribution.pdf         # Distribution of majority fractions
```

### File Specifications

#### `mixed_species_summary.csv`

**Columns:**
```csv
consensus_group,n_samples,n_species,primary_species,primary_pct,contaminating_species,contamination_count,flags,confidence_level
c15_n386,386,4,Sphyrna lewini,89.1,"S. zygaena, S. tiburo, S. mokarran",42,CONTAMINATED,high
c14_n10,10,1,Sphyrna lewini,100.0,"",0,CLEAN,high
```

**Flag Values:**
- `CLEAN`: 100% single species
- `CONTAMINATED`: <100% but >70% majority
- `AMBIGUOUS`: <70% majority (genus-level assignment)
- `SINGLETON`: Only 1-2 samples
- `REVIEW`: Unusual patterns requiring manual inspection

#### Enhanced Main Output CSV

**Add These Columns:**
```csv
# Existing columns...
processid,genotype,species,identity,...

# NEW COLUMNS:
group_majority_species,      # What most samples in this group are labeled
group_majority_pct,          # Percentage of majority
group_species_count,         # How many species in this group
matches_group_majority,      # TRUE/FALSE - does sample match group majority?
potential_misidentification, # TRUE/FALSE - likely mislabeled?
depositor_uncertainty_flag,  # TRUE/FALSE - has notes indicating uncertainty
depositor_notes,             # Full notes field
misidentification_confidence # high/medium/low/none
```

**Example Rows:**
```csv
processid,genotype,species,identity,group_majority_species,group_majority_pct,group_species_count,matches_group_majority,potential_misidentification,depositor_uncertainty_flag,depositor_notes,misidentification_confidence
ANGBF123,c15_n386,Sphyrna lewini,99.2,Sphyrna lewini,89.1,4,TRUE,FALSE,FALSE,"",none
ANGBF456,c15_n386,Sphyrna zygaena,98.5,Sphyrna lewini,89.1,4,FALSE,TRUE,FALSE,"",high
ANGBF789,c15_n386,Sphyrna lewini,97.8,Sphyrna lewini,89.1,4,TRUE,FALSE,TRUE,"potential cryptic sp",medium
```

#### `depositor_flags_summary.csv`

**Purpose:** Extract and categorize uncertainty notes

**Columns:**
```csv
processid,genotype,species_label,depositor_note,flag_category,pattern_matched
ANGBF789,c15_n386,Sphyrna lewini,"potential cryptic sp",cryptic_species,cryptic
GBSHC123,c16_n86,Sphyrna lewini,"juvenile - uncertain ID",uncertain_id,uncertain
ANGBF456,c18_n255,Sphyrna cf. lewini,"cf. indicates uncertainty",uncertain_id,cf.
```

**Flag Categories:**
- `cryptic_species`: Notes mention "cryptic", "complex"
- `uncertain_id`: Notes mention "uncertain", "tentative", "probable"
- `morphology_issues`: Notes mention "juvenile", "damaged", "incomplete"
- `taxonomy_qualifier`: Uses "cf.", "aff.", "sp."
- `mixed_sample`: Notes mention "mixed", "contamination"

**Pattern Matching Keywords:**
```python
keywords = {
    'cryptic_species': ['cryptic', 'complex', 'cf. gilberti'],
    'uncertain_id': ['uncertain', 'tentative', 'probable', 'likely', 'possibly'],
    'morphology_issues': ['juvenile', 'damaged', 'incomplete', 'poor condition'],
    'taxonomy_qualifier': ['cf.', 'aff.', 'sp.', 'spp.'],
    'mixed_sample': ['mixed', 'contamination', 'multiple']
}
```

#### `contamination_heatmap.pdf`

**Visual Specifications:**
- **Rows:** Consensus groups (Y-axis)
- **Columns:** Species (X-axis)
- **Cells:** Color intensity = log(count + 1)
- **Annotations:** Raw counts in cells
- **Highlighting:** Border around cells where species ≠ group majority
- **Color scheme:** 
  - Diagonal (correct species): Blue scale
  - Off-diagonal (contamination): Red/Orange scale
- **Marginals:** Row/column totals

**Purpose:** Instantly see contamination patterns - which species contaminate which groups.

#### `purity_distribution.pdf`

**Visual Specifications:**
- **Type:** Histogram
- **X-axis:** Majority fraction (0-100%)
- **Y-axis:** Number of groups
- **Bars colored by:**
  - Green: 100% (pure)
  - Yellow: 90-99% (minor contamination)
  - Orange: 70-89% (moderate contamination)
  - Red: <70% (ambiguous)
- **Overlay:** Threshold line at 70%

### Automated Alerts

Add to console output during analysis:

```
⚠️  QUALITY ALERT: 7 of 17 consensus groups contain multiple species
    - c15_n386: 89.1% Sphyrna lewini (42 contaminating samples)
    - c18_n255: 90.6% Sphyrna lewini (24 contaminating samples)
    
⚠️  15 samples have depositor notes indicating uncertainty
    - 8 flagged as "cryptic species"
    - 7 flagged as "uncertain ID"
    
💡 TIP: Consider running family-level analysis for quality control
    See: boldgenotyper-compare --help
```

---

## 5. Divergence Analysis Module

### Overview
Automatically calculate and report divergence statistics between and within groups.

### Implementation

Add to existing output:

```
{organism}_output/
├── divergence_analysis/                 # NEW DIRECTORY
│   ├── README.md
│   ├── pairwise_divergence_matrix.csv
│   ├── divergence_summary.csv
│   ├── within_vs_between_species.csv
│   ├── barcoding_gap.pdf
│   └── divergence_heatmap.pdf
```

### File Specifications

#### `pairwise_divergence_matrix.csv`

**Format:** Square matrix

```csv
,c15_n386_S.lewini,c18_n255_S.lewini,c14_n10_S.lewini,c16_n86_S.zygaena,c8_n98_S.mokarran
c15_n386_S.lewini,0.000,0.038,0.045,0.042,0.085
c18_n255_S.lewini,0.038,0.000,0.041,0.039,0.083
c14_n10_S.lewini,0.045,0.041,0.000,0.048,0.088
c16_n86_S.zygaena,0.042,0.039,0.048,0.000,0.079
c8_n98_S.mokarran,0.085,0.083,0.088,0.079,0.000
```

**Calculation Method:**
- Align consensus sequences with MAFFT
- Calculate uncorrected p-distance (ignoring gaps)
- `divergence = 1 - (matches / informative_positions)`

#### `divergence_summary.csv`

**Purpose:** Summary statistics on divergence

**Columns:**
```csv
comparison_type,group1,group2,species1,species2,mean_divergence,sd_divergence,min_divergence,max_divergence,n_comparisons,interpretation
within_species,c15_n386,c18_n255,S. lewini,S. lewini,0.038,0.003,0.035,0.042,10,high_divergence_within_species
within_species,c15_n386,c14_n10,S. lewini,S. lewini,0.045,0.004,0.040,0.050,10,high_divergence_within_species
between_species,c15_n386,c16_n86,S. lewini,S. zygaena,0.042,0.002,0.038,0.045,100,typical_interspecific
between_species,c15_n386,c8_n98,S. lewini,S. mokarran,0.085,0.003,0.080,0.090,100,high_interspecific
```

**Interpretation Categories:**
- `low_divergence_within_species`: <1% (typical)
- `moderate_divergence_within_species`: 1-2% (structured)
- `high_divergence_within_species`: >2% (potential cryptic species)
- `typical_interspecific`: 3-5%
- `high_interspecific`: >5%

#### `within_vs_between_species.csv`

**Purpose:** Statistical comparison of divergence distributions

**Columns:**
```csv
species,n_within_comparisons,mean_within,sd_within,max_within,n_between_comparisons,mean_between,sd_between,min_between,barcoding_gap,gap_size
Sphyrna lewini,3,0.041,0.004,0.045,200,0.063,0.021,0.038,TRUE,0.022
```

**Barcoding Gap:**
- `TRUE` if `min_between > max_within` (clear separation)
- `FALSE` if overlap exists
- `gap_size = min_between - max_within`

#### `barcoding_gap.pdf`

**Visual Specifications:**
- **Type:** Overlapping histograms or density plots
- **X-axis:** Divergence (%)
- **Y-axis:** Frequency or density
- **Colors:** 
  - Blue: Within-species comparisons
  - Red: Between-species comparisons
- **Vertical lines:** Mean divergence for each distribution
- **Shaded region:** "Barcoding gap" (if exists)
- **Annotations:** Statistics (mean, SD, min, max)

**Purpose:** Visually assess species boundaries - clear gap = good species delimitation

#### `divergence_heatmap.pdf`

**Visual Specifications:**
- **Type:** Heatmap of pairwise divergence matrix
- **Color scale:** Gradient from white (0%) to dark red (10%+)
- **Cell annotations:** Divergence values
- **Clustering:** Dendrogram on rows/columns (hierarchical clustering by divergence)
- **Species labels:** Color-coded bars on margins

**Purpose:** Quickly see relationships among all genotypes

### Console Output

Add to summary:

```
DIVERGENCE ANALYSIS
===================
Within-species divergence (S. lewini):
  - Between c15_n386 and c18_n255: 3.8% (HIGH)
  - Between c15_n386 and c14_n10: 4.5% (HIGH)
  - Between c18_n255 and c14_n10: 4.1% (HIGH)
  
⚠️  WARNING: Within-species divergence >2% suggests potential cryptic lineages
    Consider: 
    1. Nuclear marker validation
    2. Morphological examination
    3. Biogeographic context
    
Between-species divergence:
  - S. lewini vs S. zygaena: 4.2%
  - S. lewini vs S. mokarran: 8.5%
  
Barcoding gap: ABSENT (within-species max = 4.5%, between-species min = 3.8%)
```

---

## 6. Geographic Analysis Enhancement

### Overview
Better handle limited geographic metadata and provide transparency about coverage.

### Implementation

#### `geographic_coverage.csv`

**Columns:**
```csv
genotype,species,total_samples,with_lat_lon,pct_with_coords,with_ocean_basin,pct_with_basin,with_country,pct_with_country,representativeness
c15_n386,Sphyrna lewini,344,94,27.3,27,7.8,312,90.7,moderate
c18_n255,Sphyrna lewini,231,62,26.8,41,17.7,195,84.4,moderate
c14_n10,Sphyrna lewini,10,0,0.0,0,0.0,8,80.0,poor
```

**Representativeness Score:**
- `excellent`: >50% with ocean basin
- `good`: 25-50% with ocean basin
- `moderate`: 10-25% with ocean basin OR >50% with coords
- `poor`: <10% with ocean basin AND <25% with coords
- `very_poor`: <10% with any geographic data

#### Enhanced Basin Assignment

**Add column:** `basin_assignment_confidence`

```csv
processid,lat,lon,ocean_basin,distance_to_boundary_km,basin_confidence,notes
ANGBF123,-25.3,35.2,Indian Ocean,250.5,high,
ANGBF456,25.1,-80.1,Atlantic Ocean,5.2,low,near Gulf of Mexico boundary
ANGBF789,NULL,NULL,Unknown,NA,none,no coordinates available
```

**Confidence Levels:**
- `high`: >50 km from basin boundary
- `medium`: 10-50 km from boundary
- `low`: <10 km from boundary
- `ambiguous`: Equidistant from multiple basins
- `none`: No coordinates or assignment failed

#### Missing Data Report

**File:** `geographic_missing_data_report.txt`

```
Geographic Data Availability Report
====================================

Overall Coverage:
-----------------
Total Samples: 585
With Coordinates: 158 (27.0%)
With Ocean Basin: 68 (11.6%)
With Country Only: 399 (68.2%)
No Geographic Data: 28 (4.8%)

By Genotype:
------------
c15_n386 (n=344):
  ✓ Coordinates: 94 (27.3%)
  ✓ Ocean Basin: 27 (7.8%)
  ⚠ Country Only: 241 (70.1%)
  ✗ No Data: 9 (2.6%)
  → Representativeness: MODERATE

c18_n255 (n=231):
  ✓ Coordinates: 62 (26.8%)
  ✓ Ocean Basin: 41 (17.7%)
  ⚠ Country Only: 164 (71.0%)
  ✗ No Data: 5 (2.2%)
  → Representativeness: MODERATE

c14_n10 (n=10):
  ✗ Coordinates: 0 (0.0%)
  ✗ Ocean Basin: 0 (0.0%)
  ⚠ Country Only: 8 (80.0%)
  ✗ No Data: 2 (20.0%)
  → Representativeness: POOR

Limitations:
------------
1. Ocean basin assignments available for only 11.6% of samples
2. Coordinate-based visualizations rely on 27% of samples
3. Country-level data insufficient for multi-basin countries
4. Geographic interpretations should be considered preliminary

Recommendations:
----------------
1. Present ocean basin data as supplemental (Table S5)
2. Rely primarily on coordinate-based maps (Figure 2)
3. Acknowledge coverage limitations in methods/discussion
4. Consider contacting depositors for missing coordinates
5. Future sampling should prioritize GPS precision
```

### Console Warnings

```
⚠️  GEOGRAPHIC DATA WARNING
    
Only 11.6% of samples have ocean basin assignments
Only 27.0% have precise coordinates

This may limit biogeographic conclusions.
Consider:
1. Focusing interpretation on samples with coordinates
2. Presenting basin analysis as preliminary
3. Acknowledging limitations in manuscript

See: geographic_analysis/geographic_missing_data_report.txt
```

---

# MEDIUM-LOW PRIORITY: Export & Integration

## 7. Population Genetics Export Formats

### Overview
Enable seamless transfer to specialized population genetics software.

### Implementation

Add `--export-format` flag:

```bash
boldgenotyper Sphyrna_lewini.tsv --export-format arlequin,popart,dnasp
```

### Output Structure

```
{organism}_output/
└── exports/
    ├── README.md                        # Format descriptions
    ├── arlequin/
    │   ├── {organism}.arp              # Arlequin project
    │   └── populations.txt             # Population definitions
    ├── popart/
    │   ├── {organism}.nexus            # NEXUS alignment
    │   ├── {organism}.traits           # Trait file
    │   └── populations.csv             # Population mapping
    ├── dnasp/
    │   ├── {organism}.fas              # FASTA with pop labels
    │   └── {organism}.txt              # DnaSP format
    ├── structure/
    │   └── {organism}.str              # Structure format (if SNPs)
    └── generic/
        ├── alignment.fasta             # Simple FASTA
        ├── genotype_membership.csv     # Sample to genotype mapping
        └── haplotypes.csv              # Unique haplotypes
```

### File Format Specifications

#### Arlequin (.arp)

```
[Profile]
    Title="Sphyrna lewini genotypes from BOLD"
    NbSamples=3
    DataType=DNA
    GenotypicData=0
    LocusSeparator=NONE
    MissingData='?'
    
[Data]
    [[Samples]]
    
    SampleName="c15_n386"
    SampleSize=344
    SampleData= {
        1 SAMPLE001 1 ATCGATCG...
        1 SAMPLE002 1 ATCGATCG...
        ...
    }
    
    SampleName="c18_n255"
    SampleSize=231
    SampleData= {
        1 SAMPLE350 1 ATCGTTCG...
        ...
    }
```

#### PopART/POPART (.nexus)

```
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=585;
    TAXLABELS
        c15_n386_001
        c15_n386_002
        ...
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=650;
    FORMAT DATATYPE=DNA MISSING=? GAP=-;
    MATRIX
        c15_n386_001  ATCGATCGATCG...
        c15_n386_002  ATCGATCGATCG...
        ...
    ;
END;

BEGIN TRAITS;
    DIMENSIONS NTRAITS=3;
    FORMAT LABELS=YES MISSING=? SEPARATOR=COMMA;
    TRAITLABELS genotype species ocean_basin;
    MATRIX
        c15_n386_001,c15_n386,Sphyrna_lewini,Indian_Ocean
        c15_n386_002,c15_n386,Sphyrna_lewini,Indian_Ocean
        ...
    ;
END;
```

#### DnaSP (.fas with pop labels)

```
>c15_n386_001 [genotype=c15_n386;species=Sphyrna_lewini;basin=Indian_Ocean]
ATCGATCGATCGATCGATCG...
>c15_n386_002 [genotype=c15_n386;species=Sphyrna_lewini;basin=Indian_Ocean]
ATCGATCGATCGATCGATCG...
>c18_n255_001 [genotype=c18_n255;species=Sphyrna_lewini;basin=South_China_Sea]
ATCGTTCGATCGATCGATCG...
```

#### Generic CSV (`genotype_membership.csv`)

```csv
sample_id,genotype,species,lat,lon,ocean_basin,country,sequence_id,sequence_length,consensus_sequence
ANGBF11456-15,c15_n386,Sphyrna lewini,-25.3456,35.2341,Indian Ocean,South Africa,SEQ001,650,ATCGATCG...
GBSHC1234-18,c15_n386,Sphyrna lewini,-25.3456,35.2341,Indian Ocean,South Africa,SEQ001,650,ATCGATCG...
ANGBF22789-16,c18_n255,Sphyrna lewini,22.4123,114.1234,South China Sea,China,SEQ002,650,ATCGTTCG...
```

### User Benefits

- No manual reformatting required
- Preserves genotype and population structure
- Maintains metadata associations
- Enables advanced population genetics analyses
- Standardized file naming conventions

---

## 8. Metadata Enrichment Tools

### Overview
Allow users to add their own metadata or update geographic assignments.

### Command-Line Interface

```bash
boldgenotyper-enrich Sphyrnidae_annotated.csv \
  --add-metadata my_sampling_data.csv \
  --join-on processid \
  --add-ocean-basins goas_v2_updated.shp \
  --add-grouping sampling_expedition \
  --output enriched_analysis/
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--add-metadata` | Path | CSV with additional columns to merge |
| `--join-on` | String | Column name to join on (default: processid) |
| `--add-ocean-basins` | Path | Updated shapefile for basin assignment |
| `--add-grouping` | String | Column name to use for grouping (e.g., sampling_campaign) |
| `--recalculate-geography` | Flag | Re-run all geographic summaries |
| `--output` | Path | Output directory |

### Use Cases

**1. Add Field Collection Data:**
```csv
# my_sampling_data.csv
processid,sampling_date,collector,depth_m,temperature_c,life_stage
ANGBF11456-15,2018-05-12,Smith,45,18.5,adult
GBSHC1234-18,2019-03-08,Jones,32,20.1,juvenile
```

**2. Update Ocean Basins with New Shapefile:**
- User has more detailed regional boundaries
- Wants to use different basin definitions
- Has custom ecological regions

**3. Add Custom Groupings:**
- Sampling expeditions
- Time periods
- Ecological zones
- Management regions

### Output

```
enriched_analysis/
├── {organism}_enriched.csv              # Original + new columns
├── enrichment_report.txt               # Summary of changes
├── geographic_updates.csv              # Changed basin assignments
└── updated_visualizations/             # Re-generated plots with new data
    ├── distribution_map_by_expedition.pdf
    ├── distribution_bar_by_region.pdf
    └── ...
```

---

# LOW PRIORITY: Nice-to-Have

## 9. Interactive Dashboard Enhancement

### Overview
Upgrade static HTML report to interactive web dashboard.

### Features to Add

1. **Filterable Tables:**
   - Filter by genotype, species, ocean basin, identity range
   - Sort by any column
   - Export filtered subsets to CSV

2. **Interactive Maps:**
   - Zoom and pan
   - Click sample to see metadata
   - Toggle layers (basins, genotypes, countries)
   - Measure distances
   - Select region to filter

3. **Dynamic Plots:**
   - Hover for sample details
   - Click legend to show/hide genotypes
   - Adjust binning for histograms
   - Download plots as PNG/SVG

4. **Threshold Explorer:**
   - Slider to adjust identity threshold
   - Real-time update of assignment rates
   - Show flagged samples as threshold changes

5. **Contamination Inspector:**
   - Highlight mixed-species groups
   - Click to see detailed composition
   - Show potential misidentifications
   - Link to depositor notes

### Technology Stack

- **Backend:** Flask or FastAPI (Python)
- **Frontend:** Plotly Dash or Streamlit (Python)
- **Alternative:** Shiny (R) for R-focused users
- **Maps:** Leaflet.js or Plotly Express
- **Tables:** DataTables.js or AG Grid

### Implementation Priority

LOW because:
- Static HTML already provides most information
- Interactive features add complexity
- Requires local server or hosted deployment
- Nice-to-have, not essential

But if implemented:
- Significantly improves exploratory analysis
- Makes QC review much faster
- Enables non-technical collaborators to explore

---

## 10. Citation Helper

### Overview
Auto-generate properly formatted methods text and citations.

### Command-Line Interface

```bash
boldgenotyper-cite Sphyrnidae_output/ \
  --style plos \
  --output methods_section.md
```

### Output

**methods_section.md:**

````markdown
## Materials and Methods

### Data Acquisition and Sequence Analysis

COI sequence data were downloaded from the Barcode of Life Database (BOLD; 
Ratnasingham & Hebert, 2007) on [DATE]. Sequence clustering, quality control, 
and genotype assignment were performed using BOLDGenotyper v0.1.0 (Smith et 
al., in prep) with the following parameters:

- Clustering threshold: 0.015 (98.5% identity)
- Assignment threshold: 0.50 (50% identity)
- Minimum sequence length: 400 bp
- Maximum N content: 10%

Sequences were aligned with MAFFT v7.490 (Katoh & Standley, 2013) and trimmed 
with trimAl v1.4.1 (Capella-Gutiérrez et al., 2009). Pairwise genetic distances 
were calculated using uncorrected p-distance, and hierarchical clustering was 
performed using average linkage. Genotype assignment used edit distance 
calculation via the edlib algorithm (Šošić & Šikić, 2017).

Phylogenetic trees were constructed using FastTree v2.1.11 (Price et al., 2010) 
with the GTR+Gamma model. Ocean basin assignments were determined using the 
Global Oceans and Seas (GOaS) v01 shapefile (Flanders Marine Institute, 2023).

### References

Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T (2009) trimAl: a tool for 
automated alignment trimming in large-scale phylogenetic analyses. 
Bioinformatics, 25, 1972-1973.

Flanders Marine Institute (2023) Global Oceans and Seas, version 1. 
Available online at https://www.marineregions.org/

Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7: 
improvements in performance and usability. Molecular Biology and Evolution, 30, 
772-780.

Price MN, Dehal PS, Arkin AP (2010) FastTree 2 – Approximately Maximum-Likelihood 
Trees for Large Alignments. PLOS ONE, 5, e9490.

Ratnasingham S, Hebert PDN (2007) BOLD: The Barcode of Life Data System. 
Molecular Ecology Notes, 7, 355-364.

Smith ST et al. (in prep) BOLDGenotyper: Automated COI sequence analysis and 
biogeographic visualization from BOLD database downloads.

Šošić M, Šikić M (2017) Edlib: a C/C++ library for fast, exact sequence alignment 
using edit distance. Bioinformatics, 33, 1394-1395.
````

### Citation Styles Supported

- `plos` - PLOS journals
- `nature` - Nature journals
- `science` - Science format
- `apa` - APA style
- `vancouver` - Vancouver style
- `bibtex` - BibTeX format

### Implementation

Simple Python script that:
1. Reads `pipeline_parameters.json`
2. Extracts software versions from log
3. Formats according to style guide
4. Includes all relevant citations

---

# QUICK WINS: Minimal Changes, High Impact

These can be implemented immediately with minimal code changes:

## 1. Enhanced Main Output CSV

**Add these columns to existing annotated CSV:**

```python
# In taxonomy analysis phase, add:
df['group_majority_species'] = df.groupby('consensus_group')['species'].transform(
    lambda x: x.mode()[0] if len(x.mode()) > 0 else 'Unknown'
)
df['group_majority_pct'] = df.groupby('consensus_group')['species'].transform(
    lambda x: (x == x.mode()[0]).sum() / len(x) * 100 if len(x.mode()) > 0 else 0
)
df['group_species_count'] = df.groupby('consensus_group')['species'].transform('nunique')
df['matches_group_majority'] = df['species'] == df['group_majority_species']
df['potential_misidentification'] = ~df['matches_group_majority'] & (df['group_majority_pct'] > 70)
```

**Impact:** Users can immediately identify potential misidentifications without any additional tools.

## 2. Add `--export-plot-data` Flag

**Implementation:**
```python
if args.export_plot_data:
    plot_data_dir = output_dir / "plots" / "data"
    plot_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Export map data
    map_df.to_csv(plot_data_dir / "distribution_map.csv", index=False)
    
    # Export bar data
    bar_df.to_csv(plot_data_dir / "distribution_bar.csv", index=False)
    
    # Export tree
    shutil.copy(tree_file, plot_data_dir / "tree_data.nwk")
    
    # Copy plot config
    shutil.copy("plot_config_template.yaml", plot_data_dir / "../plot_config.yaml")
    
    # Copy R scripts
    shutil.copytree("scripts/plot_regeneration", plot_data_dir / "../scripts")
```

**Impact:** Immediate publication-quality plot customization without any new tools.

## 3. Add Contamination Warnings to Console

**Implementation:**
```python
# After taxonomy analysis
mixed_groups = taxonomy_df[taxonomy_df['n_species'] > 1]
if len(mixed_groups) > 0:
    print("\n⚠️  QUALITY ALERT: {len(mixed_groups)} consensus groups contain multiple species")
    for _, row in mixed_groups.head(3).iterrows():
        print(f"    - {row['consensus_group']}: {row['majority_pct']:.1f}% {row['majority_species']}")
    print("\n💡 TIP: Consider running family-level analysis for quality control")
```

**Impact:** Users immediately aware of data quality issues.

## 4. Generate `comparison_summary.csv` Automatically

**When both species and family outputs exist:**
```python
def compare_runs(species_dir, family_dir):
    """Automatically generate comparison if both runs detected."""
    species_summary = load_summary(species_dir / "reports/assignment_summary.csv")
    family_summary = load_summary(family_dir / "reports/assignment_summary.csv")
    
    comparison = {
        'metric': [...],
        'species_level': [...],
        'family_level': [...],
        'difference': [...]
    }
    
    pd.DataFrame(comparison).to_csv(family_dir / "comparison_with_species_level.csv")
```

**Impact:** Comparative analysis partially automated without new tools.

---

# IMPLEMENTATION GUIDANCE

## Phase 1: Essential (Weeks 1-4)

**Week 1-2: Comparative Analysis Module**
- Priority: Highest
- Complexity: Medium-High
- Dependencies: None
- Impact: Enables key methodology

**Tasks:**
1. Design `boldgenotyper-compare` CLI
2. Implement comparison metrics
3. Generate comparison tables (A-C)
4. Create visualization functions
5. Write auto-generated methods text
6. Add tests

**Week 3: Plot Regeneration Kit**
- Priority: High
- Complexity: Medium
- Dependencies: None
- Impact: Publication quality figures

**Tasks:**
1. Design `plot_config.yaml` schema
2. Export plot data to CSV
3. Write R regeneration scripts
4. Create example modifications
5. Write documentation (README.md)
6. Test cross-platform

**Week 4: Quick Wins**
- Priority: High
- Complexity: Low
- Dependencies: None
- Impact: Immediate improvements

**Tasks:**
1. Add columns to main CSV
2. Implement `--export-plot-data`
3. Add console warnings
4. Auto-detect and compare runs
5. Update documentation

## Phase 2: Enhanced QC (Weeks 5-8)

**Week 5-6: Contamination Detection**
- Priority: Medium-High
- Complexity: Medium
- Dependencies: Phase 1
- Impact: Automated quality control

**Week 7: Divergence Analysis**
- Priority: Medium
- Complexity: Medium
- Dependencies: None
- Impact: Species delimitation

**Week 8: Geographic Enhancement**
- Priority: Medium
- Complexity: Low-Medium
- Dependencies: None
- Impact: Transparency

## Phase 3: Advanced Features (Weeks 9-12)

**Week 9-10: Parameter Sweep**
- Priority: Medium
- Complexity: Medium-High
- Dependencies: None
- Impact: Optimization help

**Week 11: Export Formats**
- Priority: Medium-Low
- Complexity: Low-Medium
- Dependencies: None
- Impact: Broader user base

**Week 12: Metadata Enrichment**
- Priority: Low
- Complexity: Medium
- Dependencies: None
- Impact: Power users

## Phase 4: Polish (Ongoing)

**Continuous:**
- Interactive dashboard (if desired)
- Citation helper
- Documentation improvements
- User tutorials
- Additional export formats

---

## Testing Strategy

### Unit Tests
- Each output file format
- Comparison calculations
- Divergence calculations
- Plot data export

### Integration Tests
- Full pipeline with new features
- Comparative analysis workflow
- Plot regeneration from exports

### Validation Tests
- Sphyrnidae dataset (known results)
- Match published S. lewini analysis
- Verify all numbers align

### User Testing
- Beta release to marine biology labs
- Collect feedback on usability
- Iterate on interface design

---

## Documentation Needs

### README Updates
- New features overview
- Example workflows
- Updated screenshots

### Tutorials
- "Comparative Analysis Walkthrough"
- "Customizing Plots for Publication"
- "Detecting Database Contamination"
- "Parameter Optimization Guide"

### API Documentation
- New functions and classes
- File format specifications
- Configuration schemas

### Methods Templates
- Ready-to-paste methods sections
- Citation formats
- Parameter reporting

---

## Success Metrics

### Adoption Metrics
- Downloads/installations
- GitHub stars/forks
- Citations in papers
- User forum activity

### Quality Metrics
- Bug reports per release
- Documentation completeness
- Test coverage (>80% target)
- User satisfaction surveys

### Impact Metrics
- Papers using comparative analysis
- Database quality improvements
- Species discoveries enabled
- Community contributions

---

## Long-Term Vision

**BOLDGenotyper should become:**

1. **Standard Tool** for COI analysis from BOLD
2. **Quality Control Gateway** for database mining
3. **Teaching Tool** for phylogeography and barcoding
4. **Community Platform** for sharing workflows and best practices

**Key Principles:**
- Ease of use (one command = complete analysis)
- Reproducibility (all parameters logged, all outputs documented)
- Transparency (quality issues flagged automatically)
- Flexibility (easy to customize, extend, export)

---

## Questions for User Feedback

Before implementation, gather feedback on:

1. **CLI Design:** Separate tools (`boldgenotyper-compare`) or integrated flags?
2. **Default Behavior:** Should contamination detection run automatically?
3. **Output Volume:** Is the proposed output structure too complex?
4. **Plot Customization:** R scripts vs Python vs both?
5. **Export Formats:** Which population genetics tools are priorities?
6. **Documentation:** Video tutorials vs written guides vs both?

---

## Contact & Contribution

**Repository:** [GitHub URL]  
**Issues:** [GitHub Issues URL]  
**Documentation:** [Docs URL]  
**Citation:** Smith et al. (in prep)

**Contributors Welcome!**
- Feature requests
- Bug reports
- Documentation improvements
- New export formats
- Translation to other languages

---

## Document Version History

- v1.0 (2024-11-25): Initial specification based on S. lewini analysis

---

*This document is a living specification. Features and priorities may be adjusted based on user feedback, technical feasibility, and community needs.*
