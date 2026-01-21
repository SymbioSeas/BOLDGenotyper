# BOLDGenotyper: Complete Workflow Documentation for Manuscript Flow Chart

**Version:** 1.0
**Last Updated:** 2025-12-09
**Purpose:** Comprehensive documentation of data processing steps for manuscript flow chart creation

---

## Overview

BOLDGenotyper implements a **haplotype-first workflow** that processes COI barcode sequences through 10 distinct phases, from raw BOLD data to publication-ready visualizations. This document details every processing step, parameter, output, and decision point.

### Workflow Philosophy

1. **Quality-first approach**: Multi-stage quality control at input, alignment, and post-assignment stages
2. **Exact sequence variant (ESV) approach**: Preserves biological variation while filtering technical errors
3. **Geographic precision**: Marks coordinate quality to maximize genotyping coverage while maintaining geographic accuracy
4. **Comprehensive reporting**: Every decision and metric is logged and visualized

---

## PHASE 1: Data Loading and Quality Control

### 1.1: Parse BOLD TSV Metadata

**Input:**
- BOLD TSV file downloaded from boldsystems.org
- Required columns: `processid`, `nucleotides` (or `nuc`), `species`
- Recommended columns: `lat`, `lon`, `country`, `genus`, `bin_uri`

**Process:**
- Validates required columns exist
- Standardizes column names (e.g., `nuc` → `nucleotides`)
- Preserves all original BOLD metadata columns

**Output:**
- `intermediate/01_parsed_metadata.tsv` - Validated metadata

**Parameters:**
- None (standardized BOLD format)

---

### 1.2: Mark Coordinate Quality

**Input:**
- Parsed metadata with lat/lon columns

**Process:**
- **CRITICAL DESIGN DECISION**: Samples are MARKED (not filtered) for coordinate quality
- Creates quality marker columns:
  - `has_centroid_coords`: TRUE if coordinates match country/region centroids
  - `has_missing_coords`: TRUE if lat or lon is missing/NaN
  - `has_zero_coords`: TRUE if coordinates are exactly [0, 0]
  - `has_invalid_coords`: TRUE if lat/lon out of valid range (lat: -90 to 90, lon: -180 to 180)
  - `is_geographic_quality`: TRUE if all above are FALSE (suitable for basin assignment)

**Rationale:**
- Centroid coordinates can span multiple ocean basins (e.g., Mexico could be Pacific or Atlantic)
- However, sequences are still valid for genotyping
- **Solution**: Keep all samples for genotyping, but only assign ocean basins to samples with precise coordinates

**Output:**
- `intermediate/02_marked_metadata.tsv` - Metadata with quality markers
- All samples retained for downstream analysis

**Parameters:**
- `exclude_centroids`: TRUE (hardcoded) - centroid filtering ENABLED for geographic analysis
- `exclude_zero_coords`: TRUE (hardcoded)
- `exclude_missing_coords`: FALSE (hardcoded) - missing coords samples retained for genotyping

**Key Metric:**
- Number of samples with geographic-quality coordinates (logged to console)

---

### 1.3: Assign Ocean Basins

**Input:**
- Metadata with quality markers
- GOaS shapefile (`GOaS_v1_20211214/goas_v01.shp`)

**Process:**
- Filters to samples where `is_geographic_quality == TRUE`
- Performs spatial join: (lat, lon) → ocean basin polygon
- Assigns basin name (e.g., "North Atlantic Ocean", "South Pacific Ocean")
- Samples with low-quality coordinates receive "Unknown" basin
- Merges basin assignments back to full dataset

**Output:**
- `intermediate/geographic/samples_with_ocean_basins.tsv`
- All samples present, with `ocean_basin` column

**Parameters:**
- `goas_shapefile_path`: Path to shapefile (configurable)
- Can skip with `--no-geo` flag (all basins set to "Unknown")

**Behavior if GOaS not available:**
- Logs warning with download instructions
- Sets all `ocean_basin` to "Unknown"
- Pipeline continues (genotyping still works)

**Key Metrics:**
- Samples with ocean basin assigned: X
- Samples marked "Unknown": Y (centroid/missing coordinates)

---

## PHASE 2: Pre-processing Quality Control

### 2.1: Generate Initial FASTA

**Input:**
- Metadata with sequences (`nuc` column)

**Process:**
- Extracts sequences for all samples with valid (non-empty) `nuc` values
- Converts to uppercase
- Skips samples with missing/empty sequences

**Output:**
- Dictionary: `processid` → `sequence`
- Logs number of sequences loaded vs. skipped

**Parameters:**
- None (direct extraction)

**Key Metric:**
- Skipped samples (missing sequences): X

---

### 2.2: Orientation Normalization and ORF Validation

**Input:**
- Raw sequences from step 2.1

**Process:**
1. **Orientation Check** (if `enable_orientation_check == TRUE`):
   - Translates sequence in all 6 reading frames (forward + reverse complement)
   - Computes ORF (open reading frame) metrics for each:
     - Coverage: fraction of sequence covered by longest ORF without internal stop codons
     - Internal stops: number of stop codons within ORF
   - Selects orientation with best ORF (highest coverage, fewest stops)

2. **ORF Validation**:
   - Uses genetic code `mitochondrial_code` (default: 2 = vertebrate mitochondrial)
   - Calculates ORF metrics for chosen orientation
   - Flags sequences with:
     - `orf_coverage` < `orf_min_coverage` (default: 50%)
     - `orf_internal_stops` > `orf_max_internal_stops` (default: 2)

**Output:**
- `quality_control/{organism}_oriented.fasta` - Orientation-corrected sequences
- `quality_control/{organism}_orf_validation.csv` - Per-sample ORF metrics
  - Columns: `processid`, `orientation_chosen`, `orf_coverage`, `orf_internal_stops`, `is_valid_orf`

**Parameters (Adjustable):**
- `mitochondrial_code`: 2 (standard for vertebrates)
- `enable_orientation_check`: TRUE
- `orf_min_coverage`: 0.5 (50% of sequence must be in valid ORF)
- `orf_max_internal_stops`: 2 (permissive for real-world data)

**Rationale:**
- Permissive thresholds accommodate real-world BOLD data quality
- Sequences with very low ORF coverage (<30%) likely NUMTs, pseudogenes, or contamination

**Key Metrics:**
- Sequences corrected to reverse complement: X
- Sequences with invalid ORF: Y

---

### 2.3: Dynamic QC Filtering

**Input:**
- Oriented sequences
- ORF validation results

**Process:**
- Applies multi-tier filtering:

**Tier 1: Raw Sequence Length**
- **Absolute minimum**: Length ≥ `min_raw_length_abs` (default: 200 bp)
- **Median-based minimum**: Length ≥ `median_length * min_raw_length_frac_of_median` (default: 70% of median)
- Median is computed from entire dataset

**Tier 2: Ambiguous Base Content**
- N fraction ≤ `max_raw_N_fraction` (default: 5%)

**Tier 3: ORF Requirement** (if `require_valid_orf == TRUE`)
- Filters sequences where `is_valid_orf == FALSE`

**Output:**
- `quality_control/{organism}_qc_passed.tsv` - Metadata for QC-passed samples
- `quality_control/{organism}_qc_passed.fasta` - QC-passed sequences
- QC summary statistics

**Parameters (Adjustable):**
- `min_raw_length_abs`: 200 bp
- `min_raw_length_frac_of_median`: 0.7 (70%)
- `max_raw_N_fraction`: 0.05 (5%)
- `require_valid_orf`: TRUE (recommended for phylogenetics)

**Key Metrics:**
- Total input sequences: N
- Passed QC: X (X%)
- Failed QC: Y (Y%)

---

## PHASE 3: Haplotype Discovery (ESV Approach)

### 3.1: Identify Haplotypes

**Input:**
- QC-passed sequences (`quality_control/{organism}_qc_passed.fasta`)
- ORF validation data

**Process:**

**Step 3.1.1: Multiple Sequence Alignment (MSA)**
- Aligns all QC-passed sequences with MAFFT
- Algorithm: `mafft_algorithm` (default: "auto" - MAFFT chooses based on dataset size)
- Output: Aligned FASTA file

**Step 3.1.2: Gap Masking**
- Masks alignment columns with >50% gaps (`mask_gap_threshold`)
- Prevents insertion artifacts from fragmenting coverage
- Converts high-gap columns to "N" (masked)

**Step 3.1.3: Core Region Extraction**
- Identifies "core region" covered by ≥80% of sequences (`core_min_coverage`)
- Extracts this shared region from each sequence
- Ensures minimum length: `core_min_length` (default: 150 bp)
- Falls back to 100 bp with warning if 150 bp not achievable

**Step 3.1.4: Identify Exact Sequence Variants (ESVs)**
- Groups sequences by **exact identity** in core region
- Each unique core region sequence = one haplotype
- No distance threshold - true ESV approach
- Records membership: which samples belong to each haplotype

**Step 3.1.5: Singleton Error Filtering**
- **Critical QC step**: Filters singleton haplotypes (n=1 member) likely due to sequencing/PCR errors
- Calculates nearest-neighbor distance for each singleton
- **Removes** singletons with distance ≤ `min_singleton_distance` (default: 0.005 = 0.5%)
- **Rationale**:
  - Sequencing errors (Illumina ~0.1-0.3%, Sanger ~0.1%) create false singletons
  - PCR errors (~0.01-0.1% per bp) compound the issue
  - Empirical data: ~85% of singletons fall within 0.5% of true haplotypes
  - Biological rare variants (>0.5% divergence) are retained

**Step 3.1.6: Suspect Haplotype Flagging** (if `flag_suspect_haplotypes == TRUE`)
- Flags remaining singletons with distance > `max_singleton_distance` (default: 5%)
- These may represent:
  - Contamination
  - Misidentification
  - Cryptic species
  - Geographic relicts
- Flagged but NOT removed (user decision)

**Output:**
- `haplotypes/{organism}_haplotypes.fasta` - Consensus sequences for each haplotype
- `haplotypes/{organism}_haplotype_mapping.csv` - Sample-to-haplotype assignments
  - Columns: `processid`, `haplotype_id`, `haplotype_number`, `n_members`, `is_singleton`
- `haplotypes/{organism}_haplotype_stats.csv` - Per-haplotype statistics
  - Columns: `haplotype_id`, `n_members`, `is_singleton`, `is_suspect`, `nearest_neighbor_distance`

**Parameters (Adjustable):**
- `core_min_coverage`: 0.8 (80% of sequences must cover position)
- `mask_gap_threshold`: 0.5 (columns with >50% gaps masked)
- `core_min_length`: 150 bp (minimum core region length)
- `min_singleton_distance`: 0.005 (0.5% - singletons ≤ this divergence are FILTERED)
- `max_singleton_distance`: 0.05 (5% - singletons > this divergence are FLAGGED)
- `flag_suspect_haplotypes`: FALSE (default - no flagging)

**Hardcoded:**
- ESV identity: 100% (exact match in core region)
- Gap/ambiguous base handling: Ignored in distance calculations

**Key Metrics:**
- Total haplotypes identified: X (after error filtering)
- Singletons: Y
- Suspect haplotypes: Z
- Singleton error filter removed: W singletons

---

## PHASE 4: Haplotype Assignment

### 4.1: Assign Samples to Haplotypes

**Input:**
- Haplotype mapping from Phase 3
- QC-passed metadata

**Process:**
- **Direct mapping**: Uses exact assignments from Phase 3 haplotype discovery
- No distance-based re-assignment (avoids false ties)
- All samples are assigned to their ESV haplotype
- Samples that didn't pass QC (Phase 2) are NOT assigned

**Assignment Statistics:**
- `is_tie`: FALSE for all (no ambiguous assignments in ESV approach)
- `is_low_confidence`: FALSE for all (exact matches)
- `assignment_method`: "ESV_direct_mapping"

**Output:**
- `haplotype_assignments/{organism}_diagnostics.csv`
  - Columns: `processid`, `assigned_haplotype`, `is_tie`, `is_low_confidence`, `assignment_method`, `note`
- `intermediate/haplotype_assignments/{organism}_with_haplotypes.tsv`

**Parameters:**
- None (direct mapping from Phase 3)

**Key Metrics:**
- Assigned samples: X/Y (X%)
- Ties: 0 (ESV method)
- Low-confidence: 0 (ESV method)

---

## PHASE 5: Taxonomy Assignment to Haplotypes

### 5.1: Metadata-Based Taxonomy (Majority Vote)

**Input:**
- Haplotype-assigned metadata
- Haplotype membership information

**Process:**
1. **Group by haplotype**: Aggregate all samples in each haplotype
2. **Count species**: Tally species names for each haplotype
3. **Majority voting**:
   - If one species ≥ `majority_species_threshold` (default: 70%) → assign that species
   - If no majority → assign genus-level or "ambiguous"
4. **Build haplotype labels**: Format as "{species} h{N}" (e.g., "S. lewini h1")

**Output:**
- `taxonomy/{organism}_species_by_haplotype.csv`
  - Columns: `haplotype_id`, `reported_species`, `n_samples`, `fraction`
- `taxonomy/{organism}_haplotype_taxonomy.csv`
  - Columns: `haplotype_id`, `assigned_sp`, `haplotype_sp`, `assignment_level`, `assignment_notes`, `majority_fraction`

**Parameters (Adjustable):**
- `majority_species_threshold`: 0.70 (70% of samples must agree)

**Assignment Levels:**
- `species`: Clear majority (≥70%)
- `genus`: Species conflict, but genus agrees
- `ambiguous`: No clear taxonomic assignment

---

### 5.2: Sequence-Based Taxonomy (Optional, if BLAST/VSEARCH available)

**Input:**
- Haplotype consensus sequences
- Reference database (e.g., BOLD COI database)

**Process:**
1. **Sequence similarity search**: BLAST or VSEARCH against reference
2. **Threshold-based assignment**:
   - Identity ≥ `min_identity_pct` (default: 98.5%) for species-level
   - Query coverage ≥ `min_query_cov_pct` (default: 90%)
3. **Top2 conflict detection**:
   - If best hit and 2nd-best hit differ by < `top2_min_delta_pct` → ambiguous
4. **Conflict resolution**: Chooses between metadata-based and sequence-based

**Output:**
- `taxonomy/{organism}_haplotype_taxonomy_seq.csv` (if sequence-based taxonomy run)
  - Columns: `haplotype_id`, `cluster_seq_sp`, `cluster_seq_level`, `cluster_seq_best_identity`, `cluster_seq_qcov`

**Parameters (Adjustable):**
- `min_identity_pct`: 98.5% (COI species threshold)
- `min_query_cov_pct`: 90%
- `top2_min_delta_pct`: 0.5% (minimum margin for unambiguous call)
- `classifier`: "blastn" or "vsearch"

**Final Taxonomy Decision:**
- Integrates metadata-based and sequence-based assignments
- Prioritizes sequence-based if high confidence (≥98.5% identity, ≥90% coverage)
- Falls back to metadata if sequence match is poor
- Provenance tracked in `tax_provenance` column

---

## PHASE 6: Post-assignment Quality Control

### 6.1: Contamination Detection

**Input:**
- Annotated dataset with haplotype and taxonomy assignments

**Process:**
1. **Per-haplotype contamination checks**:
   - Multi-species haplotypes flagged (species diversity within haplotype)
   - Low-majority species (< 70% agreement) flagged
   - Unexpected species combinations noted

2. **Per-sample purity scores**:
   - Purity = fraction of haplotype members sharing sample's species
   - Low purity (<0.7) suggests potential contamination/mislabeling

**Output:**
- `quality_control/README.md` - QC summary report
- `quality_control/contamination_heatmap.pdf` - Species × haplotype heatmap
- `quality_control/purity_distribution.pdf` - Histogram of sample purity scores
- Enhanced `{organism}_annotated.csv` with QC columns:
  - `sample_purity`: Purity score (0-1)
  - `haplotype_is_multi_species`: Boolean flag
  - `potential_contamination`: Boolean flag

**Parameters:**
- `majority_threshold`: 0.70 (inherited from taxonomy config)

**Quality Alerts (logged to console):**
- Haplotypes with multiple species
- Samples with low purity (<0.7)
- Potential contamination cases

---

## PHASE 6.5: Species-Level Aggregation

### 6.5.1: Aggregate Samples by Species

**Input:**
- Annotated metadata with species assignments

**Process:**
- Groups samples by `assigned_sp` (from haplotype taxonomy)
- Filters to species-level assignments with confidence ≥ `min_confidence` (default: 0.7)
- Creates species-level dataset

**Output:**
- `species_analysis/species_assignments.csv`
  - All columns from annotated metadata
  - Grouped by assigned species

**Parameters:**
- `min_confidence`: 0.7 (minimum majority fraction for species assignment)

---

### 6.5.2: Calculate Species Diversity Metrics

**Input:**
- Species assignments

**Process:**
- **Per-species metrics**:
  - Number of samples
  - Number of haplotypes
  - Number of ocean basins
  - Geographic range (lat/lon bounding box)
- **Diversity indices**:
  - Haplotype richness (# unique haplotypes)
  - Haplotype evenness (Shannon index)

**Output:**
- `species_analysis/species_summary.csv`
  - Columns: `species`, `n_samples`, `n_haplotypes`, `n_basins`, `haplotype_richness`, `haplotype_evenness`
- `species_analysis/species_diversity_metrics.csv`

**Parameters:**
- None (standard diversity calculations)

---

### 6.5.3: Geographic Summary by Species

**Input:**
- Species assignments with ocean basin data

**Process:**
- Aggregates samples by species × ocean basin
- Calculates:
  - Samples per species per basin
  - Haplotypes per species per basin
  - Relative abundance within each species

**Output:**
- `species_analysis/species_geographic_summary.csv`
  - Columns: `species`, `ocean_basin`, `n_samples`, `n_haplotypes`, `relative_abundance`

---

### 6.5.4: Generate Species-Faceted Haplotype Subsets

**Input:**
- Species assignments
- Haplotype data

**Process:**
- For each species with ≥ `min_haplotypes` (default: 2):
  - Extracts all samples assigned to that species
  - Creates species-specific haplotype dataset
  - Saves to separate file for downstream analysis

**Output:**
- `species_analysis/species_facets/{Species}_haplotypes.csv` (one file per species)
  - Subset of annotated data for that species only

**Parameters:**
- `min_haplotypes`: 2 (minimum haplotypes to create species facet)

**Purpose:**
- Enables species-specific analyses (divergence, phylogeny, distribution)
- Facilitates multi-species studies with per-species detail

---

## PHASE 7: Geographic Analysis Enhancement

### 7.1: Geographic Data Quality Assessment

**Input:**
- Annotated dataset with ocean basin assignments
- GOaS spatial data (if available)

**Process:**
1. **Missing Data Report**:
   - Counts samples with missing coordinates
   - Counts samples with centroid coordinates
   - Summarizes ocean basin coverage

2. **Coordinate Precision Analysis**:
   - Calculates decimal places in lat/lon
   - Flags low-precision coordinates

3. **Basin-Level Summaries**:
   - Samples per basin
   - Haplotypes per basin
   - Species per basin (if available)

**Output:**
- `geographic_analysis/geographic_missing_data_report.txt`
- `geographic_analysis/basin_summary.csv`
- Enhanced annotated CSV with geographic metadata

**Parameters:**
- None (summary statistics only)

**Key Metrics:**
- Geographic quality samples: X/Y
- Centroid coordinates: Z
- Basins represented: N

---

## PHASE 8: Phylogenetic Analysis (Optional)

### 8.1: Build Phylogenetic Tree

**Input:**
- Haplotype consensus sequences (`haplotypes/{organism}_haplotypes.fasta`)

**Process:**

**Step 8.1.1: Filter Haplotypes**
- Filters haplotypes by:
  - Minimum sequence length: `min_consensus_length` (default: 150 bp)
  - Minimum cluster size: `min_cluster_size` (default: 3 members)
- Rationale: Small clusters often produce fragmentary consensus sequences

**Step 8.1.2: Multiple Sequence Alignment**
- Aligns filtered haplotypes with MAFFT
- Algorithm: `mafft_algorithm` (default: "auto")

**Step 8.1.3: Alignment Trimming** (if `trim_alignment == TRUE`)
- Trims gappy columns with trimAl
- Method: `trim_method` (default: "gappyout")
- Removes columns with excessive gaps

**Step 8.1.4: Tree Building**
- Constructs maximum-likelihood tree with FastTree
- Model: GTR+Gamma (hardcoded, appropriate for COI)
- Calculates bootstrap support values

**Step 8.1.5: Tree Rooting** (if outgroup specified)
- **Outgroup FASTA** (`outgroup_fasta`): External outgroup sequences
- **Outgroup Label** (`outgroup_label`): Specific tip label in tree
- **Outgroup Taxon** (`outgroup_taxon`): Species/genus name (uses LCA of matching tips)
- If no outgroup: Midpoint rooting applied if `midpoint_root == TRUE`

**Output:**
- `phylogenetic/{organism}_tree.nwk` - Newick format tree
- `phylogenetic/{organism}_tree_relabeled.nwk` - Tree with "{species} h{N}" tip labels
- `intermediate/phylogenetic/{organism}_aligned.fasta` - Alignment
- `intermediate/phylogenetic/{organism}_aligned_trimmed.fasta` - Trimmed alignment (if trimming enabled)
- `intermediate/phylogenetic/{organism}_aligned_relabeled.fasta` - Alignment with readable labels

**Parameters (Adjustable):**
- `build_tree`: FALSE (default - tree building optional)
- `min_consensus_length`: 150 bp (raise to 600 bp for full-length barcodes)
- `min_cluster_size`: 3 members (avoid singleton/doubleton consensuses)
- `trim_alignment`: TRUE
- `trim_method`: "gappyout"
- `outgroup_fasta`: None (optional)
- `outgroup_label`: None (optional)
- `outgroup_taxon`: None (optional)
- `midpoint_root`: TRUE

**Hardcoded:**
- Substitution model: GTR+Gamma
- Tree algorithm: FastTree

**Key Metrics:**
- Haplotypes in tree: X (after filtering)
- Alignment length: Y bp
- Trimmed alignment length: Z bp (if trimming)

---

### 8.2: Create Relabeled Tree

**Input:**
- Tree file (`.nwk`)
- Haplotype taxonomy (`taxonomy/{organism}_haplotype_taxonomy.csv`)

**Process:**
- Replaces haplotype IDs with readable labels from `haplotype_sp` column
- Format: "{species} h{N}" (e.g., "S. lewini h1")
- Enables publication-quality figures in external tree viewers

**Output:**
- `phylogenetic/{organism}_tree_relabeled.nwk`

**Usage:**
- Open in TreeViewer, FigTree, iTOL, or R (ape/ggtree)
- Ready for re-rooting, annotation, and figure customization

---

### 8.3: MSA Visualization (Optional)

**Input:**
- Aligned haplotype sequences
- Phylogenetic tree (for sequence ordering)

**Process:**
- Generates phylogeny-ordered multiple sequence alignment plots
- **Chunking**: Splits long alignments into `chunk_size` bp chunks (default: 100 bp)
- **Sequence selection**: Limits to `max_sequences` (default: 50) if alignment is large
- **Color scheme**: Nucleotide coloring (`color_scheme`, default: "Nucleotide")
- **Annotations**:
  - Consensus sequence bar (if `show_consensus == TRUE`)
  - Sequence logo (if `show_logo == TRUE`)
  - Conservation bar chart (if `show_conservation == TRUE`)

**Output:**
- `phylogenetic/{organism}_msa_chunk_{N}.{format}` (one file per chunk)
  - Formats: PDF, PNG (as specified in `output_formats`)

**Parameters (Adjustable):**
- `enabled`: TRUE (MSA visualization on by default)
- `chunk_size`: 100 bp
- `max_sequences`: 50
- `color_scheme`: "Nucleotide" (options: "Purine/Pyrimidine", "Clustal", "Taylor")
- `show_consensus`: TRUE
- `show_logo`: TRUE
- `show_conservation`: TRUE
- `output_formats`: ["pdf", "png"]

**Dependencies:**
- Requires `pymsaviz` package (optional)
- Gracefully skips if not installed

---

## PHASE 9: Divergence Analysis

### 9.1: Haplotype-Level Pairwise Divergence

**Input:**
- Haplotype consensus sequences
- Haplotype taxonomy

**Process:**
1. **Pairwise Distance Matrix**:
   - Calculates p-distance (uncorrected) between all haplotype pairs
   - Formula: `distance = 1 - (matches / informative_sites)`
   - Gaps and ambiguous bases (N) ignored

2. **Barcoding Gap Analysis**:
   - Separates distances into:
     - **Within-species**: Both haplotypes assigned to same species
     - **Between-species**: Haplotypes from different species
   - Identifies "barcoding gap" (minimum between-species > maximum within-species)

3. **Nearest-Neighbor Distances**:
   - For each haplotype, finds closest relative
   - Useful for identifying outliers

**Output:**
- `divergence_analysis/divergence_matrix.csv` - Full pairwise distance matrix
- `divergence_analysis/barcoding_gap.csv` - Within vs. between-species distances
- `divergence_analysis/barcoding_gap.pdf` - Histogram showing gap
- `divergence_analysis/divergence_heatmap.pdf` - Heatmap of pairwise distances

**Parameters:**
- `distance_metric_primary`: "p" (p-distance, hardcoded)
- `compute_k2p`: FALSE (Kimura 2-parameter distance, optional)

**Key Metrics (logged to console):**
- Mean within-species divergence: X%
- Mean between-species divergence: Y%
- Barcoding gap present: TRUE/FALSE
- Maximum within-species: A%
- Minimum between-species: B%

---

### 9.2: Species-Level Divergence

**Input:**
- Haplotype divergence matrix
- Species summary (from Phase 6.5)

**Process:**
1. **Species-Level Distance Matrix**:
   - Aggregates haplotype distances to species level
   - Uses mean distance between all haplotype pairs from different species

2. **Species Divergence Summary**:
   - Mean divergence between each species pair
   - Within-species divergence (if species has multiple haplotypes)

**Output:**
- `species_analysis/species_divergence_matrix.csv`
- `species_analysis/species_divergence_summary.csv`

**Parameters:**
- None (derived from haplotype divergence)

---

### 9.3: Within-Species Divergence Matrices

**Input:**
- Haplotype divergence matrix
- Species-faceted haplotype subsets (from Phase 6.5.4)

**Process:**
- For each species with ≥ `min_haplotypes` (default: 2):
  - Extracts pairwise distances between haplotypes within that species
  - Creates species-specific divergence matrix

**Output:**
- `species_analysis/species_facets/{Species}_divergence_matrix.csv` (one per species)

**Parameters:**
- `min_haplotypes`: 2

**Purpose:**
- Enables within-species population structure analysis
- Identifies intraspecific lineages

---

## PHASE 10: Visualization

### 10.1: Haplotype-Level Visualizations

#### 10.1.1: Geographic Distribution Map

**Input:**
- Annotated metadata with lat/lon
- Haplotype assignments

**Process:**
- Creates world map with:
  - Sample points color-coded by haplotype (`haplotype_sp`)
  - Point size scaled by abundance (number of samples at location)
  - Ocean basin boundaries (if GOaS available)

**Output:**
- `visualization/{organism}_distribution_map.{format}`
  - Formats: PNG (300 DPI), PDF (vector), SVG (optional)
- `visualization/{organism}_distribution_map_data.json` - Plot data for interactive HTML

**Parameters (Adjustable):**
- `map_projection`: "PlateCarree" (options: "Robinson", "Mollweide", "Mercator")
- `map_figsize`: (12, 6) inches
- `figure_dpi`: 300
- `point_alpha`: 0.8
- `point_size_range`: [3, 12]

**Color Palette:**
- Uses `color_palette` (default: "colorblind")
- First 3 haplotypes use `reference_colors` if specified:
  - Purple: #9D7ABE
  - Teal: #5AB4AC
  - Yellow: #F2CC8F

---

#### 10.1.2: Ocean Basin Abundance (Relative)

**Input:**
- Annotated metadata with ocean basin assignments

**Process:**
- Creates stacked bar chart:
  - X-axis: Ocean basin
  - Y-axis: Relative abundance (0-100%)
  - Bars stacked by haplotype

**Output:**
- `visualization/{organism}_distribution_bar.{format}`
- `visualization/{organism}_distribution_bar_data.json`

**Parameters:**
- `barplot_figsize`: (10, 6) inches
- `bar_orientation`: "vertical" or "horizontal"
- `bar_width`: 0.7

**Interpretation:**
- Shows haplotype composition within each basin
- Normalized to 100% per basin (accounts for sampling bias)

---

#### 10.1.3: Ocean Basin Abundance (Total Counts)

**Input:**
- Annotated metadata with ocean basin assignments

**Process:**
- Creates stacked bar chart:
  - X-axis: Ocean basin
  - Y-axis: Total sample count
  - Bars stacked by haplotype

**Output:**
- `visualization/{organism}_totaldistribution_bar.{format}`
- `visualization/{organism}_totaldistribution_bar_data.json`

**Interpretation:**
- Shows absolute sample counts per basin
- Reveals sampling effort differences across basins

---

### 10.2: Species-Level Visualizations

**Trigger:** Only generated if species-level analysis completed (Phase 6.5)

#### 10.2.1: Species Distribution Map

**Process:**
- Similar to haplotype map, but color-coded by `primary_species`
- Shows species-level biogeography

**Output:**
- `visualization/{organism}_species_distribution_map.{format}`

---

#### 10.2.2: Species Ocean Basin Abundance

**Process:**
- Stacked bar charts (relative and total) by species instead of haplotype

**Output:**
- `visualization/{organism}_species_distribution_bar.{format}` (relative)
- `visualization/{organism}_species_totaldistribution_bar.{format}` (total)

---

### 10.3: Species-Faceted Haplotype Visualizations

**Trigger:** Species-level analysis completed with ≥ 2 species

#### 10.3.1: Species-Faceted Distribution Maps

**Process:**
- Creates separate maps for each species
- Each map shows:
  - Only samples from that species
  - Haplotypes within that species color-coded
  - Zoomed to species' geographic range (with `map_buffer_degrees` buffer, default: 20°)
- Annotation showing samples with missing coordinates

**Output:**
- `visualization/species_facets/{Species}_haplotype_distribution_map.{format}` (one per species)

**Parameters:**
- `map_buffer_degrees`: 20.0 (degrees to extend beyond data bounds)
- `show_unknown_geography_annotation`: TRUE (shows count of samples without coordinates)
- `show_scale_bar`: TRUE (legend showing point size → sample count)

---

#### 10.3.2: Species-Faceted Basin Abundance Charts

**Process:**
- Creates separate bar charts for each species
- Shows haplotype composition within that species across basins

**Output:**
- `visualization/species_facets/{Species}_haplotype_basin_abundance.{format}` (one per species)

---

### 10.4: Haplotype-Specific Visualizations

**Trigger:** If significant haplotypes (high abundance or wide distribution)

#### 10.4.1: Per-Haplotype Distribution Maps

**Process:**
- Creates individual maps for major haplotypes
- Highlights geographic range of single haplotype
- Useful for phylogeography studies

**Output:**
- `visualization/haplotype_plots/haplotype_{ID}_distribution.{format}` (one per major haplotype)

---

## PHASE 11: Report Generation (Optional)

### 11.1: Interactive HTML Report

**Input:**
- All previous outputs
- JSON data files from visualizations

**Process:**
- Generates comprehensive HTML report with:
  - **Summary Dashboard**: Sample counts, haplotype counts, assignment rates
  - **Pipeline Parameters**: All settings used (reproducibility)
  - **Interactive Visualizations**: Plotly.js plots with filtering
    - Filter by haplotype, ocean basin, minimum sample count
    - Toggle between static and interactive views
    - Export plots as PNG/SVG
  - **Methods Section**: Publication-ready text describing workflow
  - **Data Tables**: Assignment status, taxonomy, geography
  - **Download Links**: CSV exports, filtered datasets

**Output:**
- `{organism}_summary_report.html`

**Parameters:**
- `no_report`: FALSE (default - generate report)

**Dependencies:**
- Jinja2 templates
- Plotly.js (CDN, requires internet for interactive features)

---

## Data Flow Summary

```
BOLD TSV (raw sequences)
  ↓
[PHASE 1] Mark coordinate quality → All samples retained
  ↓
[PHASE 2] QC filtering (length, N content, ORF) → QC-passed samples
  ↓
[PHASE 3] Haplotype discovery (ESV + error filtering) → Haplotypes
  ↓
[PHASE 4] Direct assignment → Samples with haplotype IDs
  ↓
[PHASE 5] Taxonomy (metadata + sequence) → Haplotypes with species labels
  ↓
[PHASE 6] QC (contamination detection) → Quality-flagged samples
  ↓
[PHASE 6.5] Species aggregation → Species-level datasets
  ↓
[PHASE 7] Geographic enhancement → Enhanced metadata
  ↓
[PHASE 8] Phylogenetic tree (optional) → Tree files
  ↓
[PHASE 9] Divergence analysis → Distance matrices
  ↓
[PHASE 10] Visualization → Publication-ready figures
  ↓
[PHASE 11] HTML report (optional) → Interactive report
```

---

## Parameter Summary Tables

### Hardcoded Parameters (Cannot be Changed via CLI)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **ESV identity threshold** | 100% (exact match in core region) | True ESV approach preserves all biological variation |
| **Gap handling in distances** | Ignored | Gaps may be alignment artifacts; focus on informative sites |
| **Ambiguous base handling** | Ignored (N excluded from distance calc) | Ns lack information; including them inflates distances |
| **Genetic code** | 2 (vertebrate mitochondrial) | Standard for COI in vertebrates |
| **Phylogenetic model** | GTR+Gamma | Appropriate for COI; allows rate heterogeneity |
| **Tree algorithm** | FastTree | Fast approximate ML; suitable for initial trees |
| **Centroid exclusion from basins** | TRUE | Centroid coords span multiple basins (geographically ambiguous) |
| **Centroid inclusion in genotyping** | TRUE | Sequence quality independent of coordinate quality |

---

### Adjustable Parameters (via Config File)

| Parameter | Default | Range | What It Controls | When to Adjust |
|-----------|---------|-------|------------------|----------------|
| **`min_raw_length_abs`** | 200 bp | 50-1000 bp | Absolute minimum sequence length | Lower for degraded samples; raise for full-length barcodes |
| **`min_raw_length_frac_of_median`** | 0.7 | 0.5-1.0 | Minimum length as fraction of median | Lower for variable-length datasets |
| **`max_raw_N_fraction`** | 0.05 | 0.0-0.2 | Maximum ambiguous base content | Raise for low-quality datasets |
| **`require_valid_orf`** | TRUE | TRUE/FALSE | Require valid ORF for QC | FALSE if accepting pseudogenes/NUMTs |
| **`orf_min_coverage`** | 0.5 | 0.3-0.9 | Minimum ORF coverage | Raise for stricter validation |
| **`orf_max_internal_stops`** | 2 | 0-5 | Maximum internal stop codons | Lower for stricter validation |
| **`core_min_coverage`** | 0.8 | 0.5-1.0 | Fraction of sequences covering position | Lower for fragmentary datasets |
| **`mask_gap_threshold`** | 0.5 | 0.2-0.8 | Max gap fraction before masking column | Lower to mask more aggressively |
| **`core_min_length`** | 150 bp | 50-500 bp | Minimum core region length | Lower for short amplicons |
| **`min_singleton_distance`** | 0.005 | 0.0-0.02 | Singleton error filter threshold | 0.0 to disable; raise to filter more aggressively |
| **`max_singleton_distance`** | 0.05 | 0.02-0.15 | Singleton suspect flag threshold | Raise for divergent taxa |
| **`flag_suspect_haplotypes`** | FALSE | TRUE/FALSE | Enable suspect flagging | TRUE to flag distant singletons |
| **`majority_species_threshold`** | 0.70 | 0.5-1.0 | Majority vote threshold | Raise for conservative taxonomy |
| **`min_identity_pct`** | 98.5% | 95-100% | Sequence similarity species threshold | Adjust based on marker/taxon |
| **`min_query_cov_pct`** | 90% | 70-100% | Minimum query coverage for taxonomy | Lower for fragmentary sequences |
| **`top2_min_delta_pct`** | 0.5% | 0.1-2.0% | Margin for unambiguous call | Raise to be more conservative |
| **`min_consensus_length`** | 150 bp | 100-1000 bp | Minimum haplotype length for tree | Raise to 600 bp for full-length analyses |
| **`min_cluster_size`** | 3 | 1-10 | Minimum haplotype members for tree | Raise to avoid fragmentary consensuses |
| **`build_tree`** | FALSE | TRUE/FALSE | Enable phylogenetic tree building | TRUE for phylogenetic analyses |
| **`trim_alignment`** | TRUE | TRUE/FALSE | Trim alignment with trimAl | FALSE to keep all columns |
| **`trim_method`** | "gappyout" | "gappyout", "strict", "automated1" | trimAl method | "strict" for aggressive trimming |
| **`map_buffer_degrees`** | 20.0 | 5-90 | Map zoom buffer (degrees) | Adjust for regional vs. global focus |
| **`figure_dpi`** | 300 | 72-600 | Raster figure resolution | 600 for publication; 150 for web |
| **`figure_format`** | ["png", "pdf", "svg"] | List of formats | Output formats | Remove "svg" to reduce file size |

---

### User-Facing Flags (CLI Only)

| Flag | Effect | Use Case |
|------|--------|----------|
| **`--no-geo`** | Skip geographic analysis (all basins = "Unknown") | Non-marine organisms, missing GOaS data |
| **`--build-tree`** | Enable phylogenetic tree building | Phylogenetic studies |
| **`--no-report`** | Skip HTML report generation | Batch processing, CI/CD pipelines |
| **`--threads N`** | Use N CPU cores | Speed up parallelizable steps |
| **`--log-level`** | Set logging verbosity | DEBUG for troubleshooting |

---

## Output Interpretation Guide

### Key Output Files for Manuscript

| File | Purpose | Include in Manuscript? |
|------|---------|------------------------|
| **`{organism}_annotated.csv`** | Complete dataset with all annotations | Supplement (data availability) |
| **`haplotypes/{organism}_haplotypes.fasta`** | Representative haplotype sequences | Supplement (GenBank/Dryad) |
| **`haplotype_assignments/{organism}_diagnostics.csv`** | Assignment quality metrics | Supplement (quality control) |
| **`taxonomy/{organism}_haplotype_taxonomy.csv`** | Haplotype species assignments | Main text (Table 1) |
| **`species_analysis/species_summary.csv`** | Per-species diversity metrics | Main text (Table 2) |
| **`divergence_analysis/barcoding_gap.csv`** | Within vs. between-species divergence | Main text (Results) |
| **`divergence_analysis/divergence_matrix.csv`** | Pairwise haplotype distances | Supplement (distance matrix) |
| **`phylogenetic/{organism}_tree_relabeled.nwk`** | Publication-ready tree | Main text (Figure 1) |
| **`visualization/{organism}_distribution_map.pdf`** | Geographic distribution | Main text (Figure 2) |
| **`visualization/{organism}_distribution_bar.pdf`** | Basin abundance | Main text (Figure 3) |
| **`quality_control/README.md`** | QC summary and alerts | Supplement (methods validation) |
| **`{organism}_pipeline_parameters.json`** | All parameters used | Supplement (reproducibility) |

---

## Decision Points for Flow Chart

### Critical User Decisions

1. **Skip geographic analysis?**
   - YES: `--no-geo` (non-marine organisms, missing GOaS)
   - NO: Proceed with basin assignment

2. **Build phylogenetic tree?**
   - YES: `--build-tree` (phylogenetic questions)
   - NO: Skip to visualization (faster, genotyping-focused)

3. **Require valid ORF?**
   - YES: `require_valid_orf = TRUE` (exclude NUMTs, pseudogenes)
   - NO: `require_valid_orf = FALSE` (accept all sequences)

4. **Enable suspect flagging?**
   - YES: `flag_suspect_haplotypes = TRUE` (flag distant singletons)
   - NO: `flag_suspect_haplotypes = FALSE` (default, no flagging)

5. **Filter singleton errors?**
   - YES: Keep `min_singleton_distance = 0.005` (remove likely errors)
   - NO: Set `min_singleton_distance = 0.0` (keep all singletons)

---

## Processing Flow Chart Elements

### Suggested Flow Chart Structure

```
┌─────────────────────────────────────────────┐
│ PHASE 1: Data Loading & Quality Marking    │
├─────────────────────────────────────────────┤
│ Input: BOLD TSV                            │
│ ↓                                           │
│ Parse metadata (required columns)          │
│ ↓                                           │
│ Mark coordinate quality                    │
│ → has_centroid_coords (flag)              │
│ → has_missing_coords (flag)               │
│ → is_geographic_quality (calculated)      │
│ ↓                                           │
│ ┌─────────────────────────────────┐        │
│ │ Decision: GOaS available?       │        │
│ │ ├─ YES → Assign ocean basins   │        │
│ │ └─ NO → All basins = "Unknown" │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Output: All samples with quality markers   │
│         and ocean basin assignments        │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 2: Pre-processing Quality Control    │
├─────────────────────────────────────────────┤
│ Input: All samples with sequences          │
│ ↓                                           │
│ Orientation normalization                  │
│ → Try forward & reverse complement        │
│ → Select best ORF (highest coverage)      │
│ ↓                                           │
│ ORF validation                             │
│ → Calculate orf_coverage                  │
│ → Count orf_internal_stops                │
│ ↓                                           │
│ Dynamic QC filtering                       │
│ → Length: ≥200bp AND ≥70% median          │
│ → N content: ≤5%                           │
│ ┌─────────────────────────────────┐        │
│ │ Decision: require_valid_orf?    │        │
│ │ ├─ TRUE → Filter invalid ORFs  │        │
│ │ └─ FALSE → Keep all sequences  │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Output: QC-passed sequences                │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 3: Haplotype Discovery (ESV)         │
├─────────────────────────────────────────────┤
│ Input: QC-passed sequences                 │
│ ↓                                           │
│ Multiple sequence alignment (MAFFT)        │
│ ↓                                           │
│ Gap masking (columns >50% gaps)            │
│ ↓                                           │
│ Core region extraction                     │
│ → Positions covered by ≥80% sequences     │
│ → Minimum core length: 150bp              │
│ ↓                                           │
│ Identify exact sequence variants (ESVs)    │
│ → Group by 100% identity in core          │
│ → Each unique sequence = 1 haplotype      │
│ ↓                                           │
│ Singleton error filtering                  │
│ → Calculate nearest-neighbor distance     │
│ → REMOVE singletons ≤0.5% divergent       │
│ ↓                                           │
│ ┌─────────────────────────────────┐        │
│ │ Decision: Flag suspects?        │        │
│ │ ├─ YES → Flag singletons >5%   │        │
│ │ └─ NO → No flagging            │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Output: High-confidence haplotypes         │
│         (after error filtering)            │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 4: Haplotype Assignment              │
├─────────────────────────────────────────────┤
│ Input: Haplotype mapping from Phase 3     │
│ ↓                                           │
│ Direct ESV mapping                         │
│ → Each sample assigned to its ESV         │
│ → No distance-based re-assignment         │
│ → No ties (exact matches only)            │
│ ↓                                           │
│ Output: All QC-passed samples assigned     │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 5: Taxonomy Assignment               │
├─────────────────────────────────────────────┤
│ Input: Haplotypes with member lists        │
│ ↓                                           │
│ Metadata-based taxonomy (Majority vote)    │
│ → Aggregate species for each haplotype    │
│ → If ≥70% agree → assign species          │
│ → Else → genus or ambiguous               │
│ ↓                                           │
│ ┌─────────────────────────────────┐        │
│ │ Optional: Sequence-based        │        │
│ │ ├─ BLAST/VSEARCH vs. reference │        │
│ │ └─ Threshold: ≥98.5% identity  │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Conflict resolution                        │
│ → Integrate metadata + sequence           │
│ → Build labels: "{species} h{N}"          │
│ ↓                                           │
│ Output: Haplotypes with species labels     │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 6: Post-assignment QC                │
├─────────────────────────────────────────────┤
│ Input: Annotated samples with taxonomy    │
│ ↓                                           │
│ Contamination detection                    │
│ → Multi-species haplotypes                │
│ → Low-majority species (<70%)             │
│ → Per-sample purity scores                │
│ ↓                                           │
│ Quality reports & alerts                   │
│ ↓                                           │
│ Output: QC-flagged samples                 │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 6.5: Species-Level Aggregation       │
├─────────────────────────────────────────────┤
│ Input: Haplotypes with species assignments │
│ ↓                                           │
│ Aggregate by species                       │
│ → Filter to species-level (≥70% majority) │
│ ↓                                           │
│ Calculate diversity metrics                │
│ → Haplotype richness per species          │
│ → Haplotype evenness (Shannon)            │
│ ↓                                           │
│ Geographic summary by species              │
│ → Samples per species per basin           │
│ ↓                                           │
│ Generate species-faceted haplotype subsets │
│ → One file per species (if ≥2 haplotypes) │
│ ↓                                           │
│ Output: Species-level datasets             │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 7: Geographic Enhancement            │
├─────────────────────────────────────────────┤
│ Input: Annotated data with basins         │
│ ↓                                           │
│ Missing data assessment                    │
│ → Count centroid coordinates              │
│ → Count missing coordinates               │
│ ↓                                           │
│ Basin-level summaries                      │
│ → Samples, haplotypes, species per basin  │
│ ↓                                           │
│ Output: Enhanced geographic metadata       │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 8: Phylogenetic Analysis (Optional)  │
├─────────────────────────────────────────────┤
│ ┌─────────────────────────────────┐        │
│ │ Decision: build_tree?           │        │
│ │ ├─ NO → Skip to Phase 9        │        │
│ │ └─ YES → Continue              │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Filter haplotypes                          │
│ → Length ≥150bp (or user threshold)       │
│ → Cluster size ≥3 members                 │
│ ↓                                           │
│ Align filtered haplotypes (MAFFT)         │
│ ↓                                           │
│ ┌─────────────────────────────────┐        │
│ │ Decision: trim_alignment?       │        │
│ │ ├─ YES → Trim with trimAl      │        │
│ │ └─ NO → Use full alignment     │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Build tree (FastTree, GTR+Gamma)          │
│ ↓                                           │
│ ┌─────────────────────────────────┐        │
│ │ Decision: Outgroup specified?   │        │
│ │ ├─ YES → Root by outgroup      │        │
│ │ └─ NO → Midpoint root (if on)  │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Relabel tree with species names            │
│ ↓                                           │
│ ┌─────────────────────────────────┐        │
│ │ Optional: MSA visualization     │        │
│ │ → Phylogeny-ordered alignment  │        │
│ │ → Chunked plots (100bp chunks) │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Output: Tree files (.nwk), alignment,      │
│         MSA plots (if enabled)             │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 9: Divergence Analysis               │
├─────────────────────────────────────────────┤
│ Input: Haplotype sequences + taxonomy     │
│ ↓                                           │
│ Pairwise p-distance matrix                 │
│ → Gaps and Ns ignored                     │
│ → Formula: 1 - (matches/informative_sites)│
│ ↓                                           │
│ Barcoding gap analysis                     │
│ → Within-species distances                │
│ → Between-species distances               │
│ ↓                                           │
│ Species-level divergence                   │
│ → Mean distance between species pairs     │
│ ↓                                           │
│ Within-species divergence matrices         │
│ → Per-species haplotype divergence        │
│ ↓                                           │
│ Output: Distance matrices, barcoding gap   │
│         plots, heatmaps                    │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 10: Visualization                    │
├─────────────────────────────────────────────┤
│ Input: All analysis outputs                │
│ ↓                                           │
│ Haplotype-level plots                      │
│ → Distribution map (global)                │
│ → Ocean basin abundance (relative & total) │
│ ↓                                           │
│ Species-level plots                        │
│ → Species distribution map                 │
│ → Species basin abundance                  │
│ ↓                                           │
│ Species-faceted haplotype plots            │
│ → Per-species haplotype distribution maps │
│ → Per-species haplotype basin abundance   │
│ → Zoomed to species range (+20° buffer)   │
│ ↓                                           │
│ Output: PNG (300 DPI), PDF (vector), SVG,  │
│         JSON (interactive HTML)            │
└─────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────┐
│ PHASE 11: Report Generation (Optional)     │
├─────────────────────────────────────────────┤
│ ┌─────────────────────────────────┐        │
│ │ Decision: no_report?            │        │
│ │ ├─ YES → Skip report           │        │
│ │ └─ NO → Generate HTML          │        │
│ └─────────────────────────────────┘        │
│ ↓                                           │
│ Interactive HTML report                    │
│ → Summary dashboard                        │
│ → Pipeline parameters                      │
│ → Interactive Plotly.js visualizations    │
│ → Methods section (publication-ready)     │
│ → Data tables                              │
│ → Download links                           │
│ ↓                                           │
│ Output: {organism}_summary_report.html     │
└─────────────────────────────────────────────┘
```

---

## Recommended Flow Chart Visualization

For your manuscript, I recommend a **multi-level flow chart**:

### Level 1: High-Level Overview (Main Figure)
- 10 phases as major boxes
- Key decision points (diamonds): GOaS available? Build tree? Require ORF?
- Inputs/outputs as rectangles
- Data flow arrows

### Level 2: Detailed Sub-Process Diagrams (Supplemental)
- One diagram per complex phase (e.g., Phase 3 haplotype discovery)
- Shows internal steps and calculations

### Level 3: Parameter Tables (Supplemental)
- Use the tables provided above
- Clearly distinguish hardcoded vs. adjustable parameters

### Visual Elements to Include:

**Color Coding:**
- 🟦 Blue: Input data
- 🟩 Green: Processing steps
- 🟨 Yellow: Decision points
- 🟧 Orange: Quality control checks
- 🟥 Red: Filtering steps (data reduction)
- 🟪 Purple: Output files

**Annotations:**
- Parameter values in boxes (e.g., "Length ≥200bp")
- Sample counts at each major step (e.g., "1000 samples → 850 QC-passed")
- Optional steps marked with dashed borders

**Legend:**
- Adjustable parameters: Dotted border
- Hardcoded parameters: Solid border
- Optional steps: Dashed border
- Critical decision points: Diamond shape

---

## Key Takeaways for Manuscript

### 1. Non-Destructive Geographic QC
- Samples with centroid/missing coordinates are **retained for genotyping**
- Only excluded from ocean basin assignment (marked as "Unknown")
- Maximizes genotyping coverage while maintaining geographic accuracy

### 2. Two-Stage Singleton QC
- **Error filtering** (0.5% threshold): Removes likely sequencing/PCR errors
- **Suspect flagging** (5% threshold, optional): Flags potential contamination
- Retains biologically meaningful rare variants

### 3. ESV Approach vs. OTU Clustering
- True ESV: 100% identity in core region
- No arbitrary distance threshold
- Followed by error filtering to remove artifacts
- Preserves intraspecific variation

### 4. Comprehensive Taxonomy Integration
- Metadata-based (majority vote) + sequence-based (BLAST/VSEARCH)
- Conflict resolution with provenance tracking
- Transparent species assignment levels (species, genus, ambiguous)

### 5. Multi-Level Analysis
- Haplotype-level: Primary unit of analysis
- Species-level: Aggregated diversity metrics
- Species-faceted: Within-species haplotype structure
- Enables both broad and fine-scale insights

---

## Glossary for Flow Chart

**ESV (Exact Sequence Variant):** Unique sequence after error filtering; differs from any other sequence by ≥1 nucleotide in core region

**Core Region:** Alignment region covered by ≥80% of sequences; shared barcode fragment

**Haplotype:** Interchangeable with ESV in this context; unique mitochondrial COI lineage

**ORF (Open Reading Frame):** Stretch of DNA that can be translated without stop codons; validates COI authenticity

**Singleton:** Haplotype with only 1 member (n=1); higher risk of sequencing error

**Centroid Coordinates:** Country or region-level coordinates (not precise GPS); unsuitable for ocean basin assignment

**Geographic Quality:** Coordinates precise enough for basin assignment (not centroid, not missing, within valid range)

**Barcoding Gap:** Separation between within-species and between-species genetic distances

**p-Distance:** Proportion of differing nucleotides between sequences (uncorrected distance)

**Majority Threshold:** Minimum fraction of samples that must agree for species assignment (default: 70%)

---

## Citation for Methodology

When describing this workflow in your manuscript:

> "COI sequences were processed using BOLDGenotyper v1.0 (Smith, 2025), which implements a haplotype-first workflow based on exact sequence variants (ESVs). Briefly, sequences underwent orientation normalization, open reading frame validation, and dynamic quality control filtering. Haplotypes were identified by extracting a core shared barcode region and grouping sequences with 100% identity, followed by singleton error filtering to remove sequences ≤0.5% divergent from nearest neighbors. Samples were assigned to haplotypes via direct ESV mapping, and taxonomy was determined by majority vote (≥70% threshold). Geographic coordinates were quality-checked to exclude country-level centroids from ocean basin assignment while retaining these samples for genotyping. [Additional details: tree building, divergence analysis, visualization as applicable.] Pipeline parameters are provided in Supplemental Table X."

---

**End of Documentation**

*This document provides a complete description of the BOLDGenotyper workflow suitable for creating manuscript flow charts and supplemental methods sections. For questions or clarifications, please refer to the source code or contact the author.*
