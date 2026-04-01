# BOLDGenotyper: Comprehensive Technical Reference

**Internal Document — For Developer Reference Only**

*Author: Steph Smith | Version: as of March 2026*

---

## Overview

BOLDGenotyper is a Python command-line pipeline that automates the discovery, geographic
characterization, and population-genetics export of COI haplotypes from BOLD database downloads.
It takes a raw BOLD TSV file as input and produces: quality-filtered sequences, a multiple
sequence alignment, ESV-based haplotype definitions, haplotype-sample mapping, geographic region
assignments, species-level summaries, divergence analyses, pop-gen export files (Arlequin,
PopART/NEXUS, DnaSP), and an interactive HTML report.

The pipeline is **linear** (each phase depends on the previous) and is invoked via two entry
points:
- `boldgenotyper` — run the full pipeline on a taxon dataset
- `boldgenotyper-sweep` — run parameter sweep to find the optimal singleton distance threshold

---

## Architecture

```
BOLD TSV file
     │
     ▼
[1] DATA INGESTION         metadata.py
     │
     ▼
[2] QUALITY CONTROL        quality_control.py
     │
     ▼
[3] ALIGNMENT              dereplication.py → MAFFT
     │
     ▼
[4] HAPLOTYPE DISCOVERY    dereplication.py (ESV)
     │
     ▼
[5] HAPLOTYPE ASSIGNMENT   haplotype_assignment.py → edlib
     │
     ▼
[6] GEOGRAPHIC ANALYSIS    geographic.py → GeoPandas
     │
     ▼
[7] TAXONOMY + DIVERGENCE  utils.py, divergence_analysis.py
     │
     ▼
[8] SPECIES ANALYSIS       species_analysis.py
     │
     ▼
[9] PHYLOGENETICS           phylogenetics.py → trimAl → FastTree
     │
     ▼
[10] POP-GEN EXPORT         popgen_export.py
     │
     ▼
[11] METADATA ANALYSIS      metadata_analysis.py
     │
     ▼
[12] VISUALIZATION          visualization.py → cartopy, matplotlib
     │
     ▼
[13] HTML REPORT            reports.py → Plotly.js
```

Configuration for every parameter lives in `boldgenotyper/config.py`.
Support utilities (FASTA I/O, ORF validation, external tool checks, taxonomy assignment)
live in `boldgenotyper/utils.py`.

---

## Module-by-Module Reference

---

### 1. Data Ingestion — `metadata.py`

#### Purpose
Parse the raw BOLD TSV download and normalize it into a clean pandas DataFrame ready for
downstream QC.

#### Key Functions

**`parse_bold_tsv(file_path, ...)`**

Opens the file with UTF-8 encoding; if that fails (common with BOLD exports containing
accented characters in collector names, institution names, etc.), falls back to latin-1
encoding. Reads the file using `pandas.read_csv` with tab separator and `low_memory=False`
to prevent dtype inference errors on large files.

After loading:
- Skips malformed rows (rows where essential columns are empty)
- Deduplicates on `processid`, keeping the first occurrence when duplicates appear
  (BOLD occasionally exports the same specimen twice in large downloads)
- Validates that the two mandatory columns are present: `processid` (unique specimen
  identifier) and `nuc` (COI sequence string)
- Logs N rows read, N after deduplication

**Why parse manually rather than using a BOLD API client?**
BOLD's web portal allows bulk downloads as TSV files and users have already downloaded their
data. Parsing TSV directly means no API key, no rate limits, and offline usability. The
pipeline is explicitly designed for the "download then analyze" workflow.

**`extract_coordinates(df, ...)`**

BOLD stores coordinates in a single `coord` column in `[lat, lon]` bracket format. This
function parses that string into separate numeric `lat` and `lon` columns.

Parsing logic:
1. Strip brackets and split on commas/whitespace
2. Attempt float conversion for each component
3. Validate ranges: latitude in −90 to +90, longitude in −180 to +180
4. Set lat/lon to NaN for rows failing validation

**`infer_coord_order(sample_coords)`**

BOLD data sometimes has lat/lon reversed. This heuristic detects the issue by checking
whether the first value falls within valid latitude range (−90 to +90) and the second within
longitude range (−180 to +180). If the pattern is reversed in a statistically dominant
fraction of records, the columns are swapped. This is a sanity check, not a universal fix
— marine datasets near ±180° meridian can trigger false positives, which is why the function
logs a warning rather than silently correcting.

**Centroid detection and exclusion**

Samples with country-level centroids (where the depositor provided only a country location,
so BOLD assigns the geographic centroid of that country as the coordinate) are identified via
the `coord_source` metadata field. These samples receive coordinates but the coordinates
represent the national centroid, not the actual collection point. For geographic analysis
(point-in-polygon ocean/basin assignment), centroid samples are excluded because placing a
shark sample at the centroid of Australia would assign it to a terrestrial region. They are
retained in all other analyses.

---

### 2. Quality Control — `quality_control.py`

#### Purpose
Filter sequences that are too short, too ambiguous, or incorrectly oriented before alignment.
All decisions are documented in the output DataFrame so users can see exactly why each
sequence passed or failed.

#### Dynamic QC Philosophy

A fixed minimum length threshold fails in heterogeneous datasets. Consider a family-level
dataset containing both ~600 bp and ~200 bp COI fragments for different genera — a fixed
300 bp cutoff would discard legitimate short sequences. BOLDGenotyper uses a **three-tier
dynamic system**:

1. **Absolute minimum (200 bp)**: Any sequence shorter than 200 bp is insufficient for
   meaningful barcode analysis regardless of the dataset. This is a hard biological floor.

2. **Adaptive length filter (≥70% of dataset median)**: Computes the median sequence length
   across all sequences after removing leading/trailing gaps, then rejects sequences below
   70% of that value. For a dataset with median 600 bp, this threshold is 420 bp. For a
   dataset with median 250 bp, it becomes 175 bp (but the absolute floor of 200 bp still
   applies). This adaptive approach handles multi-primer datasets where different amplicons
   produce different lengths legitimately.

   **Why 70%?** In our validation work, sequences at 70% of the dataset median are typically
   truncated amplicons (failed at one end) or chimeric sequences, not shorter legitimate
   barcodes. At 80%, we discard too many short-but-valid sequences from rare taxa. At 60%,
   we retain clear truncations. 70% was validated empirically across the five benchmark taxa.

3. **ORF validity filter**: Translates the sequence using vertebrate mitochondrial code (NCBI
   genetic code 2) and rejects sequences with >2 internal stop codons. The rationale for
   allowing up to 2 internal stops: BOLD sequences occasionally contain brief ambiguous
   stretches (NNNN blocks) that translate as stops but are not genuine stop codons; setting
   the threshold at 2 catches sequences with actual pseudogene contamination or nuclear
   mitochondrial DNA inserts (NUMTs) while tolerating minor ambiguity artifacts.

**`apply_orientation_normalization(sequences_dict, ...)`**

COI sequences in BOLD are deposited by many different labs and can be in either forward
or reverse complement orientation. MAFFT will still align them, but ESV dereplication
(exact sequence comparison) will completely fail to identify identical haplotypes that are
deposited in opposite orientations.

Process:
1. For each sequence, attempt ORF validation in the forward orientation
2. Attempt ORF validation in the reverse complement
3. Select the orientation producing the better ORF quality score
4. If reverse complement is better, flag the sequence as corrected and store it corrected
5. Return a dict of corrected sequences plus a DataFrame of ORF stats per sequence

**Why we normalize before alignment (not after)**: MAFFT's `--auto` mode can usually
cope with mixed-orientation input, but downstream core-region extraction and ESV comparison
require consistent orientation. Correcting before alignment ensures the aligned positions
are biologically comparable.

**ORF quality scoring (`check_orf_quality` in `utils.py`)**: Checks for the presence of a
valid start codon (ATG) in the sequence, counts internal stop codons (TAA/TAG/TGA in the
correct reading frame), and computes a coverage fraction (length of translated region /
expected COI length). Returns a composite score that is used to select the better orientation.

**`apply_dynamic_qc_filters(df, ...)`**

Applies all three filters simultaneously and annotates each row with:
- `qc_pass`: Boolean True/False
- `qc_fail_reason`: Human-readable string explaining failure ('too_short', 'below_median',
  'invalid_orf', or combinations)

Preserving fail reasons in the output allows users to audit QC decisions. This is
particularly important for publication — a reviewer might ask why certain sequences were
excluded.

**Depositor uncertainty detection**

The `UNCERTAINTY_KEYWORDS` dictionary contains terms like "cryptic", "juvenile", "cf.",
"sp.", "nr.", "aff.", "mixed", "larva", "unknown" that indicate taxonomic uncertainty in
BOLD depositor annotations. Sequences with these keywords in their species name are flagged
(not removed) to allow users to decide whether to include uncertain identifications.

---

### 3. Multiple Sequence Alignment — `dereplication.py`

#### Purpose
Align all quality-passed sequences to enable core-region extraction and ESV identification.

#### Why MAFFT?

MAFFT is the de facto standard for large-scale progressive multiple sequence alignment in
molecular ecology. Key reasons for choosing it over alternatives:

- **Speed**: MAFFT's progressive alignment algorithms (FFT-NS-1, FFT-NS-2, G-INS-1) scale
  to thousands of sequences in minutes. Clustal Omega is similar in speed; T-Coffee is much
  slower.
- **Accuracy**: Benchmarked repeatedly (Thompson et al. 2011 comparative studies) as
  top-performing for divergent sequences, which is critical for family- or order-level
  datasets.
- **Auto mode**: `--auto` selects the appropriate algorithm based on input size and
  divergence level automatically — FFT-NS-2 for large/divergent sets, G-INS-i for small/
  conserved sets. This means the same command works correctly across our five benchmark taxa
  that range from 755 to 9,244 sequences.
- **Conda/bioconda availability**: Easy installation in any environment.
- **BioPython independence**: We call MAFFT as an external subprocess, which means
  BOLDGenotyper doesn't require BioPython for alignment (only `divergence_analysis.py`
  uses BioPython for sequence I/O in the heatmap generation step).

**`run_mafft_alignment(sequences_dict, output_path, mafft_algorithm='auto', threads=1)`**

Constructs the MAFFT command with `--auto --quiet --thread N`. The `--quiet` flag suppresses
MAFFT's verbose progress output (which is appropriate for interactive use but clutters logs
in automated pipelines). Raises `AlignmentError` if MAFFT is not found in PATH, rather than
silently skipping alignment.

**Parallelization**: MAFFT uses `--thread N` where N is the user-specified thread count
(default 4). MAFFT's threading is primarily for pairwise distance computation in progressive
alignment; it doesn't parallelize perfectly but provides meaningful speedup on large datasets.

---

### 4. Core Region Extraction and ESV Haplotype Discovery — `dereplication.py`

This is the algorithmic heart of BOLDGenotyper.

#### What is an ESV?

ESV stands for **Exact Sequence Variant** (Porter & Hajibabaei 2020, *PeerJ*). An ESV is a
unique, exact sequence — no clustering at a distance threshold, no consensus across multiple
reads. For amplicon data with a single copy target per specimen (as COI barcoding is), ESVs
are the finest possible resolution: each distinct sequence is its own haplotype.

**Contrast with OTU clustering**: Traditional OTU approaches cluster sequences within X%
divergence (e.g., 97% identity = 3% divergence threshold). This collapses genuine haplotype
diversity. For example, two haplotypes differing at 1 bp out of 600 bp (0.17% divergence)
would be merged into a single OTU at 97% clustering, losing real population genetic signal.
ESVs preserve every genuine variant.

**Why ESVs are appropriate for COI haplotyping**:
1. COI is single-copy mitochondrial — no paralogs, no multi-allelic signals per specimen
2. For barcode-length COI (typically 658 bp), even 1 bp differences are biologically
   meaningful for population genetics (different haplotypes in the traditional network sense)
3. BOLD sequences represent processed, curated specimens (not raw metabarcoding reads), so
   sequencing error rates are low enough that most variants are genuine

#### Core Region Identification

**Why extract a core region at all?**

Raw COI sequences in BOLD have variable start and end positions due to:
- Different primer pairs amplifying slightly different lengths
- Sequencing quality deteriorating at sequence ends (common in Sanger sequencing)
- Some labs submitting only the "good quality" interior of a longer read

If we compare full-length sequences directly, two sequences that are biologically identical
(same haplotype) but differ only in which end positions were sequenced would appear as
different ESVs. This is a technical artifact, not biological.

The core region is the set of alignment positions covered by ≥80% of sequences after
gap-column masking. Any position covered by fewer than 80% of sequences is edge-effect
noise from primer and trimming variation.

**`identify_core_region(alignment_dict, min_coverage=0.8, gap_threshold=0.5, min_length=150)`**

Algorithm:
1. For each alignment column, count the fraction of sequences with a non-gap character
2. A column is "gap-heavy" if >50% of sequences have a gap there — these are masked as
   uninformative before coverage calculation
3. A column is in the core region if ≥80% of sequences (after gap masking) have a
   non-gap character there
4. The core region is the contiguous stretch of passing columns
5. If the core region is shorter than 150 bp, the function raises an error — 150 bp is
   the minimum length for meaningful barcode comparison

**Why 80% coverage threshold?** At 90%, we lose real data from taxa where a minority of
sequences use slightly different primers. At 70%, we include positions with too much missing
data that would make comparisons unreliable. 80% is validated across all five benchmark taxa.

**Why 50% gap-column masking threshold?** Alignment columns where more than half of sequences
have gaps typically represent insertion events in a small subset of sequences or alignment
artifacts at boundaries. These columns add noise to ESV comparison without adding signal.

**`identify_unique_haplotypes(alignment_dict, core_start, core_end, ...)`**

Algorithm:
1. For each sequence, extract the core-region substring (ungapped)
2. Group sequences by their core-region substring — all sequences with identical core
   sequences form one haplotype group
3. Sort groups by member count, descending (most abundant haplotype gets ID `h1`)
4. Name haplotypes as `h{rank}_n{member_count}`, e.g., `h1_n179` (haplotype 1, 179 members)
5. For each group, select one representative sequence (the first one, arbitrarily, since
   they are all identical in the core) to be the "consensus" FASTA sequence for downstream use
6. Record the consensus sequence, member list, member count, and rank for each haplotype

**Haplotype naming format `h{N}_n{count}`**: Rank N is assigned at haplotype discovery time
and remains stable for a given run. The `n{count}` suffix in the ID encodes the member count
directly in the filename/column so users can immediately see abundance from the ID without
joining to a separate table. This is a deliberate usability choice.

#### Pairwise Distance Calculation for Singleton Filtering

**`_compute_distance(seq1, seq2)`**

Computes p-distance = `1 − (matches / valid_sites)` where valid sites are positions where
BOTH sequences have unambiguous ACGT bases (not gaps, not N). This is the correct
uncorrected p-distance formula for COI:
- Gaps are excluded from the denominator (gaps are alignment artifacts)
- N bases are excluded from the denominator (ambiguous calls don't contribute information)
- Both sequences must have valid bases at a position for it to count

This is the standard uncorrected p-distance used in barcoding studies (e.g., "10× rule"
papers define intraspecific divergence in these terms).

**`calculate_pairwise_distances(haplotypes_dict)`**

Computes all-pairs distances in O(n²) time. Returns a condensed distance vector in the
format used by `scipy.spatial.distance.squareform`. For datasets with hundreds of haplotypes
(e.g., Carcharhiniformes with 487), this is still fast (<1 second) because:
1. The number of haplotypes is much smaller than the number of input sequences
2. Core-region sequences are already extracted so comparisons are on short, ungapped strings
3. Pure NumPy operations are used for the inner loop

#### Singleton Error Filtering

**The problem**: Sanger COI sequencing has error rates of approximately 0.1–0.5% per
position. A single sequencing error creates a false unique sequence (a false singleton
haplotype) that is nearly identical to its correct counterpart. If unfiltered, a dataset
with 1,000 samples might yield 200+ spurious singleton haplotypes from sequencing errors.

**`filter_error_singletons(haplotypes, pairwise_distances, min_singleton_distance=0.005)`**

Logic:
1. Identify singleton haplotypes (member count = 1)
2. For each singleton, find its nearest neighbor (the haplotype with smallest p-distance)
3. If nearest_neighbor_distance ≤ `min_singleton_distance`, remove the singleton — it is
   likely a sequencing error variant of the neighbor haplotype
4. If nearest_neighbor_distance > `min_singleton_distance`, retain the singleton — it
   represents genuine divergence

**Why 0.5% (0.005) as the default threshold?**
Empirical validation across 5 taxa showed that ~85% of singletons in real BOLD datasets
fall within 0.3% divergence of their nearest neighbor (1–2 bp differences in 600 bp
sequences). These are essentially certain sequencing errors. The 0.5% threshold provides
a small buffer above the 0.3% practical maximum for Sanger error while remaining well
below the 1.5% typical minimum intraspecific divergence for most organisms. This was the
lowest threshold that captured essentially all errors without reaching into genuine
intraspecific variation.

The `boldgenotyper-sweep` tool allows users to empirically determine the appropriate
threshold for their specific taxon — the default 0.5% is a reasonable starting point but
taxa with low genetic diversity (e.g., recently-expanded populations) may benefit from a
lower threshold.

#### Suspect Haplotype Flagging

**`flag_suspect_haplotypes(haplotypes, pairwise_distances, max_singleton_distance=0.05)`**

A second, separate pass that flags (but does NOT remove) singleton haplotypes that are
highly divergent from all other haplotypes:
1. Singletons where nearest_neighbor_distance > 5% are flagged as "outlier-high"
2. Singletons at or above the 95th percentile of all pairwise distances AND >10% absolute
   divergence are flagged as "outlier-extreme"
3. Singleton sequences failing ORF validation are flagged as "orf_fail"

**Why flag rather than remove?** High-divergence singletons could be:
- Contamination (most likely)
- Genuinely divergent cryptic lineages (biologically interesting!)
- NUMTs (nuclear mitochondrial DNA pseudogenes)

BOLDGenotyper cannot reliably distinguish these cases from sequence data alone. Flagging
puts the decision in the user's hands with full transparency. The HTML report shows flagged
haplotypes prominently so users can investigate.

By default, `filter_suspect_haplotypes=False` and `flag_suspect_haplotypes=False` — both
features are opt-in via CLI flags or config, because most users want to see all haplotypes
in the initial run and decide for themselves.

#### Zero-Divergence Merging

**`merge_zero_divergence_haplotypes(haplotypes, pairwise_distances)`**

After core-region extraction, two haplotypes can have 0 pairwise divergence (identical
sequences) because their variable positions lie entirely outside the core region. This
happens when:
- Two primer pairs amplify slightly different regions
- Some sequences are longer than others at both ends

These are the SAME haplotype from a biological standpoint — the core region is identical.
Union-find (disjoint sets) algorithm merges them transitively: if A and B are zero-divergent,
and B and C are zero-divergent, then A, B, and C are all merged into one haplotype
(represented by the largest-member group).

**Why Union-Find?** Transitive merging requires tracking connected components. Union-find
achieves this in near-O(n) time with path compression.

---

### 5. Haplotype Assignment — `haplotype_assignment.py`

#### Purpose

Assign each individual input sequence (each BOLD specimen) to the closest haplotype. This
maps the full dataset (thousands of sequences) onto the small set of haplotypes (tens to
hundreds).

#### Why assignment is needed

At the ESV discovery stage, sequences were grouped by identical core-region content. So in
principle, assignment is trivial: look up which haplotype group contains your processid.

In practice, assignment is a separate pass because:
1. Some sequences were excluded from haplotype discovery by QC filters but can still be
   assigned if they are sufficiently similar to a retained haplotype
2. The FASTA used for alignment may contain slight formatting differences from the metadata
   TSV, requiring processid extraction from headers
3. We want identity scores and tie detection for every sample, which requires explicit
   distance calculation

#### Identity Metric: Target-Based Identity

**`calculate_identity_with_cigar(seq1, seq2)`**

Two identity metrics are computed:

1. **target_identity** = matches / len(target_sequence)
   This is the primary assignment metric. It asks: "What fraction of the consensus (haplotype)
   sequence is represented in the sample?" This is robust to noisy 5' and 3' ends in the
   sample — if a sample has 20 extra bases of low-quality sequence at one end but perfect
   match to the consensus over the full consensus length, it gets identity = 1.0.

2. **classic_identity** = 1 − (edit_distance / max_length)
   The traditional metric, computed and stored as a diagnostic column but NOT used for
   assignment decisions.

**Why target-based?** COI barcode sequences from BOLD vary in the precise start and end
positions. If sample A is 620 bp and consensus is 600 bp, but they match perfectly over the
600 bp consensus region, the biologically correct answer is 100% identical. Classic identity
penalizes to 97% (20 extra bp out of 620). Target-based correctly gives 100%.

#### Edit Distance vs. p-Distance for Assignment

Assignment uses edit distance (via edlib or fallback Levenshtein), not p-distance.

**Why?** p-distance requires an alignment — it counts matches and mismatches at homologous
positions. For assignment, we're comparing raw sequences of slightly different lengths. Edit
distance handles this correctly without requiring an additional alignment step per comparison.
This is also much faster: edlib (C-based) computes edit distances in microseconds.

#### edlib Dependency

**`check_edlib_available()`** and graceful fallback to pure Python `levenshtein_distance()`

`edlib` is a C extension for fast global edit distance (Levenshtein distance). On a dataset
with 1,000 sequences and 50 haplotypes, naive Python Levenshtein requires 50,000 comparisons
of ~600-character strings. edlib's C implementation is ~100× faster than pure Python for
this use case.

**The fallback**: Pure Python two-row DP Levenshtein with O(n*m) time and O(min(n,m)) space
is implemented in `levenshtein_distance()`. This ensures the pipeline works even when edlib
is not installed, just more slowly.

**Why not use BioPython pairwise alignment for this?** BioPython's pairwise alignment is
designed for biological accuracy with gap penalties and substitution matrices. For assignment,
we want raw edit distance as a similarity proxy — no gap penalty tuning needed. edlib's global
NW alignment is cleaner for this purpose.

#### Tie Detection

**Tie condition**: `best_identity ≥ 0.95 AND (best_identity − runner_up_identity) < 0.003`

A tie is flagged when the best match is highly similar (≥95%) AND the second-best match is
within 0.3% identity. This captures cases where 1–2 bp differences separate the best and
second-best haplotype candidates, making the assignment genuinely ambiguous.

**Why 0.3% margin?** At 600 bp, 0.3% = ~1.8 bp. Assignment ambiguity between haplotypes
that differ by only 1–2 bp is real — a sequencing error in the sample could tip the
assignment either way. The 0.3% margin captures this zone of genuine uncertainty without
flagging the majority of unambiguous assignments (where the gap is typically >>5%).

**Why 95% minimum for tie detection?** Below 95% identity, the sequences are too divergent
for assignment uncertainty to matter — both candidates are bad matches. Tie detection is only
meaningful in the zone of high similarity.

#### Multiprocessing

Assignment uses `multiprocessing.Pool` with `partial` to parallelize across samples. Each
worker receives a processid + sequence pair and independently computes distances to all
haplotypes. Results are collected and aggregated in the main process. Thread count is
controlled by the `--threads` CLI flag.

---

### 6. Geographic Analysis — `geographic.py`

#### Purpose

Assign each sample (that has valid coordinates) to a geographic region (ocean basin,
freshwater basin, ecoregion, etc.) via point-in-polygon spatial join.

#### Core Concept: Point-in-Polygon

Every sample with valid latitude/longitude is represented as a Shapely Point geometry.
The polygon shapefile (GOaS, BasinATLAS, Ecoregions2017, or user-supplied) provides
named regions as polygon features. The spatial join (`gpd.sjoin`) assigns each point to
the polygon it falls within.

**`load_goas_data(shp_path)`**

Loads the Global Oceans and Seas (GOaS) v1 shapefile from Marine Regions. Key steps:
1. Read with GeoPandas
2. Check and reproject to WGS84 (EPSG:4326) if the shapefile is in a different CRS
3. Standardize the region name column to `basin_name`

**Why GOaS as the marine default?** GOaS is the standard geographic framework for marine
biogeography — it is used in publications like the PLOS ONE Sphyrna lewini companion paper
and is maintained by VLIZ (Flanders Marine Institute). It provides 10 named ocean regions
(North Atlantic, South Atlantic, Indian Ocean, North Pacific, South Pacific, Arctic, etc.)
that align with the biogeographic regions commonly discussed in molecular ecology.

**`create_points_geodataframe(df, lat_col='lat', lon_col='lon')`**

Creates Shapely Point objects from numeric lat/lon columns. Validates:
- Latitude in range [−90, +90]
- Longitude in range [−180, +180]
- Sets geometry to None for rows with NaN or out-of-range coordinates

Returns a GeoDataFrame in EPSG:4326.

**`assign_ocean_basins(points_gdf, goas_gdf, ...)`** and
**`assign_regions_from_shapefile(points_gdf, shapefile_path, region_field, ...)`**

These functions perform the spatial join. Steps:
1. Filter to samples with valid (non-null) geometries
2. Use `gpd.sjoin(left_df=points, right_df=polygons, how='left', predicate='within')`
3. Join is 'left' so all samples are retained (samples outside all polygons get NaN region)
4. A secondary fallback uses the BOLD `country/ocean` metadata column for samples outside
   all polygons but with a country name (e.g., open ocean samples not covered by GOaS)

**Why GeoPandas?** GeoPandas wraps Shapely's spatial operations with a familiar pandas API
and handles CRS reprojection transparently. The alternative would be PostGIS (requires a
database server) or manual polygon containment tests (much slower, requires custom code).
GeoPandas + Shapely is the standard Python geospatial stack and handles shapefiles directly
without a database.

**Why WGS84 / EPSG:4326?** BOLD coordinates are in decimal degrees WGS84. All public
shapefiles (GOaS, BasinATLAS, Ecoregions2017) are available in WGS84. No projection
mathematics needed, minimizing coordinate transformation error.

#### Freshwater Basins: BasinATLAS v10 lev03

For freshwater datasets (Salmonidae benchmark), we use BasinATLAS v10 at level 03:
- Level 03 = major global drainage basins (292 features) identified by Pfafstetter code
  (`PFAF_ID` field)
- This level captures biogeographically meaningful units (e.g., North American Pacific
  vs. Atlantic drainages) without excessive fragmentation (level 07–12 = thousands of
  small sub-basins, overkill for species-level analysis)

#### Terrestrial Ecoregions: Ecoregions2017

For terrestrial datasets (Pieridae benchmark), we use WWF Terrestrial Ecoregions of the
World (Olson et al. 2001, updated 2017). We use the `REALM` field (8 biogeographic realms:
Palearctic, Nearctic, Neotropical, Afrotropic, Indo-Malay, Australasia, Oceania, Antarctic)
rather than the `ECO_NAME` field (847 individual ecoregions) because:
- Realm is the appropriate scale for population-level biogeography of butterflies
- 847 ecoregion names produce unreadable figures; 8 realms is interpretable
- The companion Pieridae analysis focuses on inter-realm structure

---

### 7. Geographic Coverage Assessment — `geographic_enhancement.py`

#### Purpose

Quantify and report the quality of geographic data coverage per haplotype, to ensure users
understand the limitations of the geographic analysis.

**`assess_geographic_coverage(df, group_col='haplotype_id', ...)`**

For each haplotype, calculates:
- N total samples
- N samples with lat/lon coordinates
- N samples assigned to a region (excluding "Unknown")
- N samples with country data

Assigns a representativeness tier:
- **excellent**: >50% of samples have region assignment
- **good**: 25–50%
- **moderate**: 10–25%, or ≥50% have coordinates even without region assignment
- **poor**: some coordinates but <25% region coverage
- **very_poor**: almost no coordinate data

This assessment appears in the HTML report as a data quality indicator, helping users
distinguish well-supported geographic statements from speculative ones.

**Why is this a separate module?** Geographic coverage assessment was added after initial
development as a transparency feature requested during manuscript review. Separating it from
the core `geographic.py` keeps the spatial join logic clean and makes the coverage reporting
independently testable.

---

### 8. Divergence Analysis — `divergence_analysis.py`

#### Purpose

Calculate pairwise uncorrected p-distances between all haplotype consensus sequences
to assess: (1) within-species vs. between-species divergence, (2) barcoding gap presence,
(3) produce a divergence heatmap.

#### Distance Metric: Uncorrected p-Distance

`calculate_pairwise_distance(seq1, seq2)` = `1 − (matches / informative_positions)`

Informative positions exclude gaps (`-`) and ambiguous bases (`N`). Both sequences must
have valid ACGT at a position for it to count.

**Why uncorrected p-distance (not Kimura 2-parameter, not JC)?**
- **Consistency with barcoding literature**: The vast majority of COI barcoding papers
  (Hebert et al., Ward et al., Ratnasingham & Hebert) report uncorrected p-distances
  ("percent sequence divergence"). Using p-distance allows direct comparison with published
  divergence values.
- **Transparency**: p-distance is directly interpretable as "fraction of positions that differ"
  without correction model assumptions. For divergences typical in COI barcoding (<20%),
  p-distance and K2P are nearly numerically identical.
- **Appropriate scale**: Correction models (K2P, Tamura-Nei) matter most for deep
  divergences where multiple hits saturate individual sites. Within species and at family
  level, COI divergences are typically <15%, well below the saturation zone.

#### Alignment for Divergence Calculation

`align_consensus_sequences()` re-aligns the haplotype consensus FASTA using MAFFT `--auto`
before calculating pairwise distances. This creates a separate alignment from the main
pipeline alignment (which included all input sequences). Re-aligning consensus sequences
only (typically tens to a few hundred sequences vs. thousands of input sequences) is fast
and ensures the divergence matrix is calculated from a clean, haplotype-only alignment.

#### Barcoding Gap Analysis

`analyze_barcoding_gap(divergence_matrix, taxonomy_df)`

Collects all within-species pairwise divergences (pairs of haplotypes with the same species
assignment) and all between-species pairwise divergences (pairs with different species).
For each species:
- Reports mean/max within-species divergence
- Reports mean/min between-species divergence
- Tests whether a barcoding gap exists: `min(between) > max(within)`
- Reports gap size if present

The barcoding gap is a fundamental DNA barcoding concept (Hebert et al. 2003, Meyer &
Paulay 2005): species are diagnosable if the smallest interspecific divergence exceeds the
largest intraspecific divergence. BOLDGenotyper reports this empirically without making any
species delimitation claims.

#### Divergence Heatmap

`create_divergence_heatmap()` produces a clustered heatmap using seaborn's `clustermap`
with complete-linkage hierarchical clustering and adaptive figure sizing based on haplotype
count. The `mako_r` color palette is perceptually uniform and colorblind-friendly.

Cell annotations show p-distance values; diagonal cells (self-comparison, distance = 0)
show "×" to distinguish from genuine near-zero divergences.

---

### 9. Taxonomy Assignment — `utils.py`

#### `assign_consensus_taxonomy(df, group_col, species_col, majority_threshold=0.5)`

Assigns a consensus species name to each haplotype group using a majority-rule vote:
1. Count species-name strings among all samples in each haplotype group
2. If one species name appears in >50% of samples (configurable threshold), assign it
3. If no species reaches the threshold, fall back to genus-level: assign the most common
   genus name
4. If even genus cannot be determined, assign "Unknown"

The `TaxonomyConfig.majority_species_threshold` defaults to 0.70 (70%) for final taxonomy
tables, but the underlying function uses 0.50 for internal intermediate tables. This creates
a stricter gate for the final published taxonomy assignments.

**Why majority-rule rather than strict consensus?**
- BOLD species names for haplotypes are depositor-provided and often inconsistent (e.g.,
  "Sphyrna lewini" vs. "Sphyrna lewini (Griffith & Smith, 1834)")
- A haplotype shared across species boundaries (e.g., recently-diverged sister taxa) would
  have a 50/50 species split — majority rule correctly calls this ambiguous
- 70% threshold means at least 70% of members agree before a species name is assigned;
  this filters out contamination and misidentifications at the haplotype level

**`haplotype_sp` label format**: The taxonomy table creates a `haplotype_sp` column
formatted as `{species} {haplotype_id}`, e.g., "Sphyrna lewini h1_n179". This label is
used as the tip label in phylogenetic trees, making the tree immediately interpretable
without requiring a separate legend.

---

### 10. Species-Level Analysis — `species_analysis.py`

#### Purpose

Aggregate haplotype-level data to the species level for biologically meaningful summaries.

**`aggregate_samples_by_species(annotated_metadata, species_composition, min_confidence=0.7)`**

Joins sample metadata (with haplotype assignments) to the haplotype species composition
table. Filters to confident species assignments (primary_species_pct ≥ 0.7) to exclude
haplotypes with ambiguous taxonomy. Returns:
- `species_assignments`: sample-level table with species added
- `species_summary`: one row per species with N samples, N haplotypes, haplotype list

**`calculate_species_diversity(species_assignments)`**

Calculates per-species:
- **Haplotype diversity (h)**: `(n/(n−1)) × (1 − Σpᵢ²)` where pᵢ is frequency of each
  haplotype. This is the standard Nei (1987) haplotype diversity estimator. The
  correction factor `n/(n−1)` accounts for small sample sizes. This metric = probability
  that two randomly chosen samples from the species have different haplotypes.
- **Simpson's index**: `1 − Σpᵢ²` (same formula without the sampling correction).
  Included as an alternative for comparison.

**Why Nei's h?** This is the standard haplotype diversity metric in molecular ecology and
population genetics (DnaSP, Arlequin, PopART all compute it). Reporting it in BOLDGenotyper
makes the output directly comparable to any study using these standard tools.

**`generate_within_species_divergence_matrices(divergence_matrix, ...)`**

Subsets the full haplotype pairwise divergence matrix to haplotypes within each species,
enabling per-species divergence visualization. This is particularly useful for multi-species
family-level datasets (e.g., Salmonidae with 260 haplotypes across dozens of species) where
the full matrix is too large to interpret.

---

### 11. Phylogenetics — `phylogenetics.py`

#### Purpose

Build a maximum likelihood phylogenetic tree from haplotype consensus sequences to:
1. Verify ESV haplotypes cluster correctly by species/taxonomy
2. Reveal phylogeographic structure within species
3. Provide a tree for manual interpretation and figures

#### Sequence Filtering Before Tree Building

`filter_consensus_sequences(consensus_fasta, min_length=600, min_cluster_size=5)`

Two filters applied before phylogenetics:
1. **Minimum ungapped length ≥ 600 bp**: Short consensus sequences produce unreliable
   branch lengths. 600 bp = the standard "full-length barcode" length for COI-5P. Using
   the full length provides maximum phylogenetic resolution.
2. **Minimum cluster size ≥ 5**: Haplotypes with <5 members are excluded from tree
   building. Singletons and very rare haplotypes (after singleton filtering) contribute
   minimal phylogenetic signal while adding computation time. More importantly, very rare
   haplotypes are more likely to represent residual sequencing errors or contamination.

**Critical distinction**: These filters affect ONLY phylogenetics. All haplotypes
(including small clusters) are retained for assignment, geographic analysis, and pop-gen
export. The tree is exploratory/visualization; the haplotype system is the analysis result.

#### trimAl — Alignment Trimming

`run_trimal(alignment_fasta, output_fasta, method='gappyout')`

trimAl removes gap-heavy columns from the alignment before tree building.

**Why trimAl?** Phylogenetic tree building is sensitive to alignment quality. Columns
where many sequences have gaps (typically from indels in a minority of sequences or from
5'/3' end variation) produce spurious branch length estimates. trimAl's `gappyout` method
removes columns that introduce gaps without significant information content.

**Why `gappyout`?** The `gappyout` method selects a natural gap score threshold by finding
the first large drop in the sorted gap-score distribution. This automatic thresholding is
more robust across datasets than a fixed threshold (e.g., `-gt 0.1` = remove columns where
>10% have gaps). `gappyout` adapts to each dataset's gap distribution.

Alternative considered: `-automated1` (heuristic combining multiple methods). `gappyout`
was chosen for its transparent, interpretable behavior and consistent performance across
our benchmark taxa.

#### FastTree — Maximum Likelihood Tree

`run_fasttree(aligned_fasta, output_newick, model='GTR')`

Builds a maximum likelihood tree using FastTree's GTR+CAT model.

**Why FastTree over RAxML or IQ-TREE?**
- **Speed**: FastTree is 10–100× faster than RAxML or IQ-TREE for the haplotype-count
  range typical in BOLDGenotyper (20–500 sequences). RAxML with full bootstrap is
  appropriate for phylogenomics but overkill for a haplotype diagnostic tree.
- **Quality**: FastTree's CAT rate approximation is less accurate than RAxML's GAMMA model
  for precise divergence time estimation, but the primary purpose here is topology
  visualization and species-group verification, not publication-quality phylogenomics.
- **Conda availability and portability**: FastTree is a single self-contained binary.
- **GTR model**: General time-reversible is the standard substitution model for COI. It
  allows all six substitution rate parameters to vary freely, which is appropriate for
  COI's well-documented transition bias (transitions >> transversions).

**Midpoint rooting**: When no outgroup is specified (the default), the tree is midpoint-
rooted (longest branch bisected). This is the standard approach when the user doesn't
provide an outgroup. True rooting requires an outgroup; midpoint is a reasonable
approximation for visualizing overall topology.

---

### 12. Pop-Gen Export — `popgen_export.py`

#### Purpose

Export haplotype and population data in formats accepted by standard population genetics
software, allowing immediate downstream analysis without reformatting.

#### Arlequin (.arp)

`export_arlequin(haplotypes_fasta, mapping_df, output_path, ...)`

Arlequin (Excoffier & Lischer 2010) is one of the most widely used population genetics
analysis suites. BOLDGenotyper exports:
- `.arp` file: Full Arlequin project format with haplotype sequences, sample list, and
  population (geographic region) definitions
- `populations.txt`: Mapping of sample IDs to population labels

Populations are defined automatically from the geographic region assignments — samples from
each ocean basin, freshwater basin, or ecoregion form one population. This allows immediate
AMOVA (Analysis of Molecular Variance) without any manual formatting.

**Why Arlequin?** AMOVA is the standard test for hierarchical population structure in
molecular ecology. While BOLDGenotyper does not perform AMOVA internally (this would
require phylogeographic decisions about groupings that are dataset-specific), it exports
in ready-to-use Arlequin format so users can run AMOVA in 2 minutes after the pipeline
finishes.

#### PopART/NEXUS + traits

`export_popart_nexus(haplotypes_fasta, mapping_df, output_path, ...)`

PopART (Leigh & Bryant 2015) is the standard tool for haplotype network visualization
(TCS networks, median-joining networks, statistical parsimony networks).

Export includes:
- `.nexus` file: Haplotype sequences in NEXUS format with taxa named as
  `{genotype_id}_{sample_number}`
- `.traits` file: Trait matrix in PopART format mapping each haplotype to geographic
  region counts (so the network can be drawn with pie charts showing region frequencies)

#### DnaSP (.fas)

`export_dnsp_fasta()` exports a simple FASTA with haplotype sequences suitable for import
into DnaSP (Rozas et al. 2017) for nucleotide diversity (π), haplotype diversity (h),
Tajima's D, and other neutrality statistics.

**Why export to three programs rather than computing these statistics internally?**
1. Each program has established, well-tested implementations of these statistics
2. Users trust these tools and reviewers recognize them
3. Different users have different downstream needs — some want AMOVA, some want networks,
   some want neutrality tests
4. Implementing AMOVA, TCS networks, and Tajima's D correctly in Python is substantial
   extra code that duplicates existing tools
5. BOLDGenotyper's role is discovery and export; analysis is the user's domain

---

### 13. Metadata Analysis — `metadata_analysis.py`

#### Purpose

Analyze associations between haplotype identity and non-geographic metadata fields
(sex, life stage, collection date, etc.).

#### Default Fields Analyzed

`DEFAULT_CATEGORICAL_FIELDS = ['sex', 'life_stage', 'reproduction', 'country/ocean',
'country_iso', 'province/state', 'realm', 'biome', 'ecoregion', 'habitat', 'geoid']`

These fields are analyzed if present in the BOLD TSV. Fields not present are silently skipped.

#### Statistical Test: Chi-Square

For each categorical field, a chi-square test of independence tests whether haplotype
distribution differs across metadata categories (e.g., are haplotype frequencies
significantly different between male and female specimens?).

**Why chi-square?** Chi-square is the standard test for contingency table independence. It
is robust to unequal category sizes (unlike ANOVA) and requires only count data. Fisher's
exact test would be more appropriate for small cell counts, but chi-square scales to the
multi-haplotype, multi-category tables typical in this analysis. Small-cell warnings are
reported in the log.

#### Temporal Analysis

For `collection_date_start`, the module:
1. Parses dates in multiple formats (YYYY, YYYY-MM, YYYY-MM-DD, MM/DD/YYYY, etc.)
2. Extracts year
3. Plots haplotype frequency over time as a stacked bar chart
4. Tests for temporal trend using Spearman rank correlation per haplotype

**Why preserve raw values?** The module's philosophy is to present metadata as-is, without
standardization (the module docstring states this explicitly). Users know their data best.
For example, 'F', 'female', 'FEMALE', 'f' might all mean female — BOLDGenotyper does NOT
merge these. The optional `--normalize-sex` flag provides sex normalization as a deliberate
opt-in. This avoids silently merging categories that might not be equivalent.

---

### 14. Parameter Sweep — `parameter_sweep.py`

#### Purpose

Help users find the optimal `min_singleton_distance` threshold for their specific dataset
by running the pipeline at multiple thresholds and identifying the elbow point in the
haplotype count curve.

#### The Core Problem

The relationship between threshold and haplotype count is monotonic: higher threshold →
fewer singletons retained → fewer haplotypes. There is no single "correct" threshold. The
goal is to find the threshold where the curve shows diminishing returns — the elbow point —
where increasing the threshold further removes more genuine haplotypes than sequencing errors.

#### Elbow Detection Algorithm

`find_elbow_point(thresholds, n_haplotypes)`

1. Compute first derivative (slopes): `diff(n_haplotypes) / diff(thresholds)` — rate of
   haplotype loss per unit threshold increase
2. Compute second derivative (curvature): `diff(slopes) / diff(thresholds[:-1])` — how
   rapidly the rate of change itself is changing
3. Find the index of maximum absolute second derivative: this is where the curve "bends"
   most sharply (the elbow)
4. Return `thresholds[argmax(|curvature|) + 1]` — the threshold just after the bend

**Why maximum curvature (second derivative)?** Knee/elbow detection is a well-studied
problem. The "kneedle" algorithm (Satopaa et al. 2011) finds the point of maximum curvature
on a normalized curve, which is precisely what our second-derivative approach approximates.
Alternatives include distance-from-line methods, but second-derivative is simpler and
behaves well on smooth monotonic curves with a single clear bend.

**The `+1` offset**: We return the threshold one step PAST the maximum curvature, not at
it. This is because the curvature maximum marks the point where the curve is already bending,
meaning haplotypes are being lost rapidly. By stepping one threshold beyond it, we select
a threshold that has already passed the steepest loss region, retaining slightly fewer
singletons but more conservatively.

#### Acceptable Range

`generate_recommendations()` defines an "acceptable range" as all thresholds where the
resulting haplotype count is within ±20% of the elbow haplotype count:

```
range_lo = min threshold where count ≤ elbow_count × 1.2
range_hi = max threshold where count ≥ elbow_count × 0.8
```

This gives the user a window of thresholds that produce biologically similar results,
even if the exact elbow point selection is uncertain.

#### Stability Tracking

`track_group_membership(sweep_results)`

Tracks how haplotype composition changes across thresholds. For each pair of adjacent
thresholds, calculates Jaccard similarity between the set of haplotypes and their member
assignments. If Jaccard < 0.7, the transition is flagged as a composition change.

Stability scores per haplotype:
- **High**: 0 composition changes across the sweep range
- **Medium**: ≤2 composition changes
- **Low**: >3 composition changes

Stable haplotypes are confidently identified across a wide threshold range; unstable
haplotypes are at the boundary between "real variant" and "sequencing error" and deserve
scrutiny.

**Why 0.7 Jaccard threshold?** Jaccard = 0.7 means 70% overlap between adjacent threshold
states. A single member being reassigned in a 10-member haplotype gives Jaccard ≈ 0.82
(acceptable). Two members reassigned: Jaccard ≈ 0.67 (change detected). This is sensitive
enough to catch meaningful compositional shifts.

#### Sweep Plot

`create_threshold_stability_plot()` creates a 4-panel figure:
1. Haplotype count vs. threshold (with elbow marked)
2. Singletons retained vs. threshold
3. Suspect haplotypes flagged vs. threshold
4. Assignment coverage (% samples assigned) vs. threshold

This 4-panel view lets users see all relevant tradeoffs simultaneously.

---

### 15. Visualization — `visualization.py`

#### Purpose

Generate publication-quality figures: global distribution maps, regional abundance bar
charts, and faceted versions of both.

#### Color Palette Design

A custom 20-color palette approximates colorblind-friendly Okabe-Ito colors for the first
8 haplotypes, extending with seaborn's `tab20` for datasets with many haplotypes. The
primary colors match those used in the companion Sphyrna lewini PLOS ONE paper for
consistency:
- h1: Purple (#8545C1)
- h2: Teal (#10B3A5)
- h3: Yellow (#FFB031)
- h4: Blue (#3874F4)

Colors are assigned consistently across all figures in a run (same haplotype = same color
in map, bar chart, faceted plots, etc.).

#### Distribution Map

`plot_distribution_map(df, output_path, ...)`

Uses **cartopy** for cartographic projection. Samples are plotted as colored circles with
size optionally proportional to the count of samples at that location (useful for showing
where a haplotype is most abundant, not just where it occurs).

**Why cartopy?** Cartopy is the standard Python library for cartographic projections and
geographic visualization. It supports dozens of projections, handles the 180°/-180°
meridian correctly (common source of bugs in naive longitude-based plots), and integrates
with matplotlib.

Default projection: Robinson (for global maps). This projection minimizes area distortion
across most of the globe and is widely used in biogeographic publications.

#### Regional Abundance Bar Chart

`plot_ocean_basin_abundance(df, output_path, basin_column='ocean_basin', ...)`

Stacked bar chart showing haplotype proportions per geographic region. Works with any
geographic category by accepting any column name as `basin_column`. Region order:
- For GOaS ocean basins: curated biogeographic order (Arctic, North Atlantic, South
  Atlantic, Mediterranean, Indian Ocean, North Pacific, South Pacific, Southern Ocean, etc.)
- For custom regions: alphabetical

**Why stacked proportional bars (not counts)?** Regions have very different sample counts
(e.g., the North Atlantic is massively oversampled relative to the South Pacific in most
BOLD datasets). Proportional bars normalize for sampling effort, showing the distribution
of haplotypes WITHIN each region rather than absolute sample counts that reflect collection
bias.

#### Raw Data Export

Every figure exports a companion JSON file with the underlying data (e.g.,
`*_distribution_bar_data.json`). This serves two purposes:
1. Allows regeneration of figures (with different colors, fonts, etc.) without re-running
   the full pipeline
2. Enables the HTML report to use JavaScript charting with the same data

---

### 16. HTML Report — `reports.py`

#### Purpose

Generate a self-contained, interactive HTML report summarizing all analysis results for
non-specialist users (e.g., fisheries managers, quarantine officers) who may not analyze
the raw CSV/FASTA outputs.

#### Design Principles

- **Self-contained**: All JavaScript and CSS are either inline or from CDN. The HTML file
  can be emailed or shared without accompanying files (except for linked PDF/PNG figures).
- **Plotly.js for interactivity**: Tables are sortable/filterable; charts are hoverable.
  Plotly is loaded from CDN (no server needed).
- **Auto-generated Methods text**: The report includes a Methods section describing the
  pipeline settings used for the specific run. This text is generated from the config
  parameters so it is always accurate and reproducible. Users can copy-paste this section
  into their own manuscripts.

#### Report Contents (typical)

1. Executive summary: N sequences, N haplotypes, assignment coverage %
2. Haplotype table: sortable table with haplotype ID, member count, species, regions
3. Distribution map (embedded PNG)
4. Regional abundance chart (embedded PNG)
5. Divergence heatmap (embedded PNG, if tree was built)
6. Flagged samples: QC failures, suspect haplotypes, coordinate warnings
7. Parameter table: all configuration settings used
8. Auto-generated Methods section
9. Links to output files (CSV, FASTA, pop-gen exports)

**Why embed images as PNG rather than interactive Plotly?** The main publication figures
(maps, bar charts) are generated in matplotlib/cartopy because they produce better
publication-quality output than Plotly. The HTML report includes these as embedded PNGs.
Interactive Plotly is used only for the filterable tables.

---

### 17. Support Utilities — `utils.py`

#### Logging Setup

`setup_logging(log_level, log_file)` configures the `boldgenotyper` package logger with:
- Console output to stdout (for interactive use)
- Optional file output (written to `{output_dir}/{organism}_pipeline.log`)
- Format: `[YYYY-MM-DD HH:MM:SS] LEVEL: message`

Log files are essential for reproducibility — a complete log captures every parameter,
every filter decision, every external tool call, enabling exact reproduction of any run.

#### External Tool Management

`check_external_tool(tool_name, min_version)` checks PATH for MAFFT, trimAl, and FastTree
using `shutil.which`. If a tool is missing, it prints installation instructions (conda,
apt, brew). This is user-facing error handling designed to make the most common setup
problem (missing external tools) self-diagnosing.

**Tool version checking**: Some older MAFFT versions (<7.0) have different `--auto`
behavior. The version checker uses regex to extract version numbers from tool output
(trying `--version`, `-v`, etc.) and compares using tuple comparison of version components.

#### FASTA I/O

`read_fasta()` and `write_fasta()` are pure Python implementations requiring no BioPython.
This was a deliberate dependency reduction decision — BioPython is a heavy dependency (~50
MB) that many users have not installed. FASTA is a simple enough format that a custom
reader/writer is straightforward and fast.

The only place BioPython IS used is `divergence_analysis.py` for parsing aligned sequences
in divergence calculations. This could be replaced with the custom reader but was retained
for historical reasons.

#### ORF Quality Checking

`check_orf_quality(sequence, genetic_code=2)` translates the sequence using Python's
`str.translate()` with a custom translation table for NCBI genetic code 2 (vertebrate
mitochondrial code). The vertebrate mitochondrial code differs from the universal code in:
- TGA codes for tryptophan (not stop)
- ATA codes for methionine (not isoleucine)
- AGA and AGG code for stop (not arginine)

**Why code 2 as default?** The five benchmark taxa include sharks (code 2), lobsters (code
5 — invertebrate mitochondrial), butterflies (code 5), and salmon (code 2). Code 2 is the
most common for vertebrate COI barcoding. Users analyzing invertebrates should set
`--mitochondrial-code 5` (via config or CLI flag).

**ORF quality score components**:
- Has ATG start codon: +1
- N internal stop codons (in-frame): penalize
- Coverage fraction (length of unambiguous region / expected COI length): primary metric

The function returns both a score and a Boolean `is_valid` for use in QC filtering.

---

## Configuration Reference — `config.py`

All parameters are defined as Python dataclasses. The CLI maps command-line flags to these
config fields. Users can also provide a YAML config file to override defaults.

### Key Parameter Justifications

| Parameter | Default | Justification |
|---|---|---|
| `min_raw_length_abs` | 200 bp | Hard floor for meaningful barcode comparison |
| `min_raw_length_frac_of_median` | 0.70 | Adaptive threshold; removes truncated sequences while allowing shorter-amplicon datasets |
| `max_raw_N_fraction` | 0.05 | 5% ambiguous bases; typical Sanger quality cutoff |
| `require_valid_orf` | True | COI is protein-coding; ORF validation catches pseudogenes and NUMTs |
| `orf_max_internal_stops` | 2 | Allows for NNNN ambiguity artifacts while catching genuine pseudogenes |
| `core_min_coverage` | 0.80 | 80% of sequences must cover a position for it to be in the core region |
| `core_min_length` | 150 bp | Minimum viable length for barcode comparison |
| `mask_gap_threshold` | 0.50 | Columns with >50% gaps masked before core region calculation |
| `min_singleton_distance` | 0.005 (0.5%) | Removes singletons within 0.5% of nearest neighbor (sequencing error zone) |
| `max_singleton_distance` | 0.05 (5%) | Singletons >5% from nearest neighbor flagged as suspect |
| `tie_margin` | 0.003 (0.3%) | Tie flagged if best and runner-up differ by <0.3% identity |
| `tie_threshold` | 0.95 | Tie detection only active when best identity ≥ 95% |
| `majority_species_threshold` | 0.70 | 70% of members must share species name for assignment |
| `min_identity` | 0.50 | Assignment cutoff: samples with <50% identity to any haplotype are unassigned |
| `orf_min_coverage` | 0.50 | At least 50% of the expected COI length must translate successfully |
| `mafft_algorithm` | 'auto' | Automatic algorithm selection based on input size |
| `fasttree_model` | 'GTR' | Standard COI substitution model |
| `trimal_method` | 'gappyout' | Adaptive trimming; removes gap-heavy columns |
| `mitochondrial_code` | 2 | Vertebrate mitochondrial code; override to 5 for invertebrates |
| `phylo_min_length` | 600 bp | Only full-length sequences used for tree building |
| `phylo_min_cluster_size` | 5 | Minimum members for a haplotype to appear in tree |

---

## Design Decision Summary

### Why Python (not R)?

R packages (BOLDconnectR, vegan, ape, PopGenome) could theoretically cover similar ground,
but:
1. Workflow orchestration (calling MAFFT, FastTree as subprocesses) is more natural in Python
2. GeoPandas + Shapely for spatial joins has no R equivalent with comparable ease of use
3. Python's subprocess module handles external tools consistently across platforms
4. Plotly for interactive HTML reports is more mature in Python
5. The companion custom scripts (proto-BOLDGenotyper) were Python; this is a direct formalization

### Why not use BioPython throughout?

BioPython covers sequence alignment, phylogenetics, and file I/O. We use it only in
`divergence_analysis.py` because:
- BioPython is a heavy dependency (~50 MB) not universally installed
- FASTA I/O is simple enough for a custom reader (40 lines of code vs. a 50 MB import)
- Our MAFFT interface (subprocess) doesn't need BioPython's `Applications` module
- The custom `_compute_distance()` function gives us precise control over what counts
  as a valid site (excluding N and gaps consistently with barcoding conventions)

### Why no species delimitation?

BOLDGenotyper does NOT call species. This is intentional:
- Species delimitation (ASAP, ABGD, BIN system) makes claims about evolutionary entities
- Automatic species calls from COI alone are unreliable (the barcoding gap breaks down
  for many taxa)
- Users have taxon-specific knowledge about their study organisms
- BOLDGenotyper provides the haplotypes, divergences, and barcoding gap stats — the
  inputs for informed species delimitation — without making the call itself

### Why no AMOVA or neutrality tests internally?

Same logic: Arlequin (AMOVA) and DnaSP (Tajima's D, Fu's Fs) are established, peer-
reviewed implementations with well-understood statistical properties. Reimplementing
them in BOLDGenotyper adds code maintenance burden without adding value. The pop-gen
export module makes using these standard tools seamless.

### Why BOLD-specific (not iNaturalist, NCBI, iBOL directly)?

BOLD is the primary repository for COI barcoding data and provides the richest metadata
(collection dates, coordinates, life stage, specimen images) in a single download. NCBI
GenBank COI data can be imported by formatting it as a TSV with `processid` and `nuc`
columns. BOLDGenotyper is designed for BOLD's TSV format but works with any compliant TSV.

---

## External Tool Summary

| Tool | Version tested | Function | Why |
|---|---|---|---|
| MAFFT | 7.490+ | Multiple sequence alignment | Fastest+accurate for large datasets; --auto adapts to input size |
| trimAl | 1.4.1+ | Alignment trimming for phylogenetics | `gappyout` auto-thresholding; industry standard |
| FastTree | 2.1.11+ | ML phylogenetic tree | 10–100× faster than RAxML; adequate for diagnostic purposes |
| edlib | 1.3.9+ | Fast edit distance (Python C extension) | ~100× faster than pure Python Levenshtein for assignment step |
| GeoPandas | 0.14+ | Spatial join (point-in-polygon) | Standard Python geospatial stack; handles any shapefile |
| cartopy | 0.22+ | Map projections for distribution maps | Standard Python cartography; handles 180° meridian correctly |
| Plotly | 5.0+ | Interactive HTML report tables | Well-supported CDN-based interactive tables |
| scipy | 1.10+ | Hierarchical clustering for heatmaps | Standard scientific Python; `squareform`, `linkage` |
| seaborn | 0.12+ | Divergence heatmap | `clustermap` with dendrogram is simplest implementation |
| pandas | 2.0+ | Core data manipulation throughout | Standard Python data analysis library |
| numpy | 1.24+ | Numerical operations throughout | Standard Python numerical library |

---

*End of internal reference document.*
