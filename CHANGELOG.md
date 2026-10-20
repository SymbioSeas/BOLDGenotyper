# Changelog

All notable changes to BOLDGenotyper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.1] - 2026-10-20

### Fixed
- The divergence heatmap (an `--extended` output) is drawn only with `--extended`, and is
  skipped above 20,000 haplotype-by-haplotype cells, as the label heatmap is. It was drawn on
  every run and every sweep threshold: about 25 minutes per run for a family-level Pieridae
  download (~1,860 haplotypes), and at low sweep thresholds a risk of running out of memory.
  Result files are unchanged.

## [2.1.0] - 2026-10-20

Query side only: reference runs are unchanged.

### Changed
- A query that runs past the haplotype core is compared over the part the core covers.
  Haplotypes are defined on a core region, which in some datasets starts tens of bases into
  the barcode (about 54 bp in an order-level shark dataset), so a 5' mini-barcode hung over
  the core's start and got no call. Query bases beyond a haplotype's own ends
  (`query_overhang`) are now left out of the coverage that decides eligibility
  (`compared_coverage`, a new column of `query_candidates.csv`); `overlap_coverage` keeps its
  2.0.0 meaning.
- A match must compare at least `--min-overlap-length` bases (new option, default 60), so a
  sequence that touches a haplotype's end over a few bases is never eligible.

### Added
- `query_candidates.csv`: `query_overhang` and `compared_coverage`.
- Flags `compared_over_part` (the best match covers neither the whole query nor the whole
  haplotype) and `short_overlap` (no match with enough coverage compares
  `--min-overlap-length` bases).
- `min_overlap_length` in `query_provenance.json` parameters and the methods paragraph.

### Fixed
- `boldgenotyper-labels suggest --backbone worms` no longer suggests renaming species that
  WoRMS writes with their subgenus (`Holothuria (Metriatyla) albiventer`): WoRMS marks the
  plain binomial an "alternative representation" of the accepted name, which is not a
  synonym. Suggested names are given without the subgenus. (Found on a *Holothuria*
  download, where 33 of 39 suggestions were of this kind.)
- `boldgenotyper-labels suggest` retries a request that fails or times out (3 tries) instead
  of ending the run: a slow WoRMS answer stopped a 198-name run partway.
- A large run no longer runs out of memory in step 6.2: the label composition heatmap (a
  staging figure, not a run output) is skipped above 20,000 haplotype-by-label cells. A
  family-level Pieridae download (1,250 haplotypes x 584 labels) was killed there.

## [2.0.0] - 2026-10-19

Output files and columns changed throughout; see
[docs/guides/MIGRATION_v1_to_v2.md](docs/guides/MIGRATION_v1_to_v2.md) for a
file-by-file and column-by-column map. Item IDs refer to the v2.0.0 release
plan (`docs/plans/2026-09-30-v2.0.0-release-plan.md`).

### Added
- `run_provenance.json`: the full resolved configuration (including `--singleton-distance`), the values the run applied, tool versions, seeds, input md5 and row count, an optional `--download-date`, warnings and the run status (P0-1).
- One status per unique input record: `assigned`, `assigned_corrected`, `excluded_error_singleton`, `excluded_core_coverage`, `failed_qc`, `no_sequence`, with `status_note`, `nearest_haplotype` and `distance_to_nearest` (P0-2).
- Error-filtered singletons are assigned to their nearest retained haplotype when it is unique and within `--singleton-distance` (`assigned_corrected`); `--no-assign-corrected` turns this off (P0-2).
- `summary.csv` with explicit denominators and named species counts (`n_species_labels_input`, `n_species_labels_species_level`, `n_species_labels_qc_passed`, `n_species_haplotype_assigned`, `n_species_confident`) (P0-3, P0-9).
- `seq_md5`, a haplotype key that is stable across runs; FASTA headers `>H0001 seq_md5=... n_members=...` (P0-8).
- `coord_status` for every record (`ok`, `missing`, `unparseable`, `centroid`, `zero`, `out_of_range`) (P0-4).
- `--tree-min-members` (default 1) and `phylogeny/excluded_tips.txt` listing haplotypes left out of the tree with the reason (P0-7).
- A warnings summary and the run status at the end of every log; a warning when the output folder already holds a run (P0-6, P0-12).
- `--extended` for interpretive and descriptive extras (metadata cross-tabs, collection-date summaries, label-based barcoding gap, species tables, MSA plots, additional maps and bars) (P1-2 to P1-5, P1-10).
- `label_conflicts.csv` with evidence for each conflict (majority label and fraction, members, nearest haplotype carrying the record's label) (P1-6, P1-7).
- Metadata vocabulary normalization before counting (P1-5).
- The v2 run folder: `records.csv` (one row per unique record), `haplotypes.csv` (one row per haplotype), `geography/haplotype_by_region.csv`, and a generated `README.md` that defines every file and column (P2).
- `--plot-kit` to write the plot regeneration kit (P2).
- `boldgenotyper-query`: overlap alignment over the whole query, `differences`, `n_mismatches`, `n_gaps`; equivalence sets with `--max-mismatch-delta`; species support from pooled member labels with `--min-species-support`; `--min-query-coverage`; the reference run's provenance md5 in every output; `--search-records` over `records_qc_passed.fasta` (Q1-Q10).
- `boldgenotyper-sweep`: per-pair share of records whose haplotype changes between adjacent thresholds; `md5_<threshold>` columns in `group_membership_tracking.csv` (P3).
- `boldgenotyper-sweep --genetic-code`, so a sweep checks reading frames with the same code as the final run (v1 sweeps always used code 2).
- A log warning when `--geo-category` has the name of an input field (e.g. BOLD's `ecoregion`), which the region column replaces.
- `docs/guides/MIGRATION_v1_to_v2.md`.
- Taxon evidence in `haplotypes.csv`: nearest species-level haplotype and whether it lies within that taxon's spread, member BINs and the BIN's majority label, label sources (records, projects, institutions), the tree sister taxon and its support, and for haplotypes left at genus a `suggested_sp` with its `suggestion_basis`.
- Evidence taxa: a genus-level haplotype takes the species two of distance, BIN and tree name (none against) when it is the haplotype's leading label in the label rule's genus (`assignment_source` = `evidence`; `--no-evidence-taxa` reports evidence without applying it; `applied.taxon_evidence` in `run_provenance.json`). Suggestions not applied are written to `taxon_curation_suggestions.csv`, a curation file to edit and pass to `--taxon-curation`.
- `--taxon-curation FILE` on `boldgenotyper` and `boldgenotyper-query`: checked haplotype taxa applied after the label rule, matched by `seq_md5` or sequence, with `assignment_source` (`labels`/`evidence`/`curated`; a curated taxon takes precedence) in `haplotypes.csv`, `haplotype_taxon` and `assignment_source` in `label_conflicts.csv`, and the file's md5 and rows not applied in `run_provenance.json`.
- `boldgenotyper-query` output folder: `query_calls.csv` (the result to report), `query_candidates.csv`, `query_alignments.txt`, `query_provenance.json` (reference run, parameters, curation and metadata files with md5, figures, warnings), a `README.md` with a methods paragraph and every column defined, `query_report.html` (calls, one evidence card per query, figures, methods), and `query_figures/`: the queries on the reference tree (pruned to each called genus), maps per called taxon, reference records per region (`query_regions.csv`) and sampling years, each as a PDF plus one SVG per page.
- `boldgenotyper-query --query-metadata FILE` (query coordinates, dates, country and region; a query's region from its value or from the reference run's own shapefile) and `--json`.
- `boldgenotyper-query`: `call_basis` (`pooled labels` or `haplotype taxa`); a set containing a haplotype with a curated or evidence taxon is called from the haplotypes' taxa, with the flag `haplotype_taxa_disagree` when they differ; query provenance records the reference run's evidence taxa (`reference_taxon_evidence`).
- A `docs` extra (`pip install -e ".[docs]"`) and Sphinx in `environment.yml`; `sphinx-build -b html docs/source docs/build/html` builds the API pages without warnings.
- `--marker` (default `COI-5P`) and `--allow-malformed-rows` on `boldgenotyper` and `boldgenotyper-sweep`; `summary.csv` counts the rows removed while reading the input (`n_other_marker_rows`, `n_malformed_rows`, `n_duplicate_rows`), and `run_provenance.json` records them under `applied.input_rows`.
- `label_variants.csv` in every run: pairs of labels in one genus that may be one name written differently (one or two letters, a gender ending, a qualifier), with record counts and shared haplotypes. Labels are tidied as they are read (whitespace, underscores, case); names are never merged automatically.
- `--label-map FILE` on `boldgenotyper` and `boldgenotyper-sweep`: confirmed label replacements (`from`, `to`, `reason`, optional `source`) applied before labels are pooled; `records.csv` `label_mapped_from`; the map's md5 and unused rows in `run_provenance.json`.
- `boldgenotyper-labels suggest RUN_DIR`: a suggested label map from GBIF or WoRMS (spelling fixes and synonyms with the backbone and date), outside the pipeline run.
- Error-filter diagnostics: `sources.csv` (one row per source and label: records, distinct core sequences and a `haplotype_only` flag), `records.csv` `label_in_retained` for error-filtered singletons, `summary.csv` `pct_error_filtered`, `n_singleton_label_absent`, `n_haplotype_only_groups` and `n_singletons_exempt`, and log warnings when the filter may be removing real variation.
- The error filter no longer treats as errors the singletons of studies that deposit one sequence per haplotype (a source, i.e. BOLD project or GenBank accession prefix, and label with at least 10 records, at least 90% of them distinct core sequences; `haplotype.haplotype_only_min_records` and `haplotype_only_min_distinct`). In simulations with known haplotypes, filtering such studies kept 16% of their true haplotypes.

### Changed
- The HTML report is rebuilt (`boldgenotyper.report`): generated after the output layout from the final run folder, with an Overview and key findings (counts only), Records and haplotypes, Taxonomy and labels with a "Needs attention" list (haplotypes not at species level with their evidence and suggested taxa, label conflicts per haplotype and label, curated taxa, curation rows not applied), Phylogeny and Geography with genus and species menus, Sampling and metadata, Divergence (`--extended`), Methods (the statement now states the taxon rule and any curation file) and Output files (from the run README's definitions). Tables show ten rows and the true count of the rest. No external scripts; works offline.
- Figures are organized by species: each has an overview (haplotypes up to 40 haplotypes, taxa above; at most 12 colored items, the rest "Other (n)") and one page per species or genus (a multi-page PDF and one SVG per page), with a `manifest.json` per figure folder. New: `phylogeny/tree_overview` (taxa that form a clade collapsed, others listed), `tree_by_genus.pdf`, a full tree sized 0.14 in per tip; `geography/distribution_map_by_species.pdf`; with `--extended`, `region_abundance_{relative,total}` and metadata views per field, records by collection year (heatmap), sampling timeline, records by year and region, and metadata coverage. "By region" metadata views use the run's region column (country/ocean when fewer than two regions). Figures are PDF and SVG; no PNGs outside the plot kit. `VisualizationConfig` gains `max_overview_haplotypes`, `max_colored_items`, `max_overview_species`, `max_heatmap_species`.
- `boldgenotyper-query` output names: `query_summary.csv` → `query_calls.csv` and `query_results.csv` → `query_candidates.csv`, both without `#` comment lines (their content moved to `query_provenance.json`); `query_results_detailed.txt` → `query_alignments.txt`; `query_results.json` only with `--json`. The writers moved from `haplotype_query` to `query_outputs`.
- Trees are rooted: on the outgroup when one is given (`--phylo-outgroup-fasta` now roots on its sequences; before, they were only added to the tree), otherwise at the midpoint of the longest path between two tips. The configuration's `midpoint_root: true` was recorded but never applied in 1.0.0. `run_provenance.json` records the rooting (`applied.tree_rooting`).
- The tree evidence line collapses branches with support below 0.70 before reading a haplotype's sister clade, so it names a taxon only from a supported clade.
- The tree is built in Phase 5 (step 5.3), before any step reads a haplotype's taxon; Phase 8 writes the relabeled tree and alignment plots. Warning steps are `5 tree`, `5 tree reroot` and `5 taxon evidence`.
- Species aggregation counts every haplotype assigned at species level; curated and evidence taxa were left out when their label share was below 0.70. `species_analysis.aggregate_samples_by_species` loses `min_confidence`.
- Haplotype IDs are `H0001`, `H0002`, ... ordered by member count, then `seq_md5`; v1's `h1_n354` IDs carried member counts that disagreed with other files (P0-8).
- `processid` is bare in every table and FASTA header (P0-10).
- Coordinates are parsed for every record, including current BOLD portal `coord` strings; only geographic-quality coordinates are assigned a region and mapped; country- and province-centroid coordinates keep `lat`/`lon` but get no region; records without a region have an empty value instead of `Unknown` (P0-4).
- Country counts read the portal's `country/ocean` field (P0-5).
- The input is parsed once and the log is written fresh for each run (P0-12).
- The haplotype taxon (`assigned_sp`) uses the label-based `majority_fraction`: members without a species label are left out of the denominator; open-nomenclature labels are no longer assigned as species (P1-6).
- Divergence outputs describe uncorrected p-distances without taxonomic threshold interpretations (P1-1).
- `divergence/pairwise_p_distance.csv` and the per-species matrices are keyed by `haplotype_id` (display labels change when a taxon changes).
- Default metadata fields: `sex`, `life_stage`, `reproduction`, `habitat`, `country/ocean`, `province/state`. BOLD's `realm`, `biome` and `ecoregion` (which BOLD derives from the coordinates using terrestrial ecoregions), `geoid` and `country_iso` are summarized only when named in `--metadata-fields`.
- Depositor notes are searched for keywords as whole words; every hit is reported, without categories (P1-8).
- "Emergence" outputs are renamed "first sampled" and describe sampling dates (P1-4).
- Geographic coverage is reported as counts and percentages without grades (P1-9).
- Default runs write metadata coverage only; cross-tabs and collection-date summaries need `--extended`; `--no-metadata-summary` and `--no-temporal-analysis` only restrict what `--extended` adds (P1-3).
- The default run folder holds the files listed in the run README; extras move to `extended/`, exports to `exports/` (only with `--export-format`), the plot kit to `plots/` (only with `--plot-kit`) (P2).
- `report.html` (was `<org>_summary_report.html`) links figures instead of embedding them, and its Methods tab describes the haplotype workflow with the run's values (P2).
- `--export-format` has no default (v1: `all`) (P2).
- PNG figures are capped at 150 dpi and 3000 px per side; MSA chunk plots are PDF only; metadata figures are PDF and PNG (P3).
- Log lines are plain text (P3).
- The sweep's `recommendations.txt` reports the elbow with the stability table instead of an "acceptable range" (P3). With fewer than 4 thresholds it says there is no elbow (v1 named the middle value, or the only one); `parameter_sweep.find_elbow_point` then returns None.
- The conda environment installs pymsaviz for MSA plots.
- Documentation: README, API reference, parameter-sweep, popgen-export and full-dataset guides, the Sphinx sources (one page per module) and the tutorial notebooks describe 2.0.0; the benchmarking protocol notes what changed since the 1.0.0 benchmarks.
- `COMPARATIVE_ANALYSIS_GUIDE.md` compares a species-level and a family-level run by matching their `records.csv` files on `processid`; the docs no longer refer to `scripts/compare_analyses.py`, which is not distributed.
- `FULL_DATASET_COMMANDS.md` gives the shapefiles the case-study runs use (BasinATLAS level 3 `PFAF_ID` for Salmonidae, Ecoregions2017 `REALM` for Pieridae, with the region column `biogeographic_realm` so BOLD's `ecoregion` field is kept).

### Removed
- Faceted, species-faceted and per-haplotype maps and bars, country-faceted metadata figures, and PNG copies of figures (replaced by the per-species pages); the plotting functions behind them in `visualization` and `metadata_summary`.
- `VisualizationConfig` options no figure reads (`color_palette`, `reference_colors`, `figure_dpi`, `figure_format`, `map_projection`, `map_figsize`, `barplot_figsize`, `font_family`, `font_size`, `show_bootstrap_threshold`, `facet_by`, `map_buffer_degrees`, `show_unknown_geography_annotation`, `show_scale_bar`, `facet_dpi`, `facet_formats`, `save_individual_facets`, `intermediate_dpi`); a config file that sets one loads with a warning.
- Flags: `--similarity-threshold`, `--tie-margin`, `--tie-threshold` (no effect: records are mapped directly from haplotype discovery), `--no-export-plot-data`, `--no-export-format` (the plot kit and exports are opt-in), the `--no-metadata-analysis` alias (use `--no-metadata-summary`), and `boldgenotyper-query --tie-margin` (use `--max-mismatch-delta`, Q4). A removed flag stops the command with one line naming its replacement.
- Functions: `quality_control.add_contamination_columns` (use `add_label_columns`), `quality_control.detect_depositor_uncertainty` (use `find_keyword_hits`), `haplotype_query.classify_match_quality`, `divergence_analysis.analyze_barcoding_gap` (use `compute_label_based_gap`).
- `species_analysis.get_species_sample_list` and `get_species_haplotype_list` (no callers; they still filtered on the 0.70 label share).
- `metadata.assign_genotypes` and its edit-distance helpers (`levenshtein_distance`, `calculate_edit_distance`), which the pipeline did not call; `edlib` is no longer a dependency.
- The v1 identity-based assignment module `haplotype_assignment`, the `genotype_assignment` config section (a config file that still has it loads, with a warning) and the provenance `not_applied` entry; identity plots and report helpers that read v1 files.
- `<org>_annotated.csv`, the diagnostics, mapping, stats, taxonomy, species-composition and label-summary tables (folded into `records.csv` and `haplotypes.csv`) (P2).
- Duplicated outputs: plot-kit copies of the FASTA and tree, the duplicate plot config, the identity regeneration script, the generic export's consensus FASTA, `visualization/json/`, `species_assignments.csv`, the metadata coverage JSON, MSA PNGs, metadata SVGs (P2, P3).
- `potential_misidentifications.csv` confidence levels, CLEAN/CONTAMINATED/AMBIGUOUS flags, depositor-note categories (P1-6 to P1-8).
- `representativeness`, `species_geographic_range.csv`, `simpson_index` (P1-9, P1-10).
- The haplotype-consensus barcoding gap (`within_vs_between_species.csv`, `barcoding_gap.pdf`), replaced by the label-based gap with `--extended` (P1-2).
- `match_quality` in query results (Q7).

### Fixed
- Records removed by QC, core coverage or the error filter were reported as `no_sequence`, and the assignment summary reported 100% assigned (P0-2, P0-3).
- The main distribution map was not written (a gridline style error) and step failures were silent (P0-6).
- A distribution map failed to save when every plotted location stood for the same number of records (the marker sizes divided by zero); with `--extended` this dropped species maps (P0-6).
- Maps and the plot kit drew country- and province-centroid and `[0, 0]` coordinates, although the log and documentation said these are not mapped; they keep `lat`/`lon` in `records.csv` (P0-4).
- The tree silently dropped haplotypes with fewer than three members (P0-7).
- Coordinates and regions were attached to the wrong records (P0-4). `lat`, `lon` and the region were merged back by row position instead of by `processid`, so when the input contained duplicate processids (common in BOLD portal exports) every record after the first duplicate received the coordinates and region of another record, or none: in a 500-record Sphyrnidae subsample with 14 duplicates, 61 of the 101 records with `lat`/`lon` had no coordinate of their own. The geographic-quality flag also checked columns that do not exist, so records without coordinates counted as geographic quality. Region assignments, maps and region counts from 1.0.0 runs on inputs with duplicate processids are unreliable; re-run them. Country counts were zero (P0-5).
- Species counts were computed differently in different files (P0-9).
- The query tool's local alignment dropped mismatches at the ends of a query and reported 100% identity (Q2).
- Species-level divergence summaries and per-species matrices were never written (haplotype-ID mismatch).
- The report never showed the metadata tabulations; the exports README named files that are not written; the sweep's example command used a nonexistent `--input` flag.
- Datasets that collapse to a single haplotype no longer fail.
- The report's Methods and Parameters tabs read the provenance written at the start of the run, so every count read "not recorded" and a built tree was reported as not built; they now use the run's final provenance.
- Report tables said "Showing first N rows of N total" with the truncated count; they give the true count.
- With pandas 3, text fields were not stripped of surrounding whitespace on input (pandas 3 reads text as the `str` dtype); all-empty ORF columns were kept.
- Rooting on an outgroup put support values on the wrong branches (Biopython moves supports with nodes, not branches, when it reroots) and made an internal outgroup clade the root itself; rooting now keeps each support on its split and places a two-way root on the outgroup's branch. Trees rooted with `--phylo-outgroup-label` or `--phylo-outgroup-taxon` in 1.0.0 have misplaced supports along the path to the root.
- Rooted and relabeled trees (`tree.nwk`, `tree_labeled.nwk`) keep FastTree's three-decimal support values; 1.0.0's `tree_labeled.nwk` rounded them to two.
- A specimen exported with several markers kept whichever row came first, so a specimen whose 12S, 16S, ND2 or COII row came before its COI row lost its COI; rows of other markers also entered QC. Rows of the target marker are now chosen before duplicate process IDs are removed.
- A BOLD TSV with rows of the wrong field count (e.g. a block of 10,000 rows with extra columns in one portal export) was read with a fallback that silently dropped those rows and reported them as duplicate process IDs. Every row is now read or the run stops naming the first bad line; fields are read with no quoting, so a value starting with a quote character cannot swallow the following lines.

## [1.0.0] - 2025

### Added
- **`--singleton-distance` CLI flag**: Exposes `haplotype.min_singleton_distance` as a command-line parameter, closing the intended workflow of running `boldgenotyper-sweep` to identify the optimal threshold and passing it directly to the full pipeline run via `--singleton-distance`. Default: 0.005 (0.5%).
- **Custom Shapefile Support**: `--custom-shp`, `--shp-field`, `--geo-category` for geographic analysis with any polygon shapefile (freshwater basins, ecoregions, watersheds, biomes, etc.)
- **Metadata Analysis Module**: Automated haplotype-metadata association analysis with chi-square tests, coverage assessment, and temporal emergence timelines (`metadata_analysis.py`)
- **Population Genetics Export**: Export to Arlequin (.arp), PopART/NEXUS, DnaSP (.fas), and generic CSV/FASTA via `--export-format`
- **Outgroup Rooting**: Three rooting modes — `--phylo-outgroup-fasta`, `--phylo-outgroup-label`, `--phylo-outgroup-taxon`
- **Haplotype Query Tool**: `boldgenotyper-query` command for assigning new sequences to existing haplotypes
- **Parameter Sweep Tool**: `boldgenotyper-sweep` command for optimizing singleton filtering thresholds
- **Plot Customization**: YAML-configured plot regeneration from exported data (`plot_export.py`, `plot_regeneration.py`)
- **Species-Faceted Analysis**: Per-species haplotype distribution maps and bar charts
- **MSA Visualization**: Phylogeny-ordered multiple sequence alignment plots (`msa_visualization.py`)
- **Divergence Analysis**: Pairwise divergence matrices, barcoding gap analysis, within/between-species divergence
- **ESV Approach**: Exact Sequence Variant haplotype discovery (replaces clustering-based approach)
- **COI Validation**: Orientation normalization and ORF validation with configurable genetic codes
- **Dynamic QC**: Median-based adaptive quality control filtering
- **Singleton Error Filtering**: Two-stage singleton quality control with configurable thresholds

### Changed
- Geographic analysis generalized from marine-only to any organism via custom shapefiles
- Haplotype discovery uses ESV approach instead of clustering thresholds
- Version bumped to 1.0.0 for publication release
- Python requirement raised to >=3.9 (3.8 is EOL)
- FastTree is the only supported tree-building algorithm

### Removed
- PhyML and RAxML tree-building options (only FastTree is used)
- `genotype_assignment.py` (replaced by `haplotype_assignment.py`)
- `cluster_diagnostics.py` (functionality integrated elsewhere)
- `comparative_analysis.py` (standalone script at `scripts/compare_analyses.py`)
- Dead configuration classes: `DelimitationConfig`, `OutlierConfig`, `DistanceConfig`

## [0.1.0] - 2025-01-20

### Added
- Core pipeline: BOLD TSV to annotated results
- Sequence dereplication with MAFFT alignment
- Edit distance-based genotype assignment with CIGAR string parsing
- GOaS ocean basin assignment for marine organisms
- Phylogenetic tree building with FastTree
- Geographic distribution maps with Cartopy
- Interactive HTML reports with Plotly.js
- Publication-ready visualizations (PDF/SVG)
- Hierarchical configuration system with dataclasses
- Comprehensive CLI with automatic organism name detection

---

## Version Numbering

BOLDGenotyper uses [Semantic Versioning](https://semver.org/):
- **MAJOR** version for incompatible API changes
- **MINOR** version for new functionality in a backwards compatible manner
- **PATCH** version for backwards compatible bug fixes

## Release Process

1. Update version number in `setup.py`, `boldgenotyper/__init__.py`, and `CITATION.cff`
2. Update this CHANGELOG.md with release date and changes
3. Create git tag: `git tag -a v1.0.0 -m "Release version 1.0.0"`
4. Push tag: `git push origin v1.0.0`

For detailed information about any release, see the corresponding GitHub release notes at:
https://github.com/SymbioSeas/BOLDGenotyper/releases
