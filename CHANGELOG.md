# Changelog

All notable changes to BOLDGenotyper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - Unreleased

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
- A `docs` extra (`pip install -e ".[docs]"`) and Sphinx in `environment.yml`; `sphinx-build -b html docs/source docs/build/html` builds the API pages without warnings.

### Changed
- Haplotype IDs are `H0001`, `H0002`, ... ordered by member count, then `seq_md5`; v1's `h1_n354` IDs carried member counts that disagreed with other files (P0-8).
- `processid` is bare in every table and FASTA header (P0-10).
- Coordinates are parsed for every record, including current BOLD portal `coord` strings; only geographic-quality coordinates are assigned a region and mapped; country- and province-centroid coordinates keep `lat`/`lon` but get no region; records without a region have an empty value instead of `Unknown` (P0-4).
- Country counts read the portal's `country/ocean` field (P0-5).
- The input is parsed once and the log is written fresh for each run (P0-12).
- The haplotype taxon (`assigned_sp`) uses the label-based `majority_fraction`: members without a species label are left out of the denominator; open-nomenclature labels are no longer assigned as species (P1-6).
- Divergence outputs describe uncorrected p-distances without taxonomic threshold interpretations (P1-1).
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

### Deprecated
- `--similarity-threshold`, `--tie-margin` (pipeline), `--tie-threshold`: no effect in the haplotype workflow; they warn (P0-1).
- `boldgenotyper-query --tie-margin`: use `--max-mismatch-delta` (Q4).
- `--no-export-plot-data`, `--no-export-format`: no effect; the plot kit and exports are opt-in (P2).
- `quality_control.add_contamination_columns` (use `add_label_columns`), `quality_control.detect_depositor_uncertainty` (use `find_keyword_hits`), `haplotype_query.classify_match_quality`.

### Removed
- `<org>_annotated.csv`, the diagnostics, mapping, stats, taxonomy, species-composition and label-summary tables (folded into `records.csv` and `haplotypes.csv`) (P2).
- Duplicated outputs: plot-kit copies of the FASTA and tree, the duplicate plot config, the identity regeneration script, the generic export's consensus FASTA, `visualization/json/`, `species_assignments.csv`, the metadata coverage JSON, MSA PNGs, metadata SVGs (P2, P3).
- `potential_misidentifications.csv` confidence levels, CLEAN/CONTAMINATED/AMBIGUOUS flags, depositor-note categories (P1-6 to P1-8).
- `representativeness`, `species_geographic_range.csv`, `simpson_index` (P1-9, P1-10).
- The haplotype-consensus barcoding gap (`within_vs_between_species.csv`, `barcoding_gap.pdf`), replaced by the label-based gap with `--extended` (P1-2).
- `match_quality` in query results (Q7).

### Fixed
- Records removed by QC, core coverage or the error filter were reported as `no_sequence`, and the assignment summary reported 100% assigned (P0-2, P0-3).
- The main distribution map was not written (a gridline style error) and step failures were silent (P0-6).
- The tree silently dropped haplotypes with fewer than three members (P0-7).
- Most coordinates were not parsed, so most records had an `Unknown` region (P0-4); country counts were zero (P0-5).
- Species counts were computed differently in different files (P0-9).
- The query tool's local alignment dropped mismatches at the ends of a query and reported 100% identity (Q2).
- Species-level divergence summaries and per-species matrices were never written (haplotype-ID mismatch).
- The report never showed the metadata tabulations; the exports README named files that are not written; the sweep's example command used a nonexistent `--input` flag.
- Datasets that collapse to a single haplotype no longer fail.

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
