# Migrating from BOLDGenotyper 1.0.0 to 2.0.0

This guide maps every v1 output file and column to its v2 equivalent, lists the
flag changes, and shows how to rebuild v1-style tables in downstream scripts.
`<org>` is the organism name that v1 put in front of most file names; v2 file
names have no prefix except the log.

Every v2 run folder also has a generated `README.md` that defines each file and
column of that run. When this guide and a run's README differ, the README
describes what that run wrote.

Contents

1. [What changed, in brief](#1-what-changed-in-brief)
2. [Run folder layout](#2-run-folder-layout)
3. [File map](#3-file-map)
4. [Column changes](#4-column-changes)
5. [Identifiers and values](#5-identifiers-and-values)
6. [Command-line changes](#6-command-line-changes)
7. [boldgenotyper-query](#7-boldgenotyper-query)
8. [boldgenotyper-sweep](#8-boldgenotyper-sweep)
9. [Python API](#9-python-api)
10. [Recipes for downstream scripts](#10-recipes-for-downstream-scripts)

---

## 1. What changed, in brief

- **One status per record.** Every unique input record has one of six
  statuses (`assigned`, `assigned_corrected`, `excluded_error_singleton`,
  `excluded_core_coverage`, `failed_qc`, `no_sequence`). v1 reported only
  `assigned` and `no_sequence`, so records removed by QC, core coverage or the
  error filter were labeled `no_sequence`.
- **Error-filtered singletons are corrected.** Records of singleton haplotypes
  removed as likely sequencing errors are assigned to their nearest retained
  haplotype when it is unique and within `--singleton-distance`
  (`assigned_corrected`); otherwise they stay unassigned.
- **Stable haplotype keys and count-free IDs.** `seq_md5` (md5 of the core
  sequence) is stable across runs on the same data; display IDs are `H0001`,
  `H0002`, ... (v1: `h1_n354`, whose count disagreed with other files).
- **Coordinates are parsed for every record.** Only records with
  geographic-quality coordinates (`coord_status` ok) are assigned a region and
  mapped; country- and province-centroid coordinates keep `lat`/`lon` but get
  no region. A record without a region has an empty value, never `Unknown`.
  v1 attached coordinates and regions to the wrong records when the input had
  duplicate processids (see 5.4): re-run rather than reuse v1 regions.
- **Labels, not contamination.** Label conflicts are reported with their
  evidence; "contamination", "confidence" and automatic depositor categories
  are gone. Members without a species label are left out of majority
  fractions, including the rule that assigns each haplotype's taxon.
- **A smaller default folder.** `records.csv` and `haplotypes.csv` replace
  the annotated table and seven per-haplotype or per-record tables;
  interpretive and descriptive extras need `--extended`; exports and the plot
  kit are opt-in; the report links figures instead of embedding them.
- **Provenance.** `run_provenance.json` records the full configuration (including
  `--singleton-distance`), tool versions, seeds, the input md5 and the run status.

Re-running v1 data with v2 gives different haplotype IDs, statuses and
assignment counts by design. Compare haplotypes across versions by sequence
(`seq_md5`, or the sequences in `haplotypes.fasta`), not by ID.

---

## 2. Run folder layout

v1 (default run):

```
<org>_BOLDGenotyper/
├── <org>_annotated.csv, <org>_pipeline.log, <org>_pipeline_parameters.json, <org>_summary_report.html
├── haplotypes/  haplotype_assignments/  taxonomy/  quality_control/  reports/
├── divergence_analysis/  phylogenetic/  geographic_analysis/  species_analysis/
├── metadata_summary/  visualization/{pdf,svg,json,metadata}/
├── plots/        (always)
└── exports/      (always; all formats)
```

v2:

```
<org>_BOLDGenotyper/
├── README.md                    every file and column, defined
├── run_provenance.json
├── <org>_pipeline.log
├── summary.csv                  assignment summary and species counts (one row)
├── records.csv                  one row per unique input record
├── records_qc_passed.fasta      QC-passed sequences (boldgenotyper-query --search-records)
├── haplotypes.fasta
├── haplotypes.csv               one row per haplotype
├── label_conflicts.csv
├── metadata_coverage.csv
├── divergence/   haplotypes_aligned.fasta, pairwise_p_distance.csv
├── phylogeny/    tree.nwk, tree_labeled.nwk, excluded_tips.txt, tree.pdf, tree.svg   (--build-tree)
├── geography/    haplotype_by_region.csv, distribution_map.pdf, distribution_map.svg  (with coordinates)
├── report.html
├── extended/     (--extended)       metadata/, divergence/, species/, geography/, labels/
├── exports/      (--export-format)  arlequin/, dnasp/, popart/, generic/
├── plots/        (--plot-kit)       data/, scripts/, plot_config.yaml
└── intermediate/ (--keep-intermediates)
```

---

## 3. File map

| v1 file | v2 file |
|---|---|
| `<org>_annotated.csv` | `records.csv` (see 4.1) |
| `<org>_pipeline.log` | `<org>_pipeline.log` (written fresh each run; ends with a warnings summary and the run status) |
| `<org>_pipeline_parameters.json` | `run_provenance.json` |
| `<org>_summary_report.html` | `report.html` (figures linked, not embedded: open it from the run folder) |
| `haplotypes/<org>_haplotypes.fasta` | `haplotypes.fasta` (header `>H0001 seq_md5=... n_members=...`) |
| `haplotypes/<org>_haplotype_mapping.csv` | `records.csv` |
| `haplotypes/<org>_haplotype_stats.csv` | `haplotypes.csv` |
| `haplotype_assignments/<org>_diagnostics.csv` | `records.csv` |
| `haplotype_assignments/<org>_species_composition.csv` | `haplotypes.csv` |
| `taxonomy/<org>_haplotype_taxonomy.csv` | `haplotypes.csv` |
| `taxonomy/<org>_species_by_haplotype.csv` | `haplotypes.csv` `labels`, or count `records.csv` (recipe 10.3) |
| `quality_control/potential_misidentifications.csv` | `label_conflicts.csv` |
| `quality_control/mixed_species_summary.csv` | `haplotypes.csv` |
| `quality_control/depositor_flags_summary.csv` | `records.csv` `depositor_keyword_hits`; `extended/labels/depositor_keyword_hits.csv` |
| `quality_control/contamination_heatmap.pdf` | `extended/labels/label_composition_heatmap.pdf` |
| `quality_control/purity_distribution.pdf` | `extended/labels/majority_fraction_distribution.pdf` |
| `quality_control/README.md`, `divergence_analysis/README.md`, `geographic_analysis/README.md` | `README.md` |
| `reports/<org>_assignment_summary.csv` | `summary.csv` |
| `divergence_analysis/<org>_consensus_aligned.fasta` | `divergence/haplotypes_aligned.fasta` |
| `divergence_analysis/pairwise_divergence_matrix.csv` | `divergence/pairwise_p_distance.csv` |
| `divergence_analysis/divergence_heatmap.pdf` | `extended/divergence/divergence_heatmap.pdf` |
| `divergence_analysis/divergence_summary.csv` | `extended/divergence/species_pair_divergence.csv` |
| `divergence_analysis/within_vs_between_species.csv`, `barcoding_gap.pdf` | `extended/divergence/label_based_gap.csv`, `label_based_gap.pdf` (computed from record labels) |
| `phylogenetic/<org>_tree.nwk` | `phylogeny/tree.nwk` |
| `phylogenetic/<org>_tree_relabeled.nwk` | `phylogeny/tree_labeled.nwk` |
| `phylogenetic/<org>_tree.pdf`, `.svg` | `phylogeny/tree.pdf`, `.svg` |
| — | `phylogeny/excluded_tips.txt` (haplotypes left out of the tree, with the reason) |
| `phylogenetic/<org>_MSA_chunk<N>.pdf`, `.png` | `extended/divergence/msa/MSA_chunk<N>.pdf` (no PNG) |
| `geographic_analysis/geographic_coverage.csv` | `haplotypes.csv` `n_with_coords`, `n_with_region`, `n_with_country` |
| `species_analysis/species_assignments.csv` | removed (`records.csv`) |
| `species_analysis/species_summary.csv` | `extended/species/species_summary.csv` |
| `species_analysis/species_diversity.csv` | `extended/species/species_diversity.csv` |
| `species_analysis/species_basin_distribution.csv` | `extended/species/species_by_region.csv` |
| `species_analysis/species_country_distribution.csv` | `extended/species/species_by_country.csv` |
| `species_analysis/species_geographic_range.csv` | removed |
| `species_analysis/species_subsets/<Species>_haplotypes.csv` | `extended/species/per_species/<Species>_records.csv` |
| `species_analysis/species_divergence_summary.csv` | `extended/divergence/species_divergence_summary.csv` (not written by v1; see 5.8) |
| `species_analysis/species_divergence_matrices/*.csv` | `extended/divergence/within_species/*.csv` |
| `metadata_summary/<org>_metadata_summary.csv` | `metadata_coverage.csv` |
| `metadata_summary/<org>_metadata_coverage.json` | removed |
| `metadata_summary/<org>_<field>_tabulation.csv` | `extended/metadata/<field>_tabulation.csv` |
| `metadata_summary/<org>_haplotype_emergence.csv` | `extended/metadata/haplotype_first_sampled.csv` |
| `metadata_summary/<org>_temporal_trends.csv`, `_temporal_summary.csv` | `extended/metadata/temporal_trends.csv`, `temporal_summary.csv` |
| `visualization/metadata/{pdf,png}/<org>_<name>.*` | `extended/metadata/figures/<name>.pdf`, `.png` |
| `visualization/metadata/svg/*` | removed |
| `visualization/metadata/*/<org>_haplotype_emergence_timeline.*` | `extended/metadata/figures/haplotype_sampling_timeline.*` |
| `visualization/metadata/*/<org>_emergence_timeline_{species,country}_faceted.*` | `extended/metadata/figures/sampling_timeline_{species,country}_faceted.*` |
| `visualization/{pdf,svg}/<org>_distribution_map.*` | `geography/distribution_map.*` (v1 often failed to write it) |
| `visualization/{pdf,svg}/<org>_<name>.*` (bars, faceted maps, species maps and bars) | `extended/geography/<name>.*` |
| `visualization/{pdf,svg}/species_haplotype_{maps,bars}/*` | `extended/geography/per_species/*` |
| `visualization/json/*.json` | removed (the data are embedded in `report.html`) |
| `plots/data/distribution_bar_absolute.csv` | `geography/haplotype_by_region.csv`; the plot kit keeps it with `--plot-kit` |
| `plots/data/consensus_sequences.fasta`, `plots/data/tree_data.nwk` | removed (the kit reads `../haplotypes.fasta`, `../phylogeny/tree.nwk`) |
| `plots/plot_customization_config.yaml` | removed (duplicate of `plots/plot_config.yaml`) |
| `plots/scripts/regenerate_identity.py` | removed |
| `exports/generic/haplotype_consensus.fasta` | removed (`haplotypes.fasta`) |
| other `exports/*`, `plots/*` files | same paths, written only with `--export-format` / `--plot-kit` |

---

## 4. Column changes

### 4.1 `<org>_annotated.csv` → `records.csv`

`records.csv` has one row per unique input record, in input order, including
records that failed QC or have no sequence (v1's annotated table had QC-passed
records only). It carries these BOLD fields, when present in the input:
`sampleid, bin_uri, insdc_acs, marker_code, inst, phylum, class, order, family,
genus, species, subspecies, identification, identification_method,
identification_rank, identified_by, taxonomy_notes, sex, life_stage,
reproduction, notes, short_note, collection_notes, collection_date_start,
collection_date_end, country/ocean, country_iso, province/state, region, site,
coord, coord_source, coord_accuracy, depth, elev, realm, biome, ecoregion,
habitat, nuc_basecount`. Join any other input column on `processid`
(recipe 10.1).

| v1 column | v2 |
|---|---|
| `haplotype_id` | `haplotype_id` (H IDs; empty unless `status` is `assigned` or `assigned_corrected`) |
| `assigned_haplotype` | `haplotype_id` |
| `haplotype_number` | removed |
| `haplotype_sp` | `haplotype_label` |
| `assigned_sp` | `assigned_sp` (the haplotype's taxon; rule changed, see 5.6) |
| `n_members`, `is_singleton`, `assignment_level`, `assignment_notes`, `majority_fraction` | `haplotypes.csv` |
| `is_tie`, `is_low_confidence` | removed (no identity re-assignment in the haplotype workflow) |
| — | `status`, `status_note`, `nearest_haplotype`, `distance_to_nearest` |
| `qc_pass`, `qc_fail_reason` | same |
| `raw_length` | `sequence_length` |
| `raw_N_count`, `raw_N_fraction`, `qc_pass_abs`, `qc_pass_median`, `orf_valid`, `qc_pass_orf` | removed |
| `lat`, `lon` | same (now parsed for every record) |
| `has_missing_coords`, `has_centroid_coords`, `has_zero_coords`, `has_invalid_coords`, `is_geographic_quality` | `coord_status` (`is_geographic_quality` is `coord_status == "ok"`) |
| `ocean_basin` (or the `--geo-category` column) | same name; empty instead of `Unknown` |
| `basin_confidence` | `region_confidence` |
| `distance_to_boundary_km` | `distance_to_region_boundary_km` |
| `basin_notes` | removed |
| `group_majority_species`, `group_majority_pct`, `group_species_count`, `matches_group_majority` | `haplotypes.csv` `majority_label`, `majority_fraction`, `n_species_labels`; record-level `label_conflict` |
| `potential_misidentification`, `misidentification_confidence` | `label_conflict` (evidence in `label_conflicts.csv`) |
| `depositor_uncertainty_flag`, `uncertainty_category`, `uncertainty_keyword` | `depositor_keyword_hits` (every keyword matched as a whole word, `; `-joined) |
| `depositor_notes` | `notes` |
| — | `label` (the record's species, or its genus when the species is missing) |
| `nuc` and BOLD columns not listed above | not in `records.csv`; join from the input |

### 4.2 Per-haplotype tables → `haplotypes.csv`

`haplotypes.csv` columns: `haplotype_id, haplotype_label, seq_md5, length,
n_members, n_exact, n_corrected, is_singleton, assigned_sp, assignment_level,
assignment_notes, assignment_source, labels, n_labeled, n_missing_label,
n_species_labels, majority_label, n_majority, majority_fraction,
n_label_conflicts, nearest_haplotype, nearest_distance, is_suspect,
suspect_reasons, min_orf_coverage, max_internal_stops, in_tree, n_with_coords,
n_with_region, n_with_country`, then the taxon evidence columns
`nearest_species_haplotype, nearest_species_taxon, nearest_species_distance,
taxon_max_within_distance, within_taxon_spread, member_bins, bin_taxon,
label_sources, tree_sister_taxon, tree_sister_support, suggested_sp,
suggestion_basis` (new in v2; see 5.6).

| v1 file: column | v2 `haplotypes.csv` |
|---|---|
| stats: `haplotype_id` | `haplotype_id` (`h1_n354` → `H0001`; see 5.1) |
| stats: `n_members`, `is_singleton` | same; `n_members` = `n_exact` + `n_corrected` |
| stats: `nearest_neighbor_distance` | `nearest_distance`; new `nearest_haplotype` (`;`-joined when tied) |
| stats: `suspect_distance`, `suspect_distance_reason`, `suspect_orf_reason`, `is_suspect`, `suspect_reasons` | `is_suspect`, `suspect_reasons` |
| — | `seq_md5`, `length`, `n_exact`, `n_corrected`, `min_orf_coverage`, `max_internal_stops`, `in_tree` |
| taxonomy: `assigned_sp`, `assignment_level`, `assignment_notes` | same (rule changed, see 5.6) |
| taxonomy: `haplotype_sp` | `haplotype_label` |
| taxonomy: `majority_fraction` | `majority_fraction` (label-based, see 5.6) |
| species_by_haplotype: `reported_species`, `n`, `frac`, `n_in_group` | `labels` ("Genus alpha 42; Genus 3; (no species label) 2"); recipe 10.3 |
| species_composition: `n_total_samples` | `n_members` |
| species_composition: `n_species` | `n_species_labels` |
| species_composition: `primary_species`, `primary_species_pct` | `majority_label`, `majority_fraction` (a fraction, not a "100.0%" string) |
| species_composition: `species_composition`, `species_list` | `labels` |
| species_composition: `is_multi_species` | `n_species_labels > 1` |
| species_composition: `is_ambiguous` | removed (`assignment_level`, `majority_fraction`) |
| mixed_species_summary: `consensus_group` | `haplotype_id` |
| mixed_species_summary: `n_samples`, `n_species`, `primary_species`, `primary_pct` | `n_members`, `n_species_labels`, `majority_label`, `majority_fraction` |
| mixed_species_summary: `contaminating_species`, `contamination_count`, `contamination_pct`, `contamination_breakdown` | `labels`, `n_label_conflicts`; `label_conflicts.csv` |
| mixed_species_summary: `flags` (CLEAN/CONTAMINATED/AMBIGUOUS), `confidence_level` | removed |
| — | `n_labeled`, `n_missing_label`, `n_majority` |
| geographic_coverage: `genotype` | `haplotype_id` |
| geographic_coverage: `with_lat_lon` | `n_with_coords` (now members with `coord_status` ok) |
| geographic_coverage: `with_ocean_basin` | `n_with_region` |
| geographic_coverage: `with_country` | `n_with_country` (v1 reported 0: it read a column BOLD exports do not have) |
| geographic_coverage: `species`, `total_samples` | `assigned_sp`, `n_members` |
| geographic_coverage: `pct_with_coords`, `pct_with_region`, `pct_with_country`, `representativeness` | removed (divide by `n_members`) |

### 4.3 Record-level tables → `records.csv`

| v1 file: column | v2 `records.csv` |
|---|---|
| mapping: `processid`, `haplotype_id` | same |
| mapping: `haplotype_number`, `n_members`, `is_singleton` | removed; `haplotypes.csv` |
| mapping: `orf_valid`, `orf_coverage`, `internal_stops` | removed (`haplotypes.csv` `min_orf_coverage`, `max_internal_stops`) |
| diagnostics: `assigned_haplotype` | `haplotype_id` |
| diagnostics: `status` (`assigned`, `no_sequence`) | `status` (six values; see 5.3) |
| diagnostics: `is_tie`, `is_low_confidence`, `assignment_method` | removed |
| diagnostics: `note` | `status_note` |

### 4.4 `reports/<org>_assignment_summary.csv` → `summary.csv`

| v1 | v2 |
|---|---|
| `total_samples` | `n_unique` (and `n_input`, data rows before duplicate removal) |
| `samples_with_sequence` | `n_unique - n_no_sequence` |
| `n_assigned` | `n_assigned + n_assigned_corrected` |
| `n_no_sequence` | `n_no_sequence` (now only records without a sequence; v1 also counted QC, core-coverage and error-filter exclusions here) |
| `pct_assigned` | `pct_assigned` (denominator: QC-passed records, `n_qc_pass`; v1 reported 100%) |
| `pct_no_sequence` | removed |
| `n_haplotypes` | same |
| `n_species_identified` | `n_species_haplotype_assigned`; see also `n_species_labels_input`, `n_species_labels_species_level`, `n_species_labels_qc_passed`, `n_species_confident` |
| — | `n_failed_qc`, `n_qc_pass`, `n_excluded_core`, `n_core_ok`, `n_assigned_corrected`, `n_excluded_error`, `n_excluded_error_tied`, `n_excluded_error_far`, `pct_assigned_exact` |

### 4.5 Label and QC tables

`potential_misidentifications.csv` → `label_conflicts.csv`:

| v1 | v2 |
|---|---|
| `processid`, `haplotype_id` | same |
| — | `haplotype_taxon`, `assignment_source` (the haplotype's `assigned_sp` and whether it was curated) |
| `species` | `label` |
| `group_majority_species`, `group_majority_pct` | `majority_label`, `majority_fraction` (a fraction) |
| `misidentification_confidence` | removed; evidence in `n_members`, `n_labeled`, `n_majority`, `nearest_haplotype_with_label`, `distance_to_nearest_with_label` |
| `depositor_uncertainty_flag` | `depositor_keyword_hits` |
| `depositor_notes` | `depositor_note` |

Rows: v1 listed records whose species differed from the group majority by
string; v2 lists label conflicts (5.6), so genus-only labels of the majority's
genus are no longer listed, and haplotypes with tied top labels list every
species-level member.

`depositor_flags_summary.csv` → `extended/labels/depositor_keyword_hits.csv`:
`genotype` → `haplotype_id`, `species_label` → `label`, `pattern_matched` →
`keyword_hits` (every whole-word match), `flag_category` removed,
`depositor_note` unchanged.

### 4.6 Divergence tables

- `pairwise_divergence_matrix.csv` → `divergence/pairwise_p_distance.csv`:
  same values; row and column names are haplotype IDs (`H0001` instead of
  `Genus species h1_n354`); taxa are in `haplotypes.csv`. The per-species matrices
  in `extended/divergence/within_species/` are keyed the same way.
- `within_vs_between_species.csv` (`species, n_within_comparisons, mean_within,
  sd_within, max_within, n_between_comparisons, mean_between, sd_between,
  min_between, barcoding_gap, gap_size`) → `extended/divergence/label_based_gap.csv`
  (`label, n_records, n_haplotypes, max_within, min_between,
  nearest_other_label, label_based_gap`), computed from record labels with
  label-conflict records left out.
- `divergence_summary.csv` → `extended/divergence/species_pair_divergence.csv`;
  columns unchanged.

### 4.7 Metadata tables

- `<org>_metadata_summary.csv` → `metadata_coverage.csv`; columns unchanged.
  Values are normalized first (5.7), so counts can differ from v1.
- Tabulations: columns unchanged (`haplotype_sp` holds the haplotype label).
- `haplotype_emergence.csv` → `haplotype_first_sampled.csv`:
  `first_collection_date` → `first_sampled`, `first_collection_year` →
  `first_sampled_year`, `last_collection_date` → `last_sampled`,
  `last_collection_year` → `last_sampled_year`, `years_observed` →
  `years_between_first_and_last_sample`; `first_location` and `last_location`
  are empty instead of `Unknown`.
- `temporal_trends.csv` starts with a `#` comment line: read it with
  `pd.read_csv(path, comment="#")`.

### 4.8 Species tables (`--extended`)

- `species_diversity.csv`: `simpson_index` removed; new `note` column
  ("label-based, sampling-dependent").
- `species_subsets/<Species>_haplotypes.csv` (all annotated columns, grouped by
  the haplotype's primary species) → `per_species/<Species>_records.csv`
  (`records.csv` columns, grouped by `assigned_sp`, species with two or more
  haplotypes).
- `species_summary.csv`, `species_by_region.csv`, `species_by_country.csv`:
  columns unchanged.

### 4.9 Plot kit and exports

- `plots/data/*.csv` keep their v1 column names (`consensus_group` is the
  haplotype ID, `consensus_group_sp` its label) so the regeneration scripts
  work; `geography/haplotype_by_region.csv` is the default equivalent of
  `distribution_bar_absolute.csv` with `haplotype_id`, `haplotype_label`,
  `assigned_sp`, `<region>`, `n_records`, `pct_of_region`.
- `exports/` files keep their v1 formats and columns; `haplotype_id` values are
  H IDs.

---

## 5. Identifiers and values

### 5.1 Haplotype IDs

| | v1 | v2 |
|---|---|---|
| ID | `h1_n354` (count in the ID, often disagreeing with other files) | `H0001` (no count) |
| Order | by v1 cluster order | `n_members` descending, then `seq_md5` |
| Stable key | none | `seq_md5`: md5 of the uppercase, ungapped core sequence |
| FASTA header | `>h1_n354 Haplotype 1 with 422 sequences` | `>H0001 seq_md5=... n_members=...` |
| Label | `Genus species h1_n354` (`haplotype_sp`) | `Genus species H0001` (`haplotype_label`) |

H IDs are renumbered in every run. To match haplotypes between runs or
versions, use `seq_md5` (v2) or the sequences themselves.

### 5.2 processid

v2 uses the bare BOLD process ID in every table and FASTA header. v1 FASTA
headers carried an organism prefix in some files.

### 5.3 Record statuses

| status | meaning | v1 |
|---|---|---|
| `assigned` | core sequence identical to a retained haplotype | `assigned` |
| `assigned_corrected` | error-filtered singleton assigned to its nearest retained haplotype, within `--singleton-distance` | `no_sequence` |
| `excluded_error_singleton` | error-filtered singleton left unassigned (tied, too far, or `--no-assign-corrected`) | `no_sequence` |
| `excluded_core_coverage` | core region too short | `no_sequence` |
| `failed_qc` | failed length, N or reading-frame checks (`qc_fail_reason`) | `no_sequence` (and absent from the annotated table) |
| `no_sequence` | no sequence in the input | `no_sequence` (absent from the annotated table) |

### 5.4 Coordinates and regions

`coord_status` is `ok`, `missing`, `unparseable`, `centroid`, `zero` or
`out_of_range`. Only `ok` records are assigned a region and mapped.
Country- and province-centroid coordinates keep `lat`/`lon` but are not
assigned a region or mapped. A record without a region has an empty value; v1
wrote `Unknown`.

v1 merged `lat`, `lon` and the region back onto the records by row position
instead of by `processid`. When the input had duplicate processids (common in
BOLD portal exports), every record after the first duplicate received the
coordinates and region of another record, or none; v1 also treated records
without coordinates as geographic quality. v1 regions, maps and region counts
from such inputs are unreliable, so compare v2 results with a v2 re-run, not
record by record with v1.

### 5.5 Countries

Country counts read the BOLD portal's `country/ocean` field (v1 read a
`country` column that portal exports do not have and reported zero).

### 5.6 Labels and the haplotype taxon

- A record's `label` is its species, or its genus when the species is missing.
- A haplotype's `majority_label` is its most common species-level member label.
  `majority_fraction` is that label's share of the counted labels: members
  with a species-level label plus genus-only or open-nomenclature labels of
  another genus. Members without a species label (no label, or a genus-only or
  open label of the majority's genus) are left out of the denominator.
- `assigned_sp` is `majority_label` when `majority_fraction` is more than the
  majority threshold (0.70); otherwise, or on a tie, the most common genus.
  v1 counted members without a label in the denominator, so some v1 genus-level
  haplotypes are species-level in v2. Open-nomenclature labels (`Genus sp. 1`)
  are no longer assigned as a haplotype's species.
- A label conflict is a member whose species-level label differs from the
  majority label, any species-level member of a haplotype whose top labels are
  tied, or a genus-only label of another genus.
- New in v2: evidence columns in `haplotypes.csv` and a `suggested_sp` for
  haplotypes left at genus (never applied), and `--taxon-curation FILE` to apply
  checked taxa (`assignment_source` = `curated`; label conflicts of a curated
  haplotype are judged against its taxon). See the README section "Reviewing and
  Curating Haplotype Taxa".

### 5.7 Metadata values

Before counting, empty and whitespace-only values are missing (not a
category); sex `M`/`F` → `Male`/`Female` and `U` → missing; life stage
`A`/`I`/`M`/`J` → `Adult`/`Immature`/`Mature`/`Juvenile`; reproduction `S`/`A`
→ `Sexual`/`Asexual`; all-lowercase values are capitalized; BOLD GIS suffixes
(`_simplify-..._buffered-...`) are removed; whole-number codes lose `.0`.
"Emergence" outputs are now "first sampled": they describe sampling dates.

### 5.8 Outputs v1 never wrote

v1 logged the species-level divergence steps as complete but wrote neither
`species_divergence_summary.csv` nor the per-species matrices (a haplotype-ID
mismatch); v1's main distribution map also failed silently. v2 writes them.

---

## 6. Command-line changes

### `boldgenotyper`

| Flag | v1 | v2 |
|---|---|---|
| `--extended` | — | new: write `extended/` (metadata cross-tabs, collection-date summaries, label-based gap, species tables, MSA plots, extra figures) |
| `--plot-kit` | — | new: write `plots/` |
| `--no-export-plot-data` | skip `plots/` (written by default) | removed (`plots/` only with `--plot-kit`) |
| `--export-format` | default `all` | default none: `exports/` only when given |
| `--no-export-format` | skip `exports/` | removed (`exports/` only with `--export-format`) |
| `--no-metadata-summary` | skip the metadata step | only drops the cross-tabs `--extended` adds; `metadata_coverage.csv` is always written |
| `--no-metadata-analysis` | alias of `--no-metadata-summary` | removed; use `--no-metadata-summary` |
| `--metadata-fields` | default included `realm`, `biome`, `ecoregion`, `geoid`, `country_iso` | default `sex`, `life_stage`, `reproduction`, `habitat`, `country/ocean`, `province/state`; name the others to summarize them (BOLD derives `realm`, `biome` and `ecoregion` from the coordinates using terrestrial ecoregions) |
| `--no-temporal-analysis` | skip temporal outputs | only restricts what `--extended` adds |
| `--tree-min-members` | — (fixed minimum of 3 members) | new, default 1: every non-suspect haplotype is a tip |
| `--no-assign-corrected` | — | new: leave error-filtered singletons unassigned |
| `--download-date` | — | new: recorded in `run_provenance.json` |
| `--taxon-curation` | — | new: apply checked haplotype taxa (5.6) |
| `--similarity-threshold`, `--tie-margin`, `--tie-threshold` | identity re-assignment | removed (records are mapped directly from haplotype discovery) |
| `--singleton-distance` | default 0.005 | unchanged; recorded in `run_provenance.json`. Always pass it explicitly. |

Config default changed: `msa.output_formats` is `["pdf"]` (v1: PDF and PNG). The
`genotype_assignment` config section is removed; a config file that still has it loads,
with a warning.

A removed flag stops the command with one line naming what to use instead (exit status 2).

### `boldgenotyper-query`

| Flag | v2 |
|---|---|
| `--min-query-coverage` | new, default 0.95 (of the shorter of query and haplotype) |
| `--max-mismatch-delta` | new, default 0: integer differences for the equivalence set; replaces `--tie-margin` |
| `--tie-margin` | removed; use `--max-mismatch-delta` |
| `--min-species-support` | new, default 0.70 |
| `--search-records` | new: search QC-passed records (`records_qc_passed.fasta`) |
| `--taxon-curation` | new: apply a curation file to the reference run's haplotype taxa for this query |
| `--haplotypes` | accepts a v2 run folder (`haplotypes.fasta` at the root), a v1 run folder (`haplotypes/*_haplotypes.fasta`; haplotype matching only) or a FASTA file |

---

## 7. boldgenotyper-query

Queries are aligned over their whole overlap with each haplotype (global
alignment with free end gaps), so mismatches at a query's ends are counted;
v1's local alignment dropped them.

`query_summary.csv` (read with `comment="#"`; the `#` lines record the
parameters and the reference run):

| v1 | v2 |
|---|---|
| `haplotype_call` | `call_level` (`haplotype`, `haplotype_set` or `none`) |
| `assigned_haplotypes` | `equivalence_set` (`;`-joined) |
| `n_tied` | `n_equivalent` |
| `best_identity` | same |
| `species_call` | same (empty unless species support reaches `--min-species-support`) |
| `assigned_species`, `candidate_species` | `pooled_member_labels` (every member label with its count), `species_support` |
| `reason` | `flags` |
| `assign_identity`, `tie_margin` | in the `#` header lines |
| — | `query_length`, `aligned_length`, `query_coverage`, `n_mismatches`, `n_gaps`, `set_identical_over_region`, `n_exact_record_matches`, `reference_run_id` |
| — | `call_basis`: `pooled labels`, or `haplotype taxa` when a haplotype in the set has a curated taxon (no call, flag `curated_taxa_disagree`, when the set's taxa differ) |

`query_results.csv`: `match_quality` removed; new `differences`,
`n_mismatches`, `n_gaps`, `indel_in_query`, `query_coverage`,
`overlap_coverage`. `query_results.json` gains a `reference_run` block
(md5 of `run_provenance.json`, version, number of haplotypes, run date, input
md5, download date). With `--search-records`, `query_record_matches.csv`
(`query_id, processid, match_type, label, haplotype_id, status, record_length`).

---

## 8. boldgenotyper-sweep

`sweep_summary.csv`:

| v1 | v2 |
|---|---|
| `n_species_detected` | `n_species_labels_qc_passed` |
| `pct_potential_misid` | `pct_label_conflict` |
| `n_mixed_groups` | `n_haplotypes_with_label_conflict` |
| `avg_group_purity` | `mean_majority_fraction` |
| `n_samples`, `n_haplotypes`, `n_singletons`, `n_suspect`, `pct_assigned` | same (`n_samples`: QC-passed records) |

`group_membership_tracking.csv` adds `md5_<threshold>` (the `seq_md5` of each
record's haplotype). `recommendations.txt` gives the haplotype-count table with
the elbow marked and, for each pair of adjacent thresholds, the share of
records whose haplotype changes; the v1 "acceptable range" is gone.

New option: `--genetic-code` (default 2). v1 sweeps always used code 2; pass the code
of the final run.

---

## 9. Python API

| v1 | v2 |
|---|---|
| `quality_control.add_contamination_columns` | removed; use `add_label_columns` |
| `quality_control.detect_depositor_uncertainty` | removed; use `find_keyword_hits` |
| `haplotype_query.classify_match_quality` | removed (identity labels are no longer reported) |
| `divergence_analysis.analyze_barcoding_gap` | removed; `compute_label_based_gap` (`--extended`) |
| `haplotype_assignment` (identity-based assignment), `config.GenotypeAssignmentConfig` | removed: records are mapped directly from haplotype discovery |
| `visualization.plot_identity_distribution`, `plot_identity_by_status`, `plot_assignment_status`; `reports.generate_taxonomy_conflicts_report`, `generate_consensus_characterization`, `generate_sequence_quality_metrics` | removed (identity scores and the v1 files they read no longer exist) |
| `reports.generate_html_report(organism, output_dir, version)` | adds `staging_dir`, `report_path`, `links`, `provenance`; writes `report.html` |
| `utils.setup_logging(..., mode="a")` | default `mode="w"` (a fresh log per run) |
| `utils.assign_consensus_taxonomy` | label-based rule (5.6) |
| — | `layout` (`build_records_table`, `build_haplotypes_table`, `FILE_MAP`), `run_readme`, `provenance`, `record_status`, `haplotype_ids`, `taxon_labels`, `taxon_evidence`, `taxon_curation` |

---

## 10. Recipes for downstream scripts

All examples use pandas; `run` is the run folder.

### 10.1 Rebuild a v1-style annotated table

```python
import pandas as pd

records = pd.read_csv(run / "records.csv", low_memory=False)
haps = pd.read_csv(run / "haplotypes.csv")
bold = (pd.read_csv("input.tsv", sep="\t", dtype=str)     # the input named in run_provenance.json
        .drop_duplicates("processid"))                     # the pipeline keeps the first record

extra = [c for c in bold.columns if c not in records.columns]
annotated = (records[records["qc_pass"].astype(str) == "True"]            # v1 kept QC-passed records
             .merge(bold[["processid"] + extra], on="processid", how="left")
             .merge(haps[["haplotype_id", "n_members", "is_singleton", "assignment_level",
                          "assignment_notes", "majority_fraction"]], on="haplotype_id", how="left")
             .rename(columns={"haplotype_label": "haplotype_sp"}))
```

### 10.2 Find a haplotype from another run by its sequence key

```python
key = "5ddac434db8fa629d0f0e1f61dec72a2"                  # seq_md5 from another run or version
hap_id = haps.loc[haps["seq_md5"] == key, "haplotype_id"]  # empty if the sequence is not here
hap_records = records[records["haplotype_id"].isin(hap_id)]
```

### 10.3 v1 `species_by_haplotype.csv`

```python
members = records[records["status"].isin(["assigned", "assigned_corrected"])]
by_hap = (members.groupby(["haplotype_id", "label"], dropna=False).size()
          .rename("n").reset_index())
by_hap["n_in_group"] = by_hap.groupby("haplotype_id")["n"].transform("sum")
by_hap["frac"] = by_hap["n"] / by_hap["n_in_group"]
```

(v1 tallied the `species` field; `label` falls back to the genus.)

### 10.4 v1 `geographic_coverage.csv` percentages

```python
cov = haps[["haplotype_id", "assigned_sp", "n_members", "n_with_coords",
            "n_with_region", "n_with_country"]].copy()
for c in ("n_with_coords", "n_with_region", "n_with_country"):
    cov["pct_" + c[2:]] = cov[c] / cov["n_members"] * 100
```

### 10.5 Haplotype counts per region (v1 `distribution_bar_absolute.csv`)

```python
by_region = pd.read_csv(run / "geography" / "haplotype_by_region.csv")
```

### 10.6 Read a v2 temporal table

```python
trends = pd.read_csv(run / "extended" / "metadata" / "temporal_trends.csv", comment="#")
```
