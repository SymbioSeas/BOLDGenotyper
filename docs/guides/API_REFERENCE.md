# BOLDGenotyper API Reference

Version 2.0.0. The command-line tools are the supported interface; this page describes
the Python functions behind them for scripting and custom analyses. Every run folder's
`README.md` defines the output files and columns; the
[migration guide](MIGRATION_v1_to_v2.md) maps v1 outputs and functions to v2.

## Table of Contents

- [Installation](#installation)
- [Running the Pipeline](#running-the-pipeline)
- [Reading a Run Folder](#reading-a-run-folder)
- [Configuration](#configuration)
- [Module Reference](#module-reference)
- [Examples](#examples)

---

## Installation

```bash
# conda (includes MAFFT, trimAl, FastTree and pymsaviz)
conda env create -f environment.yml
conda activate boldgenotyper
pip install -e .

# pip only (install MAFFT, trimAl and FastTree separately)
pip install -e ".[msa]"
```

```python
import boldgenotyper
print(boldgenotyper.__version__)   # 2.0.0
```

---

## Running the Pipeline

### From Python via the CLI (recommended)

```python
import subprocess

subprocess.run([
    "boldgenotyper", "data/Sphyrna_lewini.tsv",
    "--output", "results/Sphyrna_lewini",
    "--singleton-distance", "0.015",
    "--build-tree", "--threads", "8",
], check=True)
```

### `cli.run_pipeline`

```python
from pathlib import Path
from boldgenotyper import cli, config

out = Path("results/Sphyrna_lewini")
cfg = config.get_default_config().update(
    output_dir=out,
    n_threads=8,
    extended=True,
    haplotype__min_singleton_distance=0.015,
    phylogenetic__build_tree=True,
)
ok = cli.run_pipeline(
    tsv_path=Path("data/Sphyrna_lewini.tsv"),
    organism="Sphyrna_lewini",
    output_dir=out,
    cfg=cfg,
)
```

`run_pipeline(tsv_path, organism, output_dir, cfg, no_report=False, skip_geo=False,
export_plot_data=False, export_popgen_formats=None, metadata_summary_enabled=False,
metadata_fields=None, normalize_sex=True, temporal_analysis=False, ...,
taxon_curation=None)` returns `True`
when the run finishes (`completed` or `completed_with_warnings` in
`run_provenance.json`). `export_plot_data` is `--plot-kit`; `export_popgen_formats` is
`--export-format`; `taxon_curation` is `--taxon-curation` (a curation file that cannot be
read or lacks required columns makes it return `False` before writing anything). Call
`utils.setup_logging(log_file=...)` first to get the run log.

---

## Reading a Run Folder

```python
import pandas as pd
from pathlib import Path

run = Path("results/Sphyrna_lewini")
summary = pd.read_csv(run / "summary.csv").iloc[0]
records = pd.read_csv(run / "records.csv", low_memory=False)   # one row per record
haps = pd.read_csv(run / "haplotypes.csv")                     # one row per haplotype

members = records[records["status"].isin(["assigned", "assigned_corrected"])]
by_region = pd.read_csv(run / "geography" / "haplotype_by_region.csv")
trends = pd.read_csv(run / "extended" / "metadata" / "temporal_trends.csv", comment="#")
```

---

## Configuration

```python
from boldgenotyper import config

cfg = config.get_default_config()
print(cfg.haplotype.min_singleton_distance)       # 0.005
print(cfg.taxonomy.majority_species_threshold)    # 0.7

# update() returns a new config; nested fields use section__field
cfg = cfg.update(haplotype__min_singleton_distance=0.015,
                 haplotype__assign_corrected=False,
                 coi_validation__mitochondrial_code=5,
                 extended=True)
```

Main sections and fields (all are recorded in `run_provenance.json`):

| Section | Fields |
|---|---|
| `qc` | `min_raw_length_abs` (200), `min_raw_length_frac_of_median` (0.7), `max_raw_N_fraction` (0.05), `require_valid_orf` (True) |
| `coi_validation` | `mitochondrial_code` (2), `orf_min_coverage` (0.5), `orf_max_internal_stops` (2) |
| `core_region` | `core_min_coverage` (0.8), `min_sequence_coverage` (0.95), `core_min_length` (150), `mask_gap_threshold` (0.5) |
| `haplotype` | `min_singleton_distance` (0.005; `--singleton-distance`), `max_singleton_distance` (0.05), `assign_corrected` (True) |
| `taxonomy` | `majority_species_threshold` (0.7) |
| `phylogenetic` | `build_tree` (False), `min_cluster_size` (1; `--tree-min-members`), `min_consensus_length` (150), `trim_alignment`, `trim_method`, outgroup options |
| `msa` | `enabled`, `chunk_size`, `max_sequences`, `output_formats` (`["pdf"]`) |
| top level | `n_threads`, `output_dir`, `keep_intermediates`, `extended` |

---

## Module Reference

Signatures are abbreviated; see each function's docstring.

### Input and coordinates — `metadata`

- `parse_bold_tsv(tsv_path, required_columns=None)` → DataFrame; duplicate processids keep the first record
- `parse_coord_column(df, coord_col="coord")` → adds float `lat`, `lon` for every parseable coordinate pair (BOLD portal `[lat, lon]` strings included)
- `mark_coordinate_quality(df)` → adds `coord_status` (`ok`, `missing`, `unparseable`, `centroid`, `zero`, `out_of_range`) and `is_geographic_quality`; removes nothing

### Quality control and labels — `quality_control`

- `apply_orientation_normalization(sequences_dict, genetic_code=2, ...)`, `apply_dynamic_qc_filters(df, sequences_dict, ...)` → oriented sequences and QC flags (`qc_pass`, `qc_fail_reason`)
- `add_label_columns(df, group_col="haplotype_id", notes_col="notes", curated_taxa=None)` → `label`, `majority_label`, `n_majority`, `n_labeled`, `n_missing_label`, `majority_fraction`, `label_conflict`, `depositor_keyword_hits`
- `summarize_haplotype_labels(df)` → one row per haplotype (labels with counts, majority, conflicts)
- `build_label_conflicts(df, haplotype_sequences=None)` → the rows of `label_conflicts.csv` (`haplotype_taxon` and `assignment_source` from `assigned_sp` and `assignment_source` when `df` has them)
- `find_keyword_hits(note)` → keywords found as whole words

### Haplotypes — `dereplication`, `haplotype_ids`, `record_status`

- `dereplication.identify_haplotypes(tsv_path, fasta_path, output_dir, min_core_coverage=0.8, min_singleton_distance=0.005, max_singleton_distance=0.05, min_sequence_coverage=0.95, assign_corrected=True, return_details=False, ...)` → `(haplotype_records, mapping, stats, [details])`; IDs are `H0001`..., FASTA headers carry `seq_md5` and `n_members`; `details` has the core length, cutoffs, excluded records and error-singleton outcomes
- `dereplication.flag_suspect_haplotypes(haplotype_records, mapping, max_singleton_distance=0.05)` → per-haplotype stats with `nearest_haplotype`, `nearest_neighbor_distance`, `is_suspect`, `suspect_reasons`
- `haplotype_ids.seq_md5(seq)` → md5 of the uppercase, ungapped sequence
- `haplotype_ids.assign_display_ids(counts)` → `{seq_md5: "H0001", ...}` by count, then md5
- `record_status.build_record_status(unique_ids, sequence_ids, qc, exact_mapping, excluded_core, error_singletons, singleton_distance)` → one status per unique record; `record_status.STATUSES`

### Taxonomy — `taxon_labels`, `taxon_evidence`, `taxon_curation`, `utils`

- `taxon_labels.label_level(label)` → `species`, `open`, `genus` or `missing`
- `taxon_labels.pool_labels(labels)` → `PooledLabels` (`text`, `leading_species`, `support`, `n_leading`, `n_species_level`, `n_compatible`, `n_missing`)
- `taxon_labels.species_count_summary(input_labels, qc_passed_labels, haplotype_taxonomy, majority_threshold)` → the `n_species_*` counts of `summary.csv`
- `utils.assign_consensus_taxonomy(df, group_col, species_col="species", genus_col="genus", majority_threshold)` → `(assign_table, species_counts)`; `assigned_sp` is the majority species when `majority_fraction` (members without a species label left out) exceeds the threshold, otherwise the genus
- `taxon_evidence.compute_evidence(taxonomy, records, sequences, tree_path=None)` → `haplotype_id` + `taxon_evidence.EVIDENCE_COLUMNS` (nearest species-level haplotype and spread, BINs, label sources, tree sister, `suggested_sp`, `suggestion_basis`); `taxon_evidence.label_source(processid, inst, insdc_acs)` → project code or GenBank accession prefix; `BIN_MAJORITY` (0.70)
- `taxon_curation.load_curation(path)` → validated curation rows (raises `CurationError`); `match_curation(curation, haplotypes, sequences)` → `(applied, not_applied)`; `apply_curation(taxonomy, applied)` → taxonomy with `assignment_source`; `curated_taxa(taxonomy)` → `{haplotype_id: taxon}` of curated haplotypes

### Run folder — `layout`, `run_readme`, `provenance`

- `layout.build_records_table(records, status, qc, taxonomy, region_col, region_confidence=None)` → `records.csv`
- `layout.build_haplotypes_table(stats, taxonomy, records, region_col, tree_tips=None, evidence=None)` → `haplotypes.csv`
- `layout.build_haplotype_by_region(records, region_col)` → `geography/haplotype_by_region.csv`
- `layout.plan_moves(staging, organism, extended=False)`, `layout.apply_moves(moves, output_dir)`, `layout.FILE_MAP`, `layout.KEY_BOLD_FIELDS`
- `run_readme.write_run_readme(output_dir, organism, region_col="ocean_basin")` → `(path, undocumented)`
- `provenance.RunProvenance`, `provenance.find_previous_run(output_dir)`

### Geography — `geographic`, `geographic_enhancement`

- `geographic.load_goas_data(shapefile_path)` → GeoDataFrame of GOaS basins
- `geographic.assign_ocean_basins(df, goas_data=None, coord_col=None, ...)` → adds `ocean_basin`
- `geographic.assign_regions_from_shapefile(df, shapefile_path, shapefile_field=None, output_column="geographic_region", ...)`
- `geographic_enhancement.enhance_geographic_analysis(df, output_dir, organism, ...)` → coverage counts and `basin_confidence`

Pass only records with `coord_status == "ok"` to the region functions, as the pipeline does.

### Phylogeny — `phylogenetics`

- `select_tree_haplotypes(records, min_members, min_length, suspect)` → `(kept_records, [(haplotype_id, reason), ...])`
- `build_phylogeny(consensus_fasta, output_prefix, threads=1, min_consensus_length=..., min_cluster_size=..., trim_alignment=True, trim_method="gappyout", outgroup_fasta=None)` → MAFFT, optional trimAl, FastTree (GTR+Gamma)
- `midpoint_root_tree(tree)`, `root_on_tips(tree, tip_names)` → rooted in place with support values kept on their splits (Biopython's own rerooting moves them); `root_on_tips` also returns the number of other tips in the outgroup's clade (0 when the outgroup is a clade)
- `midpoint_root_tree_file(tree_file)`, `reroot_tree_by_tips(tree_file, tip_names)`, `reroot_tree_by_label(tree_file, label)`, `reroot_tree_by_taxon(tree_file, taxonomy_csv, taxon)` → the same on Newick files

### Divergence — `divergence_analysis`

- `generate_divergence_analysis(consensus_fasta, taxonomy_csv, output_dir, organism, extended=False, records=None)` → aligned haplotypes and pairwise p-distances; with `extended`, the label-based gap
- `compute_label_based_gap(distance_by_id, records)` → per-label `max_within`, `min_between`, `label_based_gap`

### Metadata — `metadata_summary`

- `normalize_vocabulary(df, fields, map_sex_codes=True)` → normalized values (blank → missing, codes expanded, GIS suffixes removed)
- `run_metadata_summary(annotated_df, output_dir, organism, fields=None, normalize_sex=True, temporal_analysis=True, extended=False, cross_tabs=True, ...)` → coverage, and with `extended` the cross-tabs and collection-date summaries

### Figures — `visualization`, `plot_export`, `plot_regeneration`

- `run_figures.write_run_figures(records, run_dir, tree_path=None, haplotypes=None, region_col=None, geography=True, metadata_fields=None, extended=False, ...)` → `(manifests, warnings)`: every figure of a run, by species, from `records.csv` and the tree
- `tree_figures.write_tree_figures(tree_path, haplotypes, out_dir)`; `run_figures.write_geography_figures(records, region_col, geography_dir, extended_dir=None)`; `metadata_figures.write_metadata_figures(records, fields, region_col, out_dir)`; `metadata_figures.region_facet(records, region_col)` → `(column, caption)`
- `figure_pages.cap_items(counts, max_items=12)`, `species_groups(haplotypes)`, `haplotype_colors(haplotypes)`, `write_paged_figure(name, out_dir, overview, page, groups)`, `write_manifest(out_dir, entries)`
- `visualization.get_genotype_colors(n)` and the map helpers the figure modules share
- `plot_export.export_plots_complete(df, output_dir, organism, color_map=None)` → the `plots/` kit
- `plot_regeneration.regenerate_all_plots(base_dir)` → regenerate from a kit

PNG output is capped at 150 dpi and 3000 px per side (`utils.save_capped_png`).

### Report — `report`, `reports`

- `report.generate_html_report(run_dir, organism, provenance=None, version=None)` → `report.html` from the final run folder
- `report.data.load_run(run_dir, organism, provenance=None)` → `RunData` (tables, figure manifests, provenance); `report.insights.key_findings(run)` → `[Finding(key, text, count, link), ...]`
- `report.page`: `table(df, columns, max_rows=10, file=None)` (ten rows and the true count of the rest), `switcher`, `two_level_switcher`, `bars`, `render_page`
- `reports.compute_assignment_summary(...)`, `reports.generate_assignment_summary(...)` → `summary.csv`

### Query — `haplotype_query`, `query_assignment`, `query_outputs`, `query_overlays`, `query_report`

- `haplotype_query.align_query_to_haplotype(query_seq, query_id, haplotype_seq, haplotype_id)` → `AlignmentResult` (`differences`, `n_mismatches`, `n_gaps`, `identity_pct`, `query_coverage`, `overlap_coverage`, ...)
- `haplotype_query.query_against_haplotypes(query_fasta, haplotype_fasta, analysis_dir=None, top_n=10, min_identity=0.97, min_query_coverage=0.95, max_mismatch_delta=0, min_species_support=0.70, haplotype_taxa=None)` → `(results, metadata, verdicts)`; `haplotype_taxa` defaults to the run folder's (`load_haplotype_taxa`)
- `query_assignment.assign_query(query_id, matches, member_labels, ..., haplotype_taxa=None)` → `QueryVerdict` (`species_call`, `call_basis`, `flags`, ...)
- `query_assignment.load_haplotype_taxa(run_dir)` → `{haplotype_id: HaplotypeTaxon(taxon, level, source)}` from `haplotypes.csv`
- `query_outputs.calls_table(verdicts, reference_run=None, record_matches=None)` → `query_calls.csv` rows (`CALL_COLUMNS`); `candidates_table(results, verdicts, haplotypes=None, haplotype_taxa=None, ...)` → `query_candidates.csv` rows (`CANDIDATE_COLUMNS`); `write_calls`, `write_candidates`, `write_alignments`, `write_json`, `build_query_provenance`, `write_query_provenance`, `write_query_readme(output_dir, provenance)` → `(path, undocumented)`, `methods_sentence(provenance)`
- `query_overlays.load_reference(run_dir)` → `ReferenceRun` (records, haplotypes, provenance, region column, tree path) or None; `load_query_metadata(path, query_ids)` (raises `QueryMetadataError`); `assign_query_regions(meta, reference)`; `query_regions_table(calls, meta, reference)`; `plot_query_tree`, `plot_query_map`, `plot_query_regions`, `plot_query_timeline` → written paths; `write_query_figures(calls, meta, reference, out_dir)` → `(figures, warnings)`
- `query_report.generate_query_report(output_dir, calls, candidates, provenance, figures, alignments)` → `query_report.html`
- `query_assignment.resolve_reference(haplotypes_arg, analysis_dir_arg)` → `(fasta_path, run_dir)` for a v2 or v1 run folder or a FASTA
- `query_assignment.load_reference_run(path)` → the `reference_run` block of a v2 run folder

### Sweep — `parameter_sweep`

- `run_parameter_sweep(tsv_path, thresholds, output_dir, organism=None, threads=4, keep_intermediates=False, genetic_code=2)` → paths and the elbow threshold
- `extract_run_metrics(output_dir, organism, threshold)`, `track_group_membership(run_outputs, organism, thresholds)`, `adjacent_change_rates(tracking_df, thresholds)`

### Exports — `popgen_export`

- `export_population_genetics_formats(df, consensus_fasta_path, output_dir, organism, formats=None, pop_column="ocean_basin", haplotype_column="haplotype_id")` → `exports/`; records without a population are left out of Arlequin and DnaSP

---

## Examples

### Example 1: Query new sequences against a run

```python
from pathlib import Path
from boldgenotyper import haplotype_query, query_assignment

fasta, run = query_assignment.resolve_reference("results/Sphyrna_lewini", None)
results, metadata, verdicts = haplotype_query.query_against_haplotypes(
    query_fasta=Path("new_samples.fasta"), haplotype_fasta=fasta, analysis_dir=run)
for v in verdicts:
    print(v.query_id, v.call_level, v.equivalence_set, v.species_call, v.species_support)
```

### Example 2: Label summary of an existing run

```python
import pandas as pd

haps = pd.read_csv("results/Sphyrna_lewini/haplotypes.csv")
mixed = haps[haps["n_species_labels"] > 1]
print(mixed[["haplotype_id", "labels", "majority_label", "majority_fraction", "n_label_conflicts"]])
```

### Example 3: Batch processing

```python
import subprocess
from pathlib import Path

datasets = {"Sphyrnidae": 0.015, "Carcharhiniformes": 0.015}   # chosen with boldgenotyper-sweep
for organism, sd in datasets.items():
    subprocess.run(["boldgenotyper", f"data/{organism}.tsv", "--organism", organism,
                    "--output", f"results/{organism}", "--singleton-distance", str(sd),
                    "--build-tree", "--threads", "8"], check=True)
```

### Example 4: Regions from a custom column

```python
import pandas as pd
from boldgenotyper import popgen_export

run = "results/Sphyrna_lewini"
records = pd.read_csv(f"{run}/records.csv")
assigned = records[records["status"].isin(["assigned", "assigned_corrected"])]
popgen_export.export_population_genetics_formats(
    df=assigned, consensus_fasta_path=f"{run}/haplotypes.fasta",
    output_dir="by_country", organism="Sphyrna_lewini",
    formats=["arlequin"], pop_column="country/ocean")
```

### Comparing runs

To compare a species-level and a family-level run, match their `records.csv` files by
`processid` with pandas; see [COMPARATIVE_ANALYSIS_GUIDE.md](COMPARATIVE_ANALYSIS_GUIDE.md).
