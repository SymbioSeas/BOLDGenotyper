# Query Guide for BOLDGenotyper

**Author**: Steph Smith
**BOLDGenotyper version**: 2.3.0

---

## Table of Contents

- [Overview](#overview)
- [When to Query](#when-to-query)
- [Quick Start](#quick-start)
- [Choosing the Reference](#choosing-the-reference)
- [How a Query Is Matched](#how-a-query-is-matched)
- [The Species Call](#the-species-call)
- [Output Files](#output-files)
- [Figures and Query Metadata](#figures-and-query-metadata)
- [Curated Haplotype Taxa](#curated-haplotype-taxa)
- [Reporting a Query](#reporting-a-query)
- [Repeat Queries](#repeat-queries)
- [Troubleshooting](#troubleshooting)

---

## Overview

`boldgenotyper-query` places new COI sequences on the haplotypes of a completed
BOLDGenotyper run, without re-running the pipeline: new specimens, field samples, museum
material, or sequences from another study. For each query it reports which haplotypes the
query matches, which of them it cannot tell apart, whether their member records support a
species call, and the numbers behind the call. It also draws the queries on the reference
run's tree, maps, region counts and sampling years.

A query run takes seconds and writes to its own folder; the reference run is never changed.

---

## When to Query

- **Query** when you want to know where new sequences fall among the haplotypes of an
  existing run, and to report that against a fixed, identified reference.
- **Re-run the pipeline** when the new sequences should change the haplotypes themselves
  (for example, you are adding many records to BOLD and want them counted in haplotype
  frequencies, regions or the tree). Queries never add haplotypes.

---

## Quick Start

```bash
boldgenotyper-query \
  --query new_samples.fasta \
  --haplotypes results/Sphyrna_analysis \
  --query-metadata samples.csv \
  --output query_results/
```

- `--query`: one consensus sequence per specimen in a (multi-)FASTA, not raw forward and
  reverse reads. Gaps (`-`, `.`) and lowercase letters are cleaned automatically; a query
  that looks aligned is logged.
- `--haplotypes`: the reference run folder (recommended), or a haplotype FASTA.
- `--query-metadata`: optional; adds each query's location, region and date to the figures.

Open `query_results/README.md` first: it says which file to report and gives a methods
paragraph with this query's values.

---

## Choosing the Reference

| `--haplotypes` | Haplotype matching | Species calls | Figures |
|---|---|---|---|
| A 2.0.0 run folder (`haplotypes.fasta` at its root) | yes | yes, from the run's record labels and haplotype taxa | yes |
| A 1.0.0 run folder (`haplotypes/*_haplotypes.fasta`) | yes | no (flag `no_reference_labels`) | no |
| A haplotype FASTA | yes | no (flag `no_reference_labels`) | no |

Use the run folder whenever you have it. Every output records the reference run by the md5
of its `run_provenance.json`, so a query result can always be traced to the exact run.

---

## How a Query Is Matched

Each query is aligned to every haplotype over their whole overlap (global alignment with free
end gaps; match +1, mismatch −1, gap opening −6, gap extension −1). Differences at the
query's ends therefore count.

- `differences` = mismatches + gap columns inside the overlap.
- A match is **eligible** when it covers at least `--min-query-coverage` (0.95) of the
  shorter of query and haplotype, compares at least `--min-overlap-length` (60) bases, and
  reaches `--assign-identity` (0.97).
- Query bases beyond a haplotype's own ends (`query_overhang`) cannot be compared and are
  left out of the coverage (`compared_coverage`, 2.1.0). Haplotypes are defined on a core
  region, which in some datasets starts tens of bases into the barcode; a 5' mini-barcode
  then hangs over the core's start and is judged on the part the core covers. Such a call
  carries the flag `compared_over_part`: fewer bases were compared than the query has, so
  more haplotypes may tie.
- The **equivalence set** is every eligible haplotype within `--max-mismatch-delta` (0)
  differences of the best match: the haplotypes this query cannot tell apart.
- `call_level` is `haplotype` (one haplotype in the set), `haplotype_set` (several) or
  `none` (no eligible match).

A short query often matches several haplotypes that differ only outside its region: that
is a `haplotype_set`, not an error. The flag `set_not_identical` says the set's haplotypes
differ within the query's region (possible with `--max-mismatch-delta` above 0).

---

## The Species Call

By default the call comes from the **pooled member labels** of the set's haplotypes:
every record assigned to them, including `assigned_corrected` records. The leading species
is called when it holds at least `--min-species-support` (0.70) of the counted labels
(genus-only labels of the leading genus are compatible). `species_support` is that share.
`call_basis` is then `pooled labels`.

When a haplotype in the set has a **curated or evidence taxon** (see below), the call
comes from the haplotypes' taxa instead: the taxon they all share is called, with
`call_basis` `haplotype taxa`; if their taxa differ, or one is not at species level, there
is no call and the flag `haplotype_taxa_disagree` is set. `species_support` still reports
the pooled labels, so you can see both.

No identity value is read as "same species" or "same genus"; the numbers are reported.

**Flags** (`;`-joined in `query_calls.csv`):

| Flag | Meaning |
|---|---|
| `indel_in_query` | a gap run in the alignment whose length is not a multiple of 3 (COI codes protein: check the query for a sequencing or editing error) |
| `low_coverage` | no match covers `--min-query-coverage` of the shorter sequence |
| `short_overlap` | no match with enough coverage compares `--min-overlap-length` bases |
| `below_identity` | no match with enough coverage reaches `--assign-identity` |
| `compared_over_part` | the query runs past the haplotype core: the best match covers neither the whole query nor the whole haplotype, so the call rests on the bases the core covers |
| `set_not_identical` | the set's haplotypes differ within the query's region |
| `label_conflict_in_set` | the pooled labels name more than one species |
| `species_support_below_threshold` | the leading species is below `--min-species-support` |
| `no_species_labels` | the set's members have no species-level labels |
| `no_reference_labels` | the reference is a FASTA or a 1.0.0 run folder |
| `haplotype_taxa_disagree` | a haplotype with a curated or evidence taxon is in the set, and the set's taxa differ or one is not at species level |

---

## Output Files

| File | What it holds |
|---|---|
| `query_calls.csv` | **the result to report**: one row per query (`call_level`, `equivalence_set`, `species_call`, `call_basis`, `species_support`, `pooled_member_labels`, `differences`, `identity_pct`, `query_coverage`, `flags`, `reference_run`) |
| `query_candidates.csv` | the evidence: the top `--top-n` haplotypes per query with their numbers, whether each is eligible and in the set, and the haplotype's taxon, `assignment_source` and `n_members` |
| `query_alignments.txt` | each query's call and its best alignment |
| `query_report.html` | the calls, one evidence card per query, the figures and the methods |
| `query_figures/` | the queries placed in the reference run (see below) |
| `query_provenance.json` | the reference run (md5, version, path), query and reference files with md5, parameters, curation and metadata files, figures, warnings, files written |
| `README.md` | which file to report, a methods paragraph with this run's values, and every column |
| `query_record_matches.csv` | `--search-records`: QC-passed reference records that contain the query unchanged (including records not assigned to a haplotype) |
| `query_results.json` | `--json`: everything as one JSON document |

The CSVs have no comment lines, so spreadsheets open them directly. Other haplotype
columns (labels, evidence, region counts) are in the reference run's `haplotypes.csv`: join
on `haplotype_id`.

---

## Figures and Query Metadata

For a 2.0.0 run-folder reference, `query_figures/` holds:

- `query_tree`: the reference tree (rooted as in the reference run) pruned to each called
  genus, with each query written beside the tips of its equivalence set; queries without an
  eligible match are listed under the figure.
- `query_map`: one page per called taxon: the reference run's records of that taxon, the
  records of the queries' sets colored by haplotype, and the queries as stars.
- `query_regions` and `query_regions.csv`: reference records of each set's haplotypes per
  region, with the queries in each region.
- `query_timeline`: the sampling years of each set's haplotypes, with the query dates.

Each figure is a multi-page PDF with one SVG per page. The tree needs no metadata; the
other figures show the queries' own locations and dates when you pass `--query-metadata`:

```csv
query_id,coord,collection_date,region
sample_01,"[27.5, -82.9]",2024-05-14,
sample_02,"[29.1, -89.6]",2023-09-02,
sample_03,,2024-06-30,Gulf of Mexico
```

- `query_id` must match the FASTA IDs (the first word of each header). Rows for other IDs
  are ignored and listed in `query_provenance.json`.
- Coordinates: `lat` and `lon` columns, or `coord` as `[lat, lon]`. They are checked like
  the pipeline's (missing, zero, out of range).
- A query's **region** is its `region` value; otherwise its coordinates are placed with the
  reference run's own shapefile and field (GOaS ocean basins or the custom shapefile the run
  used), when that shapefile is still at the path the reference run recorded.
- CSV, or TSV for `.tsv`/`.txt`.

---

## Curated Haplotype Taxa

The reference run's evidence and curated taxa are used automatically (they are in its
`haplotypes.csv`, column `assignment_source`: `evidence` when the run's evidence named a
haplotype's leading label, `curated` from a `--taxon-curation` file). The methods
paragraph names how many taxa the reference run set by evidence, and the run's curation
file with its md5.

To try or apply a curation for one query run only:

```bash
boldgenotyper-query --query new_samples.fasta --haplotypes results/Sphyrna_analysis \
  --taxon-curation curation.csv --output query_results/
```

The file format and matching rules are those of the pipeline (README, "Reviewing and
Curating Haplotype Taxa"); rows match by `seq_md5` or by sequence. `query_provenance.json`
records the file's path, md5, and every row that was not applied, with the reason.
Query-time curation needs a run-folder reference.

---

## Reporting a Query

- Report `call_level`, `equivalence_set` and `species_call` with its `call_basis` from
  `query_calls.csv`, and cite the reference run by `reference_run` (the md5 of its
  `run_provenance.json`) and the BOLDGenotyper version.
- The methods paragraph in the query folder's `README.md` states this run's alignment
  scores, thresholds, equivalence delta, species support, and any curation file with its
  md5. Copy and edit it.
- Keep `query_provenance.json` with your results.

---

## Repeat Queries

Keep one reference run per reference dataset, and query the full set of your samples each
time you add some: a query run takes seconds for a handful of samples, and the calls of
earlier samples are recomputed with the same settings. Every output records the reference run's md5, so results
from different dates can be compared; a different md5 means a different reference run
(new download, settings or curation), and earlier calls should be repeated against it.

---

## Troubleshooting

**No eligible match (`call_level` none)**: check `query_coverage` and `identity_pct` of the
best candidate in `query_candidates.csv`. A short or untrimmed query may not cover 95% of
the shorter sequence (`low_coverage`); a query that overlaps the haplotype core over fewer
than `--min-overlap-length` bases is `short_overlap`; a sequence from a taxon absent from the
reference stays below `--assign-identity` (`below_identity`). `query_overhang` in
`query_candidates.csv` shows how many of the query's bases lie beyond the core.

**Queries take long**: each query is aligned to every haplotype (about 4 ms per pair), so
hundreds of queries against a family-level reference take an hour or more on one core. Pass
`--threads N` (2.3.0) to align queries in N worker processes; the outputs are the same for
any N, and `query_provenance.json` records the number.

**Several haplotypes in the set**: the query does not cover the sites where they differ.
The pooled labels then decide the species call; a longer query may separate them.

**`indel_in_query`**: look at the alignment in `query_alignments.txt`. A 1- or 2-bp gap in
a coding sequence usually means a base-calling or editing error in the query.

**No species call with a FASTA reference**: species calls need the reference run folder
(record labels); pass the folder to `--haplotypes`.

**Query regions are empty**: the reference run had no geographic analysis, the query rows
have no geographic-quality coordinates or `region`, or the reference's shapefile has moved
(the log says which).
