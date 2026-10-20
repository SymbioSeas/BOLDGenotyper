# Parameter Sweep Guide for BOLDGenotyper

**Author**: Steph Smith
**BOLDGenotyper version**: 2.2.0

---

## Table of Contents

- [Overview](#overview)
- [When to Use Parameter Sweep](#when-to-use-parameter-sweep)
- [Quick Start](#quick-start)
- [Understanding the Outputs](#understanding-the-outputs)
- [Choosing a Threshold](#choosing-a-threshold)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Publication Methods Text](#publication-methods-text)

---

## Overview

`--singleton-distance` (`haplotype.min_singleton_distance`) is the main tunable setting of
BOLDGenotyper's ESV (exact sequence variant) workflow. Singleton haplotypes (one member)
within this p-distance of another haplotype are treated as likely sequencing or PCR errors
and removed. Their records are assigned to the nearest retained haplotype when it is
unique and within the same distance (status `assigned_corrected`); otherwise they stay
unassigned (`excluded_error_singleton`).

A suitable value depends on sequencing quality, how divergent real rare variants are in
your taxon, and your study's goals. `boldgenotyper-sweep` runs the pipeline at several
values and reports how the haplotypes change, so you can choose a value from your own data.

**Always pass the chosen value explicitly** with `boldgenotyper --singleton-distance`. It is
recorded in `run_provenance.json` either way, but an explicit value makes the command
self-documenting.

### What the Sweep Reports

- **Stability**: for each pair of adjacent values, the share of records whose haplotype
  changes (compared by sequence, so renumbered IDs do not count as changes)
- The **settled threshold**, the recommendation: the smallest value from which every later
  step changes the haplotype of fewer than 1% of records (2.2.0)
- The number of haplotypes and singletons at each value, with the settled threshold and
  the **elbow** of the haplotype-count curve marked (the elbow for reference only; see
  [Choosing a Threshold](#choosing-a-threshold))
- Assignment coverage, suspect haplotypes and label conflicts at each value

---

## When to Use Parameter Sweep

### Recommended Scenarios

1. **New taxonomic groups**: no prior experience with the group
2. **Mixed datasets**: several species or genera with different divergence levels
3. **Publication**: to justify the threshold with data
4. **Uncertain data quality**: unknown level of sequencing error

### When the Default May Do

- Earlier analyses on the same organism and data source used the default (0.005) and the
  sweep showed little change around it
- Exploratory analysis only

### Computational Considerations

The sweep runs the pipeline once per value (5 by default), without geography, the report or
the tree. Use `--threads`, or test fewer values first.

---

## Quick Start

### Step 1: Run the Sweep

Default values (0.01, 0.015, 0.02, 0.03, 0.05):

```bash
boldgenotyper-sweep data/Sphyrna_lewini.tsv --output parameter_sweep/
```

Custom values (evenly spaced, as the settled threshold compares steps of the same width):

```bash
boldgenotyper-sweep data/Sphyrna_lewini.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.025,0.03 \
  --threads 8 \
  --output sweep_results/
```

Pass the `--genetic-code` of the final run (default 2, vertebrate mitochondrial; 5 for most
invertebrates), so the sweep checks the same records:

```bash
boldgenotyper-sweep data/Panulirus.tsv --genetic-code 5 --output sweep_panulirus/
```

### Step 2: Review the Results

```bash
cat parameter_sweep/recommendations.txt     # the settled threshold, stability, haplotype counts
head parameter_sweep/sweep_summary.csv      # all metrics per value
open parameter_sweep/elbow_plot.pdf
open parameter_sweep/threshold_stability.pdf
```

### Step 3: Run the Analysis with the Chosen Value

```bash
boldgenotyper data/Sphyrna_lewini.tsv \
  --singleton-distance 0.015 \
  --build-tree \
  --output final_analysis/
```

---

## Understanding the Outputs

### Directory Structure

```
parameter_sweep/
├── sweep_summary.csv              # metrics for each value tested
├── threshold_stability.pdf        # four panels: metrics against value
├── elbow_plot.pdf                 # haplotype count against value, settled threshold and elbow marked
├── group_membership_tracking.csv  # each record's haplotype at each value
├── recommendations.txt            # the settled threshold, haplotype-count and stability tables
├── README.md                      # interpretation guide
└── runs/                          # record and haplotype tables per threshold (full run folders with --keep-intermediates)
    ├── threshold_0_010/
    └── ...
```

### `sweep_summary.csv`

One row per value:

| Column | Meaning |
|---|---|
| `threshold` | the `--singleton-distance` value |
| `n_haplotypes` | retained haplotypes |
| `n_singletons` | singleton haplotypes retained |
| `n_suspect` | haplotypes flagged as suspect |
| `n_samples` | QC-passed records |
| `pct_assigned` | QC-passed records assigned to a haplotype (assigned or assigned_corrected) |
| `n_species_labels_qc_passed` | distinct species labels among QC-passed records (the same at every value) |
| `pct_label_conflict` | records with a label conflict |
| `n_haplotypes_with_label_conflict` | haplotypes with at least one label conflict |
| `mean_majority_fraction` | mean `majority_fraction` over haplotypes |

### `threshold_stability.pdf`

- **A. Haplotype count**: falls steeply, then levels off
- **B. Singleton count**: falls as the value rises
- **C. Suspect count**: haplotypes flagged as suspect
- **D. Assignment quality**: percentage assigned, and records with a label conflict

### `elbow_plot.pdf`

Haplotype count against value, with the settled threshold (green, dashed) and the elbow
(gray, dotted; the largest change in the slope of the curve) marked.

### `group_membership_tracking.csv`

One row per record:

| Column | Meaning |
|---|---|
| `t_<value>` | the record's haplotype ID at that value (`unassigned`, or `NA` when the record is not in that run) |
| `md5_<value>` | the `seq_md5` of that haplotype (empty when unassigned) |
| `n_changes` | how many times the record's co-members changed substantially between adjacent values (Jaccard similarity of co-member sets below 0.7) |
| `stability_score` | `high` (0 changes), `medium` (1-2), `low` (more than 2) |
| `mean_cluster_size`, `min_cluster_size`, `max_cluster_size` | size of the record's haplotype across values |

Haplotype IDs (`H0001`, ...) are renumbered in every run by member count, so the same ID can
mean different haplotypes at different values. Compare by `md5_<value>`.

### `recommendations.txt`

```
HAPLOTYPE COUNT
---------------
   threshold   haplotypes   singletons   assigned (% of QC-passed)
   0.010           45           20         87.3
   0.015           32            6         94.5   <- settled, elbow
   0.020           30            4         95.0
   0.030           29            2         95.1
   0.050           27            0         95.3

The settled threshold is 0.015: from there on, every step to the next
threshold changes the haplotype of fewer than 1% of records (STABILITY below).
It compares neighboring thresholds only, so it does not depend on where the
range starts.

The elbow (largest change in the slope of the haplotype-count curve) is at
0.015: 32 haplotypes, against 45 at 0.010 and 27 at 0.050.
It is shown for reference only: it depends on the range tested, because on a
curve whose drops shrink from one step to the next, the slope changes most
near the start of the range.
Here it is the second threshold tested.

The thresholds are unevenly spaced: a wider step changes more records, so
compare the shares of steps of the same width.

STABILITY
---------
Share of records whose haplotype changes between adjacent thresholds
(compared by sequence; becoming assigned or unassigned counts as a change):

   0.010 -> 0.015       41 of 912 records   4.5%
   0.015 -> 0.020        6 of 912 records   0.7%
   0.020 -> 0.030        3 of 912 records   0.3%
   0.030 -> 0.050        5 of 912 records   0.5%
```

(Numbers are illustrative.) The file ends with the command for the next step, with the
settled threshold. It does not give an "acceptable range": the choice is yours.

---

## Choosing a Threshold

### 1. The Settled Threshold

The sweep recommends the smallest value from which every later step to the next value
changes the haplotype of fewer than 1% of records (the stability table). Below it, raising
the value still moves records between haplotypes; from it on, the haplotypes are nearly
the same for nearly every record, so the result does not hinge on the exact value.

- It compares neighboring values only, so it does not move when the range starts lower or
  higher. Sweep evenly spaced values: a step twice as wide changes more records.
- It cannot be lower than the first value tested. When the first value is already settled,
  `recommendations.txt` says so: sweep lower values to see where the haplotypes settle.
- When no value settles (the last step still changes 1% or more), extend the range upward.
- The 1% cutoff is `parameter_sweep.SETTLED_MAX_PCT_CHANGED` (and the `max_pct_changed`
  argument of `run_parameter_sweep`). The stability table lists every step, so another
  cutoff can be read from it.

### 2. Stability

The stability table shows how many records change haplotype between neighboring values. A
large share between two values shows where the choice matters; a share that rises again at
a higher step (two neighbors merging) moves the settled threshold above that step.

### 3. The Elbow (for Reference)

The elbow is the largest change in the slope of the haplotype-count curve. It is reported,
but it is not the recommendation: these curves fall steeply and then level off, each step
removing fewer haplotypes than the one before, so the slope changes most near the start of
whatever range is tested. On six BOLD datasets, sweeps starting at 0.000, 0.005 and 0.010
put the elbow at 0.005, 0.010 and 0.015 in five of the six: the second value tested.
`recommendations.txt` says when the elbow is the second value tested.

Locating an elbow needs at least 4 values: with 3 there is only one change in slope, so the
middle value would always be chosen.

### 4. Assignment and Labels

`pct_assigned` rises as more error-filtered singletons are corrected into retained
haplotypes. `pct_label_conflict` rising with the value can mean that singletons carrying
other labels are being merged into their neighbors; check `label_conflicts.csv` in a run
at that value (`--keep-intermediates`).

### 5. Biology

Check the haplotypes against what you know of the taxon: too many singletons suggests
errors are being kept; very few haplotypes in a diverse dataset suggests over-filtering.

### Lenient vs. Aggressive

**Lenient (lower value, more singletons retained)**: high-quality data, rare-variant
discovery; may keep sequencing errors as haplotypes. Example: 0.005-0.010.

**Aggressive (higher value, fewer singletons retained)**: noisier data; may remove real rare
haplotypes. Example: 0.03-0.05.

### Ambiguous Results

**No threshold settles** (every step changes 1% or more of records): extend the range
upward, or check data quality.

**The share dips under 1% and rises again**: the settled threshold is above the later
rise. Read the records that change at that step (`group_membership_tracking.csv`); if
the lower value suits your question, report both in a sensitivity analysis.

---

## Advanced Usage

### Custom Ranges

```bash
# High-quality data, rare variants
boldgenotyper-sweep data/population_samples.tsv \
  --thresholds 0.001,0.003,0.005,0.007,0.01,0.015 \
  --output fine_sweep/

# Noisier data or divergent taxa
boldgenotyper-sweep data/diverse_taxa.tsv \
  --thresholds 0.03,0.05,0.07,0.10,0.15 \
  --output coarse_sweep/
```

### Keeping Full Run Folders

```bash
boldgenotyper-sweep data/samples.tsv --keep-intermediates --output sweep_with_runs/
ls sweep_with_runs/runs/threshold_0_020/     # a v2 run folder (records.csv, haplotypes.csv, ...)
```

Use them to compare haplotypes, label conflicts or records at different values.

### Records That Change

```python
import pandas as pd

t = pd.read_csv("parameter_sweep/group_membership_tracking.csv")
changed = t[t["md5_0.015"] != t["md5_0.02"]]       # haplotype changed between 0.015 and 0.02
print(changed[["processid", "t_0.015", "t_0.02"]])
```

---

## Troubleshooting

### All values give similar results

The dataset may be homogeneous, or the range too narrow. Try a broader range:

```bash
boldgenotyper-sweep data/samples.tsv --thresholds 0.005,0.02,0.05,0.10,0.15 --output broad_sweep/
```

### Many singletons at every value

Check sequence quality and the record statuses in a single run (`summary.csv`,
`records.csv` `qc_fail_reason`), and the genetic code for invertebrates (`--genetic-code 5`).

### The sweep takes too long

Use fewer values, more `--threads`, or a pilot sweep on a subsample, then refine the range.

---

## Publication Methods Text

### Example

```
The singleton error-filtering threshold (--singleton-distance) was chosen with
boldgenotyper-sweep (BOLDGenotyper v2.2.0; Smith, 2026), which ran the pipeline at
0.005 to 0.030 in steps of 0.005. The haplotype count fell from 45 to 27 across this
range. From 0.015 on, every step to the next value changed the haplotype of fewer
than 1% of records (compared by sequence), so we used 0.015: singleton haplotypes
within a p-distance of 0.015 of another haplotype were treated as likely sequencing
errors, and their records assigned to the nearest retained haplotype when it was
unique and within that distance.
```

### Short Version

```
The singleton filtering threshold (0.015) was chosen with boldgenotyper-sweep as the
smallest value from which fewer than 1% of records changed haplotype at each further
step (Supplementary Figure S1).
```

### Figure Legend for Sweep Plots

```
Figure S1: Singleton error-filtering threshold. (A) Number of haplotypes at each
--singleton-distance value; the threshold used (0.015) is marked. (B) Singleton haplotypes
retained. (C) Suspect haplotypes. (D) Percentage of QC-passed records assigned to a
haplotype, and records with a label conflict.
```

(Numbers are illustrative; use your own.)

---

## Best Practices

1. Check data quality first (one ordinary run: `summary.csv`, `records.csv`)
2. Start with the default values, then refine the range
3. Read the plots and both tables; the settled threshold is a starting point, not a rule
4. Pass the chosen value explicitly, and keep `recommendations.txt` and the plots with
   your methods

---

## Additional Resources

- Main pipeline: [README.md](../../README.md)
- Migration from 1.0.0: [MIGRATION_v1_to_v2.md](MIGRATION_v1_to_v2.md)

**ESV approach**:
- Porter, T. M., & Hajibabaei, M. (2020). Putting the I in ITS: Improved characterization of ITS diversity for the eukaryotes using exact sequence variants. Ecology and Evolution, 10(14), 7110-7127.

**Elbow method** (reported for reference):
- Thorndike, R. L. (1953). Who belongs in the family? Psychometrika, 18(4), 267-276.

---

## Support

1. **Check documentation**: this guide and README.md
2. **GitHub Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
3. **Email**: symbioseas@outlook.com

When reporting issues, include the command, `sweep_summary.csv`, and the plots.
