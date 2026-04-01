# Parameter Sweep Results — Panulirus

## What was swept

`min_singleton_distance` — the divergence cutoff below which singleton
haplotypes (single-member ESVs) are removed as likely sequencing or PCR
errors.  Higher values filter more aggressively; lower values retain more
rare variants.

## Recommended value

**0.015**  (see `recommendations.txt` for full rationale)

## Files

| File | Description |
|------|-------------|
| `sweep_summary.csv` | Per-threshold metrics table |
| `threshold_stability.pdf` | 4-panel stability visualisation |
| `elbow_plot.pdf` | Optimal-threshold detection |
| `group_membership_tracking.csv` | Sample-level haplotype stability |
| `recommendations.txt` | Automated interpretation and next steps |
| `runs/` | Per-threshold pipeline outputs |

## Reading group_membership_tracking.csv

Each row is one sample.  Columns `t_X.XXX` show which haplotype the sample
was assigned to at each tested threshold.  The key derived columns are:

- **n_changes** — number of times the sample's haplotype *composition*
  changed significantly (Jaccard similarity of cluster partners < 0.7).
  Low values (0–2) indicate stable assignment.
- **stability_score** — `high` (0 changes), `medium` (1–2), or `low` (>2).
  Samples scored `low` are worth investigating manually.
- **mean / min / max_cluster_size** — size of the sample's haplotype across
  the sweep.

> Haplotype *names* (e.g. `h1_n177`) are arbitrary and may differ between
> runs.  What matters is which other samples share the same haplotype.

## Interpretation

See `recommendations.txt`.  The recommended threshold is the elbow point
where haplotype count stops dropping sharply — balancing error removal
against retention of biologically real rare variants.

## Generated

2026-03-12 12:23:58
