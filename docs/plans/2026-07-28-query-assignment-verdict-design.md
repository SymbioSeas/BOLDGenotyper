# Query assignment verdict, species interpretation, and reference-by-run-directory — design

Date: 2026-07-28
Status: Approved, pending implementation plan

## Problem

`boldgenotyper-query` reports the top-N haplotype matches per query, but leaves
interpretation to the user. When several haplotypes tie for the best identity
(as in the real Carcharhiniformes example, where h74_n18, h110_n13, and h249_n5
all match a query at 99.61%), the summary silently reports only the first, hiding
the tie. Users cannot tell whether a query is confidently one haplotype, tied
among several, or simply not a good match for anything — nor whether it resolves
to a single species.

Separately, users must pass the exact `_haplotypes.fasta` path as the reference,
even when they have the whole previous run directory on hand.

## Goals

- Give each query a clear two-level verdict: a haplotype-level call
  (confident / tied / no confident match) and a species-level call
  (confident / ambiguous / unassigned), computed from tunable thresholds.
- Resolve the species even when the haplotype is tied, as long as the tied
  haplotypes share one species.
- Let `--haplotypes` accept a previous run directory and auto-find the reference
  FASTA (and the taxonomy metadata) from it.
- Keep the existing per-match output; add interpretation on top, not in place of.

## Non-goals

- Phylogenetic placement of queries onto a tree (Piece C) — deferred to a
  separate later exploration.
- Changing the alignment/identity computation itself.

## Design

### Verdict computation

`assign_query` receives the FULL list of matches against every haplotype, sorted
by identity descending — not the truncated top-N (that truncation applies only to
the `query_results.csv` display, so a tie member ranked beyond `--top-n` is never
missed by the verdict). `identity_pct` is on a 0-100 scale; thresholds are
fractions converted with x100.

Given `min_identity` (default 0.97) and `tie_margin` (default 0.005):

1. If there are no matches, or `best.identity_pct < min_identity * 100`:
   - `haplotype_call = "no_confident_match"`, `assigned_haplotypes = []`
   - `species_call = "unassigned"`, `assigned_species = None`
2. Otherwise:
   - `tie_floor = max(best.identity_pct - tie_margin * 100, min_identity * 100)`
   - `tie_set = [m for m in matches if m.identity_pct >= tie_floor]`
   - `haplotype_call = "confident"` if `len(tie_set) == 1` else `"tied"`
   - `assigned_haplotypes = [m.haplotype_id for m in tie_set]`

Species level (from the tie set's `assigned_sp` in the metadata, if available):

- `distinct = ` sorted unique non-null `assigned_sp` over the tie set.
- If `haplotype_call == "no_confident_match"`: `species_call = "unassigned"`.
- Else if metadata is unavailable or the tie set has no species values:
  `species_call = "unknown"`, `assigned_species = None`,
  `candidate_species = []` (query still works; species layer degrades gracefully).
- Else if `len(distinct) == 1`: `species_call = "confident"`,
  `assigned_species = distinct[0]`.
- Else: `species_call = "ambiguous"`, `candidate_species = distinct`.

`candidate_species` is empty except in the ambiguous case; `assigned_species` is
set only in the confident case (None otherwise).

`reason` is a short human-readable string, e.g. "3 haplotypes tied at 99.61%,
all Sphyrna lewini" or "best match 94.2% below 97.0% threshold".

### QueryVerdict structure

```
@dataclass(frozen=True)
class QueryVerdict:
    query_id: str
    haplotype_call: str            # confident | tied | no_confident_match
    assigned_haplotypes: List[str]
    best_identity: float           # 0-100, or 0.0 when no matches
    n_tied: int
    species_call: str              # confident | ambiguous | unassigned | unknown
    assigned_species: Optional[str]
    candidate_species: List[str]
    reason: str
```

### Output

- New `query_summary.csv`, one row per query:
  `query_id, haplotype_call, assigned_haplotypes, best_identity, n_tied,
  species_call, assigned_species, candidate_species, reason,
  assign_identity, tie_margin`
  (list columns joined with `;`; the two threshold columns record what was used).
- Keep `query_results.csv` (per-match top-N) and `query_results.json` unchanged.
- Enhance the text report's per-query summary line to state the verdict, e.g.:

  ```
  17C_001_NOAA-5317 : TIED among 3 haplotypes (h74_n18, h110_n13, h249_n5) at 99.61%
                      all Sphyrna lewini -> SPECIES: confident (Sphyrna lewini)
  ```

### Reference-by-run-directory

`--haplotypes` accepts either a FASTA file (current behavior, unchanged) or a
previous run directory:

- If the value is a directory, resolve the reference by globbing
  `<dir>/haplotypes/*_haplotypes.fasta`. Exactly one match: use it. Zero:
  error "No *_haplotypes.fasta found under <dir>/haplotypes/". More than one:
  error listing the candidates and asking the user to pass the file directly.
- If the value is a directory and `--analysis-dir` was not given, set
  `analysis_dir` to that directory so taxonomy/metadata enrichment (which the
  species verdict needs) happens automatically from the same path.

### Code structure

- New focused module `boldgenotyper/query_assignment.py`:
  `assign_query(query_id, matches, metadata, min_identity, tie_margin) ->
  QueryVerdict`. Pure and independently unit-testable, kept out of the already
  large `haplotype_query.py`.
- `haplotype_query.py` calls `assign_query` per query, writes `query_summary.csv`,
  and enriches the text report summary lines.
- `cli.py` `main_query`: add `--assign-identity` (default 0.97) and `--tie-margin`
  (default 0.005); resolve `--haplotypes` as file-or-directory and default
  `--analysis-dir` from a run directory.

### CLI flags (new)

- `--assign-identity FLOAT` (default 0.97): minimum identity (fraction) to call a
  match assignable.
- `--tie-margin FLOAT` (default 0.005): haplotypes within this fraction of the
  best match are tied with it.

## Testing

Unit tests for `assign_query`:

- Confident single (one clear best above threshold).
- Tied, all one species -> haplotype tied, species confident.
- Tied across species -> species ambiguous with candidate list.
- Best below `min_identity` -> no confident match, species unassigned.
- Missing metadata / no species values -> species unknown, haplotype call intact.
- Threshold and tie-margin edges (e.g. a match exactly at the floor is included;
  `tie_floor` never dips below `min_identity`).

Reference resolution tests:

- FASTA file path unchanged.
- Directory with exactly one `*_haplotypes.fasta` resolves and sets `analysis_dir`.
- Directory with zero or several -> clear error.

## Documentation

Update the README `boldgenotyper-query` section: document the verdict columns and
their meaning, the two threshold flags, and the run-directory form of
`--haplotypes`.
