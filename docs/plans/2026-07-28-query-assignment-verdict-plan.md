# Query Assignment Verdict Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give each `boldgenotyper-query` result a clear haplotype-level and species-level verdict, and let `--haplotypes` accept a previous run directory.

**Architecture:** A new pure module `boldgenotyper/query_assignment.py` computes a `QueryVerdict` from a query's full match list plus a species map, and resolves a reference FASTA from a file-or-directory argument. `haplotype_query.py` computes verdicts before truncating to top-N, writes `query_summary.csv`, and enriches the text report. `cli.py` adds the flags and reference resolution. Additive: existing outputs are unchanged.

**Tech Stack:** Python 3.9+, pandas (already used), pytest.

## Global Constraints

- Python floor: 3.9 (no 3.10+ syntax).
- Dependencies: add none; pandas and PyYAML already available.
- No emojis in code, logging, comments, or docs.
- `identity_pct` is on a 0-100 scale; the new threshold flags are fractions (0-1) and are multiplied by 100 for comparison.
- Additive: `query_results.csv` and `query_results.json` stay unchanged; the verdict is added alongside.

## File Structure

- Create `boldgenotyper/query_assignment.py` — `QueryVerdict`, `assign_query`, `resolve_reference`. Pure logic, no I/O beyond path globbing.
- Create `tests/test_query_assignment.py` — unit tests.
- Modify `boldgenotyper/haplotype_query.py` — compute/return verdicts; write `query_summary.csv`; enrich the text report.
- Modify `boldgenotyper/cli.py` — `main_query`: new flags, reference resolution, pass-through.
- Modify `README.md` — document the verdict, flags, and run-directory reference.

---

### Task 1: Verdict logic (`assign_query`)

**Files:**
- Create: `boldgenotyper/query_assignment.py`
- Test: `tests/test_query_assignment.py`

**Interfaces:**
- Produces: `QueryVerdict` (frozen dataclass with fields `query_id: str`, `haplotype_call: str`, `assigned_haplotypes: List[str]`, `best_identity: float`, `n_tied: int`, `species_call: str`, `assigned_species: Optional[str]`, `candidate_species: List[str]`, `reason: str`); `assign_query(query_id: str, matches: list, species_by_haplotype: Dict[str, Optional[str]], min_identity: float, tie_margin: float) -> QueryVerdict`. `matches` is the full list of match objects for one query, each with `.haplotype_id` and `.identity_pct` (0-100), sorted by identity descending. `haplotype_call` in {confident, tied, no_confident_match}; `species_call` in {confident, ambiguous, unassigned, unknown}.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_query_assignment.py
from collections import namedtuple
from boldgenotyper import query_assignment as qa

M = namedtuple("M", "haplotype_id identity_pct")


def _matches(*pairs):
    return [M(h, p) for h, p in pairs]


def test_confident_single_haplotype_and_species():
    matches = _matches(("h1", 99.6), ("h2", 92.0))
    v = qa.assign_query("q1", matches, {"h1": "Sphyrna lewini"}, 0.97, 0.005)
    assert v.haplotype_call == "confident"
    assert v.assigned_haplotypes == ["h1"]
    assert v.n_tied == 1
    assert v.species_call == "confident"
    assert v.assigned_species == "Sphyrna lewini"
    assert v.candidate_species == []


def test_tied_same_species_is_species_confident():
    matches = _matches(("h1", 99.6), ("h2", 99.6), ("h3", 99.6), ("h4", 90.0))
    sp = {"h1": "Sphyrna lewini", "h2": "Sphyrna lewini", "h3": "Sphyrna lewini"}
    v = qa.assign_query("q1", matches, sp, 0.97, 0.005)
    assert v.haplotype_call == "tied"
    assert v.assigned_haplotypes == ["h1", "h2", "h3"]
    assert v.n_tied == 3
    assert v.species_call == "confident"
    assert v.assigned_species == "Sphyrna lewini"


def test_tied_across_species_is_ambiguous():
    matches = _matches(("h1", 99.6), ("h2", 99.5))
    sp = {"h1": "Sphyrna lewini", "h2": "Sphyrna mokarran"}
    v = qa.assign_query("q1", matches, sp, 0.97, 0.005)
    assert v.haplotype_call == "tied"
    assert v.species_call == "ambiguous"
    assert v.assigned_species is None
    assert v.candidate_species == ["Sphyrna lewini", "Sphyrna mokarran"]


def test_no_confident_match_below_threshold():
    matches = _matches(("h1", 94.2), ("h2", 90.0))
    v = qa.assign_query("q1", matches, {"h1": "Sphyrna lewini"}, 0.97, 0.005)
    assert v.haplotype_call == "no_confident_match"
    assert v.assigned_haplotypes == []
    assert v.species_call == "unassigned"
    assert "94.2" in v.reason


def test_missing_species_metadata_is_unknown():
    matches = _matches(("h1", 99.6),)
    v = qa.assign_query("q1", matches, {}, 0.97, 0.005)
    assert v.haplotype_call == "confident"
    assert v.species_call == "unknown"
    assert v.assigned_species is None


def test_tie_floor_never_dips_below_min_identity():
    # best 97.1, tie_margin 0.5% -> naive floor 96.6, but min_identity 97 wins.
    matches = _matches(("h1", 97.1), ("h2", 96.8))
    v = qa.assign_query("q1", matches, {"h1": "A", "h2": "B"}, 0.97, 0.005)
    assert v.assigned_haplotypes == ["h1"]
    assert v.haplotype_call == "confident"


def test_no_matches_is_no_confident_match():
    v = qa.assign_query("q1", [], {}, 0.97, 0.005)
    assert v.haplotype_call == "no_confident_match"
    assert v.best_identity == 0.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'boldgenotyper.query_assignment'`

- [ ] **Step 3: Write minimal implementation**

```python
# boldgenotyper/query_assignment.py
"""Query assignment verdict: turn a query's matches into an interpretable call.

Given the full list of haplotype matches for one query and a species map, decide
whether the query is confidently one haplotype, tied among several, or not a
confident match to anything, and whether it resolves to a single species.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass(frozen=True)
class QueryVerdict:
    query_id: str
    haplotype_call: str
    assigned_haplotypes: List[str]
    best_identity: float
    n_tied: int
    species_call: str
    assigned_species: Optional[str]
    candidate_species: List[str]
    reason: str


def assign_query(query_id, matches, species_by_haplotype, min_identity, tie_margin):
    """Compute a QueryVerdict from a query's full, identity-sorted match list."""
    min_pct = min_identity * 100.0
    if not matches:
        return QueryVerdict(query_id, "no_confident_match", [], 0.0, 0,
                            "unassigned", None, [], "no haplotype matches")

    best = matches[0]
    if best.identity_pct < min_pct:
        return QueryVerdict(
            query_id, "no_confident_match", [], best.identity_pct, 0,
            "unassigned", None, [],
            f"best match {best.identity_pct:.2f}% below {min_pct:.1f}% threshold",
        )

    tie_floor = max(best.identity_pct - tie_margin * 100.0, min_pct)
    tie_set = [m for m in matches if m.identity_pct >= tie_floor]
    assigned_haplotypes = [m.haplotype_id for m in tie_set]
    haplotype_call = "confident" if len(tie_set) == 1 else "tied"

    tie_species = [species_by_haplotype.get(h) for h in assigned_haplotypes]
    distinct = sorted({s for s in tie_species if s})
    if not species_by_haplotype or not distinct:
        species_call, assigned_species, candidate_species = "unknown", None, []
        sp_reason = "species metadata unavailable"
    elif len(distinct) == 1:
        species_call, assigned_species, candidate_species = "confident", distinct[0], []
        sp_reason = f"all {distinct[0]}"
    else:
        species_call, assigned_species, candidate_species = "ambiguous", None, distinct
        sp_reason = f"spans {len(distinct)} species: {', '.join(distinct)}"

    if haplotype_call == "confident":
        hap_reason = f"single best match {best.identity_pct:.2f}%"
    else:
        hap_reason = f"{len(tie_set)} haplotypes tied at >= {tie_floor:.2f}%"

    return QueryVerdict(
        query_id, haplotype_call, assigned_haplotypes, best.identity_pct,
        len(tie_set), species_call, assigned_species, candidate_species,
        f"{hap_reason}; {sp_reason}",
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py -q`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/query_assignment.py tests/test_query_assignment.py
git commit -m "Add query assignment verdict logic"
```

---

### Task 2: Reference-by-run-directory resolution

**Files:**
- Modify: `boldgenotyper/query_assignment.py`
- Test: `tests/test_query_assignment.py`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `resolve_reference(haplotypes_arg, analysis_dir_arg) -> (Path, Optional[Path])`. If `haplotypes_arg` is a directory, returns the single `<dir>/haplotypes/*_haplotypes.fasta` and sets the analysis dir to that directory when `analysis_dir_arg` is None. If it is a file, returns it and `analysis_dir_arg` unchanged. Zero or several FASTAs raise `ValueError`.

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_query_assignment.py
import pytest
from pathlib import Path


def test_resolve_reference_file_unchanged(tmp_path):
    f = tmp_path / "X_haplotypes.fasta"
    f.write_text(">h1\nACGT\n")
    fasta, adir = qa.resolve_reference(f, None)
    assert fasta == f and adir is None


def test_resolve_reference_directory_finds_fasta_and_sets_analysis_dir(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    f = tmp_path / "haplotypes" / "Carcharhiniformes_haplotypes.fasta"
    f.write_text(">h1\nACGT\n")
    fasta, adir = qa.resolve_reference(tmp_path, None)
    assert fasta == f
    assert adir == tmp_path


def test_resolve_reference_directory_keeps_explicit_analysis_dir(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    (tmp_path / "haplotypes" / "X_haplotypes.fasta").write_text(">h\nA\n")
    other = tmp_path / "other"
    fasta, adir = qa.resolve_reference(tmp_path, other)
    assert adir == other


def test_resolve_reference_directory_zero_fastas_errors(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    with pytest.raises(ValueError):
        qa.resolve_reference(tmp_path, None)


def test_resolve_reference_directory_multiple_fastas_errors(tmp_path):
    (tmp_path / "haplotypes").mkdir()
    (tmp_path / "haplotypes" / "A_haplotypes.fasta").write_text(">h\nA\n")
    (tmp_path / "haplotypes" / "B_haplotypes.fasta").write_text(">h\nA\n")
    with pytest.raises(ValueError):
        qa.resolve_reference(tmp_path, None)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py -k resolve_reference -q`
Expected: FAIL with `AttributeError: module 'boldgenotyper.query_assignment' has no attribute 'resolve_reference'`

- [ ] **Step 3: Write minimal implementation**

```python
# add to boldgenotyper/query_assignment.py

def resolve_reference(haplotypes_arg, analysis_dir_arg):
    """Resolve --haplotypes (a FASTA file or a run directory) to
    (reference_fasta_path, analysis_dir)."""
    p = Path(haplotypes_arg)
    if p.is_dir():
        matches = sorted(p.glob("haplotypes/*_haplotypes.fasta"))
        if not matches:
            raise ValueError(
                f"No *_haplotypes.fasta found under {p / 'haplotypes'}. "
                "Pass the FASTA file directly with --haplotypes."
            )
        if len(matches) > 1:
            names = ", ".join(m.name for m in matches)
            raise ValueError(
                f"Multiple haplotype FASTAs under {p / 'haplotypes'} ({names}). "
                "Pass the specific file with --haplotypes."
            )
        analysis_dir = analysis_dir_arg if analysis_dir_arg is not None else p
        return matches[0], analysis_dir
    return p, analysis_dir_arg
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py -q`
Expected: PASS (12 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/query_assignment.py tests/test_query_assignment.py
git commit -m "Add reference-by-run-directory resolution"
```

---

### Task 3: Compute verdicts in the query flow and write `query_summary.csv`

**Files:**
- Modify: `boldgenotyper/haplotype_query.py`
- Modify: `tests/test_haplotype_query.py` (adapt the 2-tuple unpacking)

**Interfaces:**
- Consumes: `query_assignment.assign_query` (Task 1).
- Produces: `query_against_haplotypes(...)` now returns `(all_results, metadata, verdicts)` and accepts `min_identity: float = 0.97`, `tie_margin: float = 0.005`; `write_query_summary(verdicts, output_dir, min_identity, tie_margin) -> Path` writing `query_summary.csv`.

- [ ] **Step 1: Add the import at the top of `haplotype_query.py`**

Add near the other `from . import` lines:

```python
from . import query_assignment
```

- [ ] **Step 2: Extend `query_against_haplotypes` signature and compute verdicts**

2a. Change the signature (add two params, keep others):

```python
def query_against_haplotypes(
    query_fasta: Path,
    haplotype_fasta: Path,
    analysis_dir: Optional[Path] = None,
    top_n: int = 10,
    min_length: int = 100,
    max_length: int = 2000,
    min_identity: float = 0.97,
    tie_margin: float = 0.005,
) -> "Tuple[List[AlignmentResult], Optional[pd.DataFrame], List[query_assignment.QueryVerdict]]":
```

2b. Before the per-query loop (just after `metadata` is loaded and `all_results = []` is set up), build the species map and a verdicts list:

```python
    verdicts = []
    species_by_haplotype = {}
    if metadata is not None and 'haplotype_id' in metadata.columns and 'assigned_sp' in metadata.columns:
        for _, row in metadata.iterrows():
            sp = row['assigned_sp']
            species_by_haplotype[row['haplotype_id']] = sp if pd.notna(sp) and str(sp).strip() else None
```

2c. Inside the loop, replace:

```python
        query_results.sort(key=lambda x: x.identity_pct, reverse=True)
        all_results.extend(query_results[:top_n])
```

with:

```python
        query_results.sort(key=lambda x: x.identity_pct, reverse=True)
        verdicts.append(query_assignment.assign_query(
            query_id, query_results, species_by_haplotype, min_identity, tie_margin
        ))
        all_results.extend(query_results[:top_n])
```

2d. Change the `return all_results, metadata` line to:

```python
    return all_results, metadata, verdicts
```

- [ ] **Step 3: Add `write_query_summary`**

Add after `write_results` in `haplotype_query.py`:

```python
def write_query_summary(verdicts, output_dir, min_identity, tie_margin):
    """Write one-row-per-query verdict summary to query_summary.csv."""
    rows = []
    for v in verdicts:
        rows.append({
            'query_id': v.query_id,
            'haplotype_call': v.haplotype_call,
            'assigned_haplotypes': ';'.join(v.assigned_haplotypes),
            'best_identity': round(v.best_identity, 4),
            'n_tied': v.n_tied,
            'species_call': v.species_call,
            'assigned_species': v.assigned_species or '',
            'candidate_species': ';'.join(v.candidate_species),
            'reason': v.reason,
            'assign_identity': min_identity,
            'tie_margin': tie_margin,
        })
    summary_file = Path(output_dir) / "query_summary.csv"
    pd.DataFrame(rows).to_csv(summary_file, index=False)
    logger.info(f"Wrote query verdict summary: {summary_file}")
    return summary_file
```

- [ ] **Step 4: Adapt the existing test's unpacking**

In `tests/test_haplotype_query.py` line ~115, change:

```python
            results, metadata = haplotype_query.query_against_haplotypes(
```

to:

```python
            results, metadata, _verdicts = haplotype_query.query_against_haplotypes(
```

Also update the docstring example in `query_against_haplotypes` (around line 386) from
`>>> results, metadata = query_against_haplotypes(` to
`>>> results, metadata, verdicts = query_against_haplotypes(`.

- [ ] **Step 5: Write a test for `write_query_summary` and verify the flow**

```python
# append to tests/test_query_assignment.py
def test_write_query_summary_roundtrip(tmp_path):
    import pandas as pd
    from boldgenotyper import haplotype_query
    v = qa.QueryVerdict("q1", "tied", ["h1", "h2"], 99.6, 2,
                        "confident", "Sphyrna lewini", [], "reason text")
    out = haplotype_query.write_query_summary([v], tmp_path, 0.97, 0.005)
    df = pd.read_csv(out)
    assert list(df["query_id"]) == ["q1"]
    assert df.loc[0, "assigned_haplotypes"] == "h1;h2"
    assert df.loc[0, "species_call"] == "confident"
```

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py tests/test_haplotype_query.py -q -k "summary or query_against or partial or align"`
Expected: the summary test passes; the previously-passing haplotype_query tests still pass (the pre-existing `test_align_query_to_haplotype_partial_match` failure is unchanged — do not treat it as new).

- [ ] **Step 6: Commit**

```bash
git add boldgenotyper/haplotype_query.py tests/test_haplotype_query.py tests/test_query_assignment.py
git commit -m "Compute query verdicts and write query_summary.csv"
```

---

### Task 4: Enrich the text report with the verdict

**Files:**
- Modify: `boldgenotyper/haplotype_query.py`

**Interfaces:**
- Consumes: `QueryVerdict` list (Task 3).
- Produces: `write_results(...)` accepts `verdicts=None` and forwards to `write_text_report(...)`, which renders a verdict line per query in the summary section.

- [ ] **Step 1: Thread `verdicts` through `write_results`**

1a. Add `verdicts=None` to the `write_results` signature (end of its parameter list) and pass it into the `write_text_report(...)` call inside `write_results`, plus call `write_query_summary` when verdicts are present:

```python
    if verdicts is not None:
        write_query_summary(verdicts, output_dir, min_identity, tie_margin)
```

Note: `write_results` must also accept `min_identity` and `tie_margin` (add with defaults 0.97 / 0.005) so it can record them in the summary. Forward them from the caller.

1b. Add `verdicts=None` to `write_text_report`'s signature.

- [ ] **Step 2: Render the verdict in the report summary section**

In `write_text_report`, in the "Query Summary" section, when `verdicts` is provided, build a `{query_id: verdict}` map and replace each query's summary line with a verdict-aware line:

```python
    verdict_by_id = {v.query_id: v for v in (verdicts or [])}
```

For each query id, if a verdict exists, write (RawDescription-style, plain text):

```python
        v = verdict_by_id.get(query_id)
        if v is not None:
            if v.haplotype_call == "confident":
                hap = f"CONFIDENT: {v.assigned_haplotypes[0]} ({v.best_identity:.2f}%)"
            elif v.haplotype_call == "tied":
                hap = (f"TIED among {v.n_tied} haplotypes "
                       f"({', '.join(v.assigned_haplotypes)}) at {v.best_identity:.2f}%")
            else:
                hap = f"NO CONFIDENT MATCH (best {v.best_identity:.2f}%)"
            if v.species_call == "confident":
                sp = f"SPECIES: confident ({v.assigned_species})"
            elif v.species_call == "ambiguous":
                sp = f"SPECIES: ambiguous ({', '.join(v.candidate_species)})"
            elif v.species_call == "unknown":
                sp = "SPECIES: unknown (no metadata)"
            else:
                sp = "SPECIES: unassigned"
            f.write(f"{query_id} : {hap}\n")
            f.write(f"{' ' * len(query_id)}   {sp}\n")
            continue
```

(Keep the existing summary line as the fallback when `v is None`.)

- [ ] **Step 3: Verify report rendering with a temp run**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py -q`
Expected: PASS (all Task 1-3 tests still green).

- [ ] **Step 4: Commit**

```bash
git add boldgenotyper/haplotype_query.py
git commit -m "Show query verdict in the text report summary"
```

---

### Task 5: Wire flags and reference resolution into the CLI

**Files:**
- Modify: `boldgenotyper/cli.py`

**Interfaces:**
- Consumes: `query_assignment.resolve_reference` (Task 2); the extended `query_against_haplotypes` / `write_results` (Tasks 3-4).
- Produces: `main_query` gains `--assign-identity` (default 0.97) and `--tie-margin` (default 0.005); resolves `--haplotypes` as file-or-directory; passes thresholds and verdicts through.

- [ ] **Step 1: Add the two threshold flags to the `main_query` parser**

After the `--top-n` argument in `main_query`, add:

```python
    parser.add_argument(
        '--assign-identity',
        type=float,
        default=0.97,
        help='Minimum identity (fraction, 0-1) to call a match assignable '
             '(default: 0.97).'
    )
    parser.add_argument(
        '--tie-margin',
        type=float,
        default=0.005,
        help='Haplotypes within this fraction of the best match are tied with it '
             '(default: 0.005).'
    )
```

- [ ] **Step 2: Resolve the reference and thread everything through**

2a. Add the import inside the `try` block near `from . import haplotype_query`:

```python
        from . import haplotype_query, query_assignment
```

2b. Before the `query_against_haplotypes` call, resolve the reference:

```python
        try:
            reference_fasta, analysis_dir = query_assignment.resolve_reference(
                args.haplotypes, args.analysis_dir)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
```

2c. Change the call to use the resolved values, thresholds, and 3-tuple return:

```python
        results, metadata, verdicts = haplotype_query.query_against_haplotypes(
            query_fasta=args.query,
            haplotype_fasta=reference_fasta,
            analysis_dir=analysis_dir,
            top_n=args.top_n,
            min_length=args.min_length,
            max_length=args.max_length,
            min_identity=args.assign_identity,
            tie_margin=args.tie_margin,
        )
```

2d. Change the `write_results` call to pass verdicts and thresholds:

```python
        haplotype_query.write_results(
            results=results,
            output_dir=output_dir,
            metadata=metadata,
            haplotype_file=reference_fasta,
            analysis_dir=analysis_dir,
            verdicts=verdicts,
            min_identity=args.assign_identity,
            tie_margin=args.tie_margin,
        )
```

2e. Add the new output file to the logged output list:

```python
        logger.info(f"  - {output_dir / 'query_summary.csv'}")
```

- [ ] **Step 3: Verify end to end**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m py_compile boldgenotyper/cli.py`
Expected: no output.

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -c "import sys; sys.argv=['boldgenotyper-query','--help']; from boldgenotyper import cli; cli.main_query()" 2>&1 | grep -E "assign-identity|tie-margin"`
Expected: both flags appear in help.

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_query_assignment.py -q`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add boldgenotyper/cli.py
git commit -m "Wire query verdict flags and run-directory reference into CLI"
```

---

### Task 6: Documentation

**Files:**
- Modify: `README.md`

**Interfaces:** none (docs only).

- [ ] **Step 1: Update the `boldgenotyper-query` README section**

In the `boldgenotyper-query` documentation, add:

- A "Reading the results" subsection describing `query_summary.csv` and the two calls: `haplotype_call` (confident / tied / no_confident_match) and `species_call` (confident / ambiguous / unassigned / unknown), including that a haplotype tie can still yield a confident species when the tied haplotypes share one species.
- The two new flags `--assign-identity` (default 0.97) and `--tie-margin` (default 0.005) with one-line descriptions.
- The run-directory form of `--haplotypes`: pass a previous run directory and BOLDGenotyper finds `<dir>/haplotypes/*_haplotypes.fasta` and uses that directory for taxonomy metadata automatically. Show:
  ```bash
  boldgenotyper-query --query new_samples.fasta --haplotypes boldgenotyper_Penaeidae
  ```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "Document query verdict, thresholds, and run-directory reference"
```

---

## Self-Review

**Spec coverage:**
- Verdict computation (haplotype + species levels, thresholds) -> Task 1. Full-match-list requirement -> Task 3 (verdict computed before top-N truncation). `query_summary.csv` output -> Task 3. Enhanced text report -> Task 4. Reference-by-run-directory + auto analysis-dir -> Tasks 2, 5. CLI flags -> Task 5. Docs -> Task 6. Tests -> Tasks 1-3. All spec sections covered.

**Placeholder scan:** No TBD/TODO; code steps show complete code; integration steps quote exact old/new text.

**Type consistency:** `QueryVerdict` field names are consistent between Task 1 (definition), Task 3 (`write_query_summary`), and Task 4 (report rendering). `assign_query` and `resolve_reference` signatures match between definition (Tasks 1-2) and callers (Tasks 3, 5). `query_against_haplotypes` return becomes a 3-tuple in Task 3 and every caller (cli.py, test_haplotype_query.py, docstring) is updated in Tasks 3 and 5.

**Note for implementer:** the repo's wider `pytest tests/` does not collect cleanly (pre-existing, backlogged). Run the specific test files named in each task. The pre-existing `tests/test_haplotype_query.py::test_align_query_to_haplotype_partial_match` failure is not caused by this work.
