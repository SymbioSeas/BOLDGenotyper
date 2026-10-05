"""Query assignment verdict: turn a query's matches into numbers and a call.

Given every haplotype match for one query (from the overlap alignment) and the
reference run's member labels, decide which haplotypes the query cannot tell
apart (its equivalence set) and whether their pooled member labels support a
species call. No interpretive labels: the verdict reports numbers and flags.

Eligible match: overlap_coverage >= min_query_coverage and identity_pct >=
100 * min_identity. Ranking: eligible first, then fewer differences
(mismatches + gap columns), higher coverage, higher identity, longer
alignment, haplotype ID. Equivalence set: eligible matches within
max_mismatch_delta differences of the best.

Species call: when a haplotype in the equivalence set has a curated or
evidence taxon, the call uses the haplotypes' taxa (curated, evidence, else the
label rule's): one shared species is called (call_basis 'haplotype taxa');
otherwise no call (flag haplotype_taxa_disagree). When none has, the pooled
member labels are called at min_species_support or more (call_basis 'pooled
labels').

Flags: indel_in_query, low_coverage, below_identity, set_not_identical,
label_conflict_in_set, species_support_below_threshold, no_species_labels,
no_reference_labels, haplotype_taxa_disagree.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence

from . import taxon_labels

REVIEWED_SOURCES = ("curated", "evidence")    # haplotype taxa not from the pooled labels

CALL_LEVELS = ("haplotype", "haplotype_set", "none")


class HaplotypeTaxon(NamedTuple):
    taxon: str                            # assigned_sp
    level: str                            # species | genus | unassigned
    source: str                           # labels | curated


@dataclass(frozen=True)
class QueryVerdict:
    query_id: str
    query_length: int
    best_haplotype: Optional[str]
    best_identity: float
    aligned_length: int
    query_coverage: float
    overlap_coverage: float
    n_mismatches: int
    n_gaps: int
    differences: Optional[int]
    call_level: str                       # haplotype | haplotype_set | none
    equivalence_set: List[str]
    n_equivalent: int
    set_identical_over_region: Optional[bool]
    pooled_member_labels: str
    species_call: str                     # "" when no species is called
    species_support: Optional[float]
    flags: List[str]
    notes: str
    max_mismatch_delta: int               # delta applied to this query
    n_reference_haplotypes: int = 0
    call_basis: str = ""                  # pooled labels | haplotype taxa | "" (no call)


def is_eligible(match, min_identity: float, min_query_coverage: float) -> bool:
    return (match.overlap_coverage >= min_query_coverage
            and match.identity_pct >= min_identity * 100.0)


def rank_matches(matches: Sequence, min_identity: float, min_query_coverage: float) -> list:
    """Order matches: eligible first, then fewer differences, higher coverage,
    higher identity, longer alignment, haplotype ID."""
    return sorted(matches, key=lambda m: (
        not is_eligible(m, min_identity, min_query_coverage),
        m.differences, -m.overlap_coverage, -m.identity_pct, -m.aligned_length,
        m.haplotype_id,
    ))


def assign_query(query_id, matches, member_labels: Optional[Dict[str, List[str]]], *,
                 min_identity: float = 0.97, min_query_coverage: float = 0.95,
                 max_mismatch_delta: int = 0, min_species_support: float = 0.70,
                 haplotype_taxa: Optional[Dict[str, HaplotypeTaxon]] = None) -> QueryVerdict:
    """Compute a QueryVerdict from a query's full list of haplotype matches.

    member_labels maps haplotype_id -> labels of its member records (exact
    and corrected); None when the reference has no record labels (bare FASTA
    or a v1 run folder). haplotype_taxa maps haplotype_id -> HaplotypeTaxon;
    a curated taxon in the equivalence set makes the call from the
    haplotypes' taxa (module docstring).
    """
    ranked = rank_matches(matches, min_identity, min_query_coverage)
    eligible = [m for m in ranked if is_eligible(m, min_identity, min_query_coverage)]
    top = ranked[0] if ranked else None
    flags: List[str] = []

    def _verdict(best, call_level, eq, identical, pooled_text, species_call, support, delta, notes,
                 call_basis=""):
        return QueryVerdict(
            query_id=query_id,
            query_length=(best or top).query_length if (best or top) else 0,
            best_haplotype=best.haplotype_id if best else None,
            best_identity=best.identity_pct if best else 0.0,
            aligned_length=best.aligned_length if best else 0,
            query_coverage=best.query_coverage if best else 0.0,
            overlap_coverage=best.overlap_coverage if best else 0.0,
            n_mismatches=best.n_mismatches if best else 0,
            n_gaps=best.n_gaps if best else 0,
            differences=best.differences if best else None,
            call_level=call_level,
            equivalence_set=[m.haplotype_id for m in eq],
            n_equivalent=len(eq),
            set_identical_over_region=identical,
            pooled_member_labels=pooled_text,
            species_call=species_call,
            species_support=support,
            flags=flags,
            notes=notes,
            max_mismatch_delta=delta,
            n_reference_haplotypes=len(matches),
            call_basis=call_basis,
        )

    if not eligible:
        if top is None:
            return _verdict(None, "none", [], None, "", "", None, max_mismatch_delta,
                            "no haplotypes to match")
        if not any(m.overlap_coverage >= min_query_coverage for m in ranked):
            flags.append("low_coverage")
            why = f"no match covers {min_query_coverage:.0%} of the overlap"
        else:
            flags.append("below_identity")
            why = f"no match with enough coverage reaches {min_identity * 100:.1f}% identity"
        return _verdict(top, "none", [], None, "", "", None, max_mismatch_delta, why)

    best = eligible[0]
    if best.indel_in_query:
        flags.append("indel_in_query")
    delta = int(max_mismatch_delta)
    eq = [m for m in eligible if m.differences <= best.differences + delta]
    identical = len({m.haplotype_region for m in eq}) == 1
    if not identical:
        flags.append("set_not_identical")
    call_level = "haplotype" if len(eq) == 1 else "haplotype_set"
    notes = (f"{len(eq)} haplotype(s) within {delta} difference(s) of the best "
             f"({best.differences} differences, {best.identity_pct:.2f}% identity)")

    if member_labels is None:
        flags.append("no_reference_labels")
        return _verdict(best, call_level, eq, identical, "", "", None, delta,
                        notes + "; species support needs a v2 run folder")

    pooled = taxon_labels.pool_labels(
        label for m in eq for label in member_labels.get(m.haplotype_id, []))
    if pooled.n_distinct_species > 1:
        flags.append("label_conflict_in_set")
    support = round(pooled.support, 4) if pooled.support is not None else None

    taxa = haplotype_taxa or {}
    if any(taxa.get(m.haplotype_id) and taxa[m.haplotype_id].source in REVIEWED_SOURCES for m in eq):
        set_taxa = [taxa.get(m.haplotype_id) for m in eq]
        species = {t.taxon for t in set_taxa if t and t.level == "species"}
        if len(species) == 1 and all(t and t.level == "species" for t in set_taxa):
            return _verdict(best, call_level, eq, identical, pooled.text, species.pop(), support,
                            delta, notes + "; called from haplotype taxa (one is curated or set by evidence)",
                            call_basis="haplotype taxa")
        flags.append("haplotype_taxa_disagree")
        listed = "; ".join(f"{m.haplotype_id} {t.taxon} ({t.source})" if t else
                           f"{m.haplotype_id} no taxon" for m, t in zip(eq, set_taxa))
        return _verdict(best, call_level, eq, identical, pooled.text, "", support, delta,
                        notes + f"; haplotype taxa differ: {listed}")

    species_call = ""
    if pooled.support is None:
        flags.append("no_species_labels")
    elif pooled.leading_species and pooled.support >= min_species_support:
        species_call = pooled.leading_species
    else:
        flags.append("species_support_below_threshold")
    return _verdict(best, call_level, eq, identical, pooled.text, species_call, support,
                    delta, notes, call_basis="pooled labels" if species_call else "")


# Record tables that carry per-record statuses (v2): records.csv (v2 layout) or
# <organism>_annotated.csv (v2.0.0 development layout).
RECORD_TABLE_GLOBS = ("records.csv", "*_annotated.csv")
MEMBER_STATUSES = ("assigned", "assigned_corrected")


def load_member_labels(run_dir) -> tuple:
    """Labels of each haplotype's member records in a v2 run folder.

    Returns (labels_by_haplotype or None, reason). A record's label is its
    species label, or its genus when the species is missing. Members are the
    records with status assigned or assigned_corrected. None (with the reason)
    for a bare FASTA reference, a v1 run folder (no run_provenance.json) or a
    folder without a record table.
    """
    import pandas as pd

    if run_dir is None:
        return None, "the reference is a FASTA file, not a run folder"
    run_dir = Path(run_dir)
    if not (run_dir / "run_provenance.json").exists():
        return None, (f"{run_dir} has no run_provenance.json (a v1 run folder?); "
                      "species support needs a v2 run folder")
    df, table = load_record_table(run_dir)
    if df is None:
        return None, f"no record table with statuses ({' or '.join(RECORD_TABLE_GLOBS)}) in {run_dir}"
    df = df[df["status"].isin(MEMBER_STATUSES) & df["haplotype_id"].notna()]
    labels = {}
    for hid, lab in zip(df["haplotype_id"], df["label"]):
        labels.setdefault(hid, []).append(lab)
    return labels, f"member labels from {table.name}"


def load_haplotype_taxa(run_dir) -> Optional[Dict[str, HaplotypeTaxon]]:
    """Each haplotype's taxon from a v2 run folder's haplotypes.csv.

    source is assignment_source (labels when the column is absent: runs made
    before curation existed). None for a FASTA reference or a folder without
    run_provenance.json and haplotypes.csv.
    """
    import pandas as pd

    if run_dir is None:
        return None
    run_dir = Path(run_dir)
    table = run_dir / "haplotypes.csv"
    if not (run_dir / "run_provenance.json").exists() or not table.exists():
        return None
    df = pd.read_csv(table, dtype=str, keep_default_na=False)
    if not {"haplotype_id", "assigned_sp", "assignment_level"} <= set(df.columns):
        return None
    source = df["assignment_source"] if "assignment_source" in df.columns else ["labels"] * len(df)
    return {h: HaplotypeTaxon(t, lvl, src or "labels")
            for h, t, lvl, src in zip(df["haplotype_id"], df["assigned_sp"],
                                      df["assignment_level"], source)}


def load_record_table(run_dir) -> tuple:
    """Per-record processid, species, genus, label, haplotype_id and status of a v2 run.

    label = species, or genus when the species is missing. Returns
    (DataFrame or None, path of the table or None).
    """
    import pandas as pd

    run_dir = Path(run_dir)
    table = next((p for pattern in RECORD_TABLE_GLOBS for p in sorted(run_dir.glob(pattern))), None)
    if table is None:
        return None, None
    wanted = {"processid", "haplotype_id", "status", "species", "genus"}
    df = pd.read_csv(table, usecols=lambda c: c in wanted, dtype=str, low_memory=False)
    if not {"processid", "haplotype_id", "status"} <= set(df.columns):
        return None, table
    for col in ("species", "genus"):
        if col not in df.columns:
            df[col] = ""
    species = df["species"].fillna("").str.strip()
    df["label"] = species.where(species != "", df["genus"].fillna("").str.strip())
    return df, table


def load_reference_run(path) -> Optional[dict]:
    """The reference_run block (Q9) for a v2 run folder; None otherwise.

    provenance_md5 (md5 of run_provenance.json; the reference_run column of query_calls.csv),
    boldgenotyper_version, n_haplotypes, run_date, input_md5, download_date,
    path.
    """
    import hashlib
    import json

    if path is None:
        return None
    path = Path(path)
    prov_file = path / "run_provenance.json"
    if not path.is_dir() or not prov_file.exists():
        return None
    raw = prov_file.read_bytes()
    prov = json.loads(raw)
    applied = prov.get("applied", {})
    return {
        "provenance_md5": hashlib.md5(raw).hexdigest(),
        "boldgenotyper_version": prov.get("boldgenotyper_version"),
        "n_haplotypes": applied.get("n_haplotypes"),
        "run_date": prov.get("started_at"),
        "input_md5": prov.get("input", {}).get("md5"),
        "download_date": prov.get("input", {}).get("download_date"),
        "path": str(path),
    }


def resolve_reference(haplotypes_arg, analysis_dir_arg):
    """Resolve --haplotypes (a FASTA file or a run directory) to
    (reference_fasta_path, analysis_dir)."""
    p = Path(haplotypes_arg)
    if p.is_dir():
        if (p / "haplotypes.fasta").exists():          # v2 layout
            analysis_dir = analysis_dir_arg if analysis_dir_arg is not None else p
            return p / "haplotypes.fasta", analysis_dir
        matches = sorted(p.glob("haplotypes/*_haplotypes.fasta"))   # v1 layout
        if not matches:
            raise ValueError(
                f"No haplotypes.fasta (or haplotypes/*_haplotypes.fasta) in {p}. "
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
