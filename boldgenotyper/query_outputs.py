"""Files written by boldgenotyper-query.

query_calls.csv       one row per query: the result to report
query_candidates.csv  the top --top-n haplotypes per query, with the numbers
query_alignments.txt  human-readable alignments
query_provenance.json reference run, parameters, curation file, outputs
query_results.json    everything above as one JSON document (--json only)

The CSVs have no comment lines, so spreadsheets open them directly; what used
to be in their '#' header lines is in query_provenance.json.
"""

import hashlib
import json
import logging
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from . import __version__
from .haplotype_query import ALIGNMENT_SCORING, format_alignment_display
from .query_assignment import is_eligible

logger = logging.getLogger(__name__)

CALL_COLUMNS = [
    "query_id", "query_length", "call_level", "equivalence_set", "n_equivalent",
    "species_call", "call_basis", "species_support", "pooled_member_labels", "differences",
    "n_mismatches", "n_gaps", "identity_pct", "query_coverage", "flags",
    "n_exact_record_matches", "reference_run",
]
CANDIDATE_COLUMNS = [
    "query_id", "rank", "haplotype_id", "haplotype_label", "assigned_sp", "assignment_level",
    "assignment_source", "n_members", "in_equivalence_set", "eligible", "differences",
    "n_mismatches", "n_gaps", "indel_in_query", "identity_pct", "query_coverage",
    "overlap_coverage", "aligned_length", "query_length", "haplotype_length",
]
POOLED_LABELS_NOTE = ("pooled_member_labels and species_support use the labels of every record "
                      "assigned to the equivalence set's haplotypes, including assigned_corrected "
                      "records; genus-only labels of the leading genus are compatible. "
                      "call_basis 'haplotype taxa': a haplotype in the set has a curated or "
                      "evidence taxon, so species_call is the set's shared haplotype taxon")


def _md5(path) -> Optional[str]:
    try:
        return hashlib.md5(Path(path).read_bytes()).hexdigest()
    except OSError:
        return None


def calls_table(verdicts, reference_run: Optional[dict] = None,
                record_matches: Optional[Dict[str, int]] = None) -> pd.DataFrame:
    """One row per query (CALL_COLUMNS)."""
    run_id = (reference_run or {}).get("provenance_md5", "")
    rows = []
    for v in verdicts:
        rows.append({
            "query_id": v.query_id,
            "query_length": v.query_length,
            "call_level": v.call_level,
            "equivalence_set": ";".join(v.equivalence_set),
            "n_equivalent": v.n_equivalent,
            "species_call": v.species_call,
            "call_basis": v.call_basis,
            "species_support": v.species_support,
            "pooled_member_labels": v.pooled_member_labels,
            "differences": v.differences,
            "n_mismatches": v.n_mismatches,
            "n_gaps": v.n_gaps,
            "identity_pct": round(v.best_identity, 4),
            "query_coverage": round(v.query_coverage, 4),
            "flags": ";".join(v.flags),
            "n_exact_record_matches": (record_matches.get(v.query_id)
                                       if record_matches is not None else None),
            "reference_run": run_id,
        })
    return pd.DataFrame(rows, columns=CALL_COLUMNS)


def write_calls(calls: pd.DataFrame, output_dir) -> Path:
    path = Path(output_dir) / "query_calls.csv"
    calls.to_csv(path, index=False)
    logger.info(f"Wrote query calls: {path}")
    return path


def candidates_table(results, verdicts, haplotypes: Optional[pd.DataFrame] = None,
                     haplotype_taxa: Optional[Dict[str, Any]] = None,
                     min_identity: float = 0.97, min_query_coverage: float = 0.95) -> pd.DataFrame:
    """
    The top matches per query (CANDIDATE_COLUMNS), in rank order.

    Taxon columns come from haplotype_taxa (which carries query-time curation)
    when given, else from the reference haplotypes table.
    """
    info = (haplotypes.drop_duplicates("haplotype_id").set_index("haplotype_id")
            if haplotypes is not None and "haplotype_id" in haplotypes.columns else None)
    sets = {v.query_id: set(v.equivalence_set) for v in verdicts}
    taxa = haplotype_taxa or {}
    rows, rank = [], {}
    for r in results:
        rank[r.query_id] = rank.get(r.query_id, 0) + 1
        meta = info.loc[r.haplotype_id] if info is not None and r.haplotype_id in info.index else None
        taxon = taxa.get(r.haplotype_id)

        def field(name):
            return meta.get(name) if meta is not None and name in meta.index else None

        rows.append({
            "query_id": r.query_id,
            "rank": rank[r.query_id],
            "haplotype_id": r.haplotype_id,
            "haplotype_label": field("haplotype_label"),
            "assigned_sp": taxon.taxon if taxon else field("assigned_sp"),
            "assignment_level": taxon.level if taxon else field("assignment_level"),
            "assignment_source": taxon.source if taxon else field("assignment_source"),
            "n_members": field("n_members"),
            "in_equivalence_set": r.haplotype_id in sets.get(r.query_id, set()),
            "eligible": is_eligible(r, min_identity, min_query_coverage),
            "differences": r.differences,
            "n_mismatches": r.n_mismatches,
            "n_gaps": r.n_gaps,
            "indel_in_query": r.indel_in_query,
            "identity_pct": round(r.identity_pct, 4),
            "query_coverage": round(r.query_coverage, 4),
            "overlap_coverage": round(r.overlap_coverage, 4),
            "aligned_length": r.aligned_length,
            "query_length": r.query_length,
            "haplotype_length": r.haplotype_length,
        })
    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)


def write_candidates(candidates: pd.DataFrame, output_dir) -> Path:
    path = Path(output_dir) / "query_candidates.csv"
    candidates.to_csv(path, index=False)
    logger.info(f"Wrote query candidates: {path}")
    return path


def _haplotype_row(haplotypes, haplotype_id) -> dict:
    if haplotypes is None or "haplotype_id" not in haplotypes.columns:
        return {}
    match = haplotypes[haplotypes["haplotype_id"] == haplotype_id]
    return match.iloc[0].to_dict() if len(match) else {}


def _group_by_query(results) -> Dict[str, list]:
    grouped: Dict[str, list] = {}
    for r in results:
        grouped.setdefault(r.query_id, []).append(r)
    return grouped


def write_json(results, verdicts, output_dir, parameters: Optional[dict] = None,
               reference_run: Optional[dict] = None,
               haplotypes: Optional[pd.DataFrame] = None) -> Path:
    """query_results.json (--json): reference run, parameters, and per query
    the verdict and the top matches."""
    path = Path(output_dir) / "query_results.json"
    verdict_by_id = {v.query_id: v for v in (verdicts or [])}
    queries = []
    for qid, matches in _group_by_query(results).items():
        entry = {"query_id": qid, "query_length": matches[0].query_length, "matches": []}
        for i, r in enumerate(matches, 1):
            entry["matches"].append({
                "rank": i, "haplotype_id": r.haplotype_id, "identity_pct": r.identity_pct,
                "differences": r.differences, "n_mismatches": r.n_mismatches, "n_gaps": r.n_gaps,
                "indel_in_query": r.indel_in_query, "query_coverage": r.query_coverage,
                "overlap_coverage": r.overlap_coverage,
                "alignment": {"matches": r.matches, "length": r.aligned_length, "score": r.score,
                              "query_start": r.query_start, "query_end": r.query_end,
                              "haplotype_start": r.haplotype_start, "haplotype_end": r.haplotype_end},
                "metadata": _haplotype_row(haplotypes, r.haplotype_id),
            })
        if qid in verdict_by_id:
            entry["verdict"] = asdict(verdict_by_id[qid])
        queries.append(entry)
    n_reference = (verdicts[0].n_reference_haplotypes if verdicts
                   else len({r.haplotype_id for r in results}))
    output = {
        "query_metadata": {"timestamp": datetime.now().isoformat(), "n_queries": len(queries),
                           "n_reference_haplotypes": n_reference},
        "reference_run": reference_run,
        "parameters": parameters or {},
        "alignment": ALIGNMENT_SCORING,
        "note": POOLED_LABELS_NOTE,
        "queries": queries,
    }
    with open(path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    logger.info(f"Wrote JSON results: {path}")
    return path


def _species_line(v) -> str:
    if v.call_basis == "haplotype taxa":
        support = f"; pooled support {v.species_support:.2f}" if v.species_support is not None else ""
        return f"species {v.species_call} (haplotype taxa{support})"
    if v.species_call:
        return f"species {v.species_call} (support {v.species_support:.2f})"
    if v.species_support is not None:
        return f"no species call (support {v.species_support:.2f})"
    return "no species call"


def write_alignments(results, verdicts, output_dir, reference_run: Optional[dict] = None,
                     haplotypes: Optional[pd.DataFrame] = None,
                     haplotype_file=None) -> Path:
    """query_alignments.txt: each query's call and its best alignment."""
    path = Path(output_dir) / "query_alignments.txt"
    grouped = _group_by_query(results)
    verdict_by_id = {v.query_id: v for v in (verdicts or [])}
    n_reference = (verdicts[0].n_reference_haplotypes if verdicts
                   else len({r.haplotype_id for r in results}))

    def taxon(hid):
        row = _haplotype_row(haplotypes, hid)
        return row.get("assigned_sp") or ""

    with open(path, "w") as f:
        f.write("=" * 80 + "\nBOLDGenotyper Haplotype Query Alignments\n" + "=" * 80 + "\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Queries Processed: {len(grouped)}\n")
        f.write(f"Haplotypes Searched: {n_reference}\n")
        if reference_run:
            f.write(f"Reference Run: {reference_run.get('path')} (BOLDGenotyper "
                    f"{reference_run.get('boldgenotyper_version')}, run {reference_run.get('run_date')}, "
                    f"run_provenance.json md5 {reference_run.get('provenance_md5')})\n")
        if haplotype_file:
            f.write(f"Haplotype File: {haplotype_file}\n")
        f.write("\nQuery Summary:\n" + "-" * 80 + "\n")
        for qid, matches in grouped.items():
            v = verdict_by_id.get(qid)
            if v is None:
                best = matches[0]
                f.write(f"{qid:30s}: best {best.haplotype_id} "
                        f"({best.identity_pct:.2f}%, {best.differences} differences)\n")
                continue
            if v.call_level == "none":
                hap = f"no eligible match (best {v.best_identity:.2f}%, coverage {v.query_coverage:.2f})"
            else:
                hap = (f"{v.call_level} {';'.join(v.equivalence_set)} "
                       f"({v.differences} difference{'' if v.differences == 1 else 's'}, "
                       f"{v.best_identity:.2f}% identity, coverage {v.query_coverage:.2f})")
            pad = " " * len(qid)
            f.write(f"{qid} : {hap}\n")
            f.write(f"{pad}   {_species_line(v)}; pooled labels: {v.pooled_member_labels or '-'}\n")
            if v.flags:
                f.write(f"{pad}   flags: {', '.join(v.flags)}\n")
        f.write("\n" + "=" * 80 + "\nDetailed Results\n" + "=" * 80 + "\n\n")
        for qid, matches in grouped.items():
            best = matches[0]
            f.write(f"Query: {qid}\n")
            f.write(f"Best Match: {best.haplotype_id} (Identity: {best.identity_pct:.2f}%)\n")
            row = _haplotype_row(haplotypes, best.haplotype_id)
            if row.get("assigned_sp"):
                level = row.get("assignment_level")
                level_str = f", {level}-level assignment" if level and level != "species" else ""
                f.write(f"Species: {row['assigned_sp']} ({row.get('n_members', 'unknown')} samples{level_str})\n")
            f.write(f"\nAlignment ({best.aligned_length} bp):\n")
            f.write(format_alignment_display(best))
            f.write(f"Matches: {best.matches}/{best.aligned_length} bp ({best.identity_pct:.2f}%)\n")
            f.write(f"Mismatches: {best.n_mismatches}; gap columns: {best.n_gaps}"
                    f"{' (indel not a multiple of 3)' if best.indel_in_query else ''}; "
                    f"query coverage {best.query_coverage:.3f}\n")
            if len(matches) > 1:
                f.write("\nOther Top Matches:\n")
                for i, r in enumerate(matches[1:], 2):
                    sp = taxon(r.haplotype_id)
                    f.write(f"  {i}. {r.haplotype_id}: {r.identity_pct:.2f}%, "
                            f"{r.differences} differences{f' ({sp})' if sp else ''}\n")
            f.write("\n" + "-" * 80 + "\n\n")
    logger.info(f"Wrote alignments: {path}")
    return path


def build_query_provenance(argv, query_path, n_queries: int, reference_fasta,
                           n_reference_haplotypes: int, reference_run: Optional[dict],
                           parameters: dict, taxon_curation: Optional[dict] = None) -> dict:
    """The contents of query_provenance.json; later steps add to query_metadata,
    figures and warnings."""
    return {
        "boldgenotyper_version": __version__,
        "run_date": datetime.now().isoformat(timespec="seconds"),
        "command": [str(a) for a in argv],
        "query": {"path": str(query_path), "md5": _md5(query_path), "n_sequences": n_queries},
        "reference_fasta": {"path": str(reference_fasta), "md5": _md5(reference_fasta),
                            "n_haplotypes": n_reference_haplotypes},
        "reference_run": reference_run,
        "parameters": parameters,
        "alignment": ALIGNMENT_SCORING,
        "note": POOLED_LABELS_NOTE,
        "taxon_curation": taxon_curation,
        "reference_taxon_curation": None,
        "reference_taxon_evidence": None,
        "query_metadata": None,
        "figures": {},
        "warnings": [],
        "outputs": [],
    }


def write_query_provenance(provenance: dict, output_dir) -> Path:
    """Write query_provenance.json, listing the files present in output_dir."""
    output_dir = Path(output_dir)
    path = output_dir / "query_provenance.json"
    names = {str(p.relative_to(output_dir)) for p in output_dir.rglob("*") if p.is_file()}
    provenance["outputs"] = sorted(names | {path.name})
    with open(path, "w") as f:
        json.dump(provenance, f, indent=2, default=str)
    logger.info(f"Wrote query provenance: {path}")
    return path


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------

_CALL_DOCS = {
    "query_id": "query sequence ID (first word of the FASTA header)",
    "query_length": "query length in bp, gaps removed",
    "call_level": "haplotype (one haplotype in the equivalence set), haplotype_set (several "
                  "haplotypes the query cannot tell apart) or none (no eligible match)",
    "equivalence_set": "the eligible haplotypes within --max-mismatch-delta differences of the "
                       "best match, ';'-joined",
    "n_equivalent": "haplotypes in equivalence_set",
    "species_call": "the species called for the query; empty when none (see call_basis)",
    "call_basis": "pooled labels (the leading species of the set's pooled member labels, at "
                  "--min-species-support or more), haplotype taxa (a haplotype in the set has a "
                  "curated or evidence taxon: the taxon the set's haplotypes share), or empty "
                  "(no call)",
    "species_support": "the leading species' share of the counted pooled member labels",
    "pooled_member_labels": "labels of every member record of the set's haplotypes, with counts",
    "differences": "mismatches plus gap columns between the query and its best match, over "
                   "their overlap",
    "n_mismatches": "mismatches with the best match",
    "n_gaps": "gap columns inside the overlap with the best match",
    "identity_pct": "identity with the best match over the columns with a base in both (%)",
    "query_coverage": "share of the query's bases in the overlap with the best match",
    "flags": "';'-joined: indel_in_query (a gap run not a multiple of 3), low_coverage, "
             "below_identity, set_not_identical (the set's haplotypes differ over the query's "
             "region), label_conflict_in_set, species_support_below_threshold, "
             "no_species_labels, no_reference_labels (a FASTA reference or v1 run folder), "
             "haplotype_taxa_disagree (a haplotype with a curated or evidence taxon is in the "
             "set, and the set's taxa differ or one is not at species level)",
    "n_exact_record_matches": "QC-passed reference records that contain the query unchanged or "
                              "are contained in it (--search-records); empty otherwise",
    "reference_run": "md5 of the reference run's run_provenance.json; empty for a FASTA reference",
}
_CANDIDATE_DOCS = {
    "query_id": "query sequence ID",
    "rank": "order of the match for this query: eligible first, then fewer differences, higher "
            "coverage, higher identity, longer alignment, haplotype ID",
    "haplotype_id": "reference haplotype",
    "haplotype_label": "its label in the reference run",
    "assigned_sp": "the haplotype's taxon (with --taxon-curation, as curated for this query)",
    "assignment_level": "species or genus",
    "assignment_source": "labels, evidence or curated",
    "n_members": "records assigned to the haplotype in the reference run",
    "in_equivalence_set": "the haplotype is in the query's equivalence set",
    "eligible": "the alignment covers --min-query-coverage of the shorter sequence at "
                "--assign-identity or more",
    "differences": "mismatches plus gap columns over the overlap",
    "n_mismatches": "mismatches",
    "n_gaps": "gap columns inside the overlap",
    "indel_in_query": "a gap run whose length is not a multiple of 3",
    "identity_pct": "identity over the columns with a base in both (%)",
    "query_coverage": "share of the query's bases in the overlap",
    "overlap_coverage": "query bases in the overlap / the shorter of query and haplotype; "
                        "--min-query-coverage applies to this",
    "aligned_length": "columns with a base in both sequences",
    "query_length": "query length in bp",
    "haplotype_length": "haplotype length in bp",
}
_RECORD_MATCH_DOCS = {
    "query_id": "query sequence ID",
    "processid": "reference record",
    "match_type": "query_within_record or record_within_query",
    "label": "the record's species, or genus when the species is missing",
    "haplotype_id": "the record's haplotype; empty when it is not assigned to one",
    "status": "the record's status in the reference run",
    "record_length": "record sequence length in bp",
}

QUERY_FILE_DOCS = [
    {"pattern": "README.md", "text": "This file."},
    {"pattern": "query_calls.csv", "text": "One row per query: the result to report.",
     "columns": _CALL_DOCS},
    {"pattern": "query_candidates.csv",
     "text": "The top --top-n haplotypes per query with the numbers behind the call. Other "
             "haplotype columns: join the reference run's haplotypes.csv on haplotype_id.",
     "columns": _CANDIDATE_DOCS},
    {"pattern": "query_record_matches.csv",
     "text": "--search-records: QC-passed reference records that contain the query unchanged or "
             "are contained in it, including records not assigned to a haplotype.",
     "columns": _RECORD_MATCH_DOCS},
    {"pattern": "query_alignments.txt", "text": "Each query's call and its best alignment."},
    {"pattern": "query_provenance.json",
     "text": "The reference run (path, version, run_provenance.json md5), the query and reference "
             "files with their md5, parameters, the curation and query metadata files, figures, "
             "warnings and the files written."},
    {"pattern": "query_report.html",
     "text": "The calls, one evidence card per query, the figures and the methods, in a browser."},
    {"pattern": "query.log", "text": "The log of this query run."},
    {"pattern": "query_results.json",
     "text": "--json: the calls, candidates, reference run and parameters as one JSON document."},
    {"pattern": "query_figures/query_regions.csv",
     "text": "Reference records of each equivalence set's haplotypes per region of the reference "
             "run, with the queries in each region (query regions from --query-metadata).",
     "columns": {
         "equivalence_set": "the queries' equivalence set (query_calls.csv)",
         "species_call": "the set's species call",
         "region": "a region of the reference run's region column, or a query's own region value",
         "n_reference_records": "reference records of the set's haplotypes in the region",
         "n_queries": "queries with this set placed in the region",
         "query_ids": "those queries, ';'-joined",
     }},
    {"pattern": "query_figures/*.pdf",
     "text": "Figures placing the queries in the reference run, one page per taxon, genus or "
             "equivalence set (see Figures)."},
    {"pattern": "query_figures/*/*.svg", "text": "The same figure pages, one SVG per page."},
]


def methods_sentence(provenance: dict) -> str:
    """A methods paragraph with this query run's values."""
    p = provenance.get("parameters") or {}
    ref = provenance.get("reference_run")
    a = provenance.get("alignment") or ALIGNMENT_SCORING
    n_q = (provenance.get("query") or {}).get("n_sequences")
    fasta = provenance.get("reference_fasta") or {}
    n_h = fasta.get("n_haplotypes")
    if ref:
        target = (f"the {n_h} haplotypes of a BOLDGenotyper {ref.get('boldgenotyper_version')} "
                  f"reference run (run_provenance.json md5 {ref.get('provenance_md5')})")
    else:
        target = f"the {n_h} haplotypes in {Path(str(fasta.get('path', ''))).name}"
    subject = (f"The {n_q} query sequences were each" if n_q != 1 else "The query sequence was")
    text = (f"{subject} aligned to {target} over their whole overlap (global alignment with free "
            f"end gaps; match {a.get('match'):+d}, mismatch {a.get('mismatch')}, gap opening "
            f"{a.get('gap_open')}, gap extension {a.get('gap_extend')}). A haplotype was an eligible "
            f"match when the alignment covered at least {p.get('min_query_coverage', 0):.0%} of the "
            f"shorter sequence at {p.get('assign_identity', 0) * 100:.1f}% identity or more; the "
            f"equivalence set was every eligible haplotype within {p.get('max_mismatch_delta')} "
            f"difference(s) (mismatches plus gap columns) of the best match.")
    if ref:
        text += (f" A species was called when one species held at least "
                 f"{p.get('min_species_support', 0):.2f} of the pooled labels of the member "
                 f"records of the set's haplotypes")
        sources = []
        ref_ev = provenance.get("reference_taxon_evidence")
        if ref_ev and ref_ev.get("n_applied"):
            sources.append(f"{ref_ev.get('n_applied')} haplotype taxa set by evidence in the reference run")
        ref_cur = provenance.get("reference_taxon_curation")
        if ref_cur:
            sources.append(f"curation file md5 {ref_cur.get('md5')}, applied in the reference run")
        cur = provenance.get("taxon_curation")
        if cur:
            sources.append(f"curation file md5 {cur.get('md5')}, applied to this query run")
        if sources:
            text += (f", or, when a haplotype in the set had a curated or evidence taxon "
                     f"({'; '.join(sources)}), when the set's haplotypes shared one species")
        text += "."
    return text


def _doc_for(rel: str) -> Optional[dict]:
    from fnmatch import fnmatch
    return next((d for d in QUERY_FILE_DOCS if fnmatch(rel, d["pattern"])), None)


def write_query_readme(output_dir, provenance: dict):
    """
    Write README.md for a query run folder: which file to report, the methods
    sentence, the reference run, and every file present with its columns.
    Returns (path, undocumented files or columns).
    """
    output_dir = Path(output_dir)
    path = output_dir / "README.md"
    files = sorted({str(p.relative_to(output_dir).as_posix())
                    for p in output_dir.rglob("*") if p.is_file()} | {"README.md"})
    undocumented, sections, seen = [], [], set()
    for rel in files:
        doc = _doc_for(rel)
        if doc is None:
            undocumented.append(rel)
            continue
        columns = doc.get("columns")
        if columns and rel.endswith(".csv") and (output_dir / rel).exists():
            header = list(pd.read_csv(output_dir / rel, nrows=0).columns)
            undocumented += [f"{rel}: {c}" for c in header if c not in columns]
        if doc["pattern"] in seen:
            continue
        seen.add(doc["pattern"])
        lines = [f"### `{doc['pattern']}`", "", doc["text"], ""]
        if columns:
            lines += ["| Column | Meaning |", "|---|---|"]
            lines += [f"| `{c}` | {m} |" for c, m in columns.items()]
            lines.append("")
        sections.append("\n".join(lines))

    ref = provenance.get("reference_run")
    query = provenance.get("query") or {}
    if ref:
        ref_lines = [f"- Path: `{ref.get('path')}`",
                     f"- BOLDGenotyper {ref.get('boldgenotyper_version')}, run {ref.get('run_date')}",
                     f"- run_provenance.json md5: `{ref.get('provenance_md5')}`",
                     f"- Input md5: `{ref.get('input_md5')}`; BOLD download date: "
                     f"{ref.get('download_date') or 'not recorded'}"]
    else:
        ref_lines = [f"- A haplotype FASTA (`{(provenance.get('reference_fasta') or {}).get('path')}`), "
                     "not a v2 run folder: haplotype matching only, no species calls."]
    text = "\n".join([
        "# Query results", "",
        f"{query.get('n_sequences')} query sequence(s) from `{query.get('path')}`, "
        f"BOLDGenotyper {provenance.get('boldgenotyper_version')}, {provenance.get('run_date')}.", "",
        "## Which file to report", "",
        "`query_calls.csv`: one row per query. Report `call_level`, `equivalence_set` and "
        "`species_call` with its `call_basis`, and cite the reference run by `reference_run` "
        "(the md5 of its run_provenance.json). `query_candidates.csv` holds the evidence.", "",
        "## Methods", "", methods_sentence(provenance), "",
        "## Reference run", "", *ref_lines, "",
        "## Repeat queries", "",
        "Keep one reference run per reference dataset and query the full set of samples each "
        "time (a query run takes seconds). Every output records the reference run's md5, so "
        "results from different dates can be compared; a different md5 means a different "
        "reference.", "",
        "## Files", "", *sections,
    ])
    path.write_text(text.rstrip() + "\n")
    if undocumented:
        logger.warning(f"README does not describe: {', '.join(undocumented)}")
    return path, undocumented
