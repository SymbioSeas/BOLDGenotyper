"""Query assignment verdict: turn a query's matches into an interpretable call.

Given the full list of haplotype matches for one query and a species map, decide
whether the query is confidently one haplotype, tied among several, or not a
confident match to anything, and whether it resolves to a single species.
"""

from dataclasses import dataclass
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
