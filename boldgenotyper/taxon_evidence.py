"""Evidence for a haplotype's taxon, reported beside the label rule's call.

Evidence is never applied. A suggested taxon prompts a check of the haplotype;
a call the user accepts goes in a curation file (taxon_curation) with its reason.

Columns (EVIDENCE_COLUMNS), one row per haplotype:

    nearest_species_haplotype   nearest other haplotype assigned at species level
                                (';'-joined when several are equally near)
    nearest_species_taxon       its taxon (';'-joined when the nearest differ)
    nearest_species_distance    the p-distance to it
    taxon_max_within_distance   largest p-distance between the other species-level
                                haplotypes of nearest_species_taxon (empty when it
                                has fewer than two, or the nearest differ in taxon)
    within_taxon_spread         nearest_species_distance <= taxon_max_within_distance
    member_bins                 members' BOLD BINs with counts ('none' = no BIN)
    bin_taxon                   for the members' most common BIN: the majority
                                species-level label among the run's QC-passed records
                                in that BIN, its fraction and counts
    label_sources               per member label: records, BOLD projects (processid
                                prefix; GenBank accession prefix for mined records)
                                and institutions behind it
    tree_sister_taxon           walking from the haplotype's tip toward the root, the
                                first clade holding another species-level haplotype:
                                their taxon when they share one, else empty
    tree_sister_support         that clade's support value
    suggested_sp                for haplotypes neither at species level nor curated:
                                a taxon named by at least two of the three lines
                                (distance within the taxon's spread; BIN taxon at
                                BIN_MAJORITY or more; tree sister) and contradicted
                                by none
    suggestion_basis            which lines agree, or why there is no suggestion

Distances are dereplication's p-distance over the aligned core sequences, the
definition behind nearest_distance and the label-conflict evidence.
"""

import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .dereplication import distances_to_references
from .taxon_labels import label_level

EVIDENCE_COLUMNS = [
    "nearest_species_haplotype", "nearest_species_taxon", "nearest_species_distance",
    "taxon_max_within_distance", "within_taxon_spread", "member_bins", "bin_taxon",
    "label_sources", "tree_sister_taxon", "tree_sister_support", "suggested_sp",
    "suggestion_basis",
]
BIN_MAJORITY = 0.70
_TOL = 1e-12


def _text(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series("", index=df.index, dtype=object)
    return df[col].astype(object).where(df[col].notna(), "").astype(str).str.strip()


def _plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


def label_source(processid, inst="", insdc_acs="") -> str:
    """
    Where a record's label comes from, for counting independent sources: the
    BOLD project code (processid's leading letters), or for records mined from
    GenBank, whose project codes are shared across studies, 'GenBank:' and the
    accession's letters and first three digits.
    """
    if "genbank" in str(inst or "").lower():
        m = re.match(r"([A-Za-z_]+\d{0,3})", str(insdc_acs or "").strip())
        return f"GenBank:{m.group(1)}" if m else "GenBank"
    m = re.match(r"[A-Za-z]+", str(processid or "").strip())
    return m.group(0) if m else str(processid or "")


def _label_sources(members: pd.DataFrame) -> str:
    parts = []
    labeled = members[members["label"] != ""]
    for lab, g in labeled.groupby("label"):
        projects = {label_source(p, i, a) for p, i, a in zip(g["processid"], g["inst"], g["insdc_acs"])}
        insts = {i for i in g["inst"] if i}
        parts.append((-len(g), lab, f"{lab}: {_plural(len(g), 'record')}, "
                                    f"{_plural(len(projects), 'project')}, "
                                    f"{_plural(len(insts), 'institution')}"))
    return "; ".join(text for _, _, text in sorted(parts))


def _member_bins(members: pd.DataFrame) -> str:
    counts = Counter(b if b else "none" for b in members["bin_uri"])
    return "; ".join(f"{b} {n}" for b, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def _bin_taxon(members: pd.DataFrame, qc_passed: pd.DataFrame):
    """(text, taxon or None, fraction) for the members' most common BIN."""
    counts = Counter(b for b in members["bin_uri"] if b)
    if not counts:
        return "", None, np.nan
    top = max(counts.values())
    tied = sorted(b for b, n in counts.items() if n == top)
    if len(tied) > 1:
        return f"most common BINs tied ({'; '.join(tied)})", None, np.nan
    in_bin = qc_passed[qc_passed["bin_uri"] == tied[0]]
    species = [lab for lab in in_bin["label"] if lab and label_level(lab) == "species"]
    if not species:
        return f"{tied[0]}: no species-level labels", None, np.nan
    lab_counts = Counter(species).most_common()
    n_top = lab_counts[0][1]
    leaders = sorted(lab for lab, n in lab_counts if n == n_top)
    if len(leaders) > 1:
        return (f"labels tied ({'; '.join(leaders)}; {n_top} each of {len(species)} "
                f"QC-passed records)"), None, np.nan
    frac = n_top / len(species)
    return (f"{leaders[0]} ({frac:.2f}; {n_top} of {len(species)} QC-passed records)",
            leaders[0], frac)


def _tree_sisters(tree_path, species_taxon: Dict[str, str], ids: List[str]) -> Dict[str, tuple]:
    """haplotype_id -> (sister taxon or '', support) from the tree, when it exists."""
    if not tree_path or not Path(tree_path).exists():
        return {}
    from Bio import Phylo
    tree = Phylo.read(str(tree_path), "newick")
    parent, tips_below = {}, {}
    for clade in tree.find_clades(order="postorder"):
        names = {clade.name} if clade.is_terminal() and clade.name else set()
        for child in clade.clades:
            parent[id(child)] = clade
            names |= tips_below[id(child)]
        tips_below[id(clade)] = names

    def support(clade):
        if clade.confidence is not None:
            return float(clade.confidence)
        try:
            return float(clade.name)
        except (TypeError, ValueError):
            return np.nan

    out = {}
    wanted = set(ids)
    for tip in tree.get_terminals():
        if tip.name not in wanted:
            continue
        node = tip
        while id(node) in parent:
            node = parent[id(node)]
            others = [t for t in tips_below[id(node)] if t != tip.name and t in species_taxon]
            if others:
                taxa = {species_taxon[t] for t in others}
                out[tip.name] = ((taxa.pop(), support(node)) if len(taxa) == 1 else ("", np.nan))
                break
    return out


def _suggest(dist_taxon, bin_taxon, tree_taxon):
    lines = [(name, t) for name, t in (("distance", dist_taxon), ("BIN", bin_taxon),
                                       ("tree", tree_taxon)) if t]
    named = {t for _, t in lines}
    if len(named) == 1 and len(lines) >= 2:
        return lines[0][1], ", ".join(name for name, _ in lines) + " agree"
    if len(named) > 1:
        return "", "evidence disagrees: " + "; ".join(f"{name} {t}" for name, t in lines)
    if lines:
        return "", f"insufficient evidence ({lines[0][0]} {lines[0][1]} only)"
    return "", "insufficient evidence (no line names a taxon)"


def compute_evidence(taxonomy: pd.DataFrame, records: pd.DataFrame,
                     sequences: Dict[str, str], tree_path: Optional[Path] = None) -> pd.DataFrame:
    """
    Evidence for each haplotype's taxon: haplotype_id + EVIDENCE_COLUMNS.

    taxonomy has haplotype_id, assigned_sp, assignment_level (and
    assignment_source); records has processid, haplotype_id, label and
    optionally qc_pass, bin_uri, inst, insdc_acs (missing text columns are
    treated as empty; without qc_pass every record counts as QC-passed).
    sequences maps haplotype_id to its aligned core sequence. tree_path is a
    Newick tree with haplotype IDs as tips (no tree evidence when None or missing).
    """
    tax = taxonomy.copy()
    if "assignment_source" not in tax.columns:
        tax["assignment_source"] = "labels"
    ids = [h for h in tax["haplotype_id"] if isinstance(h, str)]
    taxon = dict(zip(tax["haplotype_id"], tax["assigned_sp"].astype(object).where(tax["assigned_sp"].notna(), "")))
    level = dict(zip(tax["haplotype_id"], tax["assignment_level"]))
    source = dict(zip(tax["haplotype_id"], tax["assignment_source"]))
    species_ids = [h for h in ids if level.get(h) == "species" and taxon.get(h) and h in sequences]
    species_taxon = {h: taxon[h] for h in species_ids}

    recs = pd.DataFrame({c: _text(records, c) for c in
                         ("processid", "haplotype_id", "label", "bin_uri", "inst", "insdc_acs")})
    if "qc_pass" in records.columns:
        passed = records["qc_pass"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])
    else:
        passed = pd.Series(True, index=records.index)
    qc_passed = recs[passed.values]
    members_by_hap = {h: g for h, g in recs[recs["haplotype_id"] != ""].groupby("haplotype_id")}
    empty_members = recs.iloc[0:0]

    ref_seqs = [sequences[h] for h in species_ids]
    spread_cache: Dict[tuple, float] = {}

    def taxon_spread(t: str, exclude: str) -> float:
        group = tuple(h for h in species_ids if species_taxon[h] == t and h != exclude)
        if len(group) < 2:
            return np.nan
        if group not in spread_cache:
            seqs = [sequences[h] for h in group]
            spread_cache[group] = float(max(distances_to_references(s, seqs).max() for s in seqs))
        return spread_cache[group]

    sisters = _tree_sisters(tree_path, species_taxon, ids)
    rows = []
    for h in ids:
        row = {"haplotype_id": h}
        others = [i for i, r in enumerate(species_ids) if r != h]
        if h in sequences and others:
            d = distances_to_references(sequences[h], [ref_seqs[i] for i in others])
            dmin = float(d.min())
            nearest = [species_ids[others[k]] for k in np.flatnonzero(d <= dmin + _TOL)]
            near_taxa = list(dict.fromkeys(species_taxon[n] for n in nearest))
            row["nearest_species_haplotype"] = ";".join(nearest)
            row["nearest_species_taxon"] = ";".join(near_taxa)
            row["nearest_species_distance"] = round(dmin, 6)
            spread = taxon_spread(near_taxa[0], h) if len(near_taxa) == 1 else np.nan
            row["taxon_max_within_distance"] = round(spread, 6) if not np.isnan(spread) else np.nan
            row["within_taxon_spread"] = (bool(dmin <= spread + _TOL) if not np.isnan(spread)
                                          else np.nan)
        else:
            row.update(nearest_species_haplotype="", nearest_species_taxon="",
                       nearest_species_distance=np.nan, taxon_max_within_distance=np.nan,
                       within_taxon_spread=np.nan)

        members = members_by_hap.get(h, empty_members)
        row["member_bins"] = _member_bins(members) if len(members) else ""
        row["bin_taxon"], bin_major, bin_frac = _bin_taxon(members, qc_passed)
        row["label_sources"] = _label_sources(members) if len(members) else ""
        row["tree_sister_taxon"], row["tree_sister_support"] = sisters.get(h, ("", np.nan))

        if level.get(h) == "species" or source.get(h) == "curated":
            row["suggested_sp"], row["suggestion_basis"] = "", ""
        else:
            dist_taxon = (row["nearest_species_taxon"] if row["within_taxon_spread"] is True else "")
            bin_line = bin_major if bin_major and bin_frac >= BIN_MAJORITY - _TOL else ""
            row["suggested_sp"], row["suggestion_basis"] = _suggest(
                dist_taxon, bin_line, row["tree_sister_taxon"])
        rows.append(row)
    out = pd.DataFrame(rows, columns=["haplotype_id"] + EVIDENCE_COLUMNS)
    out["within_taxon_spread"] = out["within_taxon_spread"].astype(object)
    return out
