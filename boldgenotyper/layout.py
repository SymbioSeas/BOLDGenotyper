"""
Assemble the v2 run folder (P2).

The pipeline steps write their files into a staging folder; this module builds
the record and haplotype tables from in-memory data and places files at their
v2 paths.

build_records_table()        records.csv: one row per unique input record
build_haplotypes_table()     haplotypes.csv: one row per haplotype
build_haplotype_by_region()  geography/haplotype_by_region.csv
plan_moves(), apply_moves()  staged files -> v2 paths, following FILE_MAP
"""

import logging
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import numpy as np
import pandas as pd

from . import quality_control, taxon_curation
from .taxon_evidence import EVIDENCE_COLUMNS

logger = logging.getLogger(__name__)

# BOLD fields carried into records.csv, in this order (those present in the input).
# Other input columns are recovered by joining on processid with the input TSV.
KEY_BOLD_FIELDS = [
    "sampleid", "bin_uri", "insdc_acs", "marker_code", "inst",
    "phylum", "class", "order", "family", "genus", "species", "subspecies",
    "identification", "identification_method", "identification_rank", "identified_by",
    "taxonomy_notes", "sex", "life_stage", "reproduction", "notes", "short_note",
    "collection_notes", "collection_date_start", "collection_date_end",
    "country/ocean", "country_iso", "province/state", "region", "site",
    "coord", "coord_source", "coord_accuracy", "depth", "elev",
    "realm", "biome", "ecoregion", "habitat", "nuc_basecount",
]

_RECORD_FRONT = [
    "processid", "status", "status_note", "haplotype_id", "haplotype_label", "assigned_sp",
    "label", "label_conflict", "nearest_haplotype", "distance_to_nearest",
    "qc_pass", "qc_fail_reason", "sequence_length", "depositor_keyword_hits",
    "lat", "lon", "coord_status",
]
_RECORD_REGION = ["region_confidence", "distance_to_region_boundary_km"]

HAPLOTYPE_COLUMNS = [
    "haplotype_id", "haplotype_label", "seq_md5", "length", "n_members", "n_exact", "n_corrected",
    "is_singleton", "assigned_sp", "assignment_level", "assignment_notes", "assignment_source",
    "labels", "n_labeled", "n_missing_label", "n_species_labels", "majority_label", "n_majority",
    "majority_fraction", "n_label_conflicts", "nearest_haplotype", "nearest_distance",
    "is_suspect", "suspect_reasons", "min_orf_coverage", "max_internal_stops", "in_tree",
    "n_with_coords", "n_with_region", "n_with_country",
] + EVIDENCE_COLUMNS


def record_columns(region_col: str, input_columns: Iterable[str]) -> List[str]:
    """Column order of records.csv for a region column and the input's columns."""
    present = set(input_columns)
    return (_RECORD_FRONT + [region_col] + _RECORD_REGION
            + [c for c in KEY_BOLD_FIELDS if c in present and c != region_col])


def _has_text(series: pd.Series) -> pd.Series:
    return series.notna() & (series.astype(str).str.strip() != "")


def build_records_table(
    records: pd.DataFrame,
    status: pd.DataFrame,
    qc: pd.DataFrame,
    taxonomy: pd.DataFrame,
    region_col: str,
    region_confidence: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    One row per unique input record, in input order (records.csv).

    Parameters
    ----------
    records : unique input records with BOLD fields and lat, lon, coord_status
        and the region column
    status : record_status.build_record_status output
    qc : processid, qc_pass, qc_fail_reason, raw_length (records with a sequence)
    taxonomy : haplotype_id, assigned_sp, haplotype_sp
    region_col : name of the region column (e.g. 'ocean_basin')
    region_confidence : optional processid, basin_confidence, distance_to_boundary_km
    """
    out = records[["processid"]].merge(status, on="processid", how="left", validate="one_to_one")
    qc_cols = [c for c in ("processid", "qc_pass", "qc_fail_reason", "raw_length") if c in qc.columns]
    out = out.merge(qc[qc_cols], on="processid", how="left", validate="one_to_one")
    out = out.rename(columns={"raw_length": "sequence_length"})
    tax = taxonomy[["haplotype_id", "assigned_sp", "haplotype_sp"]].rename(
        columns={"haplotype_sp": "haplotype_label"})
    out = out.merge(tax, on="haplotype_id", how="left", validate="many_to_one")

    labeled = quality_control.add_label_columns(
        records.assign(haplotype_id=out["haplotype_id"].values),
        group_col="haplotype_id", notes_col="notes",
        curated_taxa=taxon_curation.curated_taxa(taxonomy))
    out["label"] = labeled["label"].values
    out["label_conflict"] = labeled["label_conflict"].values
    out["depositor_keyword_hits"] = labeled["depositor_keyword_hits"].values

    for col in ("lat", "lon", "coord_status", region_col):
        out[col] = records[col].values if col in records.columns else np.nan
    out["region_confidence"] = np.nan
    out["distance_to_region_boundary_km"] = np.nan
    if region_confidence is not None and len(region_confidence):
        conf = region_confidence.set_index("processid")
        if "basin_confidence" in conf.columns:
            out["region_confidence"] = out["processid"].map(conf["basin_confidence"])
        if "distance_to_boundary_km" in conf.columns:
            out["distance_to_region_boundary_km"] = out["processid"].map(conf["distance_to_boundary_km"])

    columns = record_columns(region_col, records.columns)
    for col in columns:
        if col not in out.columns:
            out[col] = records[col].values
    return out[columns]


def build_haplotypes_table(
    stats: pd.DataFrame,
    taxonomy: pd.DataFrame,
    records: pd.DataFrame,
    region_col: str,
    tree_tips: Optional[Set[str]] = None,
    evidence: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    One row per haplotype (haplotypes.csv), from the haplotype stats, the
    taxonomy table and records.csv. tree_tips is None when no tree was built.
    evidence (taxon_evidence.compute_evidence) fills the evidence columns;
    they are empty without it.
    """
    h = stats.rename(columns={"core_length": "length",
                              "nearest_neighbor_distance": "nearest_distance"}).copy()
    tax = taxonomy.reindex(columns=["haplotype_id", "assigned_sp", "assignment_level",
                                    "assignment_notes", "assignment_source", "haplotype_sp"])
    tax = tax.rename(columns={"haplotype_sp": "haplotype_label"})
    h = h.merge(tax, on="haplotype_id", how="left", validate="one_to_one")

    members = records[records["haplotype_id"].notna()]
    summary = quality_control.summarize_haplotype_labels(members).drop(columns=["n_members"])
    h = h.merge(summary, on="haplotype_id", how="left", validate="one_to_one")

    grouped = members.groupby("haplotype_id")
    h["n_with_coords"] = h["haplotype_id"].map(
        grouped["coord_status"].apply(lambda s: int((s == "ok").sum()))).fillna(0).astype(int)
    h["n_with_region"] = h["haplotype_id"].map(
        grouped[region_col].apply(lambda s: int(_has_text(s).sum()))).fillna(0).astype(int)
    if "country/ocean" in members.columns:
        h["n_with_country"] = h["haplotype_id"].map(
            grouped["country/ocean"].apply(lambda s: int(_has_text(s).sum()))).fillna(0).astype(int)
    else:
        h["n_with_country"] = 0
    h["in_tree"] = (h["haplotype_id"].isin(tree_tips) if tree_tips is not None
                    else pd.Series(np.nan, index=h.index, dtype=object))
    if evidence is not None and len(evidence):
        ev_cols = ["haplotype_id"] + [c for c in EVIDENCE_COLUMNS if c in evidence.columns]
        h = h.merge(evidence[ev_cols], on="haplotype_id", how="left", validate="one_to_one")

    for col in HAPLOTYPE_COLUMNS:
        if col not in h.columns:
            h[col] = np.nan
    return h[HAPLOTYPE_COLUMNS].sort_values("haplotype_id").reset_index(drop=True)


def build_haplotype_by_region(records: pd.DataFrame, region_col: str) -> pd.DataFrame:
    """Records per haplotype and region, with each haplotype's share of the region."""
    members = records[records["haplotype_id"].notna() & _has_text(records[region_col])]
    counts = (members.groupby(["haplotype_id", "haplotype_label", "assigned_sp", region_col],
                              dropna=False)
              .size().rename("n_records").reset_index())
    totals = counts.groupby(region_col)["n_records"].transform("sum")
    counts["pct_of_region"] = (counts["n_records"] / totals * 100).round(2)
    return counts.sort_values([region_col, "haplotype_id"]).reset_index(drop=True)


# Staged file (glob relative to the staging folder; {org} is the organism) ->
# v2 path, and the run option that places it. The first matching entry claims a
# file. A destination ending in "/" keeps the file name without the "{org}_"
# prefix. Staged files no entry places stay in staging, which is removed with
# intermediate/ (they are duplicates or folded into records.csv/haplotypes.csv).
FILE_MAP = [
    ("haplotypes/{org}_haplotypes.fasta", "haplotypes.fasta", "default"),
    ("quality_control/label_conflicts.csv", "label_conflicts.csv", "default"),
    ("metadata_summary/{org}_metadata_summary.csv", "metadata_coverage.csv", "default"),
    ("divergence_analysis/{org}_consensus_aligned.fasta", "divergence/haplotypes_aligned.fasta", "default"),
    ("divergence_analysis/pairwise_p_distance_by_id.csv", "divergence/pairwise_p_distance.csv", "default"),
    ("phylogenetic/{org}_tree.nwk", "phylogeny/tree.nwk", "default"),
    ("phylogenetic/{org}_tree_relabeled.nwk", "phylogeny/tree_labeled.nwk", "default"),
    ("phylogenetic/{org}_excluded_tips.txt", "phylogeny/excluded_tips.txt", "default"),

    ("metadata_summary/{org}_*_tabulation.csv", "extended/metadata/", "extended"),
    ("metadata_summary/{org}_haplotype_first_sampled.csv", "extended/metadata/", "extended"),
    ("metadata_summary/{org}_temporal_*.csv", "extended/metadata/", "extended"),
    ("divergence_analysis/label_based_gap.*", "extended/divergence/", "extended"),
    ("divergence_analysis/divergence_heatmap.pdf", "extended/divergence/", "extended"),
    ("divergence_analysis/divergence_summary.csv", "extended/divergence/species_pair_divergence.csv", "extended"),
    ("species_analysis/species_divergence_summary.csv", "extended/divergence/", "extended"),
    ("species_analysis/species_divergence_matrices/*.csv", "extended/divergence/within_species/", "extended"),
    ("phylogenetic/{org}_MSA_chunk*.pdf", "extended/divergence/msa/", "extended"),
    ("species_analysis/species_summary.csv", "extended/species/", "extended"),
    ("species_analysis/species_diversity.csv", "extended/species/", "extended"),
    ("species_analysis/species_basin_distribution.csv", "extended/species/species_by_region.csv", "extended"),
    ("species_analysis/species_country_distribution.csv", "extended/species/species_by_country.csv", "extended"),
    ("quality_control/depositor_keyword_hits.csv", "extended/labels/", "extended"),
    ("quality_control/*.pdf", "extended/labels/", "extended"),
]


def plan_moves(staging: Path, organism: str, extended: bool = False) -> Dict[Path, Path]:
    """
    Map each staged file that the run places to its v2 path (relative to the
    run folder), following FILE_MAP.
    """
    staging = Path(staging)
    enabled = {"default": True, "extended": extended}
    prefix = f"{organism}_"
    claimed: Set[Path] = set()
    moves: Dict[Path, Path] = {}
    for pattern, dest, category in FILE_MAP:
        for src in sorted(staging.glob(pattern.format(org=organism))):
            if not src.is_file() or src in claimed:
                continue
            claimed.add(src)
            if not enabled[category]:
                continue
            if dest.endswith("/"):
                name = src.name[len(prefix):] if src.name.startswith(prefix) else src.name
                moves[src] = Path(dest) / name
            else:
                moves[src] = Path(dest)
    return moves


def apply_moves(moves: Dict[Path, Path], output_dir: Path) -> List[Path]:
    """Move staged files to their v2 paths; return the destination paths."""
    placed = []
    for src, rel in moves.items():
        dst = Path(output_dir) / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
        placed.append(dst)
    return placed


def write_per_species_records(records: pd.DataFrame, dest_dir: Path, min_haplotypes: int = 2) -> List[Path]:
    """
    records.csv rows of each assigned species with at least min_haplotypes
    haplotypes, one file per species (extended/species/per_species/).
    """
    members = records[records["haplotype_id"].notna() & _has_text(records["assigned_sp"])]
    written = []
    for species, group in members.groupby("assigned_sp"):
        if group["haplotype_id"].nunique() < min_haplotypes:
            continue
        dest_dir.mkdir(parents=True, exist_ok=True)
        safe = str(species).replace(" ", "_").replace("/", "_")
        path = dest_dir / f"{safe}_records.csv"
        group.to_csv(path, index=False)
        written.append(path)
    return written
