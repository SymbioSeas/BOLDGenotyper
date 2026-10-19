"""
Generate the run folder's README.md (P2): every file, and every column of
every CSV, with a definition.

FILE_DOCS maps file patterns (relative to the run folder) to a description and
the definitions of that file's columns; COLUMNS holds definitions shared by
several files. write_run_readme() returns the files and columns it could not
define, so tests can require that list to be empty.
"""

import csv
import fnmatch
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .layout import KEY_BOLD_FIELDS
from .record_status import STATUSES

logger = logging.getLogger(__name__)

REGION = "<region>"   # replaced by the run's region column (e.g. ocean_basin)

STATUS_TEXT = {
    "assigned": "core sequence identical to a retained haplotype",
    "assigned_corrected": "error-filtered singleton assigned to its nearest retained haplotype, "
                          "which is within the error-filter threshold (--singleton-distance)",
    "excluded_error_singleton": "error-filtered singleton left unassigned: equidistant from two or more "
                                "retained haplotypes, no retained haplotype within the threshold, or "
                                "correction disabled (--no-assign-corrected)",
    "excluded_core_coverage": "ungapped core region shorter than the per-sequence cutoff",
    "failed_qc": "failed pre-alignment quality control (length, N fraction or reading frame; see qc_fail_reason)",
    "no_sequence": "no sequence in the input",
}

COORD_STATUS_TEXT = {
    "ok": "parseable coordinates of geographic quality; the only records assigned a region and mapped",
    "missing": "no coordinates",
    "unparseable": "coordinates that could not be read",
    "centroid": "country or province centroid (lat/lon kept; not assigned a region or mapped)",
    "zero": "0, 0 (lat/lon kept; not assigned a region or mapped)",
    "out_of_range": "latitude or longitude out of range",
}

# Definitions shared by several files
COLUMNS = {
    "processid": "BOLD process ID (bare, as in the input)",
    "haplotype_id": "haplotype ID (H0001, H0002, ... in order of member count, then seq_md5)",
    "haplotype_label": "assigned taxon and haplotype ID, e.g. 'Genus species H0001'",
    "haplotype_sp": "assigned taxon and haplotype ID, e.g. 'Genus species H0001' (haplotype_label "
                    "in records.csv and haplotypes.csv)",
    "assigned_sp": "taxon assigned to the haplotype (see Taxonomy below)",
    "label": "the record's species label, or its genus when the species is missing",
    "seq_md5": "md5 of the haplotype's core sequence; stable across runs on the same data",
    "lat": "latitude parsed from coord (kept for centroid and zero coordinates)",
    "lon": "longitude parsed from coord",
    REGION: "region assigned from geographic-quality coordinates; empty when coord_status is not ok "
            "or the point is outside every region",
    "depositor_keyword_hits": "keywords from the depositor's notes, matched as whole words, ';'-joined",
    "majority_label": "most common species-level label among the haplotype's members; empty on a tie",
    "n_majority": "members carrying majority_label",
    "n_labeled": "members with a species or genus label",
    "majority_fraction": "n_majority / counted labels, where counted labels are the members with a "
                         "species-level label plus genus-only or open-nomenclature labels of another "
                         "genus (members without a label, and genus-only or open labels of the "
                         "majority's genus, are left out); empty when no member has a species-level label",
    "n_members": "records assigned to the haplotype (status assigned or assigned_corrected)",
    "country/ocean": "BOLD country/ocean field, copied from the input",
}

_STATUS_LIST = ", ".join(STATUSES)

FILE_DOCS = [
    # --- run folder -------------------------------------------------------
    {"pattern": "README.md", "text": "This file."},
    {"pattern": "run_provenance.json",
     "text": "Full configuration, the values the run applied (core length, cutoffs, tree tips), "
             "tool versions, random seeds, input file md5 and row count, data notes, warnings and "
             "run status. boldgenotyper-query reports its md5 as reference_run."},
    {"pattern": "*_pipeline.log", "text": "Run log; it ends with a summary of warnings and the run status."},
    {"pattern": "report.html",
     "text": "Summary report. Open it from this folder: figures are linked, not embedded."},
    {"pattern": "summary.csv",
     "text": "One row: record counts at each step, assignment percentages and species counts.",
     "columns": {
         "n_input": "data rows in the input TSV",
         "n_other_marker_rows": "rows of markers other than the target marker (--marker, "
                                "default COI-5P), removed before duplicates",
         "n_malformed_rows": "data rows whose number of fields differs from the header's, "
                             "dropped with --allow-malformed-rows (otherwise the run stops)",
         "n_duplicate_rows": "rows repeating a processid of the target marker (first kept)",
         "pct_corrected": "assigned_corrected as a percentage of all assigned records (assigned + "
                          "assigned_corrected)",
         "n_singleton_label_absent": "error-filtered singleton records whose species label no "
                                     "assigned record carries (records.csv label_in_retained False): "
                                     "the filter may have removed a species",
         "n_haplotype_only_sources": "sources in sources.csv flagged haplotype_only",
         "n_unique": "unique records of the target marker (rows of other markers and duplicate "
                     "processids removed; first kept)",
         "n_no_sequence": "records with status no_sequence",
         "n_failed_qc": "records with status failed_qc",
         "n_qc_pass": "records that passed QC (n_unique - n_no_sequence - n_failed_qc); the "
                      "denominator of the percentages",
         "n_excluded_core": "records with status excluded_core_coverage",
         "n_core_ok": "QC-passed records with enough core coverage",
         "n_assigned": "records with status assigned",
         "n_assigned_corrected": "records with status assigned_corrected",
         "n_excluded_error": "records with status excluded_error_singleton",
         "n_excluded_error_tied": "of these, equidistant from two or more retained haplotypes",
         "n_excluded_error_far": "of these, with no retained haplotype within --singleton-distance",
         "n_haplotypes": "retained haplotypes",
         "pct_assigned": "(n_assigned + n_assigned_corrected) / n_qc_pass x 100",
         "pct_assigned_exact": "n_assigned / n_qc_pass x 100",
         "n_species_labels_input": "distinct non-empty species labels among unique input records",
         "n_species_labels_species_level": "of these, labels that name a species (not open "
                                           "nomenclature such as 'sp.' or 'cf.', not genus-only)",
         "n_species_labels_qc_passed": "distinct non-empty species labels among QC-passed records",
         "n_species_haplotype_assigned": "distinct species-level assigned_sp across haplotypes",
         "n_species_confident": "distinct species-level assigned_sp among haplotypes whose "
                                "majority_fraction is at least the majority threshold",
     }},
    {"pattern": "records.csv",
     "text": "One row per unique input record, in input order, including records that failed QC "
             "or have no sequence. Other BOLD columns of the input can be joined on processid "
             "(the input's md5 is in run_provenance.json).",
     "columns": {
         "status": f"assignment status (one of {_STATUS_LIST}; see Record statuses)",
         "status_note": "why the record has its status (e.g. the QC failure or the tied haplotypes)",
         "haplotype_id": "haplotype the record is assigned to; empty unless status is assigned or "
                         "assigned_corrected",
         "nearest_haplotype": "for error-filtered singletons, the nearest retained haplotype(s)",
         "distance_to_nearest": "p-distance to nearest_haplotype",
         "label_in_retained": "for error-filtered singletons (assigned_corrected, "
                              "excluded_error_singleton) with a species-level label: True when an "
                              "assigned record carries the same label, False when none does; empty "
                              "otherwise",
         "label_mapped_from": "the deposited label, when the run's label map (--label-map) replaced "
                              "it; empty otherwise",
         "label_conflict": "True when the record's species-level label differs from its haplotype's "
                           "majority label, or the haplotype's top labels are tied (see Labels)",
         "qc_pass": "passed pre-alignment QC; empty for records without a sequence",
         "qc_fail_reason": "the QC check the record failed",
         "sequence_length": "length of the input sequence (bp)",
         "coord_status": "coordinate quality (see Coordinates)",
         "region_confidence": "high (>50 km from a region boundary), medium (10-50 km), low (<10 km) "
                              "or none (no region)",
         "distance_to_region_boundary_km": "distance from the point to the nearest region boundary",
     }},
    {"pattern": "records_qc_passed.fasta",
     "text": "Oriented sequences of the QC-passed records, headed by processid; "
             "boldgenotyper-query --search-records searches them."},
    {"pattern": "haplotypes.fasta",
     "text": "Haplotype core sequences, ungapped. Header: '>H0001 seq_md5=... n_members=...'."},
    {"pattern": "haplotypes.csv",
     "text": "One row per haplotype.",
     "columns": {
         "length": "core sequence length (bp)",
         "n_exact": "members whose core sequence is identical to the haplotype",
         "n_corrected": "members with status assigned_corrected",
         "is_singleton": "n_members == 1",
         "assignment_level": "species, genus or unassigned (see Taxonomy)",
         "assignment_notes": "how assigned_sp was chosen: the majority and its counted labels, or "
                             "why it fell back to genus; for a curated taxon, the curation "
                             "reason followed by the label rule's note; for an evidence taxon, "
                             "the lines that agree followed by the label rule's note",
         "assignment_source": "labels (the label rule), evidence (the evidence named the "
                              "haplotype's leading label; see Taxonomy) or curated (a row of the "
                              "--taxon-curation file)",
         "labels": "every member label with its count, most common first; "
                   "'(no species label)' counts members without one",
         "n_missing_label": "members without a species or genus",
         "n_species_labels": "distinct species-level labels among members",
         "n_label_conflicts": "members with label_conflict True",
         "nearest_haplotype": "nearest other haplotype(s) by p-distance on the core alignment; "
                              "';'-joined when tied",
         "nearest_distance": "p-distance to nearest_haplotype",
         "is_suspect": "flagged as suspect (distance or reading-frame checks); suspect haplotypes "
                       "are left out of the tree",
         "suspect_reasons": "why the haplotype is suspect",
         "min_orf_coverage": "lowest open-reading-frame coverage among members",
         "max_internal_stops": "most internal stop codons among members",
         "in_tree": "tip in phylogeny/tree.nwk; empty when no tree was built",
         "n_with_coords": "members with coord_status ok",
         "n_with_region": f"members with a {REGION}",
         "n_with_country": "members with a country/ocean value",
         "nearest_species_haplotype": "nearest other haplotype assigned at species level "
                                      "(';'-joined when tied); evidence, see Taxonomy",
         "nearest_species_taxon": "its assigned_sp (';'-joined when the tied haplotypes differ)",
         "nearest_species_distance": "p-distance to nearest_species_haplotype",
         "taxon_max_within_distance": "largest p-distance between the other species-level "
                                      "haplotypes of nearest_species_taxon; empty when it has "
                                      "fewer than two",
         "within_taxon_spread": "nearest_species_distance <= taxon_max_within_distance",
         "member_bins": "members' BOLD BINs with counts ('none': no BIN)",
         "bin_taxon": "for the members' most common BIN: the most common species-level label "
                      "among this run's QC-passed records in that BIN, its fraction and counts",
         "label_sources": "per member label: records, BOLD projects (processid prefix; GenBank "
                          "accession prefix for records mined from GenBank) and institutions. "
                          "Approximate: one study can span several projects",
         "tree_sister_taxon": "after collapsing branches of phylogeny/tree.nwk with support "
                              "below 0.70, walking from the haplotype's tip toward the root: "
                              "the first clade holding another species-level haplotype; their "
                              "shared taxon. Empty when they differ, when only the root holds "
                              "one, or when no tree was built",
         "tree_sister_support": "that clade's support value (FastTree SH-like local support, "
                                "0.70 or more)",
         "suggested_sp": "for haplotypes the label rule left below species level and not "
                         "curated: a taxon named by at least two evidence lines and contradicted "
                         "by none (see Taxonomy). Applied when it is the haplotype's leading "
                         "label in the label rule's genus; otherwise listed in "
                         "taxon_curation_suggestions.csv",
         "suggestion_basis": "the lines that agree, or why there is no suggestion",
     }},
    {"pattern": "sources.csv",
     "text": "One row per source of records seen by the singleton error filter (assigned, "
             "assigned_corrected, excluded_error_singleton): the BOLD project (processid prefix), "
             "or the GenBank accession prefix for records mined from GenBank. A source that "
             "deposits one sequence per haplotype makes true rare haplotypes look like errors to "
             "the filter; haplotype_only flags it (at least 10 records, at least 90% distinct "
             "sequences). Check the error filter's results for such datasets.",
     "columns": {
         "source": "BOLD project code or GenBank accession prefix",
         "n_records": "records from the source seen by the error filter",
         "n_distinct_sequences": "distinct deposited sequences among them",
         "distinct_fraction": "n_distinct_sequences / n_records",
         "n_assigned": "records with status assigned",
         "n_corrected": "records with status assigned_corrected",
         "n_excluded_error": "records with status excluded_error_singleton",
         "haplotype_only": "n_records >= 10 and distinct_fraction >= 0.90",
     }},
    {"pattern": "label_variants.csv",
     "text": "Pairs of species labels in one genus that may be one name written differently: "
             "epithets one or two letters apart, a changed gender ending, or one label that is "
             "the other plus more words. Listed only; nothing is merged. Valid species can be "
             "two letters apart and share haplotypes, so check each pair, and list the ones to "
             "merge in a label map (--label-map; see Taxonomy).",
     "columns": {
         "label_a": "the more frequent label (as deposited)",
         "label_b": "the other label (as deposited)",
         "n_a": "records with label_a",
         "n_b": "records with label_b",
         "edit_distance": "letters that differ between the epithets",
         "ending_change": "the epithets differ only in a gender ending (e.g. -us/-a)",
         "qualifier_form": "one label is the other plus more words (e.g. 'complex')",
         "n_shared_haplotypes": "haplotypes with members carrying both labels",
         "in_label_map": "either label is a 'from' in the run's label map",
     }},
    {"pattern": "label_conflicts.csv",
     "text": "One row per record whose label conflicts with its haplotype's majority label, or "
             "with its haplotype's curated or evidence taxon (see Taxonomy), with the evidence. A "
             "conflict is not a verdict on which label is right.",
     "columns": {
         "haplotype_taxon": "the haplotype's assigned_sp",
         "assignment_source": "labels, evidence or curated; unless labels, the record's label is "
                              "compared with haplotype_taxon rather than majority_label",
         "nearest_haplotype_with_label": "nearest other haplotype with at least one member carrying "
                                         "the record's label; empty when none",
         "distance_to_nearest_with_label": "p-distance from the record's haplotype to that haplotype",
         "depositor_note": "the depositor's notes field",
     }},
    {"pattern": "taxon_curation_suggestions.csv",
     "text": "Suggested taxa the run did not apply (see Taxonomy), as a curation file: delete "
             "the rows you reject, fill curator and date, and pass the file to --taxon-curation. "
             "Written only when there is such a suggestion.",
     "columns": {
         "seq_md5": "the haplotype's seq_md5 (the curation key)",
         "sequence": "its core sequence, ungapped (a key for another dataset's haplotypes)",
         "taxon": "the suggested taxon (suggested_sp in haplotypes.csv)",
         "reason": "the evidence lines that name it, and the haplotype's member labels",
         "curator": "empty: yours to fill",
         "date": "empty: yours to fill",
         "haplotype_id": "the haplotype; for reading (--taxon-curation ignores it)",
         "assigned_sp": "its taxon in this run",
         "why_not_applied": "not the leading label, leading label tied, another genus than "
                            "the label rule's, no species-level label, or evidence taxa off "
                            "(--no-evidence-taxa)",
     }},
    {"pattern": "metadata_coverage.csv",
     "text": "Coverage of BOLD metadata fields after normalization (see Metadata). Default "
             "fields: sex, life_stage, reproduction, habitat, country/ocean, province/state; "
             "others (e.g. realm, biome, ecoregion, which BOLD derives from the coordinates "
             "using terrestrial ecoregions) only when named in --metadata-fields.",
     "columns": {
         "field": "BOLD field",
         "n_total": "QC-passed records",
         "n_with_value": "records with a value",
         "pct_coverage": "n_with_value / n_total",
         "n_unique": "distinct values",
         "top_values": "most common values with counts",
         "warning": "note on the field (e.g. value variants that may be the same value)",
     }},
    # --- divergence, phylogeny, geography ----------------------------------
    {"pattern": "divergence/haplotypes_aligned.fasta",
     "text": "Haplotype sequences aligned with MAFFT; headed by haplotype label."},
    {"pattern": "divergence/pairwise_p_distance.csv",
     "text": "Uncorrected p-distance between every pair of haplotypes: differing sites / sites "
             "compared, skipping sites with a gap or N in either sequence.",
     "matrix": "Square matrix: the first column and the header are haplotype_id; taxa are "
               "in haplotypes.csv."},
    {"pattern": "phylogeny/tree.nwk",
     "text": "FastTree tree of the haplotypes, tips named by haplotype ID; support values are "
             "FastTree local support. Rooted on the outgroup when one was given, otherwise at "
             "the midpoint of the longest path between two tips (run_provenance.json, "
             "applied.tree_rooting); support values stay on their branches when rooting."},
    {"pattern": "phylogeny/tree_labeled.nwk", "text": "The same tree with tips named by haplotype label."},
    {"pattern": "phylogeny/excluded_tips.txt",
     "text": "Haplotypes left out of the tree, with the reason (suspect, fewer than "
             "--tree-min-members members, or shorter than the minimum length)."},
    {"pattern": "phylogeny/tree.*",
     "text": "The full tree: every tip, labeled by haplotype; branch values are supports of 0.70 "
             "or more."},
    {"pattern": "phylogeny/tree_overview.*",
     "text": "Tree overview: the full tree for up to 40 haplotypes; above, each taxon whose "
             "haplotypes form a clade is collapsed to one tip (haplotypes and records in "
             "parentheses), and taxa that do not form a clade are left expanded and listed."},
    {"pattern": "phylogeny/tree_by_genus.pdf",
     "text": "One page per genus with two or more haplotypes: the clade at the most recent common "
             "ancestor of the genus's haplotypes (tips of other genera inside it included)."},
    {"pattern": "phylogeny/by_genus/*", "text": "The pages of tree_by_genus.pdf, one SVG per genus."},
    {"pattern": "*/manifest.json",
     "text": "The figures in this folder: for each figure, its overview files and its per-species "
             "(or per-genus) pages."},
    {"pattern": "geography/haplotype_by_region.csv",
     "text": "Records per haplotype and region (members with a region only).",
     "columns": {
         "n_records": "members of the haplotype in the region",
         "pct_of_region": "n_records / all members in the region x 100",
     }},
    {"pattern": "geography/distribution_map.*",
     "text": "Map of records with geographic-quality coordinates, colored by haplotype for up to "
             "40 haplotypes and by taxon above (12 haplotypes or 15 taxa colored; the rest gray "
             "as 'Other'); point size is records per site."},
    {"pattern": "geography/distribution_map_by_species.pdf",
     "text": "One page per taxon with mapped records, colored by haplotype (a haplotype keeps its "
             "color on every page)."},
    {"pattern": "geography/by_species/*", "text": "The pages of distribution_map_by_species.pdf, one SVG per taxon."},
    # --- extended/ ---------------------------------------------------------
    {"pattern": "extended/metadata/*_tabulation.csv",
     "text": "--extended. Records per haplotype and value of one metadata field (records with a value).",
     "columns": {
         "value": "normalized field value",
         "n_samples": "records with this haplotype and value",
         "pct_of_haplotype": "n_samples / the haplotype's records with a value x 100",
         "pct_of_value": "n_samples / records with this value x 100",
     }},
    {"pattern": "extended/metadata/haplotype_first_sampled.csv",
     "text": "--extended. First and last collection date among each haplotype's records with a "
             "readable date. These are sampling dates, not dates of origin or loss.",
     "columns": {
         "first_sampled": "earliest collection_date_start",
         "first_sampled_year": "year of first_sampled",
         "last_sampled": "latest collection_date_start",
         "last_sampled_year": "year of last_sampled",
         "years_between_first_and_last_sample": "last_sampled_year - first_sampled_year",
         "n_total_samples": "records with a readable collection date",
         "first_location": "country/ocean of the first-sampled record; empty when missing",
         "last_location": "country/ocean of the last-sampled record; empty when missing",
     }},
    {"pattern": "extended/metadata/temporal_trends.csv",
     "text": "--extended. Records per haplotype and collection year. The first line is a '#' "
             "comment (read with comment='#'): counts by year reflect sampling.",
     "columns": {
         "time_period": "collection year",
         "n_samples": "records of the haplotype collected that year",
         "pct_of_period": "n_samples / all records collected that year x 100",
     }},
    {"pattern": "extended/metadata/temporal_summary.csv",
     "text": "--extended. Records and haplotypes per collection year.",
     "columns": {
         "year": "collection year",
         "total_samples": "records collected that year",
         "n_haplotypes": "haplotypes among them",
     }},
    {"pattern": "extended/metadata/figures/*",
     "text": "--extended. Metadata figures (PDF and SVG): per field, records by haplotype or taxon "
             "stacked by value, one page per taxon with two or more haplotypes "
             "(<field>_by_haplotype.pdf) and records by region (<field>_by_region); records by "
             "collection year (heatmap; pages by haplotype), the sampling timeline, records by "
             "year and region, and metadata coverage. 'By region' views use the run's region "
             "column, or country/ocean when the run has fewer than two regions (the caption "
             "says which)."},
    {"pattern": "extended/divergence/label_based_gap.csv",
     "text": "--extended. Barcoding gap computed from record labels (label-conflict records left out); "
             "distances between records are the p-distances between their haplotypes.",
     "columns": {
         "n_records": "records carrying the label",
         "n_haplotypes": "haplotypes those records belong to",
         "max_within": "largest p-distance between haplotypes carrying the label",
         "min_between": "smallest p-distance from those haplotypes to a haplotype carrying another label",
         "nearest_other_label": "the label at min_between",
         "label_based_gap": "True when min_between > max_within",
     }},
    {"pattern": "extended/divergence/label_based_gap.pdf", "text": "--extended. max_within against min_between per label."},
    {"pattern": "extended/divergence/divergence_heatmap.pdf", "text": "--extended. Heatmap of pairwise_p_distance.csv."},
    {"pattern": "extended/divergence/species_pair_divergence.csv",
     "text": "--extended. p-distances between haplotypes grouped by assigned_sp.",
     "columns": {
         "comparison_type": "within_species or between_species",
         "species_pair": "the species, or the two species compared",
         "mean_divergence": "mean p-distance over the haplotype pairs",
         "min_divergence": "smallest p-distance",
         "max_divergence": "largest p-distance",
     }},
    {"pattern": "extended/divergence/species_divergence_summary.csv",
     "text": "--extended. For each assigned species with two or more haplotypes, p-distances among "
             "its haplotypes and to all other haplotypes.",
     "columns": {
         "species": "assigned_sp",
         "n_haplotypes": "haplotypes assigned to the species",
         "mean_within_species": "mean p-distance among them",
         "max_within_species": "largest",
         "min_within_species": "smallest",
         "mean_between_species": "mean p-distance to haplotypes of other species",
         "max_between_species": "largest",
         "min_between_species": "smallest",
     }},
    {"pattern": "extended/divergence/within_species/within_species_divergence_summary.csv",
     "text": "--extended. Summary of the per-species matrices.",
     "columns": {
         "species": "assigned_sp",
         "n_haplotypes": "haplotypes assigned to the species",
         "n_comparisons": "haplotype pairs",
         "mean_divergence": "mean p-distance",
         "max_divergence": "largest",
         "min_divergence": "smallest",
         "std_divergence": "standard deviation",
     }},
    {"pattern": "extended/divergence/within_species/*_divergence_matrix.csv",
     "text": "--extended. p-distances among the haplotypes of one assigned species.",
     "matrix": "Square matrix: the first column and the header are haplotype_id; taxa are "
               "in haplotypes.csv."},
    {"pattern": "extended/divergence/msa/MSA_chunk*.pdf",
     "text": "--extended. Haplotype alignment in tree order, in chunks of columns."},
    {"pattern": "extended/species/species_summary.csv",
     "text": "--extended. Records and haplotypes per species, for haplotypes whose primary species "
             "carries at least 70% of members.",
     "columns": {
         "species": "primary species",
         "n_samples": "records",
         "n_haplotypes": "haplotypes",
         "haplotype_list": "their haplotype IDs",
     }},
    {"pattern": "extended/species/species_diversity.csv",
     "text": "--extended. Haplotype diversity per species label; it depends on the labels and on "
             "which specimens were sampled.",
     "columns": {
         "species": "primary species",
         "n_samples": "records",
         "n_haplotypes": "haplotypes",
         "haplotype_diversity": "Nei's h = n/(n-1) x (1 - sum of squared haplotype frequencies)",
         "dominant_haplotype": "most common haplotype",
         "dominant_haplotype_freq": "its share of the species' records",
         "note": "'label-based, sampling-dependent'",
     }},
    {"pattern": "extended/species/species_by_region.csv",
     "text": "--extended. Records per primary species and region.",
     "columns": {"primary_species": "primary species of the haplotype", "n_samples": "records"}},
    {"pattern": "extended/species/species_by_country.csv",
     "text": "--extended. Records per primary species and country/ocean.",
     "columns": {"primary_species": "primary species of the haplotype", "n_samples": "records"}},
    {"pattern": "extended/species/per_species/*_records.csv",
     "text": "--extended. The records.csv rows of one assigned species with two or more haplotypes.",
     "same_as": "records.csv"},
    {"pattern": "extended/geography/*",
     "text": "--extended. Records per region by haplotype or taxon, as shares of each region "
             "(region_abundance_relative) and counts (region_abundance_total), with one page per "
             "taxon with two or more haplotypes (PDF and SVG)."},
    {"pattern": "extended/labels/depositor_keyword_hits.csv",
     "text": "--extended. Records whose depositor notes contain a keyword (whole-word match).",
     "columns": {"depositor_note": "the depositor's notes field",
                 "keyword_hits": "keywords found, ';'-joined"}},
    {"pattern": "extended/labels/*.pdf",
     "text": "--extended. Label composition of haplotypes, and the distribution of majority_fraction."},
    # --- exports/ and plots/ -------------------------------------------------
    {"pattern": "exports/README.md", "text": "--export-format. How to use the exported files."},
    {"pattern": "exports/arlequin/*.arp", "text": "--export-format arlequin. Arlequin project file."},
    {"pattern": "exports/arlequin/populations.csv",
     "text": "--export-format arlequin. Populations (regions) in the Arlequin file; records without a "
             "region are left out of the Arlequin and DnaSP files.",
     "columns": {"population": "region", "n_samples": "records", "n_haplotypes": "haplotypes",
                 "dominant_species": "most common species label"}},
    {"pattern": "exports/dnasp/*.fas", "text": "--export-format dnasp. DnaSP FASTA, one sequence per record."},
    {"pattern": "exports/popart/*.nexus", "text": "--export-format popart. PopART NEXUS file."},
    {"pattern": "exports/popart/popart_traits.csv",
     "text": "--export-format popart. Traits table for PopART (every assigned record; records without "
             "a region are 'Unassigned').",
     "columns": {"taxon": "taxon name in the NEXUS file: haplotype ID and a running number (H0001_001)",
                 "haplotype": "haplotype ID", "species": "species label (spaces as underscores)",
                 REGION: "region, or 'Unassigned' for records without one"}},
    {"pattern": "exports/generic/haplotype_membership.csv",
     "text": "--export-format generic. Record-to-haplotype table.",
     "columns": {"species": "species label", "consensus_sequence": "the haplotype's sequence",
                 "sequence_length": "its length (bp)"}},
    {"pattern": "exports/generic/haplotype_summary.csv",
     "text": "--export-format generic. One row per haplotype.",
     "columns": {"sequence": "haplotype sequence", "sequence_length": "its length (bp)",
                 "n_samples": "members", "species": "most common species label",
                 "primary_population": "region with most members",
                 "n_populations": "regions with members"}},
    {"pattern": "plots/README.md", "text": "--plot-kit. How to regenerate the figures."},
    {"pattern": "plots/plot_config.yaml", "text": "--plot-kit. Plot settings the scripts read."},
    {"pattern": "plots/scripts/*.py",
     "text": "--plot-kit. Scripts that regenerate the map, bar charts and tree into plots/visualization/."},
    {"pattern": "plots/visualization/*", "text": "--plot-kit. Figures regenerated by the scripts."},
    {"pattern": "plots/data/distribution_bar_absolute.csv",
     "text": "--plot-kit. Data for the region bar chart (counts).",
     "columns": {"consensus_group": "haplotype ID (v1 column name, kept for the scripts)",
                 "consensus_group_sp": "haplotype label (v1 column name)",
                 "species": "assigned taxon", "n_samples": "records"}},
    {"pattern": "plots/data/distribution_bar_relative.csv",
     "text": "--plot-kit. Data for the region bar chart (shares).",
     "columns": {"consensus_group": "haplotype ID (v1 column name, kept for the scripts)",
                 "consensus_group_sp": "haplotype label (v1 column name)",
                 "species": "assigned taxon", "n_samples": "records",
                 "pct_of_basin": "n_samples / records in the region x 100"}},
    {"pattern": "plots/data/distribution_map.csv",
     "text": "--plot-kit. Data for the distribution map.",
     "columns": {"consensus_group": "haplotype ID (v1 column name, kept for the scripts)",
                 "consensus_group_sp": "haplotype label (v1 column name)",
                 "species": "assigned taxon",
                 "n_at_location": "records at the same coordinates"}},
    {"pattern": "plots/data/genotype_colors.csv",
     "text": "--plot-kit. Colors per haplotype.",
     "columns": {"genotype": "haplotype ID", "genotype_label": "haplotype label",
                 "color_hex": "color", "species": "assigned taxon", "n_samples": "members"}},
]


def _entry_for(rel: str) -> Optional[dict]:
    for entry in FILE_DOCS:
        if fnmatch.fnmatch(rel, entry["pattern"]):
            return entry
    return None


def _csv_header(path: Path) -> List[str]:
    with open(path, newline="") as f:
        for line in f:
            if not line.startswith("#"):
                return next(csv.reader([line]))
    return []


def _column_text(col: str, entry: dict, region_col: str) -> Optional[str]:
    own = entry.get("columns", {})
    if col in own:
        return own[col]
    if col == region_col:
        return own.get(REGION, COLUMNS[REGION])
    if col in COLUMNS:
        return COLUMNS[col]
    if col in KEY_BOLD_FIELDS:
        return "BOLD field, copied from the input"
    return None


def _sections(region_col: str, has_metadata_extras: bool) -> List[str]:
    lines = ["## Record statuses", "",
             "Every unique input record has one status (records.csv `status`):", ""]
    lines += [f"- `{s}`: {STATUS_TEXT[s]}" for s in STATUSES]
    lines += ["", "## Coordinates", "",
              "`coord_status` in records.csv:", ""]
    lines += [f"- `{k}`: {v}" for k, v in COORD_STATUS_TEXT.items()]
    lines += ["",
              f"Only records with `coord_status` ok are assigned a region (`{region_col}`) and mapped. "
              "Country- and province-centroid coordinates are kept in lat/lon but are not assigned a "
              "region or mapped. A record without a region has an empty value, never 'Unknown'.",
              "", "## Taxonomy and labels", "",
              "A record's label is its species, or its genus when the species is missing. A "
              "haplotype's `majority_label` is its most common species-level member label, and "
              "`majority_fraction` that label's share of the counted labels: members with a "
              "species-level label plus genus-only or open-nomenclature labels of another genus. "
              "Members without a species label (no label, or a genus-only or open label of the "
              "majority's genus) are left out of the denominator; members without any label are "
              "counted in `n_missing_label`.",
              "",
              "`assigned_sp` is `majority_label` when `majority_fraction` is more than the majority "
              "threshold (0.70 by default; run_provenance.json); otherwise, or when the top labels "
              "are tied, the most common genus (`assignment_level` genus).",
              "",
              "Labels are tidied as they are read (whitespace, underscores, the case of genus and "
              "epithet); names are not changed. label_variants.csv lists pairs of labels that may "
              "be one name written differently; nothing is merged unless a label map is given "
              "(`--label-map FILE`: from, to, reason, optional source). The map replaces deposited "
              "species labels before any label is counted, so it affects majority labels, label "
              "conflicts, evidence and queries. The species and genus columns of records.csv are "
              "the values the run used; `label_mapped_from` keeps the deposited label. "
              "run_provenance.json records the map's md5 and any rows that matched no record.",
              "",
              "With `--taxon-curation FILE`, a row of that file replaces the label rule's taxon for "
              "the haplotype it names (`assignment_source` curated); `assignment_notes` keeps the "
              "reason and the label rule's note. A row names a haplotype by `seq_md5`, or by a "
              "sequence equal to, containing or contained in its core sequence; rows that match no "
              "haplotype or several, or that disagree for one haplotype, are not applied. "
              "run_provenance.json records the file's md5 and every row not applied, with why. "
              "A curated taxon takes precedence over an evidence taxon.",
              "",
              "The evidence columns of haplotypes.csv weigh a haplotype's taxon against three "
              "lines that do not use its own majority: distance (its nearest species-level "
              "haplotype lies within that taxon's spread, `within_taxon_spread`), BIN (the "
              "members' most common BIN has one species-level label at 0.70 or more of its "
              "QC-passed records) and tree (`tree_sister_taxon`: after collapsing branches with "
              "support below 0.70, the taxon of the species-level haplotypes in the smallest "
              "clade that holds any). For a haplotype neither at species level nor curated, "
              "`suggested_sp` is a taxon at least two lines name and none contradicts. When it "
              "is the haplotype's leading species-level label and in the genus the label rule "
              "assigned, it becomes the haplotype's taxon (`assignment_source` evidence; "
              "`--no-evidence-taxa` turns this off). Other suggestions are written to "
              "taxon_curation_suggestions.csv, a curation file: delete the rows you reject and "
              "pass it to --taxon-curation. Evidence uses only the label rule's and curated "
              "taxa, never another haplotype's evidence taxon. The tree is midpoint-rooted "
              "unless an outgroup was set, so a sister clade next to the root depends on where "
              "the root falls.",
              "",
              "A label conflict is a member whose species-level label differs from the majority "
              "label, any species-level member of a haplotype whose top labels are tied, or a "
              "genus-only label of another genus. For a haplotype whose taxon is curated or set by "
              "evidence, labels are compared with that taxon instead: a species-level label conflicts when it differs "
              "from a species-level taxon, any other label when its genus differs.",
              "", "## Divergence", "",
              "Distances are uncorrected p-distances between aligned haplotype sequences. They are "
              "reported without taxonomic thresholds.",
              ""]
    if has_metadata_extras:
        lines += ["## Metadata", "",
                  "Metadata values are normalized before counting: empty and whitespace-only values "
                  "are missing; sex M/F become Male/Female and U is missing; life stage A/I/M/J "
                  "become Adult/Immature/Mature/Juvenile; reproduction S/A become Sexual/Asexual; "
                  "case variants are merged; BOLD GIS suffixes are removed.",
                  "",
                  "Counts by collection year reflect when and where specimens were sampled and "
                  "deposited in BOLD, not changes in haplotype frequency.",
                  ""]
    return lines


def write_run_readme(output_dir: Path, organism: str, region_col: str = "ocean_basin"
                     ) -> Tuple[Path, List[str]]:
    """
    Write README.md describing every file in the run folder (intermediate/
    excluded). Returns (path, undocumented), where undocumented lists files
    and 'file: column' entries without a definition.
    """
    output_dir = Path(output_dir)
    readme = output_dir / "README.md"
    rels = sorted({str(p.relative_to(output_dir)) for p in output_dir.rglob("*")
                   if p.is_file() and p.relative_to(output_dir).parts[0] != "intermediate"}
                  | {"README.md"})

    version = started = ""
    prov = output_dir / "run_provenance.json"
    if prov.exists():
        try:
            d = json.loads(prov.read_text())
            version, started = d.get("boldgenotyper_version", ""), d.get("started_at", "")
        except Exception:
            pass

    undocumented: List[str] = []
    groups: Dict[str, List[str]] = {}
    for rel in rels:
        entry = _entry_for(rel)
        if entry is None:
            undocumented.append(rel)
            continue
        groups.setdefault(entry["pattern"], []).append(rel)

    lines = [f"# {organism}: BOLDGenotyper run", ""]
    if version or started:
        lines += [f"BOLDGenotyper {version}, run started {started}.", ""]
    lines += ["This file lists every file in this folder and defines every column. "
              "Folders: `extended/` (--extended), `exports/` (--export-format), `plots/` (--plot-kit).",
              ""]
    lines += _sections(region_col, any(r.startswith("extended/metadata/") for r in rels))
    lines += ["## Files", ""]

    for entry in FILE_DOCS:
        files = groups.get(entry["pattern"])
        if not files:
            continue
        shown = ", ".join(f"`{f}`" for f in files[:12]) + (f" and {len(files) - 12} more" if len(files) > 12 else "")
        lines += [f"### {shown}", "", entry["text"], ""]
        csvs = [f for f in files if f.endswith(".csv")]
        if not csvs:
            continue
        if "matrix" in entry:
            lines += [entry["matrix"], ""]
            continue
        if "same_as" in entry:
            lines += [f"Columns as in `{entry['same_as']}`.", ""]
            continue
        columns: List[str] = []
        for f in csvs:
            for col in _csv_header(output_dir / f):
                if col not in columns:
                    columns.append(col)
        lines += ["| Column | Definition |", "|---|---|"]
        for col in columns:
            text = _column_text(col, entry, region_col)
            if text is None:
                undocumented.append(f"{entry['pattern']}: {col}")
                text = ""
            lines.append(f"| `{col}` | {text.replace(REGION, region_col)} |")
        lines.append("")

    readme.write_text("\n".join(lines) + "\n")
    if undocumented:
        logger.debug(f"README: no definition for {', '.join(undocumented)}")
    return readme, undocumented
