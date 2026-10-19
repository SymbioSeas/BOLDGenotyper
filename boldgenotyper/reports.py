"""
Assignment summary (summary.csv): record counts at each step with explicit
denominators, and species counts under named definitions. The HTML report is
built by the report package (boldgenotyper.report).
"""

import logging
from pathlib import Path
from typing import Dict, Optional
import pandas as pd

logger = logging.getLogger(__name__)


INPUT_ROW_COLUMNS = ["n_other_marker_rows", "n_malformed_rows", "n_duplicate_rows"]
SUMMARY_COLUMNS = [
    "n_input", "n_other_marker_rows", "n_malformed_rows", "n_duplicate_rows", "n_unique", "n_no_sequence", "n_failed_qc", "n_qc_pass",
    "n_excluded_core", "n_core_ok", "n_assigned", "n_assigned_corrected",
    "n_excluded_error", "n_excluded_error_tied", "n_excluded_error_far", "n_haplotypes",
    "pct_assigned", "pct_assigned_exact", "pct_error_filtered",
    "n_singleton_label_absent", "n_haplotype_only_groups", "n_singletons_exempt",
    "n_species_labels_input", "n_species_labels_species_level", "n_species_labels_qc_passed",
    "n_species_haplotype_assigned", "n_species_confident",
]
SPECIES_COUNT_COLUMNS = [c for c in SUMMARY_COLUMNS if c.startswith("n_species_")]

# Record statuses (record_status.STATUSES) and their summary columns
_STATUS_COLUMNS = {
    "assigned": "n_assigned",
    "assigned_corrected": "n_assigned_corrected",
    "excluded_error_singleton": "n_excluded_error",
    "excluded_core_coverage": "n_excluded_core",
    "failed_qc": "n_failed_qc",
    "no_sequence": "n_no_sequence",
}


def _pct_error_filtered(n_assigned: int, n_corrected: int, n_excluded: int) -> float:
    total = n_assigned + n_corrected + n_excluded
    return round((n_corrected + n_excluded) / total * 100, 2) if total else 0.0


def compute_assignment_summary(status_df, n_input, n_haplotypes, species_counts,
                               n_excluded_error_tied=0, n_excluded_error_far=0,
                               input_counts: Optional[Dict[str, int]] = None):
    """Assignment summary with explicit denominators.

    Percentages are of QC-passed records (n_qc_pass = unique records minus
    no_sequence and failed_qc). Returns a dict with SUMMARY_COLUMNS keys.
    """
    counts = status_df["status"].value_counts()

    def c(status):
        return int(counts.get(status, 0))

    n_unique = int(len(status_df))
    n_qc_pass = n_unique - c("no_sequence") - c("failed_qc")
    summary = {
        "n_input": int(n_input if n_input is not None else n_unique),
        "n_unique": n_unique,
        "n_no_sequence": c("no_sequence"),
        "n_failed_qc": c("failed_qc"),
        "n_qc_pass": n_qc_pass,
        "n_excluded_core": c("excluded_core_coverage"),
        "n_core_ok": n_qc_pass - c("excluded_core_coverage"),
        "n_assigned": c("assigned"),
        "n_assigned_corrected": c("assigned_corrected"),
        "n_excluded_error": c("excluded_error_singleton"),
        "n_excluded_error_tied": int(n_excluded_error_tied),
        "n_excluded_error_far": int(n_excluded_error_far),
        "n_haplotypes": int(n_haplotypes),
        "pct_assigned": round((c("assigned") + c("assigned_corrected")) / n_qc_pass * 100, 2) if n_qc_pass else 0.0,
        "pct_assigned_exact": round(c("assigned") / n_qc_pass * 100, 2) if n_qc_pass else 0.0,
        "pct_error_filtered": _pct_error_filtered(c("assigned"), c("assigned_corrected"),
                                                  c("excluded_error_singleton")),
        "n_singleton_label_absent": 0,      # filled by update_summary (Phase 13)
        "n_haplotype_only_groups": 0,
        "n_singletons_exempt": 0,
    }
    for key in SPECIES_COUNT_COLUMNS:
        summary[key] = int((species_counts or {}).get(key, 0))
    for key in INPUT_ROW_COLUMNS:
        summary[key] = int((input_counts or {}).get(key, 0))
    return {k: summary[k] for k in SUMMARY_COLUMNS}


def generate_assignment_summary(
    annotated_csv: str,
    diagnostics_csv: str,
    output_csv: str,
    n_input: Optional[int] = None,
    n_haplotypes: Optional[int] = None,
    species_counts: Optional[Dict[str, int]] = None,
    n_excluded_error_tied: int = 0,
    n_excluded_error_far: int = 0,
    input_counts: Optional[Dict[str, int]] = None,
) -> Dict[str, object]:
    """
    Write the assignment summary: one row with explicit denominators.

    Parameters
    ----------
    annotated_csv : str
        Path to the annotated CSV (kept for compatibility; statuses come from
        the diagnostics table)
    diagnostics_csv : str
        Path to the diagnostics CSV: one row per unique input record with
        'status' and 'assigned_haplotype'
    output_csv : str
        Path for output summary CSV
    n_input : int, optional
        Data rows in the input TSV, before duplicate removal (default: the
        number of unique records)
    n_haplotypes : int, optional
        Number of retained haplotypes (default: distinct assigned haplotypes)
    species_counts : dict, optional
        taxon_labels.species_count_summary() output
    n_excluded_error_tied, n_excluded_error_far : int, optional
        Error-singleton records left unassigned because they are equidistant
        from two or more haplotypes, or because no retained haplotype is within
        the error-filter threshold
    input_counts : dict, optional
        Rows removed while reading the input (n_other_marker_rows,
        n_malformed_rows, n_duplicate_rows)

    Returns
    -------
    dict
        The summary (SUMMARY_COLUMNS keys)
    """
    out = Path(output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)

    diag = pd.read_csv(diagnostics_csv)
    if n_haplotypes is None:
        n_haplotypes = diag["assigned_haplotype"].dropna().nunique() if "assigned_haplotype" in diag else 0

    summary = compute_assignment_summary(diag, n_input, n_haplotypes, species_counts,
                                         n_excluded_error_tied, n_excluded_error_far,
                                         input_counts)
    pd.DataFrame([summary], columns=SUMMARY_COLUMNS).to_csv(out, index=False)
    logger.info(f"Assignment summary saved: {out}")
    return summary


def update_summary(output_csv, values: Dict[str, object]) -> Dict[str, object]:
    """Set values in an existing summary.csv (columns written later in the run);
    returns the updated summary, or {} when the file does not exist."""
    out = Path(output_csv)
    if not out.exists():
        return {}
    summary = pd.read_csv(out).iloc[0].to_dict()
    summary.update(values)
    pd.DataFrame([summary], columns=SUMMARY_COLUMNS).to_csv(out, index=False)
    return {k: summary[k] for k in SUMMARY_COLUMNS}
