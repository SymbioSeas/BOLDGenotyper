"""One assignment status per unique input record.

Statuses, in order of precedence:

    no_sequence               no sequence in the input
    failed_qc                 failed pre-alignment QC (the reason is kept)
    excluded_core_coverage    ungapped core shorter than the per-sequence cutoff
    assigned                  core identical to a retained haplotype
    assigned_corrected        error-filtered singleton assigned to its nearest
                              retained haplotype, which is within the
                              error-filter threshold
    excluded_error_singleton  error-filtered singleton left unassigned: equidistant
                              from two or more retained haplotypes, no retained
                              haplotype within the threshold, or correction disabled
"""

from typing import Sequence, Set

import pandas as pd

STATUSES = (
    "assigned", "assigned_corrected", "excluded_error_singleton",
    "excluded_core_coverage", "failed_qc", "no_sequence",
)

COLUMNS = ["processid", "status", "haplotype_id", "nearest_haplotype",
           "distance_to_nearest", "status_note"]


def build_record_status(unique_ids: Sequence[str], sequence_ids: Set[str], qc: pd.DataFrame,
                        exact_mapping: pd.DataFrame, excluded_core: Sequence[str],
                        error_singletons: pd.DataFrame, singleton_distance: float) -> pd.DataFrame:
    """Return one row per unique record, in input order, with its status and a true note.

    Parameters
    ----------
    unique_ids : processids of the unique input records (after duplicate removal)
    sequence_ids : processids that have a sequence
    qc : processid, qc_pass, qc_fail_reason for records that went through QC
    exact_mapping : processid, haplotype_id for records identical to a retained haplotype
    excluded_core : processids excluded for insufficient core coverage
    error_singletons : DiscoveryDetails.error_singletons (processid, status, reason,
        haplotype_id, nearest_haplotype, distance_to_nearest, tied_haplotypes)
    singleton_distance : the error-filter threshold, quoted in notes
    """
    qc_idx = qc.set_index("processid")
    exact = dict(zip(exact_mapping["processid"], exact_mapping["haplotype_id"]))
    excluded = set(excluded_core)
    err = error_singletons.set_index("processid") if len(error_singletons) else None
    rows = []
    for pid in unique_ids:
        row = {"processid": pid, "status": None, "haplotype_id": None,
               "nearest_haplotype": None, "distance_to_nearest": float("nan"), "status_note": ""}
        if pid not in sequence_ids:
            row.update(status="no_sequence", status_note="No sequence in the input")
        elif pid in qc_idx.index and not bool(qc_idx.loc[pid, "qc_pass"]):
            reason = qc_idx.loc[pid, "qc_fail_reason"]
            reason = reason if isinstance(reason, str) and reason else "failed QC"
            row.update(status="failed_qc", status_note=f"Failed QC: {reason}")
        elif pid in excluded:
            row.update(status="excluded_core_coverage",
                       status_note="Ungapped core shorter than the per-sequence coverage cutoff")
        elif pid in exact:
            row.update(status="assigned", haplotype_id=exact[pid],
                       status_note="Core sequence identical to a retained haplotype")
        elif err is not None and pid in err.index:
            e = err.loc[pid]
            d = float(e["distance_to_nearest"])
            row["distance_to_nearest"] = d
            tied = e["tied_haplotypes"] if isinstance(e["tied_haplotypes"], str) else ""
            reason = e["reason"] if "reason" in e.index and isinstance(e["reason"], str) else (
                "tied" if tied else "correction_disabled")
            if e["status"] == "assigned_corrected":
                row.update(status="assigned_corrected", haplotype_id=e["haplotype_id"],
                           nearest_haplotype=e["nearest_haplotype"],
                           status_note=(f"Singleton within {singleton_distance} of another haplotype; "
                                        f"assigned to nearest retained haplotype {e['haplotype_id']} "
                                        f"(distance {d:.4f})"))
            elif reason == "tied":
                row.update(status="excluded_error_singleton",
                           status_note=(f"Error singleton equidistant ({d:.4f}) from "
                                        f"{', '.join(tied.split(';'))}; left unassigned"))
            elif reason == "beyond_threshold":
                row.update(status="excluded_error_singleton", nearest_haplotype=e["nearest_haplotype"],
                           status_note=(f"Error singleton; nearest retained haplotype "
                                        f"{e['nearest_haplotype']} is {d:.4f} away "
                                        f"(> {singleton_distance}); left unassigned"))
            else:
                row.update(status="excluded_error_singleton", nearest_haplotype=e["nearest_haplotype"],
                           status_note=(f"Error singleton; correction disabled "
                                        f"(nearest {e['nearest_haplotype']}, distance {d:.4f})"))
        else:
            raise ValueError(f"No status could be determined for record {pid}")
        rows.append(row)
    return pd.DataFrame(rows, columns=COLUMNS)
