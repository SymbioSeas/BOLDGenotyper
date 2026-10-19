"""Records and haplotypes: every record's status and the error filter."""

from .data import RunData
from .page import Section, bars, kv_table, num

STATUSES = [("assigned", "n_assigned"), ("assigned_corrected", "n_assigned_corrected"),
            ("excluded_error_singleton", "n_excluded_error"),
            ("excluded_core_coverage", "n_excluded_core"), ("failed_qc", "n_failed_qc"),
            ("no_sequence", "n_no_sequence")]


def build(run: RunData) -> Section:
    s = run.summary
    n_unique = int(s.get("n_unique") or 0)
    items = [(status, int(s.get(col) or 0),
              f"({100 * int(s.get(col) or 0) / n_unique:.1f}%)" if n_unique else "")
             for status, col in STATUSES]
    d = (run.config.get("haplotype") or {}).get("min_singleton_distance")
    corrected, tied, far = (int(s.get(c) or 0) for c in
                            ("n_assigned_corrected", "n_excluded_error_tied", "n_excluded_error_far"))
    html = ("<p>Every unique input record has one status (definitions in <code>README.md</code>; "
            "per record in <code>records.csv</code>).</p>\n"
            + bars(items)
            + "<h3>Error filter</h3>\n"
            + f"<p>Singleton haplotypes within p-distance {num(d, '{}')} of another haplotype are treated "
              "as likely sequencing errors; their records are assigned to the nearest retained haplotype "
              "when it is unique and within that distance.</p>\n"
            + kv_table([
                ("Assigned to the nearest haplotype (assigned_corrected)", num(corrected)),
                ("Unassigned: two or more haplotypes equally near", num(tied)),
                ("Unassigned: no haplotype within the distance", num(far)),
                ("Haplotypes retained", num(s.get("n_haplotypes"))),
                ("Rows of other markers removed", num(s.get("n_other_marker_rows") or 0)),
                ("Malformed rows dropped", num(s.get("n_malformed_rows") or 0)),
                ("Duplicate processids removed (first record kept)", num(s.get("n_duplicate_rows") or 0)),
            ]))
    return Section("records", "Records and haplotypes", html)
