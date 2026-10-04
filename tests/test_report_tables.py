"""Report tables show ten rows and say how many more are in which file (Phase 8, 12)."""

import pandas as pd

from boldgenotyper.report import page


def test_truncation_notice_counts_all_rows():
    html = page.table(pd.DataFrame({"a": range(25)}), ["a"], file="records.csv")
    assert "15 more in <code>records.csv</code> (25 in all)" in html


def test_no_notice_when_not_truncated():
    assert "more" not in page.table(pd.DataFrame({"a": range(5)}), ["a"])
