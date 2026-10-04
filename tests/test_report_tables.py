"""Report tables: truncation notices give the table's true length."""

import pandas as pd

from boldgenotyper import reports


def test_truncation_notice_counts_all_rows():
    html = reports._dataframe_to_html(pd.DataFrame({"a": range(25)}), max_rows=10)
    assert "Showing first 10 rows of 25 total" in html


def test_no_notice_when_not_truncated():
    html = reports._dataframe_to_html(pd.DataFrame({"a": range(5)}), max_rows=10)
    assert "Showing first" not in html
