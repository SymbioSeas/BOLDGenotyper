"""Text-column checks work with pandas' object and str dtypes (pandas 3 reads text as str)."""

import pandas as pd

from boldgenotyper import utils
from boldgenotyper.metadata import parse_bold_tsv


def test_is_text_series():
    assert utils.is_text_series(pd.Series(["a", None], dtype=object))
    assert utils.is_text_series(pd.Series(["a", "b"]))
    assert not utils.is_text_series(pd.Series([1.0, 2.0]))
    assert not utils.is_text_series(pd.Series([True, False]))


def test_parse_strips_whitespace(tmp_path):
    p = tmp_path / "in.tsv"
    p.write_text("processid\tnuc\tspecies\n  ID1  \t  ACGT  \t  Genus alpha \n")
    df = parse_bold_tsv(p)
    assert df.loc[0, "processid"] == "ID1"
    assert df.loc[0, "nuc"] == "ACGT"
    assert df.loc[0, "species"] == "Genus alpha"
