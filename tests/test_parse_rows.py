"""parse_bold_tsv: every row read or reported (A2); the target marker chosen before duplicates (A1)."""

import pytest

from boldgenotyper.metadata import MalformedTSVError, parse_bold_tsv

SEQ = "ACGT" * 60
HEADER = ["processid", "marker_code", "species", "notes", "nuc"]


def _tsv(tmp_path, header, rows, name="in.tsv"):
    path = tmp_path / name
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_other_marker_row_first_keeps_coi(tmp_path):
    path = _tsv(tmp_path, HEADER, [
        ["P1", "12S", "Genus alpha", "", "TTTT" * 60],
        ["P1", "COI-5P", "Genus alpha", "", SEQ],
    ])
    df = parse_bold_tsv(path)
    assert df["processid"].tolist() == ["P1"]
    assert df.loc[df.index[0], "marker_code"] == "COI-5P"
    assert df.attrs["n_other_marker_rows"] == 1
    assert df.attrs["n_duplicate_rows"] == 0


def test_specimen_with_only_other_markers_is_counted(tmp_path):
    path = _tsv(tmp_path, HEADER, [
        ["P1", "COI-5P", "Genus alpha", "", SEQ],
        ["P2", "COII", "Genus alpha", "", SEQ],
    ])
    df = parse_bold_tsv(path)
    assert "P2" not in set(df["processid"])
    assert df.attrs["n_other_marker_rows"] == 1


def test_no_marker_column_keeps_all_rows(tmp_path):
    header = ["processid", "species", "nuc"]
    path = _tsv(tmp_path, header, [["P1", "Genus alpha", SEQ], ["P2", "Genus beta", SEQ]])
    df = parse_bold_tsv(path)
    assert df["processid"].tolist() == ["P1", "P2"]
    assert df.attrs["n_other_marker_rows"] == 0


def test_true_duplicates_counted_after_marker_selection(tmp_path):
    path = _tsv(tmp_path, HEADER, [
        ["P3", "COI-5P", "Genus alpha", "", SEQ],
        ["P3", "COI-5P", "Genus alpha", "", SEQ],
    ])
    df = parse_bold_tsv(path)
    assert df["processid"].tolist() == ["P3"]
    assert df.attrs["n_duplicate_rows"] == 1


def _malformed(tmp_path):
    return _tsv(tmp_path, HEADER, [
        ["P1", "COI-5P", "Genus alpha", "", SEQ],
        ["P2", "COI-5P", "Genus alpha", "", SEQ],
        ["P3", "COI-5P", "Genus alpha", "", SEQ, "x", "y", "z"],   # file line 4
        ["P4", "COI-5P", "Genus alpha", "", SEQ],
    ])


def test_malformed_rows_stop_by_default(tmp_path):
    with pytest.raises(MalformedTSVError) as e:
        parse_bold_tsv(_malformed(tmp_path))
    assert e.value.line_numbers == [4]
    assert e.value.n_malformed == 1 and e.value.expected_fields == len(HEADER)
    assert "--allow-malformed-rows" in str(e.value)


def test_allow_malformed_rows_drops_and_counts(tmp_path):
    df = parse_bold_tsv(_malformed(tmp_path), allow_malformed_rows=True)
    assert df["processid"].tolist() == ["P1", "P2", "P4"]
    assert df.attrs["n_malformed_rows"] == 1
    assert df.attrs["malformed_lines"] == [4]
    assert df.attrs["n_rows"] == 4


def test_leading_quote_does_not_swallow_lines(tmp_path):
    path = _tsv(tmp_path, HEADER, [
        ["P1", "COI-5P", "Genus alpha", '"quoted start', SEQ],
        ["P2", "COI-5P", "Genus alpha", "", SEQ],
        ["P3", "COI-5P", "Genus alpha", "", SEQ],
    ])
    df = parse_bold_tsv(path)
    assert df["processid"].tolist() == ["P1", "P2", "P3"]
    assert df.loc[df["processid"] == "P1", "notes"].iloc[0] == '"quoted start'
