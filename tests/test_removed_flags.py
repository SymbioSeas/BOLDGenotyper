"""Flags removed in 2.0.0 fail with one line naming the replacement."""

import sys

import pytest

from boldgenotyper import cli


@pytest.mark.parametrize("flag", ["--similarity-threshold", "--tie-margin", "--tie-threshold",
                                  "--no-export-plot-data", "--no-export-format",
                                  "--no-metadata-analysis"])
def test_main_removed_flags(flag):
    msg = cli.removed_flag_error(["in.tsv", flag], "boldgenotyper")
    assert msg.startswith(f"{flag} was removed in 2.0.0")
    assert "MIGRATION_v1_to_v2.md" in msg


def test_flag_with_value_syntax():
    assert cli.removed_flag_error(["in.tsv", "--tie-margin=0.01"], "boldgenotyper")


def test_query_tie_margin_names_replacement():
    msg = cli.removed_flag_error(["--tie-margin", "0.005"], "boldgenotyper-query")
    assert "--max-mismatch-delta" in msg


def test_current_flags_pass():
    assert cli.removed_flag_error(["in.tsv", "--singleton-distance", "0.01"], "boldgenotyper") is None
    assert cli.removed_flag_error(["--max-mismatch-delta", "0"], "boldgenotyper-query") is None


def test_main_exits_2(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["boldgenotyper", "in.tsv", "--tie-threshold", "0.9"])
    with pytest.raises(SystemExit) as e:
        cli.main()
    assert e.value.code == 2
    assert "--tie-threshold was removed in 2.0.0" in capsys.readouterr().err


def test_query_returns_2(capsys):
    assert cli.main_query(["--query", "q.fasta", "--haplotypes", "h.fasta", "--tie-margin", "0.005"]) == 2
    assert "--tie-margin was removed in 2.0.0" in capsys.readouterr().err
