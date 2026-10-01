"""boldgenotyper-query end to end against a synthetic v2 run folder (requires MAFFT).

Synthetic haplotypes (tests/synthetic_data.py): A and B differ only at bases 32
and 62, so a query on A[100:320] cannot tell them apart; C and D differ from A
inside that window.
"""

import logging
import shutil

import pandas as pd
import pytest

from boldgenotyper import cli
from tests.synthetic_data import A
from tests.test_pipeline_synthetic import run_synthetic

pytestmark = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")

WINDOW = A[100:320]


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    ok, out = run_synthetic(tmp_path_factory.mktemp("query_ref"))
    assert ok
    return out


def _query(tmp_path, run_dir, records, *extra):
    fasta = tmp_path / "queries.fasta"
    fasta.write_text("".join(f">{k}\n{v}\n" for k, v in records.items()))
    out = tmp_path / "query_out"
    code = cli.main_query(["--query", str(fasta), "--haplotypes", str(run_dir),
                           "--output", str(out), *extra])
    assert code == 0
    summary = pd.read_csv(out / "query_summary.csv", comment="#", keep_default_na=False)
    return summary.set_index("query_id"), out


def test_equivalence_set_and_pooled_species(tmp_path, run_dir):
    s, out = _query(tmp_path, run_dir, {"window": WINDOW})
    row = s.loc["window"]
    assert row["call_level"] == "haplotype_set"
    assert row["equivalence_set"] == "H0001;H0002"
    assert row["set_identical_over_region"] in (True, "True")
    assert row["pooled_member_labels"] == "Genus alpha 13; Genus beta 4"
    assert row["species_call"] == "Genus alpha"
    assert round(float(row["species_support"]), 2) == 0.76
    assert row["flags"] == "label_conflict_in_set"
    assert int(row["query_length"]) == 220 and float(row["query_coverage"]) == 1.0
    results = pd.read_csv(out / "query_results.csv")
    assert "match_quality" not in results.columns


def test_gapped_input_is_cleaned(tmp_path, run_dir, caplog):
    with caplog.at_level(logging.WARNING):
        s, _ = _query(tmp_path, run_dir, {"aligned": "---" + WINDOW[:50] + "-" + WINDOW[50:] + ".."})
    assert int(s.loc["aligned", "query_length"]) == 220
    assert s.loc["aligned", "equivalence_set"] == "H0001;H0002"
    assert any("looks aligned" in r.getMessage() for r in caplog.records)


def test_indel_and_no_match_flags(tmp_path, run_dir):
    s, _ = _query(tmp_path, run_dir, {
        "deletion": WINDOW[:110] + WINDOW[111:],
        "unrelated": "ACGTTGCA" * 30,
    })
    assert "indel_in_query" in s.loc["deletion", "flags"]
    assert int(s.loc["deletion", "n_gaps"]) == 1 and int(s.loc["deletion", "n_mismatches"]) == 0
    assert s.loc["unrelated", "call_level"] == "none"
    assert s.loc["unrelated", "species_call"] == ""


def test_deprecated_tie_margin_still_runs(tmp_path, run_dir, capsys):
    s, _ = _query(tmp_path, run_dir, {"window": WINDOW}, "--tie-margin", "0.005")
    assert "deprecated" in capsys.readouterr().err
    assert s.loc["window", "call_level"] == "haplotype_set"


def test_bare_fasta_reference_has_no_species_support(tmp_path, run_dir):
    fasta = next((run_dir / "haplotypes").glob("*_haplotypes.fasta"))
    s, _ = _query(tmp_path, fasta, {"window": WINDOW})
    assert s.loc["window", "call_level"] == "haplotype_set"
    assert s.loc["window", "flags"] == "no_reference_labels"
    assert s.loc["window", "species_call"] == ""


# --- Reference reporting and record search (Q8, Q9, Q10) -------------------

import hashlib
import json

from boldgenotyper import query_assignment as qa


def test_load_reference_run(run_dir):
    ref = qa.load_reference_run(run_dir)
    md5 = hashlib.md5((run_dir / "run_provenance.json").read_bytes()).hexdigest()
    assert ref["provenance_md5"] == md5
    assert ref["n_haplotypes"] == 4
    assert ref["boldgenotyper_version"] and ref["run_date"]
    assert qa.load_reference_run(next((run_dir / "haplotypes").glob("*.fasta"))) is None


def test_reference_block_and_haplotypes_searched(tmp_path, run_dir):
    s, out = _query(tmp_path, run_dir, {"window": WINDOW}, "--top-n", "1")
    md5 = hashlib.md5((run_dir / "run_provenance.json").read_bytes()).hexdigest()
    assert s.loc["window", "reference_run_id"] == md5
    header = [l for l in (out / "query_summary.csv").read_text().splitlines() if l.startswith("#")]
    assert f"# reference_run.provenance_md5: {md5}" in header
    assert "Haplotypes Searched: 4" in (out / "query_results_detailed.txt").read_text()
    j = json.loads((out / "query_results.json").read_text())
    assert j["reference_run"]["provenance_md5"] == md5
    assert j["query_metadata"]["n_reference_haplotypes"] == 4
    assert j["parameters"]["max_mismatch_delta"] == 0
    assert j["queries"][0]["verdict"]["equivalence_set"] == ["H0001", "H0002"]


def test_pipeline_writes_qc_passed_records(run_dir):
    text = (run_dir / "records_qc_passed.fasta").read_text()
    ids = [line[1:] for line in text.splitlines() if line.startswith(">")]
    assert len(ids) == 23 and "SYN024" not in ids and "SYN022" in ids


def test_search_records_finds_unassigned_records(tmp_path, run_dir):
    s, out = _query(tmp_path, run_dir, {"window": WINDOW}, "--search-records")
    assert int(s.loc["window", "n_exact_record_matches"]) == 19
    m = pd.read_csv(out / "query_record_matches.csv").set_index("processid")
    assert m.loc["SYN022", "status"] == "excluded_error_singleton"
    assert pd.isna(m.loc["SYN022", "haplotype_id"])
    assert m.loc["SYN021", "status"] == "assigned_corrected"
    assert set(m["match_type"]) == {"query_within_record"}
