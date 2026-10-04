"""boldgenotyper-query end to end against a synthetic v2 run folder (requires MAFFT).

Synthetic haplotypes (tests/synthetic_data.py): A and B differ only at bases 32
and 62, so a query on A[100:320] cannot tell them apart; C and D differ from A
inside that window.
"""

import hashlib
import json
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
    summary = pd.read_csv(out / "query_calls.csv", keep_default_na=False)
    return summary.set_index("query_id"), out


def test_equivalence_set_and_pooled_species(tmp_path, run_dir):
    s, out = _query(tmp_path, run_dir, {"window": WINDOW})
    row = s.loc["window"]
    assert row["call_level"] == "haplotype_set"
    assert row["equivalence_set"] == "H0001;H0002"
    assert "set_not_identical" not in row["flags"]               # identical over the query
    assert row["pooled_member_labels"] == "Genus alpha 13; Genus beta 4"
    assert row["species_call"] == "Genus alpha"
    assert round(float(row["species_support"]), 2) == 0.76
    assert row["flags"] == "label_conflict_in_set"
    assert row["call_basis"] == "pooled labels"
    assert int(row["query_length"]) == 220 and float(row["query_coverage"]) == 1.0
    from boldgenotyper import query_outputs as qo
    results = pd.read_csv(out / "query_candidates.csv")
    assert list(results.columns) == qo.CANDIDATE_COLUMNS
    assert not (out / "query_results.json").exists()            # only with --json


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


def test_run_folder_curated_taxa_are_inherited(tmp_path, run_dir):
    ref = tmp_path / "ref"
    shutil.copytree(run_dir, ref)
    haps = pd.read_csv(ref / "haplotypes.csv")
    haps.loc[haps["haplotype_id"] == "H0002", "assignment_source"] = "curated"
    haps.to_csv(ref / "haplotypes.csv", index=False)
    s, _ = _query(tmp_path, ref, {"window": WINDOW})
    row = s.loc["window"]
    assert row["species_call"] == "" and row["call_basis"] == ""
    assert "curated_taxa_disagree" in row["flags"]


def test_query_time_curation(tmp_path, run_dir):
    haps = pd.read_csv(run_dir / "haplotypes.csv").set_index("haplotype_id")
    cur = tmp_path / "cur.csv"
    cur.write_text(f"seq_md5,taxon,reason\n{haps.loc['H0002', 'seq_md5']},Genus alpha,test\n"
                   f"ffff,Genus nulla,no match\n")
    s, out = _query(tmp_path, run_dir, {"window": WINDOW}, "--taxon-curation", str(cur))
    row = s.loc["window"]
    assert row["species_call"] == "Genus alpha" and row["call_basis"] == "haplotype taxa"
    prov = json.loads((out / "query_provenance.json").read_text())
    assert prov["taxon_curation"]["md5"] == hashlib.md5(cur.read_bytes()).hexdigest()
    assert prov["taxon_curation"]["n_applied"] == 1 and prov["taxon_curation"]["n_not_applied"] == 1
    assert prov["parameters"]["taxon_curation"] == str(cur)


def test_query_time_curation_needs_a_run_folder(tmp_path, run_dir):
    cur = tmp_path / "cur.csv"
    cur.write_text("seq_md5,taxon,reason\nabc,Genus alpha,test\n")
    fasta = tmp_path / "q.fasta"
    fasta.write_text(f">window\n{WINDOW}\n")
    code = cli.main_query(["--query", str(fasta), "--haplotypes", str(run_dir / "haplotypes.fasta"),
                           "--output", str(tmp_path / "o"), "--taxon-curation", str(cur)])
    assert code == 1


def test_query_readme_documents_every_file_and_column(tmp_path, run_dir):
    s, out = _query(tmp_path, run_dir, {"window": WINDOW}, "--search-records", "--json")
    readme = (out / "README.md").read_text()
    assert "## Which file to report" in readme and "query_calls.csv" in readme
    md5 = hashlib.md5((run_dir / "run_provenance.json").read_bytes()).hexdigest()
    assert md5 in readme
    from boldgenotyper import query_outputs as qo
    assert qo.write_query_readme(out, json.loads((out / "query_provenance.json").read_text()))[1] == []


def test_bare_fasta_reference_has_no_species_support(tmp_path, run_dir):
    fasta = run_dir / "haplotypes.fasta"
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
    assert qa.load_reference_run(run_dir / "haplotypes.fasta") is None


def test_reference_block_and_haplotypes_searched(tmp_path, run_dir):
    s, out = _query(tmp_path, run_dir, {"window": WINDOW}, "--top-n", "1", "--json")
    md5 = hashlib.md5((run_dir / "run_provenance.json").read_bytes()).hexdigest()
    assert s.loc["window", "reference_run"] == md5
    prov = json.loads((out / "query_provenance.json").read_text())
    assert prov["reference_run"]["provenance_md5"] == md5
    assert prov["reference_fasta"]["n_haplotypes"] == 4 and prov["query"]["n_sequences"] == 1
    assert "Haplotypes Searched: 4" in (out / "query_alignments.txt").read_text()
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


def test_v1_style_run_folder_gives_haplotype_matching_only(tmp_path, run_dir):
    v1 = tmp_path / "v1_run"
    shutil.copytree(run_dir, v1)
    (v1 / "run_provenance.json").unlink()
    s, out = _query(tmp_path, v1, {"window": WINDOW})
    row = s.loc["window"]
    assert row["call_level"] == "haplotype_set" and row["equivalence_set"] == "H0001;H0002"
    assert row["flags"] == "no_reference_labels" and row["species_call"] == ""
    assert row["reference_run"] == ""
    assert json.loads((out / "query_provenance.json").read_text())["reference_run"] is None


# --- v2 run-folder layout (P2) ------------------------------------------------

def test_resolve_reference_v2_and_v1_layouts(tmp_path, run_dir):
    assert (run_dir / "haplotypes.fasta").exists()
    fasta, analysis = qa.resolve_reference(str(run_dir), None)
    assert fasta == run_dir / "haplotypes.fasta" and analysis == run_dir
    v1 = tmp_path / "v1"
    (v1 / "haplotypes").mkdir(parents=True)
    shutil.copy(run_dir / "haplotypes.fasta", v1 / "haplotypes" / "Old_haplotypes.fasta")
    assert qa.resolve_reference(str(v1), None)[0] == v1 / "haplotypes" / "Old_haplotypes.fasta"


def test_metadata_from_v2_haplotypes_table(run_dir):
    from boldgenotyper import haplotype_query as hq
    meta = hq.load_haplotype_metadata(run_dir).set_index("haplotype_id")
    assert meta.loc["H0001", "haplotype_sp"] == "Genus alpha H0001"
    assert meta.loc["H0001", "n_members"] == 13
