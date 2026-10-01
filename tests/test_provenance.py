"""Run provenance and step warnings (P0-1, P0-6)."""

import logging

from boldgenotyper import provenance


def test_step_warnings_records_and_summarizes(caplog):
    w = provenance.StepWarnings()
    w.record("10.1 distribution map", ValueError("dash list"))
    w.record("8.3 MSA visualization", "pymsaviz not installed")
    assert len(w) == 2
    assert w.as_list()[0] == {"step": "10.1 distribution map", "message": "dash list"}
    lg = logging.getLogger("boldgenotyper.test")
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="boldgenotyper.test"):
        w.log_summary(lg, "Synth")
    lines = [r.getMessage() for r in caplog.records if r.name == "boldgenotyper.test"]
    assert lines[0] == "WARNINGS (2):"
    assert lines[-1] == "Pipeline completed with 2 warnings for Synth"


def test_step_warnings_empty_summary(caplog):
    lg = logging.getLogger("boldgenotyper.test")
    with caplog.at_level(logging.INFO, logger="boldgenotyper.test"):
        provenance.StepWarnings().log_summary(lg, "Synth")
    lines = [r.getMessage() for r in caplog.records if r.name == "boldgenotyper.test"]
    assert lines[-1] == "Pipeline completed successfully for Synth"


def test_step_warnings_captures_marked_warnings_only():
    w = provenance.StepWarnings().attach()
    vis = logging.getLogger("boldgenotyper.visualization")
    vis.warning("Failed to save distribution map: boom",
                extra={provenance.STEP_FAILURE: "10.1 distribution map"})
    vis.warning("24.7% of coordinates fall outside all regions")   # data note, unmarked
    w.detach()
    vis.warning("after detach", extra={provenance.STEP_FAILURE: "x"})
    assert w.as_list() == [{"step": "10.1 distribution map",
                            "message": "Failed to save distribution map: boom"}]


# --- RunProvenance (P0-1) ---------------------------------------------------

import hashlib
import json

from boldgenotyper import config


def test_file_md5_and_row_count(tmp_path):
    f = tmp_path / "x.tsv"
    f.write_text("processid\tnuc\na\tACGT\nb\tACGT\n")
    assert provenance.file_md5(f) == hashlib.md5(f.read_bytes()).hexdigest()
    assert provenance.count_data_rows(f) == 2


def test_tool_versions_reports_missing_tools(monkeypatch):
    monkeypatch.setattr(provenance.shutil, "which", lambda name: None)
    assert provenance.tool_versions() == {"mafft": None, "fasttree": None, "trimal": None}


def test_run_provenance_serializes_full_config(tmp_path):
    tsv = tmp_path / "in.tsv"
    tsv.write_text("processid\tnuc\na\tACGT\n")
    cfg = config.get_default_config().update(haplotype__min_singleton_distance=0.015)
    prov = provenance.RunProvenance(cfg, ["boldgenotyper", str(tsv), "--singleton-distance", "0.015"],
                                    tsv, tmp_path, download_date="2026-09-30")
    prov.applied["min_sequence_length_bp"] = 565.25
    prov.warnings.record("10.1 distribution map", "failed")
    prov.finish()
    d = json.loads(prov.write(tmp_path / "run_provenance.json").read_text())
    assert d["config"]["haplotype"]["min_singleton_distance"] == 0.015
    assert d["config"]["haplotype"]["max_singleton_distance"] == 0.05
    assert d["config"]["core_region"]["min_sequence_coverage"] == 0.95
    assert "genotype_assignment" not in d["config"]
    assert "genotype_assignment" in d["not_applied"]
    assert d["input"]["md5"] == provenance.file_md5(tsv)
    assert d["input"]["n_rows"] == 1
    assert d["input"]["download_date"] == "2026-09-30"
    assert d["command_line"].endswith("--singleton-distance 0.015")
    assert d["applied"]["min_sequence_length_bp"] == 565.25
    assert d["seeds"]["fasttree"] == provenance.FASTTREE_SEED
    assert d["warnings"] == [{"step": "10.1 distribution map", "message": "failed"}]
    assert d["status"] == "completed_with_warnings"
    assert d["boldgenotyper_version"]
