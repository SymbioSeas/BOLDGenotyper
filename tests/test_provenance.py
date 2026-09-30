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
