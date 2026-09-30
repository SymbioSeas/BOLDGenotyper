"""P0-12: each pipeline run starts a fresh log (v1 appended to the previous log)."""

import logging

from boldgenotyper import utils


def _close_package_handlers():
    lg = logging.getLogger("boldgenotyper")
    for h in list(lg.handlers):
        h.close()
        lg.removeHandler(h)


def test_setup_logging_starts_a_fresh_log(tmp_path):
    log = tmp_path / "run.log"
    utils.setup_logging(log_level="INFO", log_file=str(log))
    logging.getLogger("boldgenotyper").info("first run marker")
    utils.setup_logging(log_level="INFO", log_file=str(log))
    logging.getLogger("boldgenotyper").info("second run marker")
    _close_package_handlers()
    text = log.read_text()
    assert "first run marker" not in text
    assert "second run marker" in text
    assert text.count("Logging to file") == 1


def test_setup_logging_append_mode_is_explicit(tmp_path):
    log = tmp_path / "run.log"
    utils.setup_logging(log_file=str(log))
    logging.getLogger("boldgenotyper").info("kept")
    utils.setup_logging(log_file=str(log), mode="a")
    _close_package_handlers()
    assert "kept" in log.read_text()


def test_setup_logging_closes_replaced_handlers(tmp_path):
    utils.setup_logging(log_file=str(tmp_path / "a.log"))
    first = [h for h in logging.getLogger("boldgenotyper").handlers
             if isinstance(h, logging.FileHandler)][0]
    utils.setup_logging(log_file=str(tmp_path / "b.log"))
    _close_package_handlers()
    assert first.stream is None  # closed
