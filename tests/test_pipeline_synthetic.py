"""End-to-end checks on a small synthetic dataset (requires MAFFT).

The planted cases are described in tests/synthetic_data.py.
"""

import shutil
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from boldgenotyper import cli, config, metadata, utils
from tests.synthetic_data import EXPECTED, make_synthetic_bold_tsv

pytestmark = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")

ORGANISM = "Synth"


def run_synthetic(workdir, **cfg_updates):
    """Run the pipeline on the synthetic TSV; return (success, output_dir)."""
    workdir = Path(workdir)
    tsv = make_synthetic_bold_tsv(workdir / "Synth.tsv")
    out = workdir / "out"
    utils.setup_logging(log_level="INFO", log_file=str(out / f"{ORGANISM}_pipeline.log"))
    cfg = config.get_default_config().update(
        output_dir=out, n_threads=1, phylogenetic__build_tree=False, **cfg_updates
    )
    ok = cli.run_pipeline(
        tsv_path=tsv, organism=ORGANISM, output_dir=out, cfg=cfg,
        no_report=True, skip_geo=True, export_plot_data=False,
        export_popgen_formats=None, metadata_summary_enabled=False,
    )
    return ok, out


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory):
    ok, out = run_synthetic(tmp_path_factory.mktemp("synthetic"))
    assert ok, "pipeline failed on synthetic data"
    return out


def test_synthetic_run_completes(synthetic_run):
    assert (synthetic_run / f"{ORGANISM}_annotated.csv").exists()


def test_input_parsed_once(tmp_path):
    """P0-12 guard: the input TSV is parsed exactly once per run."""
    with mock.patch.object(metadata, "parse_bold_tsv", wraps=metadata.parse_bold_tsv) as spy:
        ok, _ = run_synthetic(tmp_path)
    assert ok
    assert spy.call_count == 1
