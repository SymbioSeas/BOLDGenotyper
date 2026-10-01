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


def test_log_ends_with_run_status(synthetic_run):
    """P0-6: the last log line says whether the run had warnings."""
    last = (synthetic_run / f"{ORGANISM}_pipeline.log").read_text().strip().splitlines()[-1]
    assert last.endswith(f"Pipeline completed successfully for {ORGANISM}") or \
        "Pipeline completed with" in last


def _diagnostics(out):
    return pd.read_csv(out / "haplotype_assignments" / f"{ORGANISM}_diagnostics.csv")


def test_statuses_match_planted_cases(synthetic_run):
    """P0-2: every unique record has its specific status."""
    diag = _diagnostics(synthetic_run).set_index("processid")
    assert len(diag) == EXPECTED["n_unique"]
    for pid, status in EXPECTED["status"].items():
        assert diag.loc[pid, "status"] == status, pid


def test_processid_is_bare_everywhere(synthetic_run):
    """P0-10: one processid format in every table."""
    mapping = pd.read_csv(synthetic_run / "haplotypes" / f"{ORGANISM}_haplotype_mapping.csv")
    assert mapping["processid"].str.match(r"^SYN\d{3}$").all()


def test_corrected_members_count_in_species_composition(synthetic_run):
    """P0-2: corrected members are pooled into their haplotype."""
    ann = pd.read_csv(synthetic_run / f"{ORGANISM}_annotated.csv")
    h1 = ann.loc[ann["haplotype_id"] == "H0001"]
    assert len(h1) == EXPECTED["members"]["A"]
    assert (h1["status"] == "assigned_corrected").sum() == 1


def test_no_assign_corrected_flag(tmp_path):
    ok, out = run_synthetic(tmp_path, haplotype__assign_corrected=False)
    assert ok
    diag = _diagnostics(out).set_index("processid")
    assert diag.loc["SYN021", "status"] == "excluded_error_singleton"


def test_assignment_summary_denominators(synthetic_run):
    """P0-3, P0-9: explicit denominators and named species counts."""
    s = pd.read_csv(synthetic_run / "reports" / f"{ORGANISM}_assignment_summary.csv").iloc[0]
    assert s["n_input"] == EXPECTED["n_input"]
    assert s["n_unique"] == EXPECTED["n_unique"]
    assert s["n_qc_pass"] == 23 and s["n_core_ok"] == 22
    assert s["pct_assigned"] < 100.0
    assert s["n_assigned_corrected"] == 1 and s["n_excluded_error"] == 1
    assert s["n_haplotypes"] == 4
    assert s["n_species_labels_input"] == 3          # Genus alpha, beta, gamma
    assert s["n_species_haplotype_assigned"] == 3


import json


def test_run_provenance_written_and_complete(synthetic_run):
    """P0-1: provenance records the full config, the values applied, and tools."""
    d = json.loads((synthetic_run / "run_provenance.json").read_text())
    assert d["status"] in ("completed", "completed_with_warnings")
    assert d["config"]["haplotype"]["min_singleton_distance"] == 0.005
    assert d["applied"]["min_sequence_length_bp"] > 0
    assert d["applied"]["core_region_length_bp"] > 0
    assert d["tools"]["mafft"]
    assert d["input"]["n_rows"] == EXPECTED["n_input"]
    assert d["data_notes"]["n_duplicate_processids"] == 1
    assert not (synthetic_run / f"{ORGANISM}_pipeline_parameters.json").exists()


def test_fresh_run_has_no_previous_run_warning(synthetic_run):
    d = json.loads((synthetic_run / "run_provenance.json").read_text())
    assert not [w for w in d["warnings"] if w["step"] == "0 output folder"]


def test_reused_output_folder_warns(tmp_path):
    ok, _ = run_synthetic(tmp_path)
    assert ok
    ok, out = run_synthetic(tmp_path)
    assert ok
    d = json.loads((out / "run_provenance.json").read_text())
    assert "0 output folder" in [w["step"] for w in d["warnings"]]
    assert d["status"] == "completed_with_warnings"
    last = (out / f"{ORGANISM}_pipeline.log").read_text().strip().splitlines()[-1]
    assert "Pipeline completed with" in last
