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


def run_synthetic(workdir, skip_geo=True, **cfg_updates):
    """Run the pipeline on the synthetic TSV; return (success, output_dir)."""
    workdir = Path(workdir)
    tsv = make_synthetic_bold_tsv(workdir / "Synth.tsv")
    out = workdir / "out"
    utils.setup_logging(log_level="INFO", log_file=str(out / f"{ORGANISM}_pipeline.log"))
    settings = dict(output_dir=out, n_threads=1, phylogenetic__build_tree=False)
    settings.update(cfg_updates)
    cfg = config.get_default_config().update(**settings)
    ok = cli.run_pipeline(
        tsv_path=tsv, organism=ORGANISM, output_dir=out, cfg=cfg,
        no_report=True, skip_geo=skip_geo, export_plot_data=False,
        export_popgen_formats=None, metadata_summary_enabled=True, temporal_analysis=True,
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


needs_fasttree = pytest.mark.skipif(shutil.which("FastTree") is None and shutil.which("fasttree") is None,
                                    reason="FastTree not installed")


@needs_fasttree
def test_tree_includes_every_non_suspect_haplotype(tmp_path):
    """P0-7: default --tree-min-members 1; exclusions listed in excluded_tips.txt."""
    from Bio import Phylo
    ok, out = run_synthetic(tmp_path, phylogenetic__build_tree=True)
    assert ok
    tips = {t.name for t in Phylo.read(out / "phylogenetic" / f"{ORGANISM}_tree.nwk", "newick").get_terminals()}
    assert tips == {"H0001", "H0002", "H0003", "H0004"}
    excluded = (out / "phylogenetic" / f"{ORGANISM}_excluded_tips.txt").read_text().splitlines()
    assert excluded == ["haplotype_id\treason"]
    d = json.loads((out / "run_provenance.json").read_text())
    assert (d["applied"]["tree_min_members"], d["applied"]["n_tree_tips"], d["applied"]["n_tree_excluded"]) == (1, 4, 0)
    assert d["warnings"] == []
    import importlib.util
    if importlib.util.find_spec("pymsaviz") is not None:   # [msa] extra; in environment.yml
        assert list((out / "phylogenetic").glob(f"{ORGANISM}_MSA_chunk*.pdf"))


@needs_fasttree
def test_tree_min_members_excludes_and_explains(tmp_path):
    ok, out = run_synthetic(tmp_path, phylogenetic__build_tree=True, phylogenetic__min_cluster_size=2)
    assert ok
    excluded = (out / "phylogenetic" / f"{ORGANISM}_excluded_tips.txt").read_text().splitlines()
    assert excluded == ["haplotype_id\treason", "H0004\tn_members 1 < 2"]


def test_ids_are_count_free_and_consistent(synthetic_run):
    """P0-8: H IDs carry no counts, and n_members agrees across files."""
    fasta = (synthetic_run / "haplotypes" / f"{ORGANISM}_haplotypes.fasta").read_text()
    headers = [line for line in fasta.splitlines() if line.startswith(">")]
    assert headers[0].startswith(">H0001 seq_md5=")
    assert headers[0].endswith(f"n_members={EXPECTED['members']['A']}")
    assert not any("_n" in h.split()[0] for h in headers)
    stats = pd.read_csv(synthetic_run / "haplotypes" / f"{ORGANISM}_haplotype_stats.csv")
    ann = pd.read_csv(synthetic_run / f"{ORGANISM}_annotated.csv")
    per_hap = ann.dropna(subset=["haplotype_id"]).groupby("haplotype_id").size()
    pd.testing.assert_series_equal(
        stats.set_index("haplotype_id")["n_members"].sort_index(), per_hap.sort_index(),
        check_names=False, check_dtype=False)
    readme = (synthetic_run / "divergence_analysis" / "README.md").read_text()
    assert "_n{" not in readme and "h1_n" not in readme


def test_one_missing_code_and_latlon_kept(synthetic_run):
    """P0-4: no 'Unknown' region code; centroid records keep lat/lon; exclusion logged."""
    ann = pd.read_csv(synthetic_run / f"{ORGANISM}_annotated.csv")
    assert "Unknown" not in set(ann["ocean_basin"].dropna())
    cent = ann["coord_status"] == "centroid"
    assert cent.any()
    assert ann.loc[cent, ["lat", "lon"]].notna().all(axis=None)
    assert not ann.loc[ann["coord_status"] == "missing", "is_geographic_quality"].any()
    log = (synthetic_run / f"{ORGANISM}_pipeline.log").read_text()
    assert "centroid coordinates; they keep lat/lon" in log


needs_goas = pytest.mark.skipif(
    not config.get_default_config().geographic.goas_shapefile_path.exists(),
    reason="GOaS shapefile not installed")


@needs_goas
def test_goas_regions_only_for_geographic_quality(tmp_path):
    """P0-4: regions only for geographic-quality records; no region is NaN, not 'Unknown'."""
    ok, out = run_synthetic(tmp_path, skip_geo=False)
    assert ok
    ann = pd.read_csv(out / f"{ORGANISM}_annotated.csv")
    assert "Unknown" not in set(ann["ocean_basin"].dropna())
    good = ann["is_geographic_quality"].astype(bool)
    assert ann.loc[~good, "ocean_basin"].isna().all()
    assert ann.loc[ann["coord"] == "[12.5, -45.25]", "ocean_basin"].notna().all()
    d = json.loads((out / "run_provenance.json").read_text())
    assert set(d["data_notes"]["coord_status"]) <= {"ok", "missing", "centroid", "zero"}
    # P0-6: every map and plot, including species-level ones, is written
    assert d["warnings"] == []
    assert (out / "visualization" / "pdf" / f"{ORGANISM}_distribution_map.pdf").exists()
    assert (out / "visualization" / "pdf" / f"{ORGANISM}_species_distribution_map.pdf").exists()
    # P0-5: country counts come from country/ocean
    coverage = pd.read_csv(next((out / "geographic_analysis").glob("*geographic_coverage.csv")))
    assert coverage["with_country"].sum() > 0


def test_divergence_readme_has_no_threshold_interpretations(synthetic_run):
    """P1-1: describe p-distance only; no taxonomic thresholds."""
    import re
    readme = (synthetic_run / "divergence_analysis" / "README.md").read_text()
    for phrase in ("cryptic", "hybridizing", "misidentification", "Typical interspecific",
                   "barcoding gap exists", "## Interpretation"):
        assert phrase.lower() not in readme.lower(), phrase
    assert not re.search(r"[<>]\s*\d+(\.\d+)?%", readme)
    assert "uncorrected p-distance" in readme


def test_label_report_replaces_contamination_outputs(synthetic_run):
    """P1-6, P1-7, P1-8: label conflicts with evidence; no contamination or confidence."""
    qc = synthetic_run / "quality_control"
    summary = pd.read_csv(qc / "haplotype_label_summary.csv").set_index("haplotype_id")
    # H0003 (C): two 'Genus alpha' and one genus-only 'Genus' member: compatible
    assert summary.loc["H0003", "n_label_conflicts"] == 0
    assert summary.loc["H0003", "majority_fraction"] == 1.0
    assert (qc / "label_conflicts.csv").exists() and (qc / "depositor_keyword_hits.csv").exists()
    assert not (qc / "potential_misidentifications.csv").exists()
    ann = pd.read_csv(synthetic_run / f"{ORGANISM}_annotated.csv")
    assert {"label", "majority_label", "label_conflict", "majority_fraction"} <= set(ann.columns)
    assert "misidentification_confidence" not in ann.columns


@pytest.fixture(scope="module")
def extended_run(tmp_path_factory):
    ok, out = run_synthetic(tmp_path_factory.mktemp("synthetic_extended"), extended=True)
    assert ok
    return out


def test_default_run_writes_no_gap_outputs(synthetic_run):
    """P1-2: no barcoding-gap outputs by default."""
    for rel in ("divergence_analysis/within_vs_between_species.csv",
                "divergence_analysis/barcoding_gap.pdf",
                "divergence_analysis/label_based_gap.csv",
                "species_analysis/species_divergence_summary.csv"):
        assert not (synthetic_run / rel).exists(), rel


def test_extended_run_writes_label_based_gap(extended_run):
    gap = pd.read_csv(extended_run / "divergence_analysis" / "label_based_gap.csv")
    assert set(gap["label"]) == {"Genus alpha", "Genus beta", "Genus gamma"}
    assert "barcoding_gap" not in gap.columns and "label_based_gap" in gap.columns
    assert not (extended_run / "divergence_analysis" / "within_vs_between_species.csv").exists()
    d = json.loads((extended_run / "run_provenance.json").read_text())
    assert d["config"]["extended"] is True and d["warnings"] == []


def test_default_metadata_output_is_coverage_only(synthetic_run):
    """P1-3, P1-4, P1-5: cross-tabs, temporal and first-sampled outputs only with --extended."""
    files = {p.name for p in (synthetic_run / "metadata_summary").iterdir()}
    assert files == {f"{ORGANISM}_metadata_coverage.json", f"{ORGANISM}_metadata_summary.csv", "README.md"}


def test_extended_metadata_outputs(extended_run):
    md = extended_run / "metadata_summary"
    sex = pd.read_csv(md / f"{ORGANISM}_sex_tabulation.csv")
    assert set(sex["value"]) == {"Male", "Female"}
    life = pd.read_csv(md / f"{ORGANISM}_life_stage_tabulation.csv")
    assert set(life["value"]) == {"Adult", "Immature"}
    first = pd.read_csv(md / f"{ORGANISM}_haplotype_first_sampled.csv")
    assert "first_sampled" in first.columns and "first_collection_date" not in first.columns
    assert not (md / f"{ORGANISM}_haplotype_emergence.csv").exists()
    assert (md / f"{ORGANISM}_temporal_trends.csv").exists()
    readme = (md / "README.md").read_text()
    assert "sampling" in readme and "emergence" not in readme.lower()
