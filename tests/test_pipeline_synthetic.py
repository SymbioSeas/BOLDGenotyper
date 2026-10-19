"""End-to-end checks on a small synthetic dataset (requires MAFFT).

The planted cases are described in tests/synthetic_data.py.
"""

import re
import shutil
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from boldgenotyper import cli, config, metadata, utils
from tests.synthetic_data import EXPECTED, make_synthetic_bold_tsv

pytestmark = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")

ORGANISM = "Synth"


def run_synthetic(workdir, skip_geo=True, metadata_summary_enabled=True, plot_kit=False,
                  export_formats=None, report=False, duplicate_at=None, shapefile_path=None,
                  shapefile_field=None, geo_category="ocean_basin", taxon_curation=None,
                  split_labels=False, edit_tsv=None, marker="COI-5P", allow_malformed_rows=False,
                  label_map=None, **cfg_updates):
    """Run the pipeline on the synthetic TSV; return (success, output_dir).

    edit_tsv(path), when given, edits the synthetic TSV before the run."""
    workdir = Path(workdir)
    tsv = make_synthetic_bold_tsv(workdir / "Synth.tsv", duplicate_at=duplicate_at,
                                  split_labels=split_labels)
    if edit_tsv is not None:
        edit_tsv(tsv)
    out = workdir / "out"
    utils.setup_logging(log_level="INFO", log_file=str(out / f"{ORGANISM}_pipeline.log"))
    settings = dict(output_dir=out, n_threads=1, phylogenetic__build_tree=False)
    settings.update(cfg_updates)
    cfg = config.get_default_config().update(**settings)
    ok = cli.run_pipeline(
        tsv_path=tsv, organism=ORGANISM, output_dir=out, cfg=cfg,
        no_report=not report, skip_geo=skip_geo, export_plot_data=plot_kit,
        export_popgen_formats=export_formats, metadata_summary_enabled=metadata_summary_enabled,
        temporal_analysis=True, shapefile_path=shapefile_path, shapefile_field=shapefile_field,
        geo_category=geo_category, taxon_curation=taxon_curation,
        marker=marker, allow_malformed_rows=allow_malformed_rows, label_map=label_map,
    )
    return ok, out


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory):
    ok, out = run_synthetic(tmp_path_factory.mktemp("synthetic"))
    assert ok, "pipeline failed on synthetic data"
    return out


DEFAULT_FILES = {
    f"{ORGANISM}_pipeline.log", "run_provenance.json", "summary.csv", "records.csv",
    "records_qc_passed.fasta", "haplotypes.fasta", "haplotypes.csv", "label_conflicts.csv",
    "label_variants.csv", "sources.csv", "metadata_coverage.csv", "divergence/haplotypes_aligned.fasta",
    "divergence/pairwise_p_distance.csv", "README.md",
}


def _files(out):
    return {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()}


def test_synthetic_run_completes(synthetic_run):
    """P2: the default run folder holds exactly the v2 default files."""
    assert _files(synthetic_run) == DEFAULT_FILES


def _records(out):
    return pd.read_csv(out / "records.csv")


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


def test_statuses_match_planted_cases(synthetic_run):
    """P0-2: every unique record has its specific status."""
    diag = _records(synthetic_run).set_index("processid")
    assert len(diag) == EXPECTED["n_unique"]
    for pid, status in EXPECTED["status"].items():
        assert diag.loc[pid, "status"] == status, pid


def test_processid_is_bare_everywhere(synthetic_run):
    """P0-10: one processid format in every table."""
    assert _records(synthetic_run)["processid"].str.match(r"^SYN\d{3}$").all()
    fasta = (synthetic_run / "records_qc_passed.fasta").read_text()
    assert all(re.match(r"^>SYN\d{3}$", l) for l in fasta.splitlines() if l.startswith(">"))


def test_corrected_members_count_in_species_composition(synthetic_run):
    """P0-2: corrected members are pooled into their haplotype."""
    ann = _records(synthetic_run)
    h1 = ann.loc[ann["haplotype_id"] == "H0001"]
    assert len(h1) == EXPECTED["members"]["A"]
    assert (h1["status"] == "assigned_corrected").sum() == 1


def test_no_assign_corrected_flag(tmp_path):
    ok, out = run_synthetic(tmp_path, haplotype__assign_corrected=False)
    assert ok
    diag = _records(out).set_index("processid")
    assert diag.loc["SYN021", "status"] == "excluded_error_singleton"


def test_assignment_summary_denominators(synthetic_run):
    """P0-3, P0-9: explicit denominators and named species counts."""
    s = pd.read_csv(synthetic_run / "summary.csv").iloc[0]
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
    assert d["applied"]["input_rows"]["n_duplicate_rows"] == 1
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
    tips = {t.name for t in Phylo.read(out / "phylogeny" / "tree.nwk", "newick").get_terminals()}
    assert tips == {"H0001", "H0002", "H0003", "H0004"}
    excluded = (out / "phylogeny" / "excluded_tips.txt").read_text().splitlines()
    assert excluded == ["haplotype_id\treason"]
    d = json.loads((out / "run_provenance.json").read_text())
    assert (d["applied"]["tree_min_members"], d["applied"]["n_tree_tips"], d["applied"]["n_tree_excluded"]) == (1, 4, 0)
    assert d["warnings"] == []
    assert {"phylogeny/tree.nwk", "phylogeny/tree_labeled.nwk", "phylogeny/excluded_tips.txt",
            "phylogeny/tree.pdf"} <= _files(out)
    assert pd.read_csv(out / "haplotypes.csv")["in_tree"].all()
    assert not (out / "extended").exists()                     # no MSA plots by default


@needs_fasttree
def test_tree_min_members_excludes_and_explains(tmp_path):
    ok, out = run_synthetic(tmp_path, phylogenetic__build_tree=True, phylogenetic__min_cluster_size=2)
    assert ok
    excluded = (out / "phylogeny" / "excluded_tips.txt").read_text().splitlines()
    assert excluded == ["haplotype_id\treason", "H0004\tn_members 1 < 2"]
    # three tips are too few for FastTree: no tree, so in_tree is empty
    assert not (out / "phylogeny" / "tree.nwk").exists()
    assert pd.read_csv(out / "haplotypes.csv")["in_tree"].isna().all()


@needs_fasttree
def test_tree_is_midpoint_rooted_by_default(tmp_path):
    from Bio import Phylo
    ok, out = run_synthetic(tmp_path / "mid", phylogenetic__build_tree=True, report=True)
    assert ok
    for name in ("tree.nwk", "tree_labeled.nwk"):
        assert len(Phylo.read(out / "phylogeny" / name, "newick").root.clades) == 2
    d = json.loads((out / "run_provenance.json").read_text())
    assert d["applied"]["tree_rooting"] == {"method": "midpoint"} and d["warnings"] == []
    assert "Midpoint" in (out / "report.html").read_text()

    ok, out = run_synthetic(tmp_path / "unrooted", phylogenetic__build_tree=True,
                            phylogenetic__midpoint_root=False)
    d = json.loads((out / "run_provenance.json").read_text())
    assert d["applied"]["tree_rooting"] == {"method": "none"}
    assert len(Phylo.read(out / "phylogeny" / "tree.nwk", "newick").root.clades) == 3


@needs_fasttree
def test_tree_is_rooted_on_an_outgroup_fasta(tmp_path):
    from Bio import Phylo
    from tests.synthetic_data import A, mutate
    outgroup = tmp_path / "outgroup.fasta"
    outgroup.write_text(">OUTGROUP1\n" + mutate(A, range(100, 160)) + "\n")
    ok, out = run_synthetic(tmp_path / "run", phylogenetic__build_tree=True,
                            phylogenetic__outgroup_fasta=str(outgroup))
    assert ok
    tree = Phylo.read(out / "phylogeny" / "tree.nwk", "newick")
    assert ["OUTGROUP1"] in [[t.name for t in c.get_terminals()] for c in tree.root.clades]
    d = json.loads((out / "run_provenance.json").read_text())
    assert d["applied"]["tree_rooting"] == {"method": "outgroup_fasta", "outgroup": str(outgroup)}


@needs_fasttree
def test_msa_plots_only_with_extended_and_pdf_only(tmp_path):
    """P2, P3: MSA chunk plots go to extended/divergence/msa/, one format."""
    import importlib.util
    if importlib.util.find_spec("pymsaviz") is None:      # [msa] extra; in environment.yml
        pytest.skip("pymsaviz not installed")
    ok, out = run_synthetic(tmp_path, phylogenetic__build_tree=True, extended=True)
    assert ok
    msa = {f for f in _files(out) if "MSA_chunk" in f}
    assert msa and all(f.startswith("extended/divergence/msa/") and f.endswith(".pdf") for f in msa)


def test_ids_are_count_free_and_consistent(synthetic_run):
    """P0-8: H IDs carry no counts, and n_members agrees across files."""
    fasta = (synthetic_run / "haplotypes.fasta").read_text()
    headers = [line for line in fasta.splitlines() if line.startswith(">")]
    assert headers[0].startswith(">H0001 seq_md5=")
    assert headers[0].endswith(f"n_members={EXPECTED['members']['A']}")
    assert not any("_n" in h.split()[0] for h in headers)
    stats = pd.read_csv(synthetic_run / "haplotypes.csv")
    ann = _records(synthetic_run)
    per_hap = ann.dropna(subset=["haplotype_id"]).groupby("haplotype_id").size()
    pd.testing.assert_series_equal(
        stats.set_index("haplotype_id")["n_members"].sort_index(), per_hap.sort_index(),
        check_names=False, check_dtype=False)
    readme = (synthetic_run / "README.md").read_text()
    assert "_n{" not in readme and "h1_n" not in readme


def test_one_missing_code_and_latlon_kept(synthetic_run):
    """P0-4: no 'Unknown' region code; centroid records keep lat/lon; exclusion logged."""
    ann = _records(synthetic_run)
    assert "Unknown" not in set(ann["ocean_basin"].dropna())
    cent = ann["coord_status"] == "centroid"
    assert cent.any()
    assert ann.loc[cent, ["lat", "lon"]].notna().all(axis=None)
    assert ann.loc[ann["coord_status"] != "ok", "ocean_basin"].isna().all()
    log = (synthetic_run / f"{ORGANISM}_pipeline.log").read_text()
    assert "centroid coordinates; they keep lat/lon" in log


needs_goas = pytest.mark.skipif(
    not config.get_default_config().geographic.goas_shapefile_path.exists(),
    reason="GOaS shapefile not installed")


@needs_goas
def test_goas_regions_only_for_geographic_quality(tmp_path):
    """P0-4: regions only for geographic-quality records; no region is NaN, not 'Unknown'."""
    ok, out = run_synthetic(tmp_path, skip_geo=False, extended=True)
    assert ok
    ann = _records(out)
    assert "Unknown" not in set(ann["ocean_basin"].dropna())
    good = ann["coord_status"] == "ok"
    assert ann.loc[~good, "ocean_basin"].isna().all()
    assert ann.loc[ann["coord"] == "[12.5, -45.25]", "ocean_basin"].notna().all()
    d = json.loads((out / "run_provenance.json").read_text())
    assert set(d["data_notes"]["coord_status"]) <= {"ok", "missing", "centroid", "zero"}
    # P0-6: every map and plot, including species-level ones, is written
    assert d["warnings"] == []
    assert {"geography/distribution_map.pdf", "geography/distribution_map.svg",
            "geography/haplotype_by_region.csv",
            "extended/geography/region_abundance_total.pdf"} <= _files(out)
    # P0-5: country counts come from country/ocean
    assert pd.read_csv(out / "haplotypes.csv")["n_with_country"].sum() > 0


def test_regions_follow_each_record_after_an_early_duplicate(tmp_path):
    """
    P0-4: lat/lon and the region belong to each record's own coordinate. A duplicate
    processid near the top leaves a gap in the row index; v1.0.0 then assigned every
    later record the coordinates and region of another record.
    """
    import geopandas as gpd
    from shapely.geometry import box

    shp = tmp_path / "halves.shp"
    gpd.GeoDataFrame({"half": ["west", "east"]},
                     geometry=[box(-180, -90, 0, 90), box(0, -90, 180, 90)],
                     crs="EPSG:4326").to_file(shp)
    ok, out = run_synthetic(tmp_path, skip_geo=False, duplicate_at=1, shapefile_path=shp,
                            shapefile_field="half", geo_category="hemisphere")
    assert ok
    rec = _records(out)
    good = rec[rec["coord_status"] == "ok"]
    assert len(good) >= 2 and set(good["hemisphere"]) == {"west", "east"}
    own = good["coord"].str.strip("[]").str.split(",", expand=True).astype(float)
    assert (own[0] == good["lat"]).all() and (own[1] == good["lon"]).all()
    assert (good["hemisphere"] == good["lon"].map(lambda x: "west" if x < 0 else "east")).all()
    assert rec.loc[rec["coord_status"] != "ok", "hemisphere"].isna().all()


def test_maps_draw_only_geographic_quality_coordinates(tmp_path):
    """P0-4: centroid and zero coordinates keep lat/lon in records.csv but are not mapped."""
    import geopandas as gpd
    from shapely.geometry import box

    shp = tmp_path / "halves.shp"
    gpd.GeoDataFrame({"half": ["west", "east"]},
                     geometry=[box(-180, -90, 0, 90), box(0, -90, 180, 90)],
                     crs="EPSG:4326").to_file(shp)
    from boldgenotyper import run_figures
    with mock.patch.object(run_figures, "_draw_map", wraps=run_figures._draw_map) as spy:
        ok, out = run_synthetic(tmp_path, skip_geo=False, shapefile_path=shp,
                                shapefile_field="half", geo_category="hemisphere")
    assert ok and spy.called
    status = _records(out).set_index("processid")["coord_status"]
    drawn = {pid for call in spy.call_args_list for pid in call.args[0]["processid"]}
    assert drawn and (status.reindex(sorted(drawn)) == "ok").all()
    assert (status == "centroid").any()           # the fixture has centroid records ...
    rec = _records(out)
    assert rec.loc[rec["coord_status"] == "centroid", "lat"].notna().all()   # ... which keep lat/lon


def test_divergence_readme_has_no_threshold_interpretations(synthetic_run):
    """P1-1: describe p-distance only; no taxonomic thresholds."""
    import re
    readme = (synthetic_run / "README.md").read_text()
    for phrase in ("cryptic", "hybridizing", "misidentification", "Typical interspecific",
                   "barcoding gap exists", "## Interpretation"):
        assert phrase.lower() not in readme.lower(), phrase
    assert not re.search(r"[<>]\s*\d+(\.\d+)?%", readme)
    assert "uncorrected p-distance" in readme


def test_label_report_replaces_contamination_outputs(synthetic_run):
    """P1-6, P1-7, P1-8: label conflicts with evidence; no contamination or confidence."""
    summary = pd.read_csv(synthetic_run / "haplotypes.csv").set_index("haplotype_id")
    # H0003 (C): two 'Genus alpha' and one genus-only 'Genus' member: compatible
    assert summary.loc["H0003", "n_label_conflicts"] == 0
    assert summary.loc["H0003", "majority_fraction"] == 1.0
    assert (synthetic_run / "label_conflicts.csv").exists()
    ann = _records(synthetic_run)
    assert {"label", "label_conflict", "depositor_keyword_hits"} <= set(ann.columns)
    assert "misidentification_confidence" not in ann.columns


@pytest.fixture(scope="module")
def extended_run(tmp_path_factory):
    ok, out = run_synthetic(tmp_path_factory.mktemp("synthetic_extended"), extended=True)
    assert ok
    return out


def test_default_run_writes_no_gap_outputs(synthetic_run):
    """P1-2: no barcoding-gap outputs by default."""
    assert not any("gap" in f or "species_divergence" in f for f in _files(synthetic_run))


def test_extended_run_writes_label_based_gap(extended_run):
    gap = pd.read_csv(extended_run / "extended" / "divergence" / "label_based_gap.csv")
    assert set(gap["label"]) == {"Genus alpha", "Genus beta", "Genus gamma"}
    assert "barcoding_gap" not in gap.columns and "label_based_gap" in gap.columns
    assert not any("within_vs_between" in f for f in _files(extended_run))
    assert "extended/labels/depositor_keyword_hits.csv" in _files(extended_run)
    d = json.loads((extended_run / "run_provenance.json").read_text())
    assert d["config"]["extended"] is True and d["warnings"] == []


def test_default_metadata_output_is_coverage_only(synthetic_run):
    """P1-3, P1-4, P1-5: cross-tabs, temporal and first-sampled outputs only with --extended."""
    assert not (synthetic_run / "extended").exists()
    summary = pd.read_csv(synthetic_run / "metadata_coverage.csv", keep_default_na=False).set_index("field")
    # blank, whitespace-only and "U" sex values are missing, not a ": N" category
    assert summary.loc["sex", "n_with_value"] == 10
    assert summary.loc["sex", "top_values"] == "Male: 5; Female: 5"
    for values in summary["top_values"]:
        assert not any(v.strip().startswith(":") for v in values.split(";"))


def test_extended_metadata_outputs(extended_run):
    md = extended_run / "extended" / "metadata"
    sex = pd.read_csv(md / "sex_tabulation.csv")
    assert set(sex["value"]) == {"Male", "Female"}
    life = pd.read_csv(md / "life_stage_tabulation.csv")
    assert set(life["value"]) == {"Adult", "Immature"}
    first = pd.read_csv(md / "haplotype_first_sampled.csv")
    assert "first_sampled" in first.columns and "first_collection_date" not in first.columns
    assert not any("emergence" in f for f in _files(extended_run))
    assert (md / "temporal_trends.csv").exists()
    readme = (extended_run / "README.md").read_text()
    assert "reflect when and where specimens were sampled" in readme and "emergence" not in readme.lower()


def test_species_diversity_only_with_extended(synthetic_run, extended_run):
    """P1-9, P1-10."""
    for run in (synthetic_run, extended_run):
        assert not any("geographic_range" in f for f in _files(run))
    div = pd.read_csv(extended_run / "extended" / "species" / "species_diversity.csv")
    assert "simpson_index" not in div.columns and set(div["note"]) == {"label-based, sampling-dependent"}


def test_no_metadata_summary_only_drops_cross_tabs(tmp_path):
    """--no-metadata-summary restricts --extended; coverage is always written (Task 3.1)."""
    ok, out = run_synthetic(tmp_path, metadata_summary_enabled=False, extended=True)
    assert ok
    files = _files(out)
    assert "metadata_coverage.csv" in files
    assert "extended/metadata/haplotype_first_sampled.csv" in files
    assert not any(f.endswith("_tabulation.csv") for f in files)


def test_sweep_reads_v2_run_folder(synthetic_run):
    """P2: the sweep reads haplotypes.csv and the QC-passed rows of records.csv."""
    from boldgenotyper import parameter_sweep
    m = parameter_sweep.extract_run_metrics(synthetic_run, ORGANISM, 0.01)
    assert m["n_haplotypes"] == 4 and m["n_samples"] == 23
    assert 0 < m["pct_assigned"] <= 100 and "pct_label_conflict" in m
    assert 0 < m["mean_majority_fraction"] <= 1
    tracking = parameter_sweep.track_group_membership({0.01: synthetic_run}, ORGANISM, [0.01])
    assert len(tracking) == 23
    md5s = set(pd.read_csv(synthetic_run / "haplotypes.csv")["seq_md5"])
    assert set(tracking["md5_0.01"]) - {""} <= md5s


def test_plot_kit_and_exports_without_duplicates(tmp_path):
    """P2: plots/ and exports/ only on request, with no copies of the FASTA or tree."""
    ok, out = run_synthetic(tmp_path, plot_kit=True, export_formats=["all"])
    assert ok
    files = _files(out)
    plots = {f for f in files if f.startswith("plots/")}
    assert {"plots/plot_config.yaml", "plots/README.md", "plots/scripts/regenerate_all.py"} <= plots
    for dup in ("plots/data/consensus_sequences.fasta", "plots/data/tree_data.nwk",
                "plots/plot_customization_config.yaml", "plots/scripts/regenerate_identity.py"):
        assert dup not in plots
    exports = {f for f in files if f.startswith("exports/")}
    assert "exports/generic/haplotype_membership.csv" in exports
    assert "exports/generic/haplotype_consensus.fasta" not in exports
    assert not [f for f in files if f.endswith(".fasta") and f not in DEFAULT_FILES]


def _undocumented(out):
    from boldgenotyper import run_readme
    return run_readme.write_run_readme(out, ORGANISM)[1]


def test_run_readme_defines_every_file_and_column(synthetic_run, extended_run, tmp_path):
    """P2: README.md lists every file and defines every CSV column."""
    assert _undocumented(synthetic_run) == []
    assert _undocumented(extended_run) == []
    ok, full = run_synthetic(tmp_path, extended=True, plot_kit=True, export_formats=["all"],
                             phylogenetic__build_tree=True)
    assert ok and _undocumented(full) == []
    readme = (full / "README.md").read_text()
    for f in _files(full):
        assert f"`{f}`" in readme or "more" in readme, f


def test_run_readme_definitions(synthetic_run):
    """P0-2, P0-9, P0-4 definitions in the README."""
    readme = (synthetic_run / "README.md").read_text()
    from boldgenotyper.record_status import STATUSES
    for status in STATUSES:
        assert f"- `{status}`:" in readme
    for col in ("n_species_labels_input", "n_species_confident", "pct_assigned"):
        assert f"| `{col}` |" in readme
    assert "Country- and province-centroid coordinates are kept in lat/lon but are not assigned a region" in readme
    assert "| `ocean_basin` |" in readme


def test_report_links_figures_and_describes_this_workflow(tmp_path):
    """P2: report.html links figures in the run folder; Methods describe the ESV workflow."""
    ok, out = run_synthetic(tmp_path, report=True, extended=True, phylogenetic__build_tree=True)
    assert ok
    html = (out / "report.html").read_text()
    assert "data:image" not in html
    assert (out / "report.html").stat().st_size < 10 * 1024 * 1024
    links = set(re.findall(r'(?:src|href)="([^"#:]+)"', html))
    assert "phylogeny/tree.svg" in links
    for link in links:
        assert (out / link).is_file(), link
    methods = html[html.index('id="section-methods"'):]
    for old in ("clustering threshold", "Similarity threshold", "Tie detection", "edlib", "Consensus genotypes"):
        assert old.lower() not in methods.lower(), old
    assert "core region" in methods and "FastTree" in methods


def test_log_has_no_glyphs(extended_run):
    """P3: plain log lines (no check marks or arrows)."""
    log = (extended_run / f"{ORGANISM}_pipeline.log").read_text()
    for glyph in ("✓", "→", "⊘", "⚠"):
        assert glyph not in log, glyph
    assert "Identity column not found" not in log


def _section(html, section_id):
    start = html.index(f'id="{section_id}"')
    end = html.find('<div class="section"', start + 1)
    return html[start:end if end != -1 else len(html)]


def test_methods_tab_reports_this_runs_values(tmp_path):
    """The Methods tab uses the run's final values, not the provenance written at startup."""
    ok, out = run_synthetic(tmp_path, report=True, phylogenetic__build_tree=True)
    assert ok
    methods = _section((out / "report.html").read_text(encoding="utf-8"), "section-methods")
    # every value is recorded except the BOLD download date (not given here)
    rows = [r for r in re.findall(r"<tr>(.*?)</tr>", methods, re.S) if "download date" not in r]
    assert rows and not any("not recorded" in r for r in rows)
    assert "No tree was built" not in methods
    assert " tips, " in methods


def test_distance_matrix_keyed_by_haplotype_id(synthetic_run):
    d = pd.read_csv(synthetic_run / "divergence" / "pairwise_p_distance.csv", index_col=0)
    ids = set(pd.read_csv(synthetic_run / "haplotypes.csv")["haplotype_id"])
    assert set(d.index) == ids and set(d.columns) == ids
    assert d.index.name == "haplotype_id"


def test_haplotypes_table_has_evidence_columns(synthetic_run):
    from boldgenotyper.taxon_evidence import EVIDENCE_COLUMNS
    haps = pd.read_csv(synthetic_run / "haplotypes.csv")
    assert set(EVIDENCE_COLUMNS) <= set(haps.columns)
    species_level = haps["assignment_level"] == "species"
    assert haps.loc[species_level, "suggested_sp"].isna().all()
    assert haps.loc[species_level, "suggestion_basis"].isna().all()
    assert haps["nearest_species_haplotype"].notna().any() and haps["label_sources"].notna().all()
    assert haps["tree_sister_taxon"].isna().all()          # no tree in this run
    d = json.loads((synthetic_run / "run_provenance.json").read_text())
    assert not [w for w in d["warnings"] if w["step"] == "5 taxon evidence"]


def test_curation_file_sets_haplotype_taxon(tmp_path):
    ok, base = run_synthetic(tmp_path / "base")
    assert ok
    target = pd.read_csv(base / "haplotypes.csv").iloc[0]
    cur = tmp_path / "curation.csv"
    cur.write_text(f"seq_md5,taxon,reason\n{target['seq_md5']},Genus curata,test\nffff,Genus nulla,no match\n")
    ok, out = run_synthetic(tmp_path / "cur", taxon_curation=cur, report=True)
    assert ok
    assert "Curated taxa" in (out / "report.html").read_text()
    row = pd.read_csv(out / "haplotypes.csv").set_index("haplotype_id").loc[target["haplotype_id"]]
    assert row["assigned_sp"] == "Genus curata" and row["assignment_source"] == "curated"
    assert row["haplotype_label"] == f"Genus curata {target['haplotype_id']}"
    applied = json.loads((out / "run_provenance.json").read_text())["applied"]["taxon_curation"]
    assert applied["n_applied"] == 1 and applied["n_not_applied"] == 1
    assert applied["md5"] and applied["not_applied"][0]["row"] == 3
    # every labeled member now differs from the curated taxon
    records = pd.read_csv(out / "records.csv")
    members = records[(records["haplotype_id"] == target["haplotype_id"]) & records["label"].notna()]
    assert len(members) and members["label_conflict"].all()
    conflicts = pd.read_csv(out / "label_conflicts.csv")
    rows = conflicts[conflicts["haplotype_id"] == target["haplotype_id"]]
    assert len(rows) == len(members)
    assert (rows["haplotype_taxon"] == "Genus curata").all() and (rows["assignment_source"] == "curated").all()


def test_invalid_curation_file_stops_the_run(tmp_path):
    cur = tmp_path / "bad.csv"
    cur.write_text("seq_md5,taxon\nabc,Genus x\n")
    ok, out = run_synthetic(tmp_path, taxon_curation=cur)
    assert not ok and not (out / "records.csv").exists()


@needs_goas
@needs_fasttree
def test_run_figures_and_manifests(tmp_path):
    ok, out = run_synthetic(tmp_path, skip_geo=False, extended=True, report=True,
                            phylogenetic__build_tree=True)
    assert ok
    files = _files(out)
    for f in ("phylogeny/tree.svg", "phylogeny/tree_overview.svg", "phylogeny/manifest.json",
              "geography/distribution_map.svg", "geography/distribution_map_by_species.pdf",
              "geography/manifest.json", "extended/geography/region_abundance_total.svg",
              "extended/metadata/figures/sex.svg", "extended/metadata/figures/records_by_year.svg",
              "extended/metadata/figures/manifest.json"):
        assert f in files, f
    assert not [f for f in files if f.endswith(".png") and not f.startswith("plots/")]
    assert not [f for f in files if "faceted" in f or "geography/per_" in f]   # old figure folders
    manifest = json.loads((out / "geography" / "manifest.json").read_text())
    for path in manifest["distribution_map"]["pages"].values():
        assert (out / "geography" / path).exists()
    html = (out / "report.html").read_text()
    assert "geography/distribution_map.svg" in html and "phylogeny/tree_overview.svg" in html
    d = json.loads((out / "run_provenance.json").read_text())
    assert d["warnings"] == []


@needs_goas
def test_default_run_writes_tree_and_map_only(tmp_path):
    ok, out = run_synthetic(tmp_path, skip_geo=False)
    assert ok and "geography/distribution_map.svg" in _files(out)
    assert not (out / "extended").exists()


@needs_fasttree
def test_tree_is_built_before_the_taxon_steps(tmp_path):
    ok, out = run_synthetic(tmp_path, phylogenetic__build_tree=True)
    assert ok
    assert (out / "phylogeny" / "tree.nwk").exists() and (out / "phylogeny" / "tree_labeled.nwk").exists()
    log = (out / f"{ORGANISM}_pipeline.log").read_text()
    assert log.index("5.3: Building phylogenetic tree") < log.index("PHASE 6:")
    assert "8.2: Creating relabeled tree" in log


def test_evidence_sets_a_split_haplotypes_taxon(tmp_path):
    ok, out = run_synthetic(tmp_path, split_labels=True)
    assert ok
    h = pd.read_csv(out / "haplotypes.csv").set_index("haplotype_id").loc["H0003"]
    assert h["assigned_sp"] == "Genus alpha" and h["assignment_level"] == "species"
    assert h["assignment_source"] == "evidence" and h["haplotype_label"] == "Genus alpha H0003"
    assert h["assignment_notes"] == ("evidence: distance, BIN agree; labels: "
                                     "majority 0.67 not above 0.70; fell back to genus")
    # evidence was computed before it was applied: the suggestion is kept
    assert h["suggested_sp"] == "Genus alpha" and h["suggestion_basis"] == "distance, BIN agree"
    rec = pd.read_csv(out / "records.csv").set_index("processid")
    assert rec.loc["SYN017", "assigned_sp"] == "Genus alpha"
    conflicts = pd.read_csv(out / "label_conflicts.csv").set_index("processid")
    assert conflicts.loc["SYN019", "haplotype_taxon"] == "Genus alpha"
    assert conflicts.loc["SYN019", "assignment_source"] == "evidence"
    applied = json.loads((out / "run_provenance.json").read_text())["applied"]["taxon_evidence"]
    assert applied == {"apply": True, "tree": False, "n_applied": 1,
                       "n_suggested_not_applied": 0, "suggestions_file": None}
    assert not (out / "taxon_curation_suggestions.csv").exists()


def test_no_evidence_taxa_writes_suggestions_that_curate(tmp_path):
    ok, out = run_synthetic(tmp_path / "off", split_labels=True, taxonomy__apply_evidence=False)
    assert ok
    h = pd.read_csv(out / "haplotypes.csv").set_index("haplotype_id").loc["H0003"]
    assert h["assignment_level"] == "genus" and h["assignment_source"] == "labels"
    sug = pd.read_csv(out / "taxon_curation_suggestions.csv")
    assert sug["haplotype_id"].tolist() == ["H0003"]
    assert sug.loc[0, "taxon"] == "Genus alpha" and sug.loc[0, "why_not_applied"] == "evidence taxa off"
    applied = json.loads((out / "run_provenance.json").read_text())["applied"]["taxon_evidence"]
    assert applied["apply"] is False and applied["n_applied"] == 0
    assert applied["suggestions_file"] == "taxon_curation_suggestions.csv"
    assert _undocumented(out) == []
    ok, cur = run_synthetic(tmp_path / "cur", split_labels=True, taxonomy__apply_evidence=False,
                            taxon_curation=out / "taxon_curation_suggestions.csv")
    assert ok
    h = pd.read_csv(cur / "haplotypes.csv").set_index("haplotype_id").loc["H0003"]
    assert h["assigned_sp"] == "Genus alpha" and h["assignment_source"] == "curated"


@needs_fasttree
def test_evidence_taxon_labels_the_tree(tmp_path):
    ok, out = run_synthetic(tmp_path, split_labels=True, phylogenetic__build_tree=True)
    assert ok
    assert "Genus alpha H0003" in (out / "phylogeny" / "tree_labeled.nwk").read_text()
    h = pd.read_csv(out / "haplotypes.csv").set_index("haplotype_id").loc["H0003"]
    assert h["assignment_source"] == "evidence" and h["suggestion_basis"] == "distance, BIN, tree agree"
    assert json.loads((out / "run_provenance.json").read_text())["applied"]["taxon_evidence"]["tree"] is True


def test_no_evidence_taxa_flag():
    assert cli._build_main_parser().parse_args(["in.tsv", "--no-evidence-taxa"]).no_evidence_taxa
    assert config.get_default_config().taxonomy.apply_evidence is True


# --- Input rows: marker selection and malformed rows (pre-tag A1, A2) -----------

def _insert_other_marker_row(tsv):
    """A 12S row for SYN002, placed before its COI row."""
    lines = tsv.read_text().splitlines()
    header = lines[0].split("\t")
    idx = next(i for i, line in enumerate(lines) if line.startswith("SYN002\t"))
    row = lines[idx].split("\t")
    row[header.index("marker_code")] = "12S"
    row[header.index("nuc")] = "TTTTAAAACCCCGGGG" * 12
    lines.insert(idx, "\t".join(row))
    tsv.write_text("\n".join(lines) + "\n")


def _append_malformed_row(tsv):
    lines = tsv.read_text().splitlines()
    lines.append(lines[2] + "\textra1\textra2")
    tsv.write_text("\n".join(lines) + "\n")
    return len(lines)          # file line number of the malformed row


def test_other_marker_row_first_keeps_coi(tmp_path):
    ok, out = run_synthetic(tmp_path, edit_tsv=_insert_other_marker_row)
    assert ok
    assert _records(out).set_index("processid").loc["SYN002", "status"] == "assigned"
    s = pd.read_csv(out / "summary.csv").iloc[0]
    assert s["n_other_marker_rows"] == 1
    assert s["n_duplicate_rows"] == 1          # the planted duplicate


def test_malformed_row_stops_run(tmp_path):
    ok, out = run_synthetic(tmp_path, edit_tsv=_append_malformed_row)
    assert ok is False
    assert not (out / "records.csv").exists()
    assert "--allow-malformed-rows" in (out / f"{ORGANISM}_pipeline.log").read_text()


def test_allow_malformed_rows_counts(tmp_path):
    lines = {}
    ok, out = run_synthetic(tmp_path, allow_malformed_rows=True,
                            edit_tsv=lambda t: lines.setdefault("n", _append_malformed_row(t)))
    assert ok
    assert pd.read_csv(out / "summary.csv").iloc[0]["n_malformed_rows"] == 1
    rows = json.loads((out / "run_provenance.json").read_text())["applied"]["input_rows"]
    assert rows["malformed_lines"] == [lines["n"]] and rows["marker"] == "COI-5P"


def test_main_parser_input_flags():
    args = cli._build_main_parser().parse_args(["in.tsv", "--marker", "COI-3P", "--allow-malformed-rows"])
    assert args.marker == "COI-3P" and args.allow_malformed_rows is True
    args = cli._build_main_parser().parse_args(["in.tsv"])
    assert args.marker == "COI-5P" and args.allow_malformed_rows is False


# --- Spelling variants (pre-tag A3) ----------------------------------------------

def _misspell_one_alpha(tsv):
    """SYN003 (haplotype A, Genus alpha) is deposited as 'Genus alpah'."""
    lines = tsv.read_text().splitlines()
    header = lines[0].split("\t")
    idx = next(i for i, line in enumerate(lines) if line.startswith("SYN003\t"))
    row = lines[idx].split("\t")
    row[header.index("species")] = "Genus alpah"
    lines[idx] = "\t".join(row)
    tsv.write_text("\n".join(lines) + "\n")


def test_label_variants_file(tmp_path):
    ok, out = run_synthetic(tmp_path, edit_tsv=_misspell_one_alpha)
    assert ok
    variants = pd.read_csv(out / "label_variants.csv")
    row = variants[(variants.label_a == "Genus alpha") & (variants.label_b == "Genus alpah")]
    assert len(row) == 1 and row.iloc[0]["n_shared_haplotypes"] == 1


def _alpah_map(tmp_path, text="from,to,reason\nGenus alpah,Genus alpha,misspelling\n"):
    path = tmp_path / "label_map.csv"
    path.write_text(text)
    return path


def test_label_map_merges_variant(tmp_path):
    ok, out = run_synthetic(tmp_path, edit_tsv=_misspell_one_alpha, label_map=_alpah_map(tmp_path))
    assert ok
    haps = pd.read_csv(out / "haplotypes.csv")
    assert not haps["labels"].fillna("").str.contains("Genus alpah").any()
    rec = _records(out).set_index("processid")
    assert rec.loc["SYN003", "label_mapped_from"] == "Genus alpah"
    assert rec.loc["SYN003", "label"] == "Genus alpha"
    applied = json.loads((out / "run_provenance.json").read_text())["applied"]["label_map"]
    assert applied["md5"] and applied["n_records_changed"] == 1 and applied["unused"] == []
    variants = pd.read_csv(out / "label_variants.csv")
    row = variants[(variants.label_a == "Genus alpha") & (variants.label_b == "Genus alpah")]
    assert len(row) == 1 and bool(row.iloc[0]["in_label_map"]) is True


def test_invalid_label_map_stops_before_output(tmp_path):
    bad = _alpah_map(tmp_path, "from,to,reason\nGenus alpah,Genus alpha,x\nGenus alpha,Genus beta,y\n")
    ok, out = run_synthetic(tmp_path, edit_tsv=_misspell_one_alpha, label_map=bad)
    assert ok is False
    assert not (out / "records.csv").exists()


# --- Error-filter diagnostics (pre-tag A4) ---------------------------------------

def _relabel_syn021(tsv):
    """SYN021 (planted assigned_corrected) gets a species found in no retained haplotype."""
    lines = tsv.read_text().splitlines()
    header = lines[0].split("\t")
    idx = next(i for i, line in enumerate(lines) if line.startswith("SYN021\t"))
    row = lines[idx].split("\t")
    row[header.index("species")] = "Genus delta"
    lines[idx] = "\t".join(row)
    tsv.write_text("\n".join(lines) + "\n")


def test_error_filter_diagnostics_in_run(tmp_path):
    ok, out = run_synthetic(tmp_path, edit_tsv=_relabel_syn021)
    assert ok
    rec = _records(out).set_index("processid")
    assert rec.loc["SYN021", "label_in_retained"] in (False, "False")
    assert rec.loc["SYN022", "label_in_retained"] in (True, "True")
    s = pd.read_csv(out / "summary.csv").iloc[0]
    assert s["n_singleton_label_absent"] == 1
    assert (out / "sources.csv").exists()
    assert "may have removed species" in (out / f"{ORGANISM}_pipeline.log").read_text()


def test_corrected_share_warning(tmp_path, monkeypatch):
    from boldgenotyper import error_filter_diagnostics
    monkeypatch.setattr(error_filter_diagnostics, "CORRECTED_SHARE_WARNING", 0.01)
    ok, out = run_synthetic(tmp_path)
    assert ok
    assert "corrected singletons" in (out / f"{ORGANISM}_pipeline.log").read_text()
