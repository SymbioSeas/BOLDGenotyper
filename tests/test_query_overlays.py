"""Query metadata, regions and figures against a synthetic reference run (requires MAFFT)."""

import json
import shutil

import pandas as pd
import pytest

from boldgenotyper import cli, query_overlays as qv
from tests.synthetic_data import A, D
from tests.test_pipeline_synthetic import run_synthetic

pytestmark = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")
needs_fasttree = pytest.mark.skipif(shutil.which("FastTree") is None and shutil.which("fasttree") is None,
                                    reason="FastTree not installed")

QUERIES = {"win": A[100:320], "gam": D[100:320], "none": "ACGTTGCA" * 30}
META = ("query_id,coord,collection_date,region\n"
        'win,"[10, -20]",2019-05-01,\n'
        'gam,"[-5, 140]",2021-07-01,\n'
        "none,,2020-01-01,\n"
        'extra,"[1, 1]",2020-01-01,\n')


def test_pruned_genus_tree_keeps_its_own_stem():
    """Pruning alone adds every branch from the old root into one stem (5.5 here)."""
    from io import StringIO
    from Bio import Phylo
    tree = Phylo.read(StringIO("(((A:1,B:1)0.9:0.5,E:3)0.8:5,(C:1,D:1)0.7:2);"), "newick")
    sub = qv._pruned_to(tree, {"A", "B"})
    assert sorted(t.name for t in sub.get_terminals()) == ["A", "B"]
    assert sub.root.branch_length == 0.5 and sub.root.confidence == 0.9
    assert tree.count_terminals() == 5


def test_timeline_marks_each_query_year_once(tmp_path):
    """Queries sharing a year get one line and a count; the caption names them by year."""
    ids = ["q1", "q2", "q3", "q4"]
    calls = pd.DataFrame({"query_id": ids, "call_level": "haplotype", "equivalence_set": "H0001",
                          "species_call": "Genus alpha"})
    meta = pd.DataFrame({"query_id": ids, "collection_year": [2020, 2020, 2021, None]})
    records = pd.DataFrame({"haplotype_id": "H0001", "collection_year": [2010, 2012, 2012]})
    ref = qv.ReferenceRun(tmp_path, records, pd.DataFrame(), {}, None, None)
    qv.plot_query_timeline(calls, meta, ref, tmp_path)
    svg = (tmp_path / "query_timeline" / "H0001.svg").read_text()
    assert svg.count("stroke-dasharray") == 2
    assert "2020: q1, q2" in svg and "2021: q3" in svg
    assert "Queries without a collection date: q4" in svg


def test_queries_at_one_site_share_a_star_and_label():
    queries = [("q1", 10, -20), ("q2", 10.0, -20.0), ("q3", 5, 5)] + [(f"x{i}", 1, 1) for i in range(4)]
    assert qv._query_sites(queries) == [(10.0, -20.0, "q1, q2"), (5.0, 5.0, "q3"),
                                        (1.0, 1.0, "x0, x1, x2\nx3")]


def _halves(path):
    import geopandas as gpd
    from shapely.geometry import box
    gpd.GeoDataFrame({"half": ["west", "east"]},
                     geometry=[box(-180, -90, 0, 90), box(0, -90, 180, 90)],
                     crs="EPSG:4326").to_file(path)
    return path


@pytest.fixture(scope="module")
def ref_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("qref")
    ok, out = run_synthetic(d, skip_geo=False, shapefile_path=_halves(d / "halves.shp"),
                            shapefile_field="half", geo_category="hemisphere", report=True,
                            phylogenetic__build_tree=True)
    assert ok
    return out


def _run_query(tmp_path, ref, *extra, meta=META):
    fasta = tmp_path / "q.fasta"
    fasta.write_text("".join(f">{k}\n{v}\n" for k, v in QUERIES.items()))
    m = tmp_path / "meta.csv"
    m.write_text(meta)
    out = tmp_path / "out"
    assert cli.main_query(["--query", str(fasta), "--haplotypes", str(ref), "--output", str(out),
                           "--query-metadata", str(m), *extra]) == 0
    return out


def test_load_query_metadata(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("Query_ID\tlat\tlon\tcollection_date\nwin\t10\t-20\t2019-05-01\nextra\t1\t1\t\n")
    meta = qv.load_query_metadata(p, ["win", "gam"])
    assert meta["query_id"].tolist() == ["win", "gam"]
    assert meta.loc[0, "coord_status"] == "ok" and meta.loc[0, "collection_year"] == 2019
    assert meta.loc[1, "coord_status"] == "missing"
    assert meta.attrs["unmatched_ids"] == ["extra"]
    p.write_text("query_id\nwin\nwin\n")
    with pytest.raises(qv.QueryMetadataError):
        qv.load_query_metadata(p, ["win"])
    p.write_text("sample\nwin\n")
    with pytest.raises(qv.QueryMetadataError):
        qv.load_query_metadata(p, ["win"])


@needs_fasttree
def test_regions_from_metadata_coordinates_or_none(tmp_path, ref_run, caplog):
    ref = qv.load_reference(ref_run)
    assert ref.region_col == "hemisphere" and ref.tree_path is not None
    p = tmp_path / "m.csv"
    p.write_text('query_id,coord,region\na,"[10, -20]",\nb,"[10, 20]",given\nc,,\n')
    meta = qv.assign_query_regions(qv.load_query_metadata(p, ["a", "b", "c"]), ref)
    assert meta[["region", "region_source"]].values.tolist() == [
        ["west", "coordinates"], ["given", "metadata"], ["", ""]]
    ref.provenance["applied"]["custom_shapefile"] = str(tmp_path / "moved.shp")
    meta = qv.assign_query_regions(qv.load_query_metadata(p, ["a"]), ref)
    assert meta.loc[0, "region"] == "" and "moved.shp" in caplog.text


@needs_fasttree
def test_query_regions_table(tmp_path, ref_run):
    out = _run_query(tmp_path, ref_run)
    t = pd.read_csv(out / "query_figures" / "query_regions.csv")
    assert list(t.columns) == qv.QUERY_REGION_COLUMNS
    win = t[(t["equivalence_set"] == "H0001;H0002") & (t["region"] == "west")].iloc[0]
    assert win["n_queries"] == 1 and win["query_ids"] == "win" and win["n_reference_records"] > 0
    prov = json.loads((out / "query_provenance.json").read_text())
    assert prov["query_metadata"]["n_matched"] == 3 and prov["query_metadata"]["unmatched_ids"] == ["extra"]


@needs_fasttree
def test_figures_with_a_tree_and_metadata(tmp_path, ref_run):
    out = _run_query(tmp_path, ref_run)
    figs = out / "query_figures"
    for name in ("query_tree", "query_map", "query_regions", "query_timeline"):
        assert (figs / f"{name}.pdf").exists(), name
        assert list((figs / name).glob("*.svg")), name
    tree_svg = (figs / "query_tree" / "Genus.svg").read_text()
    assert "win" in tree_svg and "Not placed: none" in tree_svg
    assert {p.stem for p in (figs / "query_map").glob("*.svg")} == {"Genus_alpha", "Genus_gamma"}
    prov = json.loads((out / "query_provenance.json").read_text())
    assert prov["warnings"] == [] and "query_tree" in prov["figures"]


def test_figures_without_tree_or_metadata(tmp_path):
    ok, ref = run_synthetic(tmp_path / "ref")                  # no tree, no regions
    assert ok
    fasta = tmp_path / "q.fasta"
    fasta.write_text(f">win\n{QUERIES['win']}\n")
    out = tmp_path / "out"
    assert cli.main_query(["--query", str(fasta), "--haplotypes", str(ref), "--output", str(out)]) == 0
    figs = out / "query_figures"
    assert not (figs / "query_tree.pdf").exists() and not (figs / "query_regions.pdf").exists()
    assert (figs / "query_map.pdf").exists() and (figs / "query_timeline.pdf").exists()
    prov = json.loads((out / "query_provenance.json").read_text())
    assert prov["warnings"] == []
