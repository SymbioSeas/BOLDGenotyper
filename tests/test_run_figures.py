"""Run figures from records and the tree, without a pipeline run (Phase 11)."""

import matplotlib
matplotlib.use("Agg")
import pytest

from boldgenotyper import tree_figures as tf
from tests.synthetic_data import figure_dataset


@pytest.fixture(scope="module")
def data():
    return figure_dataset()


def test_figure_dataset(data):
    records, haps, nwk = data
    assert len(haps) == 46 and set(haps["assigned_sp"]) == {
        "Alpha one", "Alpha two", "Beta three", "Beta four", "Gamma"}
    assert nwk.count("H0") == 46


def test_tree_figures(tmp_path, data):
    records, haps, nwk = data
    tree = tmp_path / "tree.nwk"
    tree.write_text(nwk)
    entries = tf.write_tree_figures(tree, haps, tmp_path)
    assert entries["tree"]["overview"] == ["tree.pdf", "tree.svg"]
    overview = (tmp_path / "tree_overview.svg").read_text()
    assert "Alpha one (12 haplotypes" in overview and "Beta four (1 haplotype" in overview
    assert "Alpha two H00" in overview                       # not a clade: left expanded
    assert "Not monophyletic: Alpha two" in overview
    assert set(entries["tree_by_genus"]["pages"]) == {"Alpha", "Beta", "Gamma"}
    beta = (tmp_path / entries["tree_by_genus"]["pages"]["Beta"]).read_text()
    assert "Alpha two" in beta                               # the genus's clade holds another genus
    assert tf.full_tree_size(46) == (10, pytest.approx(0.14 * 46 + 1.0))


from boldgenotyper import run_figures as rf


def test_geography_figures(tmp_path, data):
    records, haps, _ = data
    out = rf.write_geography_figures(records, "ocean_basin", tmp_path / "geography",
                                     tmp_path / "extended" / "geography")
    m = out["geography"]["distribution_map"]
    overview = (tmp_path / "geography" / "distribution_map.svg").read_text()
    assert "Alpha one" in overview and "H0001" not in overview          # 46 haplotypes: by taxon
    assert set(m["pages"]) == {"Alpha one", "Alpha two", "Beta three", "Beta four", "Gamma"}
    ext = out["extended/geography"]
    assert set(ext) == {"region_abundance_relative", "region_abundance_total"}
    assert set(ext["region_abundance_total"]["pages"]) == {"Alpha one", "Alpha two", "Beta three", "Gamma"}
    beta = (tmp_path / "extended" / "geography" / ext["region_abundance_total"]["pages"]["Beta three"]).read_text()
    assert "Other (3 haplotypes)" in beta                                # 14 haplotypes: 11 + Other


def test_small_run_overview_is_by_haplotype(tmp_path, data):
    records, _, _ = data
    small = records[records["assigned_sp"] == "Beta three"]
    rf.write_geography_figures(small, "ocean_basin", tmp_path, None)
    assert "Beta three H00" in (tmp_path / "distribution_map.svg").read_text()


from boldgenotyper import metadata_figures as mf


def test_region_facet(data):
    records, _, _ = data
    col, caption = mf.region_facet(records, "ocean_basin")
    n = records["ocean_basin"].notna().sum()
    assert col == "ocean_basin" and f"{n} of {len(records)} assigned records" in caption
    one = records.assign(ocean_basin="North Basin")
    assert mf.region_facet(one, "ocean_basin")[0] == "country/ocean"
    col, caption = mf.region_facet(records, None)
    assert col == "country/ocean" and "no region analysis" in caption


def test_metadata_figures(tmp_path, data):
    records, _, _ = data
    out = mf.write_metadata_figures(records, ["sex", "habitat"], "ocean_basin", tmp_path)
    assert {"sex", "sex_by_region", "habitat", "records_by_year", "sampling_timeline",
            "records_by_year_by_region", "metadata_coverage"} <= set(out)
    assert out["sex"]["pages_pdf"] == "sex_by_haplotype.pdf"
    assert set(out["records_by_year"]["pages"]) == {"Alpha one", "Alpha two", "Beta three", "Gamma"}
    assert "assigned records have an ocean_basin value" in (tmp_path / "sex_by_region.svg").read_text()
    assert not list(tmp_path.rglob("*.png"))
    none = mf.write_metadata_figures(records, ["sex"], "ocean_basin", tmp_path / "t",
                                     cross_tabs=False, temporal=False)
    assert set(none) == {"metadata_coverage"}
