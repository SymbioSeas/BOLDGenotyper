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
