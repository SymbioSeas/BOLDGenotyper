"""The HTML report (Phase 12)."""

import json
import re
import shutil

import pandas as pd
import pytest

from boldgenotyper.report import page

needs_tools = pytest.mark.skipif(shutil.which("mafft") is None or (shutil.which("FastTree") is None
                                 and shutil.which("fasttree") is None), reason="MAFFT/FastTree missing")


def test_table_truncates_with_the_true_count():
    html = page.table(pd.DataFrame({"a": range(25)}), ["a"], file="haplotypes.csv")
    assert html.count("<tr") == 11                         # header + 10 rows
    assert "15 more in <code>haplotypes.csv</code>" in html
    assert "more in" not in page.table(pd.DataFrame({"a": range(3)}), ["a"])


def test_components_escape_and_switch():
    html = page.switcher("tree", [("Genus <a>", "phylogeny/x.svg"), ("B", "phylogeny/y.svg")], "tree")
    assert "Genus &lt;a&gt;" in html and 'value="phylogeny/x.svg"' in html and "<select" in html
    assert "<select" not in page.switcher("one", [("only", "a.svg")], "x")
    two = page.two_level_switcher("f", {"sex": [("by taxon", "a.svg")], "habitat": [("by taxon", "b.svg")]}, "f")
    data = json.loads(re.search(r'<script type="application/json"[^>]*>(.*?)</script>', two, re.S).group(1))
    assert data["sex"] == [["by taxon", "a.svg"]]
    assert "width: 50.0%" in page.bars([("a", 50, ""), ("b", 100, "")])
    shell = page.render_page("T", "s", [page.Section("overview", "Overview", "<p>x</p>")])
    assert 'id="section-overview"' in shell and 'href="#section-overview"' in shell


from boldgenotyper.report import data as rdata, insights
from tests.test_pipeline_synthetic import run_synthetic
from tests.test_query_overlays import _halves


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("report_run")
    ok, out = run_synthetic(d, skip_geo=False, shapefile_path=_halves(d / "halves.shp"),
                            shapefile_field="half", geo_category="hemisphere", extended=True,
                            report=True, phylogenetic__build_tree=True)
    assert ok
    return out


@needs_tools
def test_key_findings_match_the_tables(run_dir):
    run = rdata.load_run(run_dir, "Synth")
    found = {f.key: f for f in insights.key_findings(run)}
    s = pd.read_csv(run_dir / "summary.csv").iloc[0]
    rec = pd.read_csv(run_dir / "records.csv", low_memory=False)
    haps = pd.read_csv(run_dir / "haplotypes.csv")
    assert found["assigned"].count == s["n_assigned"] + s["n_assigned_corrected"]
    assert found["mapped"].count == (rec["coord_status"] == "ok").sum()
    assert found["coordinates"].count == rec["lat"].notna().sum()
    n_conflicted = (haps["n_label_conflicts"] > 0).sum()
    assert (found["label_conflicts"].count if "label_conflicts" in found else 0) == n_conflicted
    assert "curated" not in found
    assert ("not_monophyletic" in found) == bool(run.non_monophyletic())
    assert found["dates"].text.startswith("Collection years 20")
    assert [f.key for f in insights.key_findings(run)][0] == "assigned"
