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
