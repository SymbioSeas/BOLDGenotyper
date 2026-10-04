"""query_report.html (Phase 10)."""

import re
import shutil

import pytest

from tests.test_query_overlays import QUERIES, _run_query, needs_fasttree, ref_run  # noqa: F401

pytestmark = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")


@needs_fasttree
def test_query_report(tmp_path, ref_run):
    out = _run_query(tmp_path, ref_run)
    html = (out / "query_report.html").read_text()
    assert html.count('class="query-card"') == len(QUERIES)
    assert html.count('<tr class="call-row"') == len(QUERIES)
    targets = re.findall(r'(?:src|href)="([^"#:]+)"', html)
    assert any(t.endswith(".svg") for t in targets)
    for target in targets:
        assert (out / target).exists(), target
    assert any(t.endswith("report.html") and t != "query_report.html" for t in targets)
    assert "Genus alpha" in html
    card = html[html.index('id="q-none"'):]
    assert "no eligible match" in card[:card.index("</details>")]
