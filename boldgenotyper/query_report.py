"""query_report.html: the calls, one evidence card per query, the figures and
the methods of a boldgenotyper-query run, in one self-contained page."""

import logging
import os
from html import escape
from pathlib import Path
from typing import Dict, List

import pandas as pd

from .query_outputs import methods_sentence

logger = logging.getLogger(__name__)

_CSS = """
:root { --ink: #1f2328; --muted: #59636e; --line: #d1d9e0; --accent: #c0392b; --bg: #ffffff;
        --panel: #f6f8fa; }
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: var(--ink);
       background: var(--bg); max-width: 1100px; margin: 0 auto; padding: 16px; line-height: 1.45; }
h1 { font-size: 1.6em; margin-bottom: 0.2em; }
h2 { border-bottom: 1px solid var(--line); padding-bottom: 0.2em; margin-top: 1.6em; }
.meta, .note { color: var(--muted); }
nav a { margin-right: 1em; }
table { border-collapse: collapse; width: 100%; font-size: 0.9em; margin: 0.6em 0; }
th, td { border: 1px solid var(--line); padding: 4px 6px; text-align: left; vertical-align: top; }
th { background: var(--panel); }
.wrap { overflow-x: auto; }
details.query-card { border: 1px solid var(--line); border-radius: 6px; margin: 0.6em 0;
                     padding: 0.4em 0.8em; background: var(--panel); }
details.query-card summary { cursor: pointer; font-weight: 600; }
pre { background: var(--bg); border: 1px solid var(--line); padding: 8px; overflow-x: auto;
      font-size: 0.8em; }
figure { margin: 1em 0; }
figure img { max-width: 100%; border: 1px solid var(--line); background: #fff; }
.flag { color: var(--accent); }
"""

_CALL_TABLE = ["query_id", "call_level", "equivalence_set", "species_call", "call_basis",
               "species_support", "differences", "identity_pct", "query_coverage", "flags"]
_CANDIDATE_TABLE = ["rank", "haplotype_id", "assigned_sp", "assignment_source", "n_members",
                    "in_equivalence_set", "eligible", "differences", "identity_pct",
                    "compared_coverage"]
_FIGURE_TITLES = {
    "query_tree": "Tree: queries beside their equivalence sets' haplotypes",
    "query_map": "Map: reference records and queries, by called taxon",
    "query_regions": "Regions: reference records per region, queries marked",
    "query_timeline": "Timeline: sampling years of the matched haplotypes, query dates marked",
}


def _cell(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return escape(str(value))


def _table(df: pd.DataFrame, columns: List[str], row_class: str = "") -> str:
    cols = [c for c in columns if c in df.columns]
    head = "".join(f"<th>{escape(c)}</th>" for c in cols)
    attr = f' class="{row_class}"' if row_class else ""
    body = "".join(f"<tr{attr}>" + "".join(f"<td>{_cell(r[c])}</td>" for c in cols) + "</tr>\n"
                   for _, r in df.iterrows())
    return f'<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>\n{body}</tbody></table></div>\n'


def _call_summary(r) -> str:
    if r["call_level"] == "none":
        return "no eligible match"
    species = (f"{r['species_call']} ({r['call_basis']})" if isinstance(r["species_call"], str)
               and r["species_call"] else "no species call")
    return f"{r['call_level']} {r['equivalence_set']}; {species}"


def generate_query_report(output_dir: Path, calls: pd.DataFrame, candidates: pd.DataFrame,
                          provenance: dict, figures: Dict[str, List[Path]],
                          alignments: Dict[str, str]) -> Path:
    """Write query_report.html in output_dir and return its path."""
    output_dir = Path(output_dir)
    path = output_dir / "query_report.html"
    ref = provenance.get("reference_run")
    query = provenance.get("query") or {}

    def rel(p) -> str:
        return Path(os.path.relpath(Path(p) if Path(p).is_absolute() else output_dir / p,
                                    output_dir)).as_posix()

    parts = ["<h1>Query report</h1>",
             f'<p class="meta">{_cell(query.get("n_sequences"))} query sequence(s) from '
             f'<code>{_cell(query.get("path"))}</code>; BOLDGenotyper '
             f'{_cell(provenance.get("boldgenotyper_version"))}, {_cell(provenance.get("run_date"))}</p>',
             '<nav><a href="#summary">Summary</a><a href="#calls">Calls</a>'
             '<a href="#evidence">Evidence</a><a href="#figures">Figures</a>'
             '<a href="#reference">Reference run</a><a href="#methods">Methods</a></nav>']

    levels = calls["call_level"].value_counts()
    called = calls.loc[calls["species_call"].fillna("") != "", "species_call"].value_counts()
    summary = [f"<li>{_cell(n)} with call level <code>{_cell(level)}</code></li>"
               for level, n in levels.items()]
    summary += [f"<li>{_cell(n)} called <em>{_cell(sp)}</em></li>" for sp, n in called.items()]
    parts += ['<section id="summary"><h2>Summary</h2><ul>', *summary, "</ul>",
              '<p class="note">The result to report is <code>query_calls.csv</code>; '
              '<code>README.md</code> defines every column.</p></section>']

    rows = calls.copy()
    parts += ['<section id="calls"><h2>Calls</h2>', _table(rows, _CALL_TABLE, "call-row"), "</section>"]

    cards = []
    for _, r in calls.iterrows():
        qid = str(r["query_id"])
        flags = "".join(f'<span class="flag">{escape(f)}</span> ' for f in str(r["flags"] or "").split(";")
                        if f and f != "nan")
        cand = candidates[candidates["query_id"] == r["query_id"]] if len(candidates) else candidates
        cards.append(
            f'<details class="query-card" id="q-{escape(qid)}"><summary>{escape(qid)}: '
            f"{escape(_call_summary(r))}</summary>"
            f"<p><strong>Species support</strong> {_cell(r['species_support'])}; "
            f"<strong>pooled member labels</strong> {_cell(r['pooled_member_labels']) or '-'}</p>"
            f"<p><strong>Best match</strong> {_cell(r['differences'])} difference(s), "
            f"{_cell(r['identity_pct'])}% identity, query coverage {_cell(r['query_coverage'])}"
            f"{'; flags: ' + flags if flags else ''}</p>"
            + _table(cand, _CANDIDATE_TABLE)
            + (f"<pre>{escape(alignments[qid])}</pre>" if alignments.get(qid) else "")
            + "</details>")
    parts += ['<section id="evidence"><h2>Evidence</h2>',
              '<p class="note">One card per query: the call, the candidates '
              '(<code>query_candidates.csv</code>) and the best alignment.</p>', *cards, "</section>"]

    fig_parts = []
    for name, title in _FIGURE_TITLES.items():
        paths = figures.get(name) or []
        if not paths:
            continue
        pdf = next((p for p in paths if str(p).endswith(".pdf")), None)
        fig_parts.append(f"<h3>{escape(title)}</h3>"
                         + (f'<p><a href="{escape(rel(pdf))}">All pages (PDF)</a></p>' if pdf else ""))
        for p in paths:
            if str(p).endswith(".svg"):
                fig_parts.append(f'<figure><img src="{escape(rel(p))}" alt="{escape(title)}: '
                                 f'{escape(Path(p).stem)}"><figcaption>{escape(Path(p).stem)}'
                                 f"</figcaption></figure>")
    regions_csv = (figures.get("query_regions_table") or [None])[0]
    if regions_csv:
        fig_parts.append(f'<p><a href="{escape(rel(regions_csv))}">query_regions.csv</a></p>')
    parts += ['<section id="figures"><h2>Figures</h2>',
              *(fig_parts or ['<p class="note">No figures: the reference is not a v2 run folder.</p>']),
              "</section>"]

    if ref:
        ref_report = Path(str(ref.get("path"))) / "report.html"
        link = (f'<li><a href="{escape(rel(ref_report))}">Reference run report</a></li>'
                if ref_report.exists() else "")
        ref_items = [f"<li>Path <code>{_cell(ref.get('path'))}</code></li>",
                     f"<li>BOLDGenotyper {_cell(ref.get('boldgenotyper_version'))}, run "
                     f"{_cell(ref.get('run_date'))}</li>",
                     f"<li>run_provenance.json md5 <code>{_cell(ref.get('provenance_md5'))}</code></li>",
                     f"<li>BOLD download date {_cell(ref.get('download_date')) or 'not recorded'}</li>",
                     link]
    else:
        ref_items = [f"<li>A haplotype FASTA, not a v2 run folder: "
                     f"<code>{_cell((provenance.get('reference_fasta') or {}).get('path'))}</code></li>"]
    parts += ['<section id="reference"><h2>Reference run</h2><ul>', *ref_items, "</ul></section>",
              f'<section id="methods"><h2>Methods</h2><p>{escape(methods_sentence(provenance))}</p></section>']

    html = ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            f"<title>Query report</title>\n<style>{_CSS}</style>\n</head>\n<body>\n"
            + "\n".join(parts) + "\n</body>\n</html>\n")
    path.write_text(html, encoding="utf-8")
    logger.info(f"Wrote query report: {path}")
    return path
