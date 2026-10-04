"""Output files: every file in the run folder, grouped by purpose, described
from run_readme.FILE_DOCS so the report and README.md cannot disagree."""

from collections import OrderedDict
from typing import Iterable

from .. import run_readme
from .data import RunData
from .page import Section, esc

START = {"README.md", "report.html", "summary.csv", "records.csv", "haplotypes.csv", "run_provenance.json"}
GROUPS = ["Start here", "Tables and logs", "Sequences and trees", "Figures", "Extended (--extended)",
          "Exports and plot kit", "Other"]


def _group(rel: str, pattern: str) -> str:
    if pattern in START:
        return "Start here"
    if rel.startswith("extended/"):
        return "Extended (--extended)"
    if rel.startswith(("exports/", "plots/")):
        return "Exports and plot kit"
    if rel.endswith((".pdf", ".svg")) or rel.endswith("manifest.json"):
        return "Figures"
    if rel.endswith((".fasta", ".nwk")):
        return "Sequences and trees"
    return "Tables and logs"


def build(run: RunData, also: Iterable[str] = ("report.html", "README.md")) -> Section:
    rels = {p.relative_to(run.run_dir).as_posix() for p in run.run_dir.rglob("*") if p.is_file()}
    rels |= set(also)
    rows = {g: OrderedDict() for g in GROUPS}
    for rel in sorted(rels):
        entry = run_readme._entry_for(rel)
        pattern, text = (entry["pattern"], entry["text"]) if entry else (rel, "")
        group = _group(rel, pattern) if entry else "Other"
        item = rows[group].setdefault(pattern, {"n": 0, "text": text, "example": rel})
        item["n"] += 1
    html = ('<p>Every file in this run folder, described as in <code>README.md</code> (which also defines '
            "every column).</p>\n")
    for group in GROUPS:
        if not rows[group]:
            continue
        body = ""
        for pattern, item in rows[group].items():
            name = (f'<a href="{esc(item["example"])}"><code>{esc(pattern)}</code></a>' if item["n"] == 1
                    and "*" not in pattern else f"<code>{esc(pattern)}</code> ({item['n']:,} files)")
            body += f"<tr><td>{name}</td><td>{esc(item['text'])}</td></tr>\n"
        html += f"<h3>{esc(group)}</h3>\n<div class=\"wrap\"><table><tbody>\n{body}</tbody></table></div>\n"
    return Section("files", "Output files", html)
