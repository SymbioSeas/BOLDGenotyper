"""Page shell and components of report.html: one self-contained file whose
figures are linked from the run folder (open it from that folder)."""

import json
from dataclasses import dataclass
from html import escape
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

CSS = """
:root { --ink: #1f2328; --muted: #59636e; --line: #d1d9e0; --panel: #f6f8fa; --bg: #ffffff;
        --accent: #0969da; --bar: #4c8bd6; --warn: #9a6700; }
@media (prefers-color-scheme: dark) {
  :root { --ink: #e6edf3; --muted: #9198a1; --line: #3d444d; --panel: #151b23; --bg: #0d1117;
          --accent: #4493f8; --bar: #4493f8; --warn: #d29922; }
  img.fig { background: #ffffff; }
}
* { box-sizing: border-box; }
body { margin: 0; font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: var(--ink);
       background: var(--bg); line-height: 1.5; }
header { padding: 20px 16px 8px; max-width: 1200px; margin: 0 auto; }
header h1 { margin: 0; font-size: 1.6em; }
header p { margin: 4px 0 0; color: var(--muted); }
nav { position: sticky; top: 0; z-index: 2; background: var(--bg); border-bottom: 1px solid var(--line);
      padding: 8px 16px; overflow-x: auto; white-space: nowrap; }
nav a { color: var(--accent); text-decoration: none; margin-right: 16px; font-size: 0.95em; }
main { max-width: 1200px; margin: 0 auto; padding: 0 16px 48px; }
.section { padding-top: 12px; }
.section h2 { border-bottom: 1px solid var(--line); padding-bottom: 4px; margin-top: 28px; }
.section h3 { margin-top: 22px; }
.note { color: var(--muted); font-size: 0.9em; }
.wrap { overflow-x: auto; }
table { border-collapse: collapse; margin: 8px 0; font-size: 0.9em; }
th, td { border: 1px solid var(--line); padding: 4px 8px; text-align: left; vertical-align: top; }
th { background: var(--panel); }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.more { color: var(--muted); font-size: 0.85em; margin: 2px 0 12px; }
.counts { display: flex; flex-wrap: wrap; gap: 10px; margin: 10px 0; }
.count { border: 1px solid var(--line); border-radius: 6px; padding: 8px 12px; min-width: 120px; background: var(--panel); }
.count b { display: block; font-size: 1.3em; }
.count span { color: var(--muted); font-size: 0.85em; }
ul.findings li { margin: 4px 0; }
.bars { margin: 8px 0; max-width: 760px; }
.bar-row { display: grid; grid-template-columns: minmax(120px, 220px) 1fr minmax(70px, auto); gap: 8px;
           align-items: center; margin: 3px 0; font-size: 0.9em; }
.bar-track { background: var(--panel); border-radius: 3px; height: 14px; }
.bar { display: block; background: var(--bar); height: 14px; border-radius: 3px; }
.bar-value { color: var(--muted); font-variant-numeric: tabular-nums; }
figure { margin: 10px 0 18px; }
img.fig { max-width: 100%; height: auto; border: 1px solid var(--line); border-radius: 4px; }
figcaption { font-size: 0.85em; color: var(--muted); }
.switch { margin: 6px 0; }
.switch select { font-size: 0.95em; margin-right: 8px; max-width: 100%; }
.statement { border-left: 4px solid var(--accent); background: var(--panel); padding: 10px 14px;
             font-family: Georgia, "Times New Roman", serif; }
code { font-size: 0.9em; }
"""

SCRIPT = """
document.querySelectorAll('select[data-target]').forEach(function (sel) {
  sel.addEventListener('change', function () {
    var img = document.getElementById(sel.dataset.target);
    img.src = sel.value;
    img.alt = sel.options[sel.selectedIndex].text;
    var link = document.getElementById(sel.dataset.target + '-link');
    if (link) { link.href = sel.value; }
  });
});
document.querySelectorAll('[data-two-level]').forEach(function (box) {
  var data = JSON.parse(box.querySelector('script[type="application/json"]').textContent);
  var group = box.querySelector('select.group');
  var view = box.querySelector('select.view');
  var img = box.querySelector('img');
  function show() { img.src = view.value; img.alt = view.options[view.selectedIndex].text; }
  group.addEventListener('change', function () {
    view.innerHTML = '';
    data[group.options[group.selectedIndex].dataset.group].forEach(function (pair) {
      var opt = document.createElement('option');
      opt.value = pair[1]; opt.textContent = pair[0]; view.appendChild(opt);
    });
    show();
  });
  view.addEventListener('change', show);
});
"""


@dataclass
class Section:
    id: str
    title: str
    html: str


def esc(value) -> str:
    """Escaped text; empty for None and NaN."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return escape(str(value))


def num(value, fmt: str = "{:,}") -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    try:
        return fmt.format(int(value) if float(value).is_integer() and "{:," in fmt else value)
    except (TypeError, ValueError):
        return esc(value)


def table(df: pd.DataFrame, columns: List[str], max_rows: int = 10, file: Optional[str] = None,
          headers: Optional[Dict[str, str]] = None) -> str:
    """The first max_rows rows; 'N more in <file>' with the true count."""
    cols = [c for c in columns if c in df.columns]
    if df.empty or not cols:
        return '<p class="note">None.</p>\n'
    head = "".join(f"<th>{esc((headers or {}).get(c, c))}</th>" for c in cols)
    rows = []
    for _, r in df.head(max_rows).iterrows():
        cells = []
        for c in cols:
            v = r[c]
            numeric = isinstance(v, (int, float)) and not isinstance(v, bool) and not pd.isna(v)
            cells.append(f'<td class="num">{num(v, "{:,}" if float(v).is_integer() else "{:.4g}")}</td>'
                         if numeric else f"<td>{esc(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    html = (f'<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>\n'
            + "\n".join(rows) + "\n</tbody></table></div>\n")
    if len(df) > max_rows:
        where = f" in <code>{esc(file)}</code>" if file else ""
        html += f'<p class="more">{len(df) - max_rows:,} more{where} ({len(df):,} in all).</p>\n'
    return html


def kv_table(rows: Sequence[Tuple[str, str]]) -> str:
    """Two-column table; values are HTML."""
    body = "".join(f"<tr><td><strong>{esc(k)}</strong></td><td>{v}</td></tr>\n" for k, v in rows)
    return f'<div class="wrap"><table><tbody>\n{body}</tbody></table></div>\n'


def counts(items: Sequence[Tuple[str, str]]) -> str:
    return ('<div class="counts">' + "".join(f'<div class="count"><b>{v}</b><span>{esc(k)}</span></div>'
                                             for k, v in items) + "</div>\n")


def bars(items: Sequence[Tuple[str, float, str]], unit: str = "", scale: Optional[float] = None) -> str:
    """CSS bar chart: (label, value, note); widths relative to scale (default: the largest value)."""
    if not items:
        return ""
    top = scale if scale else max((v for _, v, _ in items), default=0) or 1
    rows = []
    for label, value, note in items:
        width = 100.0 * float(value) / top if top else 0.0
        rows.append(f'<div class="bar-row"><span>{esc(label)}</span><span class="bar-track">'
                    f'<span class="bar" style="width: {width:.1f}%"></span></span>'
                    f'<span class="bar-value">{num(value, "{:,}" if float(value).is_integer() else "{:.1f}")}'
                    f'{esc(unit)}{" " + esc(note) if note else ""}</span></div>')
    return '<div class="bars">' + "\n".join(rows) + "</div>\n"


def figure(src: str, alt: str, links: Sequence[Tuple[str, str]] = ()) -> str:
    caption = " | ".join(f'<a href="{esc(href)}">{esc(text)}</a>' for text, href in links)
    return (f'<figure><img class="fig" src="{esc(src)}" alt="{esc(alt)}" loading="lazy">'
            + (f"<figcaption>{caption}</figcaption>" if caption else "") + "</figure>\n")


def switcher(switch_id: str, options: Sequence[Tuple[str, str]], alt: str,
             links: Sequence[Tuple[str, str]] = ()) -> str:
    """A figure and, with two or more options, the select that changes it."""
    if not options:
        return ""
    first_label, first_src = options[0]
    select = ""
    if len(options) > 1:
        opts = "".join(f'<option value="{esc(src)}">{esc(label)}</option>' for label, src in options)
        select = (f'<div class="switch"><label>Show: <select data-target="{esc(switch_id)}">{opts}'
                  f"</select></label></div>\n")
    caption = " | ".join(f'<a href="{esc(href)}">{esc(text)}</a>' for text, href in links)
    return (select + f'<figure><img class="fig" id="{esc(switch_id)}" src="{esc(first_src)}" '
            f'alt="{esc(alt)}: {esc(first_label)}" loading="lazy">'
            + (f"<figcaption>{caption}</figcaption>" if caption else "") + "</figure>\n")


def two_level_switcher(switch_id: str, groups: Dict[str, List[Tuple[str, str]]], alt: str) -> str:
    """Two selects: a group (e.g. a metadata field), then one of its views."""
    groups = {g: v for g, v in groups.items() if v}
    if not groups:
        return ""
    first = next(iter(groups))
    group_opts = "".join(f'<option data-group="{esc(g)}">{esc(g)}</option>' for g in groups)
    view_opts = "".join(f'<option value="{esc(src)}">{esc(label)}</option>' for label, src in groups[first])
    data = json.dumps({g: [list(p) for p in v] for g, v in groups.items()}).replace("</", "<\\/")
    return (f'<div data-two-level id="{esc(switch_id)}">'
            f'<script type="application/json">{data}</script>'
            f'<div class="switch"><select class="group">{group_opts}</select>'
            f'<select class="view">{view_opts}</select></div>'
            f'<figure><img class="fig" src="{esc(groups[first][0][1])}" alt="{esc(alt)}" loading="lazy">'
            "</figure></div>\n")


def manifest_options(run, folder: str, figure_name: str, overview_label: str) -> List[Tuple[str, str]]:
    """(label, path) for a figure's overview SVG and its pages, files that exist only."""
    entry = run.manifests.get(folder, {}).get(figure_name) or {}
    options = []
    svg = next((f for f in entry.get("overview", []) if f.endswith(".svg")), None)
    if svg and (run.run_dir / folder / svg).exists():
        options.append((overview_label, f"{folder}/{svg}"))
    for group, path in (entry.get("pages") or {}).items():
        if (run.run_dir / folder / path).exists():
            options.append((group, f"{folder}/{path}"))
    return options


def manifest_file(run, folder: str, figure_name: str, key: str) -> Optional[str]:
    """Path of a figure's overview PDF ('overview_pdf') or pages PDF ('pages_pdf'), if present."""
    entry = run.manifests.get(folder, {}).get(figure_name) or {}
    name = (next((f for f in entry.get("overview", []) if f.endswith(".pdf")), None)
            if key == "overview_pdf" else entry.get("pages_pdf"))
    return f"{folder}/{name}" if name and (run.run_dir / folder / name).exists() else None


def render_page(title: str, subtitle: str, sections: List[Section]) -> str:
    nav = "".join(f'<a href="#section-{esc(s.id)}">{esc(s.title)}</a>' for s in sections)
    body = "\n".join(f'<div class="section" id="section-{esc(s.id)}">\n<h2>{esc(s.title)}</h2>\n'
                     f"{s.html}</div>" for s in sections)
    return ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            f"<title>{esc(title)}</title>\n<style>{CSS}</style>\n</head>\n<body>\n"
            f"<header><h1>{esc(title)}</h1><p>{subtitle}</p></header>\n"
            f"<nav>{nav}</nav>\n<main>\n{body}\n</main>\n<script>{SCRIPT}</script>\n</body>\n</html>\n")
