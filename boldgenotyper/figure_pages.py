"""Helpers for figures that scale to large runs.

A figure panel colors at most a set number of items (12 by default); the rest
are drawn in gray as "Other (n)". Large figures are split into an overview and
one page per species-level (or genus-level) group: a multi-page PDF for
publication and one SVG per page for the report. Each figure folder holds a
manifest.json naming every figure's overview and pages.
"""

import json
import re
from pathlib import Path
from typing import Callable, Dict, List, Optional

import pandas as pd

OTHER_COLOR = "#9e9e9e"
UNASSIGNED = "Unassigned"


def slug(text) -> str:
    """File-name form of a group or figure name: 'Genus species' -> 'Genus_species'."""
    return re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_") or "page"


def cap_items(counts: pd.Series, max_items: int = 12, noun: str = "") -> pd.Series:
    """
    Item -> display label. With more than max_items items, the max_items - 1
    with the largest counts keep their names (ties by name) and the rest
    become 'Other (k noun)'.
    """
    labels = pd.Series(counts.index.astype(str), index=counts.index)
    if len(counts) <= max_items:
        return labels
    order = sorted(counts.index, key=lambda i: (-counts[i], str(i)))
    folded = order[max_items - 1:]
    other = f"Other ({len(folded)}{' ' + noun if noun else ''})"
    labels[folded] = other
    return labels


def article(word) -> str:
    """'an' before a vowel, else 'a' (for captions naming a column)."""
    return "an" if str(word)[:1].lower() in "aeiou" else "a"


def count(n: int, word: str, plural: str = "") -> str:
    """'1 record', '1,234 records' (plural defaults to word + 's')."""
    return f"{n:,} {word if n == 1 else (plural or word + 's')}"


def add_caption(ax, text: str) -> None:
    """A small caption between the title and the panel (a fixed 4 points above
    the panel, so tall figures keep it below the title)."""
    ax.annotate(text, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 4),
                textcoords="offset points", fontsize=7, color="#444444", ha="left", va="bottom")


def is_other(label) -> bool:
    return str(label).startswith("Other (")


def item_colors(labels: List[str]) -> Dict[str, str]:
    """Colors for labels in order; 'Other (...)' labels are gray."""
    from .visualization import get_genotype_colors

    named = [lab for lab in dict.fromkeys(labels) if not is_other(lab)]
    colors = dict(zip(named, get_genotype_colors(len(named))))
    return {lab: (OTHER_COLOR if is_other(lab) else colors[lab]) for lab in dict.fromkeys(labels)}


def species_groups(haplotypes: pd.DataFrame) -> Dict[str, List[str]]:
    """
    Group -> haplotype IDs. A group is a haplotype's assigned taxon (genus-level
    haplotypes under their genus); haplotypes by members, most first; groups
    sorted by name, haplotypes without a taxon last under 'Unassigned'.
    """
    h = haplotypes.copy()
    taxon = h["assigned_sp"].astype(object).where(h["assigned_sp"].notna(), "")
    h["_group"] = taxon.astype(str).str.strip().replace("", UNASSIGNED)
    h["_n"] = pd.to_numeric(h.get("n_members", 0), errors="coerce").fillna(0)
    groups = {}
    for name in sorted(h["_group"].unique(), key=lambda g: (g == UNASSIGNED, g)):
        g = h[h["_group"] == name].sort_values(["_n", "haplotype_id"], ascending=[False, True])
        groups[name] = g["haplotype_id"].astype(str).tolist()
    return groups


def haplotype_colors(haplotypes: pd.DataFrame, max_items: int = 12) -> Dict[str, str]:
    """Haplotype ID -> color, assigned within each group, the same on every page."""
    n = pd.to_numeric(haplotypes.set_index("haplotype_id").get("n_members"), errors="coerce").fillna(0)
    colors = {}
    for ids in species_groups(haplotypes).values():
        labels = cap_items(n.reindex(ids).fillna(0), max_items)
        palette = item_colors(list(labels.values))
        colors.update({h: palette[labels[h]] for h in ids})
    return colors


def overview_level(n_haplotypes: int, max_overview_haplotypes: int = 40) -> str:
    """'haplotype' when an overview can show each haplotype, else 'species'."""
    return "haplotype" if n_haplotypes <= max_overview_haplotypes else "species"


def _savefig(fig, path: Path) -> None:
    import matplotlib

    with matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
        fig.savefig(path, bbox_inches="tight")


def save_figure(fig, base: Path) -> List[Path]:
    """Write base.pdf and base.svg (text kept as text) and close the figure."""
    import matplotlib.pyplot as plt

    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    paths = [base.with_suffix(".pdf"), base.with_suffix(".svg")]
    try:
        for p in paths:
            _savefig(fig, p)
    finally:
        plt.close(fig)
    return paths


def write_paged_figure(name: str, out_dir, overview: Callable, page: Callable,
                       groups: Dict[str, List[str]], pages_name: Optional[str] = None,
                       pages_dir: str = "by_species") -> Dict:
    """
    Draw an overview and one page per group.

    overview() and page(group, haplotype_ids) return a matplotlib figure or
    None (nothing to draw). Writes <name>.pdf/.svg, <pages_name>.pdf (all
    pages; default <name>_by_species) and <pages_dir>/<pages_name or name>/
    <group>.svg. Returns the manifest entry, paths relative to out_dir.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    entry: Dict = {}
    fig = overview() if overview is not None else None
    if fig is not None:
        entry["overview"] = [p.name for p in save_figure(fig, out_dir / name)]

    svg_dir = out_dir / pages_dir / (pages_name or name)
    pages_name = pages_name or f"{name}_by_species"
    pages: Dict[str, str] = {}
    pdf_path = out_dir / f"{pages_name}.pdf"
    pdf = None
    # Fonts are embedded when the PDF is closed: the font type must be set
    # for the whole life of the file, not only while each page is added
    with matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
        try:
            for group, ids in groups.items():
                fig = page(group, ids)
                if fig is None:
                    continue
                if pdf is None:
                    pdf = PdfPages(pdf_path)
                try:
                    pdf.savefig(fig, bbox_inches="tight")
                    svg_dir.mkdir(parents=True, exist_ok=True)
                    svg = svg_dir / f"{slug(group)}.svg"
                    fig.savefig(svg, bbox_inches="tight")
                    pages[group] = svg.relative_to(out_dir).as_posix()
                finally:
                    plt.close(fig)
        finally:
            if pdf is not None:
                pdf.close()
    if pages:
        entry["pages_pdf"] = pdf_path.name
        entry["pages"] = pages
    return entry


def write_manifest(out_dir, entries: Dict[str, Dict]) -> Path:
    """Write (or extend) out_dir/manifest.json: figure -> overview and pages."""
    path = Path(out_dir) / "manifest.json"
    manifest = json.loads(path.read_text()) if path.exists() else {}
    manifest.update({k: v for k, v in entries.items() if v})
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


# ---------------------------------------------------------------------------
# Drawing primitives shared by run_figures and metadata_figures
# ---------------------------------------------------------------------------

def stacked_barh(table: pd.DataFrame, colors: Dict[str, str], title: str, caption: str,
                 xlabel: str, ylabel: str = "", relative: bool = False, sort_rows: bool = True):
    """Horizontal bars, one per row of table, stacked by its columns (counts)."""
    import matplotlib.pyplot as plt
    import numpy as np

    if sort_rows:
        table = table.loc[table.sum(axis=1).sort_values().index]
    values = table.div(table.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) * 100 if relative else table
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.3 * len(table) + 1.6)))
    left = np.zeros(len(table))
    for item in table.columns:
        ax.barh(table.index.astype(str), values[item].values, left=left, color=colors[item],
                edgecolor="white", linewidth=0.4, label=f"{item} ({int(table[item].sum())})")
        left += values[item].values
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="y", labelsize=7 if len(table) > 20 else 9)
    if relative:
        ax.set_xlim(0, 100)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), ncol=1 if table.shape[1] <= 16 else 2,
              frameon=False, fontsize=7, title="records in parentheses", title_fontsize=7)
    ax.set_title(title, fontsize=11, pad=16)
    add_caption(ax, caption)
    return fig


def heatmap(matrix: pd.DataFrame, title: str, caption: str, value_label: str,
            fmt: str = "{:.0f}", vmax: Optional[float] = None):
    """Rows x columns heatmap; cells annotated when the matrix is small."""
    import matplotlib.pyplot as plt

    n_rows, n_cols = matrix.shape
    fig, ax = plt.subplots(figsize=(max(6.0, 0.35 * n_cols + 3.5), max(2.5, 0.28 * n_rows + 1.8)))
    im = ax.imshow(matrix.values.astype(float), aspect="auto", cmap="Blues", vmin=0, vmax=vmax)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([str(c) for c in matrix.columns], rotation=90, fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([str(r) for r in matrix.index], fontsize=7)
    if n_rows * n_cols <= 600:
        top = vmax if vmax is not None else (matrix.values.max() if matrix.size else 0)
        for i in range(n_rows):
            for j in range(n_cols):
                v = matrix.iat[i, j]
                if v:
                    ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=6,
                            color="white" if top and v > 0.6 * top else "black")
    fig.colorbar(im, ax=ax, label=value_label, shrink=0.8)
    ax.set_title(title, fontsize=11, pad=16)
    add_caption(ax, caption)
    return fig
