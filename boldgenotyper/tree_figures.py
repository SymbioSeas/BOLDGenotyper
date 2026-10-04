"""Tree figures (phylogeny/)::

    tree.pdf/.svg            every tip, 0.14 in per tip, 6 pt labels, supports >= 0.70
    tree_overview.pdf/.svg   the full tree for up to 40 haplotypes; above, each taxon
                             whose haplotypes form a clade is collapsed to one tip
                             "Taxon (k haplotypes, n records)" and taxa that do not
                             form a clade are left expanded and listed
    tree_by_genus.pdf        one page per genus with two or more tips: the subtree at
                             the genus's most recent common ancestor (tips of other
                             genera inside it included); by_genus/tree_by_genus/<Genus>.svg

The tree is read as written (phylogeny/tree.nwk is rooted by the pipeline).
"""

import copy
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from . import figure_pages as fp

logger = logging.getLogger(__name__)

LABEL_SIZE = 6


def full_tree_size(n_tips: int) -> Tuple[float, float]:
    """Figure size (inches) for a tree with n_tips tips."""
    return 10, max(4.0, 0.14 * n_tips + 1.0)


def _haplotype_info(haplotypes: pd.DataFrame):
    h = haplotypes.drop_duplicates("haplotype_id").set_index("haplotype_id")
    ids = h.index.astype(str)
    label = (h["haplotype_label"] if "haplotype_label" in h.columns else pd.Series(ids, index=h.index))
    label = label.where(label.notna(), pd.Series(ids, index=h.index))
    taxon = h["assigned_sp"].astype(object).where(h["assigned_sp"].notna(), fp.UNASSIGNED)
    n = pd.to_numeric(h.get("n_members", 0), errors="coerce").fillna(0).astype(int)
    return dict(zip(ids, label.astype(str))), dict(zip(ids, taxon.astype(str))), dict(zip(ids, n))


def collapse_groups(tree, taxon_of: Dict[str, str], records_of: Dict[str, int]):
    """
    Copy of tree with each taxon whose tips form a clade collapsed to one tip.

    Returns (tree, taxa whose tips are not a clade). A one-tip taxon is
    relabeled "Taxon (1 haplotype, n records)".
    """
    t = copy.deepcopy(tree)
    groups: Dict[str, list] = {}
    for tip in t.get_terminals():
        groups.setdefault(taxon_of.get(tip.name, fp.UNASSIGNED), []).append(tip)
    not_clade = []
    for taxon in sorted(groups):
        tips = groups[taxon]
        n_records = sum(records_of.get(x.name, 0) for x in tips)
        if len(tips) == 1:
            tips[0].name = f"{taxon} (1 haplotype, {n_records} records)"
            continue
        names = {x.name for x in tips}
        ancestor = t.common_ancestor(tips)
        if ancestor is not t.root and {x.name for x in ancestor.get_terminals()} == names:
            ancestor.clades = []
            ancestor.name = f"{taxon} ({len(tips)} haplotypes, {n_records} records)"
        else:
            not_clade.append(taxon)
    return t, not_clade


def _draw(tree, title: str, labels: Dict[str, str], support_min: float,
          caption: Optional[str] = None):
    import matplotlib
    import matplotlib.pyplot as plt
    from Bio import Phylo

    n_tips = tree.count_terminals()
    fig, ax = plt.subplots(figsize=full_tree_size(n_tips))
    with matplotlib.rc_context({"font.size": LABEL_SIZE}):
        Phylo.draw(tree, axes=ax, do_show=False,
                   label_func=lambda c: labels.get(c.name, c.name) if c.is_terminal() else None,
                   branch_labels=lambda c: (f"{c.confidence:.2f}" if c.confidence is not None
                                            and c.confidence >= support_min else None))
    ax.set_title(title, fontsize=10)
    ax.set_yticks([])
    ax.set_ylabel("")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    if caption:
        fig.text(0.01, -0.01, caption, ha="left", va="top", fontsize=8, wrap=True)
    return fig


def write_tree_figures(tree_path, haplotypes: pd.DataFrame, out_dir,
                       max_overview_haplotypes: int = 40,
                       support_min: float = 0.70) -> Dict[str, Dict]:
    """Write tree, tree_overview and tree_by_genus; returns manifest entries."""
    from Bio import Phylo

    out_dir = Path(out_dir)
    tree = Phylo.read(str(tree_path), "newick")
    labels, taxon_of, records_of = _haplotype_info(haplotypes)
    support_note = f"Branch values: supports of {support_min:.2f} or more."
    entries: Dict[str, Dict] = {}

    paths = fp.save_figure(_draw(tree, "Haplotype tree", labels, support_min, support_note),
                           out_dir / "tree")
    entries["tree"] = {"overview": [p.name for p in paths]}

    n_tips = tree.count_terminals()
    if fp.overview_level(n_tips, max_overview_haplotypes) == "haplotype":
        overview = _draw(tree, "Haplotype tree", labels, support_min, support_note)
    else:
        collapsed, not_clade = collapse_groups(tree, taxon_of, records_of)
        caption = ("Taxa whose haplotypes form a clade are collapsed to one tip. " + support_note
                   + (f"\nNot monophyletic: {', '.join(not_clade)}" if not_clade else ""))
        overview = _draw(collapsed, f"Haplotype tree, {n_tips} haplotypes collapsed by taxon",
                         labels, support_min, caption)
    entries["tree_overview"] = {"overview": [p.name for p in fp.save_figure(overview, out_dir / "tree_overview")]}

    by_genus: Dict[str, List[str]] = {}
    for tip in tree.get_terminals():
        genus = taxon_of.get(tip.name, fp.UNASSIGNED).split()[0]
        by_genus.setdefault(genus, []).append(tip.name)
    by_genus = {g: tips for g, tips in sorted(by_genus.items()) if len(tips) >= 2}

    def page(genus, tips):
        ancestor = tree.common_ancestor(tips)
        sub = Phylo.BaseTree.Tree(root=copy.deepcopy(ancestor), rooted=True)
        others = sorted({taxon_of.get(t.name, fp.UNASSIGNED) for t in sub.get_terminals()
                         if taxon_of.get(t.name, fp.UNASSIGNED).split()[0] != genus})
        caption = (f"The clade at the most recent common ancestor of {genus}'s haplotypes. "
                   + support_note
                   + (f"\nOther taxa in this clade: {', '.join(others)}" if others else ""))
        return _draw(sub, f"{genus}", labels, support_min, caption)

    entry = fp.write_paged_figure("tree_by_genus", out_dir, None, page, by_genus,
                                  pages_name="tree_by_genus", pages_dir="by_genus")
    if entry:
        entries["tree_by_genus"] = entry
    return entries
