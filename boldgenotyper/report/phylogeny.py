"""Phylogeny: tree overview with a genus switcher, rooting, taxa that are not
a clade, and haplotypes left out of the tree."""

from .data import RunData
from .page import Section, esc, kv_table, manifest_file, manifest_options, num, switcher, table

ROOTING = {
    "midpoint": "at the midpoint of the longest path between two tips",
    "outgroup_label": "on the outgroup tip {outgroup}",
    "outgroup_taxon": "on the clade of haplotypes assigned to {outgroup}",
    "outgroup_fasta": "on the outgroup sequences in {outgroup}",
    "none": "not rooted (as written by FastTree)",
}


def build(run: RunData) -> Section:
    if not run.file("phylogeny/tree.nwk").exists():
        return Section("phylogeny", "Phylogeny",
                       '<p class="note">No tree was built in this run (--build-tree not given).</p>\n')
    applied = run.applied
    rooting = applied.get("tree_rooting") or {}
    rooted = ROOTING.get(rooting.get("method"), "not recorded").format(outgroup=esc(rooting.get("outgroup", "")))
    options = (manifest_options(run, "phylogeny", "tree_overview", "Overview")
               + manifest_options(run, "phylogeny", "tree_by_genus", "Overview"))
    links = [(text, path) for text, path in (
        ("full tree (PDF)", manifest_file(run, "phylogeny", "tree", "overview_pdf")),
        ("full tree (SVG)", "phylogeny/tree.svg" if run.file("phylogeny/tree.svg").exists() else None),
        ("tree by genus (PDF)", manifest_file(run, "phylogeny", "tree_by_genus", "pages_pdf")),
        ("tree.nwk", "phylogeny/tree.nwk"),
        ("tree_labeled.nwk", "phylogeny/tree_labeled.nwk" if run.file("phylogeny/tree_labeled.nwk").exists() else None),
    ) if path]
    html = (kv_table([("Tips", f"{num(applied.get('n_tree_tips'))} haplotypes; {num(applied.get('n_tree_excluded'))} "
                               "left out (below)"),
                      ("Rooting", rooted + "; support values stay on their branches")])
            + "<p>Overview: taxa whose haplotypes form a clade are collapsed when the run has more than "
              "40 haplotypes. Choose a genus to see the clade at its most recent common ancestor.</p>\n"
            + switcher("tree-view", options, "Tree", links)
            + '<p class="note">Open the Newick files in a tree editor such as '
              '<a href="https://treeviewer.org/">TreeViewer</a> (Bianchini &amp; Sánchez-Baracaldo, 2024) '
              "to re-root or restyle.</p>\n")
    not_mono = run.non_monophyletic()
    html += "<h3>Taxa that do not form a clade</h3>\n" + (
        f"<p>{esc(', '.join(not_mono))}</p>\n" if not_mono else '<p class="note">None.</p>\n')
    html += ("<h3>Haplotypes left out of the tree</h3>\n"
             + table(run.excluded_tips, ["haplotype_id", "reason"], file="phylogeny/excluded_tips.txt"))
    return Section("phylogeny", "Phylogeny", html)
