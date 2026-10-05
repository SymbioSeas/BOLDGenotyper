"""Tree rooting: midpoint by default, outgroups, and support values kept on their splits."""

from io import StringIO

import pytest
from Bio import Phylo

from boldgenotyper import phylogenetics as ph

# FastTree-style: trifurcating root, a distinct support on every internal split
NWK = "((A:0.1,B:0.2)0.91:0.05,((C:0.1,D:0.1)0.82:0.3,(E:0.2,F:0.9)0.73:0.4)0.64:0.2,G:0.15);"


def _read(text=NWK):
    return Phylo.read(StringIO(text), "newick")


def _splits(tree):
    """Support of each non-trivial split, keyed by the side without tip 'A'."""
    tips = frozenset(t.name for t in tree.get_terminals())
    out = {}
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        below = frozenset(t.name for t in clade.get_terminals())
        side = below if "A" not in below else tips - below
        if 1 < len(side) < len(tips) - 1:
            assert out.get(side, clade.confidence) == clade.confidence
            out[side] = clade.confidence
    return out


def test_biopython_alone_moves_supports():
    """Why rerooting goes through _reroot_keeping_support (Biopython 1.88)."""
    tree = _read()
    tree.root_at_midpoint()
    assert _splits(tree) != _splits(_read())


def test_midpoint_root_keeps_supports_on_their_splits():
    before = _splits(_read())
    tree = ph.midpoint_root_tree(_read())
    assert _splits(tree) == before
    assert len(tree.root.clades) == 2 and tree.root.confidence is None
    assert tree.total_branch_length() == pytest.approx(_read().total_branch_length())


def test_midpoint_is_halfway_along_the_longest_path():
    tree = ph.midpoint_root_tree(_read())
    depths = {t.name: d for t, d in tree.depths().items() if t.is_terminal()}
    # longest path: F to B (0.9 + 0.4 + 0.2 + 0.05 + 0.2 = 1.75)
    assert depths["F"] == pytest.approx(0.875) and depths["B"] == pytest.approx(0.875)


def test_outgroup_spanning_the_written_root():
    # {A, B, C, D, E, F} vs {G} is a split, but the written root lies inside it
    tree, n_extra = ph.root_on_tips(_read(), ["A", "B", "C", "D", "E", "F"])
    sides = sorted(sorted(t.name for t in c.get_terminals()) for c in tree.root.clades)
    assert sides == [["A", "B", "C", "D", "E", "F"], ["G"]] and n_extra == 0
    assert _splits(tree) == _splits(_read())


def test_outgroup_clade_and_non_monophyletic_outgroup():
    tree, n_extra = ph.root_on_tips(_read(), ["E", "F"])
    assert sorted(sorted(t.name for t in c.get_terminals()) for c in tree.root.clades)[1] == ["E", "F"]
    assert n_extra == 0 and _splits(tree) == _splits(_read())
    assert tree.total_branch_length() == pytest.approx(_read().total_branch_length())
    assert len(tree.root.clades) == 2
    tree, n_extra = ph.root_on_tips(_read(), ["A", "C"])        # not a clade
    assert n_extra > 0
    with pytest.raises(ValueError):
        ph.root_on_tips(_read(), ["X"])


def test_file_level_rooting(tmp_path, caplog):
    f = tmp_path / "t.nwk"
    f.write_text(NWK)
    assert ph.midpoint_root_tree_file(str(f)) == str(f)
    assert _splits(Phylo.read(f, "newick")) == _splits(_read())
    f.write_text(NWK)
    assert ph.reroot_tree_by_label(str(f), "G") == str(f)
    assert ["G"] in [[t.name for t in c.get_terminals()] for c in Phylo.read(f, "newick").root.clades]
    f.write_text(NWK)
    assert ph.reroot_tree_by_tips(str(f), ["A", "C"], description="outgroup FASTA") == str(f)
    assert any("not a clade" in r.getMessage() for r in caplog.records)
    assert ph.reroot_tree_by_label(str(f), "missing") is None


def test_rooted_and_relabeled_files_keep_three_decimal_supports(tmp_path):
    """FastTree writes supports with three decimals; Biopython's default would round 0.697 to 0.70,
    across the tree evidence threshold."""
    f = tmp_path / "t.nwk"
    f.write_text(NWK.replace("0.73", "0.697"))
    ph.midpoint_root_tree_file(str(f))
    assert 0.697 in _splits(Phylo.read(f, "newick")).values()
    f.write_text(NWK.replace("0.73", "0.697"))
    ph.reroot_tree_by_tips(str(f), ["G"])
    assert 0.697 in _splits(Phylo.read(f, "newick")).values()
    aln, taxonomy, out = tmp_path / "aln.fasta", tmp_path / "taxonomy.csv", tmp_path / "labeled.nwk"
    aln.write_text("".join(f">{t}\nACGT\n" for t in "ABCDEFG"))
    taxonomy.write_text("haplotype_id,haplotype_sp\n" + "".join(f"{t},Taxon {t}\n" for t in "ABCDEFG"))
    ph.relabel_tree_and_alignment(str(f), str(aln), str(taxonomy), str(out), str(tmp_path / "labeled.fasta"))
    assert ")0.697:" in out.read_text()
