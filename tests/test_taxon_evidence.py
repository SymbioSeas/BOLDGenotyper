"""Evidence for a haplotype's taxon (never applied)."""

import numpy as np
import pandas as pd

from boldgenotyper import taxon_evidence as te

# H1, H2: Genus alpha (spread 1 site); H3: tied labels, 1 site from H1; H4: Genus beta, far
SEQS = {"H1": "A" * 100, "H2": "A" * 99 + "C", "H3": "C" + "A" * 99, "H4": "G" * 20 + "A" * 80}
TAX = pd.DataFrame({"haplotype_id": ["H1", "H2", "H3", "H4"],
                    "assigned_sp": ["Genus alpha", "Genus alpha", "Genus", "Genus beta"],
                    "assignment_level": ["species", "species", "genus", "species"],
                    "assignment_source": ["labels"] * 4})


def _records(h3_beta_bin=""):
    rows = []
    for i in range(5):
        rows.append(("ALPHA%03d-10" % i, "H1", "Genus alpha", "BOLD:AAA0001", "Inst A"))
    rows += [("PRJA%03d-07" % i, "H3", "Genus beta", h3_beta_bin, "Inst X") for i in range(4)]
    rows += [("PRJ%s%03d-12" % (p, i), "H3", "Genus alpha", "BOLD:AAA0001", "Inst " + p)
             for i, p in enumerate("BCDE")]
    rows += [("BETA%03d-11" % i, "H4", "Genus beta", "BOLD:BBB0002", "Inst B") for i in range(3)]
    df = pd.DataFrame(rows, columns=["processid", "haplotype_id", "label", "bin_uri", "inst"])
    df["status"], df["qc_pass"], df["insdc_acs"] = "assigned", True, ""
    return df


def test_label_source():
    assert te.label_source("ESHKD057-07", "Save Our Seas") == "ESHKD"
    assert te.label_source("GBMNA1234-19", "Mined from GenBank, NCBI", "MN712345") == "GenBank:MN712"


def test_tied_haplotype_gets_a_suggestion():
    ev = te.compute_evidence(TAX, _records(), SEQS).set_index("haplotype_id")
    assert list(ev.columns) == te.EVIDENCE_COLUMNS
    h3 = ev.loc["H3"]
    assert h3["nearest_species_haplotype"] == "H1" and h3["nearest_species_taxon"] == "Genus alpha"
    assert h3["nearest_species_distance"] == 0.01 and h3["taxon_max_within_distance"] == 0.01
    assert h3["within_taxon_spread"] == True
    assert h3["member_bins"] == "BOLD:AAA0001 4; none 4"
    assert h3["bin_taxon"] == "Genus alpha (1.00; 9 of 9 QC-passed records)"
    assert "Genus beta: 4 records, 1 project, 1 institution" in h3["label_sources"]
    assert "Genus alpha: 4 records, 4 projects, 4 institutions" in h3["label_sources"]
    assert h3["suggested_sp"] == "Genus alpha" and h3["suggestion_basis"] == "distance, BIN agree"
    assert ev.loc["H1", "suggested_sp"] == "" and ev.loc["H1", "suggestion_basis"] == ""


def test_tree_sister(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text("((H3:0.01,(H1:0.01,H2:0.01)0.9:0.01)0.95:0.02,H4:0.2);")
    ev = te.compute_evidence(TAX, _records(), SEQS, tree).set_index("haplotype_id")
    assert ev.loc["H3", "tree_sister_taxon"] == "Genus alpha"
    assert ev.loc["H3", "tree_sister_support"] == 0.95
    assert ev.loc["H3", "suggestion_basis"] == "distance, BIN, tree agree"


def test_mixed_sister_clade_names_no_taxon(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text("((H3:0.01,H4:0.01)0.8:0.01,(H1:0.01,H2:0.01)0.9:0.01);")
    tax = TAX.assign(assigned_sp=["Genus alpha", "Genus beta", "Genus", "Genus beta"])
    ev = te.compute_evidence(tax, _records(), SEQS, tree).set_index("haplotype_id")
    assert ev.loc["H3", "tree_sister_taxon"] == "Genus beta"      # first clade: H3 + H4
    tree.write_text("((H3:0.01,(H1:0.01,H4:0.01)0.9:0.01)0.95:0.02,H2:0.2);")
    tax = TAX.assign(assigned_sp=["Genus alpha", "Genus alpha", "Genus", "Genus gamma"])
    ev = te.compute_evidence(tax, _records(), SEQS, tree).set_index("haplotype_id")
    assert ev.loc["H3", "tree_sister_taxon"] == ""                 # H1 and H4 differ
    assert np.isnan(ev.loc["H3", "tree_sister_support"])
    ev = te.compute_evidence(TAX, _records(), SEQS, tmp_path / "missing.nwk").set_index("haplotype_id")
    assert ev.loc["H3", "tree_sister_taxon"] == ""


def test_lines_that_disagree_give_no_suggestion():
    ev = te.compute_evidence(TAX, _records(h3_beta_bin="BOLD:BBB0002"), SEQS).set_index("haplotype_id")
    h3 = ev.loc["H3"]
    assert h3["member_bins"] == "BOLD:AAA0001 4; BOLD:BBB0002 4"
    assert h3["suggested_sp"] == ""


def test_one_line_is_insufficient():
    recs = _records()
    recs["bin_uri"] = ""
    ev = te.compute_evidence(TAX, recs, SEQS).set_index("haplotype_id")
    h3 = ev.loc["H3"]
    assert h3["member_bins"] == "none 8" and h3["bin_taxon"] == ""
    assert h3["suggested_sp"] == ""
    assert h3["suggestion_basis"] == "insufficient evidence (distance Genus alpha only)"


def test_bin_and_distance_disagree():
    recs = _records()
    recs.loc[recs["haplotype_id"] == "H3", "bin_uri"] = "BOLD:BBB0002"
    recs.loc[recs["haplotype_id"] == "H4", "bin_uri"] = "BOLD:BBB0002"
    recs.loc[recs["haplotype_id"] == "H3", "label"] = "Genus beta"
    ev = te.compute_evidence(TAX, recs, SEQS).set_index("haplotype_id")
    h3 = ev.loc["H3"]
    assert h3["suggested_sp"] == ""
    assert h3["suggestion_basis"] == "evidence disagrees: distance Genus alpha; BIN Genus beta"


def test_curated_haplotypes_get_no_suggestion_and_missing_columns_are_empty():
    tax = TAX.assign(assignment_source=["labels", "labels", "curated", "labels"])
    recs = _records()[["processid", "haplotype_id", "label"]]
    ev = te.compute_evidence(tax, recs, SEQS).set_index("haplotype_id")
    assert ev.loc["H3", "suggested_sp"] == "" and ev.loc["H3", "suggestion_basis"] == ""
    assert ev.loc["H3", "member_bins"] == "none 8"
    # H1's nearest species-level haplotype is H2; Genus alpha has no other haplotype
    assert ev.loc["H1", "nearest_species_haplotype"] == "H2"
    assert np.isnan(ev.loc["H1", "taxon_max_within_distance"])
    assert pd.isna(ev.loc["H1", "within_taxon_spread"])


def _no_bins():
    recs = _records()
    recs["bin_uri"] = ""
    return recs


def test_tree_line_collapses_unsupported_clades(tmp_path):
    tree = tmp_path / "t.nwk"
    # H3's nearest clade has support 0.00 (or none); the next holds only Genus alpha
    for newick in ("(((H3:0.01,H1:0.01)0.0:0.0,H2:0.01)0.85:0.02,H4:0.2);",
                   "(((H3:0.01,H1:0.01):0.0,H2:0.01)0.85:0.02,H4:0.2);"):
        tree.write_text(newick)
        h3 = te.compute_evidence(TAX, _no_bins(), SEQS, tree).set_index("haplotype_id").loc["H3"]
        assert h3["tree_sister_taxon"] == "Genus alpha" and h3["tree_sister_support"] == 0.85
        assert h3["suggestion_basis"] == "distance, tree agree"


def test_tree_line_names_no_taxon_without_a_supported_single_taxon_clade(tmp_path):
    tree = tmp_path / "t.nwk"
    for newick in ("(((H3:0.01,H1:0.01)0.0:0.0,H4:0.01)0.9:0.02,H2:0.2);",   # supported clade mixes taxa
                   "((H3:0.01,H1:0.01)0.3:0.0,(H2:0.01,H4:0.2)0.4:0.1);"):   # only the root holds one
        tree.write_text(newick)
        h3 = te.compute_evidence(TAX, _no_bins(), SEQS, tree).set_index("haplotype_id").loc["H3"]
        assert h3["tree_sister_taxon"] == "" and np.isnan(h3["tree_sister_support"])
        assert h3["suggestion_basis"] == ("insufficient evidence (distance Genus alpha only; "
                                          "tree: no supported clade names one taxon)")
    tree.write_text("((H1:0.01,H2:0.01)0.9:0.01,H4:0.2);")
    h3 = te.compute_evidence(TAX, _no_bins(), SEQS, tree).set_index("haplotype_id").loc["H3"]
    assert h3["suggestion_basis"] == "insufficient evidence (distance Genus alpha only; tree: not in the tree)"


def test_weakly_supported_tree_sister_does_not_count(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text("((H3:0.01,(H1:0.01,H2:0.01)0.9:0.01)0.5:0.02,H4:0.2);")
    h3 = te.compute_evidence(TAX, _no_bins(), SEQS, tree).set_index("haplotype_id").loc["H3"]
    assert h3["tree_sister_taxon"] == "" and h3["suggested_sp"] == ""
    assert h3["suggestion_basis"] == ("insufficient evidence (distance Genus alpha only; "
                                      "tree: no supported clade names one taxon)")
    tree.write_text("((H3:0.01,(H1:0.01,H2:0.01)0.9:0.01)0.7:0.02,H4:0.2);")
    ev = te.compute_evidence(TAX, _no_bins(), SEQS, tree)
    assert ev.set_index("haplotype_id").loc["H3", "suggestion_basis"] == "distance, tree agree"
