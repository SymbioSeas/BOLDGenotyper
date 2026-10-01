"""P0-9: species labels classified and counted under explicit definitions."""

import pandas as pd

from boldgenotyper import taxon_labels as tl


def test_label_levels():
    assert tl.label_level("Genus alpha") == "species"
    assert tl.label_level("Genus alpha beta") == "species"
    assert tl.label_level("Genus sp. 1") == "open"
    assert tl.label_level("Genus cf. alpha") == "open"
    assert tl.label_level("Genus aff. alpha") == "open"
    assert tl.label_level("Genus") == "genus"
    assert tl.label_level("") == "missing"
    assert tl.label_level(None) == "missing"
    assert tl.label_level(float("nan")) == "missing"
    assert tl.genus_of("Genus alpha") == "Genus"
    assert tl.genus_of("") is None


def test_species_count_summary_uses_distinct_definitions():
    input_labels = ["Genus alpha", "Genus beta", "Genus sp. 1", "Genus", "", "Genus gamma"]
    qc_labels = ["Genus alpha", "Genus beta", "Genus"]
    tax = pd.DataFrame({"assigned_sp": ["Genus alpha", "Genus beta", "Genus", "Genus alpha"],
                        "majority_fraction": [1.0, 0.6, 1.0, 0.9]})
    out = tl.species_count_summary(input_labels, qc_labels, tax, majority_threshold=0.7)
    assert out == {
        "n_species_labels_input": 5,
        "n_species_labels_species_level": 3,
        "n_species_labels_qc_passed": 3,
        "n_species_haplotype_assigned": 2,
        "n_species_confident": 1,
    }


# --- Pooled member labels (Q5, Q6) -------------------------------------------

def test_pool_labels_support_and_text():
    p = tl.pool_labels(["Genus alpha"] * 42 + ["Genus beta"] * 4 + ["Genus"] * 3 + [""] * 2)
    assert p.leading_species == "Genus alpha"
    assert p.n_leading == 42 and p.n_species_level == 46
    assert round(p.support, 2) == 0.91
    assert p.n_compatible == 3 and p.n_missing == 2 and p.n_distinct_species == 2
    assert p.text == "Genus alpha 42; Genus beta 4; Genus 3; (no species label) 2"


def test_genus_only_congeners_do_not_conflict():
    p = tl.pool_labels(["Genus alpha"] * 5 + ["Genus"] * 5 + ["Genus sp. 1"])
    assert p.support == 1.0 and p.n_compatible == 6 and p.n_distinct_species == 1


def test_genus_only_label_of_another_genus_counts_against():
    p = tl.pool_labels(["Genus alpha"] * 3 + ["Othergenus"])
    assert p.n_species_level == 4 and p.support == 0.75


def test_tied_species_have_no_leading_species():
    p = tl.pool_labels(["Genus alpha"] * 4 + ["Genus beta"] * 4)
    assert p.leading_species is None and p.support == 0.5


def test_no_species_level_labels():
    p = tl.pool_labels(["Genus", "", "Genus"])
    assert p.leading_species is None and p.support is None
    assert p.text == "Genus 2; (no species label) 1"
