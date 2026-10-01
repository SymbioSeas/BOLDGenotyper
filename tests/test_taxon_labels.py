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
