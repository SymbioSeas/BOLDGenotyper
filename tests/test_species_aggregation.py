"""Species aggregation counts every species-level haplotype, whatever its label share."""

import pandas as pd

from boldgenotyper import species_analysis as sa


def test_species_level_haplotypes_are_aggregated():
    meta = pd.DataFrame({"processid": ["a", "b", "c", "d"], "haplotype_id": ["H1", "H1", "H2", "H3"]})
    comp = pd.DataFrame({"haplotype_id": ["H1", "H2", "H3"],
                         "primary_species": ["Genus alpha", "Genus alpha", "Genus"],
                         "primary_species_pct": [1.0, 0.56, 0.5],     # H2: curated or evidence taxon
                         "is_ambiguous": [False, False, True],
                         "is_multi_species": [False, True, True]})
    _, summary = sa.aggregate_samples_by_species(meta, comp)
    assert summary.set_index("species").loc["Genus alpha", "n_samples"] == 3
    assert "Genus" not in set(summary["species"])
