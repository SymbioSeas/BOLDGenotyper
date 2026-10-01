"""Species-level divergence (--extended) finds haplotypes by their display labels.

The pipeline's divergence matrix is indexed by display labels ("Genus alpha
H0001") while haplotype_taxonomy.csv is keyed by haplotype_id; v1 matched
neither and wrote no species summary or per-species matrices.
"""

import pandas as pd
import pytest

from boldgenotyper import divergence_analysis as da
from boldgenotyper import species_analysis as sa

LABELS = ["Genus alpha H0001", "Genus alpha H0002", "Genus beta H0003"]


@pytest.fixture
def inputs(tmp_path):
    matrix = pd.DataFrame([[0, 0.01, 0.10], [0.01, 0, 0.12], [0.10, 0.12, 0]],
                          index=LABELS, columns=LABELS)
    tax = tmp_path / "haplotype_taxonomy.csv"
    pd.DataFrame({"haplotype_id": ["H0001", "H0002", "H0003"],
                  "assigned_sp": ["Genus alpha", "Genus alpha", "Genus beta"],
                  "haplotype_sp": LABELS}).to_csv(tax, index=False)
    summary = tmp_path / "species_summary.csv"
    pd.DataFrame({"species": ["Genus alpha", "Genus beta"], "n_samples": [5, 2],
                  "n_haplotypes": [2, 1]}).to_csv(summary, index=False)
    return matrix, summary, tax


def test_species_divergence_summary_uses_display_labels(tmp_path, inputs):
    matrix, summary, tax = inputs
    da.generate_species_divergence_analysis(matrix, summary, tax, tmp_path / "out")
    out = pd.read_csv(tmp_path / "out" / "species_divergence_summary.csv").set_index("species")
    assert list(out.index) == ["Genus alpha"]
    assert out.loc["Genus alpha", "max_within_species"] == pytest.approx(0.01)
    assert out.loc["Genus alpha", "min_between_species"] == pytest.approx(0.10)


def test_within_species_matrices_use_display_labels(tmp_path, inputs):
    matrix, summary, tax = inputs
    res = sa.generate_within_species_divergence_matrices(matrix, summary, tax, tmp_path / "out")
    assert list(res["matrices"]) == ["Genus alpha"]
    m = pd.read_csv(res["matrices"]["Genus alpha"], index_col=0)
    assert list(m.index) == LABELS[:2]
