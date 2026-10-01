"""P1-9, P1-10: coverage percentages without grades; one label-based diversity metric."""

import numpy as np
import pandas as pd
import pytest

from boldgenotyper import geographic_enhancement as ge
from boldgenotyper import species_analysis as sa


def _geo_records():
    return pd.DataFrame({
        "haplotype_id": ["H0001"] * 4 + ["H0002"] * 2,
        "species": ["Genus alpha"] * 4 + ["Genus beta"] * 2,
        "lat": [10.0, 11.0, np.nan, np.nan, 5.0, np.nan],
        "lon": [-50.0, -51.0, np.nan, np.nan, -40.0, np.nan],
        "ocean_basin": ["North Atlantic Ocean", np.nan, np.nan, np.nan, np.nan, np.nan],
        "country/ocean": ["Atlantis", "Atlantis", "Atlantis", np.nan, np.nan, np.nan],
    })


def test_coverage_has_percentages_and_no_grade():
    cov = ge.assess_geographic_coverage(_geo_records(), group_col="haplotype_id").set_index("genotype")
    assert "representativeness" not in cov.columns
    assert cov.loc["H0001", "pct_with_coords"] == 50.0
    assert cov.loc["H0001", "pct_with_region"] == 25.0
    assert cov.loc["H0002", "pct_with_coords"] == 50.0


def test_console_and_readme_report_percentages_only(tmp_path, capsys):
    ge.enhance_geographic_analysis(_geo_records(), tmp_path, "Synth", group_col="haplotype_id")
    out = capsys.readouterr().out
    assert "16.7%" in out and "50.0%" in out
    readme = (tmp_path / "README.md").read_text()
    for text in (out, readme):
        assert "representativeness" not in text.lower() and "poor" not in text.lower()
    assert "representativeness" not in (tmp_path / "geographic_coverage.csv").read_text()


def _assignments():
    return pd.DataFrame({
        "processid": [f"P{i}" for i in range(7)],
        "haplotype_id": ["H0001"] * 3 + ["H0002"] + ["H0003"] * 3,
        "primary_species": ["Genus alpha"] * 4 + ["Genus beta"] * 3,
        "is_ambiguous": [False] * 7,
        "ocean_basin": ["North Atlantic Ocean"] * 7,
        "country/ocean": ["Atlantis"] * 7,
        "coord": ["[10.0, -50.0]"] * 7,
    })


def test_species_diversity_has_one_metric_and_a_note(tmp_path):
    div = sa.calculate_species_diversity(_assignments(), output_dir=tmp_path)
    assert list(div.columns) == sa.SPECIES_DIVERSITY_COLUMNS
    assert sa.SPECIES_DIVERSITY_COLUMNS == ["species", "n_samples", "n_haplotypes",
                                           "haplotype_diversity", "note"]
    alpha = div.set_index("species").loc["Genus alpha"]
    # Nei's h = n/(n-1) * (1 - sum p^2) = 4/3 * (1 - (9/16 + 1/16)) = 0.5
    assert alpha.haplotype_diversity == pytest.approx(0.5)
    assert set(div["note"]) == {"label-based, sampling-dependent"}
    assert list(pd.read_csv(tmp_path / "species_diversity.csv").columns) == sa.SPECIES_DIVERSITY_COLUMNS


def test_species_geographic_summary_writes_no_range_boxes(tmp_path):
    sa.generate_species_geographic_summary(_assignments(), output_dir=tmp_path)
    assert (tmp_path / "species_basin_distribution.csv").exists()
    assert not (tmp_path / "species_geographic_range.csv").exists()
