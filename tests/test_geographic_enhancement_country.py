"""P0-5: country counts read the BOLD portal country/ocean field."""

import numpy as np
import pandas as pd

from boldgenotyper import geographic_enhancement as ge


def _df(col):
    return pd.DataFrame({"haplotype_id": ["H0001"] * 3, "species": ["Genus alpha"] * 3,
                         "lat": [1.0, np.nan, np.nan], "lon": [2.0, np.nan, np.nan],
                         "ocean_basin": ["North Atlantic Ocean", np.nan, np.nan],
                         col: ["Brazil", "", "Chile"]})


def test_country_column_prefers_portal_field():
    assert ge._country_column(_df("country/ocean")) == "country/ocean"
    assert ge._country_column(_df("country")) == "country"
    both = _df("country/ocean").assign(country=["x", "y", "z"])
    assert ge._country_column(both) == "country/ocean"
    assert ge._country_column(pd.DataFrame({"x": [1]})) is None


def test_coverage_counts_portal_country():
    cov = ge.assess_geographic_coverage(_df("country/ocean"), group_col="haplotype_id",
                                        geo_category="ocean_basin")
    assert cov.loc[0, "with_country"] == 2          # empty string is missing
    assert round(cov.loc[0, "pct_with_country"], 1) == 66.7


def test_coverage_counts_legacy_country():
    cov = ge.assess_geographic_coverage(_df("country"), group_col="haplotype_id",
                                        geo_category="ocean_basin")
    assert cov.loc[0, "with_country"] == 2
