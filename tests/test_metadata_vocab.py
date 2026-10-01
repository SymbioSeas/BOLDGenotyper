"""P1-5: metadata vocabularies are normalized before tabulation."""

import numpy as np
import pandas as pd
import pytest

from boldgenotyper import metadata_summary as ms


def _norm(field, values):
    out = ms.normalize_vocabulary(pd.DataFrame({field: values}), [field])[field]
    return [None if (isinstance(v, float) and np.isnan(v)) or v is None else v for v in out]


def test_sex_codes():
    assert _norm("sex", ["M", "F", "U", " ", "", "male", None]) == \
        ["Male", "Female", None, None, None, "Male", None]


def test_life_stage_codes_and_case():
    assert _norm("life_stage", ["A", "I", "adult", "Immature", "Immature - Not Calcified", "M", "J", "Pup"]) == \
        ["Adult", "Immature", "Adult", "Immature", "Immature - Not Calcified", "Mature", "Juvenile", "Pup"]


def test_reproduction_codes():
    assert _norm("reproduction", ["S", "A", "U", " "]) == ["Sexual", "Asexual", None, None]


def test_gis_suffixes_and_underscores():
    assert _norm("biome", ["Mangroves_simplify-0.01_buffered-0.018",
                           "Tropical_&_Subtropical_Moist_Broadleaf_Forest",
                           "Temperate_Grasslands_simplify-0.001_buffered-0.018, Mangroves_simplify-0.01_buffered-0.018"]) == \
        ["Mangroves", "Tropical & Subtropical Moist Broadleaf Forest", "Temperate Grasslands, Mangroves"]
    assert _norm("realm", ["Oceania_simplify-0_buffered-0.018", "Nearctic, Neotropic"]) == \
        ["Oceania", "Nearctic, Neotropic"]


def test_tabulation_has_no_empty_category():
    df = pd.DataFrame({"haplotype_sp": ["H1"] * 5, "reproduction": ["S", " ", "", "S", None]})
    df = ms.normalize_vocabulary(df, ["reproduction"])
    tab = ms.analyze_categorical_field(df, "reproduction", "haplotype_sp")
    assert set(tab["value"]) == {"Sexual"}
