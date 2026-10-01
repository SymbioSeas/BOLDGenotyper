"""P0-4 regression tests on current BOLD portal coordinate formats.

Fixture: tests/data/bold_portal_coords_2026-09-30.tsv holds 42 rows from a BOLD
portal TSV export (Carcharhiniformes) downloaded 2026-09-30 (public records),
trimmed to the coordinate columns. It covers each "[lat, lon]" shape in that
export (integer and decimal values, all sign combinations), country- and
province-centroid sources, other coord_source values, and empty coordinates.
"""

from pathlib import Path

import pandas as pd
import pytest

from boldgenotyper import metadata

FIXTURE = Path(__file__).parent / "data" / "bold_portal_coords_2026-09-30.tsv"


@pytest.fixture
def portal():
    return pd.read_csv(FIXTURE, sep="\t", comment="#", dtype=str, keep_default_na=False)


def test_every_nonempty_coord_is_parsed(portal):
    df = metadata.parse_coord_column(portal)
    has = portal["coord"].str.strip() != ""
    assert has.sum() == 34
    assert df.loc[has, ["lat", "lon"]].notna().all(axis=None)
    assert df.loc[~has, ["lat", "lon"]].isna().all(axis=None)


@pytest.mark.parametrize("text,expected", [
    ("[23.5863, 81.173004]", (23.5863, 81.173004)),
    ("[-32.933, 151.75]", (-32.933, 151.75)),
    ("[34.945, -75.83]", (34.945, -75.83)),
    ("[-11.524631, -54.355206]", (-11.524631, -54.355206)),
    ("[-27, 87]", (-27.0, 87.0)),
    ("[37, -8.94]", (37.0, -8.94)),
    ("[-34, 151.183]", (-34.0, 151.183)),
    ("[-35.5, 21]", (-35.5, 21.0)),
    ("[0.9, -52]", (0.9, -52.0)),
    ("[-2.43, -44]", (-2.43, -44.0)),
    ("[-9, -34.59]", (-9.0, -34.59)),
    ("[31, -79]", (31.0, -79.0)),
    ("[14, 139.1]", (14.0, 139.1)),
    ("-33.8, 151.283", (-33.8, 151.283)),   # older exports, no brackets
])
def test_coord_shapes(text, expected):
    df = metadata.parse_coord_column(pd.DataFrame({"coord": [text]}))
    assert (df.loc[0, "lat"], df.loc[0, "lon"]) == pytest.approx(expected)


def test_existing_lat_lon_columns_win():
    df = pd.DataFrame({"coord": ["[1, 2]", ""], "lat": ["10.5", ""], "lon": ["20.5", "30"]})
    out = metadata.parse_coord_column(df)
    assert out.loc[0, ["lat", "lon"]].tolist() == [10.5, 20.5]
    assert pd.isna(out.loc[1, "lat"]) and out.loc[1, "lon"] == 30.0


def test_quality_status_and_flag():
    df = pd.DataFrame({
        "coord": ["[10, 20]", "", "[-20.1, 134.2]", "[0, 0]", "[95, 10]", "[1, 1]", "garbage"],
        "coord_source": ["GPS", "", "Coordinates from country centroid", "", "",
                         "Coordinates from province centroid", ""],
    })
    out = metadata.mark_coordinate_quality(metadata.parse_coord_column(df))
    assert list(out["coord_status"]) == ["ok", "missing", "centroid", "zero",
                                         "out_of_range", "centroid", "unparseable"]
    assert list(out["is_geographic_quality"]) == [True, False, False, False, False, False, False]


def test_records_without_coordinates_are_not_geographic_quality(portal):
    """v1 marked records with no coordinates as geographic-quality."""
    out = metadata.mark_coordinate_quality(metadata.parse_coord_column(portal))
    empty = portal["coord"].str.strip() == ""
    assert not out.loc[empty, "is_geographic_quality"].any()
    assert (out.loc[empty, "coord_status"] == "missing").all()


def test_centroid_rows_keep_their_coordinates(portal):
    out = metadata.mark_coordinate_quality(metadata.parse_coord_column(portal))
    cent = out["coord_status"] == "centroid"
    assert cent.sum() == 5
    assert out.loc[cent, ["lat", "lon"]].notna().all(axis=None)
    assert not out.loc[cent, "is_geographic_quality"].any()
