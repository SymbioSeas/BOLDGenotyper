"""A --geo-category named like an input field replaces that field: say so in the log."""

from boldgenotyper import cli


def test_clash_with_bold_field_is_reported():
    message = cli.region_column_clash(["processid", "ecoregion", "realm"], "ecoregion")
    assert "ecoregion" in message
    assert "replaces" in message


def test_no_clash():
    assert cli.region_column_clash(["processid", "ecoregion"], "biogeographic_realm") is None
    assert cli.region_column_clash(["processid", "ecoregion"], "ocean_basin") is None
