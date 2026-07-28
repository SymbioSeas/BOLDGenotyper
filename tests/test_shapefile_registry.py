from pathlib import Path
import os
import pytest
import yaml as _yaml
from boldgenotyper import shapefile_registry as sr


# ---------------------------------------------------------------------------
# Task 1: path resolution
# ---------------------------------------------------------------------------

def test_resolve_absolute_path_unchanged(tmp_path):
    p = tmp_path / "abs" / "x.shp"
    assert sr._resolve_path(str(p), None, tmp_path) == p


def test_resolve_relative_against_config_dir_when_no_shapefiles_dir(tmp_path):
    got = sr._resolve_path("sub/x.shp", None, tmp_path)
    assert got == (tmp_path / "sub" / "x.shp").resolve()


def test_resolve_relative_against_shapefiles_dir_when_set(tmp_path):
    base = tmp_path / "data"
    got = sr._resolve_path("x.shp", base, tmp_path)
    assert got == (base / "x.shp").resolve()


def test_resolve_expands_user_and_env(tmp_path, monkeypatch):
    monkeypatch.setenv("SHPBASE", str(tmp_path))
    got = sr._resolve_path("$SHPBASE/x.shp", None, Path("/ignored"))
    assert got == (tmp_path / "x.shp").resolve()


# ---------------------------------------------------------------------------
# Task 2: config discovery and loading
# ---------------------------------------------------------------------------

def _write_config(tmp_path, body):
    p = tmp_path / "shapefiles.yaml"
    p.write_text(body)
    return p


def test_load_registry_missing_default_returns_empty(tmp_path, monkeypatch):
    monkeypatch.delenv(sr.ENV_VAR, raising=False)
    monkeypatch.setattr(sr, "default_config_path", lambda: tmp_path / "nope.yaml")
    reg = sr.load_registry()
    assert reg.entries == {}
    assert reg.source is None


def test_load_registry_explicit_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        sr.load_registry(tmp_path / "nope.yaml")


def test_load_registry_parses_entries_and_shapefiles_dir(tmp_path):
    cfg = _write_config(tmp_path, (
        "shapefiles_dir: data\n"
        "shapefiles:\n"
        "  Ecoregions2017:\n"
        "    path: eco/eco.shp\n"
        "    field: REALM\n"
        "    geo_category: biogeographic_realm\n"
    ))
    reg = sr.load_registry(cfg)
    assert "Ecoregions2017" in reg.entries
    assert reg.shapefiles_dir == (tmp_path / "data").resolve()
    assert reg.source == cfg


def test_load_registry_env_var_used(tmp_path, monkeypatch):
    cfg = _write_config(tmp_path, "shapefiles:\n  A:\n    path: a.shp\n")
    monkeypatch.setenv(sr.ENV_VAR, str(cfg))
    reg = sr.load_registry()
    assert "A" in reg.entries


def test_load_registry_rejects_non_mapping(tmp_path):
    cfg = tmp_path / "bad.yaml"
    cfg.write_text("- just\n- a list\n")
    with pytest.raises(ValueError):
        sr.load_registry(cfg)


# ---------------------------------------------------------------------------
# Task 3: name resolution
# ---------------------------------------------------------------------------

def test_resolve_shapefile_generic_defaults(tmp_path):
    cfg = _write_config(tmp_path, (
        "shapefiles:\n"
        "  Eco:\n"
        "    path: eco/eco.shp\n"
    ))
    reg = sr.load_registry(cfg)
    res = sr.resolve_shapefile("Eco", reg)
    assert res.path == (tmp_path / "eco" / "eco.shp").resolve()
    assert res.field is None
    assert res.geo_category == "geographic_region"
    assert res.is_goas is False


def test_resolve_shapefile_goas_routing_and_default_category(tmp_path):
    cfg = _write_config(tmp_path, (
        "shapefiles:\n"
        "  GOaS_v1_20211214:\n"
        "    path: goas/goas_v01.shp\n"
        "    field: name\n"
        "    type: goas\n"
    ))
    reg = sr.load_registry(cfg)
    res = sr.resolve_shapefile("GOaS_v1_20211214", reg)
    assert res.is_goas is True
    assert res.geo_category == "ocean_basin"
    assert res.field == "name"


def test_resolve_shapefile_unknown_name_lists_registered(tmp_path):
    cfg = _write_config(tmp_path, "shapefiles:\n  A:\n    path: a.shp\n  B:\n    path: b.shp\n")
    reg = sr.load_registry(cfg)
    with pytest.raises(ValueError) as e:
        sr.resolve_shapefile("Missing", reg)
    assert "A" in str(e.value) and "B" in str(e.value)


def test_resolve_shapefile_requires_path(tmp_path):
    cfg = _write_config(tmp_path, "shapefiles:\n  Bad:\n    field: x\n")
    reg = sr.load_registry(cfg)
    with pytest.raises(ValueError):
        sr.resolve_shapefile("Bad", reg)


# ---------------------------------------------------------------------------
# Task 4: CLI selection precedence and conflict handling
# ---------------------------------------------------------------------------

def test_resolve_selection_none_when_no_shp():
    assert sr.resolve_selection(
        shp_name=None, custom_shp=None, explicit_field=None,
        explicit_geo_category=None) is None


def test_resolve_selection_conflict_shp_and_custom(tmp_path):
    with pytest.raises(ValueError):
        sr.resolve_selection(
            shp_name="Eco", custom_shp=Path("x.shp"),
            explicit_field=None, explicit_geo_category=None,
            config_path=tmp_path / "shapefiles.yaml")


def test_resolve_selection_explicit_flags_override(tmp_path):
    cfg = _write_config(tmp_path, (
        "shapefiles:\n  Eco:\n    path: eco.shp\n    field: REALM\n"
        "    geo_category: biogeographic_realm\n"
    ))
    res = sr.resolve_selection(
        shp_name="Eco", custom_shp=None,
        explicit_field="BIOME_NAME", explicit_geo_category="biome",
        config_path=cfg)
    assert res.field == "BIOME_NAME"
    assert res.geo_category == "biome"


def test_resolve_selection_uses_registered_when_no_overrides(tmp_path):
    cfg = _write_config(tmp_path, (
        "shapefiles:\n  Eco:\n    path: eco.shp\n    field: REALM\n"
        "    geo_category: biogeographic_realm\n"
    ))
    res = sr.resolve_selection(
        shp_name="Eco", custom_shp=None,
        explicit_field=None, explicit_geo_category=None, config_path=cfg)
    assert res.field == "REALM"
    assert res.geo_category == "biogeographic_realm"


# ---------------------------------------------------------------------------
# Task 5: config template generation
# ---------------------------------------------------------------------------

def test_write_template_creates_parseable_file(tmp_path):
    target = tmp_path / "cfg" / "shapefiles.yaml"
    out = sr.write_template(target)
    assert out == target and target.exists()
    data = _yaml.safe_load(target.read_text())
    assert isinstance(data, dict) and "shapefiles" in data
    assert "GOaS_v1_20211214" in data["shapefiles"]


def test_write_template_refuses_overwrite(tmp_path):
    target = tmp_path / "shapefiles.yaml"
    target.write_text("shapefiles: {}\n")
    with pytest.raises(FileExistsError):
        sr.write_template(target)
