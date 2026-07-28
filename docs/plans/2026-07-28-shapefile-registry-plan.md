# Shapefile Registry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let users register shapefiles once by shorthand in `shapefiles.yaml` and select one with `--shp NAME`, working from any directory.

**Architecture:** A new self-contained module `boldgenotyper/shapefile_registry.py` handles config discovery, parsing, path resolution, name resolution, CLI-precedence, and template generation. `cli.py` calls it after argument parsing and populates the existing `shapefile_path`/`shapefile_field`/`geo_category` values, leaving `run_pipeline` and the geographic code untouched. The feature is purely additive.

**Tech Stack:** Python 3.9+, PyYAML (already a dependency), argparse, pytest.

## Global Constraints

- Python floor: 3.9 (no 3.10+ syntax such as `X | Y` unions or `match`).
- Dependencies: use only what is already declared (PyYAML `>=5.0.0` is available); add no new dependency.
- No emojis in any code, logging, comments, or docs.
- Cross-platform: default config path must resolve correctly on macOS, Linux, and Windows.
- Purely additive: when neither `--shp` nor a config file is used, existing behavior (GOaS default, `--custom-shp`/`--shp-field`/`--geo-category`) is unchanged.
- Follow existing config idioms (mirror `config.py`'s `create_config_template()` / `load_config_from_file()` style).

## File Structure

- Create `boldgenotyper/shapefile_registry.py` — discovery, parse, path/name resolution, CLI precedence, template. One responsibility: turning a shorthand name + flags into a resolved shapefile selection.
- Create `tests/test_shapefile_registry.py` — unit tests for the module.
- Modify `boldgenotyper/cli.py` — add `--shp`, `--shapefile-config`, `--init-shapefile-config`; wire resolution into `main()`.
- Modify `README.md` — document the registry workflow.

---

### Task 1: Data models and path resolution

**Files:**
- Create: `boldgenotyper/shapefile_registry.py`
- Test: `tests/test_shapefile_registry.py`

**Interfaces:**
- Produces: `ResolvedShapefile(name: str, path: Path, field: Optional[str], geo_category: str, is_goas: bool)` (frozen dataclass); `Registry(entries: Dict[str, Any], shapefiles_dir: Optional[Path], config_dir: Path, source: Optional[Path])` (frozen dataclass); `_resolve_path(raw: str, shapefiles_dir: Optional[Path], config_dir: Path) -> Path`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_shapefile_registry.py
from pathlib import Path
import os
import pytest
from boldgenotyper import shapefile_registry as sr


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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'boldgenotyper.shapefile_registry'`

- [ ] **Step 3: Write minimal implementation**

```python
# boldgenotyper/shapefile_registry.py
"""Shapefile registry: register shapefiles by shorthand in shapefiles.yaml.

Turns a shorthand name plus optional overrides into a resolved shapefile
selection, so users need not pass --custom-shp/--shp-field/--geo-category on
every run. Purely additive: nothing here changes default GOaS behavior.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class ResolvedShapefile:
    name: str
    path: Path
    field: Optional[str]
    geo_category: str
    is_goas: bool


@dataclass(frozen=True)
class Registry:
    entries: Dict[str, Any]
    shapefiles_dir: Optional[Path]
    config_dir: Path
    source: Optional[Path]


def _resolve_path(raw: str, shapefiles_dir: Optional[Path], config_dir: Path) -> Path:
    """Resolve a shapefile path so it works from any working directory."""
    expanded = Path(os.path.expandvars(str(raw))).expanduser()
    if expanded.is_absolute():
        return expanded
    base = shapefiles_dir if shapefiles_dir is not None else config_dir
    return (base / expanded).resolve()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/shapefile_registry.py tests/test_shapefile_registry.py
git commit -m "Add shapefile registry models and path resolution"
```

---

### Task 2: Config discovery and loading

**Files:**
- Modify: `boldgenotyper/shapefile_registry.py`
- Test: `tests/test_shapefile_registry.py`

**Interfaces:**
- Consumes: `Registry`, `_resolve_path` (Task 1).
- Produces: `default_config_path() -> Path`; `ENV_VAR = "BOLDGENOTYPER_SHAPEFILE_CONFIG"`; `load_registry(config_path: Optional[Path] = None) -> Registry`. Discovery precedence: `config_path` arg > `ENV_VAR` env var > `default_config_path()`. A user-specified path (arg or env) that is missing raises `FileNotFoundError`; a missing default returns an empty `Registry`.

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_shapefile_registry.py

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: FAIL with `AttributeError: module 'boldgenotyper.shapefile_registry' has no attribute 'ENV_VAR'`

- [ ] **Step 3: Write minimal implementation**

```python
# add near the top of boldgenotyper/shapefile_registry.py, after imports
import yaml

ENV_VAR = "BOLDGENOTYPER_SHAPEFILE_CONFIG"


def default_config_path() -> Path:
    """Return the default shapefiles.yaml location for this platform."""
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
        return base / "boldgenotyper" / "shapefiles.yaml"
    xdg = os.environ.get("XDG_CONFIG_HOME")
    base = Path(xdg) if xdg else Path.home() / ".config"
    return base / "boldgenotyper" / "shapefiles.yaml"


def _resolve_dir(raw: str, config_dir: Path) -> Path:
    d = Path(os.path.expandvars(str(raw))).expanduser()
    return d if d.is_absolute() else (config_dir / d).resolve()


def load_registry(config_path: Optional[Path] = None) -> Registry:
    """Discover and load the shapefile registry.

    Precedence: config_path arg > ENV_VAR env var > default_config_path().
    A user-specified path that is missing raises FileNotFoundError; a missing
    default returns an empty Registry so existing behavior is unaffected.
    """
    explicit = config_path is not None or bool(os.environ.get(ENV_VAR))
    if config_path is not None:
        path = Path(config_path).expanduser()
    elif os.environ.get(ENV_VAR):
        path = Path(os.environ[ENV_VAR]).expanduser()
    else:
        path = default_config_path()

    if not path.exists():
        if explicit:
            raise FileNotFoundError(f"Shapefile config not found: {path}")
        return Registry(entries={}, shapefiles_dir=None, config_dir=path.parent, source=None)

    with open(path) as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Shapefile config {path} must be a mapping at the top level.")
    entries = data.get("shapefiles", {}) or {}
    if not isinstance(entries, dict):
        raise ValueError(f"'shapefiles' in {path} must be a mapping of name to settings.")

    config_dir = path.parent
    raw_dir = data.get("shapefiles_dir")
    shp_dir = _resolve_dir(raw_dir, config_dir) if raw_dir else None
    return Registry(entries=entries, shapefiles_dir=shp_dir, config_dir=config_dir, source=path)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: PASS (9 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/shapefile_registry.py tests/test_shapefile_registry.py
git commit -m "Add shapefile config discovery and loading"
```

---

### Task 3: Name resolution

**Files:**
- Modify: `boldgenotyper/shapefile_registry.py`
- Test: `tests/test_shapefile_registry.py`

**Interfaces:**
- Consumes: `Registry`, `ResolvedShapefile`, `_resolve_path` (Tasks 1-2).
- Produces: `resolve_shapefile(name: str, registry: Registry) -> ResolvedShapefile`. Unknown name raises `ValueError` listing registered names. Missing `field` stays `None` (auto-detect downstream). `geo_category` defaults to `ocean_basin` when `type: goas`, else `geographic_region`. `is_goas` is True when `type` equals `goas` (case-insensitive).

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_shapefile_registry.py

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -k resolve_shapefile -q`
Expected: FAIL with `AttributeError: ... has no attribute 'resolve_shapefile'`

- [ ] **Step 3: Write minimal implementation**

```python
# add to boldgenotyper/shapefile_registry.py

def resolve_shapefile(name: str, registry: Registry) -> ResolvedShapefile:
    """Resolve a registered shorthand name into a ResolvedShapefile."""
    if name not in registry.entries:
        available = ", ".join(sorted(registry.entries)) or "(none registered)"
        source = registry.source or "(no config file found)"
        raise ValueError(
            f"Shapefile '{name}' not found in registry. Registered: {available}. "
            f"Config: {source}"
        )
    entry = registry.entries[name]
    if not isinstance(entry, dict) or "path" not in entry:
        raise ValueError(f"Shapefile '{name}' must define at least a 'path'.")

    is_goas = str(entry.get("type", "generic")).lower() == "goas"
    path = _resolve_path(entry["path"], registry.shapefiles_dir, registry.config_dir)
    field = entry.get("field")
    geo_category = entry.get("geo_category") or ("ocean_basin" if is_goas else "geographic_region")
    return ResolvedShapefile(name=name, path=path, field=field, geo_category=geo_category, is_goas=is_goas)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: PASS (13 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/shapefile_registry.py tests/test_shapefile_registry.py
git commit -m "Add shapefile name resolution with GOaS routing"
```

---

### Task 4: CLI selection precedence and conflict handling

**Files:**
- Modify: `boldgenotyper/shapefile_registry.py`
- Test: `tests/test_shapefile_registry.py`

**Interfaces:**
- Consumes: `load_registry`, `resolve_shapefile`, `ResolvedShapefile` (Tasks 1-3).
- Produces: `resolve_selection(*, shp_name: Optional[str], custom_shp: Optional[Path], explicit_field: Optional[str], explicit_geo_category: Optional[str], config_path: Optional[Path] = None) -> Optional[ResolvedShapefile]`. Returns `None` when `shp_name` is None (caller keeps existing behavior). Raises `ValueError` if both `shp_name` and `custom_shp` are given. Explicit `explicit_field`/`explicit_geo_category` (non-None) override the registered values.

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_shapefile_registry.py

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -k resolve_selection -q`
Expected: FAIL with `AttributeError: ... has no attribute 'resolve_selection'`

- [ ] **Step 3: Write minimal implementation**

```python
# add to boldgenotyper/shapefile_registry.py

def resolve_selection(
    *,
    shp_name: Optional[str],
    custom_shp: Optional[Path],
    explicit_field: Optional[str],
    explicit_geo_category: Optional[str],
    config_path: Optional[Path] = None,
) -> Optional[ResolvedShapefile]:
    """Resolve the --shp selection, applying explicit-flag overrides.

    Returns None when --shp was not given. Raises ValueError on --shp +
    --custom-shp conflict or an unresolvable name.
    """
    if shp_name is None:
        return None
    if custom_shp is not None:
        raise ValueError("Use either --shp NAME or --custom-shp PATH, not both.")
    registry = load_registry(config_path)
    resolved = resolve_shapefile(shp_name, registry)
    field = explicit_field if explicit_field is not None else resolved.field
    geo_category = explicit_geo_category if explicit_geo_category is not None else resolved.geo_category
    return ResolvedShapefile(
        name=resolved.name, path=resolved.path, field=field,
        geo_category=geo_category, is_goas=resolved.is_goas,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: PASS (17 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/shapefile_registry.py tests/test_shapefile_registry.py
git commit -m "Add CLI shapefile selection precedence and conflict handling"
```

---

### Task 5: Config template generation

**Files:**
- Modify: `boldgenotyper/shapefile_registry.py`
- Test: `tests/test_shapefile_registry.py`

**Interfaces:**
- Consumes: nothing from earlier tasks except module import.
- Produces: `TEMPLATE: str`; `write_template(target_path: Path) -> Path`. Creates parent dirs, refuses to overwrite an existing file (`FileExistsError`), returns the written path. `TEMPLATE` must be valid YAML that `load_registry` can parse.

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_shapefile_registry.py
import yaml as _yaml


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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -k template -q`
Expected: FAIL with `AttributeError: ... has no attribute 'write_template'`

- [ ] **Step 3: Write minimal implementation**

```python
# add to boldgenotyper/shapefile_registry.py

TEMPLATE = """# BOLDGenotyper shapefile registry
#
# Register shapefiles here once, then select one per run with: --shp NAME
# Discovery order: --shapefile-config PATH > BOLDGENOTYPER_SHAPEFILE_CONFIG env
# var > this file (~/.config/boldgenotyper/shapefiles.yaml).
#
# Optional base directory for resolving the relative 'path' values below.
# If omitted, relative paths resolve against this file's own directory.
# shapefiles_dir: ~/boldgenotyper-data/shapefiles

shapefiles:
  # Marine default (Global Oceans and Seas). Kept here so --shp GOaS_v1_20211214
  # works uniformly; GOaS is still used automatically when no --shp is given.
  GOaS_v1_20211214:
    path: GOaS_v1_20211214/goas_v01.shp
    field: name
    geo_category: ocean_basin
    type: goas

  # Terrestrial example (Ecoregions2017).
  # Ecoregions2017:
  #   path: Ecoregions2017/Ecoregions2017.shp
  #   field: REALM
  #   geo_category: biogeographic_realm

  # Freshwater example (BasinATLAS / HydroBASINS).
  # BasinATLAS_v10:
  #   path: BasinATLAS/BasinATLAS_lev03.shp
  #   field: PFAF_ID
  #   geo_category: drainage_basin
"""


def write_template(target_path: Path) -> Path:
    """Write a commented shapefiles.yaml template. Refuses to overwrite."""
    target_path = Path(target_path).expanduser()
    if target_path.exists():
        raise FileExistsError(
            f"{target_path} already exists. Edit or remove it instead of re-initialising."
        )
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(TEMPLATE)
    return target_path
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: PASS (19 passed)

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/shapefile_registry.py tests/test_shapefile_registry.py
git commit -m "Add shapefile config template generation"
```

---

### Task 6: Wire the registry into the CLI

**Files:**
- Modify: `boldgenotyper/cli.py`
- Test: `tests/test_shapefile_registry.py` (argparse + integration behavior)

**Interfaces:**
- Consumes: `shapefile_registry.resolve_selection`, `.write_template`, `.default_config_path` (Tasks 2-5).
- Produces: three new CLI flags on the main parser — `--shp` (dest `shp_name`), `--shapefile-config` (dest `shapefile_config`, `type=Path`), `--init-shapefile-config` (`store_true`, dest `init_shapefile_config`); and resolution logic in `main()`.

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_shapefile_registry.py
from boldgenotyper import cli as _cli


def test_cli_parses_shp_flags_without_prefix_collision(monkeypatch):
    # --shp must resolve to shp_name, not collide with --shp-field.
    import sys
    argv = ["prog", "data.tsv", "--shp", "Eco", "--shp-field", "REALM"]
    monkeypatch.setattr(sys, "argv", argv)
    # Build the same parser main() uses by calling parse via a helper is not
    # exposed; instead assert the option strings exist on a fresh parse.
    # We import argparse indirectly by invoking the parser construction path:
    parser = _cli._build_main_parser()
    ns = parser.parse_args(["data.tsv", "--shp", "Eco", "--shp-field", "REALM"])
    assert ns.shp_name == "Eco"
    assert ns.shapefile_field == "REALM"
```

Note: this test requires `main()`'s parser construction to be extracted into a
`_build_main_parser()` helper (Step 3) so the parser is testable without running
the pipeline.

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -k prefix_collision -q`
Expected: FAIL with `AttributeError: module 'boldgenotyper.cli' has no attribute '_build_main_parser'`

- [ ] **Step 3: Implement — extract parser, add flags, wire resolution**

3a. In `boldgenotyper/cli.py`, add the import near the other `from . import (...)` block:

```python
from . import shapefile_registry
```

3b. Refactor `main()` so the `argparse.ArgumentParser(...)` construction and all
`parser.add_argument(...)` calls move into a new module-level function
`_build_main_parser() -> argparse.ArgumentParser` that returns the parser.
`main()` then begins with:

```python
def main():
    """Main CLI entry point."""
    parser = _build_main_parser()
    args = parser.parse_args()
```

3c. Inside `_build_main_parser()`, immediately after the existing `--geo-category`
argument block, add the three new arguments:

```python
    parser.add_argument(
        '--shp',
        type=str,
        default=None,
        dest='shp_name',
        help='Shorthand name of a shapefile registered in your shapefiles.yaml '
             '(see --init-shapefile-config). Resolves its path, field, and '
             'geo-category. Explicit --shp-field/--geo-category override the '
             'registered values. Cannot be combined with --custom-shp.'
    )
    parser.add_argument(
        '--shapefile-config',
        type=Path,
        default=None,
        dest='shapefile_config',
        help='Path to a shapefile registry file. Overrides the '
             'BOLDGENOTYPER_SHAPEFILE_CONFIG env var and the default '
             '(~/.config/boldgenotyper/shapefiles.yaml).'
    )
    parser.add_argument(
        '--init-shapefile-config',
        action='store_true',
        dest='init_shapefile_config',
        help='Write a commented shapefiles.yaml template to the target config '
             'path (--shapefile-config if given, else the default location) and exit.'
    )
```

3d. In `main()`, immediately after `args = parser.parse_args()` and BEFORE the
`# Validate input file` block, add the init handling and registry resolution:

```python
    # Handle --init-shapefile-config: write template and exit.
    if args.init_shapefile_config:
        target = args.shapefile_config or shapefile_registry.default_config_path()
        try:
            written = shapefile_registry.write_template(target)
        except FileExistsError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        print(f"Wrote shapefile config template to: {written}")
        print("Edit it to register your shapefiles, then use: --shp <name>")
        return 0

    # Resolve a registered shapefile shorthand (--shp), if given.
    goas_path_override = None
    try:
        selection = shapefile_registry.resolve_selection(
            shp_name=args.shp_name,
            custom_shp=args.shapefile_path,
            explicit_field=args.shapefile_field,
            explicit_geo_category=args.geo_category,
            config_path=args.shapefile_config,
        )
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if selection is not None:
        args.geo_category = selection.geo_category
        if selection.is_goas:
            # Route through the existing GOaS ocean-basin logic; override the
            # GOaS shapefile path so the registered copy is used.
            args.shapefile_path = None
            goas_path_override = selection.path
        else:
            args.shapefile_path = selection.path
            args.shapefile_field = selection.field
```

3e. Where `cfg` is built via `cfg = cfg.update(...)`, apply the GOaS override
immediately after that block:

```python
    if goas_path_override is not None:
        cfg = cfg.update(geographic__goas_shapefile_path=goas_path_override)
```

- [ ] **Step 4: Run the new test plus the full registry suite and a CLI smoke check**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m pytest tests/test_shapefile_registry.py -q`
Expected: PASS (20 passed)

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -c "import sys; sys.argv=['boldgenotyper','--help']; from boldgenotyper import cli; cli.main()" 2>&1 | grep -E "\-\-shp|--shapefile-config|--init-shapefile-config"`
Expected: the three flags appear in help.

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -m py_compile boldgenotyper/cli.py`
Expected: no output (compiles).

- [ ] **Step 5: Commit**

```bash
git add boldgenotyper/cli.py tests/test_shapefile_registry.py
git commit -m "Wire shapefile registry into the CLI (--shp, --shapefile-config, --init-shapefile-config)"
```

---

### Task 7: Documentation

**Files:**
- Modify: `README.md`

**Interfaces:** none (docs only).

- [ ] **Step 1: Update the README shapefile section**

In `README.md`, in the "Custom Shapefiles for Non-Marine Organisms" section (around the `--custom-shp`/`--shp-field`/`--geo-category` documentation), add a new subsection *before* the explicit-flags examples titled "Registering shapefiles (recommended)". It must state:

- Run `boldgenotyper --init-shapefile-config` to create `~/.config/boldgenotyper/shapefiles.yaml`.
- Show the schema (copy the block from `docs/plans/2026-07-28-shapefile-registry-design.md` section "Schema").
- Explain discovery order: `--shapefile-config` > `BOLDGENOTYPER_SHAPEFILE_CONFIG` > default location.
- Explain path resolution: relative paths resolve against `shapefiles_dir` if set, else the config file's directory; absolute paths and `~`/env vars work.
- Show usage: `boldgenotyper data/Salmonidae.tsv --shp BasinATLAS_v10 --build-tree`.
- Note that explicit `--shp-field`/`--geo-category` override registered values, and that `--shp` and `--custom-shp` cannot be combined.
- Keep the existing explicit-flags examples, framed as the alternative when you do not want a registry.

- [ ] **Step 2: Verify the documented commands are accurate**

Run: `/Users/stephsmith/miniforge3/envs/boldgenotyper/bin/python -c "import sys; sys.argv=['boldgenotyper','--init-shapefile-config','--shapefile-config','/tmp/bg_demo.yaml']; from boldgenotyper import cli; raise SystemExit(cli.main())"`
Expected: prints "Wrote shapefile config template to: /tmp/bg_demo.yaml"; then confirm with `cat /tmp/bg_demo.yaml` that it matches the documented schema. Clean up: `rm /tmp/bg_demo.yaml`.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "Document shapefile registry workflow in README"
```

---

## Self-Review

**Spec coverage:**
- Config discovery + precedence -> Task 2. Schema/GOaS pre-registration -> Tasks 2-3, 5 (template). Path resolution -> Task 1. Invocation + precedence + conflict -> Tasks 4, 6. Code structure (new module) -> Tasks 1-5. Scaffolding command -> Tasks 5-6. Docs -> Task 7. Tests -> every task. Argparse `--shp`/`--shp-field` prefix note -> Task 6 test. All spec sections covered.

**Placeholder scan:** No TBD/TODO; every code step shows complete code; test steps show real assertions.

**Type consistency:** `ResolvedShapefile` and `Registry` field names are consistent across Tasks 1-6. `resolve_selection` keyword names (`shp_name`, `custom_shp`, `explicit_field`, `explicit_geo_category`, `config_path`) match between Task 4 definition and Task 6 caller. `dest` names (`shp_name`, `shapefile_config`, `init_shapefile_config`) match between the argparse definitions and the `main()` usage.

**Note for implementer:** Task 6 Step 3b is a mechanical refactor (moving parser construction into `_build_main_parser()`); run the existing CLI `--help` before and after to confirm output is unchanged apart from the three new flags.
