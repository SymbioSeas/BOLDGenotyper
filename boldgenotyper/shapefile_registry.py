"""Shapefile registry: register shapefiles by shorthand in shapefiles.yaml.

Turns a shorthand name plus optional overrides into a resolved shapefile
selection, so users need not pass --custom-shp/--shp-field/--geo-category on
every run. Purely additive: nothing here changes default GOaS behavior.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

ENV_VAR = "BOLDGENOTYPER_SHAPEFILE_CONFIG"


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


def default_config_path() -> Path:
    """Return the default shapefiles.yaml location for this platform."""
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
        return base / "boldgenotyper" / "shapefiles.yaml"
    xdg = os.environ.get("XDG_CONFIG_HOME")
    base = Path(xdg) if xdg else Path.home() / ".config"
    return base / "boldgenotyper" / "shapefiles.yaml"


def _resolve_path(raw: str, shapefiles_dir: Optional[Path], config_dir: Path) -> Path:
    """Resolve a shapefile path so it works from any working directory."""
    expanded = Path(os.path.expandvars(str(raw))).expanduser()
    if expanded.is_absolute():
        return expanded
    base = shapefiles_dir if shapefiles_dir is not None else config_dir
    return (base / expanded).resolve()


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
