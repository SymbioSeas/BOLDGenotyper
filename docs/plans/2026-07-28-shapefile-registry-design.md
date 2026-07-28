# Shapefile registry (shapefiles.yaml) — design

Date: 2026-07-28
Status: Approved, pending implementation plan

## Problem

Custom (non-GOaS) geographic analysis currently requires three flags on every
run: `--custom-shp <path>`, `--shp-field <field>`, and `--geo-category <label>`.
Most users work with one or a few shapefiles repeatedly, so re-supplying all
three each time is tedious and error-prone. The paths are also commonly written
relative to the repository root (e.g. `shapefiles/BasinATLAS/...`), which breaks
when BOLDGenotyper is run from a different working directory — the normal case
once it is installed and used across multiple projects.

## Goals

- Let users register a shapefile once under a shorthand name and select it with a
  single flag: `boldgenotyper data.tsv --shp Ecoregions2017`.
- Work from any working directory, across different projects, without editing the
  command per project.
- Support multiple shapefiles in one config file.
- Present one consistent shorthand mechanism across marine, freshwater, and
  terrestrial data (GOaS included).
- Be purely additive: existing behavior (GOaS default, the three explicit flags)
  is unchanged when the new flag and config are not used.

## Non-goals

- Re-routing GOaS through the generic shapefile assignment. GOaS keeps its
  existing ocean-basin code path to avoid any risk of changing marine results.
- Downloading or managing shapefile data itself (still a manual user step,
  documented in the README).
- A `--list-shapefiles` command (deferred; the unknown-name error already lists
  registered names).

## Design

### Config discovery

Resolution order (first found wins):

1. `--shapefile-config <path>` flag
2. `BOLDGENOTYPER_SHAPEFILE_CONFIG` environment variable
3. Default location: `~/.config/boldgenotyper/shapefiles.yaml`
   (`%APPDATA%\boldgenotyper\shapefiles.yaml` on Windows)

If no config is found and `--shp` is not used, behavior is unchanged.

### Schema (one file, multiple entries)

```yaml
# Optional base dir for resolving relative paths below.
# If omitted, relative paths resolve against this file's own directory.
shapefiles_dir: ~/boldgenotyper-data/shapefiles

shapefiles:
  GOaS_v1_20211214:
    path: GOaS_v1_20211214/goas_v01.shp
    field: name
    geo_category: ocean_basin
    type: goas            # routes to the marine ocean-basin logic
  Ecoregions2017:
    path: Ecoregions2017/Ecoregions2017.shp
    field: REALM
    geo_category: biogeographic_realm
  BasinATLAS_v10:
    path: BasinATLAS/BasinATLAS_lev03.shp
    field: PFAF_ID
    geo_category: drainage_basin
```

Per entry:

- `path` (required)
- `field` (optional; auto-detected via existing logic if omitted)
- `geo_category` (optional; defaults to `geographic_region`, or `ocean_basin`
  for `type: goas`)
- `type` (optional; `generic` by default, or `goas` for the marine routing)

GOaS ships pre-registered so `--shp GOaS_v1_20211214` works uniformly, while GOaS
remains the default when no `--shp` is given.

### Path resolution

For each `path`, resolve so it works from any working directory:

- Absolute paths used as-is.
- `~` and environment variables expanded.
- Relative paths resolved against `shapefiles_dir` if set, otherwise against the
  config file's own directory.

This lets users keep the config and the (large) shapefile data wherever they
like, together or apart.

### Invocation and flag precedence

- New flag `--shp NAME` looks NAME up in the registry.
- The three existing flags remain and take precedence over registered values, so
  `--shp Ecoregions2017 --geo-category realm` uses the registered path/field but
  overrides the category.
- `--shp` together with `--custom-shp` is an error (ambiguous source).
- An unknown NAME is an error that lists the registered names.

Note: `--shp` and the existing `--shp-field` share a prefix. argparse matches
exact option strings before abbreviations, so `--shp NAME` resolves to the new
flag unambiguously; the implementation must confirm this holds (and that no other
flag abbreviation regresses) rather than assume it.

### Code structure

New focused module `boldgenotyper/shapefile_registry.py`:

- `load_registry(config_path: Optional[Path] = None) -> Registry` — discovery,
  parse, and validation.
- `resolve_shapefile(name: str, registry: Registry) -> ResolvedShapefile` —
  returns `(path, field, geo_category, is_goas)`.

`cli.py` `main()` calls the resolver after argument parsing and populates
`shapefile_path` / `shapefile_field` / `geo_category` using the precedence above.
Keeping this logic in its own module keeps it out of the already-large `cli.py`
and makes it independently testable.

### Scaffolding command

`boldgenotyper --init-shapefile-config` writes a commented template (the GOaS
entry plus freshwater/terrestrial examples) to the target config path — the
`--shapefile-config` value if given, otherwise the default location — creating
parent directories as needed. It refuses to overwrite an existing file and tells
the user to edit or remove it instead. This follows the existing
`create_config_template()` pattern. Without it, users would have no easy way to
learn the schema.

## Testing

Unit tests cover:

- Discovery and precedence (flag, env var, default location).
- Path resolution: absolute, relative-to-`shapefiles_dir`, relative-to-config,
  `~` and environment-variable expansion.
- Name resolution: found, and not-found error listing registered names.
- Flag-override precedence (explicit flags beat registered values).
- `--shp` + `--custom-shp` conflict error.
- GOaS `type: goas` routing.

## Documentation

Update the README shapefile section to present the registry workflow as the
recommended path, keeping the explicit-flags workflow documented as the
alternative. Document config discovery order, schema, and `--init-shapefile-config`.
