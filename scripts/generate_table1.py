#!/usr/bin/env python3
"""
Generate Table 1 for the BOLDGenotyper MER manuscript.

Table 1: Case Study Benchmark Summary — one row per dataset (5 datasets).

Data sources:
  - output/{dataset}_full/reports/{dataset}_assignment_summary.csv
  - output/{dataset}_full/geographic_analysis/geographic_coverage.csv
  - benchmarking/parameter_sweep/{dataset}/recommendations.txt
  - output/{dataset}_full/{dataset}_pipeline.log  (for runtime)

Output:
    scripts/Table1.csv   — machine-readable
    scripts/Table1.md    — markdown formatted (for manuscript)
"""

import re
from pathlib import Path
from datetime import datetime

import pandas as pd

# ---------------------------------------------------------------------------
# Dataset metadata (fixed / documented)
# ---------------------------------------------------------------------------
DATASETS = [
    dict(
        name        = 'Sphyrnidae',
        taxon       = 'Sphyrnidae',
        rank        = 'Family',
        environment = 'Marine',
        shapefile   = 'GOaS v1',
        geo_col     = 'with_ocean_basin',
        region_type = 'ocean basins',
        validation  = 'Smith & Black (2026)',
    ),
    dict(
        name        = 'Carcharhiniformes',
        taxon       = 'Carcharhiniformes',
        rank        = 'Order',
        environment = 'Marine',
        shapefile   = 'GOaS v1',
        geo_col     = 'with_ocean_basin',
        region_type = 'ocean basins',
        validation  = '—',
    ),
    dict(
        name        = 'Panulirus',
        taxon       = 'Panulirus',
        rank        = 'Genus',
        environment = 'Marine',
        shapefile   = 'GOaS v1',
        geo_col     = 'with_ocean_basin',
        region_type = 'ocean basins',
        validation  = '—',
    ),
    dict(
        name        = 'Salmonidae',
        taxon       = 'Salmonidae',
        rank        = 'Family',
        environment = 'Freshwater',
        shapefile   = 'BasinATLAS v10 (lev. 3)',
        geo_col     = 'with_freshwater_basin',
        region_type = 'drainage basins',
        validation  = '—',
    ),
    dict(
        name        = 'Pieridae',
        taxon       = 'Pieridae',
        rank        = 'Family',
        environment = 'Terrestrial',
        shapefile   = 'WWF Ecoregions 2017',
        geo_col     = 'with_ecoregion',
        region_type = 'biogeographic realms',
        validation  = '—',
    ),
]

OUTPUT_DIR   = Path('output')
SWEEP_ROOT   = Path('benchmarking/parameter_sweep')

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def parse_runtime(name: str) -> float:
    """Return pipeline runtime in minutes from pipeline log timestamps."""
    log = OUTPUT_DIR / f'{name}_full' / f'{name}_pipeline.log'
    lines = log.read_text().splitlines()
    def ts(line):
        m = re.match(r'\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]', line)
        return datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S') if m else None
    first = next((ts(l) for l in lines if ts(l)), None)
    last  = next((ts(l) for l in reversed(lines) if ts(l)), None)
    if first and last:
        return round((last - first).total_seconds() / 60, 1)
    return float('nan')


def parse_elbow_threshold(name: str) -> tuple[float, int]:
    """Return (optimal_threshold, n_haplotypes_at_elbow)."""
    text = (SWEEP_ROOT / name / 'recommendations.txt').read_text()
    t  = float(re.search(r'PRIMARY:\s+min_singleton_distance\s*=\s*([0-9.]+)', text).group(1))
    m  = re.search(r'At [0-9.]+ \(elbow\):\s+([0-9.]+)\s+haplotypes', text)
    nh = int(float(m.group(1))) if m else None
    return t, nh


def count_unique_regions(name: str, geo_col: str) -> int:
    """Count unique occupied geographic regions from the annotated CSV."""
    ann = OUTPUT_DIR / f'{name}_full' / f'{name}_annotated.csv'
    col_map = {
        'with_ocean_basin':    'ocean_basin',
        'with_freshwater_basin': 'freshwater_basin',
        'with_ecoregion':      'ecoregion',
    }
    col = col_map.get(geo_col)
    df  = pd.read_csv(ann, low_memory=False, usecols=[col])
    vals = df[col].dropna().unique()
    # Exclude 'Unknown', 'unknown', and multi-value strings (e.g. 'Nearctic, Neotropic')
    clean = [v for v in vals if str(v) not in ('Unknown', 'unknown') and ',' not in str(v)]
    return len(clean)


# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
rows = []
for ds in DATASETS:
    name    = ds['name']
    geo_col = ds['geo_col']

    # Assignment summary
    summ = pd.read_csv(
        OUTPUT_DIR / f'{name}_full' / 'reports' / f'{name}_assignment_summary.csv'
    ).iloc[0]
    n_input     = int(summ['total_samples'])
    n_with_seq  = int(summ['samples_with_sequence'])
    n_haplotypes = int(summ['n_haplotypes'])
    pct_with_seq = round(n_with_seq / n_input * 100, 1)

    # Geographic coverage
    geo = pd.read_csv(
        OUTPUT_DIR / f'{name}_full' / 'geographic_analysis' / 'geographic_coverage.csv'
    )
    n_geo_assigned = int(geo[geo_col].sum())
    pct_geo        = round(n_geo_assigned / n_with_seq * 100, 1)

    # Unique regions
    n_regions = count_unique_regions(name, geo_col)

    # Parameter sweep
    threshold, n_hap_elbow = parse_elbow_threshold(name)

    # Runtime
    runtime_min = parse_runtime(name)

    rows.append({
        'Taxon':                   ds['taxon'],
        'Rank':                    ds['rank'],
        'Environment':             ds['environment'],
        'Input records':           n_input,
        'COI sequences':           n_with_seq,
        'COI coverage (%)':        pct_with_seq,
        'Haplotypes (full run)':   n_haplotypes,
        'Haplotypes (at elbow)':   n_hap_elbow,
        'Geo-assigned':            n_geo_assigned,
        'Geo-assigned (% of COI)': pct_geo,
        'Geographic regions':      n_regions,
        'Region type':             ds['region_type'],
        'Shapefile':               ds['shapefile'],
        'Optimal threshold':       threshold,
        'Runtime (min)':           runtime_min,
        'Validation':              ds['validation'],
    })

df = pd.DataFrame(rows)

# Sort: marine first (by N input desc), then freshwater, terrestrial
env_order = {'Marine': 0, 'Freshwater': 1, 'Terrestrial': 2}
df['_env_order'] = df['Environment'].map(env_order)
df = df.sort_values(['_env_order', 'Input records'], ascending=[True, False]).drop(columns='_env_order')
df = df.reset_index(drop=True)

# ---------------------------------------------------------------------------
# Save CSV
# ---------------------------------------------------------------------------
out_csv = Path('scripts/Table1.csv')
df.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")

# ---------------------------------------------------------------------------
# Save Markdown
# ---------------------------------------------------------------------------
# Condensed columns for the markdown table
md_df = df[[
    'Taxon', 'Rank', 'Environment',
    'Input records', 'COI sequences', 'COI coverage (%)',
    'Haplotypes (full run)', 'Haplotypes (at elbow)', 'Optimal threshold',
    'Geographic regions', 'Region type', 'Shapefile',
    'Geo-assigned', 'Geo-assigned (% of COI)',
    'Runtime (min)', 'Validation',
]].copy()

md_lines = ['# Table 1: BOLDGenotyper Case Study Summary\n']
def df_to_md(df):
    cols = df.columns.tolist()
    header = '| ' + ' | '.join(str(c) for c in cols) + ' |'
    sep    = '|' + '|'.join(['---'] * len(cols)) + '|'
    rows_  = ['| ' + ' | '.join(str(v) for v in row) + ' |' for _, row in df.iterrows()]
    return '\n'.join([header, sep] + rows_)

md_lines.append(df_to_md(md_df))
md_lines.append('\n')
md_lines.append('**Notes:**\n')
md_lines.append('- COI coverage: % of input records that contain a usable COI sequence.\n')
md_lines.append('- Haplotypes (full run): at default singleton-distance threshold (0.005).\n')
md_lines.append('- Haplotypes (at elbow): at the optimal threshold identified by parameter sweep.\n')
md_lines.append('- Geo-assigned: number of COI sequences spatially assigned to a geographic region.\n')
md_lines.append('- Geo-assigned %: proportion of COI sequences with a geographic region assignment.\n')
md_lines.append('- Panulirus: only 1/163 COI records carry geographic coordinates (BOLD data quality).\n')
md_lines.append('- Runtime on Apple M2 MacBook Pro (16 GB RAM); single-threaded pipeline.\n')
md_lines.append('- GOaS = General Ocean Atlas of Seas (10 named basins); BasinATLAS = HydroSHEDS basin atlas (level 3, 292 basins); WWF Ecoregions 2017 (REALM field, 8 realms).\n')

out_md = Path('scripts/Table1.md')
out_md.write_text('\n'.join(md_lines))
print(f"Saved: {out_md}")

# Print to console for quick check
print('\n' + '='*80)
print(md_df.to_string(index=False))
