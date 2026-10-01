# BOLDGenotyper Tutorial Notebooks

Interactive Jupyter notebooks for common BOLDGenotyper 2.0.0 workflows.

## Prerequisites

- BOLDGenotyper installed (see the main README.md)
- Jupyter Notebook or JupyterLab
- A BOLD TSV download for your taxon (see the README section "Obtaining Input Data from
  BOLD"). The repository does not ship BOLD downloads: set `INPUT_TSV` at the top of each
  notebook.

Each notebook writes its runs to `tutorial_output/` next to it, one folder per run.

## Installation

```bash
conda activate boldgenotyper
conda install -c conda-forge jupyter
```

## Launching the Tutorials

```bash
jupyter lab notebooks/
```

## Tutorial Series

### Tutorial 1: Basic Haplotype Workflow
**File**: `01_basic_genotyping_workflow.ipynb`

Run an analysis with an explicit `--singleton-distance`, then read the run folder: its
`README.md`, `summary.csv`, record statuses in `records.csv`, haplotypes in
`haplotypes.csv`, member labels and label conflicts, regions, provenance and the HTML
report.

### Tutorial 2: Choosing `--singleton-distance`
**File**: `02_parameter_sweep_optimization.ipynb`

Run `boldgenotyper-sweep`, find the elbow of the haplotype-count curve, read the stability
table (the share of records whose haplotype changes between values), run the final
analysis and document the choice.

### Tutorial 3: Labels and Label Conflicts
**File**: `03_comparative_analysis_quality_control.ipynb`

Member labels and majority fractions per haplotype, label conflicts with their evidence,
depositor-note keyword hits, and a record-by-record crosswalk between a species-level and a
family-level run.

### Tutorial 4: Custom Shapefiles for Freshwater and Terrestrial Taxa
**File**: `04_custom_shapefiles_freshwater.ipynb`
**Additional requirement**: a shapefile such as HydroBASINS (https://www.hydrosheds.org/)

Inspect a shapefile, run with `--custom-shp`, and read records per haplotype and region.

### Tutorial 5: Population Genetics Exports
**File**: `05_population_genetics_export.ipynb`

Export with `--export-format`, see which records each format includes, and compute
haplotype diversity per region from `records.csv`.

### Tutorial 6: Plot Customization for Publication
**File**: `06_plot_customization_tutorial.ipynb`
**Requires**: a run with `--plot-kit`

Edit `plots/plot_config.yaml` (colors, sizes, projections, filters by haplotype label) and
regenerate the map, bar charts and tree.

## Recommended Learning Path

- **Start here**: Tutorial 1, then Tutorial 2
- **Labels and quality control**: Tutorial 3
- **Freshwater or terrestrial taxa**: Tutorial 4
- **Population genetics**: Tutorial 5
- **Publication figures**: Tutorial 6

## Tips

- Use a new output folder for each run; a run warns if its folder already holds a run
- Start with each run folder's `README.md`: it defines every file and column
- Clear outputs before committing notebooks (Kernel → Restart & Clear Output)

## Additional Resources

- **Main README**: `../README.md`
- **Migration from 1.0.0**: `../docs/guides/MIGRATION_v1_to_v2.md`
- **API Reference**: `../docs/guides/API_REFERENCE.md`
- **Parameter Sweep Guide**: `../docs/guides/PARAMETER_SWEEP_GUIDE.md`
- **PopGen Export Guide**: `../docs/guides/POPGEN_EXPORT_GUIDE.md`
- **Custom Shapefiles**: `../docs/guides/CUSTOM_SHAPEFILE_IMPLEMENTATION.md`

### External Software
- **Arlequin**: http://cmpg.unibe.ch/software/arlequin35/
- **PopART**: http://popart.otago.ac.nz/
- **DnaSP**: http://www.ub.edu/dnasp/
- **QGIS**: https://qgis.org/

## Support

- **GitHub Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
- **Email**: symbioseas@outlook.com

**BOLDGenotyper version**: 2.0.0
