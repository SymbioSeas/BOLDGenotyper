# Testing New Plot Regeneration Features

This guide walks you through testing the new Python-based plot regeneration system that replaces the previous R-based approach.

## What Changed

### Before (R-based)
- Required separate R installation and complex package management
- Used `install_plot_deps.sh` to install R packages (sf, ggtree, etc.)
- Generated R scripts for plot customization
- Cross-platform issues (especially on Windows)

### After (Python-based)
- All dependencies managed through conda
- Everything installed with `conda env create -f environment.yml`
- Python scripts for plot customization
- Works on Mac, Linux, and Windows

## Prerequisites

Before testing, make sure you have:
1. conda or miniconda installed
2. A clean environment (or remove old environment: `conda env remove -n boldgenotyper`)

## Test 1: Clean Installation

Test that the new installation process works from scratch.

### Mac/Linux

```bash
# Navigate to repository
cd /Users/stesmith/Documents/depredation/boldgenotyper

# Option A: Use installation script
bash install.sh

# Option B: Manual installation
conda env create -f environment.yml
conda activate boldgenotyper
pip install -e .

# Verify installation
boldgenotyper --version
mafft --version
trimal --version
fasttree 2>&1 | head -1

# Check Python packages
python -c "import geopandas, cartopy, yaml, matplotlib; print('All packages OK')"
```

### Expected Results
- Environment creates without errors
- All external tools (mafft, trimal, fasttree) are accessible
- All Python packages import successfully
- `boldgenotyper --version` prints version number

### Troubleshooting
- **conda slow?** Try using mamba: `conda install mamba -c conda-forge`, then `mamba env create -f environment.yml`
- **Cartopy installation issues?** This is usually a conda solver issue. Try creating environment with `--solver=libmamba`

## Test 2: Run Pipeline with Test Data

Test that the pipeline generates the new Python regeneration scripts.

```bash
# Activate environment
conda activate boldgenotyper

# Run on test data
boldgenotyper \
  data/Sphyrnidae_test/Sphyrnidae_tsv_data.tsv \
  --organism Sphyrnidae_test \
  --output data/Sphyrnidae_test_output \
  --threshold 0.01 \
  --build-tree \
  --verbose

# Or use existing test data location
boldgenotyper \
  data/Sphyrnidae_test/Sphyrnidae_tsv_data.tsv \
  --organism Sphyrnidae \
  --output data/Sphyrnidae_test \
  --threshold 0.01 \
  --build-tree \
  --verbose
```

### Expected Results

The pipeline should complete successfully and create a `plots/` directory:

```
data/Sphyrnidae_test/plots/
├── plot_config.yaml
├── README.md
├── data/
│   ├── distribution_map.csv
│   ├── distribution_bar_relative.csv
│   ├── distribution_bar_absolute.csv
│   ├── identity_distribution.csv
│   └── genotype_colors.csv
└── scripts/
    ├── regenerate_all.py
    ├── regenerate_map.py
    ├── regenerate_bars.py
    └── regenerate_identity.py
```

### Check Files

```bash
# Check that Python scripts were created (not R scripts)
ls data/Sphyrnidae_test/plots/scripts/

# Should see:
# regenerate_all.py
# regenerate_map.py
# regenerate_bars.py
# regenerate_identity.py

# Should NOT see:
# regenerate_all.sh
# regenerate_map.R
# regenerate_bars.R
# regenerate_identity.R
# install_plot_deps.sh
```

### Verify Script Contents

```bash
# Check that Python scripts import from boldgenotyper
head -20 data/Sphyrnidae_test/plots/scripts/regenerate_map.py

# Should contain:
# #!/usr/bin/env python3
# from boldgenotyper.plot_regeneration import regenerate_distribution_map
```

## Test 3: Plot Customization

Test that plot regeneration works with custom settings.

### Step 1: Edit Configuration

```bash
cd data/Sphyrnidae_test/plots

# Open plot_config.yaml in your favorite editor
nano plot_config.yaml  # or vim, emacs, VS Code, etc.
```

### Step 2: Modify Settings

Try changing these settings:

```yaml
general:
  dpi: 600  # Increase resolution
  width_inches: 14
  height_inches: 10

map:
  projection: "mollweide"  # Change from robinson
  center_longitude: -180   # Center on Pacific
  point_alpha: 0.9

bars:
  orientation: "horizontal"  # Change from vertical

identity:
  binwidth: 0.3  # Make histogram bins smaller
  show_density: true
```

### Step 3: Regenerate Plots

```bash
# From the plots/ directory
python scripts/regenerate_all.py

# Or regenerate individual plots
python scripts/regenerate_map.py
python scripts/regenerate_bars.py
python scripts/regenerate_identity.py
```

### Expected Results

```
Regenerating all plots...

  1. Regenerating distribution map...
  ✓ Wrote data/Sphyrnidae_test/visualization/distribution_map_custom.pdf
  ✓ Wrote data/Sphyrnidae_test/visualization/distribution_map_custom.png

  2. Regenerating bar charts...
  ✓ Wrote data/Sphyrnidae_test/visualization/distribution_bar_relative_custom.pdf
  ✓ Wrote data/Sphyrnidae_test/visualization/distribution_bar_relative_custom.png
  ✓ Wrote data/Sphyrnidae_test/visualization/distribution_bar_absolute_custom.pdf
  ✓ Wrote data/Sphyrnidae_test/visualization/distribution_bar_absolute_custom.png

  3. Regenerating identity distribution...
  ✓ Wrote data/Sphyrnidae_test/visualization/identity_distribution_custom.pdf
  ✓ Wrote data/Sphyrnidae_test/visualization/identity_distribution_custom.png

✓ All plots regenerated successfully! (6 files)
  Check data/Sphyrnidae_test/visualization for updated plots
```

### Verify Output

```bash
# Check that custom plots were created
ls data/Sphyrnidae_test/visualization/*_custom.*

# View plots (Mac)
open data/Sphyrnidae_test/visualization/distribution_map_custom.pdf
open data/Sphyrnidae_test/visualization/distribution_bar_relative_custom.pdf
open data/Sphyrnidae_test/visualization/identity_distribution_custom.pdf

# View plots (Linux with PDF viewer)
evince data/Sphyrnidae_test/visualization/distribution_map_custom.pdf &
```

### Verify Changes Applied

- Map should use Mollweide projection (oval shape)
- Bar charts should be horizontal instead of vertical
- Identity histogram bins should be narrower
- All plots should have higher resolution (600 DPI)

## Test 4: Error Handling

Test that the scripts provide helpful error messages.

### Test Missing Data

```bash
cd data/Sphyrnidae_test/plots

# Temporarily move data directory
mv data data_backup

# Try to regenerate
python scripts/regenerate_all.py

# Should see clear error message:
# Error: distribution_map.csv not found in data/

# Restore data
mv data_backup data
```

### Test Invalid Configuration

```bash
# Edit plot_config.yaml and add invalid projection
nano plot_config.yaml

# Change projection to something invalid:
# map:
#   projection: "invalid_projection"

# Try to regenerate
python scripts/regenerate_map.py

# Should handle gracefully and fall back or show error

# Fix configuration
nano plot_config.yaml
```

## Test 5: Cross-Platform Compatibility (If Available)

If you have access to multiple platforms, test installation on each:

### macOS (Intel or ARM)
```bash
bash install.sh
# Run tests above
```

### Linux
```bash
bash install.sh
# Run tests above
```

### Windows (via WSL2 or native conda)
```bash
# If using WSL
bash install.sh

# If using Windows conda prompt
install.bat
```

## Test 6: Package Import

Test that the new plot_regeneration module can be imported.

```bash
conda activate boldgenotyper

python -c "from boldgenotyper import plot_regeneration; print('Import OK')"

python -c "from boldgenotyper.plot_regeneration import regenerate_all_plots; print(regenerate_all_plots.__doc__)"
```

### Expected Results
- No import errors
- Function docstring prints correctly

## Test 7: Integration with Existing Data

If you have existing output from a previous run, test regeneration on that data.

```bash
# Find an existing output directory
ls data/

# Navigate to plots directory
cd data/[your_organism]/plots

# Check if scripts exist and are Python (not R)
ls scripts/

# If they're R scripts, you can manually test by running the pipeline again
# to generate new Python scripts, or manually create them using plot_export module

# Test regeneration
python scripts/regenerate_all.py
```

## Common Issues and Solutions

### Issue: "ImportError: No module named 'geopandas'"

**Solution**: The conda environment wasn't activated correctly
```bash
conda activate boldgenotyper
python scripts/regenerate_all.py
```

### Issue: "FileNotFoundError: plot_config.yaml not found"

**Solution**: Make sure you're in the plots/ directory
```bash
cd data/[organism]/plots
python scripts/regenerate_all.py
```

### Issue: "ModuleNotFoundError: No module named 'boldgenotyper'"

**Solution**: The package isn't installed or scripts can't find it
```bash
# Reinstall package
conda activate boldgenotyper
pip install -e /path/to/boldgenotyper/repo
```

### Issue: Cartopy projection errors

**Solution**: Fall back to basic plotting (the scripts handle this automatically)
- Maps will still generate but without advanced projections
- Check that cartopy is properly installed: `python -c "import cartopy"`

## Success Criteria

✅ **Installation**
- Environment creates without errors
- All dependencies install via conda (no manual R setup needed)
- Scripts work out of the box

✅ **Pipeline Integration**
- Pipeline generates Python scripts (not R scripts)
- plot_config.yaml created correctly
- Data files exported to plots/data/

✅ **Plot Regeneration**
- All three plot types regenerate successfully
- Custom configurations apply correctly
- Output files created in correct locations

✅ **User Experience**
- Clear error messages when things go wrong
- README.md in plots/ directory provides guidance
- No need for R knowledge or installation

## Reporting Issues

If you encounter problems, please provide:

1. **Error message**: Copy the full error output
2. **Environment info**: Run `conda list` and share output
3. **OS and platform**: Mac (Intel/ARM), Linux, Windows
4. **Steps to reproduce**: Exact commands that caused the issue
5. **Files**: Share plot_config.yaml if regeneration fails

## Next Steps After Testing

Once testing is complete and everything works:

1. **Commit changes**:
   ```bash
   git add -A
   git commit -m "Replace R plot regeneration with Python implementation

   - Add plot_regeneration.py module with Python plotting functions
   - Update plot_export.py to generate Python scripts
   - Consolidate all dependencies in environment.yml
   - Remove install_plot_deps.sh (R dependencies no longer needed)
   - Add install.sh and install.bat for easy setup
   - Update README with simplified installation instructions
   - Add plot customization section to README"
   ```

2. **Update documentation** if you find any gaps

3. **Consider creating a release** if this is a major milestone

## Questions?

If you have questions about the new system or need help testing:
- Check the README.md customization section
- Look at plots/README.md after running pipeline
- Review the plot_regeneration.py module docstrings
