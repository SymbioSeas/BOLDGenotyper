# BOLDGenotyper Tutorial Notebooks

This directory contains interactive Jupyter notebooks that guide you through common BOLDGenotyper workflows.

## Overview

The tutorials are designed to be completed in sequence, with each building on concepts from previous notebooks. However, they can also be used as standalone references for specific workflows.

## Prerequisites

- BOLDGenotyper installed (see main README.md)
- Jupyter Notebook or JupyterLab
- Example datasets (included in `data/` directory)
- Basic understanding of DNA barcoding and COI sequences

## Installation

Install Jupyter in your BOLDGenotyper conda environment:

```bash
conda activate boldgenotyper
conda install -c conda-forge jupyter
```

Or with pip:

```bash
pip install jupyter
```

## Launching the Tutorials

From the repository root directory:

```bash
jupyter notebook notebooks/
```

Or with JupyterLab:

```bash
jupyter lab notebooks/
```

## Tutorial Series

### Tutorial 1: Basic Genotyping Workflow
**File**: `01_basic_genotyping_workflow.ipynb`
**Duration**: 30-45 minutes
**Dataset**: Sphyrna lewini (~600 samples)

**Topics covered**:
- Running a basic genotyping analysis
- Interpreting genotype assignments
- Exploring geographic distributions
- Understanding output files
- Navigating the HTML interactive report

**Learning objectives**:
- Understand the complete BOLDGenotyper workflow
- Interpret genotype assignment quality metrics
- Analyze ocean basin-scale genetic structure
- Navigate output directory structure

---

### Tutorial 2: Parameter Sweep for Threshold Optimization
**File**: `02_parameter_sweep_optimization.ipynb`
**Duration**: 45-60 minutes
**Dataset**: Sphyrna lewini (~600 samples)

**Topics covered**:
- Running parameter sweep analysis
- Interpreting elbow plots and stability metrics
- Detecting optimal clustering thresholds
- Comparing default vs optimized results
- Generating publication methods text

**Learning objectives**:
- Understand threshold optimization principles
- Use data-driven approaches for parameter selection
- Interpret multiple evaluation criteria
- Justify threshold choices for publication

---

### Tutorial 3: Comparative Analysis for Quality Control
**File**: `03_comparative_analysis_quality_control.ipynb`
**Duration**: 45-60 minutes
**Datasets**: Sphyrna lewini + Sphyrnidae (~600 + ~1000 samples)

**Topics covered**:
- Multi-level taxonomic comparison
- Detecting contamination and mislabeling
- Identifying cryptic diversity
- Sample reassignment recommendations
- Quantifying contamination rates

**Learning objectives**:
- Perform quality control using comparative analysis
- Interpret genotype crosswalk tables
- Make data-driven decisions about flagged samples
- Document quality control for publication

---

### Tutorial 4: Custom Shapefiles for Freshwater Organisms
**File**: `04_custom_shapefiles_freshwater.ipynb`
**Duration**: 60-90 minutes
**Dataset**: Salmonidae (~9900 samples)
**Additional requirement**: HydroBASINS shapefile (download from https://www.hydrosheds.org/)

**Topics covered**:
- Inspecting shapefile structure and attributes
- Running enrichment with custom shapefiles
- Analyzing genotype distributions by drainage basin
- Identifying basin-specific genotypes
- Extending to terrestrial organisms

**Learning objectives**:
- Use BOLDGenotyper for non-marine organisms
- Work with custom geographic classifications
- Analyze freshwater phylogeographic patterns
- Apply to conservation and biogeography questions

---

### Tutorial 5: Population Genetics Export and Integration
**File**: `05_population_genetics_export.ipynb`
**Duration**: 45-60 minutes
**Dataset**: Sphyrna lewini (~600 samples)

**Topics covered**:
- Exporting to Arlequin, PopART, and DnaSP formats
- Calculating genetic diversity metrics
- Preparing data for haplotype networks
- Setting up AMOVA analyses
- Integrating with downstream software

**Learning objectives**:
- Export BOLDGenotyper results to population genetics formats
- Perform basic diversity calculations
- Construct analysis workflows with external tools
- Generate publication-ready analyses

---

### Tutorial 6: Plot Customization for Publication
**File**: `06_plot_customization_tutorial.ipynb`
**Duration**: 30-45 minutes
**Dataset**: Any completed BOLDGenotyper analysis

**Topics covered**:
- Understanding plot_config.yaml structure
- Common YAML syntax errors and solutions
- Customizing colors, sizes, and formats
- Modifying map projections and styles
- Filtering genotypes in visualizations
- Regenerating plots with Python scripts

**Learning objectives**:
- Customize publication-ready figures without re-running analysis
- Troubleshoot common configuration errors
- Apply journal-specific formatting requirements
- Create custom color schemes and layouts
- Generate multiple plot variants efficiently

---

## Recommended Learning Path

### For Beginners
Start with tutorials in order:
1. Tutorial 1 (required foundation)
2. Tutorial 2 (parameter optimization)
3. Tutorial 6 (customizing plots)
4. Tutorial 5 (basic population genetics)

### For Advanced Users
Focus on specific needs:
- **Quality control**: Tutorials 1 + 3
- **Threshold optimization**: Tutorials 1 + 2
- **Freshwater/terrestrial organisms**: Tutorials 1 + 4
- **Population genetics**: Tutorials 1 + 5
- **Publication figures**: Tutorials 1 + 6

### For Publication Preparation
Complete all tutorials:
1. Tutorial 1 (basic workflow)
2. Tutorial 2 (justify threshold choice)
3. Tutorial 3 (demonstrate quality control)
4. Tutorial 6 (customize publication figures)
5. Tutorial 4 or 5 (organism-specific analysis)

---

## Example Datasets

All tutorials use datasets included in the `data/` directory:

| Dataset | Organism | Samples | Purpose |
|---------|----------|---------|---------|
| Sphyrna_lewini | Scalloped hammerhead shark | ~600 | Species-level analysis, quick examples |
| Sphyrnidae | Hammerhead shark family | ~1000 | Family-level analysis, comparative QC |
| Salmonidae | Salmon and trout | ~9900 | Large dataset, freshwater example |
| Cyprinodontidae | Pupfish | ~1200 | Freshwater example (alternative) |

See `data/README.md` for detailed dataset descriptions.

---

## Tips for Success

### General
- **Complete tutorials in order** - Each builds on previous concepts
- **Run all code cells** - Don't skip cells, even if they seem simple
- **Examine outputs carefully** - Understanding output structure is crucial
- **Modify parameters** - Experiment with different settings
- **Save your work** - Notebooks autosave, but manually save important results

### Technical
- **Check file paths** - Tutorial assumes you're running from repository root
- **Verify data availability** - Ensure example datasets are present
- **Monitor runtime** - Large datasets may take several minutes
- **Clear outputs** - Use "Cell > All Output > Clear" to reset notebook
- **Restart kernel** - If errors persist, restart kernel and re-run

### Troubleshooting
- **Import errors**: Verify BOLDGenotyper is installed: `conda list boldgenotyper`
- **File not found**: Check you're running from correct directory
- **Memory errors**: Reduce dataset size or use fewer threads
- **Slow execution**: Close other applications, use fewer threads

---

## Additional Resources

### Documentation
- **Main README**: `../README.md` - Complete BOLDGenotyper documentation
- **API Reference**: `../API_REFERENCE.md` - Programmatic usage
- **Parameter Sweep Guide**: `../PARAMETER_SWEEP_GUIDE.md` - Threshold optimization
- **Comparative Analysis Guide**: `../COMPARATIVE_ANALYSIS_GUIDE.md` - Quality control
- **Custom Shapefiles Guide**: `../CUSTOM_SHAPEFILES_GUIDE.md` - Non-marine organisms
- **PopGen Export Guide**: `../POPGEN_EXPORT_GUIDE.md` - Population genetics
- **Plot Customization**: See Tutorial 6 and main README section "Customizing Plots for Publication"

### Example Data
- **Dataset descriptions**: `../data/README.md`
- **Expected runtimes**: `../data/README.md`
- **Download instructions**: `../README.md` section "Obtaining Input Data from BOLD"

### External Software
- **Arlequin**: http://cmpg.unibe.ch/software/arlequin35/
- **PopART**: http://popart.otago.ac.nz/
- **DnaSP**: http://www.ub.edu/dnasp/
- **QGIS**: https://qgis.org/ (for shapefile manipulation)

---

## Contributing

Found a problem or have a suggestion for improving the tutorials?

1. **Report issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
2. **Suggest improvements**: Open an issue with the "documentation" label
3. **Submit corrections**: Fork, edit, and submit a pull request

---

## Citation

If you use these tutorials in your research or teaching, please cite:

```
Smith, S. (2025). BOLDGenotyper: Automated COI sequence genotyping and
biogeographic analysis. GitHub repository:
https://github.com/SymbioSeas/BOLDGenotyper
```

---

## Support

For questions or assistance:
- **GitHub Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
- **Email**: steph.smith@unc.edu
- **Documentation**: See main README.md

---

## Version Information

- **Tutorial Series Version**: 1.0
- **BOLDGenotyper Version**: 0.1.0
- **Last Updated**: 2025-11-25

---

## License

These tutorials are released under the same license as BOLDGenotyper. See LICENSE file for details.

---

**Happy genotyping!**
