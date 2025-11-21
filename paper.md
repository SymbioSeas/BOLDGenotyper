---
title: 'BOLDGenotyper: Automated COI genotyping and biogeographic analysis from barcode data'
tags:
  - Python
  - DNA barcoding
  - COI
  - phylogeography
  - population genetics
  - marine biology
  - genotyping
  - biogeography
authors:
  - name: Stephanie E. Smith
    orcid: 0000-0000-0000-0000
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Biology, University of North Carolina at Chapel Hill, Chapel Hill, NC 27599, USA
    index: 1
date: 20 January 2025
bibliography: paper.bib
---

# Summary

BOLDGenotyper is a comprehensive Python pipeline that automates the identification and biogeographic analysis of mitochondrial cytochrome c oxidase subunit I (COI) genotypes from the Barcode of Life Database (BOLD) [@Ratnasingham2007]. The software transforms raw BOLD sequence data into publication-ready analyses through seven integrated phases: data quality control, sequence dereplication, consensus genotype generation, sample-to-genotype assignment, taxonomic validation, phylogenetic reconstruction, and geographic visualization. BOLDGenotyper addresses a critical gap in molecular ecology by providing researchers with a standardized, reproducible workflow for analyzing intraspecific genetic variation and biogeographic patterns from the world's largest DNA barcode repository.

# Statement of need

DNA barcoding using the COI gene has become the gold standard for species identification, with the BOLD database now containing over 10 million sequences from 500,000+ species [@Ratnasingham2007; @Hebert2003]. While this wealth of data has enabled biodiversity studies and taxonomic research, the analysis of intraspecific genetic variation and population structure from BOLD data remains challenging. Researchers typically face several obstacles: (1) manual processing of large sequence datasets is time-consuming and error-prone, (2) no standardized workflow exists for identifying genotypes within species, (3) geographic coordinate quality varies widely in BOLD data, requiring careful filtering, and (4) integrating sequence, taxonomic, and geographic data for biogeographic analysis requires custom scripting.

Existing tools for COI analysis focus primarily on species delimitation [@Zhang2013; @Pons2006] or phylogenetic inference [@Stamatakis2014; @Price2010], but lack integrated workflows for intraspecific genotyping and biogeographic analysis. Population genetics software such as Arlequin [@Excoffier2010] and DnaSP [@Rozas2017] require manual sequence preparation and do not interface with BOLD data directly. Geographic analysis tools like QGIS and ArcGIS are general-purpose and not designed for molecular ecology workflows.

BOLDGenotyper fills this gap by providing an end-to-end automated pipeline specifically designed for BOLD data. The software implements rigorous quality control (coordinate filtering, sequence length thresholds, N-content filtering), sophisticated genotype assignment (edit distance-based with tie detection and confidence scoring), and publication-ready visualizations (interactive HTML reports, geographic distribution maps, phylogenetic trees). By automating these steps, BOLDGenotyper enables researchers to focus on biological interpretation rather than data wrangling.

The software has been designed with marine phylogeography in mind but is applicable to any organism with COI barcode data. A key innovation is the integration with the Global Oceans and Seas (GOaS) shapefile dataset [@MarineRegions2021], enabling automatic assignment of samples to ocean basins for biogeographic analyses. This feature is particularly valuable for studying marine species with wide distributions, where ocean basin boundaries may represent significant phylogeographic barriers.

# Design and functionality

BOLDGenotyper is implemented as a Python package with a unified command-line interface and is organized into eight core modules:

1. **Data ingestion and quality control** (`metadata.py`): Parses BOLD TSV files, validates required fields, and filters samples based on coordinate quality to exclude imprecise country-level centroids that could confound biogeographic analyses.

2. **Sequence dereplication** (`dereplication.py`): Aligns sequences with MAFFT [@Katoh2013], trims alignments with trimAl [@CapellaGutierrez2009], computes pairwise genetic distances, performs hierarchical clustering, and generates consensus sequences representing unique genotypes.

3. **Genotype assignment** (`genotype_assignment.py`): Assigns individual samples to consensus genotypes using edit distance calculations, implements tie detection for ambiguous assignments, and provides detailed diagnostic outputs for quality assessment.

4. **Taxonomic analysis** (`metadata.py`, `reports.py`): Aggregates species names within genotypes, identifies taxonomy conflicts (genotypes spanning multiple species), and generates species composition tables.

5. **Geographic analysis** (`geographic.py`): Integrates with GOaS shapefile data for ocean basin assignment, handles coordinate transformations, and provides land detection with configurable handling policies.

6. **Phylogenetic reconstruction** (`phylogenetics.py`): Optionally constructs maximum-likelihood trees using FastTree [@Price2010] with the GTR+Gamma model, generates tree visualizations, and exports Newick format for downstream analyses.

7. **Visualization** (`visualization.py`): Creates publication-quality figures including identity distribution histograms, geographic distribution maps with sample point sizing, ocean basin bar charts showing relative and total abundances, and faceted visualizations for multi-species datasets.

8. **Interactive reporting** (`reports.py`): Generates comprehensive HTML reports with interactive Plotly.js visualizations, filterable tables, complete methods sections suitable for publication, and export functionality for all plots and data.

The pipeline implements several algorithmic innovations for handling BOLD data:

- **Three-stage length filtering** prevents spurious genotypes from sequence fragments by filtering before alignment (≥400 bp), after trimming (≥300 bp ungapped), and after consensus generation (≥75% of median length).

- **Edit distance-based assignment** with CIGAR string parsing provides more accurate identity calculations than simple percent identity, properly accounting for insertions and deletions.

- **Tie detection and confidence scoring** flags ambiguous assignments where multiple genotypes match with similar identity, enabling researchers to manually review borderline cases.

- **Coordinate quality assessment** implements conservative filtering that excludes country centroids for nations bordering multiple ocean basins, prioritizing geographic accuracy over sample size.

# Example usage and validation

BOLDGenotyper has been applied to analyze population structure in the scalloped hammerhead shark (*Sphyrna lewini*), a cosmopolitan species with uncertain phylogeography. Analysis of 937 COI sequences from BOLD identified 26 distinct genotypes and revealed ocean basin-scale genetic partitioning, with distinct genotypes dominating different ocean basins. This case study demonstrates the software's utility for detecting cryptic population structure in widely distributed marine species and is being prepared for publication as a companion manuscript to this software paper.

A complete example analysis of *Sphyrna lewini* is included in the repository (`data/Sphyrna_lewini/`), providing users with reference outputs for testing and learning. The example includes the interactive HTML report, all visualizations, phylogenetic trees, and diagnostic files, demonstrating the full range of software capabilities.

# Comparison to related software

BOLDGenotyper differs from existing tools in several key aspects:

| Feature | BOLDGenotyper | BOLD Portal | TCS/PopART | Arlequin | QGIS |
|---------|---------------|-------------|------------|----------|------|
| Direct BOLD integration | ✓ | ✓ | ✗ | ✗ | ✗ |
| Automated genotyping | ✓ | ✗ | Manual | Manual | N/A |
| Geographic analysis | ✓ | Limited | ✗ | ✗ | Manual |
| Ocean basin assignment | ✓ | ✗ | ✗ | ✗ | Manual |
| Interactive reports | ✓ | ✗ | ✗ | ✗ | ✗ |
| Phylogenetic trees | ✓ | ✗ | ✗ | ✓ | ✗ |
| Publication-ready output | ✓ | Partial | Partial | Partial | ✓ |
| Reproducible workflow | ✓ | ✗ | Partial | Partial | Partial |

While the BOLD portal [@Ratnasingham2007] provides species identification and basic sequence analysis, it does not support automated genotyping or detailed biogeographic analyses. Population genetics tools like TCS [@Clement2000], PopART [@Leigh2015], and Arlequin [@Excoffier2010] require extensive manual preparation of BOLD data and do not integrate geographic analysis. BOLDGenotyper bridges these gaps by providing an end-to-end automated workflow.

# Software quality and testing

BOLDGenotyper includes 215 automated tests covering unit testing of individual functions, integration testing of multi-step workflows, and error handling for edge cases. The test suite achieves >80% code coverage for core modules and includes tests for external tool failures with graceful degradation. Comprehensive documentation includes function-level docstrings in NumPy style, a detailed README with installation instructions and usage examples, and an extensive troubleshooting guide.

The software follows best practices for scientific Python software: type hints for function signatures, modular architecture with clear separation of concerns, comprehensive logging for reproducibility, and structured configuration management using dataclasses. All code adheres to PEP 8 style guidelines enforced by Black and flake8.

# Availability and implementation

BOLDGenotyper is implemented in Python 3.8+ and distributed under the MIT license. The source code is available on GitHub at https://github.com/SymbioSeas/BOLDGenotyper and will be archived on Zenodo upon publication. The software can be installed via conda or pip and requires external tools (MAFFT, trimAl, FastTree) available through Bioconda. Comprehensive documentation, example data, and usage tutorials are provided in the GitHub repository.

# Acknowledgments

The author thanks Dr. Chelsea Black for collaboration on the *Sphyrna lewini* case study, the BOLD Systems team for maintaining the barcode database, and the BioConda community for packaging bioinformatics tools.

# References
