# Changelog

All notable changes to BOLDGenotyper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025

### Added
- **Custom Shapefile Support**: `--custom-shp`, `--shp-field`, `--geo-category` for geographic analysis with any polygon shapefile (freshwater basins, ecoregions, watersheds, biomes, etc.)
- **Metadata Analysis Module**: Automated haplotype-metadata association analysis with chi-square tests, coverage assessment, and temporal emergence timelines (`metadata_analysis.py`)
- **Population Genetics Export**: Export to Arlequin (.arp), PopART/NEXUS, DnaSP (.fas), and generic CSV/FASTA via `--export-format`
- **Outgroup Rooting**: Three rooting modes — `--phylo-outgroup-fasta`, `--phylo-outgroup-label`, `--phylo-outgroup-taxon`
- **Haplotype Query Tool**: `boldgenotyper-query` command for assigning new sequences to existing haplotypes
- **Parameter Sweep Tool**: `boldgenotyper-sweep` command for optimizing singleton filtering thresholds
- **Plot Customization**: YAML-configured plot regeneration from exported data (`plot_export.py`, `plot_regeneration.py`)
- **Species-Faceted Analysis**: Per-species haplotype distribution maps and bar charts
- **MSA Visualization**: Phylogeny-ordered multiple sequence alignment plots (`msa_visualization.py`)
- **Divergence Analysis**: Pairwise divergence matrices, barcoding gap analysis, within/between-species divergence
- **ESV Approach**: Exact Sequence Variant haplotype discovery (replaces clustering-based approach)
- **COI Validation**: Orientation normalization and ORF validation with configurable genetic codes
- **Dynamic QC**: Median-based adaptive quality control filtering
- **Singleton Error Filtering**: Two-stage singleton quality control with configurable thresholds

### Changed
- Geographic analysis generalized from marine-only to any organism via custom shapefiles
- Haplotype discovery uses ESV approach instead of clustering thresholds
- Version bumped to 1.0.0 for publication release
- Python requirement raised to >=3.9 (3.8 is EOL)
- FastTree is the only supported tree-building algorithm

### Removed
- PhyML and RAxML tree-building options (only FastTree is used)
- `genotype_assignment.py` (replaced by `haplotype_assignment.py`)
- `cluster_diagnostics.py` (functionality integrated elsewhere)
- `comparative_analysis.py` (standalone script at `scripts/compare_analyses.py`)
- Dead configuration classes: `DelimitationConfig`, `OutlierConfig`, `DistanceConfig`

## [0.1.0] - 2025-01-20

### Added
- Core pipeline: BOLD TSV to annotated results
- Sequence dereplication with MAFFT alignment
- Edit distance-based genotype assignment with CIGAR string parsing
- GOaS ocean basin assignment for marine organisms
- Phylogenetic tree building with FastTree
- Geographic distribution maps with Cartopy
- Interactive HTML reports with Plotly.js
- Publication-ready visualizations (PDF/SVG)
- Hierarchical configuration system with dataclasses
- Comprehensive CLI with automatic organism name detection

---

## Version Numbering

BOLDGenotyper uses [Semantic Versioning](https://semver.org/):
- **MAJOR** version for incompatible API changes
- **MINOR** version for new functionality in a backwards compatible manner
- **PATCH** version for backwards compatible bug fixes

## Release Process

1. Update version number in `setup.py`, `boldgenotyper/__init__.py`, and `CITATION.cff`
2. Update this CHANGELOG.md with release date and changes
3. Create git tag: `git tag -a v1.0.0 -m "Release version 1.0.0"`
4. Push tag: `git push origin v1.0.0`

For detailed information about any release, see the corresponding GitHub release notes at:
https://github.com/SymbioSeas/BOLDGenotyper/releases
