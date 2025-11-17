# Changelog

All notable changes to BOLDGenotyper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive documentation improvement plan for publication readiness
- CITATION.cff for standardized software citation
- MIT LICENSE file
- This CHANGELOG.md file

### Changed
- Documentation structure reorganized for publication

## [0.1.0] - 2025-11-17

### Added
- **Core Pipeline**: Complete end-to-end analysis pipeline from BOLD TSV to annotated results
- **Sequence Dereplication**:
  - MAFFT-based multiple sequence alignment
  - trimAl-based alignment trimming
  - Hierarchical clustering with configurable thresholds (default: 99% identity)
  - Consensus sequence generation with majority-rule voting
- **Genotype Assignment**:
  - CIGAR string-based identity calculation
  - Support for both "target_based" and "classic" identity metrics
  - Parallel processing support for large datasets
  - Comprehensive diagnostics output
- **Geographic Analysis**:
  - GOaS (Global Oceans and Seas) shapefile integration
  - Ocean basin assignment for marine organisms
  - Coordinate quality filtering (precise vs country centroid)
  - Land detection and snap-to-coast functionality
  - `--no-geo` flag for non-marine organisms or when GOaS unavailable
- **Phylogenetic Analysis** (Optional):
  - Tree building with FastTree, RAxML, or IQ-TREE
  - Tree visualization with colored tips by genotype
  - Newick format output
- **Visualization**:
  - Geographic distribution maps with Cartopy
  - Genotype abundance plots
  - Identity distribution histograms
  - Phylogenetic trees with custom styling
  - Multiple output formats (PNG, PDF, SVG)
- **Configuration System**:
  - Hierarchical configuration via dataclasses
  - YAML/JSON config file support
  - Command-line parameter overrides
- **Command-Line Interface**:
  - Single-command pipeline execution
  - Automatic organism name detection from filename
  - Flexible parameter customization
  - Comprehensive logging with configurable levels
- **Documentation**:
  - Comprehensive README with installation, usage, and troubleshooting
  - QUICKSTART guide for new users
  - CLI testing guide
  - GOaS setup documentation
  - Analysis notes and technical documentation
- **Quality Control**:
  - Coordinate precision filtering
  - Sequence length validation
  - Assignment confidence tracking
  - Comprehensive error handling

### Features by Module
- **boldgenotyper.cli**: Complete CLI interface with pipeline orchestration
- **boldgenotyper.config**: Hierarchical configuration management
- **boldgenotyper.dereplication**: Sequence clustering and consensus generation
- **boldgenotyper.genotype_assignment**: Sample-to-genotype mapping with CIGAR
- **boldgenotyper.geographic**: Spatial analysis and ocean basin assignment
- **boldgenotyper.metadata**: BOLD TSV parsing with flexible column mapping
- **boldgenotyper.phylogenetics**: Tree building and visualization
- **boldgenotyper.reports**: Summary statistics and report generation
- **boldgenotyper.utils**: Helper functions and utilities
- **boldgenotyper.visualization**: Publication-ready figure generation
- **boldgenotyper.goas_downloader**: Automated GOaS shapefile download

### Dependencies
- Python ≥3.8
- biopython ≥1.79
- pandas ≥1.3.0
- scipy ≥1.7.0
- numpy ≥1.21.0
- matplotlib ≥3.4.0
- seaborn ≥0.11.0
- Optional: geopandas, cartopy (for geographic analysis)
- External tools: MAFFT, trimAl, FastTree (optional)

### Known Limitations
- Geographic analysis currently optimized for marine organisms only
- GOaS shapefile not included in repository (must be downloaded separately)
- Large datasets (>10,000 sequences) may encounter memory limitations
- Phylogenetics module has some stub functions for planned features
- Windows support requires WSL (Windows Subsystem for Linux)

### Performance
- Small datasets (<500 samples): 2-5 minutes
- Medium datasets (500-2,000 samples): 5-15 minutes
- Large datasets (>2,000 samples): 15-60 minutes
- Parallel processing supported for genotype assignment

### Testing
- Comprehensive test suite for core modules
- Example datasets included for validation
- Integration tests for full pipeline

## [0.0.1] - Development Versions

### 2024-11-XX to 2025-11-XX
- Initial development and prototyping
- Reference scripts implementation
- Module refactoring and package structure
- Documentation iterations
- Testing and validation

---

## Version Numbering

BOLDGenotyper uses [Semantic Versioning](https://semver.org/):
- **MAJOR** version for incompatible API changes
- **MINOR** version for new functionality in a backwards compatible manner
- **PATCH** version for backwards compatible bug fixes

## Release Process

1. Update version number in `setup.py` and `boldgenotyper/__init__.py`
2. Update this CHANGELOG.md with release date and changes
3. Create git tag: `git tag -a v0.1.0 -m "Release version 0.1.0"`
4. Push tag: `git push origin v0.1.0`
5. GitHub Actions will automatically create release and publish to PyPI (when configured)
6. Zenodo will automatically create DOI for the release

## Future Releases

### [0.2.0] - Planned Features
- Terrestrial and freshwater organism support
- Additional geographic analysis modules
- Enhanced phylogenetic analysis
- Performance optimizations for very large datasets
- Docker container support
- Web interface (optional)

### [1.0.0] - Planned for Publication
- Complete API stability
- Comprehensive documentation
- Published peer-reviewed software paper
- Full test coverage
- Performance benchmarks
- Production-ready release

---

**Note**: This project is under active development. Features and API may change before version 1.0.0.

For detailed information about any release, see the corresponding GitHub release notes at:
https://github.com/SymbioSeas/BOLDGenotyper/releases
