# BOLDGenotyper - Project Instructions for Claude Code

## Project Overview

BOLDGenotyper is a Python package for automated COI sequence genotyping and biogeographic analysis from BOLD (Barcode of Life Database) data. The package enables researchers to identify genetic partitioning patterns in any organism with publicly available COI barcode sequences.

### Scientific Context
This tool was developed to support the manuscript "Ocean basin-scale genetic partitioning in *Sphyrna lewini* revealed through COI sequence analysis" and is designed to be generalizable to any taxonomic group in BOLD.

## Core Functionality

1. **Data Input**: TSV file downloaded directly from BOLD for any organism
2. **Sequence Dereplication**: Identify unique COI genotypes using MAFFT alignment and hierarchical clustering
3. **Consensus Generation**: Create consensus sequences for each genotype cluster
4. **Metadata Integration**: Map samples to genotype groups and parse geographic coordinates
5. **Geographic Filtering**: Remove samples with ambiguous locations (country centroids)
6. **Visualization**: Generate distribution maps and relative abundance plots by ocean basin
7. **Phylogenetic Analysis**: MAFFT alignment and PhyML tree construction (optional)

## Technical Requirements

### Environment
- Conda environment name: `boldgenotyper`
- Cloned from existing `depredation` environment
- Python 3.8+

### Dependencies
- **Core Python**: biopython, pandas, scipy, numpy
- **Bioinformatics Tools**: MAFFT (v7+), trimAl, PhyML
- **Optional**: edlib (for faster sequence alignment)
- **Visualization**: matplotlib, cartopy/basemap, seaborn

### File Naming Conventions
Users provide organism-specific naming, but the pipeline should:
- Parse the TSV filename to extract organism name
- Maintain consistent naming throughout: `{organism}_consensus.fasta`, `{organism}_with_genotypes.tsv`, etc.
- Handle spaces and special characters in organism names gracefully

## Input File Structure

The input TSV from BOLD contains ~86 columns including:
- `processid`: Unique sample identifier (REQUIRED)
- `nuc`: COI nucleotide sequence (REQUIRED)
- `coord`: Geographic coordinates in format `[lat, lon]` (REQUIRED for geographic analysis)
- `country/ocean`: Location metadata
- `coord_source`: Indicates if coordinates are from centroid
- Other metadata: species, date, institution, etc.

### Coordinate Filtering Rules
**EXCLUDE samples where:**
- `coord` field is empty or missing
- `coord_source` contains "centroid" or similar indicators
- Only country-level location is provided
- Coordinates are `[0, 0]` or other obvious placeholders

**This filtering is critical for accurate ocean basin assignment.**

## Analysis Workflow

### 1. Sequence Dereplication (`msa_to_consensus.py`)
- Extract all COI sequences from TSV
- Create FASTA file with headers: `{organism}_{processid}.COI-5P`
- Run MAFFT alignment
- Trim alignment with trimAl (--automated1)
- Calculate pairwise sequence distances (ignoring gaps and Ns)
- Hierarchical clustering (default threshold: 0.01 = 99% identity)
- Generate consensus sequences using majority rule
- Output naming: `consensus_c{cluster_id}_n{sample_count}`

### 2. Genotype Assignment (`consensus_group_to_metadata.py`)
- Map each sample to its best-matching consensus sequence
- Use global edit distance (edlib if available, else Levenshtein)
- Require minimum identity threshold (default: 0.90)
- Add `consensus_group` column to metadata TSV
- Generate diagnostics showing identity scores

### 3. Geographic Filtering
- Parse `coord` field from bracketed format to lat/lon floats
- Apply coordinate filtering rules (see above)
- Assign samples to ocean basins based on coordinates
- Ocean basins: North Atlantic, South Atlantic, North Pacific, South Pacific, Indian Ocean, South China Seas, Southern Ocean

### 4. Visualization
Generate publication-ready figures:
- **Global distribution map**: Points colored by genotype, sized by sample count
- **Ocean basin abundance**: Stacked bar chart showing genotype proportions
- **Phylogenetic tree**: Maximum likelihood tree with bootstrap support (if requested)
- **Sequence alignment**: Pairwise comparison highlighting differences (if requested)

### 5. Phylogenetic Analysis (Optional)
- MAFFT alignment of consensus sequences
- PhyML with GTR model and 1000 bootstrap replicates
- Tree visualization with bootstrap values

## Code Organization

```
BOLDGenotyper/
├── boldgenotyper/
│   ├── __init__.py
│   ├── core.py              # Main pipeline orchestration
│   ├── dereplication.py     # Sequence clustering and consensus
│   ├── metadata.py          # TSV parsing and genotype assignment
│   ├── geographic.py        # Coordinate filtering and ocean basin assignment
│   ├── visualization.py     # All plotting functions
│   ├── phylogenetics.py     # Tree building and analysis
│   └── utils.py             # Helper functions
├── tests/
│   └── test_*.py            # Unit tests for each module
├── examples/
│   └── example_workflow.py  # Tutorial script
├── environment.yml          # Conda environment specification
├── setup.py                 # Package installation
├── README.md                # User documentation
├── CLAUDE.md                # This file
└── LICENSE

```

## Documentation Requirements

Each module must include:
1. **Docstrings**: Comprehensive function and class documentation
2. **Type hints**: All function parameters and returns
3. **Examples**: Usage examples in docstrings
4. **Assumptions**: Clearly state any data assumptions
5. **Edge cases**: Document how edge cases are handled

### Critical Documentation Points
- Explain why coordinate filtering is necessary
- Document the clustering threshold rationale (99% identity)
- Describe consensus sequence generation method
- Explain ocean basin boundary definitions
- Note dependencies on external tools (MAFFT, PhyML)

## User Experience Goals

### For Bioinformaticians
- Clear command-line interface
- Modular design for custom workflows
- Well-documented API for programmatic use

### For Non-Technical Users (Future)
- Web app with file upload
- Automated analysis with progress indicators
- Downloadable results package (figures + data)
- No command-line experience required

## Development Priorities

1. **Phase 1** (Current): Core pipeline with CLI
   - Robust file parsing
   - Sequence analysis workflow
   - Basic visualizations
   - Comprehensive documentation

2. **Phase 2**: Enhanced features
   - Additional visualization options
   - Statistical analyses
   - Multiple organism comparison
   - Performance optimization

3. **Phase 3**: Web application
   - Flask/Streamlit interface
   - Cloud-based analysis
   - Result sharing and export
   - Interactive visualizations

## GitHub Best Practices

- **Branching**: Use feature branches (`feature/add-visualization`)
- **Commits**: Clear, descriptive commit messages
- **Issues**: Track bugs and feature requests
- **Releases**: Semantic versioning (v1.0.0, v1.1.0, etc.)
- **CI/CD**: Automated testing with GitHub Actions
- **Documentation**: Keep README and examples up-to-date

## Testing Strategy

- **Unit tests**: Each function in isolation
- **Integration tests**: Full pipeline with test data
- **Edge cases**: Empty files, malformed coordinates, missing data
- **Performance**: Test with large datasets (1000+ samples)

## Example Usage (Planned)

```bash
# Basic usage
boldgenotyper analyze \
  --input Organism_name.tsv \
  --output results/ \
  --threads 4

# With phylogenetic analysis
boldgenotyper analyze \
  --input Organism_name.tsv \
  --output results/ \
  --phylogeny \
  --outgroup Outgroup_species.fasta \
  --threads 8

# Custom clustering threshold
boldgenotyper analyze \
  --input Organism_name.tsv \
  --output results/ \
  --cluster-threshold 0.005 \
  --min-identity 0.95
```

## Notes for Claude Code

- **File paths**: Use relative paths and `pathlib` for cross-platform compatibility
- **Error handling**: Graceful failures with informative error messages
- **Logging**: Use Python logging module for progress tracking
- **Memory**: Stream large files when possible, don't load entire datasets into memory
- **External tools**: Check for tool availability before running (mafft, trimal, phyml)
- **Coordinates**: Robust parsing of various coordinate formats from BOLD
- **Visualization**: Use existing scientific color palettes, ensure figures are publication-ready

## Reference Implementation

The original analysis scripts are:
- `msa_to_consensus.py`: Sequence clustering pipeline
- `consensus_group_to_metadata.py`: Genotype assignment

These should be refactored into the modular package structure while maintaining their core functionality.

## Specific Implementation Details

### Ocean Basin Assignment
- **Reference**: GOaS_v1_20211214 (Global Oceans and Seas v1)
- **Location**: `boldgenotyper/GOaS_v1_20211214/`
- **Files**: goas_v01.shp (and associated .cpg, .dbf, .prj, .shx)
- **Method**: Use geopandas to read shapefiles and perform point-in-polygon spatial joins
- **Attribution**: Include LICENSE_GOAS_v1.txt and cite in documentation

### Phylogenetic Analysis Options
**Mode 1**: No phylogeny (default) - genotype identification only
**Mode 2**: `--phylogeny` - Build tree from consensus sequences, midpoint rooted
**Mode 3**: `--phylogeny --outgroup <fasta>` - User provides outgroup sequences

Phylogenetic analysis is optional because:
- Outgroup selection requires taxonomic expertise
- Not all users need phylogenetic trees
- Makes pipeline more flexible for exploratory analysis

### Color Palette
- Auto-assign colors based on number of genotypes detected
- Use colorblind-friendly palette (e.g., seaborn's "colorblind" or "tab10")
- For consistency with reference analysis, first 3 colors should approximate: purple (#9D7ABE), teal (#5AB4AC), yellow (#F2CC8F)
- Scale gracefully to 20+ genotypes

### Figure Output
- **Formats**: Both PNG (300 DPI) and PDF (vector)
- **Naming**: `{organism}_distribution_map.png`, `{organism}_ocean_basins.pdf`, etc.
- **Size**: Publication-ready (e.g., 10x6 inches for maps)

### Clustering Parameters
- **Default threshold**: 0.01 (99% sequence identity)
- **User configurable**: `--cluster-threshold` flag
- **Rationale**: 99% identity is standard for COI species delimitation

## Questions to Address

When developing, consider:
1. How to handle organisms with very few samples (<10)?
2. What to do when no clear genotype clusters emerge?
3. How to handle samples with identical sequences from different basins?
4. What's the best way to visualize >10 genotypes?
5. Should we include statistical tests for geographic partitioning?

## Success Criteria

The package is ready for release when:
- ✅ It successfully processes the example S. lewini dataset
- ✅ It produces identical results to the reference analysis
- ✅ Documentation allows a user to install and run without external help
- ✅ All tests pass
- ✅ Code follows PEP 8 style guidelines
- ✅ Example workflows are provided
- ✅ GitHub repository is well-organized

## Contact & Funding

- Primary Developer: Steph Smith (steph.smith@unc.edu)
- Institution: University of North Carolina, Institute of Marine Sciences
- Manuscript: "Ocean basin-scale genetic partitioning in *Sphyrna lewini* revealed through COI sequence analysis"