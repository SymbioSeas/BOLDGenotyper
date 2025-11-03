# BOLDGenoTyper Development Guide for Claude Code

## Project Overview

**BOLDGenoTyper** is a bioinformatics pipeline for identifying and analyzing COI (Cytochrome Oxidase I) genotypes from the BOLD (Barcode of Life Database) for any taxonomic group. This package enables researchers to reproduce the workflow described in the companion manuscript analyzing *Sphyrna lewini* genotypes separated by ocean basin.

## Core Objectives

1. **Reproducibility**: Enable any researcher to analyze COI genotypes for their organism of interest
2. **Accessibility**: Design for users without extensive bioinformatics experience
3. **Standardization**: Use consistent naming conventions and file structures
4. **Documentation**: Provide comprehensive documentation for every analysis step
5. **Future Web App**: Architecture should support eventual web-based interface

## Package Architecture

### Input Requirements
- **Primary Input**: TSV file downloaded directly from BOLD database for organism of interest
  - File naming convention: `Genus_species_commonname.tsv` (e.g., `Sphyrna_lewini_scallopedhammerhead.tsv`)
  - Must contain columns: `processid`, `lat`, `lon`, `nucleotides`
  - Optional metadata columns will be preserved

### Core Analysis Workflow

1. **Data Parsing** (`parse_bold_data.py`)
   - Extract organism name from filename using pattern: `Genus_species_commonname`
   - Parse TSV to extract sequences and metadata
   - Validate required columns exist
   - Generate FASTA file: `{organism}_COI.fasta`

2. **Sequence Dereplication** (`msa_to_consensus.py`)
   - Align sequences with MAFFT (--auto mode)
   - Trim alignments with trimAl (-automated1)
   - Compute pairwise distance matrix (1 - % identity, ignoring gaps/Ns)
   - Hierarchical clustering (default threshold: 0.01)
   - Generate majority-rule consensus per cluster
   - **Output**: `{organism}_consensus.fasta` with headers: `>consensus_c{cluster_id}_n{sample_count}`
   - **Assumptions**:
     - Gaps ('-') and ambiguous bases ('N') are ignored in distance calculations
     - Sequences with <90% identity to all consensus groups remain unassigned
     - Majority rule: most common base at each position becomes consensus

3. **Genotype Assignment** (`consensus_group_to_metadata.py`)
   - Map each sample (processid) to best-matching consensus sequence
   - Use global edit distance (edlib if available, else Levenshtein)
   - Calculate identity: 1 - (edit_distance / max_length)
   - Assign to consensus group if identity ≥ 90%
   - **Output**: `{organism}_with_consensus.tsv` (adds `consensus_group` column)

4. **Geographic Filtering** (`filter_coordinates.py`)
   - Remove samples without precise lat/lon coordinates
   - Exclude samples with: 
     - Missing coordinates
     - Country-level only coordinates
     - "Coordinates from country centroid" notation
   - **Rationale**: Avoid ambiguity in ocean basin assignments for countries bordering multiple basins
   - **Output**: `{organism}_filtered.tsv`

5. **Ocean Basin Assignment** (`assign_ocean_basins.py`)
   - Classify each sample by ocean basin using coordinate boundaries:
     - North Pacific Ocean
     - South Pacific Ocean
     - North Atlantic Ocean
     - South Atlantic Ocean
     - Indian Ocean
     - Southern Ocean
     - Arctic Ocean
   - **Output**: `{organism}_with_basins.tsv` (adds `ocean_basin` column)

6. **Visualization** (`generate_figures.py`)
   - **Figure 1**: Global distribution map
     - Plot sample coordinates on world map
     - Color-code by genotype
     - Size circles by sample count at location
   - **Figure 2**: Ocean basin genotype abundance
     - Stacked bar chart showing relative abundance per basin
   - Save as: `{organism}_distribution_map.png`, `{organism}_basin_abundance.png`

7. **Phylogenetic Analysis** (`run_phylogenetics.py`)
   - Input: Consensus sequences from multiple taxa (if comparative analysis desired)
   - MAFFT alignment (--auto)
   - PhyML maximum likelihood tree:
     - GTR nucleotide substitution model
     - 1,000 bootstrap replicates
   - **Output**: `{organism}_tree.newick`, `{organism}_tree.png`

8. **Sequence Alignment Visualization** (`visualize_alignment.py`)
   - Pairwise alignment of dominant genotypes
   - Highlight polymorphic sites
   - Calculate % identity
   - **Output**: `{organism}_alignment_comparison.png`

### File Structure
```
BOLDGenoTyper/
├── boldgenotyper/
│   ├── __init__.py
│   ├── parse_data.py
│   ├── dereplicate.py
│   ├── assign_genotypes.py
│   ├── filter_coords.py
│   ├── assign_basins.py
│   ├── visualize.py
│   ├── phylogenetics.py
│   └── utils.py
├── scripts/
│   └── run_full_pipeline.py
├── tests/
│   └── test_pipeline.py
├── docs/
│   ├── installation.md
│   ├── usage.md
│   ├── tutorial.md
│   └── api_reference.md
├── example_data/
│   └── Sphyrna_lewini_scallopedhammerhead.tsv
├── environment.yml
├── setup.py
├── README.md
├── LICENSE
└── .gitignore
```

## Development Environment

### Conda Environment Setup
```bash
# Clone existing 'depredation' environment to 'boldgenotyper'
conda create --name boldgenotyper --clone depredation

# Activate new environment
conda activate boldgenotyper

# Export environment for reproducibility
conda env export --no-builds > environment.yml
```

### Required Dependencies
- **Python**: ≥3.8
- **Bioinformatics tools**:
  - MAFFT v7
  - trimAl
  - PhyML
- **Python packages**:
  - biopython
  - pandas
  - numpy
  - scipy
  - matplotlib
  - cartopy (for map plotting)
  - edlib (optional, for faster alignment)

## Critical Implementation Details

### 1. Organism Name Parsing
```python
def extract_organism_from_filename(filename):
    """
    Extract genus and species from BOLD TSV filename.
    Expected format: Genus_species_commonname.tsv
    Returns: (genus, species, full_name)
    """
    base = os.path.splitext(filename)[0]
    parts = base.split('_')
    if len(parts) < 2:
        raise ValueError(f"Filename must follow pattern: Genus_species_commonname.tsv")
    genus = parts[0]
    species = parts[1]
    return genus, species, f"{genus}_{species}"
```

### 2. ProcessID Extraction
The `processid` in FASTA headers appears after the last underscore before a dot or end-of-line:
- Example: `>Sphyrna_lewini_ANGBF11456-15.COI-5P` → processid = `ANGBF11456-15`
- Use regex: `_(?P<pid>[^.\s_]+)(?:[.\s]|$)`

### 3. Coordinate Filtering Logic
```python
def should_exclude_coordinate(lat, lon, metadata_notes):
    """
    Exclude samples with:
    - Missing lat/lon (NaN, None, empty string)
    - Metadata indicating "country centroid" or "country-level" coordinates
    """
    if pd.isna(lat) or pd.isna(lon):
        return True
    if isinstance(metadata_notes, str):
        if "centroid" in metadata_notes.lower():
            return True
        if "country" in metadata_notes.lower() and not "precise" in metadata_notes.lower():
            return True
    return False
```

### 4. Ocean Basin Boundaries
Use standard oceanographic definitions. Handle edge cases:
- Mediterranean Sea → Atlantic Ocean
- Caribbean Sea → North Atlantic Ocean
- South China Sea → (separate category or Pacific Ocean - document choice)

## Testing Strategy

### Unit Tests Required
- [ ] Filename parsing with various formats
- [ ] Sequence dereplication with known sequences
- [ ] Coordinate filtering edge cases
- [ ] Ocean basin assignment accuracy

### Integration Tests
- [ ] Full pipeline run on example dataset (S. lewini)
- [ ] Verify output files created
- [ ] Check figure generation
- [ ] Validate phylogenetic tree construction

### Example Data
Include `Sphyrna_lewini_scallopedhammerhead.tsv` as test dataset with known expected outputs.

## Documentation Requirements

### README.md Must Include:
1. **Installation instructions**
   - Conda environment setup
   - Dependency installation
   - Tool availability checks (MAFFT, trimAl, PhyML)

2. **Quick Start**
   - Download data from BOLD
   - Run single command: `boldgenotyper run input.tsv`

3. **Usage Examples**
   - Basic usage
   - Custom parameters
   - Output interpretation

4. **Citation**
   - Reference to companion manuscript
   - BOLD database citation
   - Tool citations (MAFFT, PhyML, etc.)

### docs/usage.md Must Explain:
- Each analysis step in detail
- Parameter choices and defaults
- How to interpret results
- Troubleshooting common errors

### docs/tutorial.md Should:
- Walk through S. lewini example
- Explain biological interpretation
- Show how to modify for other organisms

## GitHub Best Practices

### Repository Setup
1. **Initialize with**:
   - README.md
   - LICENSE (suggest MIT or BSD)
   - .gitignore (Python, Conda, IDE-specific)

2. **Branching Strategy**:
   - `main`: stable releases
   - `dev`: development branch
   - Feature branches: `feature/description`
   - Bug fixes: `fix/description`

3. **Commit Messages**:
   - Use conventional commits format
   - Example: `feat: add ocean basin assignment module`
   - Example: `docs: update installation instructions`

4. **Version Control**:
   - Tag releases: `v0.1.0`, `v0.2.0`, etc.
   - Follow semantic versioning

### GitHub Actions (Future)
- Automated testing on push
- Documentation building
- PyPI package deployment

## Future Web Application Considerations

### Architecture Planning
- Keep core logic in separate modules from CLI
- Design functions to accept/return standard data structures (pandas DataFrames)
- Minimize side effects (file I/O separate from computation)
- Enable progress reporting for long-running tasks

### Potential Web Framework: Streamlit
```python
# Future structure for web app
def run_pipeline_with_progress(uploaded_file, progress_callback):
    """
    progress_callback: function to report progress (0-100)
    """
    # Parse data
    progress_callback(10, "Parsing BOLD data...")
    
    # Dereplicate sequences
    progress_callback(30, "Identifying unique genotypes...")
    
    # etc.
```

## Important Notes for Claude Code

1. **Always use relative paths** - never hardcode absolute paths
2. **Extract organism name from input filename** - all outputs should use this name
3. **Preserve metadata** - don't drop columns unnecessarily from TSV
4. **Document assumptions** - especially for filtering and clustering thresholds
5. **Error handling** - provide clear error messages when:
   - Required columns missing from TSV
   - External tools (MAFFT, PhyML) not found
   - Invalid coordinates
6. **Logging** - use Python logging module to track pipeline progress
7. **Reproducibility** - Set random seeds where applicable (e.g., bootstrap replicates)

## Contact and Support

**Primary Developer**: Steph Smith (@SymbioSeas)
**Email**: steph.smith@unc.edu
**Manuscript**: [Journal submission pending - include citation once published]

## References to Include in Documentation

- BOLD Database: Ratnasingham & Hebert (2007) doi:10.1111/j.1471-8286.2007.01678.x
- MAFFT: Katoh & Standley (2013) doi:10.1093/molbev/mst010
- trimAl: Capella-Gutiérrez et al. (2009) doi:10.1093/bioinformatics/btp348
- PhyML: Guindon et al. (2010) doi:10.1093/sysbio/syq010
- COI Barcoding: Hebert et al. (2003) doi:10.1098/rspb.2002.2218

---

## Development Checklist for Claude Code

### Phase 1: Core Pipeline
- [ ] Setup package structure
- [ ] Implement BOLD data parser with organism name extraction
- [ ] Port and refactor `msa_to_consensus.py` for general use
- [ ] Port and refactor `consensus_group_to_metadata.py`
- [ ] Implement coordinate filtering with documented logic
- [ ] Implement ocean basin assignment

### Phase 2: Visualization & Analysis
- [ ] Global distribution map generator
- [ ] Ocean basin abundance plotter
- [ ] Phylogenetic analysis wrapper
- [ ] Sequence alignment visualization

### Phase 3: Documentation & Testing
- [ ] Write comprehensive README
- [ ] Create usage documentation
- [ ] Write tutorial using S. lewini example
- [ ] Implement unit tests
- [ ] Add integration test with example data

### Phase 4: Packaging & Distribution
- [ ] Create setup.py for pip installation
- [ ] Finalize environment.yml
- [ ] Add CLI entry point: `boldgenotyper`
- [ ] Test installation on clean environment
- [ ] Prepare for PyPI distribution

### Phase 5: Future Enhancement (Post-Initial Release)
- [ ] Streamlit web application
- [ ] Docker containerization
- [ ] Galaxy tool wrapper
- [ ] Batch processing mode

---

**This document should be the primary reference for all development work on BOLDGenoTyper.**