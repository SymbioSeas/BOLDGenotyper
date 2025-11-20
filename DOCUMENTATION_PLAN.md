# BOLDGenotyper Documentation Update Plan for Publication

**Date**: 2025-11-17
**Status**: Draft Plan for Review
**Target Audience**: New users, researchers, developers
**Goal**: Create publication-ready, comprehensive documentation

---

## Executive Summary

BOLDGenotyper currently has **75-80% documentation coverage** with strong foundational documentation. This plan outlines specific improvements needed to achieve **90%+ coverage** and create publication-ready documentation suitable for peer-reviewed publication, new users, and developers.

### Current Strengths
- ✅ Excellent module-level docstrings (90%)
- ✅ Comprehensive function documentation (93%)
- ✅ Nearly complete type hints (94%)
- ✅ Well-written README and QUICKSTART guides (95%+)
- ✅ Good scientific rationale documentation

### Critical Gaps
- ❌ No auto-generated API documentation
- ❌ Limited tutorial/walkthrough content (40%)
- ❌ Incomplete troubleshooting reference (60%)
- ❌ Missing advanced usage examples

---

## Documentation Goals

### 1. Publication Requirements
- [ ] Comprehensive user guide suitable for methods section citation
- [ ] Clear installation and setup instructions
- [ ] Reproducible examples with test data
- [ ] Parameter documentation with biological justification
- [ ] Citation information and bibliography

### 2. New User Requirements
- [ ] Quick start guide (5-10 minutes to first result)
- [ ] Step-by-step tutorials with real data
- [ ] Common use cases documented
- [ ] Clear error messages and troubleshooting
- [ ] Video tutorials or animated GIFs (optional)

### 3. Developer Requirements
- [ ] API documentation (auto-generated)
- [ ] Code architecture overview
- [ ] Contribution guidelines
- [ ] Development setup instructions
- [ ] Testing documentation

---

## Detailed Plan by Documentation Type

## Phase 1: Core User Documentation (Priority: CRITICAL)
**Estimated Time**: 16-20 hours
**Target Completion**: Week 1

### 1.1 Update README.md
**Current Status**: 98% complete, 994 lines
**Action Items**:

- [x] Comprehensive overview ✓ (already excellent)
- [x] Installation instructions ✓ (already excellent)
- [x] GOaS setup guide ✓ (already excellent)
- [ ] **ADD**: Link to full API documentation (once created)
- [ ] **ADD**: "For Developers" section with link to CONTRIBUTING.md
- [ ] **ADD**: "How to Cite" section with BibTeX
- [ ] **ADD**: Badges (build status, coverage, version, license, DOI)
- [ ] **ENHANCE**: Add animated GIF or screenshot of output
- [ ] **REVIEW**: Ensure all example commands are tested and work
- [ ] **UPDATE**: Version history section with changelog

**New Sections to Add**:
```markdown
## For Developers
See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.

## How to Cite
If you use BOLDGenotyper in your research, please cite:

**Software Citation**:
```bibtex
@software{boldgenotyper2025,
  author = {Smith, Steph},
  title = {BOLDGenotyper: Automated COI Sequence Genotyping},
  year = {2025},
  url = {https://github.com/SymbioSeas/BOLDGenotyper},
  version = {0.1.0}
}
```

**Primary Publication**:
[To be added upon publication]

## Badges
[![Build Status](badge_url)](link)
[![Documentation](badge_url)](link)
[![License: MIT](badge_url)](link)
```

### 1.2 Enhance QUICKSTART.md
**Current Status**: 95% complete
**Action Items**:

- [ ] **ADD**: Prerequisites checklist (Python version, conda, etc.)
- [ ] **ADD**: Expected runtime for different dataset sizes
- [ ] **ADD**: Common pitfalls section
- [ ] **ADD**: Link to full tutorial for complex analyses
- [ ] **ENHANCE**: Add output interpretation guide
- [ ] **TEST**: Verify all commands work on fresh install

### 1.3 Create INSTALLATION.md (NEW)
**Target**: Comprehensive installation guide
**Contents**:

```markdown
# Installation Guide

## Table of Contents
- Quick Install (5 minutes)
- Detailed Install
- Platform-Specific Instructions
  - Linux
  - macOS
  - Windows (WSL)
- Docker Installation
- HPC/Cluster Installation
- Troubleshooting Installation
- Verifying Installation

## Quick Install
[Basic conda/pip instructions]

## Platform-Specific Instructions

### Linux
[Detailed Ubuntu, Debian, CentOS instructions]

### macOS
[Homebrew + conda instructions]

### Windows
[WSL2 setup + conda]

## Docker Installation
```bash
docker pull symbioseas/boldgenotyper:latest
docker run -v $(pwd):/data boldgenotyper data.tsv
```

## HPC/Cluster Installation
[SLURM examples, module loading]

## Troubleshooting
- Common errors during installation
- Dependency conflicts
- External tool issues
```

### 1.4 Create TUTORIAL.md (NEW)
**Target**: Step-by-step walkthrough with real data
**Contents**:

```markdown
# BOLDGenotyper Tutorial: Complete Walkthrough

## Learning Objectives
By the end of this tutorial, you will:
- Download data from BOLD database
- Run the complete pipeline
- Interpret results
- Customize parameters for your study

## Tutorial Dataset
We'll use *Sphyrna lewini* (scalloped hammerhead shark) as example.

## Step 1: Download Data from BOLD
[Detailed screenshots and instructions]

## Step 2: Prepare Your Data
[Data quality checks, file naming]

## Step 3: Run Basic Analysis
```bash
boldgenotyper data/Sphyrna_lewini.tsv
```
[Expected output, runtime, etc.]

## Step 4: Interpret Results
### Main Output File
[Column-by-column explanation of CSV]

### Genotype Assignments
[How to interpret identity scores]

### Geographic Distributions
[Reading the maps]

### Phylogenetic Trees
[If applicable]

## Step 5: Advanced Customization
[Adjusting parameters for different scenarios]

## Step 6: Publication-Ready Figures
[Exporting and formatting for papers]

## Troubleshooting Common Issues
[With this specific dataset]

## Next Steps
[Links to advanced topics]
```

### 1.5 Update TESTING_CLI.md
**Current Status**: 90% complete
**Action Items**:

- [ ] **ADD**: Integration with Python API examples
- [ ] **ADD**: Batch processing examples
- [ ] **ADD**: Custom configuration examples
- [ ] **ENHANCE**: Performance benchmarks for different dataset sizes

---

## Phase 2: API Documentation (Priority: CRITICAL)
**Estimated Time**: 12-16 hours
**Target Completion**: Week 1-2

### 2.1 Set Up Sphinx Documentation
**Action Items**:

1. **Install Sphinx and extensions**:
```bash
pip install sphinx sphinx-rtd-theme sphinx-autodoc-typehints
```

2. **Create docs/ directory structure**:
```
docs/
├── source/
│   ├── conf.py
│   ├── index.rst
│   ├── api/
│   │   ├── cli.rst
│   │   ├── config.rst
│   │   ├── dereplication.rst
│   │   ├── genotype_assignment.rst
│   │   ├── geographic.rst
│   │   ├── metadata.rst
│   │   ├── phylogenetics.rst
│   │   ├── reports.rst
│   │   ├── utils.rst
│   │   └── visualization.rst
│   ├── tutorials/
│   │   ├── quickstart.rst
│   │   ├── basic_usage.rst
│   │   └── advanced_usage.rst
│   ├── guides/
│   │   ├── installation.rst
│   │   ├── troubleshooting.rst
│   │   └── parameters.rst
│   └── development/
│       ├── contributing.rst
│       ├── architecture.rst
│       └── testing.rst
├── Makefile
└── make.bat
```

3. **Configure autodoc**:
```python
# conf.py
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx_autodoc_typehints',
]
```

4. **Generate initial documentation**:
```bash
sphinx-quickstart docs
sphinx-apidoc -o docs/source/api boldgenotyper
cd docs && make html
```

5. **Host on Read the Docs**:
- Create .readthedocs.yaml
- Connect GitHub repository
- Enable automatic builds

### 2.2 API Reference Pages
**Create comprehensive API docs for each module**:

**Example structure** (api/dereplication.rst):
```rst
Sequence Dereplication
======================

.. automodule:: boldgenotyper.dereplication
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

.. autofunction:: calculate_pairwise_distances
.. autofunction:: cluster_sequences
.. autofunction:: generate_consensus

Classes
-------

.. autoclass:: DereplicationConfig
   :members:

Examples
--------

Basic Usage
^^^^^^^^^^^

.. code-block:: python

    from boldgenotyper import dereplication

    # Load sequences
    sequences = [...]

    # Calculate distances
    distances = dereplication.calculate_pairwise_distances(sequences)

    # Cluster
    clusters = dereplication.cluster_sequences(distances, threshold=0.01)

Advanced Usage
^^^^^^^^^^^^^^

[More examples]
```

---

## Phase 3: Troubleshooting & FAQ (Priority: HIGH)
**Estimated Time**: 8-10 hours
**Target Completion**: Week 2

### 3.1 Create TROUBLESHOOTING.md (NEW)
**Contents**:

```markdown
# Troubleshooting Guide

## Table of Contents
- Installation Issues
- Data Input Issues
- Pipeline Failures
- Geographic Analysis Issues
- Visualization Issues
- Performance Issues
- Error Code Reference

## Installation Issues

### Error: "MAFFT not found"
**Symptoms**: Pipeline fails at alignment step
**Cause**: MAFFT not installed or not in PATH
**Solution**:
```bash
conda install -c bioconda mafft
# or
brew install mafft
```

### Error: "GOaS shapefile not found"
**Symptoms**: Pipeline continues but skips geographic analysis
**Cause**: GOaS shapefile not downloaded
**Solution**:
```bash
python -m boldgenotyper.goas_downloader
```
**Alternative**: Use `--no-geo` flag

[Continue for 20-30 common errors]

## Data Input Issues

### Low Assignment Rate (<30%)
**Symptoms**: Most samples unassigned to genotypes
**Diagnosis**:
1. Check diagnostics CSV for identity scores
2. Review sequence quality

**Solutions**:
- High diversity → Lower `--similarity-threshold` (try 0.80 or 0.75)
- Poor quality → Pre-filter sequences
- Multiple species → Run separately by species

[More scenarios]

## Error Code Reference

### E001: Missing Required Column
**Message**: `Column 'processid' not found in TSV`
**Cause**: Input TSV missing required column
**Fix**: Ensure TSV has columns: processid, nuc, species

[All error codes]

## Diagnostic Workflows

### My pipeline failed, now what?
1. Check log file: `{Organism}_output/{Organism}_pipeline.log`
2. Look for ERROR or CRITICAL messages
3. Find error code in reference above
4. If not found, check GitHub issues

### How do I debug low assignment rates?
[Flowchart or decision tree]
```

### 3.2 Create FAQ.md (NEW)
**Common questions**:

```markdown
# Frequently Asked Questions

## General

### What organisms can I analyze?
Any organism with COI sequences in BOLD database. Currently optimized for marine organisms; terrestrial/freshwater support coming soon.

### How much data do I need?
Minimum: 50 samples
Recommended: 200+ samples for robust genotype detection
Large-scale: 1000+ samples

### How long does analysis take?
- Small dataset (<500): 2-5 minutes
- Medium (500-2000): 5-15 minutes
- Large (>2000): 15-60 minutes

[30-40 questions]

## Technical

### What's the difference between --similarity-threshold and --min-identity?
[Detailed explanation]

### Should I use --build-tree?
Depends on your goals:
- Publication: Yes (helpful for visualization)
- Quick analysis: No (saves time)
- Large dataset: Maybe (can be slow)

[More technical questions]

## Biological

### What does 99% identity mean?
[Biological interpretation]

### How are ocean basins defined?
[GOaS explanation]

[More biological questions]
```

---

## Phase 4: Input/Output Documentation (Priority: HIGH)
**Estimated Time**: 6-8 hours
**Target Completion**: Week 2

### 4.1 Create DATA_FORMATS.md (NEW)
**Contents**:

```markdown
# Data Formats Reference

## Input Formats

### BOLD TSV Format
**File Extension**: `.tsv`
**Delimiter**: Tab
**Encoding**: UTF-8

**Required Columns**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `processid` | string | BOLD process ID (unique) | `ANGBF11456-15` |
| `nuc` or `nucleotides` | string | DNA sequence (A,C,G,T,N) | `ACGTACGT...` |
| `species` | string | Species name | `Sphyrna lewini` |

**Optional but Recommended**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `lat` | float | Latitude (-90 to 90) | `25.7617` |
| `lon` | float | Longitude (-180 to 180) | `-80.1918` |
| `coord` | string | Coordinates (lat, lon) | `25.7617, -80.1918` |
| `country` | string | Country of collection | `Mexico` |
| `coord_accuracy` | string | Coordinate precision | `precise` or `country centroid` |

**Example**:
```tsv
processid	species	nuc	lat	lon	country
ANGBF11456-15	Sphyrna lewini	ACGTACGT...	25.76	-80.19	USA
```

**Common Issues**:
- Mixed column names (`nuc` vs `nucleotides`)
- Missing coordinates
- Country centroids vs precise GPS

**Validation**:
```python
from boldgenotyper import metadata
df = metadata.parse_bold_tsv("data.tsv")
# Automatically handles column name variations
```

### FASTA Format
[For custom workflows]

### Configuration Files
[YAML/JSON format documentation]

## Output Formats

### Main Results: {Organism}_annotated.csv
**Purpose**: Complete integrated dataset

**Columns** (40+ total):

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `processid` | string | Sample identifier | `ANGBF11456-15` |
| `consensus_group` | string | Assigned genotype | `consensus_c7_n381` |
| `consensus_group_sp` | string | Genotype with species | `Sphyrna_consensus_c7` |
| `best_identity` | float | Match confidence (0-1) | `0.957` |
| `assigned_sp` | string | Species in genotype | `Sphyrna lewini` |
| `ocean_basin` | string | Ocean basin | `North Atlantic` |

[All columns documented]

### Consensus Sequences: {Organism}_consensus.fasta
**Purpose**: Representative sequence per genotype

**Header Format**:
```
>consensus_c{cluster_id}_n{sample_count}
```

**Example**:
```fasta
>consensus_c7_n381
ACGTACGTACGT...
>consensus_c12_n145
TGCATGCATGCA...
```

### Diagnostics: {Organism}_diagnostics.csv
[Detailed documentation]

[All output formats documented]
```

---

## Phase 5: Developer Documentation (Priority: MEDIUM)
**Estimated Time**: 10-12 hours
**Target Completion**: Week 3

### 5.1 Create CONTRIBUTING.md (NEW)
**Contents**:

```markdown
# Contributing to BOLDGenotyper

## Welcome!
Thank you for considering contributing to BOLDGenotyper!

## Code of Conduct
[Placeholder or link]

## How Can I Contribute?

### Reporting Bugs
[Issue template, what to include]

### Suggesting Enhancements
[Feature request template]

### Contributing Code
[Pull request process]

## Development Setup

### 1. Fork and Clone
```bash
git clone https://github.com/YOUR_USERNAME/BOLDGenotyper.git
cd BOLDGenotyper
```

### 2. Create Development Environment
```bash
conda env create -f boldgenotyper_env.yml
conda activate depredation
pip install -e ".[dev]"
```

### 3. Install Pre-commit Hooks
```bash
pip install pre-commit
pre-commit install
```

### 4. Run Tests
```bash
pytest tests/
```

## Code Style

### Python Style Guide
- Follow PEP 8
- Use Black formatter
- Maximum line length: 100
- Type hints required

### Documentation Style
- Google-style docstrings
- Examples in docstrings
- Type hints in signatures

**Example**:
```python
def calculate_distance(seq1: str, seq2: str, ignore_gaps: bool = True) -> float:
    """Calculate sequence distance.

    Args:
        seq1: First sequence (A,C,G,T,N,-)
        seq2: Second sequence
        ignore_gaps: If True, exclude gaps from distance calculation

    Returns:
        Distance between sequences (0.0 to 1.0)

    Examples:
        >>> calculate_distance("ACGT", "ACGT")
        0.0
        >>> calculate_distance("ACGT", "TGCA")
        1.0
    """
```

## Testing Guidelines
[How to write tests, coverage requirements]

## Pull Request Process
1. Create feature branch
2. Write tests
3. Update documentation
4. Run tests and linters
5. Submit PR with description

## Architecture Overview
[Link to ARCHITECTURE.md]
```

### 5.2 Create ARCHITECTURE.md (NEW)
**Contents**:

```markdown
# BOLDGenotyper Architecture

## System Overview

```
┌─────────────────────────────────────────────────────────┐
│                     BOLDGenotyper                       │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  ┌─────────┐  ┌──────────────┐  ┌─────────────────┐  │
│  │   CLI   │──│  Pipeline    │──│  Configuration  │  │
│  └─────────┘  └──────────────┘  └─────────────────┘  │
│                      │                                  │
│       ┌──────────────┼──────────────┐                 │
│       │              │              │                  │
│  ┌─────────┐  ┌────────────┐  ┌──────────┐          │
│  │Metadata │  │Dereplica-  │  │Geographic│          │
│  │ Parser  │  │   tion     │  │ Analysis │          │
│  └─────────┘  └────────────┘  └──────────┘          │
│       │              │              │                  │
│  ┌─────────┐  ┌────────────┐  ┌──────────┐          │
│  │Genotype │  │Phylo-      │  │Visualiza-│          │
│  │Assign.  │  │genetics    │  │  tion    │          │
│  └─────────┘  └────────────┘  └──────────┘          │
│                                                         │
│  ┌───────────────────────────────────────────────┐   │
│  │            External Tools (MAFFT, etc.)       │   │
│  └───────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
```

## Module Descriptions

### Core Modules

#### cli.py
**Purpose**: Command-line interface and pipeline orchestration
**Dependencies**: All other modules
**Key Functions**: `run_pipeline()`, `setup_directories()`

#### config.py
**Purpose**: Configuration management via dataclasses
**Dependencies**: None (standalone)
**Design Pattern**: Hierarchical configuration with defaults

[All modules documented]

## Data Flow

```
BOLD TSV
   │
   ├──> metadata.parse_bold_tsv()
   │       │
   │       └──> DataFrame
   │              │
   ├──> dereplication.run_dereplication()
   │       │
   │       └──> Consensus FASTA
   │              │
   ├──> genotype_assignment.assign_genotypes()
   │       │
   │       └──> Annotated TSV
   │              │
   ├──> geographic.assign_ocean_basins()
   │       │
   │       └──> Geographic TSV
   │              │
   ├──> visualization.create_maps()
   │       │
   │       └──> Figures (PNG, PDF, SVG)
   │              │
   └──> reports.generate_summary()
           │
           └──> Summary TXT
```

## Design Principles

1. **Modularity**: Each module has single responsibility
2. **Type Safety**: Extensive use of type hints
3. **Configuration**: Centralized, hierarchical config
4. **Error Handling**: Graceful degradation
5. **Logging**: Comprehensive logging at all levels

## Extension Points
[How to add new features]

## Performance Considerations
[Bottlenecks, optimization strategies]
```

### 5.3 Create TESTING.md (NEW)
**Testing documentation**:

```markdown
# Testing Guide

## Running Tests

### All Tests
```bash
pytest tests/
```

### Specific Module
```bash
pytest tests/test_dereplication.py
```

### With Coverage
```bash
pytest --cov=boldgenotyper tests/
```

## Test Organization

```
tests/
├── test_metadata.py          # Metadata parsing tests
├── test_dereplication.py     # Sequence clustering tests
├── test_genotype_assignment.py
├── test_geographic.py
├── test_phylogenetics.py
├── test_visualization.py
└── data/                     # Test datasets
    ├── Euprymna_scolopes.tsv
    └── test_sequences.fasta
```

## Writing Tests

[Guidelines, examples]

## Test Data

[Where to get, how to create]
```

---

## Phase 6: Examples and Advanced Usage (Priority: MEDIUM)
**Estimated Time**: 8-10 hours
**Target Completion**: Week 3

### 6.1 Create EXAMPLES.md (NEW)
**Use case scenarios**:

```markdown
# BOLDGenotyper Examples

## Table of Contents
- Basic Analysis
- High Diversity Taxa
- Geographic Analysis
- Phylogenetic Analysis
- Batch Processing
- Custom Configuration
- Programmatic Usage (Python API)

## Basic Analysis

### Single Species Analysis
```bash
boldgenotyper data/Sphyrna_lewini.tsv
```

## High Diversity Taxa

### Lowering Threshold for Diverse Genera
For genera like *Carcharhinus* with high genetic diversity:

```bash
boldgenotyper data/Carcharhinus.tsv \
    --similarity-threshold 0.80 \
    --min-identity 0.85
```

**When to use**:
- Assignment rate <30% with default settings
- Many species in genus
- Cryptic species complex

## Geographic Analysis

### Marine Organisms with GOaS
```bash
boldgenotyper data/Sphyrna_lewini.tsv \
    --output results/shark_geo_analysis
```

### Terrestrial Organisms (Skip Geographic)
```bash
boldgenotyper data/Terrestrial_species.tsv --no-geo
```

## Phylogenetic Analysis

### With Tree Building
```bash
boldgenotyper data/Sphyrna_lewini.tsv \
    --build-tree \
    --tree-method fasttree
```

## Batch Processing

### Multiple Species
```bash
for file in data/*.tsv; do
    boldgenotyper "$file" --output "results/$(basename $file .tsv)"
done
```

### Parallel Processing
```bash
ls data/*.tsv | parallel -j 4 boldgenotyper {} --output results/{/.}
```

## Custom Configuration

### Using Config File
Create `config.yaml`:
```yaml
dereplication:
  clustering_threshold: 0.005  # 99.5% identity

genotype_assignment:
  min_identity: 0.95

geographic:
  exclude_centroids: true

visualization:
  figure_dpi: 600
  figure_format: [png, pdf, svg]
```

Run with:
```bash
boldgenotyper data.tsv --config config.yaml
```

## Programmatic Usage (Python API)

### Basic Pipeline
```python
from pathlib import Path
from boldgenotyper import metadata, dereplication, genotype_assignment

# Load data
df = metadata.parse_bold_tsv("data.tsv")

# Extract sequences
sequences = metadata.extract_sequences(df)

# Run dereplication
consensus = dereplication.run_dereplication(
    sequences,
    threshold=0.01,
    output_dir=Path("results/consensus")
)

# Assign genotypes
df_annotated = genotype_assignment.assign_genotypes(
    df,
    sequences,
    consensus,
    min_identity=0.90
)

# Save results
df_annotated.to_csv("results/annotated.csv", index=False)
```

### Custom Workflows
[More advanced examples]
```

### 6.2 Create Jupyter Notebook Tutorials
**Target**: 2-3 interactive notebooks

1. **Tutorial 1: Basic Workflow**
   - Download data
   - Run pipeline
   - Explore results

2. **Tutorial 2: Parameter Tuning**
   - Impact of thresholds
   - Diagnostic interpretation
   - Optimization strategies

3. **Tutorial 3: Advanced Analysis**
   - Custom visualizations
   - Statistical analysis
   - Publication figures

---

## Phase 7: Publication Package (Priority: HIGH)
**Estimated Time**: 6-8 hours
**Target Completion**: Week 2-3

### 7.1 Create CITATION.cff (NEW)
**GitHub citation file**:

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
  - family-names: Smith
    given-names: Steph
    email: steph.smith@unc.edu
    affiliation: "University of North Carolina, Institute of Marine Sciences"
    orcid: "https://orcid.org/XXXX-XXXX-XXXX-XXXX"
title: "BOLDGenotyper: Automated COI Sequence Genotyping"
version: 0.1.0
date-released: 2025-XX-XX
url: "https://github.com/SymbioSeas/BOLDGenotyper"
repository-code: "https://github.com/SymbioSeas/BOLDGenotyper"
keywords:
  - bioinformatics
  - COI
  - barcoding
  - genotyping
  - phylogeography
license: MIT
```

### 7.2 Update LICENSE
**Ensure proper license**:
- Choose: MIT (recommended for scientific software)
- Add license file if missing
- Update all source file headers

### 7.3 Create CHANGELOG.md (NEW)
**Version history**:

```markdown
# Changelog

All notable changes to BOLDGenotyper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2025-XX-XX

### Added
- Initial release
- Complete pipeline from BOLD TSV to annotated results
- Sequence dereplication and consensus generation
- Geographic analysis with GOaS integration
- Phylogenetic tree building (optional)
- Publication-ready visualizations
- CLI interface
- Comprehensive documentation

### Known Issues
- Phylogenetics module has some stub functions
- Limited to marine organisms for geographic analysis
- Large datasets (>10,000 seqs) may have memory issues

## [0.0.1] - 2024-XX-XX (Development)
- Initial development version
```

### 7.4 Create DOI via Zenodo
**Steps**:
1. Link GitHub to Zenodo
2. Create release on GitHub
3. Zenodo automatically creates DOI
4. Add DOI badge to README

### 7.5 Create Software Paper (Optional)
**Target Journal**: JOSS (Journal of Open Source Software)

**Draft paper.md**:
```markdown
---
title: 'BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis'
tags:
  - Python
  - bioinformatics
  - DNA barcoding
  - COI
  - phylogeography
authors:
  - name: Steph Smith
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: University of North Carolina, Institute of Marine Sciences, USA
   index: 1
date: XX Month 2025
bibliography: paper.bib
---

# Summary

BOLDGenotyper is a Python package for automated genotyping of COI (Cytochrome Oxidase I) sequences from the BOLD (Barcode of Life Database) database. The software addresses the need for reproducible, standardized workflows in DNA barcode analysis for phylogeographic and population genetic studies.

# Statement of Need

[Why this software is needed]

# Key Features

[Core functionality]

# Usage Example

```python
from boldgenotyper import cli
cli.run_pipeline("data.tsv", ...)
```

# Acknowledgments

[Funding, collaborators]

# References
```

---

## Phase 8: Final Polish and Review (Priority: MEDIUM)
**Estimated Time**: 6-8 hours
**Target Completion**: Week 4

### 8.1 Documentation Review Checklist

#### Content Review
- [ ] All code examples tested and working
- [ ] All links functional
- [ ] No broken references
- [ ] Consistent terminology throughout
- [ ] Correct grammar and spelling
- [ ] Appropriate technical level for audience

#### Completeness Review
- [ ] All modules documented
- [ ] All functions have docstrings
- [ ] All parameters explained
- [ ] All outputs described
- [ ] All error messages documented

#### Usability Review
- [ ] Can a new user get started in <10 minutes?
- [ ] Are common questions answered?
- [ ] Are troubleshooting steps clear?
- [ ] Are examples realistic and useful?

#### Scientific Review
- [ ] All parameter defaults justified
- [ ] Citations included where appropriate
- [ ] Methods reproducible
- [ ] Assumptions clearly stated

### 8.2 Cross-Reference Check
**Ensure consistency across documents**:
- Version numbers match
- Installation instructions consistent
- Example commands work everywhere
- Parameter names consistent

### 8.3 Generate Documentation Website
**Using Sphinx or MkDocs**:

```bash
# Option 1: Sphinx
cd docs
make html
# Deploy to GitHub Pages or Read the Docs

# Option 2: MkDocs
mkdocs build
mkdocs gh-deploy
```

---

## Success Metrics

### Quantitative Goals
- [ ] 95%+ function docstring coverage
- [ ] 100% module docstring coverage
- [ ] 80%+ code comment coverage for complex logic
- [ ] All public functions have examples
- [ ] Zero broken links in documentation
- [ ] <10 minutes to first successful run (new user)

### Qualitative Goals
- [ ] New user can install and run pipeline without assistance
- [ ] Troubleshooting guide covers 90% of common issues
- [ ] Documentation suitable for Methods section in publication
- [ ] Positive feedback from beta testers
- [ ] Accepted to JOSS or similar journal (optional)

---

## Timeline and Milestones

### Week 1: Critical Documentation
- **Days 1-2**: README updates, badges, citation
- **Days 3-4**: Sphinx setup, API docs generation
- **Days 5-6**: INSTALLATION.md, TUTORIAL.md
- **Day 7**: Review and testing

### Week 2: Comprehensive Coverage
- **Days 1-2**: TROUBLESHOOTING.md, FAQ.md
- **Days 3-4**: DATA_FORMATS.md, complete API docs
- **Days 5-6**: Publication package (DOI, CITATION.cff)
- **Day 7**: Review and testing

### Week 3: Developer & Advanced
- **Days 1-2**: CONTRIBUTING.md, ARCHITECTURE.md
- **Days 3-4**: EXAMPLES.md, Jupyter notebooks
- **Days 5-6**: TESTING.md, advanced examples
- **Day 7**: Review and testing

### Week 4: Polish and Deploy
- **Days 1-2**: Comprehensive review and fixes
- **Days 3-4**: Documentation website deployment
- **Days 5-6**: Beta testing with external users
- **Day 7**: Final adjustments and release

---

## Priority Matrix

### Must Have (Critical for Publication)
1. ✅ README.md updates (badges, citation, DOI)
2. ✅ API documentation (Sphinx/MkDocs)
3. ✅ INSTALLATION.md
4. ✅ TUTORIAL.md
5. ✅ TROUBLESHOOTING.md
6. ✅ CITATION.cff
7. ✅ LICENSE

### Should Have (Important for Users)
8. ✅ FAQ.md
9. ✅ DATA_FORMATS.md
10. ✅ EXAMPLES.md
11. ✅ Jupyter notebooks
12. ✅ CHANGELOG.md

### Nice to Have (Beneficial for Community)
13. ⭐ CONTRIBUTING.md
14. ⭐ ARCHITECTURE.md
15. ⭐ TESTING.md
16. ⭐ Software paper (JOSS)
17. ⭐ Video tutorials

---

## Resource Requirements

### Tools
- [ ] Sphinx (`pip install sphinx sphinx-rtd-theme`)
- [ ] MkDocs (alternative to Sphinx)
- [ ] Jupyter (`pip install jupyter`)
- [ ] GitHub account with Pages enabled
- [ ] Zenodo account for DOI
- [ ] Read the Docs account (optional)

### Time Commitment
- **Total Estimated**: 60-80 hours
- **Per Week**: 15-20 hours
- **Duration**: 4 weeks

### Personnel
- **Primary**: Documentation lead (you)
- **Review**: 2-3 subject matter experts
- **Testing**: 3-5 new users for beta testing

---

## Maintenance Plan

### Ongoing Documentation
- [ ] Update documentation with each release
- [ ] Review and update FAQ quarterly
- [ ] Add new examples based on user questions
- [ ] Keep troubleshooting guide current

### Community Engagement
- [ ] Respond to documentation issues on GitHub
- [ ] Incorporate user feedback
- [ ] Create video tutorials based on common questions
- [ ] Host workshop or webinar (annual)

---

## Appendices

### A. Documentation Style Guide

**Tone**: Professional but friendly, approachable
**Voice**: Second person ("you can") for tutorials, third person for reference
**Technical Level**: Assume basic command-line knowledge, explain bioinformatics concepts
**Code Examples**: Always tested, include expected output
**Citations**: Use standard scientific format, include DOIs

### B. Template Library

**Function Docstring Template**:
```python
def function_name(param1: type, param2: type = default) -> return_type:
    """Short one-line summary.

    Longer description explaining what the function does, when to use it,
    and any important caveats or assumptions.

    Args:
        param1: Description of parameter 1
        param2: Description of parameter 2 (default: {default})

    Returns:
        Description of return value

    Raises:
        ExceptionType: When this exception occurs

    Examples:
        >>> function_name(arg1, arg2)
        expected_output

        >>> function_name(arg1, param2=value)
        expected_output

    References:
        Author et al. (Year). Title. Journal. DOI.
    """
```

**Tutorial Section Template**:
```markdown
## Section Title

### Learning Objectives
- Objective 1
- Objective 2

### Prerequisites
- Prerequisite 1
- Prerequisite 2

### Step-by-Step Instructions

#### Step 1: [Action]
[Explanation]

```bash
command here
```

**Expected Output**:
```
output here
```

**What This Does**:
[Explanation]

#### Step 2: [Action]
[Continue...]

### Troubleshooting
If you encounter [problem], try [solution].

### Next Steps
Now that you've completed this tutorial, you can:
- [Next topic 1]
- [Next topic 2]
```

### C. Review Checklist Template

```markdown
# Documentation Review Checklist

**Document**: [filename]
**Reviewer**: [name]
**Date**: [date]

## Content
- [ ] Technically accurate
- [ ] Complete (no missing sections)
- [ ] Current (no outdated info)
- [ ] Examples tested and working

## Clarity
- [ ] Clear and concise
- [ ] Appropriate technical level
- [ ] Well-organized structure
- [ ] Good use of headings and lists

## Consistency
- [ ] Terminology consistent
- [ ] Formatting consistent
- [ ] Links working
- [ ] Cross-references accurate

## Usability
- [ ] Meets user needs
- [ ] Easy to navigate
- [ ] Code examples copy-paste ready
- [ ] Troubleshooting helpful

## Overall Rating
- [ ] Excellent - Ready to publish
- [ ] Good - Minor revisions needed
- [ ] Fair - Significant revisions needed
- [ ] Poor - Major rewrite needed

## Comments
[Specific feedback]
```

---

## Conclusion

This comprehensive documentation plan will transform BOLDGenotyper from a well-documented research tool to a publication-ready, user-friendly software package suitable for:

1. **Publication** - Methods sections, software papers, citations
2. **New Users** - Quick starts, tutorials, troubleshooting
3. **Developers** - API reference, architecture, contribution

**Next Steps**:
1. Review this plan
2. Prioritize phases based on publication timeline
3. Assign tasks and deadlines
4. Begin with Phase 1 (Critical Documentation)

**Questions?**
- Which phases are highest priority for your publication timeline?
- Do you need help with any specific sections?
- Should we focus on user docs or developer docs first?

---

**Document Version**: 1.0
**Last Updated**: 2025-11-17
**Status**: Ready for Review
