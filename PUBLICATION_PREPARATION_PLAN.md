# BOLDGenotyper Publication Preparation Plan

**Date**: 2025-11-25
**Status**: Comprehensive Review and Action Plan
**Target**: Publication-ready documentation and codebase cleanup
**Goal**: Prepare for (1) Bioinformatics methods manuscript and (2) Biological case study manuscript

---

## Executive Summary

This document provides a comprehensive plan for updating documentation and cleaning up the boldgenotyper repository in preparation for publication. The analysis covers all modules, CLI commands, documentation files, and output structures to ensure the package is ready for new users to install, understand, and use all functions.

### Current State
- ✅ **4 CLI commands** fully functional (main, compare, sweep, enrich)
- ✅ **22 Python modules** with comprehensive functionality
- ✅ **Strong README** with 1400+ lines of documentation
- ✅ **Extensive help text** for all CLI commands
- ✅ **Multiple output formats** (HTML, CSV, PDF, PNG, population genetics exports)
- ⚠️ **18 documentation files** (some redundant/outdated)
- ⚠️ **1 redundant module** (core.py)
- ⚠️ **Filename inconsistencies** in output (spaces vs underscores)
- ⚠️ **Missing documentation** for new features (compare, sweep, enrich, custom shapefiles)

---

## Part 1: Documentation Updates

### 1.1 README.md Updates

**Status**: Strong foundation, needs expansion for new features

#### Required Updates:

1. **Installation Section** (lines 56-130)
   - ✅ Already excellent
   - ➕ ADD: Note about new CLI commands
   - ➕ ADD: Section on verifying installation of all commands

2. **GOaS Section** (lines 132-226)
   - ✅ Already comprehensive
   - ➕ UPDATE: Add note about custom shapefile support via `boldgenotyper-enrich`
   - ➕ ADD: Examples for freshwater/terrestrial organism workflows

3. **Quick Start** (lines 228-286)
   - ✅ Good examples
   - ➕ ADD: Examples for new commands:
     ```bash
     # Parameter sweep to find optimal threshold
     boldgenotyper-sweep data/Sphyrna_lewini.tsv

     # Compare species vs family analysis
     boldgenotyper-compare --species-level Sphyrna_lewini_output/ \
                           --family-level Sphyrnidae_output/

     # Enrich with custom metadata or shapefile
     boldgenotyper-enrich Sphyrnidae_annotated.csv \
                          --custom-shp freshwater_basins.shp \
                          --shp-field basin_name
     ```

4. **Pipeline Overview** (lines 370-441)
   - ✅ Comprehensive 7-phase description
   - ✅ No changes needed for main pipeline
   - ➕ ADD: New section "Advanced Workflows" after line 441:
     - Parameter sweep workflow
     - Comparative analysis workflow
     - Metadata enrichment workflow

5. **Usage Guide** (lines 443-522)
   - ✅ Good basic examples
   - ➕ ADD: Section 6.5 "Advanced Commands" with:
     - `boldgenotyper-sweep` examples and interpretation
     - `boldgenotyper-compare` for contamination detection
     - `boldgenotyper-enrich` for custom geographic regions

6. **Parameter Reference** (lines 524-608)
   - ✅ Excellent detail for main command
   - ➕ ADD: Tables for each additional command:
     - **boldgenotyper-sweep parameters**
     - **boldgenotyper-compare parameters**
     - **boldgenotyper-enrich parameters**

7. **Output Files** (lines 610-733)
   - ✅ Excellent structure documentation
   - ➕ ADD: Subsection for sweep outputs:
     ```
     parameter_sweep/
     ├── sweep_summary.csv
     ├── threshold_stability.pdf
     ├── elbow_plot.pdf
     ├── group_membership_tracking.csv
     ├── recommendations.txt
     └── runs/{threshold}/
     ```
   - ➕ ADD: Subsection for compare outputs:
     ```
     comparative_analysis/
     ├── comparison_summary.csv
     ├── genotype_crosswalk.csv
     ├── sample_reassignments.csv
     └── methods_text.md
     ```
   - ➕ ADD: Subsection for enrich outputs:
     ```
     enriched_analysis/
     ├── {organism}_enriched.csv
     ├── geographic_summary.csv
     └── custom_visualizations/
     ```
   - ➕ ADD: Subsection for popgen exports (reference POPGEN_EXPORT_GUIDE.md):
     ```
     exports/
     ├── README.md
     ├── arlequin/
     ├── popart/
     ├── dnasp/
     └── generic/
     ```

8. **Biological Context** (lines 735-808)
   - ✅ Excellent content
   - ➕ ADD: Note about using parameter sweep to empirically determine thresholds
   - ➕ ADD: Reference to `boldgenotyper-sweep` as data-driven approach

9. **Assumptions and Limitations** (lines 810-935)
   - ⚠️ NEEDS UPDATE: Line 863 states "Marine Organisms Only (v1.0)"
   - ➕ UPDATE: Mention `boldgenotyper-enrich --custom-shp` now supports any shapefile
   - ➕ ADD: Examples of freshwater/terrestrial use cases with custom shapefiles
   - ✅ Keep limitation notes but update workaround sections

10. **Troubleshooting** (lines 937-1163)
    - ✅ Comprehensive
    - ➕ ADD: Troubleshooting for new commands
    - ➕ ADD: Common issues with custom shapefiles
    - ➕ ADD: Parameter sweep interpretation issues

11. **Citation** (lines 1215-1260)
    - ✅ Good structure
    - ➕ UPDATE: Finalize DOI once available
    - ➕ UPDATE: Journal name once submitted
    - ➕ ADD: How to cite specific components (sweep, compare, enrich)

### 1.2 New Documentation Files Needed

#### 1.2.1 CUSTOM_SHAPEFILES_GUIDE.md (NEW)

**Purpose**: Comprehensive guide for using custom shapefiles for any organism type

**Contents**:
- Overview of shapefile support
- Examples for freshwater basins (Salmonidae)
- Examples for terrestrial ecoregions
- Shapefile format requirements
- Field mapping instructions
- Complete workflow examples
- Troubleshooting shapefile issues

**Priority**: HIGH (enables non-marine use cases, key differentiator)

#### 1.2.2 PARAMETER_SWEEP_GUIDE.md (NEW)

**Purpose**: Detailed guide for threshold optimization

**Contents**:
- When to use parameter sweep
- Interpreting sweep outputs
- Understanding stability metrics
- Reading elbow plots
- Making threshold decisions
- Example workflows for different taxonomic groups
- Publication-ready methods text generation

**Priority**: HIGH (key analytical workflow)

#### 1.2.3 COMPARATIVE_ANALYSIS_GUIDE.md (NEW)

**Purpose**: Guide for species vs family contamination detection

**Contents**:
- Overview of comparative analysis approach
- When to use this workflow
- Interpreting contamination signals
- Using reassignment tables
- Example: Sphyrnidae case study
- Publication methods text

**Priority**: MEDIUM (specialized use case, but important for quality control)

#### 1.2.4 API_REFERENCE.md (NEW)

**Purpose**: Programmatic usage documentation

**Contents**:
- Importing modules directly
- Using pipeline functions in Python
- Custom workflows
- Integration with other tools
- Example notebooks

**Priority**: MEDIUM (for advanced users)

### 1.3 Updates to Existing Documentation

#### POPGEN_EXPORT_GUIDE.md
- ✅ Already excellent
- ➕ ADD: Integration examples with comparative analysis
- ➕ ADD: Note about export recommendations based on sweep results

#### CHANGELOG.md
- ➕ ADD: Detailed changelog for v0.1.0 with all features
- ➕ ADD: Section for each new command
- ➕ ADD: Breaking changes (if any)

#### CONTRIBUTING.md
- ✅ Check current state
- ➕ UPDATE: Add guidelines for documentation updates
- ➕ ADD: Examples of good pull requests

#### CODE_OF_CONDUCT.md
- ✅ Verify completeness
- No changes expected

### 1.4 Documentation to Archive/Remove

**Move to `docs/archive/` directory** (not delete, for historical reference):

1. **DOCUMENTATION_PLAN.md** - Superseded by this plan
2. **PHASE1_IMPLEMENTATION_SUMMARY.md** - Internal development doc
3. **SPECIES_COMPOSITION_PLAN.md** - Internal development doc
4. **BOLDGenotyper_Enhancement_Specifications.md** - Internal development doc
5. **PUBLICATION_READY_SUMMARY.md** - Superseded by this plan
6. **ZENODO_SETUP.md** - Should be in docs/ or completed/removed
7. **ERROR_HANDLING_GUIDE.md** - Merge into CONTRIBUTING.md or remove
8. **QUICKSTART_TESTING.md** - Internal development doc
9. **TESTING_CLI.md** - Internal development doc
10. **TESTING_PLAN.md** - Internal development doc

**Consolidate**:
- **QUICKSTART.md** - Merge best parts into README Quick Start section
- **BOLDGenotyper.md** - Remove if duplicate of README, or rename to something specific

**Result**: Clean root directory with only essential user-facing documentation

### 1.5 Docstring Updates

**Action**: Audit all modules for completeness

**Priority modules** (ensure 100% coverage):
1. `cli.py` - All functions documented with examples
2. `comparative_analysis.py` - Full documentation of workflow
3. `parameter_sweep.py` - Comprehensive documentation of metrics
4. `geographic_enhancement.py` - Document custom shapefile usage
5. `metadata_enrichment.py` - Full API documentation
6. `popgen_export.py` - Format specifications

**Standard format** (NumPy style):
```python
def function_name(param1: type, param2: type) -> return_type:
    """
    Brief one-line description.

    Longer description with biological context, use cases, and
    important notes.

    Parameters
    ----------
    param1 : type
        Description of param1, including biological interpretation
    param2 : type
        Description of param2, including default behavior

    Returns
    -------
    return_type
        Description of return value

    Raises
    ------
    ExceptionType
        When this exception is raised

    Examples
    --------
    >>> result = function_name(value1, value2)
    >>> print(result)

    Notes
    -----
    - Important biological considerations
    - Performance notes for large datasets
    - Related functions

    See Also
    --------
    related_function : Related functionality
    """
```

---

## Part 2: Repository Cleanup

### 2.1 Code Cleanup

#### 2.1.1 Redundant Module: core.py

**Issue**: `boldgenotyper/core.py` exists but is not used in current CLI implementation

**Evidence**:
- `cli.py` has its own `run_pipeline()` function (lines 98+)
- `core.py` only imported in `__init__.py` and one test file (`test_core.py`)
- All pipeline orchestration now in `cli.py`

**Action**:
- **Option 1 (RECOMMENDED)**: Remove `core.py` and `test_core.py`
  - Update `__init__.py` to remove core import
  - Clean, no redundancy

- **Option 2**: Refactor `cli.py` to use `core.run_pipeline()`
  - More modular design
  - But adds complexity for minimal benefit

**Recommendation**: Remove core.py (Option 1)

**Files to modify**:
1. Delete: `boldgenotyper/core.py`
2. Delete: `tests/test_core.py`
3. Update: `boldgenotyper/__init__.py` (remove core import)

#### 2.1.2 Import Cleanup

**Action**: Verify all imports in `__init__.py` are necessary

Current `__init__.py` exports:
```python
__all__ = [
    "core",           # ← REMOVE (module being deleted)
    "dereplication",  # ✅ Keep
    "metadata",       # ✅ Keep
    "geographic",     # ✅ Keep
    "visualization",  # ✅ Keep
    "phylogenetics",  # ✅ Keep
    "utils",          # ✅ Keep
]
```

**Add to `__all__`** (for programmatic access):
```python
__all__ = [
    "dereplication",
    "metadata",
    "geographic",
    "visualization",
    "phylogenetics",
    "utils",
    "comparative_analysis",  # NEW
    "parameter_sweep",       # NEW
    "metadata_enrichment",   # NEW
    "popgen_export",         # NEW
]
```

#### 2.1.3 TODO/FIXME Comments

**Found**:
- `plot_export.py:506` - "TODO" for example modifications

**Action**:
- Review and either implement or remove TODO comments
- Document known limitations in docstrings instead

### 2.2 Output File Cleanup

#### 2.2.1 Filename Inconsistency

**Issue**: Output files have inconsistent naming with spaces vs underscores

**Evidence**:
```
data/Sphyrna_lewini/genotype_assignments/
├── Sphyrna lewini_diagnostics.csv       # ← Spaces
├── Sphyrna_lewini_diagnostics.csv       # ← Underscores (duplicate?)
├── Sphyrna lewini_identity_distribution.pdf
├── Sphyrna_lewini_identity_distribution.pdf
```

**Root cause**: Organism name extraction handling spaces inconsistently

**Action**:
1. **Standardize** on underscores throughout all file naming
2. **Update** `cli.py` organism name extraction to always replace spaces with underscores
3. **Propagate** this change to all modules that generate output filenames
4. **Test** with species names containing spaces

**Modules to update**:
- `cli.py`: `extract_organism_from_path()` - ensure underscore conversion
- `reports.py`: All file naming
- `visualization.py`: All file naming
- `phylogenetics.py`: Tree file naming
- `quality_control.py`: Report file naming
- `divergence_analysis.py`: Output file naming
- `comparative_analysis.py`: Output file naming
- `parameter_sweep.py`: Output file naming

**Testing needed**:
```bash
# Test cases
boldgenotyper "data/Species name with spaces.tsv"
boldgenotyper "data/Sphyrna_lewini.tsv"
boldgenotyper "data/Genus-species-subspecies.tsv"
```

#### 2.2.2 Directory Structure Optimization

**Current structure**: Good, but could be clearer

**Proposed improvements**:

1. **Separate intermediate files** more clearly:
   ```
   {organism}_output/
   ├── intermediate/           # Hidden from typical user view
   │   ├── 01_parsed_metadata.tsv
   │   ├── 02_filtered_metadata.tsv
   │   ├── dereplication/
   │   ├── genotype_assignments/
   │   ├── geographic/
   │   └── phylogenetic/
   ```

2. **Consider** adding `.gitkeep` files to preserve directory structure

3. **Add** top-level `RESULTS_GUIDE.md` in each output directory explaining structure

#### 2.2.3 Cleanup Empty Directories

**Action**: Implement cleanup function more aggressively

**Current**: `remove_empty_directories()` exists in `cli.py`

**Enhancement**:
- Run at end of pipeline
- Log removed directories
- Ensure it runs for all CLI commands

### 2.3 Test Data Cleanup

**Current test data** in `data/`:
```
data/
├── Cyprinodontidae/       # Test run
├── Cyprinodontidae.tsv
├── Rhizoprionodon/        # Test run
├── Rhizoprionodon.tsv
├── Salmonidae/            # Test run
├── Salmonidae.tsv
├── Sphyrna_lewini/        # Example (keep)
├── Sphyrna_lewini_scallopedhammerhead.tsv
├── Sphyrna_outgroup.tsv
├── Sphyrnidae/            # Example (keep)
├── Sphyrnidae.tsv
└── Sphyrnidae_compare/    # Example (keep)
```

**Actions**:

1. **Decide on example datasets**:
   - **Keep**: Sphyrnidae (family-level example)
   - **Keep**: Sphyrna_lewini (species-level example)
   - **Keep**: Sphyrnidae_compare (comparative analysis example)

2. **Move test runs** to `examples/` or remove:
   - Cyprinodontidae → `examples/freshwater/` (if using custom shapefile)
   - Salmonidae → `examples/freshwater/`
   - Rhizoprionodon → Remove or `examples/elasmobranchs/`

3. **Create** `data/README.md` explaining each dataset:
   ```markdown
   # Example Datasets

   ## Sphyrnidae (Family-level)
   - **Organism**: Sphyrnidae (hammerhead sharks)
   - **Samples**: ~1000
   - **Purpose**: Demonstrates family-level analysis
   - **Run time**: ~5 minutes

   ## Sphyrna_lewini (Species-level)
   - **Organism**: Sphyrna lewini (scalloped hammerhead)
   - **Samples**: ~600
   - **Purpose**: Demonstrates species-level analysis
   - **Run time**: ~3 minutes

   ## Sphyrnidae_compare
   - **Purpose**: Output from comparative analysis
   - **Demonstrates**: Contamination detection workflow
   ```

### 2.4 Remove Development Files

**Files to remove or relocate**:

1. **Root directory cleanup**:
   - Move `.DS_Store` files to `.gitignore`
   - Remove any `__pycache__` directories
   - Remove any `.pyc` files

2. **Archive internal docs** (create `docs/archive/`):
   - Move planning documents listed in section 1.4

3. **Organize docs**:
   ```
   docs/
   ├── archive/              # Historical planning docs
   ├── api/                  # API reference (future)
   ├── tutorials/            # Step-by-step guides (future)
   └── examples/             # Example notebooks (future)
   ```

---

## Part 3: Propagation of Changes

### 3.1 Filename Changes Propagation

**When standardizing filenames**, ensure changes propagate to:

1. **HTML report generation** (`reports.py`):
   - Update all file path references
   - Update plot loading (JSON data files)
   - Test interactive features

2. **Plot data export** (`plot_export.py`):
   - Ensure R scripts use correct filenames
   - Update regeneration scripts

3. **Comparative analysis** (`comparative_analysis.py`):
   - Update file loading logic
   - Ensure crosswalk tables use correct names

4. **Parameter sweep** (`parameter_sweep.py`):
   - Ensure output parsing uses correct filenames
   - Update summary report generation

5. **Tests**:
   - Update all test expectations
   - Verify fixture paths
   - Re-run full test suite

### 3.2 Testing Checklist

**Before finalizing**:

- [ ] Run full test suite: `pytest tests/ -v`
- [ ] Test main pipeline with multiple organisms
- [ ] Test all 4 CLI commands
- [ ] Test with organism names containing spaces
- [ ] Test with organism names containing hyphens
- [ ] Test with custom shapefiles (freshwater example)
- [ ] Verify all output files have consistent naming
- [ ] Check HTML reports load correctly
- [ ] Verify population genetics exports work
- [ ] Test parameter sweep workflow
- [ ] Test comparative analysis workflow
- [ ] Test metadata enrichment workflow

### 3.3 Documentation Updates Checklist

- [ ] Update README.md (all sections listed in 1.1)
- [ ] Create CUSTOM_SHAPEFILES_GUIDE.md
- [ ] Create PARAMETER_SWEEP_GUIDE.md
- [ ] Create COMPARATIVE_ANALYSIS_GUIDE.md
- [ ] Update POPGEN_EXPORT_GUIDE.md
- [ ] Update CHANGELOG.md
- [ ] Archive internal development docs
- [ ] Create data/README.md
- [ ] Add docstrings to all public functions
- [ ] Generate API documentation (Sphinx)

---

## Part 4: Priority Action Items

### Phase 1: Critical (Do First)

1. **Fix filename inconsistency** (HIGH IMPACT)
   - Standardize organism name handling
   - Update all file generation
   - Test thoroughly
   - **Estimated time**: 2-3 hours

2. **Update README with new features** (HIGH IMPACT)
   - Add examples for all 4 CLI commands
   - Document parameter sweep
   - Document comparative analysis
   - Document custom shapefiles
   - **Estimated time**: 3-4 hours

3. **Create CUSTOM_SHAPEFILES_GUIDE.md** (HIGH IMPACT)
   - Enables non-marine use cases
   - Key differentiator from v1.0 limitations
   - **Estimated time**: 2-3 hours

4. **Archive internal docs** (CLEAN UP)
   - Create docs/archive/
   - Move development docs
   - Clean root directory
   - **Estimated time**: 30 minutes

### Phase 2: Important (Do Second)

5. **Remove core.py** (CODE CLEANUP)
   - Delete redundant module
   - Update imports
   - Update tests
   - **Estimated time**: 30 minutes

6. **Create PARAMETER_SWEEP_GUIDE.md** (DOCUMENTATION)
   - Document threshold optimization
   - Interpretation guide
   - **Estimated time**: 2 hours

7. **Update all module docstrings** (DOCUMENTATION)
   - Ensure 100% coverage
   - Add examples
   - **Estimated time**: 3-4 hours

8. **Clean up test data** (ORGANIZATION)
   - Decide on examples
   - Create data/README.md
   - **Estimated time**: 1 hour

### Phase 3: Nice to Have (Do Third)

9. **Create COMPARATIVE_ANALYSIS_GUIDE.md** (DOCUMENTATION)
   - Document workflow
   - Interpretation guide
   - **Estimated time**: 2 hours

10. **Generate API documentation** (DOCUMENTATION)
    - Set up Sphinx
    - Generate HTML docs
    - **Estimated time**: 2-3 hours

11. **Create tutorial notebooks** (EXAMPLES)
    - Jupyter notebooks for common workflows
    - **Estimated time**: 4-5 hours

---

## Part 5: Publication-Specific Requirements

### 5.1 Methods Manuscript Requirements

**For bioinformatics methods paper**, ensure:

1. **Software availability**:
   - [ ] GitHub repository public and complete
   - [ ] PyPI package published (pip install boldgenotyper)
   - [ ] Conda package available (conda install -c bioconda boldgenotyper)
   - [ ] Docker image available (optional)
   - [ ] Zenodo DOI assigned

2. **Documentation completeness**:
   - [ ] Installation instructions (multiple platforms)
   - [ ] Quick start guide
   - [ ] Full parameter documentation
   - [ ] Example datasets with expected outputs
   - [ ] API reference
   - [ ] Troubleshooting guide

3. **Reproducibility**:
   - [ ] Version pinning in environment.yml
   - [ ] Example workflows with exact commands
   - [ ] Test data with checksums
   - [ ] Continuous integration (GitHub Actions)

4. **Benchmarking** (if required):
   - [ ] Performance metrics (runtime, memory)
   - [ ] Comparison with alternative tools (if applicable)
   - [ ] Scalability testing

### 5.2 Case Study Manuscript Requirements

**For Sphyrnidae biological paper**, ensure:

1. **Reproducible analysis**:
   - [ ] All commands documented
   - [ ] Parameter choices justified
   - [ ] Data availability statement
   - [ ] Supplementary methods section

2. **Methods text generation**:
   - [ ] Comparative analysis generates methods text
   - [ ] Parameter sweep generates threshold justification
   - [ ] All visualizations have figure legends

3. **Supplementary materials**:
   - [ ] Table S1: Sample metadata
   - [ ] Table S2: Genotype assignments
   - [ ] Table S3: Geographic distribution
   - [ ] Table S4: Reassignments from comparative analysis
   - [ ] Figure S1-SX: Supplementary figures

---

## Part 6: Timeline Estimate

### Week 1: Critical Fixes
- Day 1-2: Fix filename inconsistency, test thoroughly
- Day 3-4: Update README with all new features
- Day 5: Create CUSTOM_SHAPEFILES_GUIDE.md

### Week 2: Code Cleanup and Documentation
- Day 1: Remove core.py, clean up imports
- Day 2-3: Create PARAMETER_SWEEP_GUIDE.md
- Day 4-5: Update all module docstrings

### Week 3: Examples and Polish
- Day 1: Clean up test data, create data/README.md
- Day 2: Create COMPARATIVE_ANALYSIS_GUIDE.md
- Day 3: Archive internal docs, organize root directory
- Day 4-5: Testing, bug fixes, final polish

### Week 4: Publication Preparation
- Day 1-2: Generate API documentation
- Day 3: Create tutorial notebooks
- Day 4: Prepare for PyPI/Conda release
- Day 5: Final review and testing

**Total estimated time**: 3-4 weeks of focused work

---

## Part 7: Success Criteria

### Documentation Success Criteria

- [ ] New user can install and run basic analysis in <15 minutes
- [ ] All CLI commands documented with examples
- [ ] All parameters explained with biological context
- [ ] All output files documented
- [ ] Custom shapefile workflow documented with examples
- [ ] Parameter sweep interpretation guide available
- [ ] Comparative analysis workflow documented
- [ ] API reference available for programmatic use
- [ ] Root directory clean with only essential documentation
- [ ] All internal dev docs archived

### Code Success Criteria

- [ ] No redundant modules
- [ ] All output filenames consistent (underscores only)
- [ ] All imports necessary and documented
- [ ] All public functions have docstrings
- [ ] Full test coverage (>80%)
- [ ] All tests passing
- [ ] No TODO/FIXME in production code

### Repository Success Criteria

- [ ] Clean root directory structure
- [ ] Clear examples with documentation
- [ ] Organized docs/ directory
- [ ] Updated CHANGELOG
- [ ] Clear CONTRIBUTING guidelines
- [ ] All licenses and citations correct

---

## Appendix A: Files Audit

### Python Modules (boldgenotyper/)

| Module | Status | Used In | Notes |
|--------|--------|---------|-------|
| `__init__.py` | ✅ Active | Package initialization | Update exports |
| `cli.py` | ✅ Active | Main entry point | Main pipeline |
| `comparative_analysis.py` | ✅ Active | boldgenotyper-compare | Document |
| `parameter_sweep.py` | ✅ Active | boldgenotyper-sweep | Document |
| `metadata_enrichment.py` | ✅ Active | boldgenotyper-enrich | Document |
| `geographic_enhancement.py` | ✅ Active | boldgenotyper-enrich | Document custom shapefiles |
| `popgen_export.py` | ✅ Active | Main pipeline | Already documented |
| `dereplication.py` | ✅ Active | Main pipeline | Core functionality |
| `genotype_assignment.py` | ✅ Active | Main pipeline | Core functionality |
| `geographic.py` | ✅ Active | Main pipeline | Core functionality |
| `phylogenetics.py` | ✅ Active | Main pipeline | Core functionality |
| `visualization.py` | ✅ Active | Main pipeline | Core functionality |
| `reports.py` | ✅ Active | Main pipeline | Core functionality |
| `metadata.py` | ✅ Active | Main pipeline | Core functionality |
| `quality_control.py` | ✅ Active | Main pipeline | Core functionality |
| `divergence_analysis.py` | ✅ Active | Main pipeline | Core functionality |
| `cluster_diagnostics.py` | ✅ Active | Main pipeline | Core functionality |
| `plot_export.py` | ✅ Active | Main pipeline | Core functionality |
| `config.py` | ✅ Active | All modules | Configuration |
| `utils.py` | ✅ Active | All modules | Utilities |
| `goas_downloader.py` | ✅ Active | Setup | GOaS download |
| **`core.py`** | ⚠️ **REDUNDANT** | Only test_core.py | **REMOVE** |

### Documentation Files (Root)

| File | Status | Action |
|------|--------|--------|
| `README.md` | ✅ Keep | **UPDATE** (add new features) |
| `CHANGELOG.md` | ✅ Keep | **UPDATE** |
| `CODE_OF_CONDUCT.md` | ✅ Keep | Verify |
| `CONTRIBUTING.md` | ✅ Keep | Update |
| `paper.md` | ✅ Keep | Update for JOSS submission |
| `POPGEN_EXPORT_GUIDE.md` | ✅ Keep | Minor updates |
| `QUICKSTART.md` | ⚠️ Review | Merge into README or keep? |
| `DOCUMENTATION_PLAN.md` | ⚠️ Archive | Move to docs/archive/ |
| `PUBLICATION_READY_SUMMARY.md` | ⚠️ Archive | Move to docs/archive/ |
| `ZENODO_SETUP.md` | ⚠️ Archive | Complete or archive |
| `ERROR_HANDLING_GUIDE.md` | ⚠️ Archive | Merge into CONTRIBUTING |
| `QUICKSTART_TESTING.md` | ⚠️ Archive | Internal dev doc |
| `TESTING_CLI.md` | ⚠️ Archive | Internal dev doc |
| `TESTING_PLAN.md` | ⚠️ Archive | Internal dev doc |
| `PHASE1_IMPLEMENTATION_SUMMARY.md` | ⚠️ Archive | Internal dev doc |
| `SPECIES_COMPOSITION_PLAN.md` | ⚠️ Archive | Internal dev doc |
| `BOLDGenotyper_Enhancement_Specifications.md` | ⚠️ Archive | Internal dev doc |
| `BOLDGenotyper.md` | ⚠️ Review | Duplicate of README? |

### New Documentation Needed

| File | Priority | Estimated Time |
|------|----------|----------------|
| `CUSTOM_SHAPEFILES_GUIDE.md` | HIGH | 2-3 hours |
| `PARAMETER_SWEEP_GUIDE.md` | HIGH | 2 hours |
| `COMPARATIVE_ANALYSIS_GUIDE.md` | MEDIUM | 2 hours |
| `API_REFERENCE.md` | MEDIUM | 3-4 hours (with Sphinx) |
| `data/README.md` | MEDIUM | 1 hour |
| `docs/archive/README.md` | LOW | 30 min |

---

## Appendix B: Command Reference Summary

### Current CLI Commands

1. **boldgenotyper** - Main pipeline
   - Input: BOLD TSV file
   - Output: Complete genotyping analysis
   - Flags: --build-tree, --no-geo, --export-format, --export-plot-data
   - Parameters: --clustering-threshold, --similarity-threshold, --tie-margin, --tie-threshold

2. **boldgenotyper-sweep** - Parameter optimization
   - Input: BOLD TSV file
   - Output: Threshold recommendations
   - Key outputs: sweep_summary.csv, threshold_stability.pdf, recommendations.txt

3. **boldgenotyper-compare** - Contamination detection
   - Input: Two analysis directories (species-level, family-level)
   - Output: Comparison metrics, reassignment table
   - Key outputs: comparison_summary.csv, genotype_crosswalk.csv, sample_reassignments.csv

4. **boldgenotyper-enrich** - Metadata/geographic enhancement
   - Input: Annotated CSV from main pipeline
   - Options: --add-metadata, --custom-shp, --recalculate-geography
   - Output: Enriched dataset with custom annotations

---

## Next Steps

1. **Review this plan** with collaborators
2. **Prioritize** sections based on publication timeline
3. **Assign tasks** if working with team
4. **Create GitHub issues** for each major task
5. **Set milestones** for publication deadlines
6. **Begin with Phase 1** critical fixes

---

**Document prepared by**: Claude Code Analysis
**Date**: 2025-11-25
**For**: BOLDGenotyper Publication Preparation
**Contact**: steph.smith@unc.edu
