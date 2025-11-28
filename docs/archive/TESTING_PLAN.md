# BOLDGenotyper - Automated Testing Plan for Publication

**Document Version**: 1.0
**Date**: 2025-11-18
**Status**: Planning Phase
**Target**: JOSS Publication & PLOS ONE Case Study

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current Test Coverage Analysis](#current-test-coverage-analysis)
3. [Test Suite Architecture](#test-suite-architecture)
4. [Priority Test Development Plan](#priority-test-development-plan)
5. [Testing Standards & Best Practices](#testing-standards--best-practices)
6. [Documentation Requirements](#documentation-requirements)
7. [Continuous Integration Setup](#continuous-integration-setup)
8. [Timeline and Resources](#timeline-and-resources)
9. [Success Criteria for Publication](#success-criteria-for-publication)
10. [Appendix: Module Function Inventory](#appendix-module-function-inventory)

---

## Executive Summary

### Goals

Create a **minimal yet comprehensive automated test suite** for BOLDGenotyper that:
- Tests all public module functions
- Achieves 85-90% code coverage
- Meets JOSS publication standards
- Is maintainable and well-documented
- Runs efficiently in CI/CD pipelines

### Current Status

- **Total Codebase**: ~9,571 lines of Python code across 13 modules
- **Existing Tests**: 17 test files with ~100+ test cases
- **Current Coverage**: ~60-70% (estimated)
- **Test Framework**: Python `unittest` with mocking support

### Gaps to Address

**Priority 1 (Critical for JOSS):**
1. Core pipeline orchestration (`core.py`) - currently stub
2. Report generation (`reports.py`) - no tests
3. Visualization functions (`visualization.py`) - minimal tests
4. CLI integration (`cli.py`) - minimal tests
5. Configuration management (`config.py`) - partial tests

**Priority 2 (Important):**
6. Utility functions (`utils.py`) - partial tests
7. Phylogenetics (`phylogenetics.py`) - partial tests
8. GOaS downloader (`goas_downloader.py`) - no tests

### Target Metrics

- **Code Coverage**: 85-90%
- **Test Files**: Add 5-8 new comprehensive test files
- **Test Cases**: Add ~150-200 new test cases
- **Documentation**: Each test file with module-level docstring
- **CI/CD**: Automated testing on all commits and PRs

---

## Current Test Coverage Analysis

### Well-Tested Modules (75-95% coverage)

#### 1. `dereplication.py` ✅
- **Test File**: `test_dereplication.py`
- **Coverage**: ~90%
- **Test Classes**: 8 classes covering:
  - External tool detection
  - Distance calculation
  - Consensus generation
  - Clustering
  - MAFFT/trimAl integration
  - Full workflow integration
  - Edge cases
- **Status**: Excellent coverage

#### 2. `genotype_assignment.py` ✅
- **Test File**: `test_genotype_assignment.py`
- **Coverage**: ~95%
- **Test Classes**: 10 classes covering:
  - Edit distance calculation
  - Identity metrics
  - ProcessID extraction
  - Best match finding with ties
  - FASTA reading
  - Parallel processing
  - Missing data handling
  - Diagnostics generation
  - Full workflow
  - Error handling
- **Status**: Excellent coverage

#### 3. `metadata.py` ✅
- **Test File**: `test_metadata.py`
- **Coverage**: ~85%
- **Test Classes**: 9 classes covering:
  - BOLD TSV parsing
  - Column validation
  - Coordinate extraction and parsing
  - Coordinate filtering (centroids)
  - QC statistics
  - Edit distance
  - Metadata validation
  - Integration tests
  - Edge cases
- **Status**: Good coverage

#### 4. `geographic.py` ✅
- **Test File**: `test_geographic.py`
- **Coverage**: ~80%
- **Test Classes**: 9 classes covering:
  - GeoPandas availability
  - GOaS data loading
  - Point GeoDataFrame creation
  - Ocean basin assignment
  - Basin counts and statistics
  - Marine coordinate validation
  - Marine filtering
  - Edge cases
- **Status**: Good coverage

#### 5. `identity_calculation.py` ✅
- **Test File**: `test_identity_calculation.py`
- **Coverage**: ~90%
- **Test Classes**: 6 classes covering:
  - CIGAR parsing
  - CIGAR-based identity
  - Best match finding
  - Config validation
  - Method comparison
  - Edge cases
- **Status**: Excellent coverage

### Partially-Tested Modules (30-60% coverage)

#### 6. `phylogenetics.py` ⚠️
- **Test Files**: `test_phylo_viz_pipeline.py`, `test_tree_scaling.py`
- **Coverage**: ~40%
- **Current Tests**: Integration tests only
- **Missing Tests**:
  - Unit tests for `run_mafft_alignment()`
  - Unit tests for `run_fasttree()`
  - `run_phyml()` (stub function)
  - `midpoint_root_tree()` (stub function)
  - Outgroup handling
  - Tree relabeling
  - Bootstrap support handling
  - Error handling
- **Status**: Needs dedicated unit test file

#### 7. `visualization.py` ⚠️
- **Test File**: `test_viz_warnings.py`
- **Coverage**: ~15%
- **Current Tests**: Only warning tests
- **Missing Tests**: All plotting functions untested
  - `plot_distribution_map()`
  - `plot_ocean_basin_abundance()`
  - `plot_distribution_map_faceted()`
  - `plot_ocean_basin_abundance_faceted()`
  - `plot_phylogenetic_tree()`
  - `plot_identity_distribution()`
  - `plot_identity_by_status()`
  - `plot_assignment_status()`
- **Status**: Needs comprehensive mock-based tests

#### 8. `config.py` ⚠️
- **Test Coverage**: Via `test_identity_calculation.py` only
- **Coverage**: ~30%
- **Current Tests**: Only dataclass usage
- **Missing Tests**:
  - `load_config_from_file()` (YAML/JSON)
  - `load_config_from_env()` (environment variables)
  - `validate_config()` (validation logic)
  - `create_config_template()` (template generation)
  - Nested parameter updates
  - Path serialization/deserialization
  - Type conversion and validation
- **Status**: Needs dedicated test file

#### 9. `utils.py` ⚠️
- **Test Coverage**: Indirect via other tests
- **Coverage**: ~40%
- **Current Tests**: Minimal
- **Missing Tests**:
  - Logging setup and function call logging
  - External tool checking with version comparison
  - `ProgressTracker` class
  - `assign_consensus_taxonomy()`
  - `pick_final_group_taxon()`
  - Sequence validation functions
  - Sequence statistics
  - File I/O helpers
  - Time/format utilities
- **Status**: Needs dedicated test file

#### 10. `cli.py` ⚠️
- **Test File**: `test_pipeline_script.py`
- **Coverage**: ~20%
- **Current Tests**: Minimal integration only
- **Missing Tests**:
  - `extract_organism_from_path()`
  - `setup_directories()`
  - Argument parsing
  - Error handling
  - Help text generation
  - Full pipeline integration
- **Status**: Needs comprehensive CLI test file

### Untested Modules (0-10% coverage)

#### 11. `core.py` ❌
- **Test File**: None
- **Coverage**: 0%
- **Status**: **CRITICAL - Main orchestration function is a stub**
- **Issue**: `run_pipeline()` has docstring but no implementation
- **Required**:
  - Implement `run_pipeline()` function
  - Create `test_core.py` with integration tests
  - Test all pipeline phases
  - Test error handling and recovery

#### 12. `reports.py` ❌
- **Test File**: None
- **Coverage**: 0%
- **Functions to Test**:
  - `generate_taxonomy_conflicts_report()`
  - `generate_assignment_summary()`
  - `generate_consensus_characterization()`
  - `generate_sequence_quality_metrics()`
- **Required**: Create `test_reports.py` with comprehensive tests
- **Status**: **CRITICAL for JOSS - No QC report testing**

#### 13. `goas_downloader.py` ❌
- **Test File**: None
- **Coverage**: 0%
- **Functions to Test**:
  - `calculate_md5()`
  - `download_file()`
  - `extract_zip()`
  - `verify_files()`
  - `create_citation_file()`
  - `setup_goas()`
- **Required**: Create `test_goas_downloader.py` with mock-based tests
- **Status**: Lower priority but should have basic tests

---

## Test Suite Architecture

### Design Principles

1. **Modularity**: One test file per module with clear organization
2. **Mock External Dependencies**: Mock MAFFT, trimAl, FastTree, downloads, etc.
3. **Fast Execution**: Unit tests should run in <5 minutes total
4. **Comprehensive Coverage**: Cover happy paths, edge cases, and error conditions
5. **Readable**: Clear test names, good docstrings, organized test classes
6. **Maintainable**: DRY principles, reusable fixtures and helpers

### Test Organization

```
tests/
├── conftest.py                      # Pytest fixtures (to be added)
├── test_data/                       # Test data directory
│   ├── sample.fasta
│   ├── sample.tsv
│   └── ...
├── test_dereplication.py           # ✅ Existing - Good
├── test_genotype_assignment.py     # ✅ Existing - Good
├── test_metadata.py                # ✅ Existing - Good
├── test_geographic.py              # ✅ Existing - Good
├── test_identity_calculation.py    # ✅ Existing - Good
├── test_core.py                    # ❌ NEW - Critical
├── test_reports.py                 # ❌ NEW - Critical
├── test_visualization.py           # ❌ NEW - Critical
├── test_cli.py                     # ❌ NEW - Critical
├── test_config.py                  # ❌ NEW - Important
├── test_utils.py                   # ❌ NEW - Important
├── test_phylogenetics.py           # ❌ NEW - Important
├── test_goas_downloader.py         # ❌ NEW - Nice to have
└── test_integration.py             # ❌ NEW - End-to-end tests
```

### Test Categories

#### 1. Unit Tests (Primary Focus)
- Test individual functions in isolation
- Use mocking for external dependencies
- Fast execution (<1 second per test)
- Examples: distance calculation, coordinate parsing, file validation

#### 2. Integration Tests
- Test module interactions
- Use real small test data
- Moderate execution time (<10 seconds per test)
- Examples: full dereplication workflow, assignment pipeline

#### 3. End-to-End Tests
- Test complete pipeline from input to output
- Use realistic test datasets
- Slower execution (<60 seconds per test)
- Examples: full pipeline with all phases

#### 4. Mock-Based Tests (for Visualization, External Tools)
- Mock matplotlib, cartopy for visualization
- Mock subprocess calls for external tools
- Verify correct calls without actual execution
- Examples: plot functions, MAFFT/trimAl calls

---

## Priority Test Development Plan

### Phase 1: Critical Gaps (Week 1-2) - REQUIRED FOR JOSS

#### Task 1.1: Implement and Test `core.py` 🔴 CRITICAL
**Priority**: Highest
**Effort**: 8-12 hours

**Implementation Steps:**
1. **Implement `core.py:run_pipeline()` function**
   - Main orchestration logic
   - Call dereplication, assignment, geographic, visualization, phylogenetics
   - Error handling and logging
   - Progress tracking
   - Result aggregation

2. **Create `test_core.py`**
   ```python
   class TestCorePipeline(unittest.TestCase):
       """Test main pipeline orchestration."""

       def test_run_pipeline_minimal_config(self):
           """Test pipeline with minimal configuration."""

       def test_run_pipeline_full_config(self):
           """Test pipeline with all features enabled."""

       def test_run_pipeline_with_phylogeny(self):
           """Test pipeline with phylogenetic tree building."""

       def test_run_pipeline_error_handling(self):
           """Test pipeline error recovery."""

       def test_run_pipeline_creates_outputs(self):
           """Test that all expected output files are created."""

       def test_run_pipeline_with_missing_dependencies(self):
           """Test graceful degradation when optional tools missing."""
   ```

**Deliverables:**
- Implemented `core.py:run_pipeline()`
- `test_core.py` with 15-20 test cases
- Integration tests with real small dataset

#### Task 1.2: Create `test_reports.py` 🔴 CRITICAL
**Priority**: Highest
**Effort**: 6-8 hours

**Test Coverage:**
```python
class TestTaxonomyConflictsReport(unittest.TestCase):
    """Test taxonomy conflict detection and reporting."""

    def test_generate_taxonomy_conflicts_report_no_conflicts(self):
        """Test report with no conflicts."""

    def test_generate_taxonomy_conflicts_report_with_conflicts(self):
        """Test report with taxonomy mismatches."""

    def test_generate_taxonomy_conflicts_report_empty_data(self):
        """Test report with empty input."""


class TestAssignmentSummary(unittest.TestCase):
    """Test assignment summary report generation."""

    def test_generate_assignment_summary_complete_data(self):
        """Test summary with complete assignment data."""

    def test_generate_assignment_summary_with_ties(self):
        """Test summary including tie statistics."""

    def test_generate_assignment_summary_below_threshold(self):
        """Test summary with low-confidence assignments."""


class TestConsensusCharacterization(unittest.TestCase):
    """Test consensus group characterization."""

    def test_generate_consensus_characterization_complete(self):
        """Test characterization with all data."""

    def test_generate_consensus_characterization_missing_taxonomy(self):
        """Test handling of missing taxonomy."""


class TestSequenceQualityMetrics(unittest.TestCase):
    """Test sequence quality report generation."""

    def test_generate_sequence_quality_metrics_good_quality(self):
        """Test metrics for high-quality sequences."""

    def test_generate_sequence_quality_metrics_poor_quality(self):
        """Test metrics for low-quality sequences (high N content)."""

    def test_generate_sequence_quality_metrics_length_distribution(self):
        """Test sequence length analysis."""
```

**Deliverables:**
- `test_reports.py` with 20-25 test cases
- CSV output validation
- Edge case handling (empty data, missing columns)

#### Task 1.3: Create `test_visualization.py` 🔴 CRITICAL
**Priority**: Highest
**Effort**: 10-12 hours

**Test Strategy**: Mock matplotlib and cartopy to avoid image generation

```python
class TestColorPalette(unittest.TestCase):
    """Test genotype color generation."""

    def test_get_genotype_colors_small_set(self):
        """Test color generation for small genotype set."""

    def test_get_genotype_colors_large_set(self):
        """Test color generation for large genotype set."""

    def test_get_genotype_colors_colorblind_friendly(self):
        """Test that colors are distinguishable."""


class TestDistributionMap(unittest.TestCase):
    """Test global distribution map generation."""

    @patch('matplotlib.pyplot.savefig')
    @patch('cartopy.crs.PlateCarree')
    def test_plot_distribution_map_basic(self, mock_crs, mock_savefig):
        """Test basic map generation."""

    def test_plot_distribution_map_no_data(self):
        """Test map with no data points."""

    def test_plot_distribution_map_single_genotype(self):
        """Test map with single genotype."""

    def test_plot_distribution_map_no_cartopy(self):
        """Test graceful degradation without cartopy."""


class TestOceanBasinAbundance(unittest.TestCase):
    """Test ocean basin abundance plots."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_ocean_basin_abundance_basic(self, mock_savefig):
        """Test basic stacked bar chart."""

    def test_plot_ocean_basin_abundance_single_basin(self):
        """Test chart with single basin."""

    def test_plot_ocean_basin_abundance_empty_basins(self):
        """Test handling of empty basins."""


class TestFacetedPlots(unittest.TestCase):
    """Test multi-panel faceted plots."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_distribution_map_faceted(self, mock_savefig):
        """Test faceted map generation."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_ocean_basin_abundance_faceted(self, mock_savefig):
        """Test faceted abundance plots."""


class TestPhylogeneticTree(unittest.TestCase):
    """Test phylogenetic tree visualization."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_phylogenetic_tree_basic(self, mock_savefig):
        """Test basic tree plotting."""

    def test_plot_phylogenetic_tree_with_bootstrap(self):
        """Test tree with bootstrap values."""


class TestIdentityPlots(unittest.TestCase):
    """Test identity score visualization."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_identity_distribution(self, mock_savefig):
        """Test identity histogram."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_identity_by_status(self, mock_savefig):
        """Test identity box plots by status."""

    @patch('matplotlib.pyplot.savefig')
    def test_plot_assignment_status(self, mock_savefig):
        """Test assignment status charts."""
```

**Deliverables:**
- `test_visualization.py` with 25-30 test cases
- Mock-based tests (no actual image generation)
- Verify correct function calls and arguments
- Edge case handling (no data, single point, etc.)

#### Task 1.4: Create `test_cli.py` 🔴 CRITICAL
**Priority**: Highest
**Effort**: 8-10 hours

**Test Coverage:**
```python
class TestOrganismExtraction(unittest.TestCase):
    """Test organism name extraction from file paths."""

    def test_extract_organism_from_path_standard(self):
        """Test extraction from standard filename."""

    def test_extract_organism_from_path_with_underscores(self):
        """Test extraction with underscores in name."""

    def test_extract_organism_from_path_complex(self):
        """Test extraction from complex path."""


class TestDirectorySetup(unittest.TestCase):
    """Test output directory creation."""

    def test_setup_directories_creates_structure(self):
        """Test that all required directories are created."""

    def test_setup_directories_existing_dirs(self):
        """Test handling of existing directories."""


class TestCLIArgumentParsing(unittest.TestCase):
    """Test command-line argument parsing."""

    def test_cli_minimal_arguments(self):
        """Test CLI with minimal required arguments."""

    def test_cli_full_arguments(self):
        """Test CLI with all optional arguments."""

    def test_cli_config_file(self):
        """Test CLI with config file."""

    def test_cli_invalid_arguments(self):
        """Test CLI error handling for invalid args."""


class TestCLIIntegration(unittest.TestCase):
    """Test full CLI execution."""

    def test_cli_pipeline_minimal(self):
        """Test complete pipeline through CLI."""

    def test_cli_pipeline_with_phylogeny(self):
        """Test CLI with phylogeny enabled."""

    def test_cli_pipeline_error_handling(self):
        """Test CLI error reporting."""
```

**Deliverables:**
- `test_cli.py` with 20-25 test cases
- Argument parsing tests
- Integration tests with real data
- Error handling tests

---

### Phase 2: Important Gaps (Week 3) - IMPORTANT FOR QUALITY

#### Task 2.1: Create `test_config.py` 🟡
**Priority**: High
**Effort**: 6-8 hours

**Test Coverage:**
```python
class TestConfigLoading(unittest.TestCase):
    """Test configuration loading from files."""

    def test_load_config_from_yaml(self):
        """Test loading YAML config file."""

    def test_load_config_from_json(self):
        """Test loading JSON config file."""

    def test_load_config_invalid_format(self):
        """Test error handling for invalid config."""


class TestConfigEnvironment(unittest.TestCase):
    """Test environment variable override."""

    def test_load_config_from_env_basic(self):
        """Test environment variable parsing."""

    def test_load_config_from_env_nested(self):
        """Test nested parameter updates via env vars."""


class TestConfigValidation(unittest.TestCase):
    """Test configuration validation."""

    def test_validate_config_valid(self):
        """Test validation with valid config."""

    def test_validate_config_invalid_threshold(self):
        """Test validation with invalid threshold."""

    def test_validate_config_missing_files(self):
        """Test validation with missing file paths."""


class TestConfigTemplate(unittest.TestCase):
    """Test config template generation."""

    def test_create_config_template_yaml(self):
        """Test YAML template generation."""

    def test_create_config_template_json(self):
        """Test JSON template generation."""
```

**Deliverables:**
- `test_config.py` with 15-20 test cases
- YAML/JSON roundtrip tests
- Environment variable parsing tests
- Validation logic tests

#### Task 2.2: Create `test_utils.py` 🟡
**Priority**: High
**Effort**: 8-10 hours

**Test Coverage:**
```python
class TestLogging(unittest.TestCase):
    """Test logging setup and utilities."""

    def test_setup_logging_basic(self):
        """Test basic logging setup."""

    def test_log_function_call_decorator(self):
        """Test function call logging decorator."""


class TestExternalToolChecking(unittest.TestCase):
    """Test external tool detection and version checking."""

    def test_check_external_tool_found(self):
        """Test tool detection when available."""

    def test_check_external_tool_not_found(self):
        """Test tool detection when missing."""

    def test_get_tool_version(self):
        """Test version extraction."""

    def test_compare_versions(self):
        """Test semantic version comparison."""


class TestProgressTracker(unittest.TestCase):
    """Test progress tracking."""

    def test_progress_tracker_initialization(self):
        """Test tracker initialization."""

    def test_progress_tracker_update(self):
        """Test progress updates and ETA calculation."""


class TestTaxonomyAssignment(unittest.TestCase):
    """Test taxonomy consensus functions."""

    def test_assign_consensus_taxonomy_unanimous(self):
        """Test consensus with unanimous agreement."""

    def test_assign_consensus_taxonomy_majority(self):
        """Test consensus with majority vote."""

    def test_pick_final_group_taxon(self):
        """Test final taxonomy selection."""


class TestSequenceValidation(unittest.TestCase):
    """Test sequence validation functions."""

    def test_validate_sequence_valid(self):
        """Test validation of valid sequence."""

    def test_validate_sequence_too_short(self):
        """Test validation of short sequence."""

    def test_validate_sequence_high_n_content(self):
        """Test validation with high N content."""


class TestFileIO(unittest.TestCase):
    """Test file I/O helper functions."""

    def test_read_fasta(self):
        """Test FASTA reading."""

    def test_write_fasta(self):
        """Test FASTA writing."""

    def test_validate_fasta_file(self):
        """Test FASTA validation."""
```

**Deliverables:**
- `test_utils.py` with 25-30 test cases
- Comprehensive coverage of utility functions
- Mock-based external tool tests

#### Task 2.3: Expand `test_phylogenetics.py` 🟡
**Priority**: High
**Effort**: 6-8 hours

**Test Coverage:**
```python
class TestMAFFTAlignment(unittest.TestCase):
    """Test MAFFT alignment for phylogenetics."""

    @patch('subprocess.run')
    def test_run_mafft_alignment_basic(self, mock_run):
        """Test basic MAFFT alignment."""

    def test_run_mafft_alignment_error(self):
        """Test MAFFT error handling."""


class TestFastTree(unittest.TestCase):
    """Test FastTree wrapper."""

    @patch('subprocess.run')
    def test_run_fasttree_basic(self, mock_run):
        """Test basic FastTree execution."""

    @patch('subprocess.run')
    def test_run_fasttree_with_bootstrap(self, mock_run):
        """Test FastTree with bootstrap."""


class TestPhyMLStub(unittest.TestCase):
    """Test PhyML wrapper (stub function)."""

    def test_run_phyml_not_implemented(self):
        """Test that PhyML raises NotImplementedError."""


class TestTreeRooting(unittest.TestCase):
    """Test tree rooting methods."""

    def test_midpoint_root_tree(self):
        """Test midpoint rooting."""

    def test_outgroup_rooting(self):
        """Test outgroup addition and rooting."""


class TestTreeRelabeling(unittest.TestCase):
    """Test tree tip relabeling."""

    def test_relabel_tree_and_alignment(self):
        """Test relabeling with consensus_group_sp."""


class TestPhylogenyIntegration(unittest.TestCase):
    """Test complete phylogeny building workflow."""

    def test_build_phylogeny_no_tree(self):
        """Test with phylogeny disabled."""

    def test_build_phylogeny_midpoint(self):
        """Test with midpoint rooting."""

    def test_build_phylogeny_outgroup(self):
        """Test with outgroup rooting."""
```

**Deliverables:**
- Expanded phylogenetics tests with 20-25 test cases
- Mock-based FastTree/PhyML tests
- Tree manipulation tests

---

### Phase 3: Supplementary Tests (Week 4) - NICE TO HAVE

#### Task 3.1: Create `test_goas_downloader.py` 🟢
**Priority**: Medium
**Effort**: 4-6 hours

**Test Coverage:**
```python
class TestChecksumCalculation(unittest.TestCase):
    """Test MD5 checksum calculation."""

    def test_calculate_md5_known_file(self):
        """Test checksum with known file."""


class TestFileDownload(unittest.TestCase):
    """Test file download functionality."""

    @patch('urllib.request.urlopen')
    def test_download_file_success(self, mock_urlopen):
        """Test successful file download."""

    def test_download_file_network_error(self):
        """Test download error handling."""


class TestZipExtraction(unittest.TestCase):
    """Test ZIP file extraction."""

    def test_extract_zip_valid(self):
        """Test extraction of valid ZIP."""

    def test_extract_zip_invalid(self):
        """Test error handling for invalid ZIP."""


class TestGoaSSetup(unittest.TestCase):
    """Test complete GOaS setup workflow."""

    @patch('boldgenotyper.goas_downloader.download_file')
    @patch('boldgenotyper.goas_downloader.extract_zip')
    def test_setup_goas_complete(self, mock_extract, mock_download):
        """Test complete setup workflow."""
```

**Deliverables:**
- `test_goas_downloader.py` with 10-15 test cases
- Mock-based download tests
- ZIP extraction tests

#### Task 3.2: Create `test_integration.py` 🟢
**Priority**: Medium
**Effort**: 6-8 hours

**Test Coverage:**
```python
class TestEndToEndPipeline(unittest.TestCase):
    """Test complete pipeline from input to output."""

    def test_pipeline_small_dataset(self):
        """Test pipeline with small test dataset."""

    def test_pipeline_with_all_features(self):
        """Test pipeline with all features enabled."""

    def test_pipeline_output_files_exist(self):
        """Test that all expected outputs are created."""

    def test_pipeline_output_files_valid(self):
        """Test that output files are valid."""
```

**Deliverables:**
- `test_integration.py` with 8-10 end-to-end tests
- Real test dataset (small subset)
- Output validation

#### Task 3.3: Create `conftest.py` with Shared Fixtures 🟢
**Priority**: Medium
**Effort**: 4-6 hours

**Content:**
```python
"""Shared pytest fixtures for BOLDGenotyper tests."""

import pytest
import tempfile
from pathlib import Path

@pytest.fixture
def temp_dir():
    """Create temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)

@pytest.fixture
def sample_fasta():
    """Sample FASTA file for testing."""
    # Create small test FASTA
    pass

@pytest.fixture
def sample_tsv():
    """Sample BOLD TSV file for testing."""
    # Create small test TSV
    pass

@pytest.fixture
def sample_config():
    """Sample configuration for testing."""
    from boldgenotyper.config import get_default_config
    return get_default_config()
```

**Deliverables:**
- `conftest.py` with reusable fixtures
- Reduced code duplication
- Cleaner test code

---

## Testing Standards & Best Practices

### Code Style

1. **Follow PEP 8**: Use `black` for formatting
2. **Type Hints**: Include type hints in test functions where helpful
3. **Docstrings**: Every test class and function must have descriptive docstring
4. **Naming**: Use descriptive test names following pattern:
   ```python
   def test_<function>_<scenario>_<expected_result>(self):
   ```

### Test Structure

```python
def test_example_function_with_valid_input(self):
    """Test example function with valid input returns expected result."""
    # Arrange
    input_data = create_test_data()
    expected_output = "expected_value"

    # Act
    actual_output = example_function(input_data)

    # Assert
    self.assertEqual(actual_output, expected_output)
```

### Mocking Best Practices

1. **Mock External Dependencies**: Always mock subprocess calls, network requests, file I/O
2. **Use `@patch` Decorators**: Apply patches at the right level
3. **Verify Mock Calls**: Use `assert_called_once_with()` to verify behavior
4. **Don't Over-Mock**: Test real logic, only mock external boundaries

Example:
```python
@patch('subprocess.run')
def test_run_mafft_calls_correct_command(self, mock_run):
    """Test that MAFFT is called with correct arguments."""
    mock_run.return_value = Mock(returncode=0, stdout="", stderr="")

    run_mafft_alignment("input.fasta", "output.fasta")

    mock_run.assert_called_once()
    call_args = mock_run.call_args[0][0]
    self.assertIn('mafft', call_args)
```

### Test Data Management

1. **Small Test Data**: Use minimal datasets (5-10 sequences)
2. **Synthetic Data**: Create synthetic data with known properties
3. **Version Control**: Include test data in `tests/test_data/`
4. **Data Generators**: Create helper functions to generate test data

### Coverage Goals

- **Overall Coverage**: 85-90%
- **Critical Modules**: 90%+ (core, dereplication, assignment)
- **Utility Modules**: 75%+ (utils, config)
- **Visualization**: 70%+ (due to mocking complexity)

### Continuous Improvement

1. **Run Coverage Reports**: Use `pytest-cov` to generate reports
2. **Review Untested Lines**: Identify and test uncovered code
3. **Update Tests**: Keep tests in sync with code changes
4. **Refactor Tests**: Improve test quality over time

---

## Documentation Requirements

### Test File Documentation

Every test file must include:

```python
"""
Unit tests for [module name].

Tests cover:
- [Feature 1]
- [Feature 2]
- [Feature 3]
- Edge cases: [describe]
- Error handling: [describe]

Test Categories:
- Unit tests: Test individual functions
- Integration tests: Test module interactions
- Mock tests: Test external tool integration

Dependencies:
- External tools: [list if applicable]
- Test data: [describe test data requirements]
"""
```

### Test Function Documentation

Every test function must have a docstring:

```python
def test_parse_coordinates_valid_decimal_format(self):
    """
    Test parsing coordinates in decimal format (lat,lon).

    Input: "40.7128,-74.0060"
    Expected: lat=40.7128, lon=-74.0060
    """
```

### README for Tests

Create `tests/README.md`:

```markdown
# BOLDGenotyper Test Suite

## Overview

Comprehensive test suite for BOLDGenotyper package covering:
- Unit tests for all modules
- Integration tests for workflows
- End-to-end pipeline tests
- Mock-based tests for external tools

## Running Tests

### All Tests
```bash
pytest
```

### With Coverage
```bash
pytest --cov=boldgenotyper --cov-report=html
```

### Specific Module
```bash
pytest tests/test_dereplication.py
```

### Verbose Output
```bash
pytest -v
```

## Test Organization

- `test_*.py`: Test files for each module
- `test_data/`: Test data files
- `conftest.py`: Shared fixtures
- `test_integration.py`: End-to-end tests

## Writing New Tests

See TESTING_PLAN.md for guidelines and standards.
```

---

## Continuous Integration Setup

### GitHub Actions Workflow

Create `.github/workflows/tests.yml`:

```yaml
name: Tests

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main, develop ]

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        python-version: ['3.8', '3.9', '3.10', '3.11']

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}

    - name: Install system dependencies
      run: |
        if [ "$RUNNER_OS" == "Linux" ]; then
          sudo apt-get update
          sudo apt-get install -y mafft trimal fasttree
        elif [ "$RUNNER_OS" == "macOS" ]; then
          brew install mafft trimal fasttree
        fi

    - name: Install Python dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -e ".[dev,all]"

    - name: Run tests with coverage
      run: |
        pytest --cov=boldgenotyper --cov-report=xml --cov-report=term

    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml
        flags: unittests
        name: codecov-${{ matrix.os }}-py${{ matrix.python-version }}

  lint:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'

    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install flake8 black mypy

    - name: Lint with flake8
      run: |
        flake8 boldgenotyper --count --select=E9,F63,F7,F82 --show-source --statistics
        flake8 boldgenotyper --count --exit-zero --max-complexity=10 --max-line-length=100

    - name: Check formatting with black
      run: |
        black --check boldgenotyper

    - name: Type check with mypy
      run: |
        mypy boldgenotyper --ignore-missing-imports
```

### Codecov Integration

1. Sign up at https://codecov.io/
2. Add BOLDGenotyper repository
3. Add badge to README:
   ```markdown
   [![codecov](https://codecov.io/gh/SymbioSeas/BOLDGenotyper/branch/main/graph/badge.svg)](https://codecov.io/gh/SymbioSeas/BOLDGenotyper)
   ```

### Pre-commit Hooks (Optional)

Create `.pre-commit-config.yaml`:

```yaml
repos:
  - repo: https://github.com/psf/black
    rev: 23.3.0
    hooks:
      - id: black
        language_version: python3.10

  - repo: https://github.com/pycqa/flake8
    rev: 6.0.0
    hooks:
      - id: flake8
        args: [--max-line-length=100]

  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.3.0
    hooks:
      - id: mypy
        args: [--ignore-missing-imports]
```

---

## Timeline and Resources

### Phase 1: Critical Tests (Weeks 1-2)
**Effort**: 32-42 hours
**Priority**: REQUIRED FOR JOSS

| Task | Effort | Priority |
|------|--------|----------|
| Implement `core.py` and tests | 12h | 🔴 Critical |
| Create `test_reports.py` | 8h | 🔴 Critical |
| Create `test_visualization.py` | 12h | 🔴 Critical |
| Create `test_cli.py` | 10h | 🔴 Critical |
| **Total Phase 1** | **42h** | |

### Phase 2: Important Tests (Week 3)
**Effort**: 20-26 hours
**Priority**: IMPORTANT FOR QUALITY

| Task | Effort | Priority |
|------|--------|----------|
| Create `test_config.py` | 8h | 🟡 High |
| Create `test_utils.py` | 10h | 🟡 High |
| Expand `test_phylogenetics.py` | 8h | 🟡 High |
| **Total Phase 2** | **26h** | |

### Phase 3: Supplementary Tests (Week 4)
**Effort**: 14-20 hours
**Priority**: NICE TO HAVE

| Task | Effort | Priority |
|------|--------|----------|
| Create `test_goas_downloader.py` | 6h | 🟢 Medium |
| Create `test_integration.py` | 8h | 🟢 Medium |
| Create `conftest.py` | 6h | 🟢 Medium |
| **Total Phase 3** | **20h** | |

### Phase 4: Documentation & CI (Week 5)
**Effort**: 8-12 hours
**Priority**: REQUIRED FOR JOSS

| Task | Effort | Priority |
|------|--------|----------|
| Create `tests/README.md` | 2h | 🔴 Critical |
| Set up GitHub Actions | 4h | 🔴 Critical |
| Set up Codecov | 2h | 🔴 Critical |
| Generate coverage report | 2h | 🔴 Critical |
| Update documentation | 2h | 🔴 Critical |
| **Total Phase 4** | **12h** | |

### Total Effort Summary

| Phase | Effort | Deliverables |
|-------|--------|--------------|
| Phase 1 | 42h | Core, reports, viz, CLI tests |
| Phase 2 | 26h | Config, utils, phylo tests |
| Phase 3 | 20h | GOaS, integration, fixtures |
| Phase 4 | 12h | Docs, CI/CD, coverage |
| **TOTAL** | **100h** | **~200-250 new test cases** |

### Recommended Schedule

- **Week 1**: Phase 1 (Part 1) - Core + Reports
- **Week 2**: Phase 1 (Part 2) - Visualization + CLI
- **Week 3**: Phase 2 - Config + Utils + Phylo
- **Week 4**: Phase 3 - GOaS + Integration
- **Week 5**: Phase 4 - Documentation + CI/CD

---

## Success Criteria for Publication

### JOSS Requirements ✅

1. **Automated Testing** ✅
   - [ ] Test suite covers all major functions
   - [ ] Tests run automatically on CI/CD
   - [ ] Coverage report available

2. **Documentation** ✅
   - [ ] Test documentation in `tests/README.md`
   - [ ] Test docstrings for all test functions
   - [ ] Coverage badge in main README

3. **Quality Metrics** ✅
   - [ ] 85-90% code coverage achieved
   - [ ] All tests pass on multiple Python versions
   - [ ] Tests pass on Linux and macOS

4. **Maintainability** ✅
   - [ ] Clear test organization
   - [ ] Reusable fixtures and helpers
   - [ ] Mock-based tests for external dependencies

### PLOS ONE Case Study Requirements ✅

1. **Reproducibility** ✅
   - [ ] Integration tests demonstrate full workflow
   - [ ] Test data included or easily generated
   - [ ] Results can be reproduced from tests

2. **Validation** ✅
   - [ ] Tests validate scientific accuracy (distances, clustering)
   - [ ] Tests verify output formats
   - [ ] Tests check edge cases

### Publication-Ready Checklist

- [ ] **Core Pipeline**: `test_core.py` implemented and passing
- [ ] **Report Generation**: `test_reports.py` implemented and passing
- [ ] **Visualization**: `test_visualization.py` implemented and passing
- [ ] **CLI**: `test_cli.py` implemented and passing
- [ ] **Configuration**: `test_config.py` implemented and passing
- [ ] **Utilities**: `test_utils.py` implemented and passing
- [ ] **Phylogenetics**: Expanded tests implemented and passing
- [ ] **Integration**: End-to-end tests implemented and passing
- [ ] **Coverage**: 85-90% code coverage achieved
- [ ] **CI/CD**: GitHub Actions workflow set up and passing
- [ ] **Codecov**: Coverage reporting integrated
- [ ] **Documentation**: Test README and docstrings complete
- [ ] **Badge**: Coverage badge added to main README
- [ ] **Multi-platform**: Tests pass on Linux, macOS, Windows WSL
- [ ] **Multi-version**: Tests pass on Python 3.8, 3.9, 3.10, 3.11

---

## Appendix: Module Function Inventory

### Complete Function List by Module

#### `core.py` (82 lines)
- [ ] `run_pipeline()` - Main orchestration (STUB - needs implementation)

#### `dereplication.py` (885 lines)
- [x] `check_external_tools()` ✅
- [x] `dereplicate_from_fasta()` ✅
- [x] `run_mafft_alignment()` ✅
- [x] `run_trimal_trimming()` ✅
- [x] `calculate_pairwise_distances()` ✅
- [x] `_compute_distance()` ✅
- [x] `cluster_sequences()` ✅
- [x] `generate_consensus()` ✅
- [x] `dereplicate_sequences()` ✅

#### `metadata.py` (1,138 lines)
- [x] `parse_bold_tsv()` ✅
- [x] `validate_required_columns()` ✅
- [x] `_log_data_quality()` ✅
- [x] `extract_coordinates()` ✅
- [x] `infer_coord_order()` ✅
- [x] `parse_coordinates_column()` ✅
- [x] `_validate_coordinate_ranges()` ✅
- [x] `filter_by_coordinate_quality()` ✅
- [x] `get_coordinate_quality_stats()` ✅
- [x] `levenshtein_distance()` ✅
- [x] `calculate_edit_distance()` ✅
- [x] `calculate_identity()` ✅
- [x] `assign_genotypes()` ✅ (legacy location)
- [x] `_load_consensus_sequences()` ✅
- [x] `_load_raw_sequences_from_fasta()` ✅
- [x] `_get_sequences_from_tsv()` ✅
- [x] `_assign_single_sample()` ✅
- [x] `_log_assignment_stats()` ✅
- [x] `_save_diagnostics()` ✅
- [x] `validate_metadata_for_analysis()` ✅

#### `genotype_assignment.py` (1,068 lines)
- [x] `check_edlib_available()` ✅
- [x] `levenshtein_distance()` ✅
- [x] `calculate_edit_distance()` ✅
- [x] `calculate_identity()` ✅
- [x] `parse_cigar()` ✅
- [x] `calculate_identity_with_cigar()` ✅
- [x] `extract_processid_from_header()` ✅
- [x] `read_fasta_simple()` ✅
- [x] `find_best_consensus_match()` ✅
- [x] `_assignment_worker()` ✅
- [ ] `assign_species_to_sample()` (STUB)
- [ ] `assign_species_to_consensus()` (STUB)
- [x] `assign_genotypes()` ✅

#### `geographic.py` (1,156 lines)
- [x] `check_geopandas_available()` ✅
- [x] `check_cartopy_available()` ✅
- [x] `load_goas_data()` ✅
- [x] `create_points_geodataframe()` ✅
- [x] `assign_ocean_basins()` ✅
- [x] `get_basin_counts()` ✅
- [x] `get_basin_summary_stats()` ✅
- [x] `validate_marine_coordinates()` ✅
- [x] `filter_marine_samples()` ✅
- [x] `load_land_geometries()` ✅
- [x] `classify_and_snap_coordinates()` ✅

#### `phylogenetics.py` (463 lines)
- [ ] `build_phylogeny()` ⚠️ (partial)
- [ ] `run_mafft_alignment()` ⚠️ (partial)
- [ ] `run_fasttree()` ⚠️ (partial)
- [ ] `run_phyml()` (STUB)
- [ ] `midpoint_root_tree()` (STUB)
- [ ] `add_outgroup()` ⚠️
- [ ] `relabel_tree_and_alignment()` ⚠️

#### `visualization.py` (982 lines)
- [ ] `get_genotype_colors()` ❌
- [ ] `plot_distribution_map()` ❌
- [ ] `plot_ocean_basin_abundance()` ❌
- [ ] `plot_distribution_map_faceted()` ❌
- [ ] `plot_ocean_basin_abundance_faceted()` ❌
- [ ] `plot_phylogenetic_tree()` ❌
- [ ] `plot_identity_distribution()` ❌
- [ ] `plot_identity_by_status()` ❌
- [ ] `plot_assignment_status()` ❌

#### `config.py` (1,028 lines)
- [ ] `get_default_config()` ⚠️
- [ ] `load_config_from_file()` ❌
- [ ] `_load_yaml_config()` ❌
- [ ] `_load_json_config()` ❌
- [ ] `_dict_to_config()` ⚠️
- [ ] `_convert_paths_to_strings()` ❌
- [ ] `_convert_strings_to_paths()` ❌
- [ ] `load_config_from_env()` ❌
- [ ] `_parse_env_value()` ❌
- [ ] `validate_config()` ❌
- [ ] `create_config_template()` ❌

#### `utils.py` (1,327 lines)
- [ ] `setup_logging()` ⚠️
- [ ] `log_function_call()` ❌
- [ ] `check_external_tool()` ⚠️
- [ ] `get_tool_version()` ❌
- [ ] `parse_version_string()` ❌
- [ ] `compare_versions()` ❌
- [ ] `get_tool_installation_instructions()` ❌
- [ ] `create_output_directory()` ⚠️
- [ ] `safe_file_path()` ❌
- [ ] `sanitize_filename()` ❌
- [ ] `extract_organism_name()` ❌
- [ ] `read_fasta()` ⚠️
- [ ] `write_fasta()` ⚠️
- [ ] `validate_fasta_file()` ⚠️
- [ ] `validate_tsv_file()` ⚠️
- [ ] `validate_sequence()` ❌
- [ ] `calculate_sequence_identity()` ❌
- [ ] `get_sequence_stats()` ❌
- [ ] `format_elapsed_time()` ❌
- [ ] `format_file_size()` ❌
- [ ] `get_timestamp()` ❌
- [ ] `extract_processid_from_header()` ⚠️
- [ ] `ProgressTracker` class ❌
- [ ] `assign_consensus_taxonomy()` ❌
- [ ] `pick_final_group_taxon()` ❌

#### `reports.py` (415 lines)
- [ ] `generate_taxonomy_conflicts_report()` ❌
- [ ] `generate_assignment_summary()` ❌
- [ ] `generate_consensus_characterization()` ❌
- [ ] `generate_sequence_quality_metrics()` ❌

#### `goas_downloader.py` (291 lines)
- [ ] `calculate_md5()` ❌
- [ ] `download_file()` ❌
- [ ] `extract_zip()` ❌
- [ ] `verify_files()` ❌
- [ ] `create_citation_file()` ❌
- [ ] `setup_goas()` ❌

#### `cli.py` (~400-500 lines estimated)
- [ ] `extract_organism_from_path()` ⚠️
- [ ] `setup_directories()` ⚠️
- [ ] `run_pipeline()` ⚠️
- [ ] `main()` ⚠️

**Legend:**
- ✅ Well-tested (75-95% coverage)
- ⚠️ Partially tested (30-60% coverage)
- ❌ Untested or minimal coverage (0-30%)

---

## References

### Testing Best Practices
- Python Testing with pytest (Brian Okken)
- Effective Python Testing with Pytest (Lorenzo Cellentani)
- JOSS Review Checklist: https://joss.readthedocs.io/en/latest/review_checklist.html

### Tools
- pytest: https://docs.pytest.org/
- pytest-cov: https://pytest-cov.readthedocs.io/
- unittest.mock: https://docs.python.org/3/library/unittest.mock.html
- Codecov: https://codecov.io/

### Related Documentation
- `PUBLICATION_READY_SUMMARY.md` - Documentation status
- `DOCUMENTATION_PLAN.md` - Full documentation plan
- `CHANGELOG.md` - Version history
- `README.md` - Package overview

---

**Document Prepared By**: Claude Code Assistant
**Date**: 2025-11-18
**Status**: READY FOR IMPLEMENTATION
**Next Action**: Begin Phase 1 implementation
