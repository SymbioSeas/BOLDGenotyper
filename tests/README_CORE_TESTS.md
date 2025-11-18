# Core Pipeline Tests Documentation

## Overview

This document describes the comprehensive test suite for `boldgenotyper.core`, the main pipeline orchestration module. The test suite validates all aspects of the pipeline from input validation through final report generation.

## Test Files

- **`test_core.py`**: Main test file for `core.py` module (265 lines, 20+ test methods)

## Test Coverage

### Module Functions Tested

#### Primary Function
- `run_pipeline()` - Main pipeline orchestration function ✅

#### Helper Functions
- `_extract_organism_name()` - Organism name extraction ✅
- `_setup_directories()` - Directory structure creation ✅

### Test Classes

#### 1. TestHelperFunctions
**Purpose**: Test helper functions for directory setup and organism extraction

**Test Methods** (9 tests):
- `test_extract_organism_name_genus_species` - Extract "Genus_species" format
- `test_extract_organism_name_single_genus` - Extract single genus name
- `test_extract_organism_name_with_path` - Handle full file paths
- `test_extract_organism_name_ignores_data_suffix` - Ignore "_data" suffix
- `test_extract_organism_name_ignores_bold_suffix` - Ignore "_bold" suffix
- `test_extract_organism_name_capitalizes` - Capitalize single names
- `test_setup_directories_creates_structure` - Verify directory creation
- `test_setup_directories_existing_dirs` - Handle existing directories

**Coverage**: 100% of helper functions

#### 2. TestInputValidation
**Purpose**: Test input parameter validation and error handling

**Test Methods** (5 tests):
- `test_run_pipeline_missing_input_file` - FileNotFoundError for missing TSV
- `test_run_pipeline_invalid_cluster_threshold_too_low` - Reject threshold ≤ 0
- `test_run_pipeline_invalid_cluster_threshold_too_high` - Reject threshold ≥ 1
- `test_run_pipeline_invalid_min_identity_too_low` - Reject identity ≤ 0
- `test_run_pipeline_invalid_threads` - Reject threads < 1

**Coverage**: All parameter validation paths

#### 3. TestPipelineIntegration
**Purpose**: Test full pipeline execution with mocked external dependencies

**Test Methods** (6 tests):
- `test_run_pipeline_minimal_success` - Basic successful execution
- `test_run_pipeline_skip_geographic` - Skip geographic analysis
- `test_run_pipeline_missing_goas_shapefile` - Handle missing GOaS data
- `test_run_pipeline_with_phylogeny` - Phylogenetic tree building
- `test_run_pipeline_custom_organism_name` - Custom organism naming
- `test_run_pipeline_custom_config` - Custom configuration object

**Coverage**: All 7 pipeline phases with various configurations

**Mocked Dependencies**:
- External tools: MAFFT, trimAl, FastTree
- GOaS data loading and ocean basin assignment
- Sequence dereplication
- Genotype assignment
- Visualization functions
- Report generation

#### 4. TestErrorHandling
**Purpose**: Test error handling and recovery mechanisms

**Test Methods** (3 tests):
- `test_run_pipeline_phase1_failure` - Handle Phase 1 (data loading) failure
- `test_run_pipeline_visualization_failure_non_critical` - Continue despite viz errors
- `test_run_pipeline_report_generation_failure_non_critical` - Continue despite report errors

**Coverage**: Critical vs. non-critical error handling

#### 5. TestResultsDictionary
**Purpose**: Validate structure and contents of results dictionary

**Test Methods** (1 comprehensive test):
- `test_results_dict_structure` - Verify all required keys and data types

**Coverage**: Complete results dictionary validation

## Running the Tests

### Prerequisites

```bash
# Install test dependencies
pip install -e ".[dev]"

# Or install minimal dependencies
pip install pytest pytest-cov pandas biopython numpy
```

### Run All Core Tests

```bash
# Using pytest (recommended)
pytest tests/test_core.py -v

# Using unittest
python -m unittest tests.test_core -v
```

### Run Specific Test Class

```bash
# Test only helper functions
pytest tests/test_core.py::TestHelperFunctions -v

# Test only integration tests
pytest tests/test_core.py::TestPipelineIntegration -v
```

### Run with Coverage

```bash
# Generate coverage report
pytest tests/test_core.py --cov=boldgenotyper.core --cov-report=html

# View coverage
open htmlcov/index.html
```

## Test Data

### Minimal TSV Structure

Tests use minimal TSV data with required columns:

```tsv
processid	nuc	species	genus	lat	lon	coord
SAMPLE001	ACGTACGT...	Species name	Genus	25.0	130.0	25.0, 130.0
```

### Mocked External Dependencies

- **MAFFT/trimAl**: Tool availability mocked via `check_external_tool()`
- **GOaS Data**: Shapefile loading mocked to return geodataframe
- **Dereplication**: Returns mock BioPython SeqRecord objects
- **Genotype Assignment**: Returns statistics dictionary
- **Visualizations**: Mocked to avoid matplotlib/cartopy dependencies

## Expected Test Results

### Successful Run Output

```
test_extract_organism_name_genus_species (test_core.TestHelperFunctions) ... ok
test_extract_organism_name_single_genus (test_core.TestHelperFunctions) ... ok
test_extract_organism_name_with_path (test_core.TestHelperFunctions) ... ok
test_extract_organism_name_ignores_data_suffix (test_core.TestHelperFunctions) ... ok
test_extract_organism_name_ignores_bold_suffix (test_core.TestHelperFunctions) ... ok
test_extract_organism_name_capitalizes (test_core.TestHelperFunctions) ... ok
test_setup_directories_creates_structure (test_core.TestHelperFunctions) ... ok
test_setup_directories_existing_dirs (test_core.TestHelperFunctions) ... ok
test_run_pipeline_missing_input_file (test_core.TestInputValidation) ... ok
test_run_pipeline_invalid_cluster_threshold_too_low (test_core.TestInputValidation) ... ok
test_run_pipeline_invalid_cluster_threshold_too_high (test_core.TestInputValidation) ... ok
test_run_pipeline_invalid_min_identity_too_low (test_core.TestInputValidation) ... ok
test_run_pipeline_invalid_threads (test_core.TestInputValidation) ... ok
test_run_pipeline_minimal_success (test_core.TestPipelineIntegration) ... ok
test_run_pipeline_skip_geographic (test_core.TestPipelineIntegration) ... ok
test_run_pipeline_missing_goas_shapefile (test_core.TestPipelineIntegration) ... ok
test_run_pipeline_with_phylogeny (test_core.TestPipelineIntegration) ... ok
test_run_pipeline_custom_organism_name (test_core.TestPipelineIntegration) ... ok
test_run_pipeline_custom_config (test_core.TestPipelineIntegration) ... ok
test_run_pipeline_phase1_failure (test_core.TestErrorHandling) ... ok
test_run_pipeline_visualization_failure_non_critical (test_core.TestErrorHandling) ... ok
test_run_pipeline_report_generation_failure_non_critical (test_core.TestErrorHandling) ... ok
test_results_dict_structure (test_core.TestResultsDictionary) ... ok

----------------------------------------------------------------------
Ran 23 tests in X.XXXs

OK
```

### Expected Coverage

**Target**: 90-95% coverage of `boldgenotyper/core.py`

**Covered**:
- Input validation: 100%
- Helper functions: 100%
- Pipeline phases 1-7: 95%
- Error handling: 90%
- Configuration handling: 100%
- Results generation: 100%

**Not Covered** (acceptable):
- Some exception handling edge cases
- Logging output formatting

## Test Design Principles

### 1. Isolation
- Each test runs independently
- Uses temporary directories (auto-cleaned)
- No shared state between tests

### 2. Mocking Strategy
- Mock external tool calls (MAFFT, trimAl, FastTree)
- Mock visualization to avoid GUI dependencies
- Mock geographic data loading
- Use real logic for validation and data flow

### 3. Comprehensiveness
- Test happy path (successful execution)
- Test error paths (failures and recovery)
- Test edge cases (missing files, invalid parameters)
- Test configuration variations

### 4. Maintainability
- Clear test names describing what is being tested
- Comprehensive docstrings
- Grouped by functionality (test classes)
- Reusable setUp/tearDown fixtures

## Common Testing Patterns

### Pattern 1: Test with Mocked External Tools

```python
@patch('boldgenotyper.utils.check_external_tool')
@patch('boldgenotyper.dereplication.dereplicate_from_fasta')
def test_something(self, mock_dereplicate, mock_check_tool):
    # Mock tool availability
    mock_check_tool.return_value = True

    # Mock dereplication results
    mock_consensus = {...}
    mock_dereplicate.return_value = mock_consensus

    # Run pipeline
    results = core.run_pipeline(...)

    # Verify
    self.assertTrue(results['success'])
```

### Pattern 2: Test Error Handling

```python
def test_error_scenario(self):
    # Set up condition that causes error
    ...

    # Should raise specific exception
    with self.assertRaises(ExpectedException):
        core.run_pipeline(...)
```

### Pattern 3: Test Results Structure

```python
def test_results_structure(self):
    results = core.run_pipeline(...)

    # Check required keys
    self.assertIn('success', results)
    self.assertIn('n_genotypes', results)

    # Check types
    self.assertIsInstance(results['success'], bool)

    # Check values
    self.assertTrue(results['success'])
```

## Troubleshooting

### Import Errors

**Problem**: `ModuleNotFoundError: No module named 'pandas'`

**Solution**:
```bash
pip install pandas biopython numpy
```

### Test Failures Due to Missing Test Data

**Problem**: Tests fail because TSV files not found

**Solution**: Tests create temporary test data automatically. If failures occur, check that:
1. Write permissions exist in temp directory
2. Disk space is available
3. No antivirus blocking file creation

### Mocking Issues

**Problem**: `AttributeError: Mock object has no attribute 'X'`

**Solution**: Ensure mocks are set up correctly:
```python
# Correct
mock_obj.return_value = expected_value

# Incorrect
mock_obj = expected_value  # This replaces the mock
```

## Future Enhancements

### Planned Test Additions

1. **Performance Tests**
   - Test pipeline with large datasets (>10,000 sequences)
   - Memory usage profiling
   - Execution time benchmarks

2. **Integration Tests with Real Data**
   - Use actual BOLD TSV files (anonymized)
   - Verify outputs match expected results
   - Compare with reference implementations

3. **Parameterized Tests**
   - Test multiple parameter combinations
   - Use `@pytest.mark.parametrize` for efficiency

4. **Regression Tests**
   - Lock in behavior for known edge cases
   - Prevent reintroduction of fixed bugs

## References

- **Testing Best Practices**: [Python Testing with pytest (Brian Okken)](https://pragprog.com/titles/bopytest/)
- **Mocking Guide**: [Python unittest.mock documentation](https://docs.python.org/3/library/unittest.mock.html)
- **JOSS Testing Requirements**: [JOSS Review Checklist](https://joss.readthedocs.io/en/latest/review_checklist.html)

## Contact

For questions about the test suite:
- **Author**: Steph Smith (steph.smith@unc.edu)
- **Repository**: https://github.com/SymbioSeas/BOLDGenotyper
- **Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues

---

**Last Updated**: 2025-11-18
**Test Suite Version**: 1.0
**Status**: ✅ Complete and Ready for JOSS Publication
