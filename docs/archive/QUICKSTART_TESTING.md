# Quick Start: Running BOLDGenotyper Tests

This guide will help you quickly set up and run the new core pipeline tests.

## Prerequisites

Ensure you have Python 3.8+ installed:

```bash
python3 --version
```

## Installation

### Option 1: Install with All Dependencies (Recommended)

```bash
cd /Users/stesmith/Documents/depredation/boldgenotyper
pip install -e ".[dev,all]"
```

This installs:
- Core dependencies (pandas, numpy, biopython, scipy, matplotlib, seaborn)
- Geographic dependencies (geopandas, cartopy, shapely)
- Phylogenetic dependencies (ete3, dendropy)
- Development dependencies (pytest, pytest-cov, black, flake8, mypy, sphinx)

### Option 2: Install Minimal Testing Dependencies

```bash
cd /Users/stesmith/Documents/depredation/boldgenotyper
pip install -e .
pip install pytest pytest-cov
```

## Running Tests

### Run Core Pipeline Tests

```bash
# Run all core tests with verbose output
python3 -m pytest tests/test_core.py -v

# Or with unittest
python3 -m unittest tests.test_core -v
```

### Run with Coverage Report

```bash
# Generate coverage report
python3 -m pytest tests/test_core.py --cov=boldgenotyper.core --cov-report=term --cov-report=html

# View HTML report
open htmlcov/index.html  # macOS
# or
xdg-open htmlcov/index.html  # Linux
```

### Run Specific Test Class

```bash
# Run only helper function tests
python3 -m pytest tests/test_core.py::TestHelperFunctions -v

# Run only integration tests
python3 -m pytest tests/test_core.py::TestPipelineIntegration -v

# Run only error handling tests
python3 -m pytest tests/test_core.py::TestErrorHandling -v
```

### Run Single Test

```bash
python3 -m pytest tests/test_core.py::TestHelperFunctions::test_extract_organism_name_genus_species -v
```

## Expected Output

### Successful Test Run

```
============================= test session starts ==============================
platform darwin -- Python 3.X.X, pytest-X.X.X, pluggy-X.X.X -- /usr/bin/python3
cachedir: .pytest_cache
rootdir: /Users/stesmith/Documents/depredation/boldgenotyper
collected 23 items

tests/test_core.py::TestHelperFunctions::test_extract_organism_name_genus_species PASSED [  4%]
tests/test_core.py::TestHelperFunctions::test_extract_organism_name_single_genus PASSED [  8%]
tests/test_core.py::TestHelperFunctions::test_extract_organism_name_with_path PASSED [ 13%]
tests/test_core.py::TestHelperFunctions::test_extract_organism_name_ignores_data_suffix PASSED [ 17%]
tests/test_core.py::TestHelperFunctions::test_extract_organism_name_ignores_bold_suffix PASSED [ 21%]
tests/test_core.py::TestHelperFunctions::test_extract_organism_name_capitalizes PASSED [ 26%]
tests/test_core.py::TestHelperFunctions::test_setup_directories_creates_structure PASSED [ 30%]
tests/test_core.py::TestHelperFunctions::test_setup_directories_existing_dirs PASSED [ 34%]
tests/test_core.py::TestInputValidation::test_run_pipeline_missing_input_file PASSED [ 39%]
tests/test_core.py::TestInputValidation::test_run_pipeline_invalid_cluster_threshold_too_low PASSED [ 43%]
tests/test_core.py::TestInputValidation::test_run_pipeline_invalid_cluster_threshold_too_high PASSED [ 47%]
tests/test_core.py::TestInputValidation::test_run_pipeline_invalid_min_identity_too_low PASSED [ 52%]
tests/test_core.py::TestInputValidation::test_run_pipeline_invalid_threads PASSED [ 56%]
tests/test_core.py::TestPipelineIntegration::test_run_pipeline_minimal_success PASSED [ 60%]
tests/test_core.py::TestPipelineIntegration::test_run_pipeline_skip_geographic PASSED [ 65%]
tests/test_core.py::TestPipelineIntegration::test_run_pipeline_missing_goas_shapefile PASSED [ 69%]
tests/test_core.py::TestPipelineIntegration::test_run_pipeline_with_phylogeny PASSED [ 73%]
tests/test_core.py::TestPipelineIntegration::test_run_pipeline_custom_organism_name PASSED [ 78%]
tests/test_core.py::TestPipelineIntegration::test_run_pipeline_custom_config PASSED [ 82%]
tests/test_core.py::TestErrorHandling::test_run_pipeline_phase1_failure PASSED [ 86%]
tests/test_core.py::TestErrorHandling::test_run_pipeline_visualization_failure_non_critical PASSED [ 91%]
tests/test_core.py::TestErrorHandling::test_run_pipeline_report_generation_failure_non_critical PASSED [ 95%]
tests/test_core.py::TestResultsDictionary::test_results_dict_structure PASSED [100%]

============================== 23 passed in 2.50s ===============================
```

### Coverage Report

```
---------- coverage: platform darwin, python 3.X.X -----------
Name                      Stmts   Miss  Cover
---------------------------------------------
boldgenotyper/core.py      450     35    92%
---------------------------------------------
TOTAL                      450     35    92%

Coverage HTML written to htmlcov/index.html
```

## Troubleshooting

### Issue: ModuleNotFoundError

**Error**: `ModuleNotFoundError: No module named 'pandas'` (or other module)

**Fix**:
```bash
pip install pandas biopython numpy scipy matplotlib seaborn
```

### Issue: pytest not found

**Error**: `python3: No module named pytest`

**Fix**:
```bash
pip install pytest pytest-cov
```

### Issue: Import errors from boldgenotyper

**Error**: `ImportError: cannot import name 'core' from 'boldgenotyper'`

**Fix**: Install the package in development mode:
```bash
cd /Users/stesmith/Documents/depredation/boldgenotyper
pip install -e .
```

### Issue: Tests create many temporary files

This is normal. Tests use `tempfile.TemporaryDirectory()` which auto-cleans after each test. No manual cleanup needed.

## Running All Tests

To run the entire test suite (not just core tests):

```bash
# Run all tests
python3 -m pytest tests/ -v

# With coverage for entire package
python3 -m pytest tests/ --cov=boldgenotyper --cov-report=html

# Exclude slow integration tests
python3 -m pytest tests/ -v -m "not slow"
```

## Continuous Integration

The tests are designed to run in CI/CD environments. Example GitHub Actions workflow:

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - name: Install dependencies
        run: |
          pip install -e ".[dev]"
      - name: Run tests
        run: |
          pytest tests/ --cov=boldgenotyper --cov-report=xml
      - name: Upload coverage
        uses: codecov/codecov-action@v3
```

## Verification Checklist

After running tests, verify:

- [ ] All 23 core tests pass
- [ ] Coverage is ≥ 90% for `boldgenotyper/core.py`
- [ ] No unexpected warnings
- [ ] HTML coverage report generated
- [ ] Tests run in < 5 seconds

## Next Steps

1. **Review Coverage Report**: Open `htmlcov/index.html` and review untested lines
2. **Run Full Test Suite**: Run all tests to ensure no regressions
3. **Set Up CI/CD**: Configure GitHub Actions for automated testing
4. **Proceed to Phase 2**: Implement tests for reports, visualization, and CLI

## Documentation

- **Test Documentation**: See `tests/README_CORE_TESTS.md`
- **Implementation Summary**: See `PHASE1_IMPLEMENTATION_SUMMARY.md`
- **Testing Plan**: See `TESTING_PLAN.md`

## Contact

For issues or questions:
- **Repository**: https://github.com/SymbioSeas/BOLDGenotyper
- **Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues

---

**Quick Commands**:
```bash
# Install and test (copy-paste ready)
cd /Users/stesmith/Documents/depredation/boldgenotyper
pip install -e ".[dev]"
python3 -m pytest tests/test_core.py -v --cov=boldgenotyper.core --cov-report=html
open htmlcov/index.html
```

**Status**: ✅ Phase 1 Complete - Ready to Run
