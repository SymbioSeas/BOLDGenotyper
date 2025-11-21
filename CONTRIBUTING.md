# Contributing to BOLDGenotyper

Thank you for your interest in contributing to BOLDGenotyper! We welcome contributions from the community and appreciate your help in making this software better.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)
- [Pull Request Process](#pull-request-process)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing Guidelines](#testing-guidelines)
- [Documentation](#documentation)

## Code of Conduct

This project adheres to the Contributor Covenant [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to steph.smith@unc.edu.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/BOLDGenotyper.git
   cd BOLDGenotyper
   ```
3. **Set up the development environment** (see [Development Setup](#development-setup))
4. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## How to Contribute

There are many ways to contribute to BOLDGenotyper:

- **Report bugs** and request features via [GitHub Issues](https://github.com/SymbioSeas/BOLDGenotyper/issues)
- **Improve documentation** (README, docstrings, tutorials)
- **Fix bugs** or implement features
- **Add test cases** to improve coverage
- **Share your analysis** using BOLDGenotyper
- **Help other users** by answering questions in issues

## Reporting Bugs

Before creating a bug report, please:

1. **Check existing issues** to avoid duplicates
2. **Verify the bug** with the latest version from the main branch
3. **Collect information**:
   - BOLDGenotyper version (`boldgenotyper --version`)
   - Python version (`python --version`)
   - Operating system
   - Command you ran
   - Complete error message (with traceback)
   - Relevant log file excerpts

### Bug Report Template

```markdown
**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '...'
2. With input file '...'
3. See error

**Expected behavior**
What you expected to happen.

**Environment:**
- BOLDGenotyper version: [e.g., 0.1.0]
- Python version: [e.g., 3.9.7]
- OS: [e.g., macOS 12.0, Ubuntu 20.04]
- Installation method: [conda, pip, source]

**Additional context**
- Log file excerpt
- Screenshots (if applicable)
- Input data characteristics (number of sequences, file size, etc.)
```

## Suggesting Enhancements

We welcome feature requests! Please:

1. **Check existing issues** for similar requests
2. **Provide detailed context**:
   - What problem does this solve?
   - How would you use this feature?
   - Are there examples from other tools?
3. **Be specific** about the desired behavior
4. **Consider alternatives** you've thought about

### Enhancement Request Template

```markdown
**Is your feature request related to a problem?**
A clear description of what the problem is. Ex. "I'm always frustrated when..."

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
Other solutions or features you've considered.

**Use case**
How would this feature be used in practice? What research questions would it enable?

**Additional context**
Examples from other tools, scientific literature, etc.
```

## Pull Request Process

### Before Submitting

1. **Create an issue** first to discuss major changes
2. **Follow the coding standards** (see below)
3. **Add tests** for new functionality
4. **Update documentation**:
   - Docstrings for new functions/classes
   - README if user-facing changes
   - CHANGELOG.md entry
5. **Run the test suite** and ensure all tests pass:
   ```bash
   pytest tests/
   ```
6. **Check code style**:
   ```bash
   black boldgenotyper/
   flake8 boldgenotyper/
   ```

### Submitting a Pull Request

1. **Push your branch** to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

2. **Open a Pull Request** on GitHub with:
   - Clear title describing the change
   - Reference to related issue (e.g., "Fixes #123")
   - Description of changes made
   - Testing performed
   - Screenshots (if UI/visualization changes)

3. **Address review comments**:
   - Respond to all reviewer feedback
   - Make requested changes
   - Re-request review when ready

4. **Merge requirements**:
   - All tests pass
   - Code review approval from maintainer
   - No merge conflicts with main branch
   - Changelog updated (if applicable)

### Pull Request Template

```markdown
**Description**
Brief description of what this PR does.

**Related Issue**
Fixes #(issue number)

**Type of Change**
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

**Testing**
Describe the tests you ran to verify your changes:
- [ ] Existing test suite passes
- [ ] Added new tests for this feature
- [ ] Tested manually with example data

**Checklist**
- [ ] My code follows the project's coding standards
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have updated the documentation accordingly
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes
- [ ] I have updated CHANGELOG.md
```

## Development Setup

### Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/YOUR-USERNAME/BOLDGenotyper.git
cd BOLDGenotyper

# Create development environment
conda env create -f environment.yml
conda activate boldgenotyper

# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Verify installation
boldgenotyper --version
pytest --version
```

### Using pip

```bash
# Clone and navigate to repository
git clone https://github.com/YOUR-USERNAME/BOLDGenotyper.git
cd BOLDGenotyper

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -e ".[dev,all]"

# Install external tools separately
# MAFFT: https://mafft.cbrc.jp/alignment/software/
# trimAl: http://trimal.cgenomics.org/
# FastTree: http://www.microbesonline.org/fasttree/
```

### Verify Setup

```bash
# Run tests
pytest tests/ -v

# Check code style
black --check boldgenotyper/
flake8 boldgenotyper/

# Run a simple example
boldgenotyper tests/data/test_bold.tsv --no-geo --no-report
```

## Coding Standards

### Python Style Guide

- **Follow PEP 8** for Python code style
- **Use Black** for automatic formatting (line length: 88)
- **Use flake8** for linting
- **Type hints** for function signatures (Python 3.8+ compatible)
- **Docstrings** in NumPy style for all public functions

### Example Function

```python
def example_function(
    input_data: pd.DataFrame,
    threshold: float = 0.5,
    verbose: bool = False
) -> Dict[str, Any]:
    """
    Brief description of what this function does.

    Longer description with more details about the algorithm,
    assumptions, and behavior.

    Parameters
    ----------
    input_data : pd.DataFrame
        Description of this parameter
    threshold : float, optional
        Description with units and valid range (default: 0.5)
    verbose : bool, optional
        Whether to print progress (default: False)

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'result': Description of this field
        - 'metadata': Description of this field

    Raises
    ------
    ValueError
        If threshold is not between 0 and 1

    Examples
    --------
    >>> result = example_function(df, threshold=0.7)
    >>> print(result['result'])

    Notes
    -----
    Additional information about algorithm, references, etc.

    References
    ----------
    .. [1] Author et al. (2020). Paper title. Journal.
    """
    # Implementation
    pass
```

### Code Organization

- **One class/function per logical unit**
- **Keep functions focused** (do one thing well)
- **Maximum function length**: ~50-100 lines
- **Descriptive variable names** (no single-letter except indices)
- **Comments for complex logic** (not for obvious code)
- **Constants in UPPERCASE**
- **Private functions start with underscore**

### Error Handling

```python
# Do: Raise specific exceptions with helpful messages
if threshold < 0 or threshold > 1:
    raise ValueError(
        f"Threshold must be between 0 and 1, got {threshold}. "
        "This parameter controls the similarity cutoff for clustering."
    )

# Don't: Generic errors or silent failures
if threshold < 0:
    raise Exception("Bad threshold")  # Too vague
```

## Testing Guidelines

### Test Organization

```
tests/
├── test_core.py           # Core pipeline tests
├── test_dereplication.py  # Sequence clustering tests
├── test_genotype_assignment.py
├── test_visualization.py
├── test_geographic.py
├── test_phylogenetics.py
├── test_error_handling.py # Error cases
├── conftest.py            # Pytest fixtures
└── data/                  # Test data files
    ├── test_bold.tsv
    └── test_sequences.fasta
```

### Writing Tests

```python
import pytest
from boldgenotyper import core

def test_function_basic():
    """Test basic functionality with valid input."""
    result = core.some_function(valid_input)
    assert result['status'] == 'success'
    assert len(result['data']) > 0

def test_function_edge_case():
    """Test handling of edge case (empty input)."""
    result = core.some_function([])
    assert result['status'] == 'warning'
    assert 'empty' in result['message'].lower()

def test_function_error():
    """Test that invalid input raises appropriate error."""
    with pytest.raises(ValueError, match="must be positive"):
        core.some_function(-1)

@pytest.mark.parametrize("input,expected", [
    (0.0, 'zero'),
    (0.5, 'medium'),
    (1.0, 'high'),
])
def test_function_parametrized(input, expected):
    """Test multiple input scenarios."""
    result = core.some_function(input)
    assert result['category'] == expected
```

### Running Tests

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_core.py

# Run with coverage
pytest --cov=boldgenotyper --cov-report=html tests/

# Run tests matching pattern
pytest -k "genotype"

# Run with verbose output
pytest -v tests/

# Stop on first failure
pytest -x tests/
```

### Test Coverage Goals

- **Target**: >80% code coverage
- **Critical modules**: >90% coverage (core, dereplication, genotype_assignment)
- **Test types**:
  - Unit tests for individual functions
  - Integration tests for multi-step processes
  - Error handling tests for edge cases
  - Regression tests for fixed bugs

## Documentation

### Docstring Requirements

- **All public functions** must have docstrings
- **Use NumPy style** (consistent with existing code)
- **Include**:
  - Brief one-line summary
  - Detailed description
  - Parameters with types and descriptions
  - Returns with type and description
  - Raises for exceptions
  - Examples for complex functions
  - Notes for important details
  - References for scientific methods

### README Updates

Update README.md when:
- Adding new features
- Changing user-facing behavior
- Adding new parameters
- Changing installation requirements
- Adding new outputs

### CHANGELOG

Add entry to CHANGELOG.md for:
- **Features**: New functionality
- **Fixes**: Bug fixes
- **Changes**: Behavior modifications
- **Deprecations**: Features being phased out
- **Security**: Security-related changes

Format:
```markdown
## [Unreleased]

### Added
- New feature for haplotype network construction (#123)

### Fixed
- Fixed crash when input has no coordinates (#124)

### Changed
- Improved performance of distance calculation by 50% (#125)
```

## Questions?

If you have questions about contributing:

1. **Check existing issues** for similar questions
2. **Read the documentation** thoroughly
3. **Open a discussion** on GitHub Discussions (if available)
4. **Create an issue** with the "question" label
5. **Email maintainer**: steph.smith@unc.edu

## Recognition

Contributors will be:
- Listed in [CONTRIBUTORS.md](CONTRIBUTORS.md) (if created)
- Acknowledged in release notes
- Cited in academic publications if appropriate

Thank you for contributing to BOLDGenotyper! 🎉
