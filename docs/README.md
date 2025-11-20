# BOLDGenotyper Documentation

This directory contains the source files for BOLDGenotyper's documentation, built with Sphinx.

## Building the Documentation

### Prerequisites

Install documentation dependencies:

```bash
pip install sphinx sphinx-rtd-theme sphinx-autodoc-typehints
```

Or install with the dev extras:

```bash
pip install -e ".[dev]"
```

### Build HTML Documentation

From the `docs` directory:

```bash
cd docs
make html
```

The built documentation will be in `docs/build/html/`. Open `docs/build/html/index.html` in your browser.

### Build Other Formats

```bash
make latexpdf  # PDF via LaTeX
make epub      # EPUB format
make man       # Man pages
```

### Clean Build Files

```bash
make clean
```

## Documentation Structure

```
docs/
├── source/
│   ├── conf.py              # Sphinx configuration
│   ├── index.rst            # Main documentation page
│   ├── installation.rst     # Installation guide
│   ├── quickstart.rst       # Quick start tutorial
│   ├── usage.rst            # Usage guide
│   ├── troubleshooting.rst  # Troubleshooting
│   ├── api/                 # API reference
│   │   ├── cli.rst
│   │   ├── config.rst
│   │   ├── dereplication.rst
│   │   └── ...
│   ├── contributing.rst     # Contribution guidelines
│   └── changelog.rst        # Version history
├── build/                   # Generated documentation (git-ignored)
└── Makefile                 # Build commands
```

## Auto-generating API Documentation

To auto-generate API documentation from docstrings:

```bash
sphinx-apidoc -o docs/source/api boldgenotyper
```

## Hosting on Read the Docs

1. Go to [readthedocs.org](https://readthedocs.org/)
2. Import the BOLDGenotyper repository
3. Configure webhook in GitHub settings
4. Documentation will auto-build on every push

## Style Guide

- Use Google-style docstrings in code
- Follow reStructuredText (RST) syntax in .rst files
- Include code examples in docstrings
- Cross-reference other functions with ``:func:`module.function```
- Add new pages to the toctree in `index.rst`

## Contributing

To improve documentation:
1. Edit .rst files in `docs/source/`
2. Build locally to test: `make html`
3. Submit pull request with changes

## Resources

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- [Read the Docs Guide](https://docs.readthedocs.io/)
