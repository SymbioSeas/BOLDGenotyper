from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README.md for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name="boldgenotyper",
    version="1.0.0",

    # Descriptions
    description="Automated COI sequence genotyping and biogeographic analysis from BOLD database data",
    long_description=long_description,
    long_description_content_type="text/markdown",

    # Author information
    author="Steph Smith",
    author_email="symbioseas@outlook.com",

    # URLs
    url="https://github.com/SymbioSeas/BOLDGenotyper",
    project_urls={
        "Bug Reports": "https://github.com/SymbioSeas/BOLDGenotyper/issues",
        "Source": "https://github.com/SymbioSeas/BOLDGenotyper",
        "Documentation": "https://github.com/SymbioSeas/BOLDGenotyper/blob/main/README.md",
    },

    # License
    license="MIT",

    # Package discovery
    packages=find_packages(exclude=["tests", "tests.*", "docs"]),

    # Include non-Python files specified in MANIFEST.in
    include_package_data=True,

    # Python version requirement
    python_requires=">=3.9",

    # Core dependencies
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "edlib>=1.3.9",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "jinja2>=3.0.0",
        "pyyaml>=5.0.0",
        "geopandas>=0.10.0",
        "cartopy>=0.20.0",
        "shapely>=1.8.0",
        "pyproj>=3.0.0",
    ],

    # Optional dependencies
    extras_require={
        "msa": [
            "pymsaviz>=0.4.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
        ],
    },

    # Command-line interface
    entry_points={
        'console_scripts': [
            'boldgenotyper=boldgenotyper.cli:main',
            'boldgenotyper-sweep=boldgenotyper.cli:main_sweep',
            'boldgenotyper-query=boldgenotyper.cli:main_query',
        ],
    },

    # PyPI classifiers
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Natural Language :: English",
        "Typing :: Typed",
    ],

    # Keywords for PyPI search
    keywords=[
        "bioinformatics",
        "DNA barcoding",
        "COI",
        "cytochrome oxidase I",
        "genotyping",
        "haplotype",
        "phylogeography",
        "biogeography",
        "BOLD",
        "population genetics",
    ],

    zip_safe=False,
)
