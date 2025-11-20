from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README.md for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name="boldgenotyper",
    version="0.1.0",

    # Descriptions
    description="Automated COI sequence genotyping and biogeographic analysis from BOLD database data",
    long_description=long_description,
    long_description_content_type="text/markdown",

    # Author information
    author="Steph Smith",
    author_email="steph.smith@unc.edu",

    # URLs
    url="https://github.com/SymbioSeas/BOLDGenotyper",
    project_urls={
        "Bug Reports": "https://github.com/SymbioSeas/BOLDGenotyper/issues",
        "Source": "https://github.com/SymbioSeas/BOLDGenotyper",
        "Documentation": "https://github.com/SymbioSeas/BOLDGenotyper/blob/main/README.md",
        "Changelog": "https://github.com/SymbioSeas/BOLDGenotyper/blob/main/CHANGELOG.md",
    },

    # License
    license="MIT",

    # Package discovery
    packages=find_packages(exclude=["tests", "tests.*", "docs", "examples"]),

    # Include non-Python files specified in MANIFEST.in
    include_package_data=True,

    # Python version requirement
    python_requires=">=3.8",

    # Core dependencies
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "jinja2>=3.0.0",
    ],

    # Optional dependencies for specific features
    extras_require={
        "geo": [
            "geopandas>=0.10.0",
            "cartopy>=0.20.0",
            "shapely>=1.8.0",
        ],
        "phylo": [
            "ete3>=3.1.2",
            "dendropy>=4.5.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
            "sphinx>=4.5.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
        "all": [
            "geopandas>=0.10.0",
            "cartopy>=0.20.0",
            "shapely>=1.8.0",
            "ete3>=3.1.2",
            "dendropy>=4.5.0",
        ],
    },

    # Command-line interface
    entry_points={
        'console_scripts': [
            'boldgenotyper=boldgenotyper.cli:main',
        ],
    },

    # PyPI classifiers
    classifiers=[
        # Development status
        "Development Status :: 3 - Alpha",

        # Intended audience
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",

        # Topic areas
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",

        # License
        "License :: OSI Approved :: MIT License",

        # Supported Python versions
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",

        # Operating systems
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows :: Windows Subsystem for Linux",

        # Other
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
        "phylogeography",
        "biogeography",
        "BOLD",
        "population genetics",
        "sequence analysis",
        "ocean basins",
        "marine biology",
    ],

    # Minimum setuptools version
    setup_requires=["setuptools>=45.0"],

    # Specify that this package is zip-safe (or not)
    zip_safe=False,
)
