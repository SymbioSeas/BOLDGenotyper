"""
BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis

BOLDGenotyper is a Python package for automated COI sequence genotyping and
biogeographic analysis from BOLD (Barcode of Life Database) data. The package
enables researchers to identify genetic partitioning patterns in any organism
with publicly available COI barcode sequences.

Core functionality includes:
- Sequence dereplication and consensus generation
- Genotype assignment using sequence clustering
- Geographic coordinate filtering and ocean basin assignment
- Publication-ready visualization generation
- Optional phylogenetic analysis

Developed to support the manuscript: "Ocean basin-scale genetic partitioning
in Sphyrna lewini revealed through COI sequence analysis"

Author: Steph Smith (steph.smith@unc.edu)
Institution: University of North Carolina, Institute of Marine Sciences
"""

__version__ = "0.1.0"
__author__ = "Steph Smith"
__email__ = "steph.smith@unc.edu"

# Import main modules for easy access
from . import core
from . import dereplication
from . import metadata
from . import geographic
from . import visualization
from . import phylogenetics
from . import utils

__all__ = [
    "core",
    "dereplication",
    "metadata",
    "geographic",
    "visualization",
    "phylogenetics",
    "utils",
]
