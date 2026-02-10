"""
BOLDGenotyper: Automated COI Haplotype Discovery and Biogeographic Analysis

BOLDGenotyper is a Python package for automated COI haplotype discovery and
biogeographic analysis from BOLD (Barcode of Life Database) data. The package
enables researchers to identify genetic partitioning patterns in any organism
with publicly available COI barcode sequences.

Core functionality includes:
- COI sequence orientation normalization and ORF validation
- Dynamic quality control with median-based thresholds
- Haplotype discovery using exact sequence variants (ESVs)
- Core-region extraction from variable-length sequences
- Suspect haplotype flagging based on distance and ORF quality
- Haplotype assignment using sequence similarity
- Geographic coordinate filtering and ocean basin assignment
- Publication-ready visualization generation
- Optional phylogenetic analysis with outgroup rooting

Workflow:
1. Quality Control: Orientation normalization, ORF validation, dynamic filtering
2. Haplotype Discovery: Align → extract core region → identify ESVs → flag suspects
3. Haplotype Assignment: Match samples to haplotypes
4. Geographic Analysis: Filter coordinates, assign ocean basins
5. Visualization & Phylogenetics: Generate plots and optional trees

This implements the ESV (Exact Sequence Variant) approach recommended by
Porter & Hajibabaei (2020) for high-resolution COI barcoding.

Developed to support the manuscript: "Ocean basin-scale genetic partitioning
in Sphyrna lewini revealed through COI sequence analysis"

Author: Steph Smith (steph.smith@unc.edu)
Institution: University of North Carolina, Institute of Marine Sciences
"""

__version__ = "0.1.0"
__author__ = "Steph Smith"
__email__ = "steph.smith@unc.edu"

# Import main modules for easy access
from . import config
from . import utils
from . import quality_control
from . import dereplication
from . import haplotype_assignment
from . import haplotype_query
from . import cluster_diagnostics
from . import metadata
from . import geographic
from . import visualization
from . import phylogenetics
from . import comparative_analysis
from . import parameter_sweep
from . import metadata_enrichment
from . import popgen_export
from . import metadata_analysis

# Legacy alias for backward compatibility
from . import haplotype_assignment as genotype_assignment

__all__ = [
    "config",
    "utils",
    "quality_control",
    "dereplication",
    "haplotype_assignment",
    "haplotype_query",
    "cluster_diagnostics",
    "metadata",
    "geographic",
    "visualization",
    "phylogenetics",
    "comparative_analysis",
    "parameter_sweep",
    "metadata_enrichment",
    "popgen_export",
    "metadata_analysis",
    "genotype_assignment",  # Legacy alias
]
