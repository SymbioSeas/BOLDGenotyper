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
- Geographic region assignment using any polygon shapefile (ocean basins,
  freshwater basins, ecoregions, watersheds, biomes, etc.)
- Publication-ready visualization generation
- Optional phylogenetic analysis with outgroup rooting

Workflow:
1. Quality Control: Orientation normalization, ORF validation, dynamic filtering
2. Haplotype Discovery: Align → extract core region → identify ESVs → flag suspects
3. Haplotype Assignment: Match samples to haplotypes
4. Geographic Analysis: Assign samples to regions via shapefile (default: GOaS ocean basins)
5. Visualization & Phylogenetics: Generate plots and optional trees

This implements the ESV (Exact Sequence Variant) approach recommended by
Porter & Hajibabaei (2020) for high-resolution COI barcoding.

Author: Steph Smith (symbioseas@outlook.com)
Institution: University of North Carolina, Institute of Marine Sciences
"""

__version__ = "1.0.0"
__author__ = "Steph Smith"
__email__ = "symbioseas@outlook.com"

# Import main modules for easy access
from . import config
from . import utils
from . import quality_control
from . import dereplication
from . import haplotype_assignment
from . import haplotype_query
from . import metadata
from . import geographic
from . import geographic_enhancement
from . import visualization
from . import phylogenetics
from . import divergence_analysis
from . import species_analysis
from . import metadata_summary
from . import reports
from . import plot_export
from . import plot_regeneration
from . import popgen_export
from . import parameter_sweep
from . import msa_visualization

__all__ = [
    "config",
    "utils",
    "quality_control",
    "dereplication",
    "haplotype_assignment",
    "haplotype_query",
    "metadata",
    "geographic",
    "geographic_enhancement",
    "visualization",
    "phylogenetics",
    "divergence_analysis",
    "species_analysis",
    "metadata_summary",
    "reports",
    "plot_export",
    "plot_regeneration",
    "popgen_export",
    "parameter_sweep",
    "msa_visualization",
]
