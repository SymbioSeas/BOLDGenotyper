"""
Configuration Management for BOLDGenotyper

This module provides a comprehensive configuration system using dataclasses for
clean parameter management. The configuration system supports:

1. Default parameter values based on COI barcoding best practices
2. Loading configuration from YAML/JSON files
3. Environment variable overrides
4. Validation and type checking
5. Hierarchical configuration with component-specific settings

Configuration Structure:
- DereplicationConfig: Sequence clustering parameters
- GenotypeAssignmentConfig: Sample-to-genotype matching parameters
- GeographicConfig: Coordinate filtering and ocean basin assignment
- VisualizationConfig: Plot styling and output parameters
- PhylogeneticConfig: Phylogenetic tree construction parameters
- PipelineConfig: Master configuration combining all components
- TaxonomyConfig: Sequence-based species assignment thresholds and majority rules

Key Design Principles:
- Immutable configuration objects (frozen dataclasses)
- Sensible defaults based on scientific literature
- Clear documentation of parameter rationale
- Easy override mechanism for custom analyses

Example Usage:
    >>> from boldgenotyper.config import get_default_config, load_config_from_file
    >>>
    >>> # Use defaults
    >>> config = get_default_config()
    >>> print(config.dereplication.clustering_threshold)
    0.01
    >>>
    >>> # Load from file
    >>> config = load_config_from_file("my_analysis.yaml")
    >>>
    >>> # Update specific parameters
    >>> custom_config = config.update(
    ...     dereplication__clustering_threshold=0.005,
    ...     geographic__min_coordinate_quality="high"
    ... )

Author: Steph Smith (steph.smith@unc.edu)
"""

from dataclasses import dataclass, field, asdict, replace
from pathlib import Path
from typing import Optional, List, Dict, Any, Union
import os
import json
import logging

logger = logging.getLogger(__name__)

# Try to import YAML support
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    logger.debug("PyYAML not available; YAML config files not supported")


# ============================================================================
# Dereplication Configuration
# ============================================================================

@dataclass(frozen=True)
class DereplicationConfig:
    """
    Configuration for sequence dereplication and clustering.

    Parameters control how COI sequences are aligned, clustered, and converted
    to consensus sequences representing distinct genotypes.

    Attributes
    ----------
    clustering_threshold : float
        Maximum distance for sequences to cluster together (default: 0.03).
        This corresponds to 97% sequence identity threshold for
        COI-based species/genotype delimitation.

    clustering_method : str
        Hierarchical clustering linkage method (default: "average").
        Options: "average" (UPGMA), "single", "complete", "ward"

    consensus_frequency_cutoff : float
        Minimum fraction of sequences that must share a base for consensus
        calling (default: 0.7). Positions below this threshold are marked
        as ambiguous (N).

    mafft_algorithm : str
        MAFFT alignment algorithm (default: "auto").
        Options: "auto", "linsi", "ginsi", "einsi", "fftns", "fftnsi"

    trim_alignment : bool
        Whether to trim alignment with trimAl (default: True)

    trimal_method : str
        trimAl automated method (default: "automated1")

    min_sequence_length : int
        Minimum sequence length to include before alignment (default: 400 bp)

    max_n_content : float
        Maximum fraction of ambiguous bases (N) allowed (default: 0.1 = 10%)

    min_post_trim_length : int
        Minimum ungapped sequence length after trimming (default: 300 bp)
        Sequences that become too short after trimAl are removed before clustering
        to prevent fragment-based spurious genotypes.

    min_consensus_length_ratio : float
        Minimum consensus length as fraction of median consensus length (default: 0.75)
        Consensus sequences shorter than this ratio relative to the median are filtered out.
        Set to 0.0 to disable consensus length filtering.

    Notes
    -----
    The default clustering threshold of 0.03 (97% identity) allows for
    greater intraspecific variation while maintaining distinct genotype groups.

    Three-stage length filtering prevents fragment-based spurious genotypes:
    1. Pre-alignment: Remove sequences <400bp (configurable)
    2. Post-trimming: Remove sequences <300bp ungapped after trimAl
    3. Post-consensus: Remove consensus sequences <75% of median length
    """
    clustering_threshold: float = 0.03
    clustering_method: str = "average"
    consensus_frequency_cutoff: float = 0.7
    mafft_algorithm: str = "auto"
    trim_alignment: bool = True
    trimal_method: str = "automated1"
    min_sequence_length: int = 400
    max_n_content: float = 0.1
    min_post_trim_length: int = 300
    min_consensus_length_ratio: float = 0.75

    def __post_init__(self):
        """Validate configuration parameters."""
        if not 0 < self.clustering_threshold < 1:
            raise ValueError("clustering_threshold must be between 0 and 1")
        if not 0 < self.consensus_frequency_cutoff <= 1:
            raise ValueError("consensus_frequency_cutoff must be between 0 and 1")
        if self.clustering_method not in ["average", "single", "complete", "ward"]:
            raise ValueError(f"Invalid clustering_method: {self.clustering_method}")
        if self.min_sequence_length < 100:
            raise ValueError("min_sequence_length must be at least 100")
        if self.min_post_trim_length < 100:
            raise ValueError("min_post_trim_length must be at least 100")
        if not 0 <= self.min_consensus_length_ratio <= 1:
            raise ValueError("min_consensus_length_ratio must be between 0 and 1")
        if self.min_post_trim_length > self.min_sequence_length:
            logger.warning(
                f"min_post_trim_length ({self.min_post_trim_length}) is greater than "
                f"min_sequence_length ({self.min_sequence_length}). This is unusual but allowed."
            )


# ============================================================================
# Genotype Assignment Configuration
# ============================================================================

@dataclass(frozen=True)
class GenotypeAssignmentConfig:
    """
    Configuration for assigning samples to genotype groups.

    Controls how raw sequences are matched to consensus genotypes using
    edit distance calculations.

    Attributes
    ----------
    min_identity : float
        Minimum sequence identity for genotype assignment (default: 0.5).
        Samples below this threshold are flagged as unassigned.

    identity_method : str
        Method for calculating sequence identity (default: "target_based").
        Options:
        - "target_based": matches / consensus_length
          More robust to length differences and noisy 5'/3' ends.
          Recommended for COI barcoding where consensus length is canonical.
        - "classic": 1 - (edit_distance / max_length)
          Original method for backwards compatibility.
          May penalize samples with noisy ends.

    use_edlib : bool
        Prefer edlib for fast edit distance if available (default: True)

    n_threads : int
        Number of parallel threads for assignment (default: 1)

    report_ties : bool
        Report samples with ambiguous assignments (default: True)

    tie_margin : float
        Maximum identity difference between best and runner-up to call a tie (default: 0.003).
        Samples with best_identity - runner_up_identity < tie_margin are flagged for review.

    tie_threshold : float
        Minimum best identity required to consider tie detection (default: 0.95).
        Prevents flagging low-quality matches as ties.

    tie_difference_threshold : float
        DEPRECATED: Use tie_margin instead. (default: 0.01)

    Notes
    -----
    The default minimum identity of 0.5 (50%) provides a permissive threshold
    for genotype assignment, allowing for:
    - Greater intraspecific variation
    - Degraded or low-quality sequences
    - Differences between consensus and raw sequences

    Identity Method Selection:
    - Use "target_based" (default) when consensus sequences are canonical
      references and samples may have variable quality at 5'/3' ends
    - Use "classic" for backwards compatibility or when sequences are
      expected to have uniform length and quality

    Tie Detection:
    - tie_margin: Controls sensitivity of tie detection (default 0.003 = 0.3%)
    - tie_threshold: Minimum quality for considering ties (default 0.95 = 95%)
    - Ties are only called when: best_identity >= tie_threshold AND
      (best_identity - runner_up_identity) < tie_margin
    """
    min_identity: float = 0.5
    identity_method: str = "target_based"
    use_edlib: bool = True
    n_threads: int = 1
    report_ties: bool = True
    tie_margin: float = 0.003
    tie_threshold: float = 0.95
    tie_difference_threshold: float = 0.01  # Deprecated

    def __post_init__(self):
        """Validate configuration parameters."""
        if not 0 < self.min_identity <= 1:
            raise ValueError("min_identity must be between 0 and 1")
        if self.n_threads < 1:
            raise ValueError("n_threads must be at least 1")
        if self.identity_method not in ["target_based", "classic"]:
            raise ValueError(
                f"identity_method must be 'target_based' or 'classic', "
                f"got '{self.identity_method}'"
            )
            
            
# ============================================================================
# Taxonomy (Sequence-based Species Assignment) Configuration
# ============================================================================

@dataclass(frozen=True)
class TaxonomyConfig:
    """
    Configuration for sequence-based taxonomic assignment (samples & consensus).

    These parameters control BLAST/VSEARCH-style classification thresholds,
    tie-breaking behavior, and majority-vote consolidation.

    Attributes
    ----------
    min_identity_pct : float
        Minimum percent identity required to allow species-level assignment
        from sequence similarity (e.g., 98.5 for COI).

    min_query_cov_pct : float
        Minimum percent of the query aligned (coverage) to consider a hit valid.

    top2_min_delta_pct : float
        Minimum absolute percent-identity margin between the best hit and the
        best conflicting species hit to accept a species-level call.

    majority_species_threshold : float
        Fraction of members in a consensus_group that must share the same
        reported species to call a metadata majority species.

    allow_ambiguous_N_pct : float
        Maximum percent of ambiguous bases (N) tolerated in the *query* when
        making species calls. Above this, downgrade to genus or unassigned.

    classifier : str
        Which backend you intend to use for sequence classification.
        Options: "blastn", "vsearch" (currently informational).

    db_name : str
        Human-readable database tag for provenance in reports.

    locus : Optional[str]
        Locus/marker name (e.g., "COI-5P"), purely informational.
    """
    min_identity_pct: float = 98.5
    min_query_cov_pct: float = 90.0
    top2_min_delta_pct: float = 0.5
    majority_species_threshold: float = 0.70
    allow_ambiguous_N_pct: float = 2.0
    classifier: str = "blastn"
    db_name: str = "BOLD_COI_2025-10-01"
    locus: Optional[str] = "COI-5P"

    def __post_init__(self):
        if not (0 < self.min_identity_pct <= 100):
            raise ValueError("min_identity_pct must be in (0, 100]")
        if not (0 < self.min_query_cov_pct <= 100):
            raise ValueError("min_query_cov_pct must be in (0, 100]")
        if not (0 <= self.top2_min_delta_pct <= 5):
            raise ValueError("top2_min_delta_pct should be reasonable (0–5)")
        if not (0.5 <= self.majority_species_threshold <= 1.0):
            raise ValueError("majority_species_threshold must be in [0.5, 1.0]")
        if not (0 <= self.allow_ambiguous_N_pct <= 100):
            raise ValueError("allow_ambiguous_N_pct must be in [0, 100]")
        if self.classifier not in {"blastn", "vsearch"}:
            raise ValueError('classifier must be "blastn" or "vsearch"')


# ============================================================================
# Geographic Configuration
# ============================================================================

@dataclass(frozen=True)
class GeographicConfig:
    """
    Configuration for geographic data processing and ocean basin assignment.

    Controls coordinate filtering, validation, and spatial analysis using
    the Global Oceans and Seas (GOaS) dataset.

    Attributes
    ----------
    goas_shapefile_path : Optional[Path]
        Path to GOaS shapefile for ocean basin assignment.
        Default: looks for GOaS_v1_20211214/goas_v01.shp in package directory

    exclude_centroids : bool
        Exclude samples with centroid coordinates (default: True)
        Critical for accurate biogeographic analysis.

    exclude_zero_coords : bool
        Exclude coordinates at (0, 0) (default: True)

    min_coordinate_precision : Optional[int]
        Minimum decimal places for coordinates (default: None)

    validate_marine : bool
        Validate that coordinates fall in marine areas (default: False)

    ocean_basins : List[str]
        List of ocean basin names to recognize

    Notes
    -----
    Coordinate filtering is CRITICAL for accurate ocean basin assignment.
    Many BOLD records contain country-level centroids rather than actual
    collection locations, which can lead to incorrect biogeographic conclusions.

    The GOaS (Global Oceans and Seas) dataset is from:
    https://www.marineregions.org/downloads.php
    """
    goas_shapefile_path: Optional[Path] = None
    exclude_centroids: bool = True
    exclude_zero_coords: bool = True
    min_coordinate_precision: Optional[int] = None
    validate_marine: bool = False
    ocean_basins: List[str] = field(default_factory=lambda: [
        "North Atlantic Ocean",
        "South Atlantic Ocean",
        "North Pacific Ocean",
        "South Pacific Ocean",
        "Indian Ocean",
        "Southern Ocean",
        "Mediterranean Sea",
        "South China Sea",
    ])

    def __post_init__(self):
        """Validate and set default paths."""
        # Convert string path to Path object if needed
        if self.goas_shapefile_path is not None and isinstance(self.goas_shapefile_path, str):
            object.__setattr__(self, 'goas_shapefile_path', Path(self.goas_shapefile_path))

        # Set default GOaS path if not provided
        if self.goas_shapefile_path is None:
            # Look for GOaS data in package directory
            package_dir = Path(__file__).parent.parent
            default_path = package_dir / "GOaS_v1_20211214" / "goas_v01.shp"
            object.__setattr__(self, 'goas_shapefile_path', default_path)


    def validate_goas_data(self) -> bool:
        """
        Check if GOaS shapefile exists and is readable.

        Provides helpful error message if not found.
        """
        if self.goas_shapefile_path is None:
            logger.warning("No GOaS shapefile path configured")
            return False

        if not self.goas_shapefile_path.exists():
            logger.error(
                f"GOaS shapefile not found: {self.goas_shapefile_path}\n\n"
                "To download the GOaS reference data, run:\n"
                "    python setup_goas.py\n\n"
                "Or download manually from:\n"
                "    https://www.marineregions.org/downloads.php"
            )
            return False

        return True

# ============================================================================
# Visualization Configuration
# ============================================================================

@dataclass(frozen=True)
class VisualizationConfig:
    """
    Configuration for publication-ready figure generation.

    Controls plot styling, colors, output formats, and dimensions.

    Attributes
    ----------
    color_palette : str
        Matplotlib/seaborn color palette (default: "colorblind")

    reference_colors : List[str]
        Specific colors for first N genotypes to match reference analysis

    figure_dpi : int
        Resolution for raster outputs (default: 300)

    figure_format : List[str]
        Output formats (default: ["png", "pdf"])

    map_projection : str
        Cartopy projection for maps (default: "PlateCarree")

    map_figsize : tuple
        Figure size for maps in inches (default: (12, 6))

    barplot_figsize : tuple
        Figure size for bar plots in inches (default: (10, 6))

    font_family : str
        Font family for plots (default: "sans-serif")

    font_size : int
        Base font size in points (default: 11)

    show_bootstrap_threshold : int
        Minimum bootstrap value to display on trees (default: 70)

    facet_by : str
        Faceting mode for multi-panel plots (default: "species")
        Options: "species" or "genotype"
        - "species": One facet per species, multiple genotypes per facet (reduces file size)
        - "genotype": One facet per genotype (more granular)

    map_buffer_degrees : float
        Buffer margin in degrees for faceted map zoom (default: 20.0)
        Extends map bounds beyond data points for geographic context.
        Automatically constrained to valid lat/lon ranges [-90,90] × [-180,180]

    show_unknown_geography_annotation : bool
        Show count of samples without coordinates on faceted maps (default: True)

    show_scale_bar : bool
        Show scale bar indicating point size → sample count on maps (default: True)

    Notes
    -----
    Reference colors are chosen for consistency with the original analysis:
    - Purple (#9D7ABE)
    - Teal (#5AB4AC)
    - Yellow (#F2CC8F)

    Faceting by species significantly reduces output file sizes while maintaining
    all information (e.g., 2 species with 5 genotypes each = 2 facets, not 10).
    """
    color_palette: str = "colorblind"
    reference_colors: List[str] = field(default_factory=lambda: [
        "#9D7ABE",  # Purple
        "#5AB4AC",  # Teal
        "#F2CC8F",  # Yellow
    ])
    figure_dpi: int = 300
    figure_format: List[str] = field(default_factory=lambda: ["png", "pdf"])
    map_projection: str = "PlateCarree"
    map_figsize: tuple = (12, 6)
    barplot_figsize: tuple = (10, 6)
    font_family: str = "sans-serif"
    font_size: int = 11
    show_bootstrap_threshold: int = 70
    facet_by: str = "species"
    map_buffer_degrees: float = 20.0
    show_unknown_geography_annotation: bool = True
    show_scale_bar: bool = True

    def __post_init__(self):
        """Validate configuration parameters."""
        if self.figure_dpi < 72:
            raise ValueError("figure_dpi must be at least 72")
        if self.font_size < 6:
            raise ValueError("font_size must be at least 6")
        if self.facet_by not in ["species", "genotype"]:
            raise ValueError("facet_by must be 'species' or 'genotype'")
        if self.map_buffer_degrees < 0:
            raise ValueError("map_buffer_degrees must be non-negative")
        if self.map_buffer_degrees > 90:
            logger.warning(
                f"map_buffer_degrees ({self.map_buffer_degrees}) is very large. "
                "This may result in overly zoomed-out maps."
            )


# ============================================================================
# Phylogenetic Configuration
# ============================================================================

@dataclass(frozen=True)
class PhylogeneticConfig:
    """
    Configuration for phylogenetic tree construction.

    Controls alignment and tree building for optional phylogenetic analysis.

    Attributes
    ----------
    build_tree : bool
        Whether to build phylogenetic tree (default: False)

    outgroup_fasta : Optional[Path]
        Path to outgroup sequences for rooting (default: None)

    substitution_model : str
        PhyML substitution model (default: "GTR")

    n_bootstrap : int
        Number of bootstrap replicates (default: 1000)

    tree_algorithm : str
        Tree building method (default: "phyml")
        Options: "phyml", "fasttree", "raxml"

    midpoint_root : bool
        Use midpoint rooting if no outgroup (default: True)

    Notes
    -----
    Phylogenetic analysis is optional because:
    - Outgroup selection requires taxonomic expertise
    - Not all users need phylogenetic trees
    - Computational time can be significant

    Three modes of operation:
    1. No phylogeny (build_tree=False): Fastest, genotype ID only
    2. Midpoint rooting (build_tree=True, no outgroup): Visualize relationships
    3. Outgroup rooting (build_tree=True, with outgroup): Publication quality
    """
    build_tree: bool = False
    outgroup_fasta: Optional[Path] = None
    substitution_model: str = "GTR"
    n_bootstrap: int = 1000
    tree_algorithm: str = "phyml"
    midpoint_root: bool = True

    def __post_init__(self):
        """Validate configuration parameters."""
        if self.outgroup_fasta is not None and isinstance(self.outgroup_fasta, str):
            object.__setattr__(self, 'outgroup_fasta', Path(self.outgroup_fasta))
        if self.n_bootstrap < 0:
            raise ValueError("n_bootstrap must be non-negative")
        if self.tree_algorithm not in ["phyml", "fasttree", "raxml"]:
            raise ValueError(f"Invalid tree_algorithm: {self.tree_algorithm}")


# ============================================================================
# Master Pipeline Configuration
# ============================================================================

@dataclass(frozen=True)
class PipelineConfig:
    """
    Master configuration for the complete BOLDGenotyper pipeline.

    Combines all component-specific configurations and adds global settings.

    Attributes
    ----------
    dereplication : DereplicationConfig
        Sequence clustering configuration

    genotype_assignment : GenotypeAssignmentConfig
        Sample-to-genotype matching configuration

    geographic : GeographicConfig
        Coordinate filtering and ocean basin configuration

    visualization : VisualizationConfig
        Figure generation configuration
        
    taxonomy : TaxonomyConfig
        Taxonomy configuration

    phylogenetic : PhylogeneticConfig
        Phylogenetic tree construction configuration

    log_level : str
        Logging level (default: "INFO")

    n_threads : int
        Global thread count (default: 4, overrides component settings)

    output_dir : Path
        Base output directory (default: "results")

    keep_intermediates : bool
        Keep intermediate files (alignments, etc.) (default: False)

    overwrite_existing : bool
        Overwrite existing output files (default: False)
    """
    dereplication: DereplicationConfig = field(default_factory=DereplicationConfig)
    genotype_assignment: GenotypeAssignmentConfig = field(default_factory=GenotypeAssignmentConfig)
    geographic: GeographicConfig = field(default_factory=GeographicConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    taxonomy: TaxonomyConfig = field(default_factory=TaxonomyConfig)
    phylogenetic: PhylogeneticConfig = field(default_factory=PhylogeneticConfig)
    log_level: str = "INFO"
    n_threads: int = 4
    output_dir: Path = field(default_factory=lambda: Path("results"))
    keep_intermediates: bool = False
    overwrite_existing: bool = False

    def __post_init__(self):
        """Validate and normalize configuration."""
        # Convert string paths to Path objects
        if isinstance(self.output_dir, str):
            object.__setattr__(self, 'output_dir', Path(self.output_dir))

        # Validate log level
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        if self.log_level.upper() not in valid_levels:
            raise ValueError(f"log_level must be one of {valid_levels}")

        # Validate threads
        if self.n_threads < 1:
            raise ValueError("n_threads must be at least 1")

    def update(self, **kwargs) -> 'PipelineConfig':
        """
        Create a new configuration with updated values.

        Supports nested updates using double underscore notation:
        config.update(dereplication__clustering_threshold=0.005)

        Parameters
        ----------
        **kwargs
            Configuration parameters to update. Use double underscore
            for nested parameters (e.g., dereplication__clustering_threshold)

        Returns
        -------
        PipelineConfig
            New configuration object with updates

        Examples
        --------
        >>> config = get_default_config()
        >>> new_config = config.update(
        ...     n_threads=4,
        ...     dereplication__clustering_threshold=0.005,
        ...     visualization__figure_dpi=600
        ... )
        """
        # Separate top-level and nested updates
        top_level = {}
        nested = {}

        for key, value in kwargs.items():
            if '__' in key:
                component, param = key.split('__', 1)
                if component not in nested:
                    nested[component] = {}
                nested[component][param] = value
            else:
                top_level[key] = value

        # Update nested configurations
        if nested:
            for component, updates in nested.items():
                current = getattr(self, component)
                updated = replace(current, **updates)
                top_level[component] = updated

        # Create new configuration
        return replace(self, **top_level)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary.

        Returns
        -------
        Dict[str, Any]
            Configuration as nested dictionary
        """
        return asdict(self)

    def to_yaml(self, output_path: Union[str, Path]) -> None:
        """
        Save configuration to YAML file.

        Parameters
        ----------
        output_path : Union[str, Path]
            Output file path

        Raises
        ------
        ImportError
            If PyYAML is not installed
        """
        if not YAML_AVAILABLE:
            raise ImportError("PyYAML is required to save YAML config files")

        config_dict = self.to_dict()

        # Convert Path objects to strings for YAML serialization
        config_dict = _convert_paths_to_strings(config_dict)

        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

        logger.info(f"Configuration saved to {path}")

    def to_json(self, output_path: Union[str, Path]) -> None:
        """
        Save configuration to JSON file.

        Parameters
        ----------
        output_path : Union[str, Path]
            Output file path
        """
        config_dict = self.to_dict()

        # Convert Path objects to strings for JSON serialization
        config_dict = _convert_paths_to_strings(config_dict)

        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, 'w') as f:
            json.dump(config_dict, f, indent=2)

        logger.info(f"Configuration saved to {path}")


# ============================================================================
# Helper Functions
# ============================================================================

def get_default_config() -> PipelineConfig:
    """
    Get default pipeline configuration.

    Returns
    -------
    PipelineConfig
        Default configuration with recommended parameters

    Examples
    --------
    >>> config = get_default_config()
    >>> print(config.dereplication.clustering_threshold)
    0.01
    """
    return PipelineConfig()


def load_config_from_file(config_path: Union[str, Path]) -> PipelineConfig:
    """
    Load configuration from YAML or JSON file.

    Automatically detects file format based on extension.

    Parameters
    ----------
    config_path : Union[str, Path]
        Path to configuration file (.yaml, .yml, or .json)

    Returns
    -------
    PipelineConfig
        Loaded configuration

    Raises
    ------
    FileNotFoundError
        If configuration file doesn't exist
    ValueError
        If file format is not supported

    Examples
    --------
    >>> config = load_config_from_file("my_analysis.yaml")
    """
    path = Path(config_path)

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    suffix = path.suffix.lower()

    if suffix in ['.yaml', '.yml']:
        if not YAML_AVAILABLE:
            raise ImportError("PyYAML is required to load YAML config files")
        return _load_yaml_config(path)
    elif suffix == '.json':
        return _load_json_config(path)
    else:
        raise ValueError(f"Unsupported config file format: {suffix}")


def _load_yaml_config(path: Path) -> PipelineConfig:
    """Load configuration from YAML file."""
    with open(path, 'r') as f:
        config_dict = yaml.safe_load(f)

    logger.info(f"Loaded configuration from {path}")
    return _dict_to_config(config_dict)


def _load_json_config(path: Path) -> PipelineConfig:
    """Load configuration from JSON file."""
    with open(path, 'r') as f:
        config_dict = json.load(f)

    logger.info(f"Loaded configuration from {path}")
    return _dict_to_config(config_dict)


def _dict_to_config(config_dict: Dict[str, Any]) -> PipelineConfig:
    """
    Convert dictionary to PipelineConfig object.

    Handles nested configuration structures and type conversions.
    """
    # Convert string paths back to Path objects
    config_dict = _convert_strings_to_paths(config_dict)

    # Extract nested configurations
    nested_configs = {}

    if 'dereplication' in config_dict:
        nested_configs['dereplication'] = DereplicationConfig(**config_dict.pop('dereplication'))

    if 'genotype_assignment' in config_dict:
        nested_configs['genotype_assignment'] = GenotypeAssignmentConfig(
            **config_dict.pop('genotype_assignment')
        )

    if 'geographic' in config_dict:
        nested_configs['geographic'] = GeographicConfig(**config_dict.pop('geographic'))

    if 'visualization' in config_dict:
        nested_configs['visualization'] = VisualizationConfig(**config_dict.pop('visualization'))

    if 'phylogenetic' in config_dict:
        nested_configs['phylogenetic'] = PhylogeneticConfig(**config_dict.pop('phylogenetic'))
        
    if 'taxonomy' in config_dict:
        nested_configs['taxonomy'] = TaxonomyConfig(**config_dict.pop('taxonomy'))

    # Create pipeline config
    return PipelineConfig(**nested_configs, **config_dict)


def _convert_paths_to_strings(obj: Any) -> Any:
    """Recursively convert Path objects to strings for serialization."""
    if isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: _convert_paths_to_strings(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_convert_paths_to_strings(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(_convert_paths_to_strings(item) for item in obj)
    else:
        return obj


def _convert_strings_to_paths(obj: Any) -> Any:
    """Recursively convert path strings back to Path objects."""
    if isinstance(obj, dict):
        # Known path fields
        path_fields = [
            'goas_shapefile_path', 'outgroup_fasta', 'output_dir'
        ]

        result = {}
        for k, v in obj.items():
            if k in path_fields and v is not None:
                result[k] = Path(v)
            else:
                result[k] = _convert_strings_to_paths(v)
        return result
    elif isinstance(obj, list):
        return [_convert_strings_to_paths(item) for item in obj]
    else:
        return obj


def load_config_from_env() -> Dict[str, Any]:
    """
    Load configuration overrides from environment variables.

    Environment variables should be prefixed with BOLDGENOTYPER_
    and use double underscores for nesting:

    BOLDGENOTYPER_DEREPLICATION__CLUSTERING_THRESHOLD=0.005
    BOLDGENOTYPER_N_THREADS=4

    Returns
    -------
    Dict[str, Any]
        Configuration overrides from environment

    Examples
    --------
    >>> import os
    >>> os.environ['BOLDGENOTYPER_N_THREADS'] = '4'
    >>> env_config = load_config_from_env()
    >>> config = get_default_config().update(**env_config)
    """
    prefix = "BOLDGENOTYPER_"
    overrides = {}

    for key, value in os.environ.items():
        if key.startswith(prefix):
            # Remove prefix and convert to lowercase
            config_key = key[len(prefix):].lower()

            # Try to parse value
            parsed_value = _parse_env_value(value)
            overrides[config_key] = parsed_value

    if overrides:
        logger.debug(f"Loaded {len(overrides)} configuration overrides from environment")

    return overrides


def _parse_env_value(value: str) -> Any:
    """Parse environment variable value to appropriate type."""
    # Boolean
    if value.lower() in ['true', 'yes', '1']:
        return True
    if value.lower() in ['false', 'no', '0']:
        return False

    # Integer
    try:
        return int(value)
    except ValueError:
        pass

    # Float
    try:
        return float(value)
    except ValueError:
        pass

    # String
    return value


def validate_config(config: PipelineConfig) -> List[str]:
    """
    Validate configuration and return list of warnings.

    Checks for common issues like missing files, unusual parameter values, etc.

    Parameters
    ----------
    config : PipelineConfig
        Configuration to validate

    Returns
    -------
    List[str]
        List of warning messages (empty if no issues)

    Examples
    --------
    >>> config = get_default_config()
    >>> warnings = validate_config(config)
    >>> for warning in warnings:
    ...     print(f"Warning: {warning}")
    """
    warnings = []

    # Check GOaS shapefile exists
    if config.geographic.goas_shapefile_path:
        if not config.geographic.goas_shapefile_path.exists():
            warnings.append(
                f"GOaS shapefile not found: {config.geographic.goas_shapefile_path}. "
                "Ocean basin assignment will not be available."
            )

    # Check outgroup file exists
    if config.phylogenetic.build_tree and config.phylogenetic.outgroup_fasta:
        if not config.phylogenetic.outgroup_fasta.exists():
            warnings.append(
                f"Outgroup FASTA not found: {config.phylogenetic.outgroup_fasta}"
            )

    # Check for unusual clustering threshold
    if config.dereplication.clustering_threshold > 0.05:
        warnings.append(
            f"Clustering threshold ({config.dereplication.clustering_threshold}) is quite high. "
            "This may result in very broad genotype groups."
        )

    if config.dereplication.clustering_threshold < 0.001:
        warnings.append(
            f"Clustering threshold ({config.dereplication.clustering_threshold}) is very low. "
            "This may result in excessive splitting of genotypes."
        )

    # Check for low minimum identity
    if config.genotype_assignment.min_identity < 0.80:
        warnings.append(
            f"Minimum identity ({config.genotype_assignment.min_identity}) is quite low. "
            "This may result in poor genotype assignments."
        )

    # Check thread count
    if config.n_threads > os.cpu_count():
        warnings.append(
            f"Thread count ({config.n_threads}) exceeds available CPUs ({os.cpu_count()})"
        )
        
    # taxonomy sanity checks
    tx = config.taxonomy
    if tx.min_identity_pct < 97:
        warnings.append(
            f"Taxonomy min_identity_pct ({tx.min_identity_pct}) is low for COI; "
            "you may be genus-level calls promoted to species."
        )
    if tx.min_query_cov_pct < 80:
        warnings.append(
            f"Taxonomy min_query_cov_pct ({tx.min_query_cov_pct}) is low; "
            "partial hits may drive incorrect species assignments."
        )
    if tx.majority_species_threshold < 0.6:
        warnings.append(
            f"Majority treshold ({tx.majority_species_threshold}) is lenient; "
            "consider 0.7 for conservative metadata-majority species naming."
        )

    return warnings


# ============================================================================
# Configuration Templates
# ============================================================================

def create_config_template(output_path: Union[str, Path], format: str = "yaml") -> None:
    """
    Create a configuration template file with comments.

    Parameters
    ----------
    output_path : Union[str, Path]
        Output file path
    format : str
        File format: "yaml" or "json" (default: "yaml")

    Examples
    --------
    >>> create_config_template("my_config.yaml")
    """
    config = get_default_config()

    if format.lower() == "yaml":
        config.to_yaml(output_path)
    elif format.lower() == "json":
        config.to_json(output_path)
    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Created configuration template: {output_path}")
