# BOLDGenotyper API Reference

**Version**: 0.1.0
**Last Updated**: 2025-11-25

This document provides guidance for using BOLDGenotyper programmatically in Python scripts and custom workflows.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Modules](#core-modules)
- [Advanced Usage](#advanced-usage)
- [Configuration](#configuration)
- [Examples](#examples)

---

## Installation

```bash
# Install BOLDGenotyper
pip install boldgenotyper

# Or install from source with all dependencies
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper
pip install -e ".[all]"
```

---

## Quick Start

### Basic Programmatic Usage

```python
import boldgenotyper
from pathlib import Path

# Import specific modules
from boldgenotyper import metadata, dereplication, haplotype_assignment
from boldgenotyper import geographic, visualization, phylogenetics

# Check version
print(f"BOLDGenotyper version: {boldgenotyper.__version__}")

# List available modules
print(f"Available modules: {', '.join(boldgenotyper.__all__)}")
```

### Running Full Pipeline Programmatically

For most users, the CLI commands are recommended. However, you can call pipeline functions directly:

```python
from pathlib import Path
import subprocess

# Using subprocess to call CLI (recommended for full pipeline)
result = subprocess.run([
    'boldgenotyper',
    'data/Sphyrna_lewini.tsv',
    '--clustering-threshold', '0.02',
    '--build-tree',
    '--output', 'results/'
], capture_output=True, text=True)

print(result.stdout)
```

---

## Core Modules

### 1. metadata - Data Loading and Quality Control

**Purpose**: Load BOLD TSV files, validate coordinates, assign ocean basins

**Key Functions**:

```python
from boldgenotyper import metadata

# Parse BOLD TSV file
df = metadata.parse_bold_tsv('data/samples.tsv')
print(f"Loaded {len(df)} samples")

# Filter by coordinate quality
filtered_df = metadata.filter_coordinates(df, exclude_centroids=True)
print(f"Retained {len(filtered_df)} samples with high-quality coordinates")

# Extract sequences
sequences = metadata.extract_sequences(df, min_length=400, max_n_content=0.10)
print(f"Extracted {len(sequences)} valid sequences")
```

**Main Functions**:
- `parse_bold_tsv(tsv_path)`: Load BOLD TSV file
- `filter_coordinates(df, exclude_centroids=True)`: Filter low-quality coordinates
- `extract_sequences(df, min_length=400, max_n_content=0.10)`: Extract and filter sequences
- `validate_required_columns(df)`: Check for required columns

---

### 2. dereplication - Sequence Clustering and Consensus Generation

**Purpose**: Cluster sequences, generate consensus genotypes

**Key Functions**:

```python
from boldgenotyper import dereplication
from pathlib import Path

# Cluster sequences
consensus_fasta, cluster_info = dereplication.dereplicate_sequences(
    sequences=sequences,  # BioPython SeqRecord list
    output_dir=Path('output/dereplication/'),
    clustering_threshold=0.03,  # 97% identity
    threads=8
)

print(f"Generated {len(cluster_info)} consensus genotypes")
print(f"Consensus sequences: {consensus_fasta}")
```

**Main Functions**:
- `dereplicate_sequences(sequences, output_dir, clustering_threshold, threads)`: Main clustering function
- `align_sequences(fasta_path, output_path, threads)`: MAFFT alignment
- `trim_alignment(alignment_path, output_path)`: trimAl trimming
- `calculate_distance_matrix(alignment_path)`: Pairwise distances
- `hierarchical_clustering(distance_matrix, threshold)`: Cluster sequences
- `generate_consensus(alignment, cluster_assignments)`: Create consensus sequences

---

### 3. haplotype_assignment - Sample-to-Haplotype Matching

**Purpose**: Assign samples to haplotype sequences using edit distance

**Key Functions**:

```python
from boldgenotyper import haplotype_assignment

# Assign samples to haplotypes
stats = haplotype_assignment.assign_haplotypes(
    metadata_path='data/samples.tsv',
    fasta_path='output/sequences.fasta',
    consensus_path='output/haplotypes.fasta',
    output_path='output/samples_with_haplotypes.tsv',
    min_identity=0.50,
    tie_margin=0.003,
    tie_min_identity=0.95,
    n_processes=8
)

print(f"Assigned {stats['assigned']}/{stats['total']} samples")
```

**Main Functions**:
- `assign_genotypes(sample_sequences, consensus_fasta, similarity_threshold, tie_margin, tie_threshold, threads)`: Assign all samples
- `calculate_edit_distance(seq1, seq2)`: Compute edit distance between two sequences
- `detect_ties(identities, tie_margin, tie_threshold)`: Flag ambiguous assignments
- `generate_diagnostics(assignments)`: Create diagnostics table

---

### 4. geographic - Ocean Basin Assignment

**Purpose**: Assign samples to ocean basins using GOaS shapefile

**Key Functions**:

```python
from boldgenotyper import geographic

# Assign ocean basins
df_with_basins = geographic.assign_ocean_basins(
    df=df_with_coords,
    goas_shapefile='shapefiles/GOaS_v1_20211214/goas_v01.shp'
)

print(f"Unique basins: {df_with_basins['ocean_basin'].unique()}")
```

**Main Functions**:
- `assign_ocean_basins(df, goas_shapefile)`: Assign samples to basins
- `load_goas_shapefile(shapefile_path)`: Load GOaS data
- `point_to_basin(lat, lon, goas_gdf)`: Get basin for single coordinate
- `filter_valid_coordinates(df)`: Remove invalid coordinates

---

### 5. phylogenetics - Tree Construction

**Purpose**: Build phylogenetic trees with FastTree

**Key Functions**:

```python
from boldgenotyper import phylogenetics

# Build tree
tree_path = phylogenetics.build_tree(
    consensus_fasta='output/consensus.fasta',
    output_dir=Path('output/phylogenetic/'),
    model='GTR+G',  # GTR+Gamma model
    threads=8
)

print(f"Tree saved to: {tree_path}")

# Visualize tree
phylogenetics.plot_tree(
    tree_path=tree_path,
    output_path='output/tree_plot.pdf',
    tip_labels=True
)
```

**Main Functions**:
- `build_tree(consensus_fasta, output_dir, model, threads)`: Construct ML tree with FastTree
- `plot_tree(tree_path, output_path, tip_labels)`: Generate tree visualization
- `relabel_tree(tree_path, label_mapping)`: Replace tip labels with readable names

---

### 6. visualization - Geographic and Distribution Plots

**Purpose**: Generate publication-ready visualizations

**Key Functions**:

```python
from boldgenotyper import visualization

# Geographic distribution map
visualization.plot_distribution_map(
    df=annotated_df,
    output_path='output/distribution_map.pdf',
    genotype_column='consensus_group',
    basin_column='ocean_basin'
)

# Abundance bar chart by basin
visualization.plot_distribution_bars(
    df=annotated_df,
    output_path='output/distribution_bars.pdf',
    genotype_column='consensus_group',
    basin_column='ocean_basin',
    relative=True  # Relative abundance (normalized)
)

# Identity distribution histogram
visualization.plot_identity_distribution(
    df=annotated_df,
    output_path='output/identity_dist.pdf',
    identity_column='identity'
)
```

**Main Functions**:
- `plot_distribution_map(df, output_path, genotype_column, basin_column)`: Global map
- `plot_distribution_bars(df, output_path, genotype_column, basin_column, relative)`: Bar charts
- `plot_identity_distribution(df, output_path, identity_column)`: Identity histogram
- `plot_faceted_maps(df, output_dir, facet_by)`: Multiple maps (one per species/genotype)

---

### 7. compare_analyses.py - Multi-Level Comparison (Standalone Script)

**Purpose**: Compare species vs. family analyses for quality control. This utility
lives in `scripts/compare_analyses.py` and is run directly rather than imported.

**Usage**:

```bash
python scripts/compare_analyses.py \
    --species-level Sphyrna_lewini_output/ \
    --family-level Sphyrnidae_output/ \
    --output comparative_analysis/ \
    --generate-reassignment-table \
    --majority-threshold 0.70
```

**Main Functions** (importable for programmatic use):

```python
from scripts import compare_analyses  # or sys.path-based import

results = compare_analyses.compare_analyses(
    species_path=Path('Sphyrna_lewini_output/'),
    family_path=Path('Sphyrnidae_output/'),
    output_dir=Path('comparison/'),
    generate_reassignment_table=True,
    majority_threshold=0.70
)
# results keys: comparison_summary, genotype_crosswalk, sample_reassignments,
#               methods_text, readme (all Path objects)
```

- `compare_analyses(species_path, family_path, output_dir, generate_reassignment_table, majority_threshold)`: Full comparison workflow
- `load_analysis_results(analysis_path)`: Load annotated CSV from a completed run
- `generate_comparison_summary(species_df, family_df, ...)`: High-level metric table
- `generate_genotype_crosswalk(species_df, family_df, ...)`: Group-level mapping table
- `generate_sample_reassignments(species_df, family_df, ...)`: Sample-level reassignment table
- `generate_methods_text(species_meta, family_meta, summary_df, output_path, majority_threshold)`: Auto-generate methods text

---

### 8. parameter_sweep - Threshold Optimization

**Purpose**: Test multiple thresholds to find optimal value

**Key Functions**:

```python
from boldgenotyper import parameter_sweep

# Run sweep
sweep_results = parameter_sweep.run_parameter_sweep(
    tsv_path=Path('data/samples.tsv'),
    thresholds=[0.01, 0.015, 0.02, 0.03, 0.05],
    output_dir=Path('parameter_sweep/'),
    threads=8,
    keep_intermediates=False
)

print(f"Recommended threshold: {sweep_results['recommended_threshold']}")
print(f"Elbow point: {sweep_results['elbow_point']}")
```

**Main Functions**:
- `run_parameter_sweep(tsv_path, thresholds, output_dir, threads, keep_intermediates)`: Run full sweep
- `calculate_sweep_metrics(threshold_results)`: Compile metrics across thresholds
- `detect_elbow_point(metrics_df)`: Identify optimal threshold
- `plot_stability(metrics_df, output_path)`: Generate stability plots

---

### 9. popgen_export - Population Genetics Formats

**Purpose**: Export genotypes to external software formats

**Key Functions**:

```python
from boldgenotyper import popgen_export

# Export to multiple formats
export_results = popgen_export.export_population_genetics_formats(
    df=annotated_df,
    consensus_fasta_path=Path('output/consensus.fasta'),
    output_dir=Path('exports/'),
    organism='Sphyrna_lewini',
    formats=['arlequin', 'popart', 'dnasp', 'generic']
)

print(f"Exported formats: {list(export_results.keys())}")
```

**Main Functions**:
- `export_population_genetics_formats(df, consensus_fasta_path, output_dir, organism, formats)`: Export all formats
- `export_arlequin(df, consensus_fasta, output_path)`: Arlequin .arp file
- `export_popart(df, consensus_fasta, output_path)`: PopART NEXUS file
- `export_dnasp(df, consensus_fasta, output_path)`: DnaSP FASTA file
- `export_generic(df, consensus_fasta, output_dir)`: Generic CSV/FASTA files

---

## Configuration

### Using Configuration Objects

```python
from boldgenotyper import config

# Load default configuration
cfg = config.get_default_config()

# Access settings
print(f"Default clustering threshold: {cfg.dereplication.clustering_threshold}")
print(f"Default similarity threshold: {cfg.genotype_assignment.min_identity}")
print(f"Default threads: {cfg.n_threads}")

# Modify settings
cfg_custom = cfg.update(
    dereplication__clustering_threshold=0.02,
    genotype_assignment__min_identity=0.85,
    n_threads=16
)

print(f"Updated clustering threshold: {cfg_custom.dereplication.clustering_threshold}")
```

### Configuration Parameters

**Dereplication**:
- `clustering_threshold`: Maximum distance for clustering (default: 0.03)
- `min_sequence_length`: Minimum sequence length (default: 400)
- `max_n_content`: Maximum N content proportion (default: 0.10)
- `consensus_frequency_cutoff`: Consensus base frequency (default: 0.70)

**Genotype Assignment**:
- `min_identity`: Minimum similarity for assignment (default: 0.50)
- `tie_margin`: Identity difference to call tie (default: 0.003)
- `tie_threshold`: Minimum identity for tie detection (default: 0.95)

**Geographic**:
- `goas_shapefile_path`: Path to GOaS shapefile
- `exclude_centroids`: Exclude country centroids (default: True)

**General**:
- `n_threads`: Number of parallel threads (default: 4)

---

## Examples

### Example 1: Custom Haplotype Pipeline

```python
from boldgenotyper import metadata, dereplication, haplotype_assignment
from pathlib import Path

# 1. Load data
df = metadata.parse_bold_tsv('data/samples.tsv')

# 2. Identify haplotypes using ESV approach
haplotypes = dereplication.identify_haplotypes(
    tsv_path='data/samples.tsv',
    output_dir=Path('custom_output/'),
    n_threads=16
)

# 3. Assign samples to haplotypes
stats = haplotype_assignment.assign_haplotypes(
    metadata_path='data/samples.tsv',
    fasta_path='custom_output/sequences.fasta',
    consensus_path='custom_output/haplotypes.fasta',
    output_path='custom_output/samples_with_haplotypes.tsv',
    min_identity=0.80,
    n_processes=16
)

# 4. Load results
import pandas as pd
final_df = pd.read_csv('custom_output/samples_with_haplotypes.tsv', sep='\t')

# 5. Save results
final_df.to_csv('custom_output/annotated.csv', index=False)
print(f"Saved {len(final_df)} annotated samples")
```

### Example 2: Geographic Analysis with Custom Shapefile

```python
from boldgenotyper import metadata, geographic
import geopandas as gpd
from pathlib import Path

# Load data
df = metadata.parse_bold_tsv('data/Salmonidae.tsv')

# Load custom shapefile (freshwater basins)
basins_gdf = gpd.read_file('data/hydrobasins_lev04.shp')

# Assign samples to basins
def assign_to_basin(row, gdf):
    try:
        from shapely.geometry import Point
        point = Point(row['lon'], row['lat'])
        matches = gdf[gdf.contains(point)]
        if len(matches) > 0:
            return matches.iloc[0]['MAIN_BAS']
        return 'Unknown'
    except:
        return 'Unknown'

df['freshwater_basin'] = df.apply(lambda row: assign_to_basin(row, basins_gdf), axis=1)

# Summarize by basin
basin_summary = df.groupby('freshwater_basin').agg({
    'processid': 'count',
    'species': lambda x: x.nunique()
}).rename(columns={'processid': 'n_samples', 'species': 'n_species'})

print(basin_summary)
```

### Example 3: Batch Processing Multiple Datasets

```python
from pathlib import Path
import subprocess

# List of datasets to process
datasets = [
    'data/Sphyrna_lewini.tsv',
    'data/Sphyrna_tiburo.tsv',
    'data/Sphyrna_zygaena.tsv'
]

# Process each dataset
for dataset_path in datasets:
    dataset = Path(dataset_path)
    organism = dataset.stem
    output_dir = Path(f'batch_results/{organism}/')

    print(f"Processing {organism}...")

    result = subprocess.run([
        'boldgenotyper',
        str(dataset),
        '--clustering-threshold', '0.02',
        '--build-tree',
        '--export-format', 'all',
        '--output', str(output_dir),
        '--threads', '8'
    ], capture_output=True, text=True)

    if result.returncode == 0:
        print(f"  ✓ {organism} complete")
    else:
        print(f"  ✗ {organism} failed: {result.stderr}")
```

### Example 4: Custom Visualization

```python
from boldgenotyper import visualization
import pandas as pd
import matplotlib.pyplot as plt

# Load annotated data
df = pd.read_csv('output/annotated.csv')

# Filter to assigned samples only
df_assigned = df[df['assignment_status'] == 'assigned'].copy()

# Create custom plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Identity distribution
ax = axes[0, 0]
df_assigned['identity'].hist(bins=50, ax=ax, edgecolor='black')
ax.set_xlabel('Sequence Identity')
ax.set_ylabel('Frequency')
ax.set_title('Identity Distribution')

# Genotypes by ocean basin
ax = axes[0, 1]
basin_counts = df_assigned.groupby(['ocean_basin', 'consensus_group']).size().unstack(fill_value=0)
basin_counts.plot(kind='bar', stacked=True, ax=ax, legend=False)
ax.set_xlabel('Ocean Basin')
ax.set_ylabel('Sample Count')
ax.set_title('Genotypes by Basin')

# Sample counts by genotype
ax = axes[1, 0]
genotype_counts = df_assigned['consensus_group'].value_counts()
genotype_counts.plot(kind='barh', ax=ax)
ax.set_xlabel('Sample Count')
ax.set_ylabel('Genotype')
ax.set_title('Sample Distribution')

# Geographic scatter
ax = axes[1, 1]
basins = df_assigned['ocean_basin'].unique()
for basin in basins:
    basin_df = df_assigned[df_assigned['ocean_basin'] == basin]
    ax.scatter(basin_df['lon'], basin_df['lat'], label=basin, alpha=0.6)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Geographic Distribution')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig('custom_visualization.pdf', dpi=300, bbox_inches='tight')
print("Custom visualization saved")
```

### Example 5: Integration with External Tools

```python
from boldgenotyper import popgen_export
import pandas as pd
import subprocess

# 1. Run BOLDGenotyper pipeline
subprocess.run([
    'boldgenotyper',
    'data/samples.tsv',
    '--clustering-threshold', '0.02',
    '--output', 'boldgen_output/'
])

# 2. Load results
df = pd.read_csv('boldgen_output/annotated.csv')

# 3. Export for PopART (haplotype network)
export_results = popgen_export.export_population_genetics_formats(
    df=df,
    consensus_fasta_path='boldgen_output/consensus.fasta',
    output_dir='popgen_exports/',
    organism='samples',
    formats=['popart']
)

print(f"PopART file: {export_results['popart']['nexus']}")

# 4. (Optional) Open in PopART
# Manual step: Open the .nexus file in PopART software

# 5. Export for Arlequin (population statistics)
export_results = popgen_export.export_population_genetics_formats(
    df=df,
    consensus_fasta_path='boldgen_output/consensus.fasta',
    output_dir='popgen_exports/',
    organism='samples',
    formats=['arlequin']
)

print(f"Arlequin file: {export_results['arlequin']['arp']}")
```

---

## Error Handling

### Common Errors and Solutions

```python
from boldgenotyper import metadata

try:
    df = metadata.parse_bold_tsv('data/samples.tsv')
except FileNotFoundError:
    print("Error: TSV file not found. Check the path.")
except ValueError as e:
    print(f"Error parsing TSV: {e}")
    print("Check that required columns are present: processid, nucleotides, species")

try:
    sequences = metadata.extract_sequences(df)
except KeyError as e:
    print(f"Missing required column: {e}")
except Exception as e:
    print(f"Error extracting sequences: {e}")
```

### Logging

```python
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('boldgenotyper')

# Now BOLDGenotyper will output detailed progress
from boldgenotyper import dereplication
# ... dereplication functions will log progress
```

---

## Additional Resources

### Full Documentation

- **README.md**: Main documentation and CLI usage
- **CUSTOM_SHAPEFILES_GUIDE.md**: Custom geographic analysis
- **PARAMETER_SWEEP_GUIDE.md**: Threshold optimization
- **COMPARATIVE_ANALYSIS_GUIDE.md**: Quality control workflow
- **POPGEN_EXPORT_GUIDE.md**: Population genetics exports

### Source Code

Explore module source code for advanced usage:
- GitHub: https://github.com/SymbioSeas/BOLDGenotyper
- Module source: `boldgenotyper/` directory

### Support

- **Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
- **Email**: symbioseas@outlook.com

---

**Last Updated**: 2025-11-25
**BOLDGenotyper Version**: 0.1.0
