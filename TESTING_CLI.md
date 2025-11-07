# Testing the BOLDGenotyper CLI

## Installation

### Step 1: Install the package in development mode

From the `boldgenotyper` directory, run:

```bash
pip install -e .
```

This installs the package and registers the `boldgenotyper` command.

### Step 2: Verify installation

Check that the command is available:

```bash
boldgenotyper --version
```

You should see: `BOLDGenotyper 1.0.0`

### Step 3: View help

```bash
boldgenotyper --help
```

## Testing with Your Own Data

### Basic Usage

The simplest way to run the pipeline:

```bash
boldgenotyper /path/to/your/data.tsv
```

This will:
- Automatically extract the organism name from the filename
- Create an output directory named `{Organism}_output` in the current directory
- Use default parameters (90% similarity threshold, 4 threads, no tree building)

### Examples

#### 1. Basic run with Euprymna data

```bash
boldgenotyper tests/data/Euprymna_scolopes.tsv
```

Expected output directory: `Euprymna_output/`

#### 2. Specify organism name explicitly

If your TSV filename doesn't contain the organism name:

```bash
boldgenotyper tests/data/samples.tsv --organism Euprymna
```

#### 3. Specify custom output directory

```bash
boldgenotyper tests/data/Carcharhinus.tsv --output results/shark_analysis
```

#### 4. Adjust similarity threshold for highly diverse taxa

For genera with high genetic diversity (like Carcharhinus), lower the threshold:

```bash
boldgenotyper tests/data/Carcharhinus.tsv --similarity-threshold 0.80
```

#### 5. Enable phylogenetic tree building

Requires MAFFT and FastTree to be installed:

```bash
boldgenotyper tests/data/Euprymna.tsv --build-tree
```

#### 6. Use more threads for faster processing

```bash
boldgenotyper tests/data/Carcharhinus.tsv --threads 8
```

#### 7. Increase logging verbosity for debugging

```bash
boldgenotyper tests/data/Euprymna.tsv --log-level DEBUG
```

#### 8. Complete example with all options

```bash
boldgenotyper tests/data/Carcharhinus.tsv \
    --organism Carcharhinus \
    --output results/Carcharhinus_analysis_v2 \
    --similarity-threshold 0.80 \
    --threads 8 \
    --build-tree \
    --log-level INFO
```

## Output Directory Structure

After running the pipeline, you'll find an organized directory structure:

```
{Organism}_output/
├── {Organism}_annotated.csv              # Main results file
├── {Organism}_pipeline.log               # Complete log file
├── intermediate/                         # Temporary processing files
│   ├── 01_parsed_metadata.tsv
│   ├── 02_filtered_metadata.tsv
│   ├── {Organism}.fasta
│   └── dereplication/
│       ├── {Organism}_aligned.fasta
│       ├── {Organism}_trimmed.fasta
│       └── {Organism}_consensus.fasta
├── consensus_sequences/                  # Consensus genotype sequences
│   └── {Organism}_consensus.fasta
├── genotype_assignments/                 # Assignment results
│   ├── {Organism}_with_genotypes.tsv
│   └── {Organism}_diagnostics.csv
├── taxonomy/                            # Taxonomy assignment results
│   ├── {Organism}_consensus_taxonomy.csv
│   └── {Organism}_species_by_consensus.csv
├── geographic/                          # Geographic analysis
│   └── samples_with_ocean_basins.tsv
├── phylogenetic/                        # Phylogenetic trees (if --build-tree)
│   ├── {Organism}_tree.nwk
│   ├── {Organism}_aligned.fasta
│   └── {Organism}_trimmed.fasta
├── visualization/                       # Plots and figures
│   ├── {Organism}_distribution_map.png
│   ├── {Organism}_distribution_map.pdf
│   ├── {Organism}_identity_distribution.png
│   ├── {Organism}_tree.png              # (if --build-tree)
│   └── {Organism}_tree.pdf
└── reports/                            # Summary reports
    └── {Organism}_pipeline_summary.txt
```

## Key Output Files

### 1. Main Results
- **`{Organism}_annotated.csv`**: Complete dataset with genotype assignments, taxonomy, and geography

### 2. Consensus Sequences
- **`consensus_sequences/{Organism}_consensus.fasta`**: Unique genotype sequences

### 3. Assignments
- **`genotype_assignments/{Organism}_with_genotypes.tsv`**: Sample-to-genotype mapping
- **`genotype_assignments/{Organism}_diagnostics.csv`**: Assignment quality metrics

### 4. Taxonomy
- **`taxonomy/{Organism}_consensus_taxonomy.csv`**: Species assignments per genotype
- **`taxonomy/{Organism}_species_by_consensus.csv`**: Species composition within each genotype

### 5. Visualizations
- **`visualization/{Organism}_distribution_map.png/pdf`**: Geographic distribution of genotypes
- **`visualization/{Organism}_identity_distribution.png/pdf`**: Assignment quality histogram
- **`visualization/{Organism}_tree.png/pdf`**: Phylogenetic tree (if enabled)

## Troubleshooting

### Command not found

If you get `boldgenotyper: command not found`:

1. Make sure you installed with `pip install -e .`
2. Check your Python environment is activated
3. Try reinstalling: `pip uninstall boldgenotyper && pip install -e .`

### Missing dependencies

If you get import errors:

```bash
pip install -r requirements.txt
```

Or install dependencies manually:

```bash
pip install biopython pandas scipy numpy matplotlib seaborn
```

### External tools not found

For phylogenetic analysis, you need:

- **MAFFT**: `conda install -c bioconda mafft` or `brew install mafft`
- **FastTree**: `conda install -c bioconda fasttree` or `brew install fasttree`
- **trimAl**: `conda install -c bioconda trimal` or `brew install trimal`

For geographic analysis:

- **Geopandas**: `conda install -c conda-forge geopandas`
- **Cartopy**: `conda install -c conda-forge cartopy`

### Pipeline fails partway through

Check the log file in the output directory:

```bash
less {Organism}_output/{Organism}_pipeline.log
```

The pipeline is designed to continue even if non-critical steps fail (e.g., visualization, tree building).

## Comparing with Test Scripts

The CLI integrates functionality from:

1. **test_pipeline_script.py** → Phases 1-4 (data loading through taxonomy)
2. **test_phylo_viz_pipeline.py** → Phase 5 (phylogenetics and visualization)
3. **test_full_pipeline.py** → Phase 6-7 (additional visualizations and reports)

You can verify CLI output matches the test scripts by comparing key files.

## Performance Notes

- **Small datasets (<500 samples)**: ~2-5 minutes
- **Medium datasets (500-2,000 samples)**: ~5-15 minutes
- **Large datasets (>2,000 samples)**: ~15-60 minutes

Most time is spent on:
1. Sequence alignment (MAFFT)
2. Genotype assignment (pairwise comparisons)
3. Phylogenetic tree building (if enabled)

Use `--threads` to speed up parallel steps.

## Next Steps

After running the pipeline:

1. **Review the annotated CSV**: Main results file with all data integrated
2. **Check assignment rates**: Look at diagnostics and identity distribution plots
3. **Examine taxonomy**: Review species assignments in the taxonomy directory
4. **Visualize results**: Check distribution maps and phylogenetic trees
5. **Read the summary**: Pipeline summary report has key statistics

## Advanced Usage

For more control over parameters, you can:

1. Create a custom configuration YAML file
2. Modify source code in `boldgenotyper/cli.py`
3. Use the individual modules programmatically in Python scripts

See the full documentation for details.
