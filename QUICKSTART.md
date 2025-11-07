# BOLDGenotyper Quick Start Guide

## Installation (Already Done!)

The package is installed and the `boldgenotyper` command is ready to use.

## Testing with Your Own Dataset

### Step 1: Prepare Your Data

You need a BOLD TSV file containing:
- Sample metadata
- COI sequences in the `nuc` column
- Geographic coordinates (optional, in `coord` column)
- Species identifications in `species` column

### Step 2: Run the Pipeline

**Basic command:**
```bash
boldgenotyper /path/to/your/data.tsv
```

**Example with your Euprymna data (if you have it):**
```bash
boldgenotyper tests/data/Euprymna_scolopes.tsv
```

**Example with Carcharhinus (lower threshold for high diversity):**
```bash
boldgenotyper tests/data/Carcharhinus.tsv --similarity-threshold 0.80 --threads 8
```

### Step 3: Check the Results

After the pipeline completes, you'll find:

```
Euprymna_output/                        # or YourOrganism_output/
├── Euprymna_annotated.csv              # ← START HERE! Main results
├── Euprymna_pipeline.log               # Complete log
├── consensus_sequences/                 # Genotype sequences
│   └── Euprymna_consensus.fasta
├── genotype_assignments/               # Assignment details
│   ├── Euprymna_with_genotypes.tsv
│   └── Euprymna_diagnostics.csv
├── taxonomy/                           # Species per genotype
│   └── Euprymna_consensus_taxonomy.csv
├── visualization/                      # Plots
│   ├── Euprymna_distribution_map.png
│   ├── Euprymna_identity_distribution.png
│   └── Euprymna_tree.png              # (if --build-tree used)
└── reports/                            # Summary statistics
    └── Euprymna_pipeline_summary.txt
```

## Common Use Cases

### 1. Quick analysis with defaults
```bash
boldgenotyper mydata.tsv
```

### 2. Highly diverse genus (lower similarity threshold)
```bash
boldgenotyper mydata.tsv --similarity-threshold 0.80
```

### 3. Custom output location
```bash
boldgenotyper mydata.tsv --output results/analysis_2025
```

### 4. With phylogenetic tree
```bash
boldgenotyper mydata.tsv --build-tree
```
*Requires: MAFFT and FastTree installed*

### 5. Full custom run
```bash
boldgenotyper mydata.tsv \
    --organism "Genus species" \
    --output results/my_analysis \
    --similarity-threshold 0.85 \
    --threads 8 \
    --build-tree \
    --log-level DEBUG
```

## What to Look For in Results

### 1. Main Results File: `{Organism}_annotated.csv`
Open in Excel or your favorite data viewer. Key columns:
- `consensus_group_sp`: Genotype assignment with species name
- `assigned_sp`: Species assigned to the genotype
- `lat`, `lon`, `ocean_basin`: Geographic data
- `best_identity`: How well the sample matched its genotype

### 2. Assignment Success Rate
Check the log file or terminal output for:
```
✓ Assigned XXX/YYYY samples to genotypes
✓ Assignment rate: XX.X%
```

Low assignment rates (<30%) suggest:
- High genetic diversity (try lowering `--similarity-threshold`)
- Poor sequence quality
- Multiple cryptic species

### 3. Visualizations
- **Identity distribution**: Shows assignment confidence
- **Distribution map**: Geographic spread of genotypes
- **Phylogenetic tree**: Evolutionary relationships

## Troubleshooting

### "No such file or directory"
Make sure the TSV path is correct:
```bash
ls tests/data/  # Check your data directory
```

### "Pipeline failed"
1. Check the log file: `less {Organism}_output/{Organism}_pipeline.log`
2. Try with `--log-level DEBUG` for more details
3. Verify your TSV has required columns: `processid`, `nuc`, `species`

### Low assignment rates
This is often **biologically expected** for diverse genera. Try:
- Lowering threshold: `--similarity-threshold 0.80` or `0.75`
- Check sequence quality in diagnostics file
- Review consensus taxonomy to see if genotypes make sense

### Missing visualizations
Some plots require:
- Geographic coordinates for distribution maps
- MAFFT/FastTree for phylogenetic trees
- Geopandas/Cartopy for map plotting

The pipeline will continue even if these fail.

## Performance Tips

- Use `--threads 8` (or more) for faster processing
- Small datasets (<500 samples): ~2-5 minutes
- Large datasets (>2000 samples): ~15-60 minutes
- Tree building adds 5-20 minutes depending on size

## Next Steps

1. **Explore `{Organism}_annotated.csv`** - This has everything!
2. **Check `visualization/` folder** - Look at the plots
3. **Review `taxonomy/` folder** - See species assignments per genotype
4. **Read `reports/` folder** - Pipeline summary statistics

## Getting Help

- View all options: `boldgenotyper --help`
- Full documentation: See `TESTING_CLI.md`
- Report issues: Check log file first, then contact developer

## Example Session

```bash
# Navigate to your data directory
cd /path/to/your/bold/data

# Run the pipeline (organism name inferred from filename)
boldgenotyper Euprymna_scolopes.tsv --threads 8

# Pipeline runs... (5-10 minutes)

# Check results
ls Euprymna_output/
open Euprymna_output/Euprymna_annotated.csv           # Main results
open Euprymna_output/visualization/                    # View plots
less Euprymna_output/Euprymna_pipeline.log            # Check log

# Success! 🎉
```

## Tips for Your Own Data

1. **File naming**: Name your TSV like `Genus_species.tsv` for auto-detection
2. **Data quality**: More samples = better genotype resolution
3. **Thresholds**: Start with default (0.90), adjust if needed
4. **Rerun easily**: Just change the output directory
5. **Compare runs**: Use different output directories for different parameters

Happy genotyping! 🧬
