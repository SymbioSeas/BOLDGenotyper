# Population Genetics Export Formats - User Guide

## Overview

BOLDGenotyper now supports exporting genotyping results to multiple population genetics software formats, enabling seamless integration with downstream analysis tools.

## Quick Start

### Basic Usage

Export to all formats:
```bash
boldgenotyper Sphyrnidae.tsv --export-format all
```

Export to specific formats:
```bash
boldgenotyper Sphyrnidae.tsv --export-format arlequin popart dnasp
```

### Available Formats

1. **Arlequin** (`arlequin`) - For Arlequin software
2. **PopART/NEXUS** (`popart`) - For PopART haplotype network analysis
3. **DnaSP** (`dnasp`) - For DNA sequence polymorphism analysis
4. **Generic** (`generic`) - CSV and FASTA formats for general use
5. **All** (`all`) - Export all formats at once

## Output Structure

When you run the pipeline with `--export-format`, exports are saved to:

```
{organism}_output/
└── exports/
    ├── README.md                    # Format descriptions and usage instructions
    ├── arlequin/
    │   ├── {organism}.arp          # Arlequin project file
    │   └── populations.txt         # Population definitions
    ├── popart/
    │   ├── {organism}.nexus        # NEXUS alignment with traits
    │   └── populations.csv         # Population mapping
    ├── dnasp/
    │   └── {organism}.fas          # FASTA with population labels
    └── generic/
        ├── alignment.fasta         # Consensus sequences
        ├── genotype_membership.csv # Sample-to-genotype mapping
        └── haplotypes.csv          # Haplotype summary table
```

## Format Details

### 1. Arlequin Format (.arp)

**Software:** Arlequin (http://cmpg.unibe.ch/software/arlequin35/)
**Use for:** Population genetics analysis, AMOVA, genetic diversity, population structure

**Files:**
- `{organism}.arp` - Main project file
- `populations.txt` - Population metadata

**How to use:**
1. Open Arlequin
2. File → Open → Select the .arp file
3. Configure your analysis
4. Execute → Run computations

### 2. PopART/NEXUS Format

**Software:** PopART (http://popart.otago.ac.nz/)
**Use for:** Haplotype network visualization, phylogeography, TCS networks

**Files:**
- `{organism}.nexus` - NEXUS alignment with trait data
- `populations.csv` - Population assignments

**How to use:**
1. Open PopART
2. File → Open → Select the .nexus file
3. Choose network method (TCS, Median-Joining, etc.)
4. Compute network

### 3. DnaSP Format (.fas)

**Software:** DnaSP (http://www.ub.edu/dnasp/)
**Use for:** Nucleotide diversity, neutrality tests, haplotype analysis

**Files:**
- `{organism}.fas` - FASTA with metadata in headers

**Header format:**
```
>sample_id [genotype=c1_n10;species=Sphyrna_lewini;basin=Atlantic;processid=SAMPLE001]
```

**How to use:**
1. Open DnaSP
2. File → Open Data File → Select the .fas file
3. Parse population assignments from headers
4. Run analyses

### 4. Generic Formats

**Use for:** Custom analyses, R/Python scripts, any software accepting standard formats

**Files:**

1. **alignment.fasta** - Consensus sequences for all genotypes
   ```
   >c1_n10 n=10 species=Sphyrna_lewini
   ATCGATCGATCG...
   ```

2. **genotype_membership.csv** - Complete sample metadata
   ```csv
   processid,consensus_group,species,lat,lon,ocean_basin,consensus_sequence,sequence_length
   SAMPLE001,c1_n10,Sphyrna lewini,-25.3,35.2,Indian Ocean,ATCG...,650
   ```

3. **haplotypes.csv** - Haplotype summary statistics
   ```csv
   haplotype_id,sequence,sequence_length,n_samples,species,primary_basin,n_basins
   c1_n10,ATCG...,650,386,Sphyrna lewini,Indian Ocean,3
   ```

## Example Workflows

### Workflow 1: Haplotype Network in PopART

```bash
# Run pipeline with PopART export
boldgenotyper Sphyrnidae.tsv --export-format popart

# Output: Sphyrnidae_output/exports/popart/Sphyrnidae.nexus
```

Then in PopART:
1. File → Open → Select `Sphyrnidae.nexus`
2. Network → TCS Network
3. Color by trait: `ocean_basin`
4. Generate network

### Workflow 2: Population Structure in Arlequin

```bash
# Run pipeline with Arlequin export
boldgenotyper Sphyrnidae.tsv --export-format arlequin

# Output: Sphyrnidae_output/exports/arlequin/Sphyrnidae.arp
```

Then in Arlequin:
1. File → Open → Select `Sphyrnidae.arp`
2. Settings → AMOVA
3. Define population structure
4. Execute → AMOVA

### Workflow 3: Custom Analysis in R

```bash
# Export generic formats
boldgenotyper Sphyrnidae.tsv --export-format generic
```

Then in R:
```r
# Load genotype membership
library(tidyverse)
membership <- read_csv("Sphyrnidae_output/exports/generic/genotype_membership.csv")

# Load haplotypes
haplotypes <- read_csv("Sphyrnidae_output/exports/generic/haplotypes.csv")

# Analyze geographic distribution
membership %>%
  group_by(consensus_group, ocean_basin) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = ocean_basin, y = n, fill = consensus_group)) +
  geom_col()
```

### Workflow 4: Export Multiple Formats

```bash
# Export all formats for comprehensive analysis
boldgenotyper Sphyrnidae.tsv \
  --export-format all \
  --clustering-threshold 0.015 \
  --threads 8
```

## Integration with Full Pipeline

The export formats integrate seamlessly with the full pipeline:

```bash
boldgenotyper Sphyrnidae.tsv \
  --clustering-threshold 0.015 \
  --similarity-threshold 0.50 \
  --build-tree \
  --export-format all \
  --export-plot-data \
  --threads 8
```

This will:
1. Run complete genotyping pipeline
2. Generate all visualizations
3. Export population genetics formats
4. Export plot data for custom figures
5. Generate HTML report

## File Naming Conventions

All exported files use the organism name from your input file:

- Input: `Sphyrna_lewini.tsv` → Organism: `Sphyrna_lewini`
- Exports: `Sphyrna_lewini.arp`, `Sphyrna_lewini.nexus`, etc.

## Data Included in Exports

### Genotype Information
- Consensus sequences for each genotype
- Sample assignments to genotypes
- Genotype labels (e.g., `c15_n386`)

### Geographic Data (if available)
- Ocean basin assignments
- Latitude/longitude coordinates
- Country information

### Taxonomic Data
- Species assignments
- Sample identifiers (processid)

### Metadata
- Sequence lengths
- Sample counts per genotype
- Population structure

## Tips and Best Practices

1. **Format Selection**: Choose formats based on your analysis goals
   - Phylogeography → PopART
   - Population structure → Arlequin
   - Diversity metrics → DnaSP
   - Custom analyses → Generic

2. **Data Quality**: Exports only include samples successfully assigned to genotypes
   - Unassigned samples are excluded
   - Check assignment diagnostics before exporting

3. **Downstream Analysis**: Always check the generated README.md for format-specific details

4. **Large Datasets**: Consider exporting only needed formats to save time
   ```bash
   # Only export what you need
   boldgenotyper large_dataset.tsv --export-format generic
   ```

5. **Combining with Other Tools**: Generic formats are ideal for custom pipelines
   ```bash
   # Export generic, then process with custom scripts
   boldgenotyper data.tsv --export-format generic
   python my_analysis.py exports/generic/genotype_membership.csv
   ```

## Troubleshooting

### Issue: No sequences in export

**Cause**: No samples were assigned to genotypes
**Solution**: Check genotype assignment diagnostics, adjust similarity threshold

### Issue: Missing geographic data

**Cause**: GOaS shapefile not available
**Solution**: Download manually from https://www.marineregions.org/downloads.php (search for `GOaS_v1_20211214.zip`, registration form required), extract to `shapefiles/GOaS_v1_20211214/`, then verify with `python -m boldgenotyper.goas_downloader`. Alternatively, use `--no-geo` to skip geographic analysis.

### Issue: Format-specific software won't open file

**Cause**: Incorrect file format or path
**Solution**: Check README.md for format specifications, verify software version

## Advanced Usage

### Custom Population Grouping

By default, exports group samples by consensus genotype. For custom grouping:

```python
from boldgenotyper import popgen_export
import pandas as pd

# Load annotated data
df = pd.read_csv("Sphyrnidae_output/Sphyrnidae_annotated.csv")

# Export grouped by ocean basin instead
results = popgen_export.export_population_genetics_formats(
    df=df,
    consensus_fasta_path="consensus.fasta",
    output_dir="custom_exports",
    organism="Sphyrnidae",
    formats=['arlequin'],
    group_by='ocean_basin'  # Custom grouping
)
```

## Citation

If you use these exports in a publication, please cite:

1. **BOLDGenotyper:**
   Smith et al. (in prep) BOLDGenotyper: Automated genotyping and quality control for DNA barcode data.

2. **BOLD Database:**
   Ratnasingham S, Hebert PDN (2007) BOLD: The Barcode of Life Data System.
   Molecular Ecology Notes 7: 355-364.

3. **Software-specific citations:**
   - Arlequin: Excoffier & Lischer (2010)
   - PopART: Leigh & Bryant (2015)
   - DnaSP: Rozas et al. (2017)

## Support

For questions or issues:
- GitHub: https://github.com/yourusername/boldgenotyper
- Documentation: [BOLDGenotyper docs]
- Email: symbioseas@outlook.com

---

**Version:** 1.0.0
**Last Updated:** 2025-11-25
**Module:** `boldgenotyper.popgen_export`
