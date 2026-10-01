# Population Genetics Export Formats - User Guide

## Overview

BOLDGenotyper exports haplotype assignments to population genetics software formats
with `--export-format`. Nothing is exported unless you ask for it.

Populations are the regions of the run's region column (`--geo-category`; `ocean_basin`
with GOaS, or the field of a custom shapefile). Because haplotypes are exact sequence
variants of the core region, every record of a haplotype has the same core sequence, so
the haplotype's sequence is used for each of its records.

## Quick Start

Export to all formats:
```bash
boldgenotyper Sphyrnidae.tsv --singleton-distance 0.015 --export-format all
```

Export to specific formats:
```bash
boldgenotyper Sphyrnidae.tsv --singleton-distance 0.015 --export-format arlequin popart dnasp
```

### Available Formats

1. **Arlequin** (`arlequin`) - For Arlequin software
2. **PopART/NEXUS** (`popart`) - For PopART haplotype networks
3. **DnaSP** (`dnasp`) - For DNA sequence polymorphism analysis
4. **Generic** (`generic`) - CSV tables for general use
5. **All** (`all`) - All of the above

## Which Records Are Exported

Exports include the records assigned to a haplotype (status `assigned` or
`assigned_corrected` in `records.csv`). A record has a region only when its coordinates
are of geographic quality (`coord_status` ok); records with missing, centroid or zero
coordinates, or outside every region, have no region.

| Format | Records without a region |
|---|---|
| Arlequin | left out (every sample needs a population) |
| DnaSP | left out (populations are read from the sequence names) |
| PopART | kept, with the region trait `Unassigned` |
| Generic | kept, with an empty region |

The log reports how many records each format left out. If many records have no region,
check `coord_status` in `records.csv` before interpreting population comparisons.

## Output Structure

```
{organism}_output/
└── exports/
    ├── README.md                    # format descriptions and usage
    ├── arlequin/
    │   ├── {organism}.arp           # Arlequin project file
    │   └── populations.csv          # populations: records, haplotypes, most common species
    ├── popart/
    │   ├── {organism}.nexus         # NEXUS alignment with traits
    │   └── popart_traits.csv        # traits per taxon
    ├── dnasp/
    │   └── {organism}.fas           # FASTA named by population
    └── generic/
        ├── haplotype_membership.csv # record-to-haplotype table
        └── haplotype_summary.csv    # one row per haplotype
```

The haplotype sequences themselves are in the run's `haplotypes.fasta`. The run's
`README.md` defines every column of these files.

## Format Details

### 1. Arlequin Format (.arp)

**Software:** Arlequin (http://cmpg.unibe.ch/software/arlequin35/)
**Use for:** AMOVA, genetic diversity, population structure

**Files:**
- `{organism}.arp` - project file, one sample group per region
- `populations.csv` - `population`, `n_samples`, `n_haplotypes`, `dominant_species`

**How to use:**
1. Open Arlequin
2. File → Open → Select the .arp file
3. Configure your analysis
4. Execute → Run computations

### 2. PopART/NEXUS Format

**Software:** PopART (http://popart.otago.ac.nz/)
**Use for:** Haplotype networks, phylogeography

**Files:**
- `{organism}.nexus` - alignment with one taxon per record and trait data
- `popart_traits.csv` - `taxon` (haplotype ID and a running number, e.g. `H0001_001`),
  `haplotype`, `species`, and the region (`Unassigned` when none)

**How to use:**
1. Open PopART
2. File → Open → Select the .nexus file
3. Choose network method (TCS, Median-Joining, etc.)
4. Color by the region trait

### 3. DnaSP Format (.fas)

**Software:** DnaSP (http://www.ub.edu/dnasp/)
**Use for:** Nucleotide diversity, neutrality tests, haplotype analysis

**Files:**
- `{organism}.fas` - one sequence per record with a region

**Header format** (`{population}|{processid}`, then metadata):
```
>Indian_Ocean|EMADA222-09 [haplotype=H0001;species=Sphyrna_lewini;ocean_basin=Indian_Ocean;processid=EMADA222-09]
```

**How to use:**
1. Open DnaSP
2. File → Open Data File → Select the .fas file
3. Define Sequence Sets → Auto-detect by name (reads the population before `|`)
4. Run analyses

### 4. Generic Formats

**Use for:** R/Python scripts and other software

1. **haplotype_membership.csv** - one row per assigned record
   ```csv
   processid,haplotype_id,species,lat,lon,ocean_basin,country/ocean,consensus_sequence,sequence_length
   SAMPLE001,H0001,Sphyrna lewini,-25.3,35.2,Indian Ocean,Mozambique,ATCG...,650
   ```

2. **haplotype_summary.csv** - one row per haplotype
   ```csv
   haplotype_id,sequence,sequence_length,n_samples,species,primary_population,n_populations
   H0001,ATCG...,650,386,Sphyrna lewini,Indian Ocean,3
   ```

## Example Workflows

### Haplotype Network in PopART

```bash
boldgenotyper Sphyrnidae.tsv --singleton-distance 0.015 --export-format popart
# Output: Sphyrnidae_output/exports/popart/Sphyrnidae.nexus
```

Then in PopART: File → Open → `Sphyrnidae.nexus`; Network → TCS Network; color by
`ocean_basin`.

### Population Structure in Arlequin

```bash
boldgenotyper Sphyrnidae.tsv --singleton-distance 0.015 --export-format arlequin
# Output: Sphyrnidae_output/exports/arlequin/Sphyrnidae.arp
```

Then in Arlequin: File → Open → `Sphyrnidae.arp`; Settings → AMOVA; Execute.

### Custom Analysis in R

The run's own tables need no export:

```r
library(tidyverse)
records <- read_csv("Sphyrnidae_output/records.csv")

records %>%
  filter(status %in% c("assigned", "assigned_corrected"), !is.na(ocean_basin)) %>%
  count(haplotype_id, ocean_basin) %>%
  ggplot(aes(x = ocean_basin, y = n, fill = haplotype_id)) +
  geom_col()
```

(`geography/haplotype_by_region.csv` already holds these counts.)

## Tips and Best Practices

1. **Format selection**
   - Phylogeography and networks → PopART
   - Population structure → Arlequin
   - Diversity statistics → DnaSP
   - Custom analyses → Generic, or the run's `records.csv`

2. **Records without a region** are left out of Arlequin and DnaSP; see
   [Which Records Are Exported](#which-records-are-exported).

3. **Corrected records**: `assigned_corrected` records (error-filtered singletons
   assigned to their nearest haplotype) are exported with that haplotype's sequence.
   Use `--no-assign-corrected` to leave them out.

4. **Format details**: see `exports/README.md` in the run folder.

## Troubleshooting

### Issue: No sequences in export

**Cause**: no records were assigned to a haplotype.
**Solution**: check `summary.csv` and the statuses in `records.csv`.

### Issue: Arlequin or DnaSP file is empty or small

**Cause**: few records have a region (missing or centroid coordinates, or no shapefile).
**Solution**: check `coord_status` in `records.csv` and the warnings at the end of the log.
For GOaS, download from https://www.marineregions.org/downloads.php
(`GOaS_v1_20211214.zip`, registration form required), extract to
`shapefiles/GOaS_v1_20211214/`, and verify with `python -m boldgenotyper.goas_downloader`.

### Issue: Format-specific software won't open file

**Cause**: incorrect file format or path.
**Solution**: check `exports/README.md`; verify the software version.

## Advanced Usage

### Custom Population Grouping

Populations come from the region column. To export by another column, for example
country, call the exporter on the run's `records.csv`:

```python
from boldgenotyper import popgen_export
import pandas as pd

run = "Sphyrnidae_output"
records = pd.read_csv(f"{run}/records.csv")
records = records[records["status"].isin(["assigned", "assigned_corrected"])]

popgen_export.export_population_genetics_formats(
    df=records,
    consensus_fasta_path=f"{run}/haplotypes.fasta",
    output_dir="custom_exports",
    organism="Sphyrnidae",
    formats=["arlequin", "dnasp"],
    pop_column="country/ocean",
)
```

## Citation

If you use these exports in a publication, please cite:

1. **BOLDGenotyper:** see [CITATION.cff](../../CITATION.cff).

2. **BOLD Database:**
   Ratnasingham S, Hebert PDN (2007) BOLD: The Barcode of Life Data System.
   Molecular Ecology Notes 7: 355-364.

3. **Software-specific citations:**
   - Arlequin: Excoffier & Lischer (2010)
   - PopART: Leigh & Bryant (2015)
   - DnaSP: Rozas et al. (2017)

## Support

- GitHub: https://github.com/SymbioSeas/BOLDGenotyper/issues
- Email: symbioseas@outlook.com

---

**BOLDGenotyper version:** 2.0.0
**Module:** `boldgenotyper.popgen_export`
