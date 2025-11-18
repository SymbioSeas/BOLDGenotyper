---
title: 'BOLDGenotyper: A reproducible pipeline for COI dereplication, genotyping, and global biogeographic analysis using BOLD reference data'
tags:
  - Python
  - bioinformatics
  - DNA barcoding
  - COI
  - BOLD
  - genotyping
  - biogeography
authors:
  - name: Steph Smith
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: "1"
  - name: Chelsea Black
    orcid: 0000-0000-0000-0000
    affiliation: "1"
affiliations:
 - name: University of North Carolina at Chapel Hill, Institute of Marine Sciences, USA
   index: 1
date: 2025-XX-XX
bibliography: paper.bib
---

# Summary

DNA barcoding using the mitochondrial cytochrome oxidase I (COI) gene has become a foundational tool for species identification, population genetics, and biodiversity research. With more than two million COI sequences now publicly available, large-scale comparative analyses increasingly rely on robust, transparent workflows that can parse raw BOLD downloads, delineate unique haplotypes, and quantify biogeographic structure across global datasets. However, researchers often resort to ad hoc scripts or manual data cleaning steps, which complicate reproducibility and hinder cross-study comparison.

**BOLDGenotyper** is a Python-based pipeline for automated COI dereplication, consensus sequence generation, genotype assignment, ocean-basin classification, and global biogeographic visualization using sequence records directly exported from the Barcode of Life Database (BOLD). The workflow standardizes COI meta-analysis across taxa and environments, providing an entirely reproducible end-to-end process from raw BOLD metadata to publication-ready figures and statistical outputs.

BOLDGenotyper supports any organism represented in BOLD and has already been applied in population and conservation genetic analyses, including global COI assessments of the scalloped hammerhead shark *Sphyrna lewini* (Smith & Black, in review). The pipeline offers modular components, enabling users to run only genotyping, only phylogeny, or full biogeography depending on project needs.

# Statement of need

Despite the widespread use of COI barcoding, there is no widely adopted open-source workflow for **dereplicating sequences**, **assigning haplotype/genotype clusters**, and **quantifying global spatial patterns** using BOLD data. Many studies rely on custom scripts, proprietary tools, or inconsistent data-cleaning approaches that are difficult to reproduce and rarely portable across organisms.

BOLDGenotyper addresses this need by providing:

1. **Standardized, transparent COI processing**, including clustering, consensus calling, and sequence identity scoring.  
2. **Marine-oriented biogeographic classification**, using the Global Oceans and Seas (GOaS) shapefile to assign samples to ocean basins.  
3. **Publication-quality geographic visualizations**, including genotype maps, basin-level summaries, and territory polygons.  
4. **Statistical outputs**, including contingency tables and chi-square tests for genotype–basin associations.  
5. **Organism-agnostic design**, enabling analyses across marine, freshwater, and terrestrial taxa (with optional `--no-geo` mode).  
6. **Reproducible, documented code**, with optional diagnostics and parallelization.

The tool is intended for researchers working in biodiversity monitoring, fisheries management, conservation genetics, environmental DNA, and population genomics—fields where COI remains a key marker and reproducibility is increasingly prioritized.

# Functionality

BOLDGenotyper implements the following main modules:

### **1. BOLD TSV parsing**
Converts standard BOLD downloads into FASTA and structured metadata tables, validating required columns and formatting record identifiers in a standardized way.

### **2. Sequence dereplication and consensus clustering**
- Aligns sequences with MAFFT [@mafft]  
- Trims alignments using trimAl [@trimal]  
- Computes pairwise identity using informative sites only  
- Performs hierarchical clustering (default 99% identity threshold)  
- Generates consensus sequences using majority-rule aggregation  

Consensus sequences are labeled with genotype IDs (e.g., `c7_n381`) that reflect both cluster identity and sample count.

### **3. Genotype assignment**
Each input sequence is assigned to its closest consensus cluster via global alignment using the `edlib` library. Assignment confidence, runner-up hits, and optional diagnostics are recorded.

### **4. Geographic filtering and ocean-basin assignment**
- Uses the Global Oceans and Seas shapefile (IHO v3)  
- Filters ambiguous or country-centroid coordinates  
- Detects land points and supports multiple land policies  
- Assigns coordinates to ocean basins and major marine regions  

### **5. Visualization and statistical analysis**
Generates:
- Global genotype maps (PDF/SVG)  
- ocean-basin stacked barplots  
- standardized residuals heatmaps  
- territory polygons based on sample clusters  

Paired R scripts provide chi-square tests, residual analyses, and effect-size measures.

### **6. Reproducibility**
The pipeline is fully modular and callable via command-line scripts or individual Python modules, enabling integration into larger analyses or workflow managers.

# Example usage

After downloading a BOLD TSV file (`Genus_species_commonname.tsv`), the typical workflow is:

```bash
python BOLD_tsv_to_fasta.py data/Sphyrna_lewini.tsv
python msa_to_consensus.py data/
python consensus_group_to_metadata.py \
    --metadata data/Sphyrna_lewini.tsv \
    --fasta data/Sphyrna_lewini.fasta \
    --consensus data/consensus/Sphyrna_lewini_consensus.fasta \
    --out data/Sphyrna_lewini_with_consensus.tsv

python plot_shark_genotypes_ocean_basins_complex.py \
    --csv data/Sphyrna_lewini_with_consensus.tsv \
    --outdir results/