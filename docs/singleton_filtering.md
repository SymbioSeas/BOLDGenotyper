# Singleton Error Filtering in BOLDGenotyper

## Overview

BOLDGenotyper implements two-stage singleton quality control to remove sequencing/PCR error-derived haplotypes while preserving genuine biological variation.

## Rationale

### The Problem: Singleton Inflation

Exact Sequence Variant (ESV) approaches identify all unique sequences, including those created by technical errors. Sequencing and PCR errors create false singleton haplotypes that differ by 1-3 bp from true haplotypes:

**Sequencing errors:**
- Illumina: ~0.1-0.3% error rate (Schirmer et al. 2015)
- Sanger: ~0.1% error rate  
- In 650 bp sequence: ~1-2 errors per read

**PCR errors:**
- Taq polymerase: ~0.01-0.1% error rate per bp
- Accumulated across PCR cycles

**Empirical evidence (Sphyrnidae dataset):**
- 198 total haplotypes identified
- 156 singletons (78.8%)
- **85.9% of singletons differ by ≤0.3% from nearest neighbor**
- Median singleton divergence: 0.0%
- 14x more haplotypes than BOLD BINs (2% clustering)

This pattern strongly indicates error-driven singleton inflation, not biological diversity.

## Two-Stage Quality Control

### Stage 1: Error Filtering (NEW)
**Applied:** After haplotype identification, before suspect flagging  
**Default:** ENABLED (`min_singleton_distance = 0.005`)  
**Action:** Removes singletons ≤0.5% divergent from nearest neighbor

**Threshold selection:**
- 0.5% divergence = ~3 bp difference in 650 bp sequence
- Conservative: above typical sequencing error rate
- Empirical validation: removes ~85% of singletons
- Preserves rare variants and geographic structure

**Impact on Sphyrnidae dataset:**
- Before filtering: 198 haplotypes (156 singletons)
- After filtering: 58 haplotypes (22 singletons)
- 70% reduction, closer to 14 BINs (4.1 haplotypes/BIN)

### Stage 2: Suspect Flagging
**Applied:** After error filtering  
**Default:** DISABLED (`flag_suspect_haplotypes = False`)  
**Action:** Flags (but doesn't remove) singletons >5% divergent

**Purpose:** Identify potential contamination or misidentification

## Configuration

### Using Default Settings (Recommended)

```python
# Default configuration automatically applies 0.5% filtering
boldgenotyper data/sequences.tsv
```

### Custom Thresholds

```python
from boldgenotyper.config import HaplotypeConfig

# More conservative filtering (1% threshold)
config = HaplotypeConfig(min_singleton_distance=0.01)

# Disable filtering (keep all singletons)
config = HaplatypeConfig(min_singleton_distance=0.0)

# Enable suspect flagging
config = HaplotypeConfig(
    min_singleton_distance=0.005,  # Error filter at 0.5%
    max_singleton_distance=0.05,   # Suspect flag at 5%
    flag_suspect_haplotypes=True
)
```

### Threshold Selection Guide

| Threshold | % Singletons Removed | Use Case |
|-----------|---------------------|----------|
| 0.003 (0.3%) | ~85% | Very conservative, removes likely errors only |
| **0.005 (0.5%)** | **~90%** | **DEFAULT - recommended balance** |
| 0.01 (1%) | ~93% | Aggressive filtering, may remove rare variants |
| 0.02 (2%) | ~95% | Very aggressive, approaches BIN clustering |

## Validation Approaches

### 1. Compare to BOLD BINs
BINs represent 2% clustering. ESV haplotypes should be 2-5x BIN count, not 10-15x.

### 2. Check Singleton Distribution
Plot singleton divergence histogram. Bimodal distribution indicates:
- Peak at 0-0.3%: sequencing errors
- Peak at 1-5%: rare variants, geographic structure

### 3. Geographic/Species Patterns
Real singletons cluster by geography or species. Error singletons appear randomly distributed.

### 4. Phylogenetic Validation
Real haplotypes cluster with species/populations. Error singletons appear as long branches scattered across the tree.

## Implementation Details

### Filtering Algorithm

```python
def filter_error_singletons(haplotype_records, haplotype_mapping, min_singleton_distance):
    """
    For each singleton haplotype:
    1. Calculate distance to nearest neighbor
    2. If distance <= min_singleton_distance:
       - Remove haplotype from records
       - Remove samples from mapping
    3. Else: keep haplotype
    
    Multi-member haplotypes always retained (never filtered).
    """
```

### Pipeline Integration

```
Input FASTA
    ↓
MAFFT alignment
    ↓  
Core region extraction
    ↓
Identify unique haplotypes (ESVs)
    ↓
**Filter error singletons** ← NEW STEP
    ↓
Flag suspect haplotypes (optional)
    ↓
Output filtered haplotypes
```

## References

**Sequencing error rates:**
- Schirmer et al. (2015). Insight into biases and sequencing errors for amplicon sequencing with the Illumina MiSeq platform. *Nucleic Acids Research*, 43(6), e37.

**PCR and ESV approaches:**
- Edgar (2013). UPARSE: highly accurate OTU sequences from microbial amplicon reads. *Nature Methods*, 10(10), 996-998.

**ESV best practices:**
- Porter & Hajibabaei (2020). Scaling up: A guide to high-throughput genomic approaches for biodiversity analysis. *Molecular Ecology*.

## FAQ

**Q: Why not just use clustering (BINs)?**  
A: Clustering at 2% loses fine-scale variation important for population genetics and species delimitation. Error filtering + ESVs preserves biological signal while removing technical noise.

**Q: Can I disable filtering?**  
A: Yes, set `min_singleton_distance=0.0`. However, this is not recommended for most applications due to singleton inflation.

**Q: How do I know if filtering is too aggressive?**  
A: Check if your haplotype count drops below your BIN count, or if you lose expected geographic/species structure. Empirically, 0.5% is conservative.

**Q: What about multi-member haplotypes with errors?**  
A: Errors creating identical sequences in multiple samples is extremely unlikely. Multi-member haplotypes are almost certainly real.

**Q: Should I also filter by ORF validation?**  
A: Yes! ORF filtering (enabled by default) removes NUMTs and contamination. Singleton error filtering is complementary and addresses a different artifact.
