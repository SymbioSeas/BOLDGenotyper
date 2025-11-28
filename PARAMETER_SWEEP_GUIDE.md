# Parameter Sweep Guide for BOLDGenotyper

**Author**: Steph Smith
**Last Updated**: 2025-11-25
**Version**: 1.0.0

---

## Table of Contents

- [Overview](#overview)
- [When to Use Parameter Sweep](#when-to-use-parameter-sweep)
- [Quick Start](#quick-start)
- [Understanding the Outputs](#understanding-the-outputs)
- [Interpreting Results](#interpreting-results)
- [Threshold Selection Strategy](#threshold-selection-strategy)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Publication Methods Text](#publication-methods-text)

---

## Overview

The clustering threshold is the most critical parameter in BOLDGenotyper, determining how genetically similar sequences must be to group into the same consensus genotype. Selecting an appropriate threshold is challenging because:

1. **Taxonomic variation**: Optimal thresholds vary by taxonomic group
2. **Marker variation**: COI divergence rates differ across taxa
3. **Study goals**: Fine-scale population studies require different thresholds than species delimitation
4. **Data quality**: Sequencing errors and alignment artifacts affect clustering

The `boldgenotyper-sweep` command addresses this challenge through **systematic empirical testing** of multiple threshold values to identify the optimal clustering parameter for your specific dataset.

### What Parameter Sweep Does

Parameter sweep runs the complete genotyping pipeline multiple times, each with a different clustering threshold, then analyzes how key metrics change across thresholds to identify:

- **Stability point**: Where grouping structure stabilizes (elbow point)
- **Assignment efficiency**: Threshold maximizing successful genotype assignments
- **Clustering coherence**: Consistency of sample groupings across threshold ranges
- **Identity scores**: Relationship between threshold and sequence identity

### Advantages Over Manual Selection

- **Data-driven**: Based on your actual sequences, not literature estimates
- **Reproducible**: Documented, systematic approach
- **Comprehensive**: Tests multiple thresholds simultaneously
- **Visual**: Clear plots showing optimal threshold
- **Automated**: No manual threshold testing required

---

## When to Use Parameter Sweep

### Recommended Scenarios

1. **Novel taxonomic groups**: No published COI threshold recommendations available
2. **Mixed datasets**: Multiple species or genera with different divergence rates
3. **Publication preparation**: Need to justify threshold choice with data
4. **Uncertain data quality**: Unknown level of sequencing error or variation
5. **Fine-scale studies**: Population genetics requiring precise threshold
6. **First-time analysis**: No prior experience with your organism group

### When Manual Selection is Sufficient

- Literature provides well-established thresholds for your exact taxon
- Pilot analysis with standard threshold (0.03) yields reasonable results
- Limited computational resources (sweep requires multiple pipeline runs)
- Exploratory analysis only (not for publication)

### Computational Considerations

Parameter sweep is **computationally intensive** because it runs the full pipeline for each threshold:

- 5 thresholds tested (default) = 5 complete pipeline runs
- Recommended: Run on computing cluster or overnight on workstation
- Option: Use `--threads` to parallelize each run
- Option: Test fewer thresholds initially (e.g., 3 values)

---

## Quick Start

### Step 1: Run Parameter Sweep

Basic usage with default thresholds (0.01, 0.015, 0.02, 0.03, 0.05):

```bash
boldgenotyper-sweep data/Sphyrna_lewini.tsv \
  --output parameter_sweep/
```

Custom threshold range:

```bash
boldgenotyper-sweep data/Sphyrna_lewini.tsv \
  --thresholds 0.005,0.01,0.015,0.02,0.025,0.03 \
  --threads 8 \
  --output sweep_results/
```

### Step 2: Review Recommendations

```bash
# Quick summary
cat parameter_sweep/recommendations.txt

# Detailed metrics
head parameter_sweep/sweep_summary.csv

# Visual analysis
open parameter_sweep/threshold_stability.pdf
open parameter_sweep/elbow_plot.pdf
```

### Step 3: Run Final Analysis with Optimal Threshold

```bash
# Use recommended threshold from sweep
boldgenotyper data/Sphyrna_lewini.tsv \
  --clustering-threshold 0.020 \
  --build-tree \
  --output final_analysis/
```

---

## Understanding the Outputs

### Directory Structure

```
parameter_sweep/
├── sweep_summary.csv              # Metrics for all thresholds tested
├── threshold_stability.pdf        # Multi-panel visualization
├── elbow_plot.pdf                # Optimal threshold detection
├── group_membership_tracking.csv  # Sample clustering stability
├── recommendations.txt            # Interpretation and suggested threshold
├── README.md                      # Auto-generated interpretation guide
└── runs/                          # Individual run outputs (if --keep-intermediates)
    ├── threshold_0_010/
    ├── threshold_0_015/
    ├── threshold_0_020/
    └── ...
```

### Key Output Files

#### 1. `sweep_summary.csv`

Summary metrics for each threshold tested:

```csv
threshold,n_consensus_groups,n_assigned,pct_assigned,mean_identity,median_identity,n_ties,n_low_confidence,n_singletons
0.010,45,823,87.3,96.5,97.2,12,34,8
0.015,32,891,94.5,95.8,96.5,8,28,4
0.020,24,903,95.8,95.2,96.1,5,22,2
0.030,18,897,95.1,94.5,95.6,6,25,1
0.050,12,856,90.8,93.1,94.8,9,31,0
```

**Key columns**:
- `n_consensus_groups`: Number of genotypes identified
- `n_assigned`: Samples successfully assigned
- `pct_assigned`: Assignment success rate
- `mean_identity`: Average sequence identity to assigned genotype
- `n_singletons`: Groups with only 1 sample (potential oversplitting)

#### 2. `threshold_stability.pdf`

Multi-panel figure showing:
- **Panel A**: Number of consensus groups vs. threshold (look for elbow)
- **Panel B**: Assignment rate vs. threshold (maximize this)
- **Panel C**: Mean identity vs. threshold (should be high and stable)
- **Panel D**: Singleton groups vs. threshold (minimize these)

**How to interpret**: The optimal threshold is often where:
- Panel A shows an elbow (curve flattens)
- Panel B is near maximum
- Panel C remains high (>90%)
- Panel D is low (<5% of groups)

#### 3. `elbow_plot.pdf`

Focused plot for elbow detection showing number of groups vs. threshold with:
- Smoothed curve
- Potential elbow point marked
- Rate of change indicators

**Elbow point**: Where the curve transitions from steep to gradual, indicating diminishing returns from stricter thresholds.

#### 4. `group_membership_tracking.csv`

Tracks how individual samples cluster across thresholds using **biologically meaningful** similarity metrics:

```csv
processid,t_0.010,t_0.015,t_0.020,t_0.030,n_changes,stability_score,mean_cluster_size
SAMPLE001,c15,c8,c5,c3,1,high,156.3
SAMPLE002,c15,c8,c5,c3,1,high,156.3
SAMPLE003,c42,c38,c22,c8,3,low,12.5
```

**Important columns**:
- `t_X.XXX`: Consensus group at each threshold (group IDs are arbitrary)
- `n_changes`: Number of times cluster **composition** changed significantly
- `stability_score`: Overall stability (high/medium/low)
- `mean_cluster_size`: Average cluster size across thresholds

**Note**: Group names (c15, c8, etc.) are arbitrary labels that change between thresholds. What matters is **which other samples cluster together**, not the group name. Changes are counted based on Jaccard similarity of cluster members, not group IDs.

#### 5. `recommendations.txt`

Auto-generated interpretation with:
- Recommended threshold and justification
- Data-driven rationale based on elbow detection
- Warnings about potential issues
- Suggested next steps

Example:

```
PARAMETER SWEEP RECOMMENDATIONS
================================

Recommended Threshold: 0.020
Confidence: High

Rationale:
- Elbow point detected at threshold 0.020
- Assignment rate: 95.8% (near maximum)
- Mean identity: 95.2% (stable and high)
- Singleton groups: 2 (minimal oversplitting)
- Clustering stability: 89% of samples show stable grouping

This threshold balances genotype resolution with assignment efficiency
and minimizes artifactual splitting due to sequencing variation.

Interpretation:
At this threshold, sequences with >98% identity are grouped together,
consistent with intraspecific variation in COI sequences for this
taxonomic group.

Next Steps:
1. Run full analysis with --clustering-threshold 0.020
2. Verify genotypes align with taxonomic expectations
3. Check group_membership_tracking.csv for unstable samples (n_changes > 2)
```

---

## Interpreting Results

### Strategy 1: Elbow Point Detection

The **elbow method** identifies where adding more groups provides diminishing returns.

**Visual inspection** (`threshold_stability.pdf`, Panel A):
1. Look for the point where the curve transitions from steep to gradual
2. This is your **optimal threshold**

**Characteristics of a good elbow**:
- Clear transition point (not gradual throughout)
- Stable assignment rate beyond the elbow
- High mean identity maintained beyond the elbow

**Example interpretation**:

```
Thresholds: 0.01, 0.015, 0.02, 0.03, 0.05
Groups:      45,  32,    24,  18,  12

Analysis:
0.01 → 0.015: Drops 13 groups (steep)
0.015 → 0.02: Drops 8 groups (moderate)
0.02 → 0.03:  Drops 6 groups (gradual)
0.03 → 0.05:  Drops 6 groups (gradual)

Elbow: ~0.02 (transition from steep to gradual)
```

### Strategy 2: Assignment Efficiency Maximization

Choose the threshold that maximizes successful genotype assignments while maintaining high identity scores.

**Decision criteria**:
1. Assignment rate >90%
2. Mean identity >90%
3. Low proportion of singletons (<5%)
4. Minimal ties and low-confidence assignments

**Example**:

```csv
threshold,pct_assigned,mean_identity,pct_singletons
0.010,87.3,96.5,17.8  # Too strict: many singletons
0.015,94.5,95.8,12.5  # Good assignment, some singletons
0.020,95.8,95.2,8.3   # ← OPTIMAL: Best balance
0.030,95.1,94.5,5.6   # Lower identity
0.050,90.8,93.1,0.0   # Too permissive: low identity
```

### Strategy 3: Clustering Stability Analysis

Examine `group_membership_tracking.csv` to identify samples with unstable clustering.

**Stability categories**:
- **High (n_changes = 0)**: Samples always cluster with the same partners (ideal)
- **Medium (n_changes = 1-2)**: Minor composition changes (acceptable)
- **Low (n_changes > 2)**: Unstable clustering (problem samples)

**Interpretation**:

```bash
# Count stability levels
cut -d',' -f7 group_membership_tracking.csv | sort | uniq -c

# Output:
# 756 high
#  89 medium
#  18 low
```

If >90% samples have high stability, the threshold range is appropriate.

**Problem samples** (n_changes > 2) may indicate:
1. Intermediate genotypes (genuine biological signal)
2. Low-quality sequences
3. Contamination or chimeric sequences
4. Alignment artifacts

**Action**: Review low-stability samples individually:

```bash
# Extract low-stability samples
awk -F',' '$7 == "low" {print $1}' group_membership_tracking.csv > problem_samples.txt

# Check these in diagnostics file
# Consider excluding from final analysis
```

### Strategy 4: Biological Validation

Compare results to taxonomic expectations:

**Expected patterns**:
- **Single species**: 1-10 genotypes depending on population structure
- **Closely related species**: Genotypes should separate by species at higher thresholds
- **Genus-level**: Many genotypes, species-specific clustering

**Red flags**:
- Too many genotypes (oversplitting): Consider higher threshold
- Single genotype for multi-species dataset: Threshold too permissive
- Taxonomically mixed genotypes: Re-evaluate species IDs or threshold

---

## Threshold Selection Strategy

### Decision Tree

```
START
  |
  ├─> Is there a clear elbow in the plot?
  |     YES → Use elbow threshold
  |     NO  → Continue to next criterion
  |
  ├─> Is assignment rate maximized at one threshold?
  |     YES → Use that threshold (if identity remains high)
  |     NO  → Continue to next criterion
  |
  ├─> Are >90% of samples stable (high/medium)?
  |     YES → Use threshold at stability transition
  |     NO  → Examine data quality issues
  |
  └─> Do results match taxonomic expectations?
        YES → Finalize threshold choice
        NO  → Re-evaluate or consult literature
```

### Conservative vs. Permissive Strategies

**Conservative (lower threshold, more groups)**:
- **Use when**: Population genetics, fine-scale structure, within-species variation
- **Advantages**: Captures subtle variation, avoids lumping distinct haplotypes
- **Disadvantages**: May oversplit due to sequencing errors
- **Example**: Threshold 0.01 for salmon population structure

**Permissive (higher threshold, fewer groups)**:
- **Use when**: Species delimitation, contamination screening, multi-species datasets
- **Advantages**: Robust to sequencing errors, clear genotype separation
- **Disadvantages**: May lump distinct haplotypes
- **Example**: Threshold 0.05 for genus-level diversity

### Handling Ambiguous Results

**Scenario 1: No clear elbow**

```
Thresholds: 0.01  0.015  0.02  0.03  0.05
Groups:     48    42     38    32    25

Problem: Gradual decline, no obvious elbow
```

**Solution**:
- Focus on assignment efficiency and identity scores
- Choose middle value (0.02) as conservative compromise
- Check if all samples belong to same species (gradual pattern expected)

**Scenario 2: Multiple potential elbows**

```
Thresholds: 0.01  0.015  0.02  0.03  0.05
Groups:     52    38     36    22    12

Analysis: Possible elbows at 0.015 and 0.03
```

**Solution**:
- First elbow (0.015): Fine-scale population structure
- Second elbow (0.03): Broader genotype grouping
- Choose based on research question:
  - Population genetics → 0.015
  - Species delimitation → 0.03

**Scenario 3: Conflicting signals**

```
Elbow suggests: 0.02
Max assignment: 0.03
Best identity: 0.015
```

**Solution**:
- Prioritize biological relevance
- Check taxonomic composition of groups
- Run final analysis with 2-3 candidates and compare manually
- Choose threshold where genotypes align with taxonomic units

---

## Advanced Usage

### Custom Threshold Ranges

Test very fine-scale resolution:

```bash
# Fine-grained population genetics
boldgenotyper-sweep data/population_samples.tsv \
  --thresholds 0.001,0.003,0.005,0.007,0.01,0.015 \
  --output fine_sweep/
```

Test broad taxonomic groups:

```bash
# Genus or family level
boldgenotyper-sweep data/diverse_taxa.tsv \
  --thresholds 0.03,0.05,0.07,0.10,0.15 \
  --output coarse_sweep/
```

### Keeping Intermediate Files

Retain full pipeline outputs for each threshold:

```bash
boldgenotyper-sweep data/samples.tsv \
  --keep-intermediates \
  --output sweep_with_runs/

# Access individual run outputs
ls sweep_with_runs/runs/threshold_0_020/
# Each contains complete pipeline outputs (phylogenetic trees, visualizations, etc.)
```

**Use cases**:
- Compare phylogenetic trees across thresholds
- Examine geographic patterns at different resolutions
- Extract specific genotype sequences from each run

### Parallel Execution

Speed up sweep with parallel threads:

```bash
# Use 16 cores for faster processing
boldgenotyper-sweep data/large_dataset.tsv \
  --threads 16 \
  --output sweep_parallel/

# Note: Each threshold run uses threads independently
# Total parallelization = --threads parameter
```

### Integrating Sweep Results with Main Analysis

Recommended workflow:

```bash
# 1. Run sweep to determine optimal threshold
boldgenotyper-sweep data/samples.tsv --output sweep/

# 2. Review recommendations
cat sweep/recommendations.txt

# 3. Run main analysis with optimal threshold
OPTIMAL_THRESHOLD=$(grep "Recommended" sweep/recommendations.txt | awk '{print $NF}')

boldgenotyper data/samples.tsv \
  --clustering-threshold $OPTIMAL_THRESHOLD \
  --build-tree \
  --export-format all \
  --output final_analysis/

# 4. Document threshold choice in methods
# Reference sweep results as justification
```

---

## Troubleshooting

### Problem: All thresholds yield very similar results

**Symptoms**:
- Number of groups changes by <10% across threshold range
- Flat lines in stability plots

**Possible causes**:
1. **Dataset too homogeneous**: All sequences are very similar (single population)
2. **Threshold range too narrow**: Not spanning the critical transition zone
3. **Data quality issues**: High sequencing error obscuring real variation

**Solutions**:

```bash
# Check sequence identity distribution
python3 -c "
import pandas as pd
from Bio import SeqIO, AlignIO
from scipy.spatial.distance import pdist, squareform
# Calculate pairwise distances
# Look for bimodal distribution (within vs. between groups)
"

# Try broader threshold range
boldgenotyper-sweep data/samples.tsv \
  --thresholds 0.005,0.02,0.05,0.10,0.15 \
  --output broad_sweep/
```

### Problem: Too many singletons at all thresholds

**Symptoms**:
- High proportion of 1-sample groups across all thresholds
- Many samples flagged as low-stability

**Possible causes**:
1. **Poor sequence quality**: Errors preventing clustering
2. **Highly diverse dataset**: Genuinely many unique genotypes
3. **Alignment issues**: Sequences not properly aligned

**Solutions**:

```bash
# Check for sequence length variation
python3 -c "
import pandas as pd
df = pd.read_csv('data/samples.tsv', sep='\t')
print(df['nucleotides'].str.len().describe())
"

# Filter short or poor-quality sequences
# Re-run pipeline with stricter quality filters
# Check alignment quality in intermediate files
```

### Problem: Recommendations file suggests "inconclusive"

**Meaning**: Automated analysis couldn't identify a clear optimal threshold

**Causes**:
- No clear elbow point
- Conflicting optimization criteria
- Unusual data structure

**Manual interpretation steps**:

1. **Visual inspection**: Examine plots yourself
2. **Literature comparison**: Check published thresholds for your taxon
3. **Biological validation**: Run 2-3 candidates and compare taxonomic coherence
4. **Conservative choice**: When in doubt, use middle value

### Problem: Sweep takes too long

**Typical runtimes** (1000 sequences, 5 thresholds):
- Single thread: 2-4 hours
- 8 threads: 30-60 minutes
- 16 threads: 15-30 minutes

**Solutions**:

```bash
# 1. Reduce number of thresholds
boldgenotyper-sweep data/samples.tsv \
  --thresholds 0.01,0.02,0.03 \  # Only 3 instead of 5
  --output quick_sweep/

# 2. Increase threads
boldgenotyper-sweep data/samples.tsv \
  --threads 16 \
  --output faster_sweep/

# 3. Subsample dataset initially
head -500 data/samples.tsv > data/samples_subset.tsv
boldgenotyper-sweep data/samples_subset.tsv --output pilot_sweep/
# Then run full sweep with refined threshold range
```

---

## Publication Methods Text

### Example Methods Section

```
Clustering threshold optimization was performed using the boldgenotyper-sweep
command, which systematically tested thresholds of 0.01, 0.015, 0.02, 0.03,
and 0.05. For each threshold, the complete genotyping pipeline was executed,
and metrics including number of consensus groups, assignment success rate,
mean sequence identity, and clustering stability were calculated.

The optimal threshold was determined using the elbow method (Smith, 2025),
identifying the point where grouping structure stabilized while maintaining
high assignment efficiency (>95%) and sequence identity (>90%). Based on this
analysis, a clustering threshold of 0.020 was selected, corresponding to
sequences with ≥98% identity being grouped into consensus genotypes.

Clustering stability was assessed by tracking sample groupings across all
tested thresholds using Jaccard similarity of cluster composition. A total of
89.3% of samples showed high stability (consistent clustering partners across
thresholds), indicating robust genotype delineation.

References:
Smith, S. (2025). BOLDGenotyper: Automated COI sequence genotyping and
biogeographic analysis. https://github.com/SymbioSeas/BOLDGenotyper
```

### Short Version (for supplementary methods)

```
Clustering threshold (0.020) was determined empirically using parameter sweep
analysis (boldgenotyper-sweep), which tested five threshold values and
identified the elbow point balancing genotype resolution and assignment
efficiency (see Supplementary Figure S1).
```

### Figure Legend for Sweep Plots

```
Figure S1: Parameter sweep analysis for clustering threshold optimization.
(A) Number of consensus genotypes identified at each threshold value. The
elbow point at 0.020 (arrow) indicates optimal threshold balancing resolution
and stability. (B) Assignment success rate showing near-maximum efficiency
(95.8%) at threshold 0.020. (C) Mean sequence identity to assigned genotypes
remains high (>95%) across threshold range. (D) Singleton groups (genotypes
with single sample) minimized at threshold 0.020, indicating appropriate
grouping without oversplitting. Vertical dashed line indicates selected
threshold.
```

---

## Best Practices

### Before Running Sweep

1. **Check data quality**: Review sequence lengths, N content, coordinate quality
2. **Estimate runtime**: Calculate based on dataset size and available threads
3. **Choose threshold range**: Start with defaults, refine if needed
4. **Document settings**: Record command used for reproducibility

### During Sweep

1. **Monitor progress**: Check log files periodically
2. **Watch for errors**: Ensure all threshold runs complete successfully
3. **Resource monitoring**: Verify sufficient disk space and memory

### After Sweep

1. **Review all outputs**: Don't rely solely on automated recommendations
2. **Visual inspection**: Always examine plots yourself
3. **Biological validation**: Verify genotypes match taxonomic expectations
4. **Document decision**: Save recommendations.txt and plots for methods section
5. **Test edge cases**: Run final analysis and spot-check problem samples

---

## Additional Resources

### Related Documentation

- Main pipeline: [README.md](README.md)
- Custom shapefiles: [CUSTOM_SHAPEFILES_GUIDE.md](CUSTOM_SHAPEFILES_GUIDE.md)
- Comparative analysis: See README.md "Advanced Workflows" section

### Relevant Literature

**Elbow method for clustering**:
- Thorndike, R. L. (1953). Who belongs in the family? Psychometrika, 18(4), 267-276.

**COI barcoding thresholds**:
- Hebert, P. D. N., Cywinska, A., Ball, S. L., & deWaard, J. R. (2003). Biological identifications through DNA barcodes. Proceedings of the Royal Society B, 270(1512), 313-321.
- Ward, R. D., Zemlak, T. S., Innes, B. H., Last, P. R., & Hebert, P. D. N. (2005). DNA barcoding Australia's fish species. Philosophical Transactions of the Royal Society B, 360(1462), 1847-1857.

**Clustering stability**:
- Ben-Hur, A., Elisseeff, A., & Guyon, I. (2002). A stability based method for discovering structure in clustered data. Pacific Symposium on Biocomputing, 7, 6-17.

---

## Support

For questions or issues with parameter sweep:

1. **Check documentation**: Review this guide and README.md
2. **Examine outputs**: recommendations.txt often provides specific guidance
3. **GitHub Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
4. **Email**: steph.smith@unc.edu

When reporting issues, include:
- Command used
- Summary statistics from sweep_summary.csv
- Description of unexpected results
- Plots if visual interpretation is unclear

---

**Created**: 2025-11-25
**Documentation version**: 1.0.0
**BOLDGenotyper version**: 0.1.0+

For main BOLDGenotyper documentation, see [README.md](README.md).
