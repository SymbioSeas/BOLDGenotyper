# Comparative Analysis Guide for BOLDGenotyper

**Author**: Steph Smith
**Last Updated**: 2025-11-25
**Version**: 1.0.0

---

## Table of Contents

- [Overview](#overview)
- [When to Use Comparative Analysis](#when-to-use-comparative-analysis)
- [Quick Start](#quick-start)
- [Understanding the Workflow](#understanding-the-workflow)
- [Interpreting Results](#interpreting-results)
- [Use Case Examples](#use-case-examples)
- [Troubleshooting](#troubleshooting)
- [Publication Methods Text](#publication-methods-text)

---

## Overview

Comparative analysis in BOLDGenotyper compares species-level and family-level (or higher taxonomic level) genotyping results to detect:

1. **Contamination**: Samples labeled as one species but clustering with another
2. **Mislabeling**: Taxonomic misidentifications in BOLD data
3. **Cryptic species**: Distinct genetic lineages within nominal species
4. **Quality control**: Validation of species-level assignments

### The Concept

When you analyze sequences at different taxonomic levels, you can identify samples that behave unexpectedly:

- **Species-level analysis**: Groups samples within a single species (e.g., *Sphyrna lewini*)
- **Family-level analysis**: Groups samples across multiple species (e.g., all Sphyrnidae)

If a sample labeled as Species A consistently clusters with Species B samples in the family-level analysis, this suggests:
- Contamination (wrong tissue sampled)
- Mislabeling (taxonomic error)
- Introgression (hybridization)

### Why This Matters

**For Data Quality**:
- BOLD data quality varies by contributor
- Contamination can occur during tissue sampling or lab processing
- Taxonomic expertise varies across research groups

**For Biological Discovery**:
- May reveal cryptic species
- Identifies gene flow between species
- Validates taxonomic classifications

**For Publication**:
- Strengthens confidence in results
- Provides quality control documentation
- May identify samples to exclude or re-examine

---

## When to Use Comparative Analysis

### Recommended Scenarios

1. **Multi-species studies**: Analyzing a family or genus with multiple species
2. **Quality control**: Preparing data for publication
3. **Contamination screening**: Suspect samples with unusual distributions
4. **Taxonomic validation**: Verifying species identifications
5. **Cryptic diversity**: Searching for undescribed lineages
6. **Mixed datasets**: BOLD downloads with questionable taxonomy

### Prerequisites

To run comparative analysis, you need:

1. **Species-level analysis completed**: At least one species analyzed individually
2. **Family/genus-level analysis completed**: Higher taxonomic group including that species
3. **Overlapping samples**: Same specimens present in both datasets
4. **Consistent thresholds**: Ideally use same or similar clustering thresholds

### What You Can Detect

**Definitive signals**:
- Sample labeled Species A, clusters 100% with Species B in family analysis
- Clear contamination or mislabeling

**Ambiguous signals**:
- Sample labeled Species A, clusters 70% with Species B
- Could be contamination, mislabeling, or genuine hybrid

**Expected patterns**:
- Most samples cluster with their labeled species (>95%)
- Occasional low-confidence assignments (acceptable)

---

## Quick Start

### Step 1: Run Both Analyses

Run species-level and family-level analyses separately:

```bash
# Species-level analysis
boldgenotyper data/Sphyrna_lewini.tsv \
  --clustering-threshold 0.02 \
  --output Sphyrna_lewini_output/

# Family-level analysis
boldgenotyper data/Sphyrnidae.tsv \
  --clustering-threshold 0.03 \
  --output Sphyrnidae_output/
```

**Important**: Use appropriate thresholds for each level (species-level may need stricter threshold).

### Step 2: Run Comparative Analysis

```bash
python scripts/compare_analyses.py \
  --species-level Sphyrna_lewini_output/ \
  --family-level Sphyrnidae_output/ \
  --output comparative_analysis/
```

With detailed reassignment table:

```bash
python scripts/compare_analyses.py \
  --species-level Sphyrna_lewini_output/ \
  --family-level Sphyrnidae_output/ \
  --generate-reassignment-table \
  --output comparative_analysis/
```

### Step 3: Review Results

```bash
# Quick summary
cat comparative_analysis/comparison_summary.csv

# Genotype relationships
cat comparative_analysis/genotype_crosswalk.csv

# Sample-level details (if generated)
cat comparative_analysis/sample_reassignments.csv

# Publication methods text
cat comparative_analysis/methods_text.md
```

---

## Understanding the Workflow

### What Comparative Analysis Does

1. **Load both datasets**: Species-level and family-level annotated CSVs
2. **Match samples**: Identify samples present in both analyses
3. **Compare assignments**: Check if samples cluster consistently
4. **Analyze taxonomy**: Determine majority species for each family-level genotype
5. **Flag discrepancies**: Identify samples with unexpected assignments
6. **Generate reports**: Summary statistics, crosswalk tables, reassignment recommendations

### Input Requirements

**Species-level directory** must contain:
- `{organism}_annotated.csv`: Main results with genotype assignments
- Optional: `genotype_assignments/{organism}_diagnostics.csv` for identity scores
- Optional: `taxonomy/{organism}_consensus_taxonomy.csv` for species composition

**Family-level directory** must contain:
- Same structure as species-level
- Should include samples from the species analysis

**Overlapping samples**: The comparison only analyzes samples present in both datasets.

### Algorithm Details

**Genotype Crosswalk**:
1. For each species-level genotype, find all samples assigned to it
2. For each sample, check its family-level genotype assignment
3. Calculate percentage of samples from each species-level genotype that map to each family-level genotype
4. Identify dominant mapping pattern

**Contamination Detection**:
1. For each sample, compare species label to majority species in family-level genotype
2. If mismatch occurs and identity is high, flag as potential contamination
3. Calculate confidence based on:
   - Percentage match to alternative species
   - Identity score at family level
   - Presence of multiple species in family genotype

---

## Interpreting Results

### Output Files

#### 1. `comparison_summary.csv`

High-level metrics comparing the two analyses:

```csv
metric,species_level,family_level
total_samples,617,995
consensus_groups,10,17
species_detected,1,7
mean_identity,95.30,93.64
groups_with_100pct_majority,8,12
mixed_species_groups,0,2
potential_misidentifications,0,3
```

**Key metrics**:
- **consensus_groups**: Number of genotypes at each level (family should have more)
- **species_detected**: Number of distinct species (family should have more)
- **mixed_species_groups**: Family-level genotypes with multiple species (expected)
- **potential_misidentifications**: Samples flagged for review

**Interpretation**:
- If `potential_misidentifications` = 0: Clean dataset
- If `potential_misidentifications` > 0: Review flagged samples
- If `mixed_species_groups` at species level > 0: Possible contamination or cryptic diversity

#### 2. `genotype_crosswalk.csv`

Shows how species-level genotypes map to family-level genotypes:

```csv
species_genotype,n_samples,family_genotype,n_in_family,pct_overlap,dominant_species_family,species_match
c1_n156,156,c15_n386,156,100.0,Sphyrna lewini,TRUE
c2_n98,98,c15_n386,98,100.0,Sphyrna lewini,TRUE
c3_n45,45,c8_n67,45,100.0,Sphyrna lewini,TRUE
c4_n12,12,c18_n15,10,83.3,Sphyrna zygaena,FALSE
```

**Key columns**:
- `species_genotype`: Genotype from species-level analysis
- `family_genotype`: Corresponding family-level genotype
- `pct_overlap`: Percentage of species samples in this family genotype
- `dominant_species_family`: Majority species in the family genotype
- `species_match`: Do species names match? (TRUE = expected, FALSE = investigate)

**Red flags**:
- `species_match = FALSE`: Species-level genotype clusters with different species at family level
- `pct_overlap < 80%`: Samples from one species genotype split across multiple family genotypes
- Multiple species genotypes → same family genotype with low overlap: Potential oversplitting at species level

#### 3. `sample_reassignments.csv` (if `--generate-reassignment-table`)

Sample-level details for all flagged specimens:

```csv
processid,species_label,species_genotype,family_genotype,family_dominant_species,family_species_pct,identity_family,recommendation,confidence
SAMPLE12345,Sphyrna lewini,c4_n12,c18_n15,Sphyrna zygaena,93.3,0.96,Review: potential mislabel,High
SAMPLE67890,Sphyrna lewini,c7_n8,c22_n10,Sphyrna tiburo,80.0,0.92,Review: possible contamination,Medium
```

**Key columns**:
- `species_label`: Original BOLD species name
- `family_dominant_species`: Species that dominates the family-level genotype
- `family_species_pct`: Percentage of that family genotype from the dominant species
- `recommendation`: Automated interpretation
- `confidence`: High/Medium/Low based on evidence strength

**Recommendations**:
- **"Retain"**: Sample behaves as expected
- **"Review: potential mislabel"**: Strong evidence for taxonomic error
- **"Review: possible contamination"**: Moderate evidence for sample mix-up
- **"Exclude"**: Strong evidence for problem (contamination + low identity)

#### 4. `methods_text.md`

Auto-generated methods section ready for publication:

```markdown
## Comparative Analysis Methods

To assess data quality and detect potential contamination, we performed
comparative analysis between species-level (Sphyrna lewini, n=617) and
family-level (Sphyrnidae, n=995) genotyping using the python scripts/compare_analyses.py
command (Smith, 2025).

[... detailed methods with parameter values ...]

Results: Of 617 samples analyzed at both taxonomic levels, 614 (99.5%)
clustered consistently with their species designation. Three samples showed
discordant clustering patterns and were flagged for further review...
```

### Decision Framework

**For each flagged sample, consider**:

1. **Evidence strength**:
   - High: >90% family genotype from alternative species, identity >95%
   - Medium: 70-90% alternative species, identity >90%
   - Low: <70% or identity <90%

2. **Biological plausibility**:
   - Are the two species sympatric (occur together)?
   - Is contamination plausible in sampling context?
   - Could this be hybridization or introgression?

3. **Action**:
   - **High confidence contamination**: Exclude from analysis
   - **Medium confidence**: Review original metadata, consider excluding
   - **Low confidence or biological interest**: Investigate further, possibly retain

---

## Use Case Examples

### Example 1: Clean Dataset (No Issues)

**Setup**:
```bash
boldgenotyper data/Sphyrna_lewini.tsv --output species/
boldgenotyper data/Sphyrnidae.tsv --output family/
python scripts/compare_analyses.py --species-level species/ --family-level family/ --output compare/
```

**Results** (`comparison_summary.csv`):
```csv
metric,species_level,family_level
potential_misidentifications,0,0
```

**Interpretation**: All samples cluster as expected. No contamination detected. Dataset is clean.

**Next steps**: Proceed with analysis confidently. Mention quality control in methods.

---

### Example 2: Contamination Detected

**Results** (`sample_reassignments.csv`):
```csv
processid,species_label,family_dominant_species,family_species_pct,confidence
GBSHK001,Sphyrna lewini,Sphyrna tiburo,95.2,High
GBSHK002,Sphyrna lewini,Sphyrna tiburo,95.2,High
```

**Interpretation**: Two samples labeled as *S. lewini* cluster strongly with *S. tiburo* (bonnethead shark). High confidence contamination.

**Investigation**:
1. Check collection metadata: Were these samples collected alongside *S. tiburo*?
2. Review original photos/vouchers if available
3. Check for field notes mentioning possible mix-ups

**Action**: Exclude GBSHK001 and GBSHK002 from *S. lewini* analysis. Document in methods:

```
Two samples (GBSHK001, GBSHK002) initially labeled as S. lewini were excluded
based on comparative analysis indicating consistent clustering with S. tiburo
genotypes (>95% at family level, sequence identity >96%), suggesting
contamination or mislabeling.
```

---

### Example 3: Cryptic Diversity Discovery

**Results** (`genotype_crosswalk.csv`):
```csv
species_genotype,family_genotype,dominant_species_family,species_match
c1_n234,c5_n234,Sphyrna lewini,TRUE
c2_n45,c18_n45,Sphyrna lewini,FALSE → Actually Sphyrna gilberti
```

**Interpretation**: One species-level genotype (c2) clusters separately at family level and shows distinct molecular characteristics. Potential cryptic species.

**Next steps**:
1. Examine geographic distribution of c2 samples
2. Review morphological characteristics
3. Consider phylogenetic analysis with nuclear markers
4. Consult taxonomic literature on *S. gilberti* (cryptic species)

**Publication**: This becomes a biological finding, not a quality control issue.

---

### Example 4: Hybrid Detection

**Results**:
```csv
processid,species_label,family_dominant_species,family_species_pct,identity_family
HYBRID01,Sphyrna lewini,Sphyrna lewini,52.0,0.94
```

**Interpretation**: Sample shows intermediate clustering. Family genotype contains 52% *S. lewini* and 48% *S. tiburo* (hypothetically).

**Investigation**:
1. Check if species ranges overlap (sympatry required for hybridization)
2. Review morphology for intermediate characteristics
3. Consider nuclear DNA analysis
4. Look for other potential hybrids

**Action**: Flag as "potential hybrid" rather than contamination. Retain in dataset but annotate.

---

### Example 5: Multi-Species Dataset Quality Control

**Setup**: Family-level dataset with 7 species

**Results** (`comparison_summary.csv`):
```csv
species,samples,misidentifications
Sphyrna lewini,617,2
Sphyrna tiburo,145,0
Sphyrna zygaena,98,1
Sphyrna mokarran,76,0
```

**Interpretation**: Different species show different error rates. *S. lewini* and *S. zygaena* have issues.

**Pattern analysis**:
- Check if problems cluster by geographic region (suggests one collection)
- Check if problems cluster by processid prefix (suggests one contributor/lab)
- Check collection dates for temporal patterns

**Action**: Exclude problematic samples. Consider contacting BOLD curators about systematic issues.

---

## Troubleshooting

### Problem: "No overlapping samples found"

**Cause**: Species-level samples not present in family-level dataset

**Solution**:
```bash
# Check which samples are in each dataset
cut -d',' -f1 species/Sphyrna_lewini_annotated.csv | sort > species_ids.txt
cut -d',' -f1 family/Sphyrnidae_annotated.csv | sort > family_ids.txt
comm -12 species_ids.txt family_ids.txt | wc -l

# If count = 0, datasets don't overlap
# Verify you're using correct files
```

**Common mistakes**:
- Different BOLD downloads at different times
- Species TSV doesn't include all species samples
- Filtered out samples at different stages

---

### Problem: All samples flagged as misidentified

**Cause**: Likely using incompatible analyses (different thresholds, different markers)

**Solution**:
```bash
# Check thresholds used
grep "clustering" species/Sphyrna_lewini_pipeline_parameters.json
grep "clustering" family/Sphyrnidae_pipeline_parameters.json

# If very different (e.g., 0.01 vs 0.05), re-run with similar thresholds
boldgenotyper data/Sphyrna_lewini.tsv --clustering-threshold 0.03 --output species_v2/
boldgenotyper data/Sphyrnidae.tsv --clustering-threshold 0.03 --output family_v2/
```

---

### Problem: Unexpected "mixed species groups" at species level

**Symptom**: comparison_summary.csv shows mixed_species_groups > 0 for species-level analysis

**Interpretation**: This shouldn't happen for single-species dataset

**Possible causes**:
1. **Contamination**: Some samples are actually from other species
2. **Dataset contamination**: Species TSV contains multiple species
3. **Taxonomy errors**: BOLD records have wrong species labels

**Investigation**:
```bash
# Check species diversity in your "species-level" dataset
cut -d',' -f{species_column} species/Sphyrna_lewini_annotated.csv | sort | uniq -c

# If multiple species present, you have a multi-species dataset
# Run comparative analysis to identify which samples are contaminants
```

---

### Problem: Low identity scores in family analysis

**Symptom**: Many samples with identity <90% in family-level analysis

**Cause**: Family-level analysis groups divergent species, reducing overall identity

**Interpretation**: This is **normal and expected**. Different species have different COI sequences.

**Action**: Focus on **dominant species** patterns, not absolute identity. A sample can have 85% identity to family-level genotype but still clearly cluster with its species if that's the highest identity among options.

---

### Problem: Genotype crosswalk shows many-to-one mapping

**Symptom**: Multiple species genotypes → single family genotype

**Interpretation**: Could indicate:
1. **Oversplitting** at species level (threshold too strict)
2. **Normal population structure**: Multiple haplotypes within species
3. **Expected pattern**: Family genotype represents the species clade

**Analysis**:
```bash
# Check how many species genotypes map to each family genotype
awk -F',' 'NR>1 {print $3}' genotype_crosswalk.csv | sort | uniq -c | sort -rn

# If one family genotype contains ALL species genotypes, this is expected
# If family genotypes are fragmented, investigate
```

---

## Advanced Analysis

### Quantifying Contamination Rates

Calculate overall dataset quality:

```python
import pandas as pd

# Load reassignment table
df = pd.read_csv('comparative_analysis/sample_reassignments.csv')

# Count issues by confidence
confidence_counts = df['confidence'].value_counts()
total_samples = len(df)

contamination_rate = (df['recommendation'].str.contains('Review|Exclude').sum() / total_samples) * 100

print(f"Contamination rate: {contamination_rate:.2f}%")
print(f"High confidence issues: {confidence_counts.get('High', 0)}")
print(f"Medium confidence issues: {confidence_counts.get('Medium', 0)}")
```

### Geographic Patterns in Contamination

Check if contamination clusters geographically:

```python
# Load full annotated data
species_df = pd.read_csv('species/Sphyrna_lewini_annotated.csv')
reassign_df = pd.read_csv('comparative_analysis/sample_reassignments.csv')

# Merge to get coordinates
merged = species_df.merge(reassign_df, on='processid')

# Group by country or ocean basin
contamination_by_region = merged.groupby('country')['recommendation'].apply(
    lambda x: (x.str.contains('Review|Exclude').sum() / len(x)) * 100
)

print(contamination_by_region.sort_values(ascending=False))
```

### Temporal Patterns

Check if contamination varies by collection period:

```python
# If you have collection dates
merged['year'] = pd.to_datetime(merged['collection_date']).dt.year

contamination_by_year = merged.groupby('year')['recommendation'].apply(
    lambda x: (x.str.contains('Review|Exclude').sum() / len(x)) * 100
)

# Plot trend over time
import matplotlib.pyplot as plt
contamination_by_year.plot(kind='line', marker='o')
plt.ylabel('Contamination rate (%)')
plt.xlabel('Collection year')
plt.title('Temporal pattern in sample contamination')
plt.savefig('contamination_trend.pdf')
```

---

## Publication Methods Text

### Full Methods Section

```
Data Quality Assessment

To assess data quality and detect potential sample contamination or
mislabeling, we performed comparative taxonomic analysis using the
compare_analyses.py utility (Smith, 2025). This approach compares
haplotype assignment results obtained at different taxonomic levels to identify
samples that cluster inconsistently with their taxonomic designation.

We first analyzed COI sequences from Sphyrna lewini (n=617) independently
using the ESV haplotype discovery approach, which generated 10 consensus
haplotypes. Subsequently, we analyzed a broader dataset containing all
Sphyrnidae species (n=995), generating 17 consensus haplotypes.

For the 617 samples present in both analyses, we compared species-level
and family-level genotype assignments. Samples were flagged as potentially
problematic if they consistently clustered with a different species in the
family-level analysis with high sequence identity (>90%) and strong
majority assignment (>80% of co-clustering samples from the alternative
species).

Results: Of 617 S. lewini samples analyzed, 614 (99.5%) showed consistent
clustering patterns at both taxonomic levels, supporting their species
designation. Three samples (0.5%) showed evidence of contamination or
mislabeling, clustering primarily with S. tiburo (n=2, >95% S. tiburo in
family genotype) or S. zygaena (n=1, 87% S. zygaena). These samples were
excluded from downstream analyses. No evidence of cryptic diversity or
hybridization was detected.

This quality control step provides confidence that our dataset accurately
represents S. lewini genetic diversity and is not confounded by
taxonomic errors or sample contamination.

References:
Smith, S. (2025). BOLDGenotyper: Automated COI sequence genotyping and
biogeographic analysis. https://github.com/SymbioSeas/BOLDGenotyper
```

### Supplementary Table Legend

```
Table S4: Sample-level reassignment recommendations from comparative
analysis. For each sample, we report the species designation from BOLD
records, the genotype assignment at species level, the genotype assignment
at family level, the dominant species within the family-level genotype,
the percentage of that genotype composed of the dominant species, the
sequence identity to the family-level consensus, and an automated
recommendation based on these metrics. Samples flagged for review show
evidence of contamination or mislabeling based on clustering with
alternative species at the family level.
```

---

## Best Practices

### Before Running Comparison

1. **Use similar thresholds**: Species and family thresholds should be within 0.01-0.02 of each other
2. **Complete both analyses**: Ensure both species and family runs completed successfully
3. **Check for overlap**: Verify species samples are actually present in family dataset
4. **Document parameters**: Save parameters used for both analyses

### During Interpretation

1. **Start with summary**: Review comparison_summary.csv first
2. **Check crosswalk**: Understand how genotypes relate
3. **Prioritize high confidence**: Focus on clear signals first
4. **Consider biology**: Not all discrepancies are errors (hybrids, cryptic species)
5. **Verify with metadata**: Cross-check flagged samples with collection data

### For Publication

1. **Document exclusions**: Clearly state which samples were removed and why
2. **Report statistics**: Contamination rate, number flagged, confidence levels
3. **Provide supplementary table**: sample_reassignments.csv as Table S4
4. **Include methods text**: Use auto-generated methods as starting point
5. **Archive original data**: Retain both cleaned and uncleaned datasets

---

## Additional Resources

### Related Documentation

- Main pipeline: [README.md](README.md)
- Parameter optimization: [PARAMETER_SWEEP_GUIDE.md](PARAMETER_SWEEP_GUIDE.md)
- Custom shapefiles: [CUSTOM_SHAPEFILES_GUIDE.md](CUSTOM_SHAPEFILES_GUIDE.md)

### Relevant Literature

**DNA barcode contamination**:
- Huemer, P., et al. (2014). Testing DNA barcode performance in 1000 species: a European perspective. *Molecular Ecology Resources*, 14(6), 1176-1191.

**Taxonomic validation**:
- Collins, R. A., & Cruickshank, R. H. (2013). The seven deadly sins of DNA barcoding. *Molecular Ecology Resources*, 13(6), 969-975.

**Quality control**:
- Pentinsaari, M., Salmela, H., Mutanen, M., & Roslin, T. (2016). Molecular evolution of a widely-adopted taxonomic marker (COI) across the animal tree of life. *Scientific Reports*, 6, 35275.

---

## Support

For questions about comparative analysis:

1. **Check this guide**: Review interpretation sections
2. **Examine outputs**: methods_text.md often provides context-specific guidance
3. **GitHub Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues
4. **Email**: symbioseas@outlook.com

When reporting issues, include:
- Command used
- Comparison summary statistics
- Unexpected patterns observed
- Research context (species, geographic scope)

---

**Created**: 2025-11-25
**Documentation version**: 1.0.0
**BOLDGenotyper version**: 0.1.0+

For main BOLDGenotyper documentation, see [README.md](README.md).
