# BOLDGenotyper Flow Chart Comparison

This document compares the two flow chart versions to help you choose the right one for your manuscript.

---

## 📊 Quick Comparison

| Feature | **Simplified** (Main Text) | **Detailed** (Supplemental) |
|---------|---------------------------|----------------------------|
| **Purpose** | High-level workflow overview | Complete step-by-step process |
| **File** | `flowchart_main_text_simplified.mmd` | `flowchart_supplemental_detailed.mmd` |
| **Phases Shown** | 6 major steps | 11 phases with sub-steps |
| **Decision Points** | 3 key decisions | 10 decision points |
| **Complexity** | Simple, easy to follow | Comprehensive, detailed |
| **Target Audience** | General readers, reviewers | Methods enthusiasts, reproducibility |
| **Recommended Use** | Main text Figure 1 or 2 | Supplemental Figure S1 |
| **Typical Dimensions** | ~10 inch wide × ~13 inch tall | ~12 inch wide × ~30 inch tall |
| **File Size (SVG)** | ~50-100 KB | ~150-300 KB |
| **File Size (PNG 300 DPI)** | ~500 KB - 1 MB | ~2-4 MB |

---

## 🎯 Simplified Version (Main Text Figure)

### **Use this version when:**
- ✅ Submitting to journal main text
- ✅ Presenting workflow in talks/posters
- ✅ Target audience needs high-level understanding
- ✅ Space is limited (1-page figure)
- ✅ Want to emphasize key processing steps

### **What's Included:**

```
1. Quality Control & Filtering
   ├─ Orientation normalization
   ├─ ORF validation
   ├─ Length & N content filters
   └─ Coordinate quality marking

2. Haplotype Discovery (ESV Approach)
   ├─ MAFFT alignment
   ├─ Core region extraction
   ├─ 100% identity grouping
   └─ Singleton error filtering

3. Haplotype Assignment & Taxonomy
   ├─ Direct ESV mapping
   ├─ Majority vote species
   └─ Build species labels

4. Geographic Analysis
   ├─ Geographic quality check (DECISION)
   └─ Ocean basin assignment (GOaS)

5. Phylogenetic Analysis (OPTIONAL)
   ├─ Build tree decision (DECISION)
   ├─ MAFFT + FastTree GTR+Γ
   └─ Divergence matrix

6. Visualization & Reports
   ├─ Distribution maps
   ├─ Basin abundance charts
   └─ Interactive HTML
```

### **Sample Counts Shown:**
- Total input samples (N)
- QC-passed samples (M, ~85%)
- Haplotypes identified (K, after error filtering)
- Geographic quality samples (X with basins, Y as Unknown)

### **Decision Points:**
1. **QC Passed?** → Yes (continue) / No (excluded)
2. **Geographic Quality Coords?** → High (basin assignment) / Low (marked Unknown)
3. **Build Tree?** → Yes (phylogenetic analysis) / No (skip to visualization)

---

## 🔬 Detailed Version (Supplemental Figure)

### **Use this version when:**
- ✅ Including in supplemental materials
- ✅ Reviewers request detailed methods
- ✅ Want to document every processing step
- ✅ Emphasizing reproducibility
- ✅ Technical audience needs full workflow

### **What's Included:**

```
PHASE 1: Data Loading & Quality Marking
├─ 1.1: Parse BOLD TSV (validate required columns)
├─ 1.2: Mark Coordinate Quality (flags: centroid, missing, invalid)
└─ 1.3: Assign Ocean Basins (GOaS spatial join OR --no-geo)

PHASE 2: Pre-processing QC
├─ 2.1: Extract Sequences (processid → sequence dict)
├─ 2.2: Orientation Normalization (try forward + rev comp, select best ORF)
├─ 2.3: ORF Validation (coverage ≥50%, internal stops ≤2)
├─ 2.4: Dynamic QC Filtering (length, N content, median-based)
└─ 2.5: ORF Requirement Decision (require_valid_orf flag)

PHASE 3: Haplotype Discovery
├─ 3.1: Multiple Sequence Alignment (MAFFT auto)
├─ 3.2: Gap Masking (columns >50% gaps)
├─ 3.3: Core Region Extraction (≥80% coverage, ≥150bp)
├─ 3.4: Identify ESVs (group by 100% identity in core)
├─ 3.5: Singleton Error Filter (remove ≤0.5% divergent)
└─ 3.6: Suspect Flagging Decision (flag >5% divergent singletons)

PHASE 4: Haplotype Assignment
└─ 4.1: Direct ESV Mapping (exact assignments, no ties)

PHASE 5: Taxonomy Assignment
├─ 5.1: Metadata-Based Taxonomy (majority vote ≥70%)
├─ 5.2: Sequence-Based Taxonomy Decision (BLAST/VSEARCH)
│   ├─ BLAST/VSEARCH vs reference (≥98.5% identity, ≥90% coverage)
│   └─ Skip sequence taxonomy
└─ 5.3: Conflict Resolution (integrate metadata + sequence, build labels)

PHASE 6: Post-assignment QC
├─ 6.1: Contamination Detection (multi-species haplotypes, purity scores)
└─ 6.2: Quality Reports (contamination heatmap, purity distribution)

PHASE 6.5: Species-Level Aggregation
├─ 6.5.1: Aggregate by Species (filter to species-level ≥70%)
├─ 6.5.2: Diversity Metrics (richness, evenness)
├─ 6.5.3: Geographic Summary (samples per species per basin)
└─ 6.5.4: Species-Faceted Subsets (one file per species if ≥2 haplotypes)

PHASE 7: Geographic Analysis Enhancement
└─ 7.1: Geographic Quality Assessment (missing data report, precision, summaries)

PHASE 8: Phylogenetic Analysis (OPTIONAL)
├─ Build Tree Decision
├─ 8.1: Filter Haplotypes (length ≥150bp, cluster size ≥3)
├─ 8.2: Align Filtered Haplotypes (MAFFT)
├─ 8.3: Trim Alignment Decision (trimAl gappyout OR full alignment)
├─ 8.4: Build Tree (FastTree GTR+Gamma, bootstrap)
├─ 8.5: Tree Rooting Decision (outgroup OR midpoint)
├─ 8.6: Relabel Tree (species h# labels)
└─ 8.7: MSA Visualization Decision (phylogeny-ordered plots)

PHASE 9: Divergence Analysis
├─ 9.1: Pairwise p-Distance Matrix (gaps/Ns ignored)
├─ 9.2: Barcoding Gap Analysis (within vs between-species)
├─ 9.3: Species-Level Divergence (mean distance between pairs)
└─ 9.4: Within-Species Divergence (per-species matrices)

PHASE 10: Visualization
├─ 10.1: Haplotype-Level Plots (maps, basin abundance)
├─ 10.2: Species-Level Plots (species maps, species abundance)
└─ 10.3: Species-Faceted Haplotype Plots (per-species maps + charts)

PHASE 11: Report Generation (OPTIONAL)
├─ No Report Decision
└─ 11.1: Interactive HTML Report (dashboard, parameters, plots, methods, tables)
```

### **Parameters Shown:**
- All adjustable thresholds with values
- All decision points with options
- Genetic code (2 = vertebrate mito)
- Distance formula (1 - matches/informative_sites)
- Trimming method (gappyout)
- Resolution specifications (300 DPI, PNG/PDF/SVG)

### **Decision Points (10 total):**
1. GOaS Available? (Phase 1)
2. require_valid_orf? (Phase 2)
3. flag_suspect_haplotypes? (Phase 3)
4. Sequence-based taxonomy available? (Phase 5)
5. build_tree enabled? (Phase 8)
6. trim_alignment? (Phase 8)
7. Outgroup specified? (Phase 8)
8. MSA viz enabled? (Phase 8)
9. no_report flag? (Phase 11)
10. Multiple geographic quality decisions (Phase 4)

---

## 🎨 Visual Differences

### Simplified:
- **Layout:** Vertical flow, single column
- **Boxes:** Large, grouped by major phase
- **Labels:** High-level descriptions
- **Colors:** 5 main colors (input, process, decision, exclude, output)
- **Annotations:** Sample counts at key stages (N, M, K)
- **Space:** Compact, fits on one page

### Detailed:
- **Layout:** Vertical flow, wider boxes for sub-steps
- **Boxes:** Smaller, one per sub-step (60+ boxes)
- **Labels:** Detailed step descriptions with parameters
- **Colors:** 6 colors (adds "skip" style for optional paths)
- **Annotations:** All parameters, sample counts, percentages
- **Space:** Multi-page tall (vertical scroll)

---

## 📝 Figure Captions (Templates)

### Simplified Version:

> **Figure 1. BOLDGenotyper workflow overview.** COI sequences from the Barcode of Life Database (BOLD) undergo quality control including orientation normalization, open reading frame (ORF) validation, and dynamic filtering (Phase 1). Haplotypes are identified using an exact sequence variant (ESV) approach: sequences are aligned, a core shared region is extracted, and sequences with 100% identity are grouped, followed by singleton error filtering to remove sequences ≤0.5% divergent (Phase 2). Samples are assigned to haplotypes via direct ESV mapping, and taxonomy is determined by majority vote (≥70% threshold) (Phase 3). Coordinates are quality-checked, and samples with precise GPS coordinates are assigned to ocean basins via spatial join with the Global Oceans and Seas (GOaS) dataset; samples with centroid or missing coordinates are marked as "Unknown" but retained for genotyping (Phase 4). Optional phylogenetic analysis constructs maximum-likelihood trees using FastTree with the GTR+Γ model and calculates pairwise divergence matrices (Phase 5). Visualizations include global distribution maps, ocean basin abundance charts (relative and total), species-level plots, and species-faceted haplotype maps (Phase 6). Output files include annotated datasets, haplotype sequences, phylogenetic trees, and interactive HTML reports. Blue boxes: input data; green boxes: processing steps; yellow diamonds: decision points; coral dashed boxes: excluded data; purple boxes: output data. Sample counts are provided at key stages (N = total input, M = QC-passed, K = haplotypes identified).

### Detailed Version:

> **Figure S1. Complete BOLDGenotyper workflow with all processing steps and decision points.** This supplemental figure provides a comprehensive view of the BOLDGenotyper pipeline, detailing all 11 phases, sub-steps, parameters, and decision points. **Phase 1** loads BOLD TSV data, marks coordinate quality (centroid, missing, invalid), and assigns ocean basins using GOaS spatial data. **Phase 2** applies pre-processing quality control: orientation normalization selects the best ORF orientation, ORF validation checks coverage and internal stops, and dynamic filtering applies length and N content thresholds. **Phase 3** discovers haplotypes using the ESV approach: sequences are aligned with MAFFT, a core region covering ≥80% of sequences is extracted, exact sequence variants (ESVs) are identified by 100% identity grouping, and singleton errors ≤0.5% divergent are filtered. **Phase 4** assigns samples to haplotypes via direct ESV mapping. **Phase 5** assigns taxonomy using metadata-based majority voting (≥70%) with optional sequence-based validation (BLAST/VSEARCH at ≥98.5% identity). **Phase 6** performs post-assignment quality control, detecting contamination via multi-species haplotypes and purity scores. **Phase 6.5** aggregates samples by species, calculates diversity metrics (haplotype richness, evenness), and generates species-faceted subsets. **Phase 7** enhances geographic analysis with missing data reports and basin summaries. **Phase 8** (optional) builds phylogenetic trees using MAFFT alignment, trimAl trimming, and FastTree GTR+Γ model, with optional outgroup rooting and MSA visualization. **Phase 9** calculates pairwise divergence matrices (p-distance), barcoding gap analysis, species-level divergence, and within-species haplotype divergence. **Phase 10** generates visualizations including haplotype and species distribution maps, ocean basin abundance charts, and species-faceted plots. **Phase 11** (optional) creates an interactive HTML report with summary dashboard, pipeline parameters, and downloadable data. All adjustable parameters are shown with default values; decision points (yellow diamonds) indicate user-configurable options. Color coding as in Figure 1. See main text Methods for detailed parameter descriptions and BOLDGenotyper_Workflow_Documentation.md for complete specifications.

---

## 🚀 Quick Usage Guide

### Option 1: Render Online (No Installation)
1. Visit https://mermaid.live
2. Copy contents of `.mmd` file
3. Paste into editor
4. Download as SVG or PNG

### Option 2: Render Locally (Best Quality)
```bash
# Install Mermaid CLI once
npm install -g @mermaid-js/mermaid-cli

# Navigate to docs directory
cd /Users/stesmith/Documents/depredation/boldgenotyper/docs/

# Run automated rendering script
bash render_flowcharts.sh

# Output:
#   flowchart_main_text_simplified.svg
#   flowchart_main_text_simplified.png
#   flowchart_supplemental_detailed.svg
#   flowchart_supplemental_detailed.png
```

---

## 📐 Recommended Journal Formatting

| Journal Type | Recommended Version | Format | Notes |
|--------------|---------------------|--------|-------|
| **High-impact general** (Nature, Science) | Simplified | SVG or PNG | Clear, concise for broad audience |
| **Ecology/Evolution** (Mol Ecol, Evol Appl) | Simplified main + Detailed supp | Both | Methods detail appreciated |
| **Bioinformatics** (Bioinformatics, BMC Bioinf) | Detailed main | SVG | Technical audience expects detail |
| **Methods journals** (Mol Ecol Resour, MethodsX) | Both required | SVG + PNG | Workflow documentation critical |
| **Marine biology** (Mar Biol, Mar Ecol Prog Ser) | Simplified | PNG or PDF | Standard workflow diagram |

---

## 💡 Tips for Customization

### Color Palette:
Both diagrams use a **colorblind-friendly palette** tested with Color Oracle.

If your journal requires grayscale:
- Input (blue) → Light gray with thick border
- Process (green) → Medium gray
- Decision (yellow) → White with thick border
- Exclude (coral) → Dark gray with dashed border
- Output (purple) → White with double border

### Box Sizes:
- **Simplified:** Boxes auto-size based on text
- **Detailed:** May need manual adjustment in Inkscape/Illustrator for very long parameter descriptions

### Additional Annotations:
You can add:
- Sample count annotations (e.g., "N=1000 → M=850")
- Parameter value callouts (e.g., "Default: 70%")
- Reference citations (e.g., "GTR+Γ [Price et al. 2010]")

---

## ✅ Recommendation Summary

**For most manuscripts:**
- **Main text:** Use simplified version
- **Supplemental:** Use detailed version
- **Format:** SVG for both (universal, scalable, small file size)

**If journal requires single format:**
- Use simplified version only
- Include link to GitHub repository for full workflow documentation

**If page limits are strict:**
- Combine simplified workflow + parameter table (from documentation)
- Reference detailed flow chart in online repository

---

**Questions?** See `FLOWCHART_RENDERING_GUIDE.md` for complete rendering instructions or open an issue at https://github.com/SymbioSeas/BOLDGenotyper/issues
