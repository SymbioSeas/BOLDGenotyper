# BOLDGenotyper Manuscript Preparation Plan
# Comprehensive Strategy for Publication Readiness

**Date**: 2026-01-28
**Purpose**: Strategic plan to prepare BOLDGenotyper for peer-reviewed publication
**Target Journals**: Molecular Ecology Resources (primary), Bioinformatics, BMC Bioinformatics
**Goal**: Publication-ready tool with compelling case studies demonstrating value

---

## Executive Summary

This plan addresses the critical question: **What scope should BOLDGenotyper emphasize for publication, and what work is needed to get there?**

### Recommended Positioning

**Primary Scope**: COI-based haplotype identification and phylogeographic analysis pipeline with **universal applicability** (marine, freshwater, terrestrial) but **demonstrated value** through marine case studies.

**Rationale**:
1. Your existing code already supports universal use (custom shapefiles, --no-geo flag)
2. Marine examples showcase unique value (ocean basin phylogeography)
3. Positioning as "marine-only" would artificially limit adoption and citations
4. Freshwater/terrestrial examples demonstrate generalizability without requiring equal development effort

### Core Value Proposition

**"BOLDGenotyper is the first automated pipeline integrating COI sequence clustering, haplotype assignment, phylogeographic analysis, and population genetics export into a single reproducible workflow."**

What sets it apart:
- **Integration** over fragmentation (8 manual steps → 1 command)
- **Flexibility** across habitats (ocean basins, drainage basins, ecoregions)
- **Quality control** via multi-level comparative analysis
- **Reproducibility** through comprehensive documentation and testing

---

## Part 1: Scope Definition and Positioning

### 1.1 What BOLDGenotyper CAN Do (Emphasize These)

#### Core Capabilities (Universal, All Organisms)
✅ **Automated haplotype discovery** from COI sequences
✅ **Edit distance-based assignment** with tie detection
✅ **Parameter optimization** via data-driven sweep analysis
✅ **Quality control** via species vs. family comparative analysis
✅ **Phylogenetic analysis** with FastTree integration
✅ **Population genetics exports** (Arlequin, POPART, DNASP)
✅ **Interactive HTML reports** with publication methods text
✅ **Metadata analysis** (temporal, categorical, geographic coverage)

#### Geographic Capabilities (Habitat-Dependent)
✅ **Ocean basin classification** (marine organisms via GOaS shapefile)
✅ **Custom geographic regions** (any shapefile: drainage basins, ecoregions, protected areas)
✅ **No-geography mode** (pure haplotyping without spatial component)

### 1.2 What BOLDGenotyper CANNOT Do (Acknowledge Limitations)

❌ **Multi-locus analysis** - Currently COI-only (future: concatenated alignments)
❌ **Haplotype networks** - Exports to PopART but doesn't generate TCS/MJ networks internally
❌ **Advanced population genetics** - Data preparation only; Arlequin/DnaSP perform Fst, AMOVA
❌ **Species delimitation** - Intraspecific focus; not a BOLD BIN replacement
❌ **Real-time sequencing** - Batch processing only (future: streaming for Nanopore)
❌ **Metabarcoding** - Single-species focus, not community ecology

### 1.3 Positioning Strategy: Three Modules, One Pipeline

**Frame BOLDGenotyper as three integrated modules with universal applicability:**

| Module | Function | Applicability | Unique Value |
|--------|----------|---------------|--------------|
| **Haplotyping** | Sequence clustering, consensus generation, assignment | Universal | Edit distance ties, parameter sweep, quality filtering |
| **Metadata Analysis** | Temporal/categorical patterns, data coverage | Universal | Automated visualization of sampling bias, geographic/temporal gaps |
| **Phylogeography** | Geographic classification, basin-specific analyses | Habitat-specific | Ocean basins (marine), custom shapefiles (any), or skipped |

**Key Message**: "All three modules work together, but phylogeography is *optional* and *customizable* based on organism and research question."

### 1.4 Marine vs. Terrestrial: The Strategy

**Primary Framing** (60% emphasis):
- "Designed for phylogeographic analysis with special support for marine organisms via ocean basin classification"
- Case studies: 3 marine examples (Sphyrnidae, Carcharhinus, Rhinopristiformes)
- Showcase ocean basin phylogeography (unique value proposition)

**Generalizability Demonstration** (40% emphasis):
- "Applicable to any organism with COI data via custom shapefiles or geography-free mode"
- Case studies: 1 freshwater (Salmonidae), 1 terrestrial (if available)
- Document custom shapefile workflow thoroughly

**Benefits of This Approach**:
1. Doesn't oversell terrestrial capabilities we haven't extensively tested
2. Doesn't undersell by calling it "marine-only" when it's not
3. Invites adoption across communities while being honest about development focus
4. Positions ocean basin phylogeography as competitive advantage, not limitation

---

## Part 2: Critical Gaps and Requirements for Publication

### 2.1 Essential Prerequisites (Must Complete Before Submission)

#### **Gap 1: Robust Benchmarking and Validation**

**Problem**: No quantitative comparison to manual workflows or existing tools.

**Required Work**:
1. **Runtime benchmarking** across dataset sizes (100, 500, 1000, 5000, 10000 sequences)
   - Measure total pipeline time
   - Break down by phase (QC, alignment, clustering, assignment, phylogenetics)
   - Compare to estimated manual workflow time (literature-based + expert opinion)

2. **Accuracy validation** of haplotype assignments
   - Use Sphyrnidae dataset with manual curation as "gold standard"
   - Calculate precision/recall for genotype assignments
   - Quantify how many ambiguous/tied assignments flag correctly

3. **Parameter sweep validation**
   - Test on 3-5 diverse datasets (different taxonomic groups, sample sizes)
   - Compare sweep-recommended thresholds to literature-based thresholds
   - Demonstrate elbow detection successfully identifies optimal values

4. **Comparative analysis validation**
   - Quantify contamination detection rate (true positives vs. false positives)
   - Use known mislabeled samples if available
   - Document reassignment success rate

**Deliverable**: Benchmarking section for manuscript with figures showing runtime, accuracy, and scalability.

#### **Gap 2: Large-Scale Case Studies**

**Problem**: Need compelling demonstrations of BOLDGenotyper's value beyond basic functionality.

**Required Work**:

**Case Study 1: Sphyrnidae (Family-Level Marine)**
- 1,000+ sequences, multiple species
- Ocean basin phylogeography
- Comparative analysis (species vs. family)
- Parameter sweep for threshold optimization
- **Status**: Data available, needs complete analysis

**Case Study 2: Carcharhinus (Large Marine Dataset)**
- 5,000+ sequences
- Scalability demonstration
- Multiple species within genus
- **Status**: Data available, needs analysis

**Case Study 3: Rhinopristiformes (Conservation Application)**
- Endangered group (wedgefishes, guitarfishes)
- Geographic structure relevant to management units
- **Status**: Data available, needs analysis

**Case Study 4: Salmonidae (Freshwater Example)**
- 10,000+ sequences
- Custom shapefile workflow (drainage basins)
- Demonstrates non-marine applicability
- **Status**: Data available, custom shapefile needed

**Case Study 5: Terrestrial Example (Optional but Valuable)**
- Candidate: Lepidoptera or Coleoptera with good BOLD coverage
- Ecoregion-based analysis
- Demonstrates full generalizability
- **Status**: Need to identify suitable dataset

**Deliverable**: 3-5 case study vignettes in manuscript + supplementary materials with complete results.

#### **Gap 3: Terrestrial/Freshwater Robustness Testing**

**Problem**: Pipeline developed primarily with marine examples; uncertain if edge cases exist for terrestrial organisms.

**Required Testing**:
1. **Coordinate handling** - Test with terrestrial lat/lon (not ocean)
   - Verify land detection doesn't break pipeline
   - Confirm --no-geo mode works seamlessly

2. **Custom shapefile workflow** - End-to-end test with:
   - Freshwater drainage basins (HydroBASINS or similar)
   - Terrestrial ecoregions (WWF ecoregions or similar)
   - Custom conservation areas

3. **Edge case testing**:
   - Species with highly fragmented distributions
   - High-latitude organisms (Arctic/Antarctic)
   - Island endemics

**Deliverable**: Confidence that pipeline is robust for any organism, documented in test suite and tutorials.

#### **Gap 4: Publication-Ready Documentation**

**Problem**: Existing docs are excellent but need final polish and case study integration.

**Required Work**:
1. **Update README** with case study examples
2. **Create 2-3 detailed tutorials**:
   - Tutorial 1: Marine organism with ocean basins (Sphyrnidae)
   - Tutorial 2: Freshwater organism with custom shapefile (Salmonidae)
   - Tutorial 3: Terrestrial organism or geography-free mode
3. **Methods text generation** - Auto-generate publication methods for each case study
4. **Troubleshooting expansion** - Add sections based on terrestrial testing

**Deliverable**: Documentation that enables independent use by reviewers and early adopters.

### 2.2 Important Enhancements (Should Complete If Possible)

#### **Enhancement 1: Haplotype Network Export**

**Rationale**: PopART is standard for COI haplotype networks; seamless integration would increase adoption.

**Implementation**:
- Generate PopART-compatible frequency tables (already done via popgen export)
- Create trait files for geographic/temporal groupings
- Document workflow in tutorial

**Effort**: Low (mostly documentation, export code exists)

#### **Enhancement 2: Automated Methods Text Refinement**

**Rationale**: Publication-ready methods text is major selling point.

**Implementation**:
- Refine methods templates based on case studies
- Add species-specific details (sample sizes, date ranges)
- Generate complete Methods + Results paragraph

**Effort**: Medium (template improvements)

#### **Enhancement 3: Performance Optimization for Large Datasets**

**Rationale**: Demonstrate BOLDGenotyper can handle family-level analyses (10,000+ sequences).

**Implementation**:
- Profile bottlenecks in Carcharhinus/Salmonidae runs
- Optimize MAFFT parameters for large datasets
- Implement chunking for edit distance calculations if needed

**Effort**: Medium-High (profiling + optimization)

### 2.3 Future Directions (Explicitly Defer to Future Work)

Document these as planned enhancements, not current limitations:

🔮 **Multi-locus support** - Concatenated alignments for nuclear + mitochondrial
🔮 **Haplotype network generation** - Internal TCS/median-joining implementation
🔮 **Population genetics statistics** - Fst, AMOVA, diversity indices calculated internally
🔮 **Web interface** - Galaxy integration or standalone web app
🔮 **Machine learning threshold optimization** - Adaptive parameter selection
🔮 **Real-time analysis** - Integration with Oxford Nanopore streaming

---

## Part 3: Manuscript Structure and Content

### 3.1 Target Journal: Molecular Ecology Resources

**Why This Journal**:
- Perfect audience (population geneticists, phylogeographers)
- Software/methods focused
- High impact in field (IF: 5.8)
- Published similar tools (PopART, STRUCTURE)

**Format**: Application Note (4-6 pages) or Full Article (if extensive case studies)

### 3.2 Manuscript Outline

**Title**: "BOLDGenotyper: An automated pipeline for COI haplotype identification and phylogeographic analysis"

**Abstract** (250 words):
- Problem: Gap between species ID (BOLD) and population genetics (manual workflows)
- Solution: Automated pipeline integrating clustering, assignment, geography, export
- Innovation: Data-driven optimization, flexible geography, multi-level QC
- Validation: 3-5 case studies across marine/freshwater demonstrating utility
- Availability: Open source, conda/pip installable

**Introduction** (~1200 words):
1. DNA barcoding success and BOLD database growth
2. Phylogeographic applications of COI data
3. Current workflow: fragmented, manual, error-prone (cite 5-7 tools each doing one step)
4. Need for integrated pipeline
5. BOLDGenotyper overview and positioning

**Methods** (~2000 words):
1. Pipeline architecture (flowchart)
2. Module 1: Haplotyping (dereplication, clustering, assignment)
3. Module 2: Metadata analysis (temporal, categorical, coverage)
4. Module 3: Phylogeography (ocean basins, custom shapefiles, --no-geo)
5. Parameter optimization (sweep algorithm)
6. Quality control (comparative analysis)
7. Exports and visualization
8. Implementation (Python, dependencies, testing)

**Results** (~1500 words):
1. **Case Study 1**: Sphyrnidae ocean basin phylogeography
2. **Case Study 2**: Carcharhinus scalability demonstration
3. **Case Study 3**: Salmonidae custom shapefile workflow
4. **Benchmarking**: Runtime, accuracy, comparison to manual
5. **Validation**: Parameter sweep and QC effectiveness

**Discussion** (~1200 words):
1. Fills critical gap in phylogeography toolkit
2. Advantages over fragmented workflows (time savings, reproducibility)
3. Universal applicability (marine to terrestrial)
4. Quality control innovation (contamination detection)
5. Limitations and future directions
6. Potential for broad adoption

**Availability**:
- GitHub, Zenodo DOI, conda/pip instructions
- Comprehensive documentation and tutorials
- Example datasets

**Supplementary Material**:
- Tutorial 1: Marine workflow
- Tutorial 2: Freshwater workflow with custom shapefile
- Tutorial 3: Comparative analysis for QC
- Extended case study results
- Benchmark data and scripts

### 3.3 Key Figures

**Figure 1: Pipeline Overview**
- Flowchart: Input → Modules → Outputs
- Highlight decision points (--build-tree, --no-geo, --export-format)
- Show integration with external tools (MAFFT, FastTree, PopART)

**Figure 2: Sphyrnidae Case Study**
- Panel A: Global distribution map colored by ocean basin
- Panel B: Haplotype distribution across basins (bar chart)
- Panel C: Phylogenetic tree with ocean basin annotations
- Panel D: Parameter sweep results with optimal threshold

**Figure 3: Comparative Analysis Workflow**
- Species vs. family genotype crosswalk
- Example contamination detection
- Reassignment recommendations

**Figure 4: Performance Benchmarks**
- Panel A: Runtime vs. dataset size (100-10,000 sequences)
- Panel B: Assignment accuracy (precision/recall)
- Panel C: Comparison to manual workflow time

**Figure 5: Universal Applicability**
- Panel A: Marine (ocean basins)
- Panel B: Freshwater (drainage basins)
- Panel C: Terrestrial (ecoregions) or no-geography mode
- Same pipeline, different shapefiles/modes

---

## Part 4: Implementation Timeline

### Phase 1: Foundation (Weeks 1-2)

**Week 1: Benchmarking and Validation**
- [ ] Runtime benchmarking across dataset sizes
- [ ] Accuracy validation with Sphyrnidae manual curation
- [ ] Parameter sweep validation on multiple datasets
- [ ] Comparative analysis validation (contamination detection)

**Week 2: Terrestrial/Freshwater Testing**
- [ ] Test coordinate handling with terrestrial species
- [ ] Create custom shapefile workflow (drainage basins for Salmonidae)
- [ ] Test ecoregion-based analysis (terrestrial example)
- [ ] Document edge cases and add tests

### Phase 2: Case Studies (Weeks 3-5)

**Week 3: Marine Case Studies**
- [ ] Complete Sphyrnidae analysis (family-level, ocean basins)
- [ ] Complete Carcharhinus analysis (scalability, large dataset)
- [ ] Complete Rhinopristiformes analysis (conservation application)

**Week 4: Non-Marine Case Studies**
- [ ] Complete Salmonidae analysis (freshwater, drainage basins)
- [ ] Identify and analyze terrestrial example (optional but valuable)
- [ ] Generate publication-quality figures for all case studies

**Week 5: Case Study Synthesis**
- [ ] Write Methods text for each case study
- [ ] Write Results vignettes
- [ ] Create supplementary materials
- [ ] Generate all figures

### Phase 3: Documentation and Repository Polish (Week 6)

- [ ] Update README with case study examples
- [ ] Create Tutorial 1: Marine workflow (Sphyrnidae)
- [ ] Create Tutorial 2: Freshwater workflow (Salmonidae)
- [ ] Create Tutorial 3: Terrestrial/no-geo workflow
- [ ] Update troubleshooting based on testing
- [ ] Archive outdated planning documents
- [ ] Clean up test data (keep examples, remove old runs)

### Phase 4: Manuscript Writing (Weeks 7-9)

**Week 7: Methods and Results**
- [ ] Write Methods section (pipeline description)
- [ ] Write Results section (case studies)
- [ ] Create all main figures

**Week 8: Introduction and Discussion**
- [ ] Write Introduction (position in literature)
- [ ] Write Discussion (implications, limitations, future)
- [ ] Create supplementary materials

**Week 9: Finalization**
- [ ] Write Abstract
- [ ] Compile References
- [ ] Format for target journal
- [ ] Internal review and revision

### Phase 5: Pre-Submission (Week 10)

- [ ] Submit to Zenodo for DOI
- [ ] Test installation on clean systems (Mac, Linux, Windows)
- [ ] Verify all documentation links
- [ ] Co-author review
- [ ] Final revisions
- [ ] Submit to journal

**Total Timeline**: 10 weeks to submission

---

## Part 5: Scope Decisions and Trade-offs

### 5.1 What to Include in v0.1.0 for Publication

**Must Have** (In Current Code, Just Need Testing/Docs):
✅ Core haplotyping workflow
✅ Ocean basin classification
✅ Custom shapefile support
✅ Parameter sweep
✅ Comparative analysis
✅ Population genetics exports
✅ Interactive HTML reports

**Should Have** (Minor Additions):
🔧 Enhanced methods text generation
🔧 Tutorial notebooks for each use case
🔧 Performance optimizations for large datasets

**Won't Have** (Defer to v0.2.0+):
❌ Multi-locus support
❌ Internal haplotype network generation
❌ Fst/AMOVA calculations
❌ Web interface
❌ Machine learning threshold selection

### 5.2 Marine vs. Terrestrial: Final Decision

**Decision**: Position as **universal COI phylogeography tool** with **demonstrated value in marine systems** via ocean basin analyses.

**Manuscript Language**:
- Introduction: "BOLDGenotyper is applicable to any organism with COI barcode data..."
- Methods: "The pipeline includes optional geographic analysis with built-in support for marine organisms (ocean basins via GOaS) and customizable support for freshwater (drainage basins) and terrestrial (ecoregions) organisms via user-provided shapefiles."
- Results: "We demonstrate BOLDGenotyper's utility with marine case studies... and generalizability with freshwater example..."
- Discussion: "While our case studies focus on marine organisms, the modular design and custom shapefile support enable application to any habitat."

**Why This Works**:
1. Honest about development focus (marine)
2. Doesn't artificially limit scope
3. Demonstrates generalizability with real examples
4. Invites broader adoption
5. Sets up future papers on terrestrial applications by other groups (citations!)

### 5.3 What Can Fail and Pipeline Still Provides Value?

**Scenario 1: User has terrestrial organism, no shapefile**
- ✅ Pipeline runs with `--no-geo` flag
- ✅ Still get: haplotyping, assignment, phylogenetics, exports, metadata analysis
- ✅ Value: Complete genotyping + QC + exports for population genetics

**Scenario 2: Ocean basins don't make biological sense**
- ✅ User provides custom shapefile (e.g., Arctic vs. Antarctic)
- ✅ Pipeline adapts geographic analysis to user's regions
- ✅ Value: Flexible phylogeography tool, not rigid ocean-only

**Scenario 3: User only wants haplotype IDs, no geography**
- ✅ `--no-geo` mode provides exactly this
- ✅ Value: Fast, automated genotyping with tie detection and QC

**Key Insight**: The three-module design (haplotyping, metadata, phylogeography) means failure of one doesn't break the others. This is a strength to emphasize!

---

## Part 6: Success Metrics

### 6.1 Publication Success

**Acceptance Criteria**:
- [ ] Manuscript accepted in Molecular Ecology Resources (primary) or equivalent
- [ ] Zenodo DOI assigned
- [ ] Software version 0.1.0 released
- [ ] Installation tested on Mac/Linux/Windows

**Timeline**: Submission by March 2026, acceptance by June 2026

### 6.2 Community Adoption Metrics (Year 1)

**Software Usage**:
- Target: 50+ GitHub stars
- Target: 100+ conda downloads/month
- Target: 10+ independent users (based on issues/questions)

**Citations**:
- Target: 10-20 citations in first year
- Target: Use in 3+ independent studies

### 6.3 Documentation Success

**User Experience**:
- [ ] New user can install and run basic analysis in <15 minutes
- [ ] All use cases documented with tutorials
- [ ] Troubleshooting guide covers 90% of common issues

---

## Part 7: Specific To-Do Items Before Submission

### 7.1 Code Changes

**Priority 1 (Essential)**:
- [ ] Fix any filename inconsistencies (spaces vs. underscores)
- [ ] Ensure all CLI commands generate consistent outputs
- [ ] Test --no-geo mode thoroughly
- [ ] Test custom shapefile workflow end-to-end
- [ ] Add logging for parameter choices (for reproducibility)

**Priority 2 (Important)**:
- [ ] Optimize runtime for >5000 sequences (profile bottlenecks)
- [ ] Add progress bars for long-running steps
- [ ] Improve error messages for common mistakes
- [ ] Generate benchmark scripts for reproducibility

### 7.2 Documentation Changes

**Priority 1 (Essential)**:
- [ ] README: Add case study examples in Quick Start
- [ ] Create Marine Tutorial (Sphyrnidae)
- [ ] Create Freshwater Tutorial (Salmonidae + custom shapefile)
- [ ] Update CHANGELOG for v0.1.0 release
- [ ] Create CITATION.cff file

**Priority 2 (Important)**:
- [ ] API_REFERENCE.md: Document programmatic usage
- [ ] Create comparison table vs. other tools
- [ ] Add FAQ section to README
- [ ] Create video tutorial (optional but valuable)

### 7.3 Testing and Validation

**Priority 1 (Essential)**:
- [ ] Run full test suite and achieve >80% coverage
- [ ] Test on fresh conda environments (Mac/Linux)
- [ ] Test each tutorial start-to-finish
- [ ] Benchmark runtime on standardized datasets

**Priority 2 (Important)**:
- [ ] Test on Windows (if possible)
- [ ] Stress test with 10,000+ sequences
- [ ] Test with problematic datasets (many missing values, poor coordinates, etc.)

### 7.4 Case Study Checklist

For **each** case study, generate:
- [ ] Complete annotated output CSV
- [ ] All publication figures (PNG 300dpi + PDF)
- [ ] Interactive HTML report
- [ ] Methods text paragraph
- [ ] Supplementary tables (sample metadata, haplotype assignments, geographic distribution)
- [ ] Interpretation notes (biological findings)

---

## Part 8: Risk Assessment and Mitigation

### 8.1 Potential Reviewer Concerns

**Concern 1**: "How is this different from BOLD BINs?"
- **Mitigation**: Clear explanation in Introduction - BINs are species-level, BOLDGenotyper is intraspecific
- **Evidence**: Case studies showing within-species haplotype structure

**Concern 2**: "Why not just use existing tools (VSEARCH + custom scripts)?"
- **Mitigation**: Benchmarking showing time savings, reproducibility advantages
- **Evidence**: Survey of recent phylogeography papers showing fragmented workflows

**Concern 3**: "Limited to marine organisms?"
- **Mitigation**: Freshwater case study + terrestrial example
- **Evidence**: Custom shapefile documentation, --no-geo mode

**Concern 4**: "Not enough novel methodology?"
- **Mitigation**: Emphasize integration, parameter sweep, comparative analysis as innovations
- **Evidence**: No existing tool combines all these features

**Concern 5**: "Insufficient validation?"
- **Mitigation**: Comprehensive benchmarking section
- **Evidence**: Accuracy metrics, multiple case studies, comparison to manual curation

### 8.2 Technical Risks

**Risk 1**: Terrestrial testing reveals major bugs
- **Probability**: Low (code is habitat-agnostic)
- **Impact**: Medium (would delay submission)
- **Mitigation**: Test early in Phase 1

**Risk 2**: Large datasets (>10,000) fail or are too slow
- **Probability**: Medium (Salmonidae has 10,000+)
- **Impact**: Medium (would need optimization)
- **Mitigation**: Profile and optimize in Phase 2

**Risk 3**: Custom shapefile workflow is too complicated for users
- **Probability**: Medium
- **Impact**: Low (doesn't break core functionality)
- **Mitigation**: Excellent tutorial, provide example shapefiles

### 8.3 Timeline Risks

**Risk 1**: Case studies take longer than expected
- **Probability**: High (always do)
- **Impact**: High (could push timeline)
- **Mitigation**: Focus on 3 case studies (Sphyrnidae, Carcharhinus, Salmonidae) as minimum, others optional

**Risk 2**: Co-author review takes longer than 1 week
- **Probability**: Medium
- **Impact**: Medium
- **Mitigation**: Start writing early, involve co-authors in planning

---

## Part 9: Long-Term Vision

### 9.1 How This Publication Sets Up Future Work

**Paper 1** (This Manuscript): **Methods Paper**
- Tool description and validation
- Molecular Ecology Resources or Bioinformatics
- Citation baseline for users

**Paper 2** (Concurrent or Following): **Biological Paper**
- Sphyrnidae phylogeography findings
- Marine Biology, Molecular Ecology, or Conservation Genetics
- Demonstrates scientific value of tool

**Paper 3+** (Future): **Applications**
- Other groups use BOLDGenotyper for their systems
- Citations accumulate
- Tool becomes community standard

### 9.2 Community Building

**Short Term** (Months 1-6):
- Respond quickly to GitHub issues
- Present at conferences (e.g., Evolution, Molecular Ecology conferences)
- Engage with BOLD Systems team

**Medium Term** (Months 6-12):
- Workshop at conferences
- Integration with Galaxy or other platforms
- Build user community via Google Group or Slack

**Long Term** (Years 2-3):
- Contribute to BOLD v6 development?
- Collaborate on multi-locus extensions
- Establish as standard in phylogeography courses

---

## Part 10: Action Plan Summary

### Immediate Actions (Next 2 Weeks)

1. **Decide on scope framing** (marine-focused universal vs. truly universal)
   - Recommendation: Universal with marine demonstrated value

2. **Identify case study datasets** (3 minimum, 5 ideal)
   - Required: Sphyrnidae, Carcharhinus, Salmonidae
   - Optional: Rhinopristiformes, terrestrial example

3. **Begin benchmarking** (runtime, accuracy, scalability)
   - Set up standardized benchmark protocol
   - Run on existing datasets

4. **Test terrestrial/freshwater robustness**
   - Create custom shapefile for Salmonidae
   - Test --no-geo mode thoroughly

### Next Month (Weeks 3-6)

1. **Complete all case studies**
2. **Finish benchmarking and validation**
3. **Create tutorials**
4. **Polish documentation**

### Following Month (Weeks 7-10)

1. **Write manuscript**
2. **Generate all figures**
3. **Internal review**
4. **Submit**

---

## Part 11: Questions for Discussion

Before proceeding with implementation, please consider:

1. **Scope**: Do you agree with "universal tool with marine demonstrated value" positioning, or prefer "marine-first with future terrestrial"?

2. **Case Studies**: Which 3-5 case studies should we prioritize? Do you have terrestrial datasets available?

3. **Target Journal**: Is Molecular Ecology Resources the right target, or prefer Bioinformatics/BMC Bioinformatics?

4. **Timeline**: Is 10 weeks realistic given other commitments, or should we plan for 15-20 weeks?

5. **Co-authors**: Who should be involved in the manuscript? Need to coordinate with them early.

6. **Funding acknowledgments**: Any grants to acknowledge?

7. **Companion biological paper**: Should Sphyrnidae biological findings be in this paper or separate?

---

## Conclusion

BOLDGenotyper is in excellent shape for publication. The code is mature, well-documented, and tested. The main gaps are:

1. **Benchmarking** - Quantitative validation needed
2. **Case studies** - Demonstrate value across use cases
3. **Terrestrial testing** - Ensure robustness beyond marine
4. **Manuscript writing** - Tell the story compellingly

With focused effort over 10-15 weeks, BOLDGenotyper will be ready for a high-impact publication that establishes it as the standard tool for COI phylogeography.

The recommended positioning as a **universal COI phylogeography pipeline with demonstrated marine applications** balances honesty about development focus with openness to broad adoption. This approach maximizes potential impact and citations while setting realistic expectations.

---

**Next Step**: Review this plan, make decisions on key questions, then proceed with implementation phase by phase.

**Document Author**: Claude (Sonnet 4.5)
**Date**: 2026-01-28
**For**: Steph Smith, UNC Institute of Marine Sciences
