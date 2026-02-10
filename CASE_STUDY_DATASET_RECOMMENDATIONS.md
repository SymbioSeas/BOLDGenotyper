# BOLDGenotyper Case Study Dataset Recommendations

**Purpose**: Select optimal organisms for demonstrating BOLDGenotyper capabilities
**Target**: 5 diverse case studies spanning marine, freshwater, and terrestrial habitats
**Date**: 2026-01-28

---

## Selection Criteria

For each case study, we need datasets that:

1. **Have sufficient BOLD coverage** (>500 sequences minimum, >1,000 ideal)
2. **Span appropriate geographic ranges** (demonstrate phylogeographic patterns)
3. **Have published literature** (for validation and comparison)
4. **Represent different taxonomic groups** (demonstrate generalizability)
5. **Test different BOLDGenotyper features** (each case study highlights specific capabilities)

---

## Recommended Datasets

### Case Study 1: Sphyrnidae (Marine, Family-Level) ✓ CONFIRMED

**Organism**: Sphyrnidae (hammerhead sharks)
**Primary Species**: *Sphyrna lewini* (scalloped hammerhead)
**Sample Size**: ~1,400 sequences (family-level)
**Geographic Range**: Cosmopolitan, all major ocean basins

**Purpose**:
- Demonstrate ocean basin phylogeography
- Family-level analysis with multiple species
- Parameter sweep optimization
- Comparative analysis (species vs. family)
- **Companion biological manuscript** submitted to PLOS ONE

**Expected Patterns**:
- Ocean basin-scale structure in *S. lewini*
- Multiple species with different distributions
- Well-studied system with published phylogeography (Duncan et al. 2006)

**BOLDGenotyper Features Highlighted**:
- ✓ Ocean basin classification
- ✓ Multi-species family-level analysis
- ✓ Parameter sweep
- ✓ Comparative analysis for contamination detection
- ✓ Complete workflow demonstration

**BOLD Search**: `Sphyrnidae` (family level, COI marker)

---

### Case Study 2: Large Elasmobranch Dataset (Marine, Scalability)

**Recommended Organisms** (choose one with best coverage):

**Option A: Carcharhinidae (Requiem Sharks)** - RECOMMENDED
- **Sample Size Estimate**: 8,000-10,000 sequences
- **Geographic Range**: Cosmopolitan, multiple ocean basins
- **Species**: Multiple (*Carcharhinus* genus well-represented)
- **Purpose**: Demonstrate scalability, runtime with large datasets
- **Literature**: Extensive phylogeographic literature for multiple species
- **BOLD Search**: `Carcharhinidae` (family level, COI marker)

**Option B: Rajidae (Skates)**
- **Sample Size Estimate**: 3,000-5,000 sequences
- **Geographic Range**: Global, temperate to polar waters
- **Purpose**: Alternative if Carcharhinidae too large or problematic
- **BOLD Search**: `Rajidae` (family level, COI marker)

**BOLDGenotyper Features Highlighted**:
- ✓ Scalability (>5,000 sequences)
- ✓ Runtime benchmarking
- ✓ Memory efficiency
- ✓ Multi-species haplotype assignment
- ✓ Large phylogenetic tree construction

**Recommendation**: **Carcharhinidae** - Likely to have 8,000+ sequences, widely studied, demonstrates pipeline handles very large datasets efficiently.

---

### Case Study 3: Marine Invertebrate (Accuracy Validation)

**Recommended Organisms** (choose one with best published phylogeography):

**Option A: Panulirus (Spiny Lobsters)** - RECOMMENDED
- **Sample Size Estimate**: 1,000-2,000 sequences
- **Geographic Range**: Circumtropical, multiple ocean basins
- **Species**: Multiple species with distinct distributions
- **Published Phylogeography**: YES - extensive literature on *P. argus* (Caribbean), *P. interruptus* (Eastern Pacific)
- **Purpose**: Validate haplotype assignments against published studies
- **BOLD Search**: `Panulirus` (genus level, COI marker)

**Key Literature for Validation**:
- Silberman et al. (1994) - Caribbean spiny lobster phylogeography
- Perez-Enriquez et al. (2001) - Pacific spiny lobster population structure
- Multiple recent studies on population connectivity

**Option B: Mytilus (Mussels)**
- **Sample Size Estimate**: 2,000-3,000 sequences
- **Geographic Range**: Global, all ocean basins
- **Species**: *M. edulis*, *M. galloprovincialis*, *M. trossulus*
- **Published Phylogeography**: YES - extremely well-studied model system
- **Purpose**: Validate against extensive published data
- **BOLD Search**: `Mytilus` (genus level, COI marker)

**Option C: Portunidae (Swimming Crabs)**
- **Sample Size Estimate**: 1,500-2,500 sequences
- **Geographic Range**: Widespread
- **Species**: Multiple genera (*Portunus*, *Callinectes*, *Scylla*)
- **Published Phylogeography**: YES - particularly *Callinectes sapidus*
- **BOLD Search**: `Portunidae` (family level, COI marker)

**BOLDGenotyper Features Highlighted**:
- ✓ Accuracy validation (compare to published haplotypes)
- ✓ Demonstrates agreement with manual analyses
- ✓ Ocean basin phylogeography in invertebrates
- ✓ Benchmark against literature

**Recommendation**: **Panulirus** - Moderate size (~1,500 sequences), excellent published phylogeography for multiple species, demonstrates accuracy validation clearly.

**Validation Protocol**:
1. Run BOLDGenotyper on Panulirus dataset
2. Compare identified haplotypes to published studies
3. Calculate concordance (how many BOLDGenotyper haplotypes match literature?)
4. Quantify assignment accuracy (precision/recall)
5. Document any discrepancies and investigate causes

---

### Case Study 4: Salmonidae (Freshwater, Drainage Basins) ✓ CONFIRMED

**Organism**: Salmonidae (salmon, trout, char)
**Primary Taxa**: *Oncorhynchus* (Pacific salmon/trout), *Salmo* (Atlantic salmon/trout)
**Sample Size**: 10,000+ sequences (family level)
**Geographic Range**: Northern hemisphere freshwater and anadromous

**Purpose**:
- Demonstrate freshwater organism analysis
- Custom shapefile workflow (drainage basins)
- Very large dataset scalability
- Well-studied phylogeography

**Expected Patterns**:
- Drainage basin structure (Columbia River, Fraser River, Great Lakes, European rivers)
- Postglacial colonization patterns
- Anadromous vs. freshwater resident populations

**BOLDGenotyper Features Highlighted**:
- ✓ Custom shapefile integration (HydroSHEDS or HydroBASINS)
- ✓ Freshwater organism analysis
- ✓ Very large dataset (>10,000 sequences)
- ✓ Alternative to ocean basin classification

**Custom Shapefile**:
- **HydroBASINS** (https://www.hydrosheds.org/products/hydrobasins)
  - Global drainage basin boundaries at multiple scales
  - Level 4-6 basins appropriate for Salmonidae
  - Free download, well-documented

- **Alternative**: USGS HUC (Hydrologic Unit Codes) for North America

**BOLD Search**: `Salmonidae` (family level, COI marker)

**Tutorial Value**: This will be the primary example for custom shapefile documentation.

---

### Case Study 5: Terrestrial Arthropod (Terrestrial Demonstration)

**Recommended Organisms**:

**Option A: Pieridae (White and Sulfur Butterflies)** - RECOMMENDED
- **Sample Size Estimate**: 3,000-5,000 sequences
- **Geographic Range**: Global, all continents except Antarctica
- **Species**: Multiple widespread genera (*Pieris*, *Colias*, *Eurema*)
- **Published Phylogeography**: YES - *Pieris rapae* (cabbage white) extremely well-studied
- **Purpose**: Terrestrial demonstration, ecoregion analysis
- **BOLD Coverage**: Excellent - butterflies heavily barcoded
- **BOLD Search**: `Pieridae` (family level, COI marker)

**Why Pieridae**:
- Excellent BOLD coverage (DNA barcoding flagship group)
- Well-studied phylogeography and biogeography
- Global distribution demonstrates terrestrial patterns
- *Pieris rapae* alone has extensive invasion genetics literature
- Manageable size (~3,000 sequences estimated)

**Option B: Formicidae (Ants)**
- **Sample Size Estimate**: Variable, 2,000-4,000 sequences
- **Geographic Range**: Global
- **Published Phylogeography**: YES - multiple species
- **BOLD Coverage**: Good to excellent
- **BOLD Search**: `Formicidae` (family level, COI marker)

**Option C: Carabidae (Ground Beetles)**
- **Sample Size Estimate**: 4,000-6,000 sequences
- **Geographic Range**: Global
- **Published Phylogeography**: YES - extensive
- **BOLD Coverage**: Very good
- **BOLD Search**: `Carabidae` (family level, COI marker)

**Option D: Lycaenidae (Blues, Coppers, Hairstreaks)**
- **Sample Size Estimate**: 2,000-4,000 sequences
- **Geographic Range**: Global
- **Published Phylogeography**: YES - particularly blues
- **BOLD Coverage**: Excellent (butterfly DNA barcoding well-developed)
- **BOLD Search**: `Lycaenidae` (family level, COI marker)

**BOLDGenotyper Features Highlighted**:
- ✓ Terrestrial organism analysis
- ✓ Custom shapefile (WWF terrestrial ecoregions)
- ✓ --no-geo mode option (if ecoregions not biologically relevant)
- ✓ Demonstrates complete generalizability

**Custom Shapefile**:
- **WWF Terrestrial Ecoregions** (https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world)
  - 867 ecoregions worldwide
  - Biologically meaningful boundaries
  - Free download, widely used
  - Appropriate scale for butterfly/insect phylogeography

**Recommendation**: **Pieridae** - Excellent BOLD coverage, well-studied phylogeography, manageable size, clear demonstration of terrestrial applicability.

---

## Alternative Consideration: Marine/Terrestrial Spanning Group

**Question**: Should Case Study 5 be exclusively terrestrial, or include a group spanning marine and terrestrial?

**Exclusively Terrestrial (Recommended)**:

**Pros**:
- Clear habitat separation (marine → freshwater → terrestrial)
- Demonstrates pipeline truly works outside ocean contexts
- Simpler narrative for manuscript
- Stronger demonstration of universality

**Cons**:
- Doesn't showcase flexibility of switching shapefiles within one analysis

**Marine/Terrestrial Spanning Group** (e.g., Grapsidae - shore crabs, or Decapoda subset):

**Pros**:
- Showcases custom shapefile flexibility
- Interesting biological scenario (coastal/terrestrial transition)
- Could use both ecoregions AND ocean basins in one analysis

**Cons**:
- More complex to analyze and present
- Might confuse the narrative
- Requires two shapefiles for one organism
- Less clear demonstration of "terrestrial works"

**Recommendation**: Use **exclusively terrestrial** for Case Study 5 (Pieridae). If you want to showcase marine/terrestrial spanning, consider adding a **6th supplementary example** rather than replacing terrestrial demonstration. The five core case studies should have clear habitat separation.

---

## Summary Table: Recommended Final Dataset List

| # | Organism | Habitat | Size | Purpose | BOLDGenotyper Features | Custom Shapefile |
|---|----------|---------|------|---------|----------------------|------------------|
| 1 | **Sphyrnidae** | Marine | ~1,400 | Ocean basin phylogeography, companion paper | All features, full workflow | GOaS (ocean basins) |
| 2 | **Carcharhinidae** | Marine | 8,000-10,000 | Scalability, runtime | Large dataset handling | GOaS (ocean basins) |
| 3 | **Panulirus** | Marine | ~1,500 | Accuracy validation | Agreement with literature | GOaS (ocean basins) |
| 4 | **Salmonidae** | Freshwater | 10,000+ | Custom shapefile, freshwater | Drainage basin analysis | HydroBASINS |
| 5 | **Pieridae** | Terrestrial | 3,000-5,000 | Terrestrial demonstration | Ecoregion analysis or --no-geo | WWF Ecoregions |

**Habitat Coverage**: 3 marine, 1 freshwater, 1 terrestrial
**Size Range**: 1,400 to 10,000+ sequences (demonstrates scalability)
**Geographic Features**: Ocean basins (3x), drainage basins (1x), ecoregions (1x)
**Total Sequences**: ~24,000-28,000 across all case studies

---

## Data Download Protocol

For each organism, download from BOLD:

### Search Parameters
1. Navigate to https://boldsystems.org/
2. Click "Public Data Portal"
3. Search by taxonomic name (family or genus level)
4. Select "COI-5P" marker (standard barcode region)
5. Download "Combined TSV"

### Required Fields (verify in downloaded files)
**Essential**:
- processid (sample identifier)
- nucleotides (COI sequence)
- species_name
- lat (latitude)
- lon (longitude)

**Highly Recommended**:
- country
- province_state
- coord_source
- identification_method
- bin_uri (BIN assignment)
- institution_storing
- collectors
- collection_date

### Quality Pre-Screening (Before Full Analysis)
Run quick checks on downloaded data:

```bash
# Check sample sizes
wc -l *.tsv

# Check sequence length distribution
cut -f X *.tsv | awk '{print length}' | sort | uniq -c

# Check coordinate availability
cut -f lat,lon *.tsv | grep -v "^$" | wc -l

# Check species diversity
cut -f species_name *.tsv | sort | uniq -c | sort -rn | head -20
```

---

## Expected Dataset Characteristics

After download, verify each dataset meets criteria:

### Sphyrnidae
- **Expect**: >1,000 sequences, 5-8 species, global distribution
- **Key Species**: *Sphyrna lewini* (majority), *S. mokarran*, *S. zygaena*
- **Coordinates**: 60-80% with lat/lon, mixed quality (centroid vs. precise)
- **Sequence Length**: 400-700 bp typical for COI

### Carcharhinidae
- **Expect**: >8,000 sequences, 20+ species, global distribution
- **Key Genera**: *Carcharhinus*, *Prionace*, *Triaenodon*, *Galeocerdo*
- **Coordinates**: 50-70% with lat/lon
- **Sequence Length**: 400-650 bp

### Panulirus
- **Expect**: 1,000-2,000 sequences, 5-8 species, circumtropical
- **Key Species**: *P. argus*, *P. interruptus*, *P. cygnus*, *P. ornatus*
- **Coordinates**: 70-85% with lat/lon (marine inverts better georeferenced)
- **Sequence Length**: 400-700 bp

### Salmonidae
- **Expect**: >10,000 sequences, 30+ species, Northern Hemisphere
- **Key Genera**: *Oncorhynchus*, *Salmo*, *Salvelinus*
- **Coordinates**: 60-80% with lat/lon, freshwater + anadromous
- **Sequence Length**: 400-700 bp

### Pieridae
- **Expect**: 3,000-5,000 sequences, 50+ species, global
- **Key Genera**: *Pieris*, *Colias*, *Eurema*, *Pontia*
- **Coordinates**: 75-90% with lat/lon (butterflies well-documented)
- **Sequence Length**: 400-700 bp (very consistent for Lepidoptera)

---

## Validation Literature (To Gather)

For each case study, compile key phylogeographic papers for comparison:

### Sphyrnidae
- Duncan et al. (2006) - Global *S. lewini* phylogeography
- Daly-Engel et al. (2012) - Sex-biased dispersal
- Your submitted PLOS ONE manuscript

### Carcharhinidae
- Multiple species-specific studies
- Focus on *Carcharhinus* genus reviews

### Panulirus
- Silberman et al. (1994) - *P. argus* Caribbean
- Perez-Enriquez et al. (2001) - *P. interruptus* Pacific
- Recent connectivity studies

### Salmonidae
- Numerous studies - focus on postglacial colonization
- Drainage basin structure papers

### Pieridae
- *Pieris rapae* invasion genetics
- Biogeographic reviews of Pieridae

---

## Next Steps After Download

1. **Download all 5 datasets from BOLD** using protocol above
2. **Verify data quality** (sample sizes, coordinate coverage, sequence lengths)
3. **Organize datasets** in `data/` directory:
   ```
   data/
   ├── Sphyrnidae/
   ├── Carcharhinidae/
   ├── Panulirus/
   ├── Salmonidae/
   └── Pieridae/
   ```
4. **Obtain custom shapefiles**:
   - HydroBASINS for Salmonidae
   - WWF Ecoregions for Pieridae
5. **Quick test runs** (small parameter sweep) to verify pipeline compatibility
6. **Begin benchmarking protocol** once all data ready

---

## Questions to Resolve

Before downloading, please confirm:

1. **Carcharhinidae vs. Rajidae** for Case Study 2? (Carcharhinidae recommended for size)
2. **Panulirus vs. Mytilus** for Case Study 3? (Panulirus recommended for phylogeography)
3. **Pieridae vs. Formicidae** for Case Study 5? (Pieridae recommended for coverage)
4. **Additional supplementary example** for marine/terrestrial spanning? Or keep to 5 core?

---

**Document prepared for dataset selection**
**Next step**: Download datasets from BOLD, organize, and proceed to benchmarking protocol
**Date**: 2026-01-28
