# CLI Migration Implementation Plan

## Phase Reorganization

### OLD STRUCTURE:
```
Phase 1: Data Loading and Quality Control
  - Parse TSV
  - Mark coordinate quality
  - Assign ocean basins

Phase 2: Sequence Dereplication and Consensus Generation
  - Generate FASTA
  - dereplicate_from_fasta() [clustering-based]

Phase 3: Genotype Assignment
  - assign_genotypes()

Phase 4: Taxonomy Assignment
  - assign_consensus_taxonomy()

Phase 4.5: Quality Control and Contamination Detection
  - add_contamination_columns()
  - generate_quality_control_report()

Phase 4.7: Geographic Analysis Enhancement

Phase 5: Phylogenetic Analysis

Phase 5.5: Divergence Analysis

Phase 6: Visualization
```

### NEW STRUCTURE:
```
Phase 1: Data Loading & Coordinate Quality
  - Parse TSV
  - Mark coordinate quality
  - Assign ocean basins
  [KEEP AS-IS]

Phase 2: Pre-processing QC (NEW)
  - Generate initial FASTA from raw sequences
  - Orientation normalization (apply_orientation_normalization)
  - Dynamic QC filtering (apply_dynamic_qc_filters)
  - Save QC-passed FASTA for haplotype discovery

Phase 3: Haplotype Discovery (REPLACES OLD PHASE 2)
  - identify_haplotypes() [NEW haplotype-first workflow]
    * Aligns sequences
    * Extracts core region
    * Identifies unique haplotypes (ESVs)
    * Flags suspect haplotypes
  - Outputs: haplotypes.fasta, haplotype_mapping.csv, haplotype_stats.csv

Phase 4: Haplotype Assignment (REPLACES OLD PHASE 3)
  - assign_haplotypes() [renamed from assign_genotypes]
  - Use haplotypes.fasta instead of consensus.fasta

Phase 5: Taxonomy Assignment (OLD PHASE 4)
  - assign_consensus_taxonomy() → update to use haplotype_id column
  - Update all references: consensus_group → haplotype_id

Phase 6: Post-assignment QC (OLD PHASE 4.5)
  - Contamination detection
  - Update to use haplotype_id

Phase 7: Geographic Enhancement (OLD PHASE 4.7)
  [KEEP AS-IS]

Phase 8: Phylogenetic Analysis (OLD PHASE 5)
  - Update to use haplotypes.fasta
  - Update terminology in outputs

Phase 9: Divergence Analysis (OLD PHASE 5.5)
  - Update to use haplotypes.fasta

Phase 10: Visualization (OLD PHASE 6)
  - Update to use haplotype terminology
```

## Key File Path Changes

### Inputs Created:
- `intermediate/qc/{organism}_oriented.fasta` - After orientation normalization
- `intermediate/qc/{organism}_qc_passed.fasta` - After dynamic QC
- `intermediate/qc/orf_validation_results.csv` - ORF stats per sample
- `intermediate/qc/qc_filter_results.csv` - QC stats

### Haplotype Discovery Outputs:
- `haplotypes/{organism}_haplotypes.fasta` - Main haplotype sequences
- `haplotypes/{organism}_haplotype_mapping.csv` - Sample → haplotype mapping
- `haplotypes/{organism}_haplotype_stats.csv` - Haplotype stats with flags
- `intermediate/haplotype_discovery/{organism}_aligned.fasta` - Full alignment
- `intermediate/haplotype_discovery/{organism}_core_region.fasta` - Core alignment (diagnostic)

### Assignment Outputs:
- `haplotype_assignments/{organism}_with_haplotypes.tsv` - Annotated metadata
- `haplotype_assignments/{organism}_diagnostics.csv` - Assignment diagnostics

## Column Name Changes

Throughout the pipeline:
- `consensus_group` → `haplotype_id`
- `consensus_group_sp` → `haplotype_sp`
- `consensus_group_short` → `haplotype_short`
- `n_assigned_to_consensus` → `n_assigned_to_haplotype`

## Critical Integration Points

### Downstream Module Updates Needed:
1. **taxonomy assignment** - Use `haplotype_id` instead of `consensus_group`
2. **quality_control** - Update group_col parameter to `haplotype_id`
3. **divergence_analysis** - Read from `haplotypes/{organism}_haplotypes.fasta`
4. **phylogenetics** - Read from `haplotypes/{organism}_haplotypes.fasta`
5. **visualization** - Update all plots to use `haplotype_id` column
6. **reports** - Update all references to use haplotype terminology

## Error Handling Strategy

Permissive approach:
- Log warnings for failed sequences
- Continue with valid sequences
- Track and report statistics:
  * N sequences loaded
  * N failed orientation normalization
  * N failed ORF validation
  * N failed dynamic QC
  * N passed to haplotype discovery
  * N haplotypes identified
  * N haplotypes flagged as suspect

## CLI Arguments Update

DEPRECATE:
- `--clustering-threshold`

KEEP/ADD:
- `--min-core-coverage` (default: 0.8)
- `--max-singleton-distance` (default: 0.05)
- `--build-tree` (existing)
- `--no-geo` (existing)
- `--threads` (existing)

Use config.py defaults for all other parameters.
