# MSA Visualization Implementation Summary

## Executive Summary

**Recommendation: Use Python (`pymsaviz`) instead of R integration**

**Rationale:**
- ✅ Maintains Python-only ecosystem
- ✅ Simpler installation and maintenance
- ✅ Better integration with existing BioPython code
- ✅ Equivalent functionality to ggmsa
- ✅ Publication-quality output

## Implementation Phases

### Phase 1: Core Implementation (Week 1)
**Goal:** Basic MSA visualization working

**Tasks:**
1. ☐ Add `pymsaviz>=0.4.0` to requirements.txt
2. ☐ Create `boldgenotyper/msa_visualization.py` module
3. ☐ Implement core functions:
   - `create_phylo_ordered_msa()`
   - `chunk_alignment()`
   - `order_by_tree()`
4. ☐ Add `MSAConfig` to config.py
5. ☐ Write comprehensive docstrings

**Deliverables:**
- Working MSA plot generation
- Phylogeny-ordered sequences
- Chunked display for long alignments

### Phase 2: Integration (Week 2)
**Goal:** Integrated into phylogenetic workflow

**Tasks:**
1. ☐ Modify `phylogenetics.py` to call MSA visualization
2. ☐ Add error handling and graceful degradation
3. ☐ Implement adaptive sizing logic
4. ☐ Add logging and progress indicators
5. ☐ Test on Sphyrnidae dataset

**Deliverables:**
- MSA plots generated automatically after tree building
- Robust error handling
- Multiple output formats (PDF, PNG)

### Phase 3: Testing & Documentation (Week 3)
**Goal:** Production-ready with complete documentation

**Tasks:**
1. ☐ Write unit tests (90%+ coverage)
2. ☐ Write integration tests
3. ☐ Update README with MSA documentation
4. ☐ Create usage examples
5. ☐ Document configuration options
6. ☐ Add troubleshooting guide

**Deliverables:**
- Comprehensive test suite
- Complete user documentation
- Developer documentation

### Phase 4: Enhancement (Optional, Later)
**Goal:** Advanced features

**Tasks:**
1. ☐ Interactive HTML output
2. ☐ SNP/variant annotation
3. ☐ Codon position highlighting
4. ☐ Reference sequence comparison
5. ☐ Custom color schemes

## File Changes Required

### New Files:
```
boldgenotyper/msa_visualization.py  (NEW - ~300 lines)
tests/test_msa_visualization.py     (NEW - ~200 lines)
docs/msa_visualization.md           (NEW - documentation)
```

### Modified Files:
```
boldgenotyper/config.py             (ADD MSAConfig class)
boldgenotyper/phylogenetics.py      (ADD MSA call)
requirements.txt                     (ADD pymsaviz)
README.md                           (ADD MSA section)
```

## Configuration Example

```python
# config.py
from dataclasses import dataclass, field
from typing import List

@dataclass(frozen=True)
class MSAConfig:
    """
    Configuration for MSA visualization.
    
    Attributes
    ----------
    enabled : bool
        Whether to generate MSA plots (default: True)
    chunk_size : int
        Number of alignment positions per chunk (default: 100)
    max_sequences : int
        Maximum sequences to display (default: 50)
    color_scheme : str
        Nucleotide color scheme (default: "Chemistry_NT")
    show_consensus : bool
        Show consensus sequence bar (default: True)
    show_logo : bool
        Show sequence logo (default: True)
    show_conservation : bool
        Show conservation bar chart (default: True)
    output_formats : List[str]
        Output file formats (default: ["pdf", "png"])
    """
    enabled: bool = True
    chunk_size: int = 100
    max_sequences: int = 50
    color_scheme: str = "Chemistry_NT"
    show_consensus: bool = True
    show_logo: bool = True
    show_conservation: bool = True
    output_formats: List[str] = field(default_factory=lambda: ["pdf", "png"])
```

## Core Function Skeleton

```python
# msa_visualization.py
from pathlib import Path
from typing import List, Optional, Tuple
import logging

from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from pymsaviz import MsaViz

logger = logging.getLogger(__name__)


def create_phylo_ordered_msa(
    alignment_path: Path,
    tree_path: Path,
    output_dir: Path,
    organism: str,
    chunk_size: int = 100,
    max_sequences: int = 50,
    color_scheme: str = "Chemistry_NT",
    show_consensus: bool = True,
    show_logo: bool = True,
    output_formats: List[str] = None
) -> List[Path]:
    """
    Create phylogeny-ordered MSA visualization with chunking.
    
    Generates multiple MSA plots for long alignments, with sequences
    ordered by phylogenetic tree topology.
    
    Parameters
    ----------
    alignment_path : Path
        Path to aligned FASTA file
    tree_path : Path
        Path to phylogenetic tree (Newick)
    output_dir : Path
        Output directory for plots
    organism : str
        Organism name for file naming
    chunk_size : int, optional
        Alignment positions per chunk (default: 100)
    max_sequences : int, optional
        Maximum sequences to display (default: 50)
    color_scheme : str, optional
        Color scheme for nucleotides (default: "Chemistry_NT")
    show_consensus : bool, optional
        Show consensus bar (default: True)
    show_logo : bool, optional
        Show sequence logo (default: True)
    output_formats : List[str], optional
        Output formats (default: ["pdf", "png"])
        
    Returns
    -------
    List[Path]
        Paths to generated MSA plot files
        
    Notes
    -----
    For alignments >100 sequences, the function will downsample to
    max_sequences phylogenetically diverse representatives.
    
    For alignments >chunk_size bp, multiple plots are generated:
    - {organism}_MSA_chunk1.pdf (positions 1-100)
    - {organism}_MSA_chunk2.pdf (positions 101-200)
    - etc.
    
    Examples
    --------
    >>> plots = create_phylo_ordered_msa(
    ...     alignment_path="aligned.fasta",
    ...     tree_path="tree.nwk",
    ...     output_dir="output/",
    ...     organism="Sphyrnidae"
    ... )
    >>> len(plots)
    7  # For 650 bp alignment with 100 bp chunks
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if output_formats is None:
        output_formats = ["pdf", "png"]
    
    logger.info(f"Generating MSA visualization for {organism}")
    
    # Load alignment
    alignment = AlignIO.read(alignment_path, "fasta")
    logger.info(f"  Loaded alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} bp")
    
    # Order by tree
    ordered_alignment = order_by_tree(alignment, tree_path)
    
    # Downsample if needed
    if len(ordered_alignment) > max_sequences:
        ordered_alignment = downsample_alignment(
            ordered_alignment, 
            max_sequences,
            tree_path
        )
        logger.info(f"  Downsampled to {len(ordered_alignment)} sequences")
    
    # Create chunks
    chunks = chunk_alignment(ordered_alignment, chunk_size)
    logger.info(f"  Creating {len(chunks)} chunks")
    
    # Generate plots for each chunk
    output_files = []
    for i, (start, end, chunk_records) in enumerate(chunks, 1):
        for fmt in output_formats:
            output_file = output_dir / f"{organism}_MSA_chunk{i}.{fmt}"
            
            # Create MSA plot
            mv = MsaViz(
                chunk_records,
                color_scheme=color_scheme,
                show_label=True,
                show_seq_char=True,
                wrap_length=None
            )
            
            # Add features
            if show_consensus:
                mv.add_consensus_bar(color=color_scheme)
            if show_logo:
                mv.add_seq_logo(color_scheme=color_scheme)
            
            # Save
            mv.savefig(output_file, dpi=300 if fmt == "png" else None)
            output_files.append(output_file)
            logger.info(f"  Saved chunk {i} ({start}-{end} bp): {output_file}")
    
    return output_files


def order_by_tree(
    alignment: MultipleSeqAlignment,
    tree_path: Path
) -> MultipleSeqAlignment:
    """
    Reorder alignment sequences to match phylogenetic tree tip order.
    
    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Input alignment
    tree_path : Path
        Path to phylogenetic tree
        
    Returns
    -------
    MultipleSeqAlignment
        Reordered alignment
    """
    # Read tree
    tree = Phylo.read(tree_path, "newick")
    
    # Get tip order (leaf nodes in tree traversal order)
    tip_order = [tip.name.replace("'", "") for tip in tree.get_terminals()]
    
    # Create mapping
    seq_dict = {rec.id: rec for rec in alignment}
    
    # Reorder
    ordered_records = []
    for tip in tip_order:
        if tip in seq_dict:
            ordered_records.append(seq_dict[tip])
    
    # Add any sequences not in tree (shouldn't happen but safe)
    remaining = [rec for rec in alignment if rec.id not in tip_order]
    ordered_records.extend(remaining)
    
    return MultipleSeqAlignment(ordered_records)


def chunk_alignment(
    alignment: MultipleSeqAlignment,
    chunk_size: int
) -> List[Tuple[int, int, List]]:
    """
    Split alignment into chunks for visualization.
    
    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Full alignment
    chunk_size : int
        Number of positions per chunk
        
    Returns
    -------
    List[Tuple[int, int, List]]
        List of (start, end, records) tuples
    """
    aln_length = alignment.get_alignment_length()
    chunks = []
    
    for start in range(0, aln_length, chunk_size):
        end = min(start + chunk_size, aln_length)
        
        # Extract chunk from each sequence
        chunk_records = [rec[start:end] for rec in alignment]
        
        chunks.append((start + 1, end, chunk_records))  # 1-indexed
    
    return chunks


def downsample_alignment(
    alignment: MultipleSeqAlignment,
    max_sequences: int,
    tree_path: Path
) -> MultipleSeqAlignment:
    """
    Downsample alignment to max_sequences phylogenetically diverse representatives.
    
    Uses tree to select phylogenetically dispersed sequences rather than
    random sampling.
    
    Parameters
    ----------
    alignment : MultipleSeqAlignment
        Full alignment
    max_sequences : int
        Target number of sequences
    tree_path : Path
        Phylogenetic tree for guided sampling
        
    Returns
    -------
    MultipleSeqAlignment
        Downsampled alignment
    """
    # TODO: Implement phylogenetically-informed downsampling
    # For now, just take first max_sequences
    logger.warning(
        f"Downsampling from {len(alignment)} to {max_sequences} sequences. "
        "Using simple truncation for now."
    )
    return alignment[:max_sequences]
```

## Integration Code

```python
# Add to phylogenetics.py

def run_phylogenetic_analysis(
    # ... existing parameters ...
    cfg: Config = None
):
    # ... existing tree building code ...
    
    # Add MSA visualization
    if cfg.msa.enabled:
        logger.info("Generating MSA visualization...")
        
        try:
            from boldgenotyper import msa_visualization
            
            msa_plots = msa_visualization.create_phylo_ordered_msa(
                alignment_path=aligned_fasta,
                tree_path=output_tree,
                output_dir=output_dir,
                organism=organism,
                chunk_size=cfg.msa.chunk_size,
                max_sequences=cfg.msa.max_sequences,
                color_scheme=cfg.msa.color_scheme,
                show_consensus=cfg.msa.show_consensus,
                show_logo=cfg.msa.show_logo,
                output_formats=cfg.msa.output_formats
            )
            
            logger.info(f"  Generated {len(msa_plots)} MSA plots")
            
        except ImportError as e:
            logger.warning(
                f"MSA visualization skipped: pymsaviz not installed. "
                f"Install with: pip install pymsaviz"
            )
        except Exception as e:
            logger.warning(f"MSA visualization failed: {e}")
            logger.debug("Full traceback:", exc_info=True)
```

## Testing Checklist

- [ ] Install pymsaviz: `pip install pymsaviz`
- [ ] Test with small alignment (5 seqs, 100 bp)
- [ ] Test with medium alignment (20 seqs, 650 bp)
- [ ] Test with large alignment (100 seqs, 1000 bp)
- [ ] Test phylo-ordering matches tree
- [ ] Test chunking produces correct ranges
- [ ] Test with missing pymsaviz (graceful degradation)
- [ ] Test with missing tree file
- [ ] Verify output quality (colors, consensus, logo)
- [ ] Check memory usage on large datasets
- [ ] Validate all output formats (PDF, PNG)

## Expected Outputs

For a typical COI dataset (50 haplotypes, 650 bp):

```
phylogenetic/
├── Sphyrnidae_MSA_chunk1.pdf    (~200 KB, positions 1-100)
├── Sphyrnidae_MSA_chunk1.png    (~500 KB, positions 1-100)
├── Sphyrnidae_MSA_chunk2.pdf    (~200 KB, positions 101-200)
├── Sphyrnidae_MSA_chunk2.png    (~500 KB, positions 101-200)
├── Sphyrnidae_MSA_chunk3.pdf
├── Sphyrnidae_MSA_chunk3.png
├── Sphyrnidae_MSA_chunk4.pdf
├── Sphyrnidae_MSA_chunk4.png
├── Sphyrnidae_MSA_chunk5.pdf
├── Sphyrnidae_MSA_chunk5.png
├── Sphyrnidae_MSA_chunk6.pdf
├── Sphyrnidae_MSA_chunk6.png
├── Sphyrnidae_MSA_chunk7.pdf    (~200 KB, positions 601-650)
└── Sphyrnidae_MSA_chunk7.png    (~500 KB, positions 601-650)

Total: ~9.8 MB for complete MSA visualization
```

## Next Steps

1. **Review this plan** - Ensure it meets your needs
2. **Install pymsaviz** - `pip install pymsaviz` and test basic functionality
3. **Create module skeleton** - Start with `msa_visualization.py`
4. **Implement core function** - Get basic MSA plot working
5. **Add phylo-ordering** - Match your R code's tree-based ordering
6. **Add chunking** - Split long alignments
7. **Integrate** - Add to phylogenetic workflow
8. **Test** - Run on Sphyrnidae dataset
9. **Document** - Update README and docstrings
10. **Iterate** - Refine based on results

## Questions to Consider

1. **Default chunk size** - Is 100 bp appropriate for COI? (current: 92 bp in R code)
2. **Max sequences** - What's reasonable? 50? 100?
3. **Output formats** - PDF + PNG? Add SVG?
4. **Consensus features** - Which are most useful?
5. **Color scheme** - Chemistry_NT only, or offer options?

Let me know if you'd like me to proceed with implementation!
