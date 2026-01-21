# MSA Visualization Integration Plan for BOLDGenotyper

## 1. Package Selection: pymsaviz

**Why pymsaviz:**
- Purpose-built for MSA visualization
- Similar functionality to ggmsa
- Supports consensus sequence, sequence logo, and position markers
- Works with BioPython alignments (already used in BOLDGenotyper)
- Good documentation and active maintenance
- Publication-quality output

**Alternative considered:**
- Custom matplotlib solution - too much development overhead
- plotly - good but more complex for static images
- rpy2 + ggmsa - adds R dependency complexity

## 2. Implementation Location

**New module:** `boldgenotyper/msa_visualization.py`

**Integration points:**
1. Phylogenetic workflow (after tree generation)
2. Haplotype discovery (core region visualization)
3. Quality control (contamination sequence comparison)

## 3. Key Features to Implement

### 3.1 Phylogeny-Ordered MSA
- Read aligned FASTA from phylogenetic workflow
- Order sequences by phylogenetic tree (like R chunk)
- Display in tree tip order for biological interpretation

### 3.2 Chunked Visualization
- Handle large alignments (500+ bp) by creating multiple panels
- Default chunk size: 100 bp (configurable)
- Auto-generate multiple output files (pt1, pt2, pt3, etc.)

### 3.3 Customization Options
- Color scheme (Chemistry_NT for nucleotides)
- Consensus sequence display
- Sequence logo
- Conservation bar chart
- Configurable width and height

### 3.4 Adaptive Sizing
- Small datasets (<20 sequences): show all
- Medium datasets (20-50): compress spacing
- Large datasets (>50): sample or paginate

## 4. File Structure

```
boldgenotyper/
├── msa_visualization.py (NEW)
│   ├── create_msa_plot()
│   ├── create_phylo_ordered_msa()
│   ├── chunk_alignment()
│   └── get_adaptive_sizing()
├── phylogenetics.py (MODIFY)
│   └── Call MSA visualization after tree building
├── config.py (MODIFY)
│   └── Add MSAConfig class
└── visualization.py (EXISTING - could merge here)
```

## 5. Configuration

Add to `config.py`:

```python
@dataclass(frozen=True)
class MSAConfig:
    """Configuration for MSA visualization."""
    
    enabled: bool = True
    chunk_size: int = 100
    color_scheme: str = "Chemistry_NT"
    show_consensus: bool = True
    show_logo: bool = True
    show_conservation: bool = True
    max_sequences: int = 50
    char_width: float = 0.7
    output_formats: List[str] = field(default_factory=lambda: ["pdf", "png"])
```

## 6. Output Files

Location: `<output_dir>/phylogenetic/`

Files generated:
- `{organism}_MSA_chunk1.pdf` - First 100 bp
- `{organism}_MSA_chunk2.pdf` - Next 100 bp
- `{organism}_MSA_chunk1.png` - High-res PNG
- `{organism}_MSA_full_thumbnail.pdf` - Overview of full alignment

## 7. Dependencies

Add to requirements.txt:
```
pymsaviz>=0.4.0
```

## 8. Code Organization

### 8.1 Main Function

```python
def create_phylo_ordered_msa(
    alignment_path: Path,
    tree_path: Path,
    output_dir: Path,
    organism: str,
    chunk_size: int = 100,
    config: MSAConfig = None
) -> List[Path]:
    """
    Create phylogeny-ordered MSA visualization.
    
    Generates chunked MSA plots with sequences ordered by phylogenetic
    tree topology. Includes consensus sequence, sequence logo, and
    conservation metrics.
    
    Parameters
    ----------
    alignment_path : Path
        Path to aligned FASTA file
    tree_path : Path
        Path to phylogenetic tree (Newick format)
    output_dir : Path
        Output directory for MSA plots
    organism : str
        Organism name for file naming
    chunk_size : int, optional
        Number of alignment positions per chunk (default: 100)
    config : MSAConfig, optional
        Configuration object
        
    Returns
    -------
    List[Path]
        Paths to generated MSA plot files
    """
```

### 8.2 Helper Functions

```python
def chunk_alignment(alignment, chunk_size):
    """Split alignment into chunks for visualization."""
    
def order_by_tree(alignment, tree):
    """Reorder alignment sequences to match tree tip order."""
    
def get_adaptive_sizing(n_sequences, alignment_length):
    """Calculate optimal figure dimensions."""
    
def add_sequence_annotations(plot, alignment, positions):
    """Add SNP positions, codon markers, etc."""
```

## 9. Error Handling

**Robust handling for:**
- Missing tree file (skip phylo-ordering, use default order)
- Large number of sequences (>100) - sample or warn
- Very long alignments (>1000 bp) - auto-chunk more aggressively
- Memory constraints - downsample resolution
- Missing pymsaviz - graceful degradation with warning

## 10. Integration with Phylogenetic Workflow

Modify `phylogenetics.py`:

```python
def run_phylogenetic_analysis(...):
    # ... existing tree building code ...
    
    # Add MSA visualization
    if cfg.msa.enabled:
        logger.info("Generating MSA visualization...")
        try:
            msa_plots = msa_visualization.create_phylo_ordered_msa(
                alignment_path=aligned_fasta,
                tree_path=tree_file,
                output_dir=output_dir,
                organism=organism,
                config=cfg.msa
            )
            logger.info(f"  Generated {len(msa_plots)} MSA plots")
        except ImportError:
            logger.warning("pymsaviz not installed - skipping MSA plots")
        except Exception as e:
            logger.warning(f"MSA visualization failed: {e}")
```

## 11. Documentation

### 11.1 Module Docstring

```python
"""
MSA Visualization Module

Generates publication-quality multiple sequence alignment visualizations
with phylogenetic ordering, consensus sequences, and conservation metrics.

Features:
- Phylogeny-ordered sequences (matches tree topology)
- Chunked display for long alignments
- Nucleotide coloring (Chemistry_NT scheme)
- Consensus sequence and sequence logo
- Conservation bar charts
- Adaptive sizing for different dataset sizes

The module integrates with the phylogenetic workflow to create MSA plots
ordered by tree topology, making it easy to identify conserved regions
and lineage-specific mutations.
"""
```

### 11.2 README Section

Add to main README.md:

```markdown
### MSA Visualization

BOLDGenotyper automatically generates multiple sequence alignment (MSA)
visualizations ordered by phylogenetic tree topology. These plots help
identify:

- Conserved regions across haplotypes
- Lineage-specific mutations
- Sequence quality issues
- Codon structure

**Output files:**
- `phylogenetic/{organism}_MSA_chunk*.pdf` - Chunked MSA views
- Color-coded nucleotides (Chemistry_NT scheme)
- Consensus sequence and sequence logo
- Conservation metrics

**Configuration:**
```yaml
msa:
  enabled: true
  chunk_size: 100  # bp per chunk
  show_consensus: true
  show_logo: true
```

## 12. Testing Strategy

### 12.1 Unit Tests

```python
def test_chunk_alignment():
    """Test alignment chunking logic."""
    
def test_order_by_tree():
    """Test phylogenetic ordering."""
    
def test_adaptive_sizing():
    """Test figure size calculation."""
```

### 12.2 Integration Tests

```python
def test_create_msa_plot_small():
    """Test with small alignment (5 seqs, 100 bp)."""
    
def test_create_msa_plot_large():
    """Test with large alignment (50 seqs, 650 bp)."""
    
def test_missing_dependencies():
    """Test graceful degradation without pymsaviz."""
```

### 12.3 Visual Tests

- Generate test plots with known alignments
- Manually verify: ordering, colors, consensus, logo
- Compare with ggmsa output for validation

## 13. Performance Considerations

**Memory:**
- Limit to 100 sequences max by default
- Sample if >100 (take phylogenetically diverse subset)
- Clear matplotlib figures after saving

**Speed:**
- pymsaviz is fast (<1 sec for 50 seqs × 650 bp)
- Chunking is efficient
- Only bottleneck is I/O (writing multiple files)

**Disk:**
- Each chunk ~100-500 KB (PDF)
- 650 bp alignment → ~7 chunks → ~3.5 MB total
- Reasonable for modern systems

## 14. Future Enhancements

**Phase 1 (initial):**
- Basic MSA plots ordered by tree
- Chunked display
- Consensus and logo

**Phase 2 (later):**
- Interactive HTML output (plotly)
- SNP annotation overlay
- Codon position highlighting
- Reference sequence comparison
- Export to common formats (MEGA, Jalview)

## 15. Migration Path

**Step 1:** Implement core module
**Step 2:** Add to phylogenetic workflow as optional
**Step 3:** Test on Sphyrnidae dataset
**Step 4:** Enable by default in config
**Step 5:** Document in user guide
**Step 6:** Collect feedback and iterate

## 16. Backward Compatibility

- Feature is opt-in via config (enabled by default)
- Missing pymsaviz → warning, not error
- No changes to existing outputs
- Pipeline runs normally if MSA generation fails
