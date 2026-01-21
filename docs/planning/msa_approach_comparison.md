# MSA Visualization Approach Comparison

## Option Comparison Table

| Feature | pymsaviz (Python) | ggmsa (R via rpy2) | Custom matplotlib |
|---------|-------------------|---------------------|-------------------|
| **Installation** | `pip install pymsaviz` | Requires R + rpy2 + BioConductor | No extra deps |
| **Learning curve** | Low (BioPython-like) | Medium (R knowledge needed) | High (build from scratch) |
| **Maintenance** | Single language | Two languages | High maintenance |
| **Features** | Consensus, logo, markers | Rich features, ggplot2 | Full control, more work |
| **Performance** | Fast | Fast (but subprocess overhead) | Fast |
| **Output quality** | Publication-ready | Publication-ready | Depends on effort |
| **Integration** | Native Python objects | Conversion overhead | Native |
| **Dependencies** | matplotlib, BioPython | R, rpy2, Biostrings, ggmsa | Just matplotlib |

## Recommended: pymsaviz

### Pros:
✅ Designed specifically for MSA visualization
✅ Very similar API to ggmsa
✅ Works with BioPython (already used in BOLDGenotyper)
✅ Active development and good documentation
✅ No language barrier for contributors
✅ Easier CI/CD (no R installation needed)

### Cons:
⚠ Slightly fewer styling options than ggmsa
⚠ Newer package (but mature enough for production)

## Code Example Comparison

### Your R Code:
```r
dss <- Biostrings::readDNAStringSet(aln_path)
dss <- dss[tip_order]

p_msa <- ggmsa(dss,
               char_width = 0.7,
               start = 1,
               end = 92,
               color = "Chemistry_NT",
               consensus_views = TRUE,
               disagreement = TRUE) +
  geom_seqlogo(color = "Chemistry_NT") +
  geom_msaBar()
```

### Equivalent Python (pymsaviz):
```python
from pymsaviz import MsaViz
from Bio import AlignIO, Phylo

# Read alignment
alignment = AlignIO.read(aln_path, "fasta")

# Read tree and get tip order
tree = Phylo.read(tree_path, "newick")
tip_order = [tip.name for tip in tree.get_terminals()]

# Reorder alignment by tree
ordered_records = sorted(
    alignment,
    key=lambda x: tip_order.index(x.id) if x.id in tip_order else 999
)

# Create MSA visualization
mv = MsaViz(
    ordered_records,
    start=1,
    end=92,
    color_scheme="Chemistry_NT",
    wrap_length=None,
    show_label=True,
    show_seq_char=True
)

# Add consensus and conservation
mv.add_consensus_bar(color="Chemistry_NT")
mv.add_seq_logo(color_scheme="Chemistry_NT")

# Save
mv.savefig("output.pdf", dpi=300)
```

## Alternative: BioViz (Another Python Option)

```python
# Alternative using a different package
from bioviz import MSAVisualizer

viz = MSAVisualizer(alignment)
viz.set_color_scheme("nucleotide")
viz.add_consensus()
viz.add_logo()
viz.save("output.pdf")
```

## Why Not R Integration?

### rpy2 Approach Issues:
❌ Requires R installation on all systems
❌ BioConductor setup complexity
❌ Version compatibility issues (R, rpy2, packages)
❌ Harder to debug across language boundary
❌ More complex CI/CD pipeline
❌ Contributors need R knowledge

### When R Might Be Worth It:
- If you already have R infrastructure
- If you need very specific ggmsa features
- If your users already use R
- If you're OK with the added complexity

**For BOLDGenotyper: Not worth it**
- Python-first ecosystem
- Installation simplicity is important
- Contributors are Python-focused
- pymsaviz provides equivalent functionality

## Decision: Go with pymsaviz

**Final recommendation:**
1. Use `pymsaviz` for MSA visualization
2. Implement as optional feature (graceful degradation)
3. Keep R code in your analysis scripts for custom work
4. Document both approaches for users who want R output
