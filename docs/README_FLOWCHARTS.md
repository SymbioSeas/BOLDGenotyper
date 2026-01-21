# BOLDGenotyper Manuscript Flow Charts

This directory contains comprehensive documentation and flow charts for the BOLDGenotyper workflow, designed for manuscript preparation.

---

## 📁 Files in This Directory

### 📊 **Flow Chart Files (Mermaid Format)**

| File | Purpose | Use For |
|------|---------|---------|
| `flowchart_main_text_simplified.mmd` | High-level workflow overview (6 steps) | **Main text Figure 1 or 2** |
| `flowchart_supplemental_detailed.mmd` | Complete workflow (11 phases, all sub-steps) | **Supplemental Figure S1** |

### 📖 **Documentation Files**

| File | Purpose |
|------|---------|
| `BOLDGenotyper_Workflow_Documentation.md` | **50+ page comprehensive workflow documentation**<br/>• All 11 phases with inputs/outputs<br/>• Every processing step explained<br/>• All parameters (hardcoded vs adjustable)<br/>• Decision points clearly marked<br/>• Output interpretation guides<br/>• Publication-ready citation text |
| `FLOWCHART_RENDERING_GUIDE.md` | **Complete guide to converting Mermaid → SVG/PNG**<br/>• 5 rendering methods (online, CLI, VS Code, Python, GitHub)<br/>• Publication settings (300 DPI, color profiles)<br/>• Troubleshooting tips<br/>• Post-processing with Inkscape/Illustrator |
| `FLOWCHART_COMPARISON.md` | **Side-by-side comparison of simplified vs detailed**<br/>• When to use each version<br/>• What's included in each<br/>• Sample figure captions<br/>• Journal-specific recommendations |
| `README_FLOWCHARTS.md` | **This file** - Quick start guide |

### 🛠️ **Utility Scripts**

| File | Purpose |
|------|---------|
| `render_flowcharts.sh` | **Automated rendering script**<br/>Generates all SVG and PNG files with one command |

---

## 🚀 Quick Start

### Option 1: Render Online (Easiest, No Installation)

1. **Open Mermaid Live Editor:**
   ```
   https://mermaid.live
   ```

2. **For simplified version (main text):**
   - Copy contents of `flowchart_main_text_simplified.mmd`
   - Paste into editor
   - Click **"Actions"** → **"Download SVG"**

3. **For detailed version (supplemental):**
   - Copy contents of `flowchart_supplemental_detailed.mmd`
   - Paste into editor
   - Click **"Actions"** → **"Download SVG"**

4. **Done!** You now have publication-ready SVG files.

---

### Option 2: Render Locally (Best for High-Resolution)

**Prerequisites:**
```bash
# Install Node.js (if not already installed)
# macOS: brew install node
# Linux: sudo apt install nodejs npm
# Windows: Download from https://nodejs.org/

# Install Mermaid CLI
npm install -g @mermaid-js/mermaid-cli
```

**Render all diagrams:**
```bash
cd /Users/stesmith/Documents/depredation/boldgenotyper/docs/

# Run automated script
bash render_flowcharts.sh
```

**Output files:**
- `flowchart_main_text_simplified.svg` (scalable vector)
- `flowchart_main_text_simplified.png` (3000×4000 px, 300 DPI equivalent)
- `flowchart_supplemental_detailed.svg` (scalable vector)
- `flowchart_supplemental_detailed.png` (3500×8000 px, 300 DPI equivalent)

---

## 📋 What Each Flow Chart Includes

### **Simplified Version** (Main Text)
```
✓ 6 major processing phases
✓ 3 key decision points
✓ Sample count flow (N → M → K)
✓ Color-coded by step type
✓ ~10" wide × ~13" tall

Phases:
1. Quality Control & Filtering
2. Haplotype Discovery (ESV)
3. Assignment & Taxonomy
4. Geographic Analysis
5. Phylogenetic (optional)
6. Visualization & Reports
```

### **Detailed Version** (Supplemental)
```
✓ 11 complete phases
✓ 60+ processing sub-steps
✓ 10 decision points
✓ All parameter values shown
✓ ~12" wide × ~30" tall

All Phases:
1. Data Loading & Quality Marking
2. Pre-processing QC
3. Haplotype Discovery
4. Haplotype Assignment
5. Taxonomy Assignment
6. Post-assignment QC
6.5. Species-Level Aggregation
7. Geographic Enhancement
8. Phylogenetic Analysis
9. Divergence Analysis
10. Visualization
11. Report Generation
```

---

## 🎯 Which Version Should I Use?

| Scenario | Recommendation |
|----------|---------------|
| **Journal main text figure** | Simplified version |
| **Supplemental materials** | Detailed version |
| **Presentation/poster** | Simplified version |
| **Methods-focused paper** | Both versions |
| **GitHub README** | Simplified version (rendered automatically) |
| **Documentation** | Detailed version |

**Most common:** Use **simplified in main text** + **detailed in supplement**

See `FLOWCHART_COMPARISON.md` for detailed recommendations by journal type.

---

## 📐 Publication-Ready Settings

### **For Main Text (Simplified):**
- **Format:** SVG (preferred) or PNG
- **Resolution (PNG):** 3000×4000 px (~10"×13" at 300 DPI)
- **Background:** Transparent or white
- **File size:** ~50-100 KB (SVG), ~500 KB - 1 MB (PNG)

### **For Supplement (Detailed):**
- **Format:** SVG (highly recommended due to complexity)
- **Resolution (PNG):** 3500×8000 px (~12"×27" at 300 DPI)
- **Background:** Transparent or white
- **File size:** ~150-300 KB (SVG), ~2-4 MB (PNG)

---

## 🎨 Color Scheme

Both diagrams use a **publication-friendly, colorblind-safe** palette:

| Color | Use | Hex Code |
|-------|-----|----------|
| Light Blue | Input data | `#bbdefb` |
| Light Green | Processing steps | `#c8e6c9` |
| Light Yellow | Decision points | `#fff9c4` |
| Light Coral | Excluded data (dashed) | `#ffccbc` |
| Light Purple | Output data | `#e1bee7` |
| Light Cyan | Count/metric boxes | `#b3e5fc` |

**Tested with:** Color Oracle (colorblindness simulator) ✅

---

## 📝 Sample Figure Captions

### For Simplified Version (Main Text):

> **Figure 1. BOLDGenotyper workflow overview.** COI sequences undergo quality control (orientation normalization, ORF validation, dynamic filtering), haplotype discovery via exact sequence variant (ESV) approach with singleton error filtering, assignment and taxonomy determination by majority vote, geographic analysis with ocean basin assignment for high-quality coordinates, optional phylogenetic tree construction, and visualization. Blue: input; green: processing; yellow: decisions; coral dashed: excluded; purple: output. N = total samples, M = QC-passed, K = haplotypes.

### For Detailed Version (Supplemental):

> **Figure S1. Complete BOLDGenotyper workflow.** Detailed view of all 11 pipeline phases showing sub-steps, parameters, and decision points. Phases: (1) Data loading and coordinate quality marking, (2) Pre-processing QC, (3) Haplotype discovery with ESV approach, (4) Direct haplotype assignment, (5) Taxonomy assignment, (6) Post-assignment QC, (6.5) Species aggregation, (7) Geographic enhancement, (8) Phylogenetic analysis, (9) Divergence analysis, (10) Visualization, (11) Report generation. Yellow diamonds indicate user-configurable decisions. See Methods for parameter descriptions.

---

## 🔧 Customization

### Edit the Mermaid Code:

The `.mmd` files are plain text and can be edited in any text editor.

**Common edits:**
- Change box labels: Edit text between `[` and `]`
- Add/remove steps: Add/remove boxes and arrows
- Adjust colors: Modify `classDef` statements at bottom
- Change layout: Adjust `flowchart TD` (top-down) to `flowchart LR` (left-right)

### Post-Process in Vector Editor:

After rendering to SVG:
1. Open in **Inkscape** (free) or **Adobe Illustrator**
2. Fine-tune box positions, arrow curves
3. Adjust text size/font
4. Add annotations, callouts, legends
5. Export to PDF, EPS, or high-res PNG

---

## 📚 Additional Resources

| Resource | Location |
|----------|----------|
| **Complete workflow documentation** | `BOLDGenotyper_Workflow_Documentation.md` |
| **Rendering guide (5 methods)** | `FLOWCHART_RENDERING_GUIDE.md` |
| **Version comparison** | `FLOWCHART_COMPARISON.md` |
| **Mermaid documentation** | https://mermaid.js.org/ |
| **Mermaid Live Editor** | https://mermaid.live |
| **Color Oracle (colorblind simulator)** | https://colororacle.org/ |

---

## ✅ Manuscript Submission Checklist

Before submitting flow charts to journal:

- [ ] Rendered to correct format (SVG or PNG ≥300 DPI)
- [ ] All text legible at final print size
- [ ] Colors distinguishable in grayscale (if journal prints B&W)
- [ ] File size within journal limits
- [ ] Figure caption written
- [ ] Figure cited in manuscript text
- [ ] Meets journal-specific figure requirements
- [ ] Supplemental figure numbering correct

---

## 🆘 Troubleshooting

### "Command not found: mmdc"
```bash
npm install -g @mermaid-js/mermaid-cli
```

### "Diagram too large to render"
```bash
export NODE_OPTIONS="--max-old-space-size=4096"
mmdc -i input.mmd -o output.svg
```

### "Text is cut off"
Increase canvas size:
```bash
mmdc -i input.mmd -o output.png -w 4000 -H 6000
```

### "Need higher resolution"
```bash
mmdc -i input.mmd -o output.png -w 6000 -H 8000
```

For more help, see `FLOWCHART_RENDERING_GUIDE.md` or open an issue.

---

## 🎓 Citation

If you use these flow charts in your manuscript, please cite:

**Software:**
> Smith, S. (2025). BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis. Version 1.0. https://github.com/SymbioSeas/BOLDGenotyper

**Flow Chart Tool:**
> Mermaid. https://mermaid.js.org/

---

## 📧 Questions or Issues?

- **Documentation questions:** See `BOLDGenotyper_Workflow_Documentation.md`
- **Rendering issues:** See `FLOWCHART_RENDERING_GUIDE.md`
- **Bug reports:** https://github.com/SymbioSeas/BOLDGenotyper/issues
- **Email:** steph.smith@unc.edu

---

**Ready to render?** Run `bash render_flowcharts.sh` or visit https://mermaid.live 🚀
