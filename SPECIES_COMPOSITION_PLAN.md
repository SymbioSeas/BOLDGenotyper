# Species Composition Tracking: Implementation Plan

## Problem Statement

When taxonomy assignment reverts to genus level or cannot be confidently assigned to species, it's because reference sequences from BOLD were reported under multiple conflicting species names. Users need to easily see which reported species contributed to each consensus group to:
- Identify potential misidentifications
- Understand taxonomic complexity in their data
- Make informed decisions about genotype interpretation

**Example Use Case**: A consensus group "Rhizoprionodon c13_n73" might contain samples reported as both *R. porosus* and *R. terranovae*. Geographic distribution analysis suggests some *R. terranovae* reports may be misidentifications.

## Current State

### What Works
- Species information is preserved from BOLD TSV input through the pipeline
- Output file `{organism}_with_genotypes.tsv` contains both `consensus_group` and original `species` columns
- `reports.py` has existing code (lines 74-80) that computes species composition for conflict reports

### What's Missing
- No explicit species composition summary output file
- Species composition not prominently displayed in HTML report
- No way to quickly identify which consensus groups have multiple species
- No flagging of taxonomically ambiguous consensus groups

## Proposed Solution

### Phase 1: Data Collection & Analysis

**Location**: `boldgenotyper/genotype_assignment.py`

Add new function after genotype assignment:

```python
def compute_species_composition(
    assignments_with_metadata: pd.DataFrame,
    min_abundance_pct: float = 1.0
) -> pd.DataFrame:
    """
    Compute species composition for each consensus group.

    Args:
        assignments_with_metadata: DataFrame with consensus_group and species columns
        min_abundance_pct: Only report species above this percentage (default 1%)

    Returns:
        DataFrame with columns:
            - consensus_group
            - n_total_samples
            - n_species
            - species_list (all species, comma-separated)
            - primary_species (most abundant species)
            - primary_species_pct
            - species_composition (detailed breakdown with counts and percentages)
            - is_multi_species (True if >1 species with >min_abundance_pct)
            - is_ambiguous (True if no species has >70% abundance)
    """

    # Filter to only assigned samples
    assigned = assignments_with_metadata[
        assignments_with_metadata['consensus_group'].notna()
    ].copy()

    # Group by consensus_group and compute composition
    composition_data = []

    for group, group_df in assigned.groupby('consensus_group'):
        species_counts = group_df['species'].value_counts()
        total_samples = len(group_df)

        # Filter species by min_abundance
        species_pcts = (species_counts / total_samples * 100)
        significant_species = species_pcts[species_pcts >= min_abundance_pct]

        # Build detailed composition string
        composition_parts = [
            f"{sp}: {species_counts[sp]} ({species_pcts[sp]:.1f}%)"
            for sp in significant_species.index
        ]
        composition_str = "; ".join(composition_parts)

        # Compute flags
        is_multi_species = len(significant_species) > 1
        primary_species = species_counts.index[0]
        primary_pct = species_pcts.iloc[0]
        is_ambiguous = primary_pct < 70.0

        composition_data.append({
            'consensus_group': group,
            'n_total_samples': total_samples,
            'n_species': len(significant_species),
            'species_list': ", ".join(significant_species.index),
            'primary_species': primary_species,
            'primary_species_pct': f"{primary_pct:.1f}%",
            'species_composition': composition_str,
            'is_multi_species': is_multi_species,
            'is_ambiguous': is_ambiguous
        })

    return pd.DataFrame(composition_data).sort_values('n_total_samples', ascending=False)
```

**Integration Point**: Call this function in `assign_genotypes()` after line 1076, where the output TSV is written.

### Phase 2: CSV Output

**Location**: `boldgenotyper/genotype_assignment.py`

Add to `assign_genotypes()` function:

```python
# After writing {organism}_with_genotypes.tsv (line 1076)

# Compute species composition
species_comp = compute_species_composition(
    merged_output,
    min_abundance_pct=1.0  # Could be added to config
)

# Write species composition file
composition_path = output_dir / f"{organism}_species_composition.csv"
species_comp.to_csv(composition_path, index=False)
logger.info(f"Species composition written to {composition_path}")

# Flag multi-species groups
multi_species_groups = species_comp[species_comp['is_multi_species']]
if not multi_species_groups.empty:
    logger.warning(
        f"Found {len(multi_species_groups)} consensus groups with multiple species. "
        f"See {composition_path} for details."
    )
```

**New Output File**: `{organism}_species_composition.csv`

Example content:
```csv
consensus_group,n_total_samples,n_species,species_list,primary_species,primary_species_pct,species_composition,is_multi_species,is_ambiguous
Rhizoprionodon c13_n73,73,2,"Rhizoprionodon porosus, Rhizoprionodon terranovae",Rhizoprionodon porosus,67.1%,"Rhizoprionodon porosus: 49 (67.1%); Rhizoprionodon terranovae: 24 (32.9%)",True,True
Sphyrna c1_n842,842,1,Sphyrna lewini,Sphyrna lewini,100.0%,Sphyrna lewini: 842 (100.0%),False,False
Carcharhinus c7_n156,156,3,"Carcharhinus limbatus, Carcharhinus amblyrhynchoides, Carcharhinus tilstoni",Carcharhinus limbatus,82.1%,"Carcharhinus limbatus: 128 (82.1%); Carcharhinus amblyrhynchoides: 18 (11.5%); Carcharhinus tilstoni: 10 (6.4%)",True,False
```

### Phase 3: HTML Report Enhancement

**Location**: `boldgenotyper/reports.py`

#### 3.1 Add Species Composition Section

Add new section after "Genotype Assignments" (around line 2750):

```python
def _build_species_composition_section(
    results_dir: Path,
    organism: str
) -> str:
    """Build HTML section showing species composition of consensus groups."""

    html = '<div class="tab-content" id="species-composition">\n'
    html += '<h2>Species Composition of Consensus Groups</h2>\n'

    # Load species composition file
    comp_file = results_dir / f"{organism}_species_composition.csv"
    if not comp_file.exists():
        html += '<p class="alert alert-warning">Species composition analysis not available.</p>\n'
        html += '</div>\n'
        return html

    comp_df = pd.read_csv(comp_file)

    # Summary statistics
    total_groups = len(comp_df)
    multi_species = comp_df['is_multi_species'].sum()
    ambiguous = comp_df['is_ambiguous'].sum()

    html += '<div class="stats-summary">\n'
    html += f'<div class="stat-box"><strong>{total_groups}</strong><br>Total Consensus Groups</div>\n'
    html += f'<div class="stat-box warning"><strong>{multi_species}</strong><br>Multi-Species Groups</div>\n'
    html += f'<div class="stat-box alert"><strong>{ambiguous}</strong><br>Ambiguous Groups (&lt;70% primary)</div>\n'
    html += '</div>\n'

    # Filter options
    html += '<div class="control-group">\n'
    html += '<label>Filter:</label>\n'
    html += '<select id="comp-filter" onchange="filterCompositionTable()">\n'
    html += '<option value="all">All Groups</option>\n'
    html += '<option value="multi">Multi-Species Only</option>\n'
    html += '<option value="ambiguous">Ambiguous Only</option>\n'
    html += '</select>\n'
    html += '</div>\n'

    # Interactive table
    html += '<div class="table-container">\n'
    html += '<table id="composition-table" class="data-table sortable">\n'
    html += '<thead><tr>\n'
    html += '<th onclick="sortTable(0)">Consensus Group ▼</th>\n'
    html += '<th onclick="sortTable(1)">Total Samples</th>\n'
    html += '<th onclick="sortTable(2)"># Species</th>\n'
    html += '<th>Primary Species</th>\n'
    html += '<th>Primary %</th>\n'
    html += '<th>Species Composition</th>\n'
    html += '<th>Flags</th>\n'
    html += '</tr></thead>\n'
    html += '<tbody>\n'

    for _, row in comp_df.iterrows():
        # Apply row highlighting for multi-species/ambiguous
        row_class = ""
        if row['is_ambiguous']:
            row_class = "class='alert-row'"
        elif row['is_multi_species']:
            row_class = "class='warning-row'"

        flags = []
        if row['is_multi_species']:
            flags.append('<span class="badge badge-warning">Multi-species</span>')
        if row['is_ambiguous']:
            flags.append('<span class="badge badge-alert">Ambiguous</span>')
        flags_html = " ".join(flags) if flags else "-"

        html += f'<tr {row_class} data-multi="{row["is_multi_species"]}" data-ambiguous="{row["is_ambiguous"]}">\n'
        html += f'<td><strong>{row["consensus_group"]}</strong></td>\n'
        html += f'<td>{row["n_total_samples"]}</td>\n'
        html += f'<td>{row["n_species"]}</td>\n'
        html += f'<td><em>{row["primary_species"]}</em></td>\n'
        html += f'<td>{row["primary_species_pct"]}</td>\n'
        html += f'<td class="composition-cell">{row["species_composition"]}</td>\n'
        html += f'<td>{flags_html}</td>\n'
        html += '</tr>\n'

    html += '</tbody>\n'
    html += '</table>\n'
    html += '</div>\n'

    html += '</div>\n'
    return html
```

#### 3.2 Add to Navigation Tabs

Update tab navigation to include "Species Composition" tab after "Genotype Assignments".

#### 3.3 Add CSS Styling

```css
.stats-summary {
    display: flex;
    gap: 20px;
    margin-bottom: 30px;
}

.stat-box {
    flex: 1;
    padding: 20px;
    background: #f8f9fa;
    border-radius: 8px;
    text-align: center;
    border-left: 4px solid #007bff;
}

.stat-box.warning {
    border-left-color: #ff9800;
    background: #fff3e0;
}

.stat-box.alert {
    border-left-color: #f44336;
    background: #ffebee;
}

.alert-row {
    background-color: #ffebee !important;
}

.warning-row {
    background-color: #fff3e0 !important;
}

.badge {
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 0.85em;
    font-weight: 600;
}

.badge-warning {
    background: #ff9800;
    color: white;
}

.badge-alert {
    background: #f44336;
    color: white;
}

.composition-cell {
    font-family: 'Courier New', monospace;
    font-size: 0.9em;
}

.data-table {
    width: 100%;
    border-collapse: collapse;
}

.data-table th {
    background: #007bff;
    color: white;
    padding: 12px;
    text-align: left;
    cursor: pointer;
}

.data-table td {
    padding: 10px;
    border-bottom: 1px solid #ddd;
}

.data-table tbody tr:hover {
    background: #e3f2fd;
}
```

#### 3.4 Add JavaScript for Filtering

```javascript
function filterCompositionTable() {
    const filter = document.getElementById('comp-filter').value;
    const table = document.getElementById('composition-table');
    const rows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr');

    for (let row of rows) {
        const isMulti = row.getAttribute('data-multi') === 'True';
        const isAmbiguous = row.getAttribute('data-ambiguous') === 'True';

        let show = false;
        if (filter === 'all') {
            show = true;
        } else if (filter === 'multi' && isMulti) {
            show = true;
        } else if (filter === 'ambiguous' && isAmbiguous) {
            show = true;
        }

        row.style.display = show ? '' : 'none';
    }
}
```

### Phase 4: Additional Enhancements

#### 4.1 Add to Genotype Assignments Table

Modify the existing genotype assignments table in the HTML report to include a "Species Composition" column that shows a summary for each consensus group.

#### 4.2 Integrate with Geographic Distribution

In the distribution map visualization, add ability to:
- Color points by reported species within a consensus group
- Toggle between genotype view and species view
- Highlight geographic clusters of specific species

#### 4.3 Export Functionality

Add buttons to export:
- Filtered species composition table (CSV)
- Species composition summary (JSON for custom analysis)

## Testing Plan

1. **Test with single-species groups**: Verify correct calculation when all samples are same species
2. **Test with multi-species groups**: Verify percentages sum to 100%
3. **Test with missing species data**: Handle NaN/empty species names gracefully
4. **Test with rare species**: Verify min_abundance_pct filtering works correctly
5. **Test HTML rendering**: Verify table sorting, filtering, and styling work properly

## Implementation Order

1. ✅ **Step 1**: Add `compute_species_composition()` function to genotype_assignment.py
2. ✅ **Step 2**: Integrate function call in `assign_genotypes()` and write CSV output
3. ✅ **Step 3**: Add species composition section to HTML report (reports.py)
4. ✅ **Step 4**: Add CSS styling for new section
5. ✅ **Step 5**: Add JavaScript for table filtering
6. ✅ **Step 6**: Test with real data (Rhizoprionodon example)
7. ✅ **Step 7**: Add navigation tab integration
8. ⏸️ **Step 8** (Optional): Integrate with geographic distribution visualization

## File Modifications Required

```
boldgenotyper/
├── genotype_assignment.py
│   └── Add: compute_species_composition() function (after line 823)
│   └── Modify: assign_genotypes() to call function and write CSV (after line 1076)
│
├── reports.py
│   └── Add: _build_species_composition_section() function (around line 2750)
│   └── Add: CSS styling for species composition (in <style> block)
│   └── Add: JavaScript filtering function (in <script> block)
│   └── Modify: Main report builder to include new tab
│
└── config.py (optional)
    └── Add: min_species_abundance_pct parameter to PipelineConfig
```

## Expected Output Example

For your Rhizoprionodon example, the user would see:

**In CSV** (`Rhizoprionodon_species_composition.csv`):
```
consensus_group,n_total_samples,n_species,species_list,primary_species,primary_species_pct,species_composition,is_multi_species,is_ambiguous
Rhizoprionodon c13_n73,73,2,"Rhizoprionodon porosus, Rhizoprionodon terranovae",Rhizoprionodon porosus,67.1%,"Rhizoprionodon porosus: 49 (67.1%); Rhizoprionodon terranovae: 24 (32.9%)",True,True
```

**In HTML Report**:
- Species Composition tab shows this group highlighted in amber/red
- Clear indication this is an ambiguous multi-species group
- Easy to see 67.1% reported as *R. porosus*, 32.9% as *R. terranovae*
- User can cross-reference with geographic distribution to investigate hypothesis about misidentifications

## Benefits

1. **Transparency**: Users can see exactly which species contributed to each consensus group
2. **Quality Control**: Quickly identify potentially problematic consensus groups
3. **Research Insights**: Support hypotheses about misidentifications or cryptic species
4. **Decision Support**: Help users decide whether to accept genus-level assignments or investigate further
5. **Documentation**: Automated tracking for reproducibility and publication

## Future Enhancements

- Add BOLD link-outs for each species in composition
- Statistical tests for geographic clustering by species within consensus groups
- Integration with external taxonomy databases (WoRMS, NCBI) for validation
- Automated flagging of known problematic species complexes
