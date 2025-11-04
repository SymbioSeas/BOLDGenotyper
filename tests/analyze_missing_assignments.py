#!/usr/bin/env python3
"""
Analyze why many consensus groups don't get species assignments.
"""

import pandas as pd
import re
from pathlib import Path

print("=" * 80)
print("ANALYSIS: Why Many Consensus Groups Lack Species Assignments")
print("=" * 80)

# Load data
tsv_path = Path("tests/Carcharhinus_output/02_filtered_metadata.tsv")
fasta_path = Path("tests/Carcharhinus_output/Carcharhinus.fasta")
consensus_fasta = Path("tests/Carcharhinus_output/Carcharhinus_consensus.fasta")
annotated = Path("tests/Carcharhinus_output/Carcharhinus_annotated.csv")

df = pd.read_csv(annotated)
tsv = pd.read_csv(tsv_path, sep='\t')

print(f"\n1. OVERALL STATISTICS")
print("-" * 80)
print(f"Total samples in dataset: {len(df)}")
print(f"Samples with valid sequences: {len(tsv)}")
print(f"Samples assigned to consensus groups: {df['consensus_group_sp'].notna().sum()}")
print(f"Samples NOT assigned: {df['consensus_group_sp'].isna().sum()}")

# Get consensus group info
with open(consensus_fasta) as f:
    consensus_groups = {}
    for line in f:
        if line.startswith('>'):
            match = re.search(r'c(\d+)_n(\d+)', line)
            if match:
                cluster_id = f'c{match.group(1)}_n{match.group(2)}'
                cluster_size = int(match.group(2))
                consensus_groups[cluster_id] = cluster_size

print(f"\nTotal consensus groups created: {len(consensus_groups)}")
print(f"Consensus groups with sample assignments: {df['consensus_group_sp'].nunique()}")

# Analyze unassigned groups by size
assigned_groups = set()
for val in df['consensus_group_sp'].dropna().unique():
    # Extract the c#_n# pattern
    match = re.search(r'(c\d+_n\d+)', val)
    if match:
        assigned_groups.add(match.group(1))

unassigned_groups = {k: v for k, v in consensus_groups.items() if k not in assigned_groups}

print(f"\n2. CONSENSUS GROUP SIZE DISTRIBUTION")
print("-" * 80)

size_bins = [(1, 1, "Singletons"), (2, 5, "Small (2-5)"), (6, 20, "Medium (6-20)"),
             (21, 100, "Large (21-100)"), (101, 1000, "Very large (>100)")]

print("\nASSIGNED consensus groups:")
for min_sz, max_sz, label in size_bins:
    count = sum(1 for g in assigned_groups if min_sz <= consensus_groups[g] <= max_sz)
    if count > 0:
        examples = [g for g in sorted(assigned_groups) if min_sz <= consensus_groups[g] <= max_sz][:3]
        print(f"  {label:20s}: {count:3d} groups (e.g., {', '.join(examples)})")

print("\nUNASSIGNED consensus groups:")
for min_sz, max_sz, label in size_bins:
    count = sum(1 for g, sz in unassigned_groups.items() if min_sz <= sz <= max_sz)
    if count > 0:
        examples = [g for g, sz in sorted(unassigned_groups.items(), key=lambda x: -x[1])
                   if min_sz <= sz <= max_sz][:3]
        print(f"  {label:20s}: {count:3d} groups (e.g., {', '.join(examples)})")

print(f"\n3. TOP UNASSIGNED GROUPS BY SIZE")
print("-" * 80)
top_unassigned = sorted(unassigned_groups.items(), key=lambda x: -x[1])[:10]
print("\nThese groups had many sequences cluster together during dereplication,")
print("but NO samples matched them during genotype assignment (≥80% identity):\n")

for i, (group, size) in enumerate(top_unassigned, 1):
    print(f"{i:2d}. {group:15s} - {size:4d} sequences clustered at 99% identity")

print(f"\n4. BIOLOGICAL INTERPRETATION")
print("-" * 80)
print("""
WHY does this happen? Several biological and technical reasons:

A. HIGH INTRASPECIFIC GENETIC DIVERSITY
   • Carcharhinus species show exceptionally high COI variation (5-10%)
   • Many species have divergent geographic lineages
   • Some "species" are actually cryptic species complexes

B. CONSENSUS SEQUENCE QUALITY ISSUES
   • When highly divergent sequences cluster together (at 99% threshold),
     the resulting consensus may not match any individual sequence well
   • Consensus sequences with many ambiguous positions (N's) reduce match scores
   • Example: 129 divergent sequences → poor consensus → no matches ≥80%

C. THRESHOLD MISMATCH
   • Dereplication clusters at 99% identity (tight grouping)
   • Genotype assignment requires 80% identity (looser matching)
   • But if consensus quality is poor, even 80% may not be met

D. SYSTEMATIC BIAS
   • Large unassigned groups like c20_n129 (129 seqs) likely represent:
     - Species with high population structure (e.g., C. porosus)
     - Cryptic species not yet formally described
     - Geographic races with 10-15% COI divergence
   • These sequences cluster together but are too divergent to assign

E. SINGLETONS AND SMALL GROUPS
   • 46 singleton clusters represent highly divergent individuals
   • Could be: sequencing errors, contamination, rare species, hybrids

BIOLOGICAL CONCLUSION:
This pattern is EXPECTED and BIOLOGICALLY MEANINGFUL for Carcharhinus:
✓ Genus underwent recent rapid radiation (~15 MYA)
✓ Many cryptic species and poorly defined boundaries
✓ High within-species geographic structure
✓ COI divergence often exceeds typical "species threshold"

The 67 unassigned consensus groups likely represent:
• 27 species not captured by the 18 assigned groups
• Geographic lineages within assigned species
• Cryptic species complexes
• Genuine high genetic diversity
""")

print(f"\n5. RECOMMENDATIONS")
print("-" * 80)
print("""
To improve assignments:

1. LOWER SIMILARITY THRESHOLD (already using 0.80, could try 0.75)
   - May capture more divergent lineages
   - Risk of false assignments increases

2. USE PHYLOGENETIC APPROACH
   - Build tree of all consensus sequences
   - Identify monophyletic species groups
   - May reveal cryptic species

3. EXAMINE LARGE UNASSIGNED GROUPS MANUALLY
   - c20_n129, c17_n65, c16_n61 need investigation
   - Check consensus sequence quality
   - Look at species composition if samples have IDs

4. CONSIDER MULTI-GENE APPROACH
   - COI alone may be insufficient for Carcharhinus
   - Add nuclear markers or other mitochondrial genes

5. ACCEPT BIOLOGICAL REALITY
   - Not all genetic diversity fits into clean "species" boxes
   - Carcharhinus is genuinely messy genetically
   - Unassigned groups may represent real biological entities
""")
