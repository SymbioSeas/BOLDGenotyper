#!/usr/bin/env python3
"""
Quick test to verify consensus_group_sp labels are correctly formatted.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd

# Test with Euprymna data if available
euprymna_annotated = Path("tests/Euprymna_output/Euprymna_annotated.csv")

if not euprymna_annotated.exists():
    print(f"❌ Test data not found: {euprymna_annotated}")
    print("Please run the pipeline first to generate test data.")
    sys.exit(1)

print("=" * 70)
print("Testing consensus_group_sp Label Formation")
print("=" * 70)

# Load data
df = pd.read_csv(euprymna_annotated)

# Check what columns we have
print(f"\nDataframe columns: {list(df.columns)}")
print(f"Total samples: {len(df)}")

# Analyze consensus_group_sp values
if 'consensus_group_sp' in df.columns:
    print(f"\nAnalyzing consensus_group_sp values:")
    print("-" * 70)

    assigned = df[df['consensus_group_sp'].notna()]
    print(f"Samples with consensus_group_sp: {len(assigned)}")

    # Check label format
    labels = assigned['consensus_group_sp'].unique()
    print(f"\nUnique consensus_group_sp labels: {len(labels)}")

    # Categorize labels
    with_species = []
    without_species = []

    for label in labels:
        # If it starts with "c" followed by digits, it's missing the species
        if str(label).strip().startswith('c') and str(label)[1].isdigit():
            without_species.append(label)
        else:
            with_species.append(label)

    print(f"\n✓ Labels WITH species name: {len(with_species)}")
    for label in sorted(with_species)[:10]:
        print(f"  • {label}")
    if len(with_species) > 10:
        print(f"  ... and {len(with_species) - 10} more")

    print(f"\n❌ Labels WITHOUT species name (PROBLEMATIC): {len(without_species)}")
    for label in sorted(without_species)[:10]:
        # Find the corresponding assigned_sp value
        sample_row = assigned[assigned['consensus_group_sp'] == label].iloc[0]
        assigned_sp = sample_row.get('assigned_sp', 'N/A')
        majority_frac = sample_row.get('majority_fraction', 'N/A')
        print(f"  • {label:20s} (assigned_sp={assigned_sp}, majority_frac={majority_frac})")
    if len(without_species) > 10:
        print(f"  ... and {len(without_species) - 10} more")

    # Summary
    print(f"\n{'=' * 70}")
    if len(without_species) == 0:
        print("✅ SUCCESS: All labels have species names!")
    else:
        print(f"⚠️  ISSUE: {len(without_species)} labels are missing species names")
        print(f"   This might be expected if:")
        print(f"   • majority_fraction < 0.70 (not enough consensus)")
        print(f"   • No cluster_seq data available")
        print(f"   • assigned_sp is genuinely empty (no species data)")
    print("=" * 70)

else:
    print("❌ Column 'consensus_group_sp' not found in dataframe!")
