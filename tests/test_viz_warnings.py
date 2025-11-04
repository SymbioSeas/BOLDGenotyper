#!/usr/bin/env python3
"""
Quick test to verify distribution map warning messages for Carcharhinus data.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from boldgenotyper import visualization, utils
import pandas as pd

# Setup logging
logger = utils.setup_logging(log_level="INFO")

print("=" * 70)
print("Testing Distribution Map Warning Messages")
print("=" * 70)

# Load the annotated Carcharhinus data
annotated_file = Path("tests/Carcharhinus_output/Carcharhinus_annotated.csv")
print(f"\nLoading data from: {annotated_file}")

df = pd.read_csv(annotated_file)
print(f"Loaded {len(df)} samples")

# Show data availability
total = len(df)
with_coords = (df['lat'].notna() & df['lon'].notna()).sum()
with_genotype = df['consensus_group_sp'].notna().sum()
with_both = ((df['lat'].notna() & df['lon'].notna()) & df['consensus_group_sp'].notna()).sum()

print(f"\nData availability:")
print(f"  Total samples: {total}")
print(f"  With coordinates: {with_coords} ({with_coords/total*100:.1f}%)")
print(f"  With genotype: {with_genotype} ({with_genotype/total*100:.1f}%)")
print(f"  With BOTH: {with_both} ({with_both/total*100:.1f}%)")

# Test plot_distribution_map (should trigger warning)
print(f"\n{'=' * 70}")
print("Test 1: plot_distribution_map()")
print("=" * 70)

output_path = Path("tests/Carcharhinus_output/test_distribution_map.png")
visualization.plot_distribution_map(
    df=df,
    output_path=output_path,
    genotype_column='consensus_group_sp',
    latitude_col='lat',
    longitude_col='lon'
)

# Test plot_distribution_map_faceted (should trigger warning)
print(f"\n{'=' * 70}")
print("Test 2: plot_distribution_map_faceted()")
print("=" * 70)

output_path_faceted = Path("tests/Carcharhinus_output/test_distribution_map_faceted.png")
visualization.plot_distribution_map_faceted(
    df=df,
    output_path=output_path_faceted,
    genotype_column='consensus_group_sp',
    species_column='species',
    latitude_col='lat',
    longitude_col='lon'
)

print(f"\n{'=' * 70}")
print("✓ Test completed")
print("=" * 70)
