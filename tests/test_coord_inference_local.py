#!/usr/bin/env python
"""Quick test of coordinate order inference."""

import pandas as pd
from boldgenotyper.metadata import infer_coord_order, parse_coordinates_column

# Test data with known lat/lon order
df_latlon = pd.DataFrame({
    'sample': ['S1', 'S2', 'S3'],
    'coord': ['[34.5, -76.2]', '[-33.9, 151.2]', '[41.2, -8.9]']
})

# Test data with lon/lat order (swapped)
df_lonlat = pd.DataFrame({
    'sample': ['S1', 'S2', 'S3'],
    'coord': ['[-76.2, 34.5]', '[151.2, -33.9]', '[-8.9, 41.2]']
})

print("Testing coordinate order inference...\n")

# Test 1: lat/lon format
print("Test 1: Known lat/lon format")
df1 = parse_coordinates_column(df_latlon, coord_column='coord', infer_order=True)
print(f"  Detected order and parsed coordinates:")
print(df1[['latitude', 'longitude']].head())
print(f"  ✓ First sample: lat={df1['latitude'].iloc[0]}, lon={df1['longitude'].iloc[0]}")

# Test 2: lon/lat format (should auto-detect and swap)
print("\nTest 2: Known lon/lat format (should auto-detect and swap)")
df2 = parse_coordinates_column(df_lonlat, coord_column='coord', infer_order=True)
print(f"  Detected order and parsed coordinates:")
print(df2[['latitude', 'longitude']].head())
print(f"  ✓ First sample: lat={df2['latitude'].iloc[0]}, lon={df2['longitude'].iloc[0]}")

# Both should result in same coordinates (after auto-detection)
if abs(df1['latitude'].iloc[0] - df2['latitude'].iloc[0]) < 0.01:
    print("\n✅ SUCCESS: Auto-detection correctly handled both formats!")
else:
    print("\n⚠ WARNING: Results differ - check implementation")

print("\n" + "="*60)