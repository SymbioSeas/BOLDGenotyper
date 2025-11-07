#!/usr/bin/env python
"""Quick test of snap-to-coast functionality."""

import pandas as pd
from boldgenotyper.geographic import classify_and_snap_coordinates

print("Testing snap-to-coast...\n")

# Test data with known inland point
df = pd.DataFrame({
    'sample': ['CentralAustralia', 'Pacific', 'CoastalPoint'],
    'latitude': [-25.0, 0.0, -33.0],
    'longitude': [133.0, -150.0, 151.0]  # Central Australia is ~400km inland
})

print("Original coordinates:")
print(df[['sample', 'latitude', 'longitude']])

print("\nSnapping land points to coast...")
df_snapped = classify_and_snap_coordinates(
    df,
    policy='snap',  # Snap land points to coast
    resolution='110m'
)

print("\nAfter snapping:")
print(df_snapped[['sample', 'latitude', 'longitude', 'location_type']])

# Calculate distances for snapped points
snapped_samples = df_snapped[df_snapped['location_type'] == 'snapped_to_coast']

if len(snapped_samples) > 0:
    print(f"\n✓ Snapped {len(snapped_samples)} sample(s) to coast:")
    
    for idx in snapped_samples.index:
        orig_lat = df.loc[idx, 'latitude']
        orig_lon = df.loc[idx, 'longitude']
        new_lat = df_snapped.loc[idx, 'latitude']
        new_lon = df_snapped.loc[idx, 'longitude']
        
        print(f"  {df.loc[idx, 'sample']}:")
        print(f"    Before: [{orig_lat:.2f}°, {orig_lon:.2f}°]")
        print(f"    After:  [{new_lat:.2f}°, {new_lon:.2f}°]")
        
        # Rough distance calculation (degrees)
        dist_deg = ((new_lat - orig_lat)**2 + (new_lon - orig_lon)**2)**0.5
        dist_km = dist_deg * 111  # Rough conversion
        print(f"    Distance: ~{dist_km:.0f} km")
else:
    print("\n⚠ No land points detected (resolution might be too coarse)")

print("\n" + "="*60)