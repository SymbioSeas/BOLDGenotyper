#!/usr/bin/env python
"""Test Priority 1 features with your own BOLD data."""

import pandas as pd
import sys
from boldgenotyper.metadata import parse_bold_tsv, parse_coordinates_column
from boldgenotyper.geographic import classify_and_snap_coordinates, load_land_geometries

# Configuration
TSV_FILE = "tests/data/Crassostrea.tsv"  # Change to your file name
COORD_COLUMN = "coord"     # BOLD coordinate column name (usually "coord")

print("="*70)
print("Testing Priority 1 Features with Your Data")
print("="*70)

# Step 1: Load your data
print(f"\nStep 1: Loading {TSV_FILE}...")
try:
    df = parse_bold_tsv(TSV_FILE)
    print(f"  ✓ Loaded {len(df)} samples")
    print(f"  ✓ Columns: {len(df.columns)}")
    
    # Show sample of data
    if 'processid' in df.columns:
        print(f"  ✓ First processid: {df['processid'].iloc[0]}")
    if COORD_COLUMN in df.columns:
        print(f"  ✓ First coordinate: {df[COORD_COLUMN].iloc[0]}")
except Exception as e:
    print(f"  ✗ Error loading file: {e}")
    sys.exit(1)

# Step 2: Parse coordinates with automatic order inference
print(f"\nStep 2: Parsing coordinates with automatic order inference...")
try:
    df = parse_coordinates_column(
        df,
        coord_column=COORD_COLUMN,
        infer_order=True  # Automatically detect lat/lon vs lon/lat
    )
    
    n_valid = df[['latitude', 'longitude']].notna().all(axis=1).sum()
    print(f"  ✓ Parsed {n_valid}/{len(df)} coordinates successfully")
    print(f"  ✓ Latitude range: {df['latitude'].min():.2f} to {df['latitude'].max():.2f}")
    print(f"  ✓ Longitude range: {df['longitude'].min():.2f} to {df['longitude'].max():.2f}")
    
    # Check if ranges make sense
    if df['latitude'].abs().max() > 90:
        print("  ⚠ WARNING: Latitude values > 90° detected - coordinates may be swapped!")
    if df['longitude'].abs().max() > 180:
        print("  ⚠ WARNING: Longitude values > 180° detected - may need normalization")
        
except Exception as e:
    print(f"  ✗ Error parsing coordinates: {e}")
    sys.exit(1)

# Step 3: Classify coordinates (mark policy first to inspect)
print(f"\nStep 3: Classifying coordinates as land/ocean...")
print("  (First run will download Natural Earth data ~1-10 MB)")

try:
    # Load land geometry once
    land = load_land_geometries(resolution="110m")
    
    # Classify with 'mark' policy to inspect
    df_marked = classify_and_snap_coordinates(
        df,
        land_geometry=land,
        lat_col='latitude',
        lon_col='longitude',
        policy='mark'  # Just mark, don't modify yet
    )
    
    # Show classification results
    location_counts = df_marked['location_type'].value_counts()
    print(f"\n  Location type distribution:")
    for loc_type, count in location_counts.items():
        pct = count / len(df_marked) * 100
        print(f"    {loc_type}: {count} ({pct:.1f}%)")
    
    # Identify land points
    land_samples = df_marked[df_marked['location_type'] == 'country_centroid']
    n_land = len(land_samples)
    
    if n_land > 0:
        print(f"\n  ✓ Detected {n_land} samples with land coordinates")
        print(f"    These are likely country/region centroids from BOLD")
        print(f"\n  Sample land coordinates:")
        cols_to_show = ['latitude', 'longitude']
        if 'country' in land_samples.columns:
            cols_to_show.insert(0, 'country')
        if 'processid' in land_samples.columns:
            cols_to_show.insert(0, 'processid')
        print(land_samples[cols_to_show].head(5))
    else:
        print(f"\n  ✓ No land coordinates detected - all samples in ocean")
        
except Exception as e:
    print(f"  ✗ Error classifying coordinates: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 4: Snap to coast (if land points exist)
if n_land > 0:
    print(f"\nStep 4: Snapping {n_land} land coordinates to nearest coastline...")
    
    try:
        df_snapped = classify_and_snap_coordinates(
            df,
            land_geometry=land,  # Reuse loaded geometry
            lat_col='latitude',
            lon_col='longitude',
            policy='snap'  # Snap land points to coast
        )
        
        # Show results
        location_counts_after = df_snapped['location_type'].value_counts()
        print(f"\n  Location type distribution after snapping:")
        for loc_type, count in location_counts_after.items():
            pct = count / len(df_snapped) * 100
            print(f"    {loc_type}: {count} ({pct:.1f}%)")
        
        # Calculate snap distances
        snapped_samples = df_snapped[df_snapped['location_type'] == 'snapped_to_coast']
        
        if len(snapped_samples) > 0:
            print(f"\n  ✓ Successfully snapped {len(snapped_samples)} samples to coast")
            
            # Show a few examples
            print(f"\n  Example snapped coordinates:")
            for idx in snapped_samples.index[:3]:  # First 3 examples
                orig_lat = df_marked.loc[idx, 'latitude']
                orig_lon = df_marked.loc[idx, 'longitude']
                new_lat = df_snapped.loc[idx, 'latitude']
                new_lon = df_snapped.loc[idx, 'longitude']
                
                # Rough distance
                dist_deg = ((new_lat - orig_lat)**2 + (new_lon - orig_lon)**2)**0.5
                dist_km = dist_deg * 111
                
                sample_id = df_snapped.loc[idx, 'processid'] if 'processid' in df_snapped.columns else idx
                print(f"    {sample_id}:")
                print(f"      Before: [{orig_lat:.3f}°, {orig_lon:.3f}°]")
                print(f"      After:  [{new_lat:.3f}°, {new_lon:.3f}°]")
                print(f"      Snap distance: ~{dist_km:.0f} km")
        
        # Save results
        output_file = TSV_FILE.replace('.tsv', '_processed.tsv')
        df_snapped.to_csv(output_file, sep='\t', index=False)
        print(f"\n  ✓ Saved processed data to: {output_file}")
        
    except Exception as e:
        print(f"  ✗ Error snapping coordinates: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
else:
    print(f"\nStep 4: No land coordinates to snap - skipping")
    
    # Save marked data
    output_file = TSV_FILE.replace('.tsv', '_processed.tsv')
    df_marked.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved processed data to: {output_file}")

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"Total samples: {len(df)}")
print(f"Valid coordinates: {n_valid}")
if n_land > 0:
    print(f"Land coordinates detected: {n_land} ({n_land/len(df)*100:.1f}%)")
    print(f"Snapped to coast: {len(snapped_samples) if 'snapped_samples' in locals() else 0}")
else:
    print(f"Land coordinates detected: 0 (all samples in ocean)")

print(f"\nProcessed data saved to: {output_file}")
print("\n✅ All tests completed successfully!")
print("="*70)