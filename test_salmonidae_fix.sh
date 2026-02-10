#!/bin/bash
# Test custom shapefile coordinate parsing fix

echo "Testing Salmonidae with custom shapefile coordinate parsing..."
echo "=============================================="

# Remove old output
rm -rf output/Salmonidae_full_test

# Run with fixed code
boldgenotyper data/Salmonidae.tsv \
  --output output/Salmonidae_full_test \
  --organism Salmonidae \
  --build-tree \
  --custom-shp hybas_pour_lev07_v1_shp/hybas_pour_lev07_v1.shp \
  --shp-field HYBAS_ID \
  --geo-category freshwater_basin \
  --threads 8

echo ""
echo "=============================================="
echo "Checking results..."
echo "=============================================="

# Check coordinates
python3 << 'EOF'
import pandas as pd

df = pd.read_csv('output/Salmonidae_full_test/Salmonidae_annotated.csv', low_memory=False)
print(f"\nTotal samples: {len(df)}")
print(f"Samples with lat: {df['lat'].notna().sum()} ({df['lat'].notna().sum()/len(df)*100:.1f}%)")
print(f"Samples with lon: {df['lon'].notna().sum()} ({df['lon'].notna().sum()/len(df)*100:.1f}%)")
print(f"Samples with freshwater_basin: {df['freshwater_basin'].notna().sum()} ({df['freshwater_basin'].notna().sum()/len(df)*100:.1f}%)")
print(f"\nUnique basins: {df['freshwater_basin'].value_counts().head(10)}")
print(f"\nSample data:")
print(df[['processid', 'lat', 'lon', 'freshwater_basin']].head(10))
EOF

# Check visualizations
echo ""
echo "Visualization directory contents:"
ls -lh output/Salmonidae_full_test/visualization/

echo ""
echo "Done!"
