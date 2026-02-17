# scripts/validate_recommended_thresholds.py
import pandas as pd

# Load summary
df = pd.read_csv("benchmarking/parameter_sweep/sweep_summary_all_datasets.csv")

# Add literature comparison (you'll need to research these)
literature_thresholds = {
    'Sphyrnidae': 0.025,  # From your manual analysis or literature
    'Panulirus': 0.020,   # From Silberman et al. or similar
    'Salmonidae': 0.015,  # From salmonid phylogeography papers
    'Pieridae': 0.030     # From butterfly DNA barcoding literature
    'Carcharhiniformes': 0.015, # From XXX
}

df['literature_threshold'] = df['dataset'].map(literature_thresholds)
df['threshold_difference'] = abs(df['recommended_threshold'] - df['literature_threshold'])
df['threshold_ratio'] = df['recommended_threshold'] / df['literature_threshold']

print("Recommended vs. Literature Thresholds:")
print(df[['dataset', 'recommended_threshold', 'literature_threshold',
         'threshold_difference', 'threshold_ratio']].to_string(index=False))

# Check if recommendations are reasonable (within 50% of literature values)
reasonable = (df['threshold_ratio'] >= 0.5) & (df['threshold_ratio'] <= 1.5)
print(f"\n{reasonable.sum()}/{len(df)} datasets have reasonable threshold recommendations")

df.to_csv("benchmarking/parameter_sweep/threshold_validation.csv", index=False)
