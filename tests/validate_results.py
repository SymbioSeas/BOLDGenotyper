#!/usr/bin/env python3
"""
Validate pipeline results and create summary visualizations.

Usage:
    python tests/validate_results.py
"""

import pandas as pd
from pathlib import Path
import sys

OUTPUT_DIR = Path("tests/output")

def validate_results():
    """Validate pipeline outputs."""
    
    print("="*70)
    print("  BOLDGenotyper Results Validation")
    print("="*70)
    
    errors = []
    warnings = []
    
    # Check files exist
    print("\n1. Checking output files...")
    required_files = [
        "01_parsed_metadata.tsv",
        "02_filtered_metadata.tsv",
    ]
    
    for filename in required_files:
        filepath = OUTPUT_DIR / filename
        if filepath.exists():
            print(f"  ✓ {filename}")
        else:
            errors.append(f"Missing required file: {filename}")
            print(f"  ✗ {filename} - MISSING")
    
    # Load main results
    print("\n2. Loading results...")
    try:
        # Find the genotypes file
        genotype_files = list(OUTPUT_DIR.glob("*_with_genotypes.tsv"))
        if not genotype_files:
            errors.append("No genotype assignment file found")
            print("  ✗ No genotypes file found")
            return errors, warnings
        
        genotypes_file = genotype_files[0]
        df = pd.read_csv(genotypes_file, sep='\t')
        print(f"  ✓ Loaded {len(df)} samples from {genotypes_file.name}")
        
    except Exception as e:
        errors.append(f"Failed to load results: {e}")
        return errors, warnings
    
    # Validate data quality
    print("\n3. Validating data quality...")
    
    # Check for required columns
    required_cols = ['processid', 'consensus_group']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        errors.append(f"Missing columns: {missing_cols}")
    else:
        print(f"  ✓ All required columns present")
    
    # Check genotype assignments
    if 'consensus_group' in df.columns:
        assigned = df['consensus_group'].notna().sum()
        unassigned = df['consensus_group'].isna().sum()
        
        print(f"\n  Genotype assignments:")
        print(f"    Assigned: {assigned} ({assigned/len(df)*100:.1f}%)")
        print(f"    Unassigned: {unassigned} ({unassigned/len(df)*100:.1f}%)")
        
        if unassigned > len(df) * 0.3:
            warnings.append(f"High unassignment rate: {unassigned/len(df)*100:.1f}%")
        
        # Show genotype distribution
        if assigned > 0:
            genotype_counts = df['consensus_group'].value_counts()
            print(f"\n  Unique genotypes: {len(genotype_counts)}")
            print(f"  Genotype distribution:")
            for genotype, count in genotype_counts.head(10).items():
                if pd.notna(genotype):
                    print(f"    {genotype}: {count} samples")
    
    # Check ocean basin assignments
    if 'ocean_basin' in df.columns:
        basins = df['ocean_basin'].value_counts()
        print(f"\n  Ocean basins: {len(basins)} distinct")
        if len(basins) == 0:
            warnings.append("No ocean basin assignments found")
        else:
            print(f"  Basin distribution:")
            for basin, count in basins.head(5).items():
                print(f"    {basin}: {count} samples")
    
    # Check diagnostics file
    print("\n4. Checking diagnostics...")
    diag_files = list(OUTPUT_DIR.glob("*_diagnostics.csv"))
    if diag_files:
        diag = pd.read_csv(diag_files[0])
        print(f"  ✓ Found diagnostics ({len(diag)} records)")
        
        if 'identity' in diag.columns:
            mean_identity = diag['identity'].mean()
            min_identity = diag['identity'].min()
            
            print(f"    Mean identity: {mean_identity:.3f}")
            print(f"    Min identity: {min_identity:.3f}")
            
            if min_identity < 0.80:
                warnings.append(f"Low minimum identity: {min_identity:.3f}")
            
            # Check for ties
            if 'runner_up_identity' in diag.columns:
                diag['diff'] = diag['identity'] - diag['runner_up_identity']
                close_calls = (diag['diff'] < 0.02).sum()
                if close_calls > 0:
                    print(f"    ⚠ {close_calls} samples with ambiguous assignments (close ties)")
    else:
        warnings.append("No diagnostics file found")
    
    # Check consensus sequences
    print("\n5. Checking consensus sequences...")
    fasta_files = list(OUTPUT_DIR.glob("*_consensus.fasta"))
    if fasta_files:
        with open(fasta_files[0]) as f:
            consensus_seqs = [line for line in f if line.startswith('>')]
        print(f"  ✓ Found {len(consensus_seqs)} consensus sequences")
        
        # Show consensus IDs
        print(f"    Consensus genotypes:")
        for header in consensus_seqs[:10]:
            print(f"      {header.strip()[1:]}")  # Remove '>' and newline
    else:
        errors.append("No consensus FASTA file found")
    
    # Summary
    print("\n" + "="*70)
    print("  Validation Summary")
    print("="*70)
    
    if errors:
        print(f"\n✗ ERRORS ({len(errors)}):")
        for i, error in enumerate(errors, 1):
            print(f"  {i}. {error}")
    else:
        print("\n✓ No errors found")
    
    if warnings:
        print(f"\n⚠ WARNINGS ({len(warnings)}):")
        for i, warning in enumerate(warnings, 1):
            print(f"  {i}. {warning}")
    else:
        print("\n✓ No warnings")
    
    if not errors and not warnings:
        print("\n🎉 All validations passed! Your pipeline is working correctly.")
    
    print("\n" + "="*70)
    
    return errors, warnings


if __name__ == '__main__':
    errors, warnings = validate_results()
    
    # Exit with error code if there are errors
    sys.exit(1 if errors else 0)
