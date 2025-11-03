#!/usr/bin/env python3
"""
Manual pipeline test script for BOLDGenotyper.

This script runs the complete pipeline step-by-step on test data,
with checkpoints and validation at each stage.

Usage:
    python tests/test_pipeline_manual.py

Author: Steph Smith
"""

import sys
from pathlib import Path

# Add package to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from boldgenotyper import utils, config, metadata, geographic
from boldgenotyper import dereplication, genotype_assignment
import pandas as pd
import logging

# Setup logging
logger = utils.setup_logging(log_level="INFO", log_file="tests/test_run.log")

# Paths
TEST_DATA_DIR = Path("tests/data")
TEST_TSV = TEST_DATA_DIR / "Sphyrnidae_test.tsv"
OUTPUT_DIR = Path("tests/output")
OUTPUT_DIR.mkdir(exist_ok=True)

def print_section(title):
    """Print a section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)

def check_file_exists(filepath, description):
    """Check if file exists and print status."""
    if filepath.exists():
        size = utils.format_file_size(filepath.stat().st_size)
        print(f"✓ {description}: {filepath} ({size})")
        return True
    else:
        print(f"✗ {description} not found: {filepath}")
        return False

def print_dataframe_info(df, name):
    """Print DataFrame information."""
    print(f"\n{name}:")
    print(f"  Rows: {len(df)}")
    print(f"  Columns: {len(df.columns)}")
    print(f"  Columns: {', '.join(df.columns[:10])}" + ("..." if len(df.columns) > 10 else ""))
    if len(df) > 0:
        print(f"  First row sample:")
        for col in ['processid', 'species', 'coord', 'nuc'][:4]:
            if col in df.columns:
                val = df[col].iloc[0]
                if isinstance(val, str) and len(val) > 50:
                    val = val[:50] + "..."
                print(f"    {col}: {val}")

def main():
    """Run the complete pipeline test."""
    
    print_section("BOLDGenotyper Pipeline Test")
    print(f"Test data: {TEST_TSV}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # ========================================================================
    # STEP 1: Load Configuration
    # ========================================================================
    print_section("STEP 1: Load Configuration")
    
    try:
        cfg = config.get_default_config()
        cfg = cfg.update(
            output_dir=OUTPUT_DIR,
            n_threads=2,
            keep_intermediates=True,
            log_level="DEBUG"
        )
        print("✓ Configuration loaded")
        print(f"  Clustering threshold: {cfg.dereplication.clustering_threshold}")
        print(f"  Min identity: {cfg.genotype_assignment.min_identity}")
        print(f"  Exclude centroids: {cfg.geographic.exclude_centroids}")
        
        # Validate configuration
        warnings = config.validate_config(cfg)
        if warnings:
            print(f"\n⚠ Configuration warnings ({len(warnings)}):")
            for warning in warnings[:3]:
                print(f"  - {warning}")
        
    except Exception as e:
        print(f"✗ Configuration failed: {e}")
        return False
    
    # ========================================================================
    # STEP 2: Parse Metadata
    # ========================================================================
    print_section("STEP 2: Parse BOLD Metadata")
    
    try:
        # Check input file
        if not check_file_exists(TEST_TSV, "Input TSV"):
            return False
        
        # Parse TSV
        print("\nParsing TSV...")
        df = metadata.parse_bold_tsv(TEST_TSV)
        print_dataframe_info(df, "Parsed data")
        
        # Validate required columns
        required = ['processid', 'nuc', 'coord']
        missing = [col for col in required if col not in df.columns]
        if missing:
            print(f"✗ Missing required columns: {missing}")
            return False
        print(f"✓ All required columns present")
        
        # Save parsed data
        parsed_tsv = OUTPUT_DIR / "01_parsed_metadata.tsv"
        df.to_csv(parsed_tsv, sep='\t', index=False)
        print(f"✓ Saved parsed data: {parsed_tsv}")
        
    except Exception as e:
        print(f"✗ Metadata parsing failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # ========================================================================
    # STEP 3: Filter Coordinates
    # ========================================================================
    print_section("STEP 3: Filter Coordinates")
    
    try:
        print(f"Samples before filtering: {len(df)}")
        
        # Apply coordinate filters
        df_filtered = metadata.filter_by_coordinate_quality(df, cfg.geographic)
        
        print(f"Samples after filtering: {len(df_filtered)}")
        print(f"Removed: {len(df) - len(df_filtered)} samples")
        
        if len(df_filtered) == 0:
            print("✗ No samples remaining after filtering!")
            print("\nCheck your test data for:")
            print("  - Valid coordinates (not empty)")
            print("  - Non-centroid coordinates")
            print("  - Non-zero coordinates")
            return False
        
        # Parse coordinates
        print("\nExtracting lat/lon coordinates...")
        coords_data = []
        for _, row in df_filtered.iterrows():
            try:
                lat, lon = metadata.extract_coordinates(row['coord'])
                coords_data.append({
                    'processid': row['processid'],
                    'latitude': lat,
                    'longitude': lon
                })
            except Exception as e:
                print(f"  Warning: Could not parse coord for {row['processid']}: {e}")
        
        df_coords = pd.DataFrame(coords_data)
        print(f"✓ Parsed coordinates for {len(df_coords)} samples")
        
        # Show coordinate range
        if len(df_coords) > 0:
            print(f"\nCoordinate ranges:")
            print(f"  Latitude: {df_coords['latitude'].min():.2f} to {df_coords['latitude'].max():.2f}")
            print(f"  Longitude: {df_coords['longitude'].min():.2f} to {df_coords['longitude'].max():.2f}")
        
        # Save filtered data
        filtered_tsv = OUTPUT_DIR / "02_filtered_metadata.tsv"
        df_filtered.to_csv(filtered_tsv, sep='\t', index=False)
        print(f"\n✓ Saved filtered data: {filtered_tsv}")
        
    except Exception as e:
        print(f"✗ Coordinate filtering failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # ========================================================================
    # STEP 4: Assign Ocean Basins
    # ========================================================================
    print_section("STEP 4: Assign Ocean Basins")
    
    try:
        # Check GOaS data
        goas_path = cfg.geographic.goas_shapefile_path
        if not goas_path.exists():
            print(f"✗ GOaS shapefile not found: {goas_path}")
            print("\nTo download GOaS data, run:")
            print("  python setup_goas.py")
            print("\nSkipping ocean basin assignment...")
            df_with_basins = df_filtered.copy()
            df_with_basins['ocean_basin'] = 'Unknown'
        else:
            print(f"✓ GOaS shapefile found: {goas_path}")
            
            # Load GOaS data
            print("Loading GOaS data...")
            goas_data = geographic.load_goas_data(goas_path)
            print(f"✓ Loaded {len(goas_data)} ocean regions")
            
            # Assign basins
            print("\nAssigning ocean basins...")
            df_with_basins = geographic.assign_ocean_basins(
                df_filtered,
                goas_data=goas_data,
                coord_col='coord'
            )
            
            # Show basin distribution
            basin_counts = df_with_basins['ocean_basin'].value_counts()
            print(f"\n✓ Ocean basin assignments:")
            for basin, count in basin_counts.items():
                print(f"  {basin}: {count} samples")
        
        # Save data with basins
        basins_tsv = OUTPUT_DIR / "03_with_ocean_basins.tsv"
        df_with_basins.to_csv(basins_tsv, sep='\t', index=False)
        print(f"\n✓ Saved data with basins: {basins_tsv}")
        
    except Exception as e:
        print(f"✗ Ocean basin assignment failed: {e}")
        import traceback
        traceback.print_exc()
        # Continue anyway - not critical for testing other components
        df_with_basins = df_filtered.copy()
        df_with_basins['ocean_basin'] = 'Unknown'
    
    # ========================================================================
    # STEP 5: Generate FASTA
    # ========================================================================
    print_section("STEP 5: Generate FASTA File")
    
    try:
        # Extract organism name
        organism = utils.extract_organism_name(TEST_TSV)
        print(f"Organism name: {organism}")
        
        # Create FASTA records
        fasta_records = []
        for _, row in df_filtered.iterrows():
            header = f"{organism}_{row['processid']}.COI-5P"
            sequence = row['nuc']
            
            # Validate sequence
            is_valid, reason = utils.validate_sequence(sequence, min_length=100)
            if is_valid:
                fasta_records.append((header, sequence))
            else:
                print(f"  Warning: Skipping {row['processid']}: {reason}")
        
        print(f"\n✓ Created {len(fasta_records)} FASTA records")
        
        # Show sequence statistics
        if len(fasta_records) > 0:
            seq_lengths = [len(seq) for _, seq in fasta_records]
            print(f"\nSequence length statistics:")
            print(f"  Min: {min(seq_lengths)} bp")
            print(f"  Max: {max(seq_lengths)} bp")
            print(f"  Mean: {sum(seq_lengths)/len(seq_lengths):.1f} bp")
        
        # Write FASTA
        fasta_path = OUTPUT_DIR / f"{organism}.fasta"
        utils.write_fasta(fasta_records, fasta_path)
        print(f"\n✓ Saved FASTA: {fasta_path}")
        
    except Exception as e:
        print(f"✗ FASTA generation failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # ========================================================================
    # STEP 6: Dereplicate Sequences
    # ========================================================================
    print_section("STEP 6: Dereplicate Sequences (Clustering)")
    
    try:
        print(f"Input: {len(fasta_records)} sequences")
        print(f"Clustering threshold: {cfg.dereplication.clustering_threshold} (99% identity)")
        
        # Check external tools
        if not utils.check_external_tool("mafft"):
            print("✗ MAFFT not found. Please install mafft.")
            return False
        if not utils.check_external_tool("trimal"):
            print("✗ trimAl not found. Please install trimal.")
            return False
        
        # Run dereplication
        print("\nRunning sequence alignment and clustering...")
        print("This may take a few minutes depending on dataset size...")
        
        consensus_records = dereplication.dereplicate_from_fasta(
            input_fasta=str(fasta_path),
            output_dir=str(OUTPUT_DIR / "dereplication"),
            threshold=cfg.dereplication.clustering_threshold,
            frequency_cutoff=cfg.dereplication.consensus_frequency_cutoff
        )
        
        print(f"\n✓ Identified {len(consensus_records)} unique genotypes")
        
        # Show genotype information
        print("\nGenotype clusters:")
        # Iterate over items (id, SeqRecord)
        for header, rec in sorted(consensus_records.items(), key=lambda x: x[0]):
            print(f" {header}: {len(rec.seq)} bp")
            
        # Save consensus sequences
        consensus_path = OUTPUT_DIR / f"{organism}_consensus.fasta"
        utils.write_fasta(
            [(rec.id, str(rec.seq)) for rec in consensus_records.values()],
            consensus_path
        )
        print(f"\n✓ Saved consensus sequences: {consensus_path}")
        
    except Exception as e:
        print(f"✗ Dereplication failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # ========================================================================
    # STEP 7: Assign Genotypes
    # ========================================================================
    print_section("STEP 7: Assign Samples to Genotypes")
    
    try:
        print(f"Assigning {len(df_filtered)} samples to {len(consensus_records)} genotypes...")
        print(f"Minimum identity threshold: {cfg.genotype_assignment.min_identity}")
    
        # Define output paths for this step
        annotated_tsv   = OUTPUT_DIR / f"{organism}_with_genotypes.tsv"
        diagnostics_csv = OUTPUT_DIR / f"{organism}_diagnostics.csv"
    
        # Call the library API (paths in, files out; returns stats dict)
        stats = genotype_assignment.assign_genotypes(
            metadata_path=str(filtered_tsv),         # <-- path to the filtered TSV you already saved in Step 3
            fasta_path=str(fasta_path),              # raw reads FASTA from Step 5
            consensus_path=str(consensus_path),      # consensus FASTA from Step 6
            output_path=str(annotated_tsv),          # where to write updated metadata with 'consensus_group'
            min_identity=cfg.genotype_assignment.min_identity,
            n_processes=cfg.n_threads,
            diagnostics_path=str(diagnostics_csv)    # set None to skip diagnostics file
        )
    
        # Load the outputs for downstream use in Step 8
        df_with_genotypes = pd.read_csv(annotated_tsv, sep="\t")
        diagnostics = pd.read_csv(diagnostics_csv) if diagnostics_csv.exists() else pd.DataFrame()
    
        # Show genotype distribution
        genotype_counts = df_with_genotypes['consensus_group'].value_counts(dropna=True)
        print(f"\n✓ Genotype assignments:")
        for genotype, count in genotype_counts.items():
            print(f"  {genotype}: {count} samples")
    
        # Unassigned summary (below threshold or no sequence)
        unassigned = df_with_genotypes['consensus_group'].isna().sum()
        if unassigned > 0:
            print(f"\n⚠ Warning: {unassigned} samples unassigned (below identity threshold or no sequence)")
    
        print("\nAssignment stats:", stats)
        print(f"✓ Saved annotated metadata: {annotated_tsv}")
        if diagnostics_csv.exists():
            print(f"✓ Saved diagnostics: {diagnostics_csv}")
    
    except Exception as e:
        print(f"✗ Genotype assignment failed: {e}")
        import traceback; traceback.print_exc()
        return False
    
        
    # ========================================================================
    # STEP 8: Generate Summary
    # ========================================================================
    print_section("STEP 8: Summary Statistics")
    
    try:
        # Combine all information
        df_final = df_with_genotypes.merge(
            df_with_basins[['processid', 'ocean_basin']],
            on='processid',
            how='left'
        )
        
        # Overall summary
        print(f"\n{'='*70}")
        print(f"Pipeline Summary for {organism}")
        print(f"{'='*70}")
        print(f"\nInput:")
        print(f"  Total samples in TSV: {len(df)}")
        print(f"  Samples with valid coordinates: {len(df_filtered)}")
        print(f"  Unique genotypes identified: {len(consensus_records)}")
        
        print(f"\nGenotype Distribution:")
        for genotype, count in genotype_counts.items():
            if pd.notna(genotype):
                percent = (count / len(df_filtered)) * 100
                print(f"  {genotype}: {count} ({percent:.1f}%)")
        
        if 'ocean_basin' in df_final.columns:
            print(f"\nOcean Basin Distribution:")
            for basin, count in df_final['ocean_basin'].value_counts().items():
                print(f"  {basin}: {count} samples")
        
        print(f"\nOutput Files:")
        output_files = [
            "01_parsed_metadata.tsv",
            "02_filtered_metadata.tsv",
            "03_with_ocean_basins.tsv",
            f"{organism}.fasta",
            f"{organism}_consensus.fasta",
            f"{organism}_with_genotypes.tsv",
            f"{organism}_diagnostics.csv",
        ]
        for filename in output_files:
            filepath = OUTPUT_DIR / filename
            if filepath.exists():
                print(f"  ✓ {filepath}")
        
        print(f"\n{'='*70}")
        print("✓ Pipeline test completed successfully!")
        print(f"{'='*70}\n")
        
        return True
        
    except Exception as e:
        print(f"✗ Summary generation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
