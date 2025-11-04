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
import os
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
TEST_TSV = TEST_DATA_DIR / "Euprymna.tsv"
OUTPUT_DIR = Path("tests/Euprymna_output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Adjustable similarity threshold (use 0.80 for highly diverse taxa like Carcharhinus)
SIMILARITY_THRESHOLD = float(os.environ.get("BG_SIMILARITY_THRESHOLD", "0.80"))

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

        # Update configuration including adjustable similarity threshold
        cfg = cfg.update(
            genotype_assignment__min_identity=SIMILARITY_THRESHOLD,
            output_dir=OUTPUT_DIR,
            n_threads=2,
            keep_intermediates=True,
            log_level="DEBUG"
        )
        print("✓ Configuration loaded")
        print(f"  Clustering threshold: {cfg.dereplication.clustering_threshold}")
        print(f"  Min identity (similarity threshold): {cfg.genotype_assignment.min_identity}")
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
                coord_col="coord"
            )
            
            nn = int(df_with_basins[['lat','lon']].notna().all(axis=1).sum()) if {'lat','lon'}.issubset(df_with_basins.columns) else 0
            print(f"✓ Parsed coordinate rows carried into basins: {nn}")
            unk = int((df_with_basins.get('ocean_basin','Unknown') == 'Unknown').sum())
            print(f" Unknown basin count: {unk} / {len(df_with_basins)}")
            
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
        print(f"Clustering threshold: {cfg.dereplication.clustering_threshold:.3f} (~{cfg.dereplication.clustering_threshold*100:.0f}% identity)")
        
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
    # STEP 7.5: Assign taxonomy to consensus groups (majority rule)
    # ========================================================================
    print_section("STEP 7.5: Assign Taxonomy to Consensus Groups")
    
    try:
        # df_with_genotypes should already include at least: processid, consensus_group, species (and ideally genus)
        majority_thr = cfg.taxonomy.majority_species_threshold
        assign_table, species_counts = utils.assign_consensus_taxonomy(
            df_with_genotypes,
            group_col="consensus_group",
            species_col="species",
            genus_col="genus",               # if not present, helper derives from species
            majority_threshold=majority_thr
        )
        
        required = {"consensus_group", "assigned_sp", "assignment_level", "assignment_notes"}
        missing = required - set(assign_table.columns)
        if missing:
            raise RuntimeError(f"assign_consensus_taxonomy returned missing columns: {missing}\n"f"Got columns: {list(assign_table.columns)}")
    
        # Build 'consensus_group_sp' like "Sphyrna zygaena c7_n46" or "Sphyrna c9_n12"
        def _strip_prefix(s: str) -> str:
            return s.replace("consensus_", "") if isinstance(s, str) else s
    
        assign_table["consensus_group_short"] = assign_table["consensus_group"].map(_strip_prefix)
    
        def _join_label(row):
            sp = row["assigned_sp"]
            short = row["consensus_group_short"]
            if not sp:
                return short  # fall back to bare group if totally unassigned
            return f"{sp} {short}"
    
        assign_table["consensus_group_sp"] = assign_table.apply(_join_label, axis=1)
    
        # Save per-group species tallies and final assignment
        species_counts_out = OUTPUT_DIR / f"{organism}_species_by_consensus.csv"
        assign_table_out   = OUTPUT_DIR / f"{organism}_consensus_taxonomy.csv"
        species_counts.to_csv(species_counts_out, index=False)
        assign_table.to_csv(assign_table_out, index=False)
    
        # Attach 'assigned_sp' and 'consensus_group_sp' to every row by group
        df_with_genotypes = df_with_genotypes.merge(
            assign_table[["consensus_group", "assigned_sp", "consensus_group_sp", "assignment_level", "assignment_notes", "majority_fraction"]],
            on="consensus_group",
            how="left",
            validate="many_to_one"
        )
        
        # Load and merge cluster-seq assignments
        cluster_seq_path = OUTPUT_DIR / f"{organism}_consensus_taxonomy_seq.csv"
        if cluster_seq_path.exists():
            cluster_seq_df = pd.read_csv(cluster_seq_path)
            # expected columns: consensus_group, cluster_seq_sp, cluster_seq_level,
            #					cluster_seq_best_identity, cluster_seq_qcov
            df_with_genotypes = df_with_genotypes.merge(
                cluster_seq_df, on="consensus_group", how="left", validate="many_to_one"
            )
        
        def _final_label(row):
            final_sp, final_level, prov = utils.pick_final_group_taxon(
                cluster_sp=row.get("cluster_seq_sp", ""),
                cluster_level=row.get("cluster_seq_level", ""),
                cluster_id=row.get("cluster_seq_best_identity", 0.0),
                cluster_qcov=row.get("cluster_seq_qcov", 0.0),
                majority_sp=row.get("assigned_sp", ""),
                majority_level=row.get("assignment_level", ""),
                majority_frac=row.get("majority_fraction", 0.0),
                cfg_taxonomy=cfg.taxonomy,
            )
            short = row["consensus_group"].replace("consensus_", "") if isinstance(row.get("consensus_group"), str) else ""
            # Build the display label
            label = f"{final_sp} {short}".strip() if final_sp else short
            return pd.Series({
                "final_group_sp": final_sp,
                "final_group_level": final_level,
                "tax_provenance": prov,
                "consensus_group_sp": label
            })
            
        final_cols = df_with_genotypes.apply(_final_label, axis=1)
        for c in final_cols.columns:
            df_with_genotypes[c] = final_cols[c]
            
        df_with_genotypes = df_with_genotypes.loc[:, ~df_with_genotypes.columns.duplicated(keep='last')]
        
    
        # Flag taxonomy conflicts at the processid level
        # conflict if assigned at species-level and original species != assigned_sp;
        # if assigned at genus-level, conflict if original genus != assigned_sp
        if "genus" not in df_with_genotypes.columns:
            df_with_genotypes["genus"] = df_with_genotypes["species"].map(lambda s: s.split(" ")[0] if isinstance(s, str) and s.strip() else "")
    
        def _conflict(row):
            lvl = row.get("assignment_level", "")
            if lvl == "species":
                return (isinstance(row.get("species",""), str) 
                        and row["species"].strip() != row.get("assigned_sp",""))
            elif lvl == "genus":
                return (isinstance(row.get("genus",""), str)
                        and row["genus"].strip() != row.get("assigned_sp",""))
            return False
    
        df_with_genotypes["taxonomy_conflict"] = df_with_genotypes.apply(_conflict, axis=1)
    
        # Update on-disk TSV so later steps read the enriched columns
        annotated_tsv = OUTPUT_DIR / f"{organism}_with_genotypes.tsv"
        df_with_genotypes.to_csv(annotated_tsv, sep="\t", index=False)
        print(f"✓ Saved updated genotypes with taxonomy: {annotated_tsv}")
        print(f"✓ Species tallies by group: {species_counts_out}")
        print(f"✓ Chosen assignments per group: {assign_table_out}")
    
    except Exception as e:
        print(f"✗ Taxonomy assignment failed: {e}")
        import traceback; traceback.print_exc()
        return False

    # ========================================================================
    # STEP 8: Generate Summary
    # ========================================================================
    print_section("STEP 8: Summary Statistics")
    
    try:
        # Build one final annotated table with geography + genotype
        merged_out        = OUTPUT_DIR / f"{organism}_annotated.csv"
        geno_path         = OUTPUT_DIR / f"{organism}_with_genotypes.tsv"   # already written in Step 7
        geo_path          = OUTPUT_DIR / "03_with_ocean_basins.tsv"         # already written in Step 4
        crosstab_out      = OUTPUT_DIR / f"{organism}_genotype_by_basin.csv"
        crosstab_sp_out   = OUTPUT_DIR / f"{organism}_genotypeSP_by_basin.csv"
        basin_summary_out = OUTPUT_DIR / f"{organism}_basin_summary.csv"
    
        # Re-load (robust if someone reruns Step 8 standalone)
        geno_df = df_with_genotypes if 'df_with_genotypes' in locals() else pd.read_csv(geno_path, sep="\t")
        geo_df  = df_with_basins    if 'df_with_basins'    in locals() else pd.read_csv(geo_path,  sep="\t")
    
        # Only append the geography columns we need
        geo_keep = [c for c in ['processid', 'lat', 'lon', 'ocean_basin'] if c in geo_df.columns]
        
        # Merge first
        df_final = geno_df.merge(geo_df[geo_keep], on='processid', how='left', validate='one_to_one')
        # Drop duplicate-named columns
        df_final = df_final.loc[:, ~df_final.columns.duplicated(keep='last')]
        
        # Bring new taxonomy columns to the front
        front_more = [c for c in [
            "assigned_sp", "assignment_level", "assignment_notes", "taxonomy_conflict", "consensus_group_sp"
        ] if c in df_final.columns]
        
        # For specific column ordering, do it here
        cols_front = [c for c in ["processid"] + front_more if c in df_final.columns]
        cols_geo   = [c for c in ["lat","lon","ocean_basin"] if c in df_final.columns]
        other_cols = [c for c in df_final.columns if c not in cols_front + cols_geo]
        df_final   = df_final[cols_front + other_cols + cols_geo]
        
        # Save the final unified annotated CSV
        df_final.to_csv(merged_out, index=False)
        
        # --------- Summary stats ----------
        print(f"\n{'='*70}")
        print(f"Pipeline Summary for {organism}")
        print(f"{'='*70}")
        print(f"\nInput:")
        print(f"  Total samples in TSV: {len(df)}")
        print(f"  Samples with valid coordinates after filtering: {len(df_filtered)}")
        print(f"  Unique genotypes identified: {len(consensus_records)}")
    
        # Genotype distribution (from Step 7)
        print(f"\nGenotype Distribution:")
        for genotype, count in genotype_counts.items():
            if pd.notna(genotype):
                percent = (count / len(df_filtered)) * 100
                print(f"  {genotype}: {count} ({percent:.1f}%)")
    
        # Coordinate & basin coverage
        has_lat = 'lat' in df_final.columns
        has_lon = 'lon' in df_final.columns
        has_basin = 'ocean_basin' in df_final.columns
    
        if has_lat and has_lon:
            n_coords = int((df_final['lat'].notna() & df_final['lon'].notna()).sum())
            print(f"\nCoordinate Coverage:")
            print(f"  Parsed lat/lon for: {n_coords} / {len(df_final)} "
                  f"({(n_coords/len(df_final))*100:.1f}%)")
            if n_coords > 0:
                lat_min, lat_max = df_final['lat'].min(), df_final['lat'].max()
                lon_min, lon_max = df_final['lon'].min(), df_final['lon'].max()
                print(f"  Latitude range:  {lat_min:.2f} to {lat_max:.2f}")
                print(f"  Longitude range: {lon_min:.2f} to {lon_max:.2f}")
    
        if has_basin:
            print(f"\nOcean Basin Distribution:")
            basin_counts = df_final['ocean_basin'].fillna('Unknown').value_counts()
            for basin, count in basin_counts.items():
                print(f"  {basin}: {count} samples")
    
        # Genotype × Basin crosstab (printed and saved)
        if has_basin:
            ct = pd.crosstab(
                df_final['consensus_group'].fillna('Unassigned'),
                df_final['ocean_basin'].fillna('Unknown')
            )
            # Sort for readable display: rows by total desc, columns by total desc
            ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index,
                        ct.sum(axis=0).sort_values(ascending=False).index]
    
            print(f"\nGenotype × Ocean Basin (top 10 rows):")
            # print only top 10 rows to keep logs readable
            print(ct.head(10).to_string())
    
            ct.to_csv(crosstab_out)
            
        # GenotypeSP x Basin crosstab (species-augmented labels)
        if has_basin and 'consensus_group_sp' in df_final.columns:
            ct_sp = pd.crosstab(
                df_final['consensus_group_sp'].fillna('Unassigned'),
                df_final['ocean_basin'].fillna('Unknown')
            )
            ct_sp = ct_sp.loc[ct_sp.sum(axis=1).sort_values(ascending=False).index, ct_sp.sum(axis=0).sort_values(ascending=False).index]
            print(f"\nGenotypeSP x Ocean Basin (top 10 rows):")
            print(ct_sp.head(10).to_string())
            ct_sp.to_csv(crosstab_sp_out)
    
        # Basin summary: counts, assigned rate, mean identity (if available)
        if has_basin:
            cols = ['processid', 'ocean_basin', 'consensus_group']
            if 'best_identity' in df_final.columns:
                cols.append('best_identity')
            bs = df_final[cols].copy()
            bs['assigned'] = bs['consensus_group'].notna()
            grp = bs.groupby(bs['ocean_basin'].fillna('Unknown'), dropna=False)
    
            summary = grp.agg(
                total=('processid', 'count'),
                assigned=('assigned', 'sum')
            )
            summary['assigned_rate'] = (summary['assigned'] / summary['total']).round(3)
            if 'best_identity' in bs.columns:
                # mean identity computed on assigned rows only
                summary['mean_best_identity'] = grp.apply(
                    lambda g: g.loc[g['assigned'], 'best_identity'].mean()
                ).round(4)
    
            summary = summary.sort_values(['total', 'assigned'], ascending=False)
            summary.to_csv(basin_summary_out)
            print(f"\nPer-basin assignment summary (saved): {basin_summary_out}")
    
        # --------- Output file manifest ----------
        print(f"\nOutput Files:")
        output_files = [
            "01_parsed_metadata.tsv",
            "02_filtered_metadata.tsv",
            "03_with_ocean_basins.tsv",
            f"{organism}.fasta",
            f"{organism}_consensus.fasta",
            f"{organism}_with_genotypes.tsv",
            f"{organism}_diagnostics.csv",
            f"{organism}_species_by_consensus.csv",
            f"{organism}_consensus_taxonomy.csv",
            f"{organism}_annotated.csv",          # NEW unified file
            f"{organism}_genotype_by_basin.csv",  # NEW crosstab
            f"{organism}_genotypeSP_by_basin.csv",
            f"{organism}_basin_summary.csv",      # NEW per-basin summary
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
