"""
Population Genetics Export Module for BOLDGenotyper

This module exports genotyping results to various population genetics software formats,
enabling seamless integration with downstream analysis tools.

Supported formats:
- Arlequin (.arp): Population genetics and genomics analysis
- PopART/NEXUS: Haplotype network visualization and phylogeography
- DnaSP (.fas): DNA sequence polymorphism analysis
- Structure (.str): Genetic structure analysis (if SNPs available)
- Generic: CSV and FASTA formats for general use

Author: Steph Smith (steph.smith@unc.edu)
"""

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
from datetime import datetime
import logging

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def export_arlequin(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    group_by: str = 'consensus_group'
) -> Path:
    """
    Export data in Arlequin format (.arp).

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    consensus_seqs : dict
        Dictionary mapping consensus group IDs to sequences
    output_dir : Path
        Output directory for Arlequin files
    organism : str
        Organism name
    group_by : str
        Column to use for population grouping (default: 'consensus_group')

    Returns
    -------
    Path
        Path to created .arp file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    arp_file = output_dir / f"{organism}.arp"

    # Filter to samples with assignments
    df_assigned = df[df[group_by].notna()].copy()

    if len(df_assigned) == 0:
        logger.warning("No samples with genotype assignments to export")
        return arp_file

    # Get unique populations/genotypes
    populations = sorted(df_assigned[group_by].unique())
    n_populations = len(populations)

    # Start building Arlequin file
    content = f"""[Profile]
    Title="{organism} genotypes from BOLD - exported by BOLDGenotyper"
    NbSamples={n_populations}
    DataType=DNA
    GenotypicData=0
    LocusSeparator=NONE
    MissingData='?'

[Data]
    [[Samples]]

"""

    # Add each population
    for pop in populations:
        pop_df = df_assigned[df_assigned[group_by] == pop]
        sample_size = len(pop_df)

        content += f"""    SampleName="{pop}"
    SampleSize={sample_size}
    SampleData= {{
"""

        # Add sequences for this population
        for idx, row in pop_df.iterrows():
            sample_id = row['processid']
            # Use the consensus sequence if available, otherwise try to get individual sequence
            if pop in consensus_seqs:
                seq = consensus_seqs[pop]
            elif 'nuc' in row and pd.notna(row['nuc']):
                seq = str(row['nuc'])
            else:
                continue

            # Arlequin format: frequency sample_name haplotype_id sequence
            content += f"        1 {sample_id} 1 {seq}\n"

        content += "    }\n    \n"

    # Write file
    with open(arp_file, 'w') as f:
        f.write(content)

    logger.info(f"  ✓ Arlequin format exported: {arp_file}")

    # Also create a populations definition file
    pop_file = output_dir / "populations.txt"
    pop_mapping = df_assigned.groupby(group_by).agg({
        'processid': 'count',
        'species': lambda x: x.mode()[0] if len(x.mode()) > 0 else 'Unknown'
    }).reset_index()
    pop_mapping.columns = ['population', 'n_samples', 'species']

    if 'ocean_basin' in df_assigned.columns:
        basin_mapping = df_assigned.groupby(group_by)['ocean_basin'].agg(
            lambda x: x.mode()[0] if len(x.mode()) > 0 else 'Unknown'
        )
        pop_mapping = pop_mapping.merge(
            basin_mapping.rename('ocean_basin'),
            left_on='population',
            right_index=True,
            how='left'
        )

    pop_mapping.to_csv(pop_file, index=False)

    return arp_file


def export_popart_nexus(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    group_by: str = 'consensus_group'
) -> Tuple[Path, Path]:
    """
    Export data in PopART/NEXUS format.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    consensus_seqs : dict
        Dictionary mapping consensus group IDs to sequences
    output_dir : Path
        Output directory
    organism : str
        Organism name
    group_by : str
        Column to use for population grouping

    Returns
    -------
    tuple
        Paths to .nexus file and .traits file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    nexus_file = output_dir / f"{organism}.nexus"
    traits_file = output_dir / f"{organism}.traits"

    # Filter to samples with assignments
    df_assigned = df[df[group_by].notna()].copy()

    if len(df_assigned) == 0:
        logger.warning("No samples with genotype assignments to export")
        return nexus_file, traits_file

    # Build taxa list with unique sample IDs
    taxa = []
    sequences = []
    traits = []

    sample_counter = {}  # Track sample number within each genotype

    for idx, row in df_assigned.iterrows():
        genotype = row[group_by]
        sample_id = row['processid']

        # Create unique taxon name: genotype_samplenum
        if genotype not in sample_counter:
            sample_counter[genotype] = 0
        sample_counter[genotype] += 1

        taxon_name = f"{genotype}_{sample_counter[genotype]:03d}"
        taxa.append(taxon_name)

        # Get sequence - use consensus for this genotype
        if genotype in consensus_seqs:
            seq = consensus_seqs[genotype]
        elif 'nuc' in row and pd.notna(row['nuc']):
            seq = str(row['nuc'])
        else:
            seq = '?' * 650  # placeholder

        sequences.append(seq)

        # Build traits
        species = row.get('species', 'Unknown') if pd.notna(row.get('species')) else 'Unknown'
        species = species.replace(' ', '_')

        basin = row.get('ocean_basin', 'Unknown') if pd.notna(row.get('ocean_basin')) else 'Unknown'
        basin = basin.replace(' ', '_')

        traits.append((taxon_name, genotype, species, basin))

    # Get sequence length
    seq_length = max(len(s) for s in sequences) if sequences else 650

    # Write NEXUS file
    nexus_content = f"""#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX={len(taxa)};
    TAXLABELS
        {' '.join(taxa)}
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR={seq_length};
    FORMAT DATATYPE=DNA MISSING=? GAP=-;
    MATRIX
"""

    for taxon, seq in zip(taxa, sequences):
        # Pad sequences to same length
        seq_padded = seq.ljust(seq_length, '?')
        nexus_content += f"        {taxon}  {seq_padded}\n"

    nexus_content += """    ;
END;

BEGIN TRAITS;
    DIMENSIONS NTRAITS=3;
    FORMAT LABELS=YES MISSING=? SEPARATOR=COMMA;
    TRAITLABELS genotype species ocean_basin;
    MATRIX
"""

    for taxon, genotype, species, basin in traits:
        nexus_content += f"        {taxon},{genotype},{species},{basin}\n"

    nexus_content += """    ;
END;
"""

    with open(nexus_file, 'w') as f:
        f.write(nexus_content)

    logger.info(f"  ✓ PopART NEXUS format exported: {nexus_file}")

    # Also create a populations CSV file
    pop_csv = output_dir / "populations.csv"
    pop_df = pd.DataFrame(traits, columns=['sample', 'genotype', 'species', 'ocean_basin'])
    pop_df.to_csv(pop_csv, index=False)

    return nexus_file, traits_file


def export_dnasp(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    group_by: str = 'consensus_group'
) -> Path:
    """
    Export data in DnaSP format (FASTA with pop labels).

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    consensus_seqs : dict
        Dictionary mapping consensus group IDs to sequences
    output_dir : Path
        Output directory
    organism : str
        Organism name
    group_by : str
        Column to use for population grouping

    Returns
    -------
    Path
        Path to created .fas file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fas_file = output_dir / f"{organism}.fas"

    # Filter to samples with assignments
    df_assigned = df[df[group_by].notna()].copy()

    if len(df_assigned) == 0:
        logger.warning("No samples with genotype assignments to export")
        return fas_file

    records = []
    sample_counter = {}

    for idx, row in df_assigned.iterrows():
        genotype = row[group_by]
        sample_id = row['processid']

        # Create unique sample name
        if genotype not in sample_counter:
            sample_counter[genotype] = 0
        sample_counter[genotype] += 1
        unique_id = f"{genotype}_{sample_counter[genotype]:03d}"

        # Get sequence
        if genotype in consensus_seqs:
            seq = consensus_seqs[genotype]
        elif 'nuc' in row and pd.notna(row['nuc']):
            seq = str(row['nuc'])
        else:
            continue

        # Build metadata string
        species = row.get('species', 'Unknown') if pd.notna(row.get('species')) else 'Unknown'
        species = species.replace(' ', '_')

        basin = row.get('ocean_basin', 'Unknown') if pd.notna(row.get('ocean_basin')) else 'Unknown'
        basin = basin.replace(' ', '_')

        # DnaSP format: >sample_id [key=value;key=value;...]
        description = f"[genotype={genotype};species={species};basin={basin};processid={sample_id}]"

        record = SeqRecord(
            Seq(seq),
            id=unique_id,
            description=description
        )
        records.append(record)

    # Write FASTA file
    SeqIO.write(records, fas_file, 'fasta')

    logger.info(f"  ✓ DnaSP format exported: {fas_file}")

    return fas_file


def export_generic(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    group_by: str = 'consensus_group'
) -> Dict[str, Path]:
    """
    Export generic formats: CSV membership table, FASTA alignment, haplotypes table.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    consensus_seqs : dict
        Dictionary mapping consensus group IDs to sequences
    output_dir : Path
        Output directory
    organism : str
        Organism name
    group_by : str
        Column to use for population grouping

    Returns
    -------
    dict
        Paths to created files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {}

    # Filter to samples with assignments
    df_assigned = df[df[group_by].notna()].copy()

    if len(df_assigned) == 0:
        logger.warning("No samples with genotype assignments to export")
        return files

    # 1. Genotype membership table
    membership_file = output_dir / "genotype_membership.csv"

    membership_cols = ['processid', group_by]
    if 'species' in df_assigned.columns:
        membership_cols.append('species')
    if 'lat' in df_assigned.columns and 'lon' in df_assigned.columns:
        membership_cols.extend(['lat', 'lon'])
    if 'ocean_basin' in df_assigned.columns:
        membership_cols.append('ocean_basin')
    if 'country' in df_assigned.columns:
        membership_cols.append('country')

    membership_df = df_assigned[membership_cols].copy()

    # Add consensus sequence if available
    membership_df['consensus_sequence'] = membership_df[group_by].map(consensus_seqs)
    membership_df['sequence_length'] = membership_df['consensus_sequence'].apply(
        lambda x: len(x) if pd.notna(x) else 0
    )

    membership_df.to_csv(membership_file, index=False)
    files['membership'] = membership_file

    # 2. Alignment FASTA (consensus sequences)
    alignment_file = output_dir / "alignment.fasta"
    records = []

    for genotype, seq in consensus_seqs.items():
        # Get metadata for this genotype
        genotype_df = df_assigned[df_assigned[group_by] == genotype]
        if len(genotype_df) == 0:
            continue

        n_samples = len(genotype_df)
        species = genotype_df['species'].mode()[0] if 'species' in genotype_df.columns and len(genotype_df['species'].mode()) > 0 else 'Unknown'

        record = SeqRecord(
            Seq(seq),
            id=genotype,
            description=f"n={n_samples} species={species}"
        )
        records.append(record)

    SeqIO.write(records, alignment_file, 'fasta')
    files['alignment'] = alignment_file

    # 3. Haplotypes table
    haplotypes_file = output_dir / "haplotypes.csv"

    haplotypes = []
    for genotype, seq in consensus_seqs.items():
        genotype_df = df_assigned[df_assigned[group_by] == genotype]
        if len(genotype_df) == 0:
            continue

        haplotype_info = {
            'haplotype_id': genotype,
            'sequence': seq,
            'sequence_length': len(seq),
            'n_samples': len(genotype_df),
            'species': genotype_df['species'].mode()[0] if 'species' in genotype_df.columns and len(genotype_df['species'].mode()) > 0 else 'Unknown'
        }

        if 'ocean_basin' in genotype_df.columns:
            basins = genotype_df['ocean_basin'].value_counts()
            haplotype_info['primary_basin'] = basins.index[0]
            haplotype_info['n_basins'] = len(basins)

        haplotypes.append(haplotype_info)

    haplotypes_df = pd.DataFrame(haplotypes)
    haplotypes_df.to_csv(haplotypes_file, index=False)
    files['haplotypes'] = haplotypes_file

    logger.info(f"  ✓ Generic formats exported: {len(files)} files")

    return files


def create_export_readme(
    output_dir: Path,
    organism: str,
    exported_formats: List[str],
    files: Dict[str, Any]
) -> Path:
    """
    Create README.md explaining exported formats.

    Parameters
    ----------
    output_dir : Path
        Base export directory
    organism : str
        Organism name
    exported_formats : list
        List of format names that were exported
    files : dict
        Dictionary of exported file paths

    Returns
    -------
    Path
        Path to created README
    """
    readme_path = output_dir / "README.md"

    content = f"""# Population Genetics Exports for {organism}

## Overview

This directory contains genotyping results exported to various population genetics
software formats for downstream analysis. All exports were generated by
BOLDGenotyper on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.

## Exported Formats

"""

    if 'arlequin' in exported_formats:
        content += """### Arlequin Format

**Software:** Arlequin (http://cmpg.unibe.ch/software/arlequin35/)
**Purpose:** Population genetics and genomics analysis
**Files:**
- `arlequin/{organism}.arp` - Main Arlequin project file
- `arlequin/populations.txt` - Population definitions and metadata

**Usage:**
1. Open Arlequin
2. File → Open → Select the .arp file
3. Configure analysis parameters
4. Execute → Run computations

"""

    if 'popart' in exported_formats:
        content += """### PopART/NEXUS Format

**Software:** PopART (http://popart.otago.ac.nz/)
**Purpose:** Haplotype network visualization and phylogeography
**Files:**
- `popart/{organism}.nexus` - NEXUS format alignment with trait data
- `popart/populations.csv` - Population mapping table

**Usage:**
1. Open PopART
2. File → Open → Select the .nexus file
3. Choose network method (TCS, Median-Joining, etc.)
4. Adjust parameters and compute network

"""

    if 'dnasp' in exported_formats:
        content += """### DnaSP Format

**Software:** DnaSP (http://www.ub.edu/dnasp/)
**Purpose:** DNA sequence polymorphism analysis
**Files:**
- `dnasp/{organism}.fas` - FASTA format with population labels in headers

**Usage:**
1. Open DnaSP
2. File → Open Data File → Select the .fas file
3. Configure population assignments from header labels
4. Run analyses (diversity, neutrality tests, etc.)

"""

    if 'generic' in exported_formats:
        content += """### Generic Formats

**Purpose:** General use and custom analyses
**Files:**
- `generic/alignment.fasta` - Consensus sequences for all genotypes
- `generic/genotype_membership.csv` - Sample-to-genotype mapping with metadata
- `generic/haplotypes.csv` - Summary table of unique haplotypes

**Usage:**
These files can be imported into any software that accepts standard formats:
- FASTA alignments: MEGA, Geneious, PhyloSuite, etc.
- CSV tables: R, Python, Excel, etc.

"""

    content += """## File Format Details

### Genotype Naming

Genotypes are named using the pattern: `c{cluster_id}_n{sample_count}`

Example: `c15_n386` = cluster 15, containing 386 samples

### Metadata Fields

- **genotype**: Consensus group identifier
- **species**: Taxonomic species assignment (if available)
- **ocean_basin**: Geographic region (if available)
- **processid**: BOLD database sample identifier

### Sequences

- **Type:** COI barcode sequences (mitochondrial DNA)
- **Length:** Variable, typically 400-700 bp
- **Missing data:** Represented as `?` or `N` depending on format
- **Gaps:** Represented as `-`

## Citation

If you use these data in a publication, please cite:

1. **BOLDGenotyper:**
   Smith et al. (in prep) BOLDGenotyper: Automated genotyping and quality
   control for DNA barcode data.

2. **BOLD Database:**
   Ratnasingham S, Hebert PDN (2007) BOLD: The Barcode of Life Data System.
   Molecular Ecology Notes 7: 355-364.

3. **Original Data Source:**
   [Your publication or data citation here]

## Support

For questions about these exports or BOLDGenotyper:
- GitHub: https://github.com/yourusername/boldgenotyper
- Email: support@example.com

For format-specific software support, consult the respective software documentation.
"""

    with open(readme_path, 'w') as f:
        f.write(content)

    logger.info(f"  ✓ Export README created: {readme_path}")

    return readme_path


def export_population_genetics_formats(
    df: pd.DataFrame,
    consensus_fasta_path: Path,
    output_dir: Path,
    organism: str,
    formats: Optional[List[str]] = None,
    group_by: str = 'consensus_group'
) -> Dict[str, Any]:
    """
    Master function to export to multiple population genetics formats.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with genotype assignments
    consensus_fasta_path : Path
        Path to consensus sequences FASTA file
    output_dir : Path
        Base output directory (will create subdirectories for each format)
    organism : str
        Organism name
    formats : list, optional
        List of formats to export. Options: 'arlequin', 'popart', 'dnasp', 'generic', 'all'
        Default is ['all']
    group_by : str
        Column to use for population grouping (default: 'consensus_group')

    Returns
    -------
    dict
        Summary of exported files and formats
    """
    if formats is None or 'all' in formats:
        formats = ['arlequin', 'popart', 'dnasp', 'generic']

    logger.info("Exporting population genetics formats...")

    # Create exports directory
    exports_dir = output_dir / "exports"
    exports_dir.mkdir(parents=True, exist_ok=True)

    results = {'formats': formats, 'files': {}}

    # Load consensus sequences
    consensus_seqs = {}
    if consensus_fasta_path.exists():
        for record in SeqIO.parse(consensus_fasta_path, 'fasta'):
            consensus_seqs[record.id] = str(record.seq)
    else:
        logger.warning(f"Consensus FASTA not found: {consensus_fasta_path}")
        logger.warning("Will attempt to use sequences from dataframe")

    # Export each requested format
    if 'arlequin' in formats:
        arlequin_dir = exports_dir / "arlequin"
        arp_file = export_arlequin(df, consensus_seqs, arlequin_dir, organism, group_by)
        results['files']['arlequin'] = arp_file

    if 'popart' in formats:
        popart_dir = exports_dir / "popart"
        nexus_file, traits_file = export_popart_nexus(df, consensus_seqs, popart_dir, organism, group_by)
        results['files']['popart_nexus'] = nexus_file
        results['files']['popart_traits'] = traits_file

    if 'dnasp' in formats:
        dnasp_dir = exports_dir / "dnasp"
        fas_file = export_dnasp(df, consensus_seqs, dnasp_dir, organism, group_by)
        results['files']['dnasp'] = fas_file

    if 'generic' in formats:
        generic_dir = exports_dir / "generic"
        generic_files = export_generic(df, consensus_seqs, generic_dir, organism, group_by)
        results['files']['generic'] = generic_files

    # Create README
    readme_path = create_export_readme(exports_dir, organism, formats, results['files'])
    results['files']['readme'] = readme_path

    logger.info(f"✓ Population genetics exports complete: {exports_dir}")

    return results
