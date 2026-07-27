"""
Population Genetics Export Module for BOLDGenotyper

Exports BOLDGenotyper haplotype assignments to standard population-genetics
software formats. Populations are defined by a user-specified geographic
region column (``pop_column``), which by default is ``geo_category`` (e.g.,
'ocean_basin', 'PFAF_ID', 'REALM') — i.e., the same shapefile field used in
Phase 5 to assign samples to populations. Samples without a population
assignment are excluded from formats that require structured populations
(Arlequin, DnaSP); they are retained for haplotype-network visualization
(PopART) where region is encoded as a coloring trait rather than a
population definition.

Because BOLDGenotyper dereplicates ESVs at 100% identity, samples assigned
to a given haplotype have identical trimmed core sequences. The haplotype
consensus sequence is therefore equivalent to each member sample's trimmed
sequence, and is used as the per-sample sequence in all exports.

Supported formats:
- Arlequin (.arp): Population genetics and genomics analysis
- PopART/NEXUS: Haplotype network visualization and phylogeography
- DnaSP (.fas): DNA sequence polymorphism analysis
- Generic: CSV and FASTA formats for general use

Author: Steph Smith (symbioseas@outlook.com)
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


# Sentinel strings used by geographic.py for samples that could not be assigned.
# These are treated as equivalent to a missing population in all formats that
# require a structured population partition (Arlequin, DnaSP). Case-insensitive.
_UNASSIGNED_SENTINELS = {"unknown", "unassigned", "na", "n/a", "none", ""}


def _has_population(value: Any) -> bool:
    """Return True iff value is a real population label (not null/sentinel)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() not in _UNASSIGNED_SENTINELS


def export_arlequin(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    pop_column: str = 'ocean_basin',
    haplotype_column: str = 'haplotype_id'
) -> Path:
    """
    Export data in Arlequin format (.arp).

    Each population in the .arp file corresponds to one unique value of
    ``pop_column`` (e.g., one ocean basin, drainage basin, or ecoregion).
    Samples lacking a value in ``pop_column`` are excluded.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with haplotype assignments
    consensus_seqs : dict
        Dictionary mapping haplotype IDs to consensus sequences
    output_dir : Path
        Output directory for Arlequin files
    organism : str
        Organism name
    pop_column : str
        Column whose unique values define populations (default: 'ocean_basin').
        Typically the same column used in Phase 5 geographic assignment.
    haplotype_column : str
        Column holding each sample's haplotype assignment, used to look up
        the sequence in ``consensus_seqs`` (default: 'haplotype_id').

    Returns
    -------
    Path
        Path to created .arp file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    arp_file = output_dir / f"{organism}.arp"

    if pop_column not in df.columns:
        logger.warning(
            f"Arlequin export skipped: population column '{pop_column}' not found "
            f"in dataframe (available: geographic assignment may have been disabled)"
        )
        return arp_file

    # Keep samples with BOTH a haplotype assignment and a real population label
    # (excluding sentinels like "Unknown" used for samples that couldn't be assigned).
    pop_mask = df[pop_column].map(_has_population)
    df_eligible = df[df[haplotype_column].notna() & pop_mask].copy()

    n_total = int(df[haplotype_column].notna().sum())
    n_dropped = n_total - len(df_eligible)
    if n_dropped > 0:
        logger.info(
            f"  Arlequin: {n_dropped} haplotype-assigned samples lack an assigned '{pop_column}' "
            f"and were excluded from .arp export"
        )

    if len(df_eligible) == 0:
        logger.warning(
            f"No samples have both a haplotype assignment and a '{pop_column}' value; "
            f"Arlequin .arp will not be written."
        )
        return arp_file

    populations = sorted(df_eligible[pop_column].astype(str).unique())
    n_populations = len(populations)

    # Arlequin .arp header
    content = (
        f"[Profile]\n"
        f"    Title=\"{organism} populations exported by BOLDGenotyper "
        f"(populations = {pop_column})\"\n"
        f"    NbSamples={n_populations}\n"
        f"    DataType=DNA\n"
        f"    GenotypicData=0\n"
        f"    LocusSeparator=NONE\n"
        f"    MissingData='?'\n"
        f"\n"
        f"[Data]\n"
        f"    [[Samples]]\n"
        f"\n"
    )

    # One SampleData block per population
    skipped_no_seq = 0
    for pop in populations:
        pop_df = df_eligible[df_eligible[pop_column].astype(str) == pop]

        rows = []
        for _, row in pop_df.iterrows():
            sample_id = str(row['processid']).replace(' ', '_')
            hap = row[haplotype_column]
            seq = consensus_seqs.get(hap)
            if seq is None:
                skipped_no_seq += 1
                continue
            # Arlequin sequence-data syntax: sample_name  count  sequence
            rows.append(f"        {sample_id} 1 {seq}")

        if not rows:
            continue

        pop_safe = str(pop).replace(' ', '_').replace('"', '')
        content += (
            f"    SampleName=\"{pop_safe}\"\n"
            f"    SampleSize={len(rows)}\n"
            f"    SampleData= {{\n"
        )
        content += "\n".join(rows) + "\n"
        content += "    }\n\n"

    # Structure block: one group containing all populations (enables AMOVA out of the box)
    content += "[[Structure]]\n"
    content += f"    StructureName=\"All populations\"\n"
    content += f"    NbGroups=1\n"
    content += f"    Group={{\n"
    for pop in populations:
        pop_safe = str(pop).replace(' ', '_').replace('"', '')
        content += f"        \"{pop_safe}\"\n"
    content += "    }\n"

    with open(arp_file, 'w') as f:
        f.write(content)

    logger.info(
        f"  ✓ Arlequin .arp exported: {arp_file} "
        f"({n_populations} populations, {len(df_eligible)} samples"
        + (f", {skipped_no_seq} samples skipped for missing consensus)" if skipped_no_seq else ")")
    )

    # Companion populations.csv with per-population summary (samples, haplotypes, dominant species)
    pop_file = output_dir / "populations.csv"
    summary_rows = []
    for pop in populations:
        pop_df = df_eligible[df_eligible[pop_column].astype(str) == pop]
        n_samples = len(pop_df)
        n_haps = pop_df[haplotype_column].nunique()
        if 'species' in pop_df.columns:
            mode_sp = pop_df['species'].mode()
            dom_species = mode_sp.iloc[0] if len(mode_sp) > 0 else 'Unknown'
        else:
            dom_species = 'Unknown'
        summary_rows.append({
            'population': pop,
            'n_samples': n_samples,
            'n_haplotypes': n_haps,
            'dominant_species': dom_species,
        })
    pd.DataFrame(summary_rows).to_csv(pop_file, index=False)

    return arp_file


def export_popart_nexus(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    pop_column: str = 'ocean_basin',
    haplotype_column: str = 'haplotype_id'
) -> Tuple[Path, Path]:
    """
    Export data in PopART/NEXUS format for haplotype-network construction.

    All samples with a haplotype assignment are included as taxa. Geographic
    region is encoded as a trait (used by PopART to color network nodes)
    rather than as a hard population partition; samples lacking a
    ``pop_column`` value are retained with region='Unassigned' so the network
    still shows them.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with haplotype assignments
    consensus_seqs : dict
        Dictionary mapping haplotype IDs to consensus sequences
    output_dir : Path
        Output directory
    organism : str
        Organism name
    pop_column : str
        Column whose values are used as the region trait (default: 'ocean_basin')
    haplotype_column : str
        Column holding each sample's haplotype assignment (default: 'haplotype_id')

    Returns
    -------
    tuple
        Paths to .nexus file and .traits file (.traits is currently empty;
        traits are written inline in the .nexus)
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    nexus_file = output_dir / f"{organism}.nexus"
    traits_file = output_dir / f"{organism}.traits"

    df_assigned = df[df[haplotype_column].notna()].copy()

    if len(df_assigned) == 0:
        logger.warning("No samples with haplotype assignments to export")
        return nexus_file, traits_file

    has_pop = pop_column in df_assigned.columns

    taxa: List[str] = []
    sequences: List[str] = []
    traits: List[Tuple[str, str, str, str]] = []
    sample_counter: Dict[str, int] = {}
    skipped_no_seq = 0

    for _, row in df_assigned.iterrows():
        hap = row[haplotype_column]
        if hap not in consensus_seqs:
            skipped_no_seq += 1
            continue

        sample_counter[hap] = sample_counter.get(hap, 0) + 1
        taxon_name = f"{hap}_{sample_counter[hap]:03d}"
        taxa.append(taxon_name)
        sequences.append(consensus_seqs[hap])

        species_val = row.get('species')
        species = (str(species_val).replace(' ', '_')
                   if pd.notna(species_val) else 'Unknown')

        if has_pop and _has_population(row.get(pop_column)):
            region = str(row[pop_column]).replace(' ', '_')
        else:
            region = 'Unassigned'

        traits.append((taxon_name, str(hap), species, region))

    if not taxa:
        logger.warning("No taxa produced for PopART NEXUS export.")
        return nexus_file, traits_file

    seq_length = max(len(s) for s in sequences)

    nexus_content = "#NEXUS\n\n"
    nexus_content += "BEGIN TAXA;\n"
    nexus_content += f"    DIMENSIONS NTAX={len(taxa)};\n"
    nexus_content += "    TAXLABELS\n"
    nexus_content += "        " + " ".join(taxa) + "\n"
    nexus_content += "    ;\nEND;\n\n"

    nexus_content += "BEGIN CHARACTERS;\n"
    nexus_content += f"    DIMENSIONS NCHAR={seq_length};\n"
    nexus_content += "    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
    nexus_content += "    MATRIX\n"
    for taxon, seq in zip(taxa, sequences):
        seq_padded = seq.ljust(seq_length, '?')
        nexus_content += f"        {taxon}  {seq_padded}\n"
    nexus_content += "    ;\nEND;\n\n"

    nexus_content += "BEGIN TRAITS;\n"
    nexus_content += "    DIMENSIONS NTRAITS=3;\n"
    nexus_content += "    FORMAT LABELS=YES MISSING=? SEPARATOR=COMMA;\n"
    nexus_content += f"    TRAITLABELS haplotype species {pop_column};\n"
    nexus_content += "    MATRIX\n"
    for taxon, hap, species, region in traits:
        nexus_content += f"        {taxon},{hap},{species},{region}\n"
    nexus_content += "    ;\nEND;\n"

    with open(nexus_file, 'w') as f:
        f.write(nexus_content)

    msg = f"  ✓ PopART NEXUS exported: {nexus_file} ({len(taxa)} taxa)"
    if skipped_no_seq:
        msg += f"; {skipped_no_seq} samples skipped for missing consensus"
    logger.info(msg)

    # Companion CSV (one row per sample/taxon)
    pop_csv = output_dir / "popart_traits.csv"
    pd.DataFrame(traits, columns=['taxon', 'haplotype', 'species', pop_column]).to_csv(
        pop_csv, index=False
    )

    return nexus_file, traits_file


def export_dnasp(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    pop_column: str = 'ocean_basin',
    haplotype_column: str = 'haplotype_id'
) -> Path:
    """
    Export data in DnaSP-compatible FASTA format with population labels.

    FASTA record IDs use the convention ``{population}|{processid}`` so DnaSP's
    "Define Sequence Sets" → "Auto-detect by name" reads populations directly
    from the header. Samples lacking a value in ``pop_column`` are excluded.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with haplotype assignments
    consensus_seqs : dict
        Dictionary mapping haplotype IDs to consensus sequences
    output_dir : Path
        Output directory
    organism : str
        Organism name
    pop_column : str
        Column whose unique values define populations (default: 'ocean_basin')
    haplotype_column : str
        Column holding each sample's haplotype assignment (default: 'haplotype_id')

    Returns
    -------
    Path
        Path to created .fas file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fas_file = output_dir / f"{organism}.fas"

    if pop_column not in df.columns:
        logger.warning(
            f"DnaSP export skipped: population column '{pop_column}' not found in dataframe."
        )
        return fas_file

    pop_mask = df[pop_column].map(_has_population)
    df_eligible = df[df[haplotype_column].notna() & pop_mask].copy()
    n_total = int(df[haplotype_column].notna().sum())
    n_dropped = n_total - len(df_eligible)
    if n_dropped > 0:
        logger.info(
            f"  DnaSP: {n_dropped} haplotype-assigned samples lack an assigned '{pop_column}' "
            f"and were excluded from .fas export"
        )

    if len(df_eligible) == 0:
        logger.warning(
            f"No samples have both a haplotype assignment and a '{pop_column}' value; "
            f"DnaSP .fas will not be written."
        )
        return fas_file

    records = []
    skipped_no_seq = 0

    for _, row in df_eligible.iterrows():
        hap = row[haplotype_column]
        seq = consensus_seqs.get(hap)
        if seq is None:
            skipped_no_seq += 1
            continue

        sample_id = str(row['processid']).replace(' ', '_')
        pop = str(row[pop_column]).replace(' ', '_')
        species_val = row.get('species')
        species = (str(species_val).replace(' ', '_')
                   if pd.notna(species_val) else 'Unknown')

        record = SeqRecord(
            Seq(seq),
            id=f"{pop}|{sample_id}",
            description=(
                f"[haplotype={hap};species={species};{pop_column}={pop};processid={sample_id}]"
            ),
        )
        records.append(record)

    SeqIO.write(records, fas_file, 'fasta')

    msg = (
        f"  ✓ DnaSP .fas exported: {fas_file} "
        f"({len(records)} samples across {df_eligible[pop_column].nunique()} populations)"
    )
    if skipped_no_seq:
        msg += f"; {skipped_no_seq} skipped for missing consensus"
    logger.info(msg)

    return fas_file


def export_generic(
    df: pd.DataFrame,
    consensus_seqs: Dict[str, str],
    output_dir: Path,
    organism: str,
    pop_column: str = 'ocean_basin',
    haplotype_column: str = 'haplotype_id'
) -> Dict[str, Path]:
    """
    Export generic formats: membership CSV, consensus FASTA, haplotype summary.

    Includes all samples with a haplotype assignment regardless of whether they
    have a population (``pop_column``) value. Useful for users who want to
    pull haplotype assignments into custom scripts.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with haplotype assignments
    consensus_seqs : dict
        Dictionary mapping haplotype IDs to consensus sequences
    output_dir : Path
        Output directory
    organism : str
        Organism name
    pop_column : str
        Column whose values are reported as the population label
    haplotype_column : str
        Column holding each sample's haplotype assignment

    Returns
    -------
    dict
        Paths to created files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    files: Dict[str, Path] = {}

    df_assigned = df[df[haplotype_column].notna()].copy()
    if len(df_assigned) == 0:
        logger.warning("No samples with haplotype assignments to export")
        return files

    # 1. Sample → haplotype membership table
    membership_file = output_dir / "haplotype_membership.csv"
    membership_cols = ['processid', haplotype_column]
    if 'species' in df_assigned.columns:
        membership_cols.append('species')
    if 'lat' in df_assigned.columns and 'lon' in df_assigned.columns:
        membership_cols.extend(['lat', 'lon'])
    if pop_column in df_assigned.columns:
        membership_cols.append(pop_column)
    if 'country/ocean' in df_assigned.columns:
        membership_cols.append('country/ocean')

    membership_df = df_assigned[membership_cols].copy()
    membership_df['consensus_sequence'] = membership_df[haplotype_column].map(consensus_seqs)
    membership_df['sequence_length'] = membership_df['consensus_sequence'].apply(
        lambda x: len(x) if isinstance(x, str) else 0
    )
    membership_df.to_csv(membership_file, index=False)
    files['membership'] = membership_file

    # 2. Consensus FASTA (one record per haplotype)
    alignment_file = output_dir / "haplotype_consensus.fasta"
    records = []
    for hap, seq in consensus_seqs.items():
        sub = df_assigned[df_assigned[haplotype_column] == hap]
        if len(sub) == 0:
            continue
        species_mode = sub['species'].mode() if 'species' in sub.columns else pd.Series(dtype=object)
        species = species_mode.iloc[0] if len(species_mode) else 'Unknown'
        records.append(SeqRecord(
            Seq(seq), id=str(hap),
            description=f"n={len(sub)} species={species}"
        ))
    SeqIO.write(records, alignment_file, 'fasta')
    files['alignment'] = alignment_file

    # 3. Haplotype summary table (one row per haplotype)
    haplotypes_file = output_dir / "haplotype_summary.csv"
    rows = []
    for hap, seq in consensus_seqs.items():
        sub = df_assigned[df_assigned[haplotype_column] == hap]
        if len(sub) == 0:
            continue
        species_mode = sub['species'].mode() if 'species' in sub.columns else pd.Series(dtype=object)
        info = {
            'haplotype_id': hap,
            'sequence': seq,
            'sequence_length': len(seq),
            'n_samples': len(sub),
            'species': species_mode.iloc[0] if len(species_mode) else 'Unknown',
        }
        if pop_column in sub.columns:
            counts = sub[pop_column].dropna().value_counts()
            info['primary_population'] = counts.index[0] if len(counts) else 'Unassigned'
            info['n_populations'] = len(counts)
        rows.append(info)
    pd.DataFrame(rows).to_csv(haplotypes_file, index=False)
    files['haplotypes'] = haplotypes_file

    logger.info(f"  ✓ Generic exports written: {len(files)} files")
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
- **geographic region**: Geographic classification (ocean basin, drainage basin, ecoregion, etc.)
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
- GitHub: https://github.com/SymbioSeas/BOLDGenotyper
- Email: symbioseas@outlook.com

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
    pop_column: str = 'ocean_basin',
    haplotype_column: str = 'haplotype_id'
) -> Dict[str, Any]:
    """
    Master function to export to multiple population genetics formats.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated dataframe with haplotype assignments
    consensus_fasta_path : Path
        Path to consensus sequences FASTA file (Phase 2 output)
    output_dir : Path
        Base output directory (will create ``exports/`` and per-format subdirs)
    organism : str
        Organism name
    formats : list, optional
        List of formats to export. Options: 'arlequin', 'popart', 'dnasp',
        'generic', 'all'. Default is ['all'].
    pop_column : str
        Column whose unique values define populations in Arlequin/DnaSP, and
        the region trait for PopART. Default 'ocean_basin'; typically the
        same column used in Phase 5 geographic assignment.
    haplotype_column : str
        Column holding each sample's haplotype assignment (default: 'haplotype_id')

    Returns
    -------
    dict
        Summary of exported files and formats
    """
    if formats is None or 'all' in formats:
        formats = ['arlequin', 'popart', 'dnasp', 'generic']

    logger.info(
        f"Exporting population genetics formats (populations defined by '{pop_column}')..."
    )

    exports_dir = output_dir / "exports"
    exports_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Any] = {'formats': formats, 'pop_column': pop_column, 'files': {}}

    # Load consensus sequences
    consensus_seqs: Dict[str, str] = {}
    cfp = Path(consensus_fasta_path)
    if cfp.exists():
        for record in SeqIO.parse(str(cfp), 'fasta'):
            consensus_seqs[record.id] = str(record.seq)
    else:
        logger.warning(f"Consensus FASTA not found: {cfp}; popgen export may be empty.")

    if 'arlequin' in formats:
        results['files']['arlequin'] = export_arlequin(
            df, consensus_seqs, exports_dir / "arlequin", organism,
            pop_column=pop_column, haplotype_column=haplotype_column,
        )

    if 'popart' in formats:
        nexus_file, traits_file = export_popart_nexus(
            df, consensus_seqs, exports_dir / "popart", organism,
            pop_column=pop_column, haplotype_column=haplotype_column,
        )
        results['files']['popart_nexus'] = nexus_file
        results['files']['popart_traits'] = traits_file

    if 'dnasp' in formats:
        results['files']['dnasp'] = export_dnasp(
            df, consensus_seqs, exports_dir / "dnasp", organism,
            pop_column=pop_column, haplotype_column=haplotype_column,
        )

    if 'generic' in formats:
        results['files']['generic'] = export_generic(
            df, consensus_seqs, exports_dir / "generic", organism,
            pop_column=pop_column, haplotype_column=haplotype_column,
        )

    readme_path = create_export_readme(exports_dir, organism, formats, results['files'])
    results['files']['readme'] = readme_path

    logger.info(f"✓ Population genetics exports complete: {exports_dir}")
    return results
