Sequence Dereplication
======================

The ``dereplication`` module provides functions for sequence clustering and consensus generation.

.. automodule:: boldgenotyper.dereplication
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The sequence dereplication process identifies unique genotypes by:

1. **Aligning sequences** with MAFFT
2. **Trimming alignments** with trimAl
3. **Calculating pairwise distances** between sequences
4. **Hierarchical clustering** at a specified threshold (default: 99% identity)
5. **Generating consensus sequences** for each cluster

Key Functions
-------------

.. autofunction:: boldgenotyper.dereplication.calculate_pairwise_distances

.. autofunction:: boldgenotyper.dereplication.cluster_sequences

.. autofunction:: boldgenotyper.dereplication.generate_consensus

.. autofunction:: boldgenotyper.dereplication.run_dereplication

Configuration
-------------

.. autoclass:: boldgenotyper.config.DereplicationConfig
   :members:
   :undoc-members:

Examples
--------

Basic Usage
^^^^^^^^^^^

.. code-block:: python

    from pathlib import Path
    from boldgenotyper import dereplication
    from boldgenotyper.config import DereplicationConfig

    # Configure dereplication
    config = DereplicationConfig(
        clustering_threshold=0.01,  # 99% identity
        consensus_frequency_cutoff=0.7,
        min_sequence_length=400
    )

    # Run dereplication
    consensus_path = dereplication.run_dereplication(
        fasta_path=Path("sequences.fasta"),
        output_dir=Path("output/"),
        config=config
    )

Advanced Usage
^^^^^^^^^^^^^^

.. code-block:: python

    from boldgenotyper import dereplication, utils
    from Bio import SeqIO

    # Load sequences
    sequences = list(SeqIO.parse("sequences.fasta", "fasta"))

    # Calculate pairwise distances
    distances = dereplication.calculate_pairwise_distances(
        sequences,
        ignore_gaps=True,
        ignore_ambiguous=True
    )

    # Cluster sequences
    clusters = dereplication.cluster_sequences(
        distances,
        threshold=0.01,
        method='average'
    )

    # Generate consensus for each cluster
    consensus_sequences = []
    for cluster_id, member_indices in clusters.items():
        members = [sequences[i] for i in member_indices]
        consensus = dereplication.generate_consensus(
            members,
            frequency_cutoff=0.7
        )
        consensus_sequences.append(consensus)

Scientific Background
---------------------

The 99% identity threshold is based on standard COI barcoding practices:

- **Hebert et al. (2003)**: Established COI as DNA barcode marker
- **Meier et al. (2006)**: 98-99% identity for intraspecific variation
- **BOLD Standards**: 99% used for Barcode Index Numbers (BINs)

References
----------

* Hebert, P.D.N., et al. (2003). Biological identifications through DNA barcodes.
  *Proceedings of the Royal Society B*, 270(1512), 313-321.

* Ratnasingham, S. & Hebert, P.D.N. (2013). A DNA-based registry for all animal
  species: The Barcode Index Number (BIN) system. *PLoS ONE*, 8(7), e66213.

See Also
--------

* :mod:`boldgenotyper.genotype_assignment` - Assign samples to consensus genotypes
* :mod:`boldgenotyper.config` - Configuration management
* :mod:`boldgenotyper.utils` - Utility functions
