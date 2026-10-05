BOLDGenotyper Documentation
============================

**COI haplotype discovery and biogeographic analysis from BOLD data**

BOLDGenotyper identifies COI haplotypes in Barcode of Life Data System (BOLD) records
for any taxonomic group: quality control, haplotype discovery as exact sequence variants
of a core alignment region, one status for every record, taxonomic labels with
label-conflict evidence, optional phylogenetics, region assignment from
geographic-quality coordinates, and figures. Every run folder includes a ``README.md``
that defines each output file and column.

.. image:: https://img.shields.io/badge/python-3.9%2B-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python Version

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://img.shields.io/badge/version-2.0.0-green.svg
   :target: https://github.com/SymbioSeas/BOLDGenotyper/releases
   :alt: Version

Quick Links
-----------

* `GitHub Repository <https://github.com/SymbioSeas/BOLDGenotyper>`_
* `README (installation, usage, outputs) <https://github.com/SymbioSeas/BOLDGenotyper/blob/main/README.md>`_
* `API guide <https://github.com/SymbioSeas/BOLDGenotyper/blob/main/docs/guides/API_REFERENCE.md>`_
* `Migration from 1.0.0 <https://github.com/SymbioSeas/BOLDGenotyper/blob/main/docs/guides/MIGRATION_v1_to_v2.md>`_
* `Changelog <https://github.com/SymbioSeas/BOLDGenotyper/blob/main/CHANGELOG.md>`_
* `Issue Tracker <https://github.com/SymbioSeas/BOLDGenotyper/issues>`_

Quick Start
-----------

Installation::

    git clone https://github.com/SymbioSeas/BOLDGenotyper.git
    cd BOLDGenotyper
    conda env create -f environment.yml
    conda activate boldgenotyper
    pip install -e .

Basic usage (choose ``--singleton-distance`` with ``boldgenotyper-sweep``)::

    boldgenotyper data/Sphyrna_lewini.tsv --singleton-distance 0.015 --build-tree

API Reference
-------------

.. toctree::
   :maxdepth: 1

   api/cli
   api/config
   api/metadata
   api/quality_control
   api/dereplication
   api/haplotype_ids
   api/record_status
   api/taxon_labels
   api/taxon_evidence
   api/taxon_curation
   api/layout
   api/run_readme
   api/provenance
   api/geographic
   api/geographic_enhancement
   api/shapefile_registry
   api/goas_downloader
   api/phylogenetics
   api/msa_visualization
   api/divergence_analysis
   api/species_analysis
   api/metadata_summary
   api/visualization
   api/figure_pages
   api/tree_figures
   api/run_figures
   api/metadata_figures
   api/plot_export
   api/plot_regeneration
   api/haplotype_query
   api/query_assignment
   api/query_outputs
   api/query_overlays
   api/query_report
   api/parameter_sweep
   api/popgen_export
   api/reports
   api/report
   api/utils

Citation
--------

See ``CITATION.cff`` in the repository.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
