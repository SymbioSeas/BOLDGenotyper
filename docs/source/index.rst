BOLDGenotyper Documentation
============================

**Automated COI sequence genotyping and biogeographic analysis from BOLD database data**

BOLDGenotyper is a Python package that enables researchers to identify and analyze COI (Cytochrome Oxidase I) genotypes from the BOLD (Barcode of Life Database) for any taxonomic group. The pipeline performs sequence dereplication, genotype assignment, geographic filtering, ocean basin classification, and visualization of biogeographic patterns.

.. image:: https://img.shields.io/badge/python-3.8%2B-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python Version

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://img.shields.io/badge/version-0.1.0-green.svg
   :target: https://github.com/SymbioSeas/BOLDGenotyper/releases
   :alt: Version

Quick Links
-----------

* `GitHub Repository <https://github.com/SymbioSeas/BOLDGenotyper>`_
* `Issue Tracker <https://github.com/SymbioSeas/BOLDGenotyper/issues>`_
* `Changelog <https://github.com/SymbioSeas/BOLDGenotyper/blob/main/CHANGELOG.md>`_

Quick Start
-----------

Installation::

    git clone https://github.com/SymbioSeas/BOLDGenotyper.git
    cd BOLDGenotyper
    conda env create -f boldgenotyper_env.yml
    conda activate depredation

Basic Usage::

    boldgenotyper data/Sphyrna_lewini.tsv

For detailed installation and usage instructions, see the :doc:`installation` guide and :doc:`quickstart`.

Key Features
------------

* **Automated genotyping**: Identifies unique COI genotypes through sequence clustering
* **Geographic analysis**: Filters coordinates and assigns samples to ocean basins
* **Visualization**: Creates publication-ready maps and figures
* **Phylogenetic analysis**: Optional tree building and visualization
* **Modular design**: Flexible pipeline with customizable parameters

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   usage
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/cli
   api/config
   api/dereplication
   api/genotype_assignment
   api/geographic
   api/metadata
   api/phylogenetics
   api/reports
   api/utils
   api/visualization

.. toctree::
   :maxdepth: 2
   :caption: Development

   contributing
   changelog

Citation
--------

If you use BOLDGenotyper in your research, please cite:

.. code-block:: bibtex

   @software{boldgenotyper2025,
     author = {Smith, Steph},
     title = {BOLDGenotyper: Automated COI Sequence Genotyping},
     year = {2025},
     version = {0.1.0},
     url = {https://github.com/SymbioSeas/BOLDGenotyper},
     doi = {10.5281/zenodo.XXXXXXX}
   }

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
