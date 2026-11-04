phenocor
============

**phenocor** is a Python library for analyzing phenotype correlations in
`GA4GH Phenopacket <https://www.ga4gh.org/product/phenopackets/>`_
cohorts using
`Human Phenotype Ontology <https://hpo.jax.org/>`_ (HPO) annotations.

It provides tools for building analysis-ready phenotype datasets,
accounting for the HPO hierarchy, and identifying pairwise associations
between phenotypic features.



Overview
--------

Phenotypic features can show patterns of association across individuals. 
Characterizing these relationships can help identify how phenotypes are 
related within Phenopacket cohorts.

**phenocor enables:**

- Construction of structured phenotype datasets from Phenopackets
- Hierarchy-aware processing of HPO annotations
- Pairwise correlation analysis of HPO terms
- Interactive visualization of correlation heatmaps



Quick Example
-------------

.. code-block:: python

   from pathlib import Path
   from google.protobuf.json_format import Parse
   from phenopackets import Phenopacket
   from phenocor import (
      PhenotypeDatasetBuilder,
      HPOCorrelationAnalyzer,
   )

   # Load phenopackets
   phenopacket_dir = Path("path/to/your/phenopackets/")

   phenopackets = []
   for file_path in sorted(phenopacket_dir.glob("*.json")):
      with open(file_path, "r", encoding="utf-8") as f:
         data: str = f.read()
         phenopacket: Phenopacket = Parse(data, Phenopacket())
         phenopackets.append(phenopacket)

   # Build dataset
   dataset = PhenotypeDatasetBuilder(phenopackets).build()

   # Run correlation analysis
   analyzer = HPOCorrelationAnalyzer(dataset)
   results = analyzer.compute_correlation_matrix()
   results.correlation_results.head()

This example shows the basic workflow:

**load phenopackets → build dataset → compute correlations**



.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorial

.. toctree::
   :maxdepth: 1
   :caption: Usage

   usage/dataset
   usage/correlation

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api/phenocor