Build Dataset
=============

The ``PhenotypeDatasetBuilder`` converts Phenopackets conforming to the
`GA4GH Phenopacket Schema <https://phenopacket-schema.readthedocs.io/en/latest/>`_
into an analysis-ready ``PhenotypeDataset`` for phenotype correlation analysis.

The builder extracts HPO annotations as observed, excluded, or unreported states. 
Observed states propagate to ancestor terms, and excluded states propagate to 
descendant terms, restricted to terms represented in the dataset. Unreported states 
remain unknown unless an observed or excluded state can be inferred through HPO 
propagation.

This example builds a single dataset from all Phenopackets in the directory. 
For disease-specific analyses, select individuals with the relevant disease 
before building the dataset.

Each input Phenopacket should represent a different individual.
Phenopacket IDs define matrix row identities and must be unique.
Different Phenopacket IDs are not automatically deduplicated using
subject IDs.

If the same term is explicitly recorded as both observed and excluded,
observed status takes precedence and a warning is emitted. The same
rule applies when multiple input IDs resolve to the same canonical
HPO term.

During hierarchy propagation, explicit annotations are preserved when
they conflict with inferred annotations. If both observed and excluded
states are inferred for an unknown value, observed status takes
precedence and a warning is emitted.


Core usage
----------

Phenopackets can be loaded from JSON files using the official
``phenopackets`` protobuf classes and ``google.protobuf`` utilities.


.. code-block:: python

    from pathlib import Path
    from google.protobuf.json_format import Parse
    from phenopackets import Phenopacket
    from phenocor import PhenotypeDatasetBuilder

    # Locate your local phenopacket directory
    phenopacket_dir = Path("path/to/your/phenopackets/")

    # Iterate and parse JSON files into formal Phenopacket objects
    phenopackets = []
    for file_path in sorted(phenopacket_dir.glob("*.json")):
        with open(file_path, "r", encoding="utf-8") as f:
            data = f.read()
            phenopacket = Parse(data, Phenopacket())
            phenopackets.append(phenopacket)

    # Initialize the builder and construct the analysis-ready dataset
    builder = PhenotypeDatasetBuilder(phenopackets)
    dataset = builder.build()


Dataset overview
^^^^^^^^^^^^^^^^

The resulting ``PhenotypeDataset`` contains four main components:

- **phenotype_matrix** — individual-level HPO data.
- **hpo_labels** — HPO term labels.
- **relationship_mask** — excludes identical terms and ancestor–descendant term pairs from correlation analysis.
- **phenopackets** — the original phenopacket objects.

Individual and HPO identifiers are available through:

.. code-block:: python

    dataset.individual_ids
    dataset.hpo_ids

Associated PubMed and disease identifiers can be retrieved with:

.. code-block:: python

    dataset.get_pmids()
    dataset.get_diseases()


HPO configuration
-----------------

By default, the latest available HPO release is used.
For reproducible analyses, explicitly select an HPO release or supply
a versioned local ontology file. Record the ontology release together
with the phenopacket data version.

A specific HPO release can be selected for reproducible analyses:

.. code-block:: python

    builder = PhenotypeDatasetBuilder(
        phenopackets,
        hpo_release="v2023-10-09",
    )

Alternatively, a local HPO ontology file can be supplied:

.. code-block:: python

    builder = PhenotypeDatasetBuilder(
        phenopackets,
        hpo_file="path/to/hp.json",
    )


Next steps
----------

After constructing the dataset, see :doc:`correlation` to perform phenotype–phenotype correlation analysis.
