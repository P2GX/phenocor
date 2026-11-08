Correlation Analysis
====================

Pairwise associations between HPO terms can be analyzed using
``HPOCorrelationAnalyzer``.

For each pair, individuals with an unknown state for either term are
excluded. Self-pairs and ancestor-descendant pairs are excluded.
A remaining pair is tested only if it meets ``min_individuals`` and
both phenotypes vary within the pair-specific effective sample.

The analysis reports the Phi coefficient and a two-sided Fisher's exact
test p-value. Benjamini-Hochberg correction is applied across all valid
pairs tested in this dataset. Changing ``min_individuals`` can change
the set of tests and therefore the adjusted p-values.


Core usage
----------

.. code-block:: python

    from phenocor import HPOCorrelationAnalyzer

    analyzer = HPOCorrelationAnalyzer(dataset)

    result = analyzer.compute_correlation_matrix(
        min_individuals=20,
        n_jobs=-1,
    )

    result.correlation_results.head()

Result columns
--------------

``HPO_A``, ``HPO_B``
    Canonical HPO IDs of the two features.

``correlation``
    Phi coefficient.

``p_value``
    Two-sided Fisher's exact test p-value.

``adj_p_value``
    Benjamini-Hochberg-adjusted p-value.

``n(A:E/B:E)``
    Number of individuals with both features excluded.

``n(A:E/B:O)``
    Number of individuals with A excluded and B observed.

``n(A:O/B:E)``
    Number of individuals with A observed and B excluded.

``n(A:O/B:O)``
    Number of individuals with both features observed.

``n_individuals``
    Number of individuals with known states for both features after
    propagation. This equals the sum of the four contingency counts.

``HPO_A_label``, ``HPO_B_label``
    Feature labels, included when available.

``pmids``, ``n_pmids``
    Unique source PubMed IDs and their count, included when
    ``include_pmids=True``. IDs are aggregated across individuals
    contributing to the pair. They do not represent independent
    replications or studies confirming the association.

The Phi coefficient ranges from ``-1`` to ``1``:

- positive values indicate that two phenotypes are observed together more often than expected under independence,
- negative values indicate that two phenotypes are observed together less often than expected under independence,
- values near zero indicate little or no association.


Key parameters
--------------

``min_individuals``
    Minimum number of individuals with known states for both HPO terms after 
    propagation required for a pair to be tested. The default is ``20``.

``n_jobs``
    Number of parallel jobs used for pairwise calculations. Set to ``-1``
    to use all available CPU cores.


Significant correlations
------------------------

Significant correlations can be retrieved using an adjusted p-value
threshold:

.. code-block:: python

    significant = result.get_significant_correlations()

By default, correlations with adjusted p-values below ``0.05`` are returned.

A different significance threshold can be specified:

.. code-block:: python

    significant = result.get_significant_correlations(
        adj_pval_threshold=0.01,
    )

Correlations can also be filtered by their absolute Phi coefficient:

.. code-block:: python

    significant = result.get_significant_correlations(
        adj_pval_threshold=0.05,
        min_abs_correlation=0.3,
    )


Save results
------------

By default, ``save_correlation_results()`` saves correlations with adjusted
p-values below ``0.05``:

.. code-block:: python

    result.save_correlation_results(
        output_file="correlation_results.csv",
    )

Additional filtering by correlation strength can be applied:

.. code-block:: python

    result.save_correlation_results(
        output_file="correlation_results.csv",
        adj_pval_threshold=0.05,
        min_abs_correlation=0.3,
    )

To save all tested correlations without significance filtering:

.. code-block:: python

    result.save_correlation_results(
        output_file="correlation_results.csv",
        adj_pval_threshold=None,
    )

Results can be saved as ``.csv`` or ``.xlsx`` files.


Visualization
-------------

Correlation results can be visualized as a heatmap:

.. code-block:: python

    fig = result.plot_correlation_heatmap_with_significance()

The displayed correlations can be filtered by statistical significance and
correlation strength:

.. code-block:: python

    fig = result.plot_correlation_heatmap_with_significance(
        adj_pval_threshold=0.05,
        min_abs_correlation=0.1,
    )

``min_abs_correlation`` sets an optional minimum absolute Phi coefficient.
This filter does not determine statistical significance.

``adj_pval_threshold`` sets the adjusted p-value threshold.

A custom subtitle can also be provided:

.. code-block:: python

    fig = result.plot_correlation_heatmap_with_significance(
        title="Phenotype correlations",
    )

The heatmap can be saved directly:

.. code-block:: python

    result.save_correlation_heatmap(
        output_file="correlation_heatmap.html",
    )

The same filtering options can be applied when saving:

.. code-block:: python

    result.save_correlation_heatmap(
        output_file="correlation_heatmap.html",
        adj_pval_threshold=0.05,
        min_abs_correlation=0.1,
    )

The output format is determined automatically from the file extension.

- ``.html`` saves an interactive Plotly figure with zooming, panning,
  and hover tooltips.
- ``.svg``, ``.png``, ``.pdf``, and ``.jpg`` save static figures.

For static image formats, ``width``, ``height``, and ``scale`` can be used
to control the output size and resolution.

Interpretation and limitations
------------------------------

Associations describe phenotype annotation states among individuals
with known states for both features in the selected cohort.

Unreported features remain unknown unless their states can be inferred
through hierarchy propagation. Individuals with known states may differ
systematically from those with unknown states.

Associations may be influenced by reporting practices, publication
source, age, and cohort composition. They do not establish causal
relationships between phenotypes.

Static image export uses Kaleido and requires a compatible browser
installation. HTML export does not require a browser installation.