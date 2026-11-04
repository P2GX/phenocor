import numpy as np
import pandas as pd
import pytest

from phenocor.analysis.hpo_correlation_analyzer import HPOCorrelationAnalyzer
from phenocor.core import PhenotypeDataset


@pytest.fixture
def large_dataset():
    n_samples = 60

    matrix = pd.DataFrame(
        {
            "HP:0001": [1] * 30 + [0] * 30,
            "HP:0002": [0] * 10 + [1] * 30 + [0] * 20,
            "HP:0003": [0] * 20 + [1] * 30 + [0] * 10,
            "HP:0004": [0] * 15 + [1] * 30 + [0] * 15,
        },
        index=[f"Patient_{i}" for i in range(n_samples)],
        dtype=float,
    )

    # Synthetic terms are treated as unrelated; self-pairs are excluded.
    mask = pd.DataFrame(
        0.0,
        index=matrix.columns,
        columns=matrix.columns,
    )
    np.fill_diagonal(mask.values, np.nan)

    return PhenotypeDataset(
        phenotype_matrix=matrix,
        hpo_labels={c: c for c in matrix.columns},
        relationship_mask=mask,
        phenopackets=[],
    )


@pytest.fixture
def strongly_correlated_dataset():
    n_samples = 80

    matrix = pd.DataFrame(
        {
            "HP:0001": [1] * 40 + [0] * 40,
            "HP:0002": [1] * 40 + [0] * 40,
            "HP:0003": [0] * 40 + [1] * 40,
            "HP:0004": [1, 0] * 40,
        },
        index=[f"Patient_{i}" for i in range(n_samples)],
        dtype=float,
    )

    mask = pd.DataFrame(
        0.0,
        index=matrix.columns,
        columns=matrix.columns,
    )
    np.fill_diagonal(mask.values, np.nan)

    return PhenotypeDataset(
        phenotype_matrix=matrix,
        hpo_labels={c: c for c in matrix.columns},
        relationship_mask=mask,
        phenopackets=[],
    )


def test_initialization(large_dataset):
    analyzer = HPOCorrelationAnalyzer(large_dataset)

    assert analyzer.hpo_matrix.shape == (60, 4)
    assert analyzer.n_hpo_terms == 4


def test_compute_correlation_matrix(large_dataset):
    analyzer = HPOCorrelationAnalyzer(large_dataset)

    results = analyzer.compute_correlation_matrix(
        min_individuals=10,
        n_jobs=1,
        include_pmids=False,
    )

    table = results.correlation_results

    assert isinstance(table, pd.DataFrame)
    assert len(table) == 6

    expected_columns = {
        "HPO_A",
        "HPO_B",
        "correlation",
        "p_value",
        "adj_p_value",
        "n(A:E/B:E)",
        "n(A:E/B:O)",
        "n(A:O/B:E)",
        "n(A:O/B:O)",
        "n_individuals",
    }
    assert expected_columns.issubset(table.columns)
    assert (table["n_individuals"] == 60).all()

    count_columns = [
        "n(A:E/B:E)",
        "n(A:E/B:O)",
        "n(A:O/B:E)",
        "n(A:O/B:O)",
    ]
    assert (
        table[count_columns].sum(axis=1) == table["n_individuals"]
    ).all()


def test_compute_correlation_matrix_happy_path(
    strongly_correlated_dataset,
):
    analyzer = HPOCorrelationAnalyzer(strongly_correlated_dataset)

    results = analyzer.compute_correlation_matrix(
        min_individuals=10,
        n_jobs=1,
        include_pmids=False,
    )

    table = results.correlation_results

    assert len(table) == 6

    expected_correlations = {
        frozenset(("HP:0001", "HP:0002")): 1.0,
        frozenset(("HP:0001", "HP:0003")): -1.0,
        frozenset(("HP:0002", "HP:0003")): -1.0,
        frozenset(("HP:0001", "HP:0004")): 0.0,
        frozenset(("HP:0002", "HP:0004")): 0.0,
        frozenset(("HP:0003", "HP:0004")): 0.0,
    }

    for _, row in table.iterrows():
        pair = frozenset((row["HPO_A"], row["HPO_B"]))
        assert row["correlation"] == pytest.approx(
            expected_correlations[pair]
        )


def test_filter_weak_correlations(strongly_correlated_dataset):
    analyzer = HPOCorrelationAnalyzer(strongly_correlated_dataset)

    results = analyzer.compute_correlation_matrix(
        min_individuals=10,
        n_jobs=1,
        include_pmids=False,
    )

    filtered = results.get_significant_correlations(
        adj_pval_threshold=0.05,
        min_abs_correlation=0.9,
    )

    expected_pairs = {
        frozenset(("HP:0001", "HP:0002")),
        frozenset(("HP:0001", "HP:0003")),
        frozenset(("HP:0002", "HP:0003")),
    }
    actual_pairs = {
        frozenset((row["HPO_A"], row["HPO_B"]))
        for _, row in filtered.iterrows()
    }

    assert actual_pairs == expected_pairs
    assert len(filtered) == 3
    assert np.allclose(filtered["correlation"].abs(), 1.0)
    assert (filtered["adj_p_value"] < 0.05).all()