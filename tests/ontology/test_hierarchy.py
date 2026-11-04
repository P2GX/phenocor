import pathlib

import numpy as np
import pandas as pd
import pytest

from phenocor.ontology import HPOHierarchyEngine


TEST_DIR = pathlib.Path(__file__).parent.parent.resolve()
HP_JSON_FILE = TEST_DIR / "data" / "hp.json"


@pytest.fixture
def hpo_engine():
    return HPOHierarchyEngine(hpo_file=str(HP_JSON_FILE))


@pytest.fixture
def toy_terms():
    return ["HP:0001250", "HP:0020219", "HP:0012759"]


def test_propagate_hierarchy_basic(hpo_engine):
    """
    Test observed upward propagation and excluded downward propagation.
    """
    df = pd.DataFrame(
        {
            "HP:0020219": [1, 1, np.nan],
            "HP:0001250": [np.nan, 1, 0],
            "HP:0012759": [1, np.nan, 1],
        },
        index=["Patient_1", "Patient_2", "Patient_3"],
    )

    propagated = hpo_engine.propagate(df)

    assert isinstance(propagated, pd.DataFrame)
    assert list(propagated.index) == ["Patient_1", "Patient_2", "Patient_3"]

    # Patient_1: child observed -> ancestor observed
    assert propagated.loc["Patient_1", "HP:0020219"] == 1
    assert propagated.loc["Patient_1", "HP:0001250"] == 1
    assert propagated.loc["Patient_1", "HP:0012759"] == 1

    # Patient_2: explicitly observed terms remain observed
    assert propagated.loc["Patient_2", "HP:0020219"] == 1
    assert propagated.loc["Patient_2", "HP:0001250"] == 1

    # Patient_3: excluded ancestor -> excluded descendant
    assert propagated.loc["Patient_3", "HP:0001250"] == 0
    assert propagated.loc["Patient_3", "HP:0020219"] == 0
    assert propagated.loc["Patient_3", "HP:0012759"] == 1


def test_propagate_removes_invalid_terms(hpo_engine):
    """
    Invalid HPO terms should be filtered out.
    """
    df = pd.DataFrame(
        {
            "HP:0020219": [1, np.nan],
            "HP:9999999": [1, 0],
        },
        index=["Patient_1", "Patient_2"],
    )

    propagated = hpo_engine.propagate(df)

    assert "HP:0020219" in propagated.columns
    assert "HP:9999999" not in propagated.columns


def test_propagate_preserves_index(hpo_engine):
    """
    Row index should be preserved after propagation.
    """
    df = pd.DataFrame(
        {
            "HP:0020219": [1, 0],
            "HP:0001250": [np.nan, np.nan],
        },
        index=["A", "B"],
    )

    propagated = hpo_engine.propagate(df)

    assert list(propagated.index) == ["A", "B"]


def test_propagate_duplicate_columns_are_merged(hpo_engine):
    """
    Duplicate canonical columns should be merged with max().
    """
    values = np.array(
        [
            [1, np.nan],
            [0, 1],
        ]
    )
    df = pd.DataFrame(
        values,
        index=["Patient_1", "Patient_2"],
        columns=["HP:0020219", "HP:0020219"],
    )

    propagated = hpo_engine.propagate(df)

    assert list(propagated.columns).count("HP:0020219") == 1
    assert propagated.loc["Patient_1", "HP:0020219"] == 1
    assert propagated.loc["Patient_2", "HP:0020219"] == 1


def test_build_relationship_mask(hpo_engine, toy_terms):
    """
    Ancestor/descendant/self pairs should be NaN; unrelated pairs should be 0.
    """
    terms = list(toy_terms)

    # optional but harmless: prepare terms for term manager cache
    hpo_engine.propagate(
    pd.DataFrame([[np.nan] * len(terms)], columns=terms)
    )

    mask = hpo_engine.build_relationship_mask(terms)

    assert isinstance(mask, pd.DataFrame)
    assert mask.shape == (len(terms), len(terms))
    assert list(mask.index) == terms
    assert list(mask.columns) == terms

    # diagonal must be NaN
    assert all(pd.isna(mask.loc[t, t]) for t in terms)

    # HP:0001250 is ancestor of HP:0020219
    assert pd.isna(mask.loc["HP:0001250", "HP:0020219"])
    assert pd.isna(mask.loc["HP:0020219", "HP:0001250"])

    # if HP:0012759 is unrelated to the other two in your test ontology
    assert mask.loc["HP:0001250", "HP:0012759"] == 0
    assert mask.loc["HP:0012759", "HP:0001250"] == 0


def test_build_relationship_mask_is_symmetric(hpo_engine, toy_terms):
    """
    Relationship mask should be symmetric.
    """
    terms = list(toy_terms)
    mask = hpo_engine.build_relationship_mask(terms)

    pd.testing.assert_frame_equal(mask, mask.T)


def test_get_labels_returns_dict(hpo_engine):
    df = pd.DataFrame({"HP:0001250": [1]})
    hpo_engine.propagate(df)

    labels = hpo_engine.get_labels()

    assert isinstance(labels, dict)
    assert "HP:0001250" in labels
    assert labels["HP:0001250"]

