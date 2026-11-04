from typing import IO, Any
from dataclasses import dataclass
import logging

import numpy as np
import pandas as pd
import phenopackets as ppkt

from .dataset import PhenotypeDataset
from ..ontology import HPOHierarchyEngine

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class _HpoObservation:
    """
    Store HPO observations for an individual.

    Attributes
    ----------
    individual_id : str
        Unique identifier for the individual.

    observed_terms : frozenset[str]
        HPO term IDs explicitly observed in the individual.

    excluded_terms : frozenset[str]
        HPO term IDs explicitly excluded in the individual.
    """

    individual_id: str
    observed_terms: frozenset[str]
    excluded_terms: frozenset[str]


class PhenotypeDatasetBuilder:
    """
    Build an analysis-ready ``PhenotypeDataset`` from Phenopackets.

    The builder extracts observed and excluded HPO annotations, constructs
    the phenotype status matrix, propagates annotations through the HPO
    hierarchy, and constructs the ontology relationship mask.
    """
    
    def __init__(
        self,
        phenopackets: list[ppkt.Phenopacket],
        *,
        hpo_file: str | IO[Any] | None = None,
        hpo_release: str | None = None,
    ) -> None:
        """
        Parameters
        ----------
        phenopackets : list[ppkt.Phenopacket]
            Input Phenopackets.

        hpo_file : str | IO[Any] | None, optional
            Local path, URL, or file-like object pointing to an HPO
            ontology file.

        hpo_release : str | None, optional
            HPO release identifier. If ``None`` and ``hpo_file`` is not
            provided, the latest available release is loaded.
        """
        self._validate_phenopackets(phenopackets)

        self.phenopackets = phenopackets
        self._hpo_engine = HPOHierarchyEngine(
            hpo_file=hpo_file,
            release=hpo_release,
        )

    def build(self) -> PhenotypeDataset:
        """
        Build an analysis-ready phenotype dataset.

        Returns
        -------
        PhenotypeDataset
            Dataset containing the propagated phenotype status matrix,
            HPO term labels, ontology relationship mask, and source
            Phenopackets.

        Raises
        ------
        ValueError
            If no HPO terms are found in the input Phenopackets or no valid
            HPO terms remain after ontology processing.
        """
        observations = self._parse_hpo_observations(self.phenopackets)
        raw_matrix = self._build_raw_hpo_matrix(observations)

        if raw_matrix.shape[1] == 0:
            raise ValueError(
                "No HPO terms found in phenopackets. "
                "Check `phenotypic_features`."
            )

        propagated_matrix = self._hpo_engine.propagate(raw_matrix)

        if propagated_matrix.shape[1] == 0:
            raise ValueError(
                "No valid HPO terms remain after ontology processing."
            )

        relationship_mask = (
            self._hpo_engine.build_relationship_mask(
                propagated_matrix.columns
            )
        )

        return PhenotypeDataset(
            phenotype_matrix=propagated_matrix,
            hpo_labels=self._hpo_engine.get_labels(),
            relationship_mask=relationship_mask,
            phenopackets=self.phenopackets,
        )

    @staticmethod
    def _validate_phenopackets(phenopackets: list[ppkt.Phenopacket]) -> None:
        """
        Validate input Phenopackets.

        The input must be non-empty, contain only Phenopacket objects,
        and contain non-empty unique Phenopacket IDs.

        Parameters
        ----------
        phenopackets : list[ppkt.Phenopacket]
            Phenopackets to validate.

        Raises
        ------
        TypeError
            If any element is not a Phenopacket.

        ValueError
            If the input is empty or Phenopacket IDs are missing or
            duplicated.
        """
        if not phenopackets:
            raise ValueError("Phenopackets cannot be empty")

        if not all(isinstance(p, ppkt.Phenopacket) for p in phenopackets):
            raise TypeError("`phenopackets` must contain only Phenopacket objects.")
        
        seen: set[str] = set()

        for phenopacket in phenopackets:
            individual_id = getattr(phenopacket, "id", None)

            if not individual_id:
                raise ValueError("Every phenopacket must have a non-empty `id`.")

            individual_id = str(individual_id)

            if individual_id in seen:
                raise ValueError(
                    f"Duplicate phenopacket id detected: {individual_id!r}. "
                    "Phenopacket ids define dataset row identities and must be unique."
                )

            seen.add(individual_id)

    @staticmethod
    def _parse_hpo_observations(phenopackets: list[ppkt.Phenopacket]) -> list[_HpoObservation]:
        """
        Extract observed and excluded HPO terms from Phenopackets.

        If the same HPO term is explicitly recorded as both observed and
        excluded for an individual, observed status takes precedence and
        a warning is emitted.

        Parameters
        ----------
        phenopackets : list[ppkt.Phenopacket]
            Phenopackets from which HPO annotations are extracted.

        Returns
        -------
        list[_HpoObservation]
            HPO observations for each individual.
        """
        observations: list[_HpoObservation] = []

        for phenopacket in phenopackets:
            individual_id = str(phenopacket.id)

            observed_terms: set[str] = set()
            excluded_terms: set[str] = set()

            for feature in phenopacket.phenotypic_features:
                term_id = getattr(getattr(feature, "type", None), "id", None)

                if not term_id:
                    continue

                term_id = str(term_id)

                if getattr(feature, "excluded", False):
                    excluded_terms.add(term_id)
                else:
                    observed_terms.add(term_id)

            conflicts = observed_terms & excluded_terms
            if conflicts:
                logger.warning(
                    "Individual %s has conflicting HPO annotations. "
                    "Observed takes precedence. Conflicts: %s",
                    individual_id,
                    sorted(conflicts),
                )
                excluded_terms -= conflicts

            observations.append(
                _HpoObservation(
                    individual_id=individual_id,
                    observed_terms=frozenset(observed_terms),
                    excluded_terms=frozenset(excluded_terms),
                )
            )

        return observations


    @staticmethod
    def _build_raw_hpo_matrix(observations: list[_HpoObservation]) -> pd.DataFrame:
        """
        Construct the raw phenotype status matrix.

        Parameters
        ----------
        observations : list[_HpoObservation]
            HPO observations for each individual.

        Returns
        -------
        pd.DataFrame
            Phenotype status matrix with individuals as rows and HPO term
            IDs as columns. Values are ``1`` for observed, ``0`` for
            excluded, and ``NaN`` for unknown.
        """
        hpo_ids = sorted(
            set().union(
                *(obs.observed_terms | obs.excluded_terms for obs in observations)
            )
        )

        individual_ids = [obs.individual_id for obs in observations]

        matrix = np.full(
            shape=(len(observations), len(hpo_ids)),
            fill_value=np.nan,
            dtype=float,
        )

        hpo_id_to_column = {hpo_id: col for col, hpo_id in enumerate(hpo_ids)}

        for row, observation in enumerate(observations):
            for hpo_id in observation.excluded_terms:
                matrix[row, hpo_id_to_column[hpo_id]] = 0.0

            for hpo_id in observation.observed_terms:
                matrix[row, hpo_id_to_column[hpo_id]] = 1.0

        return pd.DataFrame(
            matrix,
            index=pd.Index(individual_ids, name="individual_id"),
            columns=pd.Index(hpo_ids, name="hpo_id"),
        )
    
