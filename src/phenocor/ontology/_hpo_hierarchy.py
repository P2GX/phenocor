from typing import IO, Any
from collections.abc import Sequence
from collections import defaultdict
import logging

import numpy as np
import pandas as pd
import hpotk

from ._term_manager import HPOTermManager


logger = logging.getLogger(__name__)


class HPOHierarchyEngine:
    """
    Apply HPO hierarchy-aware operations to phenotype status matrices.

    The engine provides two main operations:

    - propagation of observed and excluded phenotype annotations through 
      the HPO hierarchy;
    - construction of a pairwise mask for excluding ontologically related 
      HPO terms from downstream correlation analysis.

    Phenotype status matrices are expected to contain:
    
        1 = observed
        0 = excluded
        NaN = unknown

    Before propagation, HPO term IDs are resolved to canonical primary IDs. 
    Invalid HPO terms are removed, and columns resolving to the same 
    canonical HPO term are merged.

    Notes 
    ----- 
    Propagation is restricted to HPO terms represented in the input matrix. 
    Observed annotations propagate to represented ancestors, while excluded 
    annotations propagate to represented descendants. 
    
    Explicit annotations are preserved when they conflict with propagated 
    annotations.
    """

    def __init__(
        self,
        hpo_file: str | IO[Any] | None = None,
        release: str | None = None,
    ) -> None:
        """ 
        Parameters 
        ---------- 
        hpo_file : str | IO[Any] | None, optional 
            Local path, URL, or file-like object pointing to an HPO 
            ontology file. 
        
        release : str | None, optional  
            HPO release identifier. If ``None`` and ``hpo_file`` is not 
            provided, the latest available release is loaded. 
        """
        self._term_manager = HPOTermManager(hpo_file=hpo_file, release=release)

    @property
    def hpo(self) -> hpotk.MinimalOntology:
        """Return the underlying HPO ontology."""
        return self._term_manager.hpo

    def propagate(
        self, 
        matrix: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Propagate phenotype annotations through the ontology hierarchy.

        Observed annotations (``1``) propagate upward to represented 
        ancestor terms. Excluded annotations (``0``) propagate downward 
        to represented descendant terms.

        Before propagation, HPO term IDs are resolved to canonical primary 
        IDs. Invalid terms are removed. If multiple input columns resolve 
        to the same canonical HPO term, they are merged. When conflicting 
        values occur during this merge, observed status (``1``) takes 
        precedence and a warning is emitted.

        Propagation only fills previously unknown values (``NaN``). 
        Therefore, explicit observed or excluded annotations are never 
        overwritten by propagated values. Conflicts between explicit and 
        propagated annotations are reported as warnings.

        When both observed and excluded states are inferred for an unknown 
        value, observed status takes precedence and a warning is emitted.

        Parameters
        ----------
        matrix : pd.DataFrame
            Phenotype status matrix with individuals as rows and HPO term 
            IDs as columns. Values must be ``1`` (observed), ``0`` 
            (excluded), or ``NaN`` (unknown).

        Returns
        -------
        pd.DataFrame
            Phenotype status matrix with canonical HPO term IDs and 
            hierarchy-propagated annotations.
        """
        if matrix.empty:
            return matrix.copy()

        input_hpo_ids = set(matrix.columns)
        canonical_hpo_ids = self._term_manager.prepare_terms(input_hpo_ids)

        self._check_mapping_conflicts(matrix)

        id_mapping = self._term_manager.get_id_mapping()

        matrix = matrix.rename(columns=id_mapping)

        # Multiple input IDs may resolve to the same canonical HPO term. 
        # max() gives precedence to observed (1) over excluded (0).
        matrix = matrix.T.groupby(level=0).max().T
        matrix = matrix.loc[:, matrix.columns.isin(canonical_hpo_ids)]
        hpo_ids = set(matrix.columns)

        ancestors = {t: self._term_manager.get_ancestors(t) for t in hpo_ids}
        descendants = {t: self._term_manager.get_descendants(t) for t in hpo_ids}

        propagated_observed: dict[str, pd.Series] = {}
        propagated_excluded: dict[str, pd.Series] = {}

        for hpo_id in hpo_ids:     
            observed_mask = matrix[hpo_id] == 1
            if observed_mask.any():
                for ancestor_id in ancestors.get(hpo_id, []):
                    if ancestor_id not in hpo_ids:
                        continue
                    if ancestor_id in propagated_observed:
                        propagated_observed[ancestor_id] = propagated_observed[ancestor_id] | observed_mask
                    else:
                        propagated_observed[ancestor_id] = observed_mask.copy()

            excluded_mask = matrix[hpo_id] == 0

            if excluded_mask.any():
                for descendant_id in descendants.get(hpo_id, []):
                    if descendant_id not in hpo_ids:
                        continue
                    if descendant_id in propagated_excluded:
                        propagated_excluded[descendant_id] = propagated_excluded[descendant_id] | excluded_mask
                    else:
                        propagated_excluded[descendant_id] = excluded_mask.copy()

        for hpo_id in hpo_ids:
            observed_propagation = propagated_observed.get(hpo_id, None)
            excluded_propagation = propagated_excluded.get(hpo_id, None)
            
            explicit_observed = matrix[hpo_id] == 1
            explicit_excluded = matrix[hpo_id] == 0

            if (
                observed_propagation is not None
                and excluded_propagation is not None
            ):
                conflicts = (
                    observed_propagation
                    & excluded_propagation
                    & matrix[hpo_id].isna()
                )
                if conflicts.any():
                    logger.warning(
                        "HPO propagation conflict: %d individual(s) "
                        "have both OBSERVED (1) and EXCLUDED (0) "
                        "inferred for %r. Preserving OBSERVED (1).",
                        int(conflicts.sum()),
                        hpo_id,
                    )

            if observed_propagation is not None:
                conflicts = observed_propagation & explicit_excluded
                if conflicts.any():
                    logger.warning(
                        "HPO propagation conflict: %d individual(s) " 
                        "have a descendant term implying OBSERVED (1) " 
                        "for %r, but the term is explicitly marked as " 
                        "EXCLUDED (0). Preserving the explicit value.",
                        int(conflicts.sum()), hpo_id
                    )
                update_mask = observed_propagation & matrix[hpo_id].isna()
                matrix.loc[update_mask, hpo_id] = 1

            if excluded_propagation is not None:
                conflicts = excluded_propagation & explicit_observed
                if conflicts.any():
                    logger.warning(
                        "HPO propagation conflict: %d individual(s) " 
                        "have an ancestor term implying EXCLUDED (0) " 
                        "for %r, but the term is explicitly marked as " 
                        "OBSERVED (1). Preserving the explicit value.",
                        int(conflicts.sum()), hpo_id
                    )
                update_mask = excluded_propagation & matrix[hpo_id].isna()
                matrix.loc[update_mask, hpo_id] = 0

        return matrix

    def _check_mapping_conflicts(
        self, 
        matrix: pd.DataFrame
    ) -> None:
        """
        Detect conflicts before merging HPO IDs to canonical terms.

        Multiple input HPO IDs may resolve to the same canonical term. 
        This method checks whether such columns contain both observed 
        (``1``) and excluded (``0``) values for the same individual.

        When a conflict is detected, a warning is emitted. The subsequent 
        merge performed in :meth:`propagate` gives precedence to observed 
        status (``1``).

        Parameters
        ----------
        matrix : pd.DataFrame
            Phenotype status matrix before HPO ID canonicalization.
        """
        id_mapping = self._term_manager.get_id_mapping()
        
        target_to_sources: dict[str, list[str]] = defaultdict(list)

        for source_id, target_id in id_mapping.items():
            if source_id in matrix.columns:
                target_to_sources[target_id].append(source_id)
                
        multi_mappings = {
            target_id: source_ids for target_id, source_ids in target_to_sources.items() if len(source_ids) > 1
        }
        
        if not multi_mappings:
            return  

        for target_id, source_ids in multi_mappings.items():
            values = matrix[source_ids]
            
            has_observed = (values == 1).any(axis=1)
            has_excluded = (values == 0).any(axis=1)

            conflict_mask = has_observed & has_excluded
            
            if conflict_mask.any():
                conflict_individual_ids = matrix.index[conflict_mask].tolist()
                
                logger.warning(
                    "HPO ID mapping conflict: input terms %s all resolve " 
                    "to canonical term %r, but individual(s) %s contain " 
                    "both OBSERVED (1) and EXCLUDED (0) values. " 
                    "The canonicalized matrix will preserve OBSERVED (1).",
                    source_ids,
                    target_id,
                    conflict_individual_ids if len(conflict_individual_ids) <= 5 else f"{conflict_individual_ids[:5]}... (Total: {len(conflict_individual_ids)})"
                )
    
    def build_relationship_mask(
        self, 
        hpo_ids: Sequence[str]
    ) -> pd.DataFrame:
        """
        Build a pairwise mask for ontologically related HPO terms.

        Self-pairs and ancestor-descendant pairs are marked as ``NaN``. 
        Unrelated HPO term pairs are marked as ``0``.

        The mask is intended to exclude structurally related HPO term
        pairs from downstream correlation analysis.

        This method is intended for internal use after ``propagate()``
        has prepared the ontology caches for the supplied terms.

        Parameters
        ----------
        hpo_ids : Sequence[str]
            Canonical HPO term IDs, typically corresponding to columns of 
            a propagated phenotype matrix. Input order is preserved.

        Returns
        -------
        pd.DataFrame
            Square matrix indexed and columned by ``hpo_ids``: 
            
            - ``NaN`` indicates a self-pair or ancestor-descendant pair; 
            - ``0`` indicates a pair eligible for downstream evaluation
              with respect to ontology structure.
        """
        hpo_ids = list(hpo_ids)
        hpo_id_to_index = {hpo_id: index for index, hpo_id in enumerate(hpo_ids)}
        n_hpo_terms = len(hpo_ids)
        mask = np.zeros((n_hpo_terms, n_hpo_terms), dtype=float)

        for index, hpo_id in enumerate(hpo_ids):
            related_hpo_ids = self._term_manager.get_ancestors(hpo_id) | self._term_manager.get_descendants(hpo_id)
            related_indices = [hpo_id_to_index[t] for t in related_hpo_ids if t in hpo_id_to_index]
            mask[index, related_indices] = np.nan
            mask[related_indices, index] = np.nan
            mask[index, index] = np.nan

        return pd.DataFrame(mask, index=hpo_ids, columns=hpo_ids)

    def get_labels(self) -> dict[str, str]:
        """Return labels for HPO terms prepared by this engine."""
        return self._term_manager.get_labels()
