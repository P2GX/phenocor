from typing import IO, Any
import logging

import hpotk

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _load_ontology(
    file: str | IO[Any] | None = None,
    release: str | None = None
) -> hpotk.MinimalOntology:
    """
    Load a minimal Human Phenotype Ontology (HPO) instance.

    The ontology is loaded from ``file`` when provided. Otherwise, it is 
    loaded from the configured ontology store using ``release``. If neither 
    is provided, the latest available HPO release is loaded.

    Parameters
    ----------
    file : str | IO[Any] | None, optional
        Local path, URL, or file-like object pointing to an HPO ontology file, such as ``hp.json``.

    release : str | None, optional
        HPO release identifier, for example ``"v2023-10-09"``. 
        Used only when ``file`` is ``None``.

    Returns
    -------
    hpotk.MinimalOntology
        Loaded HPO ontology.
    """
    if file is not None:
        ontology = hpotk.load_minimal_ontology(file)
        logger.info("Loaded HPO ontology from file: %s", file)
        return ontology
    
    store = hpotk.configure_ontology_store()
    
    if release is not None:
        ontology = store.load_minimal_hpo(release=release)
        logger.info("Loaded HPO ontology from release: %s", release)
        return ontology

    ontology = store.load_minimal_hpo()

    logger.info("Loaded HPO ontology from the latest available release.")

    return ontology


class HPOTermManager:
    """
    Resolve HPO term IDs and cache ontology information. 
    
    This internal helper provides access to the loaded HPO ontology and 
    maintains caches for canonical term IDs, labels, ancestors, and descendants. 
    
    Term IDs supplied to :meth:`prepare_terms` are resolved to their canonical 
    primary identifiers. Invalid term IDs are skipped with a warning.
    
    Notes 
    ----- 
    This class is intended for internal use by ``HPOHierarchyEngine``. 
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
            Local path, URL, or file-like object pointing to an HPO ontology file.

        release : str | None, optional
            HPO release identifier. If ``None`` and ``hpo_file`` is not
            provided, the latest available release is loaded.
        """
        self._hpo = _load_ontology(file=hpo_file, release=release)

        self._ancestor_cache: dict[str, set[str]] = {}
        self._descendant_cache: dict[str, set[str]] = {}
        self._label_cache: dict[str, str] = {}
        self._id_mapping_cache: dict[str, str] = {}

    @property
    def hpo(self) -> hpotk.MinimalOntology:
        """Return the underlying HPO ontology."""
        return self._hpo

    def resolve_term_id(
        self, 
        term_id: str,
    ) -> str:
        """
        Resolve an HPO term ID to its canonical primary identifier.

        Parameters
        ----------
        term_id : str
            HPO term ID, for example ``"HP:0001250"``.

        Returns
        -------
        str
            Canonical primary HPO term ID.

        Raises
        ------
        ValueError
            If ``term_id`` cannot be resolved in the loaded ontology.
        """
        term = self._hpo.get_term(term_id=term_id)

        if term is None:
            raise ValueError(f"HPO term {term_id!r} was not found in the ontology.")
        
        return term.identifier.value

    def prepare_terms(
        self, 
        terms: set[str],
    ) -> set[str]:
        """
        Resolve HPO term IDs and cache their ontology information.
        
        Each input term is resolved to its canonical primary identifier. 
        For valid terms, the method caches:

        - original-to-canonical ID mappings, 
        - HPO labels, 
        - ancestor IDs, 
        - descendant IDs.

        Invalid terms are skipped with a warning.

        Parameters
        ----------
        terms : set[str]
            HPO term IDs to prepare.

        Returns
        -------
        set[str]
            Canonical HPO term IDs corresponding to valid input terms.
        """
        canonical_hpo_ids: set[str] = set()

        for hpo_id in terms:
            try:
                canonical_hpo_id = self.resolve_term_id(hpo_id)

            except ValueError as exc:
                logger.warning(
                    "Skipping invalid HPO term %r: %s",
                    hpo_id,
                    exc,
                )
                continue

            if canonical_hpo_id != hpo_id:
                logger.info("Resolved term %r to canonical ID %r", hpo_id, canonical_hpo_id)

            canonical_hpo_ids.add(canonical_hpo_id)

            self._id_mapping_cache[hpo_id] = canonical_hpo_id

            if canonical_hpo_id not in self._label_cache:
                term = self._hpo.get_term(term_id=canonical_hpo_id)
                self._label_cache[canonical_hpo_id] = term.name if term is not None else ""

        uncached_hpo_ids = canonical_hpo_ids - self._ancestor_cache.keys()

        for hpo_id in uncached_hpo_ids:
            try:
                self._ancestor_cache[hpo_id] = {
                    ancestor.value for ancestor in self._hpo.graph.get_ancestors(hpo_id)
                }
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to fetch ancestors for HPO term {hpo_id!r}."
                ) from exc

            try:
                self._descendant_cache[hpo_id] = {
                    descendant.value for descendant in self._hpo.graph.get_descendants(hpo_id)
                }
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to fetch descendants for HPO term {hpo_id!r}."
                ) from exc

        return canonical_hpo_ids

    def get_ancestors(self, hpo_id: str) -> set[str]:
        """Return cached ancestor terms for a canonical HPO term."""
        return self._ancestor_cache.get(hpo_id, set())
    
    def get_descendants(self, hpo_id: str) -> set[str]:
        """Return cached descendant terms for a canonical HPO term."""
        return self._descendant_cache.get(hpo_id, set())

    def get_labels(self) -> dict[str, str]:
        """Return cached labels for canonical HPO terms."""
        return dict(self._label_cache)

    def get_id_mapping(self) -> dict[str, str]:
        """Return cached mappings from original to canonical HPO terms."""
        return dict(self._id_mapping_cache)