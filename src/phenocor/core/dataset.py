from dataclasses import dataclass, field

import pandas as pd
import phenopackets as ppkt


@dataclass
class PhenotypeDataset:
    """
    Store an analysis-ready phenotype dataset.

    The dataset contains a phenotype status matrix together with HPO term 
    labels, a pairwise ontology relationship mask, and the source Phenopackets.

    Phenotype matrix values follow the convention:

        1 = observed
        0 = excluded
        NaN = unknown

    Attributes 
    ---------- 
    phenotype_matrix : pd.DataFrame 
        Phenotype status matrix with individuals as rows and canonical HPO term 
        IDs as columns. 

    hpo_labels : dict[str, str] 
        Mapping from canonical HPO term IDs to their labels.

    relationship_mask : pd.DataFrame | None 
        Pairwise HPO relationship mask aligned to ``phenotype_matrix``. 
        Self-pairs and ancestor-descendant pairs are represented by ``NaN``; 
        unrelated pairs are represented by ``0``. 

    phenopackets : list[ppkt.Phenopacket] 
        Source Phenopackets corresponding to the individuals in the phenotype matrix.
    """

    phenotype_matrix: pd.DataFrame
    hpo_labels: dict[str, str] = field(default_factory=dict)
    relationship_mask: pd.DataFrame | None = None
    phenopackets: list[ppkt.Phenopacket] = field(default_factory=list)
        
    @property
    def individual_ids(self) -> pd.Index:
        """Return individual IDs in the phenotype matrix."""
        return self.phenotype_matrix.index

    @property
    def hpo_ids(self) -> pd.Index:
        """Return HPO term IDs in the phenotype matrix."""
        return self.phenotype_matrix.columns

    def get_pmids(self) -> pd.Series:
        """
        Return PubMed IDs associated with each individual. 
        
        PubMed IDs are extracted from Phenopacket metadata external 
        references whose identifiers begin with ``"PMID:"``.

        Returns
        -------
        pd.Series
            Series indexed by individual ID. Each value is a sorted list 
            of unique PubMed IDs without the ``"PMID:"`` prefix.
        """
        pmids_by_individual: dict[str, list[str]] = {}

        for phenopacket in self.phenopackets:
            individual_id = str(phenopacket.id)
            pmids_by_individual[individual_id] =  self._extract_pmids(phenopacket)

        return pd.Series(
            pmids_by_individual, 
            name="pmids", 
            dtype="object"
        ).reindex(self.individual_ids)

    def get_diseases(self) -> pd.Series:
        """
        Return observed disease IDs associated with each individual.

        Diseases explicitly marked as excluded are ignored.

        Returns
        -------
        pd.Series
            Series indexed by individual ID. Each value is a sorted list 
            of unique disease IDs.
        """
        diseases_by_individual: dict[str, list[str]] = {}

        for phenopacket in self.phenopackets:
            individual_id = str(phenopacket.id)
            disease_ids: set[str] = set()

            for disease in getattr(phenopacket, "diseases", []):
                if bool(getattr(disease, "excluded", False)):
                    continue

                term = getattr(disease, "term", None)

                disease_id = getattr(term, "id", None)

                if disease_id:
                    disease_ids.add(str(disease_id))

            diseases_by_individual[individual_id] = sorted(disease_ids)

        return pd.Series(
            diseases_by_individual,
            name="diseases",
            dtype="object",
        ).reindex(self.individual_ids)
        
    @staticmethod
    def _extract_pmids(phenopacket: ppkt.Phenopacket) -> list[str]:
        """
        Extract PubMed IDs from a Phenopacket metadata.

        Parameters
        ----------
        phenopacket : ppkt.Phenopacket
            Phenopacket from which PubMed IDs are extracted.  

        Returns
        -------
        list[str]
            Sorted unique PubMed IDs without the ``"PMID:"`` prefix.
        """ 
        meta_data = getattr(phenopacket, "meta_data", None)

        if meta_data is None:
            return []

        pmids: set[str] = set()
        
        for ref in getattr(meta_data, "external_references", []):
            ref_id = getattr(ref, "id", None)
            if isinstance(ref_id, str) and ref_id.startswith("PMID:"):
                pmids.add(ref_id.removeprefix("PMID:"))

        return sorted(pmids)
