# phenocor

**phenocor** is a Python library for analyzing pairwise relationships between Human Phenotype Ontology (HPO) features in [GA4GH Phenopacket](https://www.ga4gh.org/product/phenopackets/) cohorts.

It provides **Pairwise association analysis**, which measures associations between binary phenotype annotations using the phi coefficient and Fisher's exact test.


---

## Installation

```bash
pip install phenocor
```


## Features

* Construction of phenotype matrices from GA4GH Phenopackets, with separate representation of observed, excluded, and unreported HPO annotations
* Ontology-aware propagation of phenotype annotations
* Pairwise phenotype association analysis using the phi coefficient, Fisher's exact test, and Benjamini–Hochberg correction, with exclusion of ancestor–descendant HPO term pairs
* Tabular results and interactive heatmaps, including contingency counts, effective sample sizes, adjusted p-values, and publication provenance where available


## Quickstart

```python
from pathlib import Path
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket
from phenocor import (
    PhenotypeDatasetBuilder,
    HPOCorrelationAnalyzer,
)

# Load phenopackets
phenopacket_dir = Path("path/to/your/fbn1_phenopackets/")

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
```

For complete workflows, visualization options, and API details, see the [Documentation](https://phenocor.readthedocs.io/).



