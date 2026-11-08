import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import scipy.stats
from joblib import Parallel, delayed
from scipy.sparse import coo_matrix, triu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from ..core import PhenotypeDataset


logger = logging.getLogger(__name__)


class CorrelationResult:
    """
    Store and summarize pairwise HPO correlation results.

    This class contains the long-format pairwise association results,
    the symmetric Phi coefficient matrix, the corresponding raw p-value
    matrix, and optional HPO term labels.

    Statistical calculations are performed by ``HPOCorrelationAnalyzer``.
    This class provides utilities for filtering, exporting, and visualizing
    the resulting associations.

    Parameters
    ----------
    correlation_results : pd.DataFrame
        Long-format table containing pairwise HPO association results.

    coef_matrix : pd.DataFrame
        Symmetric matrix of Phi coefficients.

    pval_matrix : pd.DataFrame
        Symmetric matrix of raw Fisher's exact test p-values.

    hpo_labels : dict[str, str] | None, optional
        Mapping from HPO term IDs to human-readable labels.
    """

    def __init__(
        self, 
        correlation_results:pd.DataFrame, 
        coef_matrix: pd.DataFrame, 
        pval_matrix: pd.DataFrame, 
        hpo_labels: dict[str, str] | None = None
    ) -> None:
        self.correlation_results = correlation_results      
        self.coef_matrix = coef_matrix           
        self.pval_matrix = pval_matrix         
        self.hpo_labels = hpo_labels or {}

    def _filter_results(
        self,
        *,
        adj_pval_threshold: float | None = 0.05,
        min_abs_correlation: float | None = None,
    ) -> pd.DataFrame:
        """
        Filter pairwise correlation results by statistical significance
        and correlation strength.

        Parameters
        ----------
        adj_pval_threshold : float | None, default=0.05
            Adjusted p-value threshold. If ``None``, no significance
            filtering is applied.

        min_abs_correlation : float | None, optional
            Minimum absolute Phi coefficient. If ``None``, no
            correlation-strength filtering is applied.

        Returns
        -------
        pd.DataFrame
            Filtered pairwise correlation results.

        Raises
        ------
        ValueError
            If ``adj_pval_threshold`` or ``min_abs_correlation`` is outside its
            valid range.
        """
        if adj_pval_threshold is not None and not 0.0 < adj_pval_threshold <= 1.0:
            raise ValueError(
                "`adj_pval_threshold` must be greater than 0 and at most 1."
            )

        if (
            min_abs_correlation is not None
            and not 0.0 <= min_abs_correlation <= 1.0
        ):
            raise ValueError(
                "`min_abs_correlation` must be between 0 and 1."
            )

        results = self.correlation_results.copy()

        if adj_pval_threshold is not None:
            results = results[
                results["adj_p_value"] < adj_pval_threshold
            ]

        if min_abs_correlation is not None:
            results = results[
                results["correlation"].abs()
                >= min_abs_correlation
            ]

        return results.reset_index(drop=True)

    def get_significant_correlations(
        self,
        adj_pval_threshold: float = 0.05,
        *,
        min_abs_correlation: float | None = None,
    ) -> pd.DataFrame:
        """
        Return statistically significant pairwise correlations.

        Parameters
        ----------
        adj_pval_threshold : float, default=0.05
            Adjusted p-value significance threshold.

        min_abs_correlation : float | None, optional
            Minimum absolute Phi coefficient.

        Returns
        -------
        pd.DataFrame
            Significant pairwise correlations satisfying the specified
            filters.
        """
        return self._filter_results(
            adj_pval_threshold=adj_pval_threshold,
            min_abs_correlation=min_abs_correlation,
        )

    def save_correlation_results(
        self,
        output_file: str | Path = "correlation_results.csv",
        *,
        adj_pval_threshold: float | None = 0.05,
        min_abs_correlation: float | None = None,
    ) -> None:
        """
        Save pairwise correlation results.

        By default, only correlations with adjusted p-values below 0.05
        are saved. If ``adj_pval_threshold`` is ``None``, no significance filtering
        is applied.

        Parameters
        ----------
        output_file : str | Path, default="correlation_results.csv"
            Output file path. Supported formats are ``.csv`` and
            ``.xlsx``.

        adj_pval_threshold : float | None, default=0.05
            Adjusted p-value threshold. If ``None``, no significance
            filtering is applied.

        min_abs_correlation : float | None, optional
            Minimum absolute Phi coefficient. If ``None``, no
            correlation-strength filtering is applied.

        Raises
        ------
        ValueError
            If the output format is unsupported.
        """
        results = self._filter_results(
            adj_pval_threshold=adj_pval_threshold,
            min_abs_correlation=min_abs_correlation,
        )

        output_path = Path(output_file)
        extension = output_path.suffix.lower()

        if extension == ".csv":
            results.to_csv(
                output_path,
                index=False,
            )

        elif extension == ".xlsx":
            results.to_excel(
                output_path,
                index=False,
            )

        else:
            raise ValueError(
                "Unsupported output format. "
                "Use '.csv' or '.xlsx'."
            )

    @staticmethod
    def _format_hpo_pair( 
        hpo_id: str, 
        label: str | None
    ) -> str:
        """Format an HPO term for display."""
        if label:
            return f"{label} ({hpo_id})"
        return hpo_id
    
    @staticmethod
    def _format_pmids_for_tooltip(
        pmids: str | list[str] | None,
        max_pmids: int = 5,
    ) -> str:
        """Format PMID values for hover text."""
        if pmids is None or pmids == "":
            return "None"

        if isinstance(pmids, str):
            pmid_list = [p.strip() for p in pmids.split(";") if p.strip()]
        else:
            pmid_list = [str(p).strip() for p in pmids if str(p).strip()]

        if not pmid_list:
            return "None"

        if len(pmid_list) <= max_pmids:
            return ", ".join(pmid_list)

        shown = ", ".join(pmid_list[:max_pmids])
        remaining = len(pmid_list) - max_pmids
        return f"{shown} ... (+{remaining} more)"

    def plot_correlation_heatmap_with_significance(
        self,
        *,
        adj_pval_threshold: float | None = 0.05,
        min_abs_correlation: float | None = 0.1,
        title: str | None = None,
    ) -> go.Figure:
        """
        Create an interactive heatmap of pairwise HPO correlations.

        By default, only correlations with adjusted p-values below 0.05
        are displayed. Correlations can additionally be filtered by
        absolute Phi coefficient.

        Parameters
        ----------
        adj_pval_threshold : float | None, default=0.05
            Adjusted p-value threshold. If ``None``, no significance
            filtering is applied.

        min_abs_correlation : float | None, default=0.1
            Minimum absolute Phi coefficient. If ``None``, no
            correlation-strength filtering is applied.

        title : str | None, optional
            Optional title displayed for the heatmap.

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive correlation heatmap.

        Raises
        ------
        ValueError
            If no correlations remain after filtering.
        """
        results = self._filter_results(
            adj_pval_threshold=adj_pval_threshold,
            min_abs_correlation=min_abs_correlation,
        )

        if results.empty:
            raise ValueError(
                "No correlations remain after filtering."
            )


        hpo_ids = pd.Index(
            pd.unique(
                pd.concat(
                    [
                        results["HPO_A"],
                        results["HPO_B"],
                    ],
                    ignore_index=True,
                )
            )
        )

        coef_matrix = self.coef_matrix.loc[
            hpo_ids,
            hpo_ids,
        ].copy()

        pval_matrix = self.pval_matrix.loc[
            hpo_ids,
            hpo_ids,
        ].copy()

        included_pairs = {
            frozenset((row["HPO_A"], row["HPO_B"]))
            for _, row in results.iterrows()
        }

        for hpo_a in hpo_ids:
            for hpo_b in hpo_ids:
                if (
                    hpo_a == hpo_b
                    or frozenset((hpo_a, hpo_b))
                    not in included_pairs
                ):
                    coef_matrix.loc[hpo_a, hpo_b] = np.nan
                    pval_matrix.loc[hpo_a, hpo_b] = np.nan

        n_rows, n_cols = coef_matrix.shape
        cell_size = 60  # Base pixel size per cell
        max_dim = max(n_rows, n_cols)
        fig_size = min(1200, max_dim * cell_size)  # Cap total figure size to avoid excessive width

        title_fontsize = max(14 + max_dim // 2, 28)
        label_fontsize = max(8, 12 - max_dim // 8)
        annot_fontsize = max(6, 12 - max_dim // 8)

        triangle_mask = pd.DataFrame(
            np.tril(np.ones(coef_matrix.shape, dtype=bool), k=0),
            index=coef_matrix.index,
            columns=coef_matrix.columns
        )
        coef_matrix = coef_matrix.where(triangle_mask)
        pval_matrix = pval_matrix.where(triangle_mask)

        nan_bg = pd.DataFrame(np.nan, index=coef_matrix.index, columns=coef_matrix.columns)
        nan_bg[triangle_mask & coef_matrix.isna()] = 2

        text_matrix = np.where(
            np.isnan(coef_matrix.values),
            "",
            coef_matrix.round(2).astype(str)
        )

        result_lookup: dict[tuple[str, str], dict[str, object]] = {}
        for _, hpo_a in results.iterrows():
            forward = {
                "Coefficient": hpo_a["correlation"],
                "P_value": hpo_a["p_value"],
                "P_value_corrected": hpo_a.get("adj_p_value", None), 
                "Count_00": hpo_a["n(A:E/B:E)"],
                "Count_01": hpo_a["n(A:E/B:O)"],
                "Count_10": hpo_a["n(A:O/B:E)"],
                "Count_11": hpo_a["n(A:O/B:O)"],
                "n_individuals": hpo_a["n_individuals"],
            }

            backward = {
                "Coefficient": hpo_a["correlation"],
                "P_value": hpo_a["p_value"],
                "P_value_corrected": hpo_a.get("adj_p_value", None),
                "Count_00": hpo_a["n(A:E/B:E)"],
                "Count_01": hpo_a["n(A:O/B:E)"],  # swapped
                "Count_10": hpo_a["n(A:E/B:O)"],  # swapped
                "Count_11": hpo_a["n(A:O/B:O)"],
                "n_individuals": hpo_a["n_individuals"],
            }

            if "n_pmids" in hpo_a.index:
                forward["n_pmids"] = hpo_a["n_pmids"]
                forward["pmids"] = hpo_a.get("pmids", "")
                backward["n_pmids"] = hpo_a["n_pmids"]
                backward["pmids"] = hpo_a.get("pmids", "")

            result_lookup[(hpo_a["HPO_A"], hpo_a["HPO_B"])] = forward
            result_lookup[(hpo_a["HPO_B"], hpo_a["HPO_A"])] = backward

        hover_text:list[list[str]] = []

        for i, hpo_a in enumerate(coef_matrix.index):
            hover_row: list[str] = []
            for j, hpo_b in enumerate(coef_matrix.columns):
                correlation = coef_matrix.iloc[i, j]
                p_value = pval_matrix.iloc[i, j]

                display_a = self._format_hpo_pair(hpo_a, self.hpo_labels.get(hpo_a))
                display_b = self._format_hpo_pair(hpo_b, self.hpo_labels.get(hpo_b))

                if not triangle_mask.iloc[i, j] or np.isnan(correlation):
                    hover_row.append("")
                    continue

                stats = result_lookup.get((hpo_a, hpo_b), {})
                pmid_block = ""

                if "n_pmids" in stats:
                    pmid_text = self._format_pmids_for_tooltip(
                        stats.get("pmids", ""),
                        max_pmids=4,
                    )
                    pmid_block = (
                        f"<b>N_PMIDs</b>: {int(stats.get('n_pmids', 0))}<br>"
                        f"<b>PMIDs</b>: {pmid_text}"
                    )
                hover_row.append(
                    f"<b>HPO_A</b>: {display_a}<br>"
                    f"<b>HPO_B</b>: {display_b}<br>"
                    f"<b>Corr</b>: {correlation:.2f}<br><b>p-val</b>: {p_value:.6f}<br>"
                    f"<b>adj_p_val</b>: {stats.get('P_value_corrected', np.nan):.6f}<br>"
                    f"<b>Counts(A/B)</b>: "
                    f"<b>E/E</b>: {stats.get('Count_00', 0)}, "
                    f"<b>E/O</b>: {stats.get('Count_01', 0)}, "
                    f"<b>O/E</b>: {stats.get('Count_10', 0)}, "
                    f"<b>O/O</b>: {stats.get('Count_11', 0)}<br>"
                    f"<b>Total_individuals</b>: {stats.get('n_individuals', 0)}<br>"
                    f"{pmid_block}"
                )
            hover_text.append(hover_row)
        
        display_labels = {
            hpo_id: self.hpo_labels.get(
                hpo_id,
                hpo_id,
            )
            for hpo_id in hpo_ids
        }

        coef_matrix = coef_matrix.rename(
            index=display_labels,
            columns=display_labels,
        )

        nan_bg = nan_bg.rename(
            index=display_labels,
            columns=display_labels,
        )

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=nan_bg.to_numpy(),
            x=coef_matrix.columns,
            y=coef_matrix.index,
            colorscale=[[0, "#eef4fb"], [1, "#eef4fb"]],
            showscale=False,
            hoverinfo="skip",
            xgap=1,
            ygap=1,
        ))
        fig.add_trace(go.Heatmap(
                z=coef_matrix.to_numpy(),
                x=coef_matrix.columns,
                y=coef_matrix.index,
                colorscale=[ 
                [0.00, "#203864"],   # navy
                [0.50, "#F7F4ED"],   # ivory
                [1.00, "#7A1F3D"]    # wine
                ],
                zmin=-1,
                zmax=1,
                zmid=0,
                text=text_matrix,
                texttemplate=f"<span style='font-size:{annot_fontsize}px'>%{{text}}</span>",
                hovertext=hover_text,
                hoverinfo="text",
                colorbar=dict(title="Corr.", len=0.8, thickness=title_fontsize),
                xgap=1,
                ygap=1,
            ))
        
        max_ylabel_len = max(len(str(lbl)) for lbl in coef_matrix.index) if not coef_matrix.empty else 10
        left_margin = 60 + max_ylabel_len * label_fontsize

        clean_subtitle = title.strip() if title and title.strip() else ""

        main_title = "<b>Phi Coefficient Matrix for HPO Pairwise Associations</b>"
        
        full_title = f"{main_title}<br><span style='font-size:0.8em'>{clean_subtitle}</span>" if clean_subtitle else main_title

        fig.update_layout(
            title=dict(
                text=full_title,
                x=0.5,
                xanchor="center",
                yanchor="top",
                font=dict(
                    size=min(title_fontsize, 24),
                    family="Arial"
                )
            ),
            xaxis=dict(
                tickangle=90,
                tickfont=dict(size=label_fontsize),
            ),
            yaxis=dict(
                tickfont=dict(size=label_fontsize),
                scaleanchor="x",
                scaleratio=1
            ),
            width=fig_size + left_margin,
            height=fig_size + left_margin,
            plot_bgcolor="white",
            paper_bgcolor="white"
        )
        fig.update_yaxes(autorange="reversed")
        self.fig = fig
        return fig
    
    def save_correlation_heatmap(
        self,
        output_file: str | Path = "correlation_heatmap.html",
        *,
        adj_pval_threshold: float | None = 0.05,
        min_abs_correlation: float | None = None,
        title: str | None = None,
        width: int | None = None,
        height: int | None = None,
        scale: float = 2.0,
    ) -> None:
        """
        Save the correlation heatmap as HTML or a static image.

        Parameters
        ----------
        output_file : str | Path, default="correlation_heatmap.html"
            Output path. Supported extensions are ``.html``, ``.png``,
            ``.svg``, ``.pdf``, and ``.jpg``.

        adj_pval_threshold : float | None, default=0.05
            Adjusted p-value threshold. If ``None``, no significance
            filtering is applied.

        min_abs_correlation : float | None, optional
            Minimum absolute Phi coefficient.

        title : str | None, optional
            Optional subtitle displayed below the main heatmap title.

        width : int | None, optional
            Export width in pixels for static image formats.

        height : int | None, optional
            Export height in pixels for static image formats.

        scale : float, default=2.0
            Resolution multiplier for raster image formats.

        Raises
        ------
        ValueError
            If the output format is unsupported.
        """
        output_path = Path(output_file)
        extension = output_path.suffix.lower()

        supported_extensions = {
            ".html",
            ".png",
            ".svg",
            ".pdf",
            ".jpg",
        }

        if extension not in supported_extensions:
            raise ValueError(
                "Unsupported output format. "
                "Use one of: .html, .png, .svg, .pdf, .jpg."
            )

        figure = self.plot_correlation_heatmap_with_significance(
            adj_pval_threshold=adj_pval_threshold,
            min_abs_correlation=min_abs_correlation,
            title=title,
        )

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        if extension == ".html":
            figure.write_html(
                str(output_path)
            )

        else:
            figure.write_image(
                str(output_path),
                width=width,
                height=height,
                scale=scale,
            )

class HPOCorrelationAnalyzer:
    """
    Analyze pairwise associations between HPO terms.

    Pairwise associations are quantified using the Phi coefficient.
    Statistical significance is assessed using Fisher's exact test,
    followed by Benjamini-Hochberg correction for multiple testing.

    Only ontology-unrelated HPO term pairs with sufficient pairwise
    non-missing observations are evaluated.
    """

    def __init__(
        self,  
        dataset: PhenotypeDataset, 
    ) -> None:
        """
        Parameters
        ----------
        dataset : PhenotypeDataset
            Analysis-ready phenotype dataset.
        """
        if not isinstance(dataset, PhenotypeDataset):
            raise TypeError("`dataset` must be a `PhenotypeDataset` instance.")

        if dataset.relationship_mask is None:
            raise ValueError(
                "`dataset.relationship_mask` is required for ontology-aware "
                "analysis. Build the dataset with `PhenotypeDatasetBuilder`, "
                "or provide a relationship mask aligned to "
                "`dataset.phenotype_matrix`."
            )
 
        mask = dataset.relationship_mask
        if not (
            mask.index.equals(dataset.hpo_ids)
            and mask.columns.equals(dataset.hpo_ids)
        ):
            raise ValueError(
                "`dataset.relationship_mask` must be indexed and columned by "
                "the same HPO IDs, in the same order, as "
                "`dataset.phenotype_matrix`."
            )

        self.dataset= dataset
        self.hpo_matrix = dataset.phenotype_matrix
        self.hpo_ids = dataset.hpo_ids
        self.hpo_labels = self.dataset.hpo_labels
        self.individual_ids = self.dataset.individual_ids
        self.relationship_mask = dataset.relationship_mask.to_numpy(
            copy=True
        )

        self.n_hpo_terms = len(self.hpo_ids)

    @staticmethod
    def _calculate_stats(
        count_00: int,
        count_01: int,
        count_10: int,
        count_11: int,
        ) -> tuple[float, float]:
        """
        Compute the Phi coefficient and Fisher's exact test p-value.

        The contingency table is defined as::

                  B=0   B=1
            A=0   00    01
            A=1   10    11

        where ``0`` denotes excluded and ``1`` denotes observed.

        Parameters
        ----------
        count_00 : int
            Number of individuals with both HPO terms excluded.

        count_01 : int
            Number of individuals with HPO A excluded and HPO B observed.

        count_10 : int
            Number of individuals with HPO A observed and HPO B excluded.

        count_11 : int
            Number of individuals with both HPO terms observed.

        Returns
        -------
        tuple[float, float]
            Phi coefficient and two-sided Fisher's exact test p-value.
            If either phenotype has no variation, both values are ``NaN``.
        """
        denominator = np.sqrt(
            (count_00 + count_01)
            * (count_10 + count_11)
            * (count_00 + count_10)
            * (count_01 + count_11)
        )

        if denominator == 0:
            return np.nan, np.nan

        numerator = (
            count_00 * count_11
            - count_01 * count_10
        )

        phi = numerator / denominator

        contingency_table = np.array(
            [
                [count_00, count_01],
                [count_10, count_11],
            ],
            dtype=int,
        )

        p_value = scipy.stats.fisher_exact(
            contingency_table,
            alternative="two-sided",
        ).pvalue

        return float(phi), float(p_value)

    def _calculate_pairwise_correlation(
        self,
        col_a: int,
        col_b: int, 
        pmids_by_individual: pd.Series | None = None,
    ) -> tuple[int, int, float, float, dict[str,Any]]:
        """
        Calculate the association between two HPO terms.

        Only individuals with non-missing values for both terms contribute
        to the contingency table.

        Parameters
        ----------
        col_a : int
            Column index of HPO term A.

        col_b : int
            Column index of HPO term B.

        pmids_by_individual : pd.Series | None, optional
            PubMed IDs indexed by individual ID. If ``None``, PubMed IDs
            are not included.

        Returns
        -------
        tuple[int, int, float, float, dict[str, Any]]
            Column indices, Phi coefficient, Fisher's exact test p-value,
            and pairwise count information.
        """
        matrix = self.hpo_matrix.values
        valid_mask = (~np.isnan(matrix[:, col_a])) & (~np.isnan(matrix[:, col_b]))
        col_a_values = matrix[valid_mask, col_a]
        col_b_values = matrix[valid_mask, col_b]

        count_11 = np.sum((col_a_values == 1) & (col_b_values == 1))
        count_10 = np.sum((col_a_values == 1) & (col_b_values == 0))
        count_01 = np.sum((col_a_values == 0) & (col_b_values == 1))
        count_00 = np.sum((col_a_values == 0) & (col_b_values == 0))
        n_individuals = len(col_a_values)

        counts: dict[str, Any] = {
            "00": count_00,
            "01": count_01,
            "10": count_10,
            "11": count_11,
            "N": n_individuals,
            "n_pmids": np.nan,
            "pmids": [],
        }

        # Phi is undefined if either variable has no variation
        # within the pair-specific effective sample.
        if (
            np.unique(col_a_values).size < 2
            or np.unique(col_b_values).size < 2
        ):
            return (
                col_a,
                col_b,
                np.nan,
                np.nan,
                counts,
            )
        
        try:
            phi, p_value = self._calculate_stats(
                count_00=count_00,
                count_01=count_01,
                count_10=count_10,
                count_11=count_11,
            )

            if pmids_by_individual is not None:
                individual_ids = self.individual_ids[valid_mask]

                pmids_list = (
                    pmids_by_individual
                    .loc[individual_ids]
                    .to_numpy()
                )

                all_pmids = sorted(
                    {
                        str(pmid)
                        for pmids in pmids_list
                        if pmids is not None
                        for pmid in pmids
                        if pd.notna(pmid)
                    }
                )

                counts["n_pmids"] = len(all_pmids)
                counts["pmids"] = all_pmids

            return (
                col_a,
                col_b,
                phi,
                p_value,
                counts,
            )

        except Exception as exc:
            logger.error(
                "Error calculating correlation for columns %d and %d: %s",
                col_a,
                col_b,
                exc,
            )
            return col_a, col_b, np.nan, np.nan, counts
        
    def compute_correlation_matrix(
        self,  
        min_individuals: int = 20,
        n_jobs: int = -1,
        include_pmids: bool = True
    ) -> CorrelationResult:
        """
        Calculate pairwise associations between HPO terms.

        Self-pairs and ancestor-descendant pairs are excluded using the
        ontology relationship mask. A remaining pair is tested only if at
        least ``min_individuals`` individuals have non-missing values for
        both HPO terms.

        Parameters
        ----------
        min_individuals : int, default=20
            Minimum number of individuals with non-missing values for both
            HPO terms required to test a pairwise association. This is a
            pair-specific effective sample size threshold, not a minimum
            cohort size.

        n_jobs : int, default=-1
            Number of parallel jobs. ``-1`` uses all available CPUs.

        include_pmids : bool, default=True
            If ``True``, aggregate unique PubMed IDs from individuals
            contributing to each tested pair.

        Returns
        -------
        CorrelationResult
            Pairwise association results together with symmetric Phi
            coefficient and raw p-value matrices.

        Raises
        ------
        TypeError
            If ``min_individuals`` is not an integer.

        ValueError
            If ``min_individuals`` is less than 1.
        """
        if not isinstance(min_individuals, int):
            raise TypeError(
                "`min_individuals` must be an integer."
            )

        if min_individuals < 1:
            raise ValueError(
                "`min_individuals` must be at least 1."
            )
        matrix = self.hpo_matrix.to_numpy()

        non_missing = ~np.isnan(matrix)

        pairwise_counts = (
            non_missing.T.astype(int)
            @ non_missing.astype(int)
        )

        upper_triangle = triu(
            coo_matrix(pairwise_counts),
            k=1,
        )

        row_indices = upper_triangle.row
        col_indices = upper_triangle.col
        effective_sample_sizes = upper_triangle.data

        ontology_values = self.relationship_mask[
            row_indices,
            col_indices,
        ]

        ontology_eligible = ~np.isnan(
            ontology_values
        )

        n_pairs_after_ontology_filter = int(
            np.sum(ontology_eligible)
        )

        eligible_indices = np.where(
            ontology_eligible
            & (effective_sample_sizes >= min_individuals)
        )[0]

        candidate_rows = row_indices[eligible_indices]
        candidate_cols = col_indices[eligible_indices]

        pairs = list(
            zip(
                candidate_rows,
                candidate_cols,
            )
        )

        if not pairs:
            logger.warning(
                "[Pairwise association: no eligible pairs]\n"
                "Pairs after ontology filtering: %d\n"
                "Pairs meeting effective N >= %d: 0\n"
                "No pairs were tested.",
                n_pairs_after_ontology_filter,
                min_individuals,
            )
            empty_results = pd.DataFrame(
                columns=[
                    "HPO_A",
                    "HPO_B",
                    "correlation",
                    "p_value",
                    "adj_p_value",
                ]
            )

            empty_matrix = pd.DataFrame(
                index=self.hpo_ids,
                columns=self.hpo_ids,
                dtype=float,
            )

            return CorrelationResult(
                correlation_results=empty_results,
                coef_matrix=empty_matrix,
                pval_matrix=empty_matrix.copy(),
                hpo_labels=self.hpo_labels,
            )

        pmids_by_individual = (
            self.dataset.get_pmids()
            if include_pmids
            else None
        )
        
        results = Parallel(n_jobs=n_jobs)(
            delayed(
                self._calculate_pairwise_correlation
            )(
                col_a,
                col_b,
                pmids_by_individual=pmids_by_individual,
            )
            for col_a, col_b in tqdm(
                pairs,
                desc="Calculating pairwise correlations",
            )
        )
        
        coef_matrix = np.full((self.n_hpo_terms, self.n_hpo_terms), np.nan)
        pvalue_matrix = np.full((self.n_hpo_terms, self.n_hpo_terms), np.nan)

        result_rows: list[dict[str, Any]] = []
        for (col_a, col_b, phi, p_value, counts) in results:
            if np.isnan(phi):
                continue

            coef_matrix[col_a, col_b] = phi
            coef_matrix[col_b, col_a] = phi

            pvalue_matrix[col_a, col_b] = p_value
            pvalue_matrix[col_b, col_a] = p_value

            hpo_a = self.hpo_ids[col_a]
            hpo_b = self.hpo_ids[col_b]

            row_data: dict[str, Any] = {
                "HPO_A": hpo_a,
                "HPO_B": hpo_b,
                "correlation": phi,
                "p_value": p_value,
                "n(A:E/B:E)": counts["00"],
                "n(A:E/B:O)": counts["01"],
                "n(A:O/B:E)": counts["10"],
                "n(A:O/B:O)": counts["11"],
                "n_individuals": counts["N"],
            }

            hpo_a_label = self.hpo_labels.get(hpo_a)
            hpo_b_label = self.hpo_labels.get(hpo_b)

            if hpo_a_label:
                row_data["HPO_A_label"] = hpo_a_label

            if hpo_b_label:
                row_data["HPO_B_label"] = hpo_b_label

            if include_pmids:
                row_data["n_pmids"] = counts["n_pmids"]
                row_data["pmids"] = ";".join(
                    counts["pmids"]
                )

            result_rows.append(row_data) 

        if not result_rows:
            logger.warning(
                "[Pairwise association: no valid results]\n"
                "Pairs evaluated: %d\n"
                "Pairs producing a valid Phi coefficient: 0\n"
                "For every evaluated pair, at least one HPO term had "
                "no variation within the pair-specific effective sample.",
                len(pairs),
            )
            empty_results = pd.DataFrame(
                columns=[
                    "HPO_A",
                    "HPO_B",
                    "correlation",
                    "p_value",
                    "adj_p_value",
                ]
            )

            empty_matrix = pd.DataFrame(
                index=self.hpo_ids,
                columns=self.hpo_ids,
                dtype=float,
            )

            return CorrelationResult(
                correlation_results=empty_results,
                coef_matrix=empty_matrix,
                pval_matrix=empty_matrix.copy(),
                hpo_labels=self.hpo_labels,
            )
        
        correlation_results = pd.DataFrame(
            result_rows
        )

        p_values = correlation_results[
            "p_value"
        ].to_numpy()

        _, adjusted_p_values, _, _ = multipletests(
            p_values,
            method="fdr_bh",
        )

        p_value_column = int(
            correlation_results.columns.get_loc(
                "p_value"
            )
        )

        correlation_results.insert(
            p_value_column + 1,
            "adj_p_value",
            adjusted_p_values,
        )

        correlation_results.sort_values(
            by="adj_p_value",
            ascending=True,
            inplace=True,
        )

        correlation_results.reset_index(
            drop=True,
            inplace=True,
        )

        valid_hpo_mask = ~np.isnan(
            coef_matrix
        ).all(axis=0)

        valid_hpo_ids = self.hpo_ids[
            valid_hpo_mask
        ]

        coef_df = pd.DataFrame(
            coef_matrix[
                np.ix_(
                    valid_hpo_mask,
                    valid_hpo_mask,
                )
            ],
            index=valid_hpo_ids,
            columns=valid_hpo_ids,
        )

        pval_df = pd.DataFrame(
            pvalue_matrix[
                np.ix_(
                    valid_hpo_mask,
                    valid_hpo_mask,
                )
            ],
            index=valid_hpo_ids,
            columns=valid_hpo_ids,
        )

        return CorrelationResult(
            correlation_results=correlation_results,
            coef_matrix=coef_df,
            pval_matrix=pval_df,
            hpo_labels=self.hpo_labels,
        )
