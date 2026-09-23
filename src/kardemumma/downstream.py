"""
Downstream analysis: differential expression testing and pathway/gene-set
enrichment analysis on protein/peptide quantification data.

:class:`DownStream` follows the same SDRF-driven convention as
:mod:`kardemumma.prm` plotting functions (e.g. ``plot_peptide_concentration_by_group``):
it takes a peptide/protein-level quantification DataFrame plus the SDRF table, and
compares two groups of an experimental variable column — by default the last
``factor value[...]`` column in the SDRF, since that is where SDRF convention
places the primary comparison variable. Pathway enrichment tests a list of IDs
of interest against user-supplied pathway/gene-set definitions using a
hypergeometric (Fisher's exact) test.
"""

import logging
from typing import Dict, Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .prm import _build_concentration_plot_df, _resolve_group_col

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Differential Expression
# ---------------------------------------------------------------------------


class DownStream:
    """
    Differential expression analysis between two groups of an SDRF experimental variable.

    Takes the wide-format output of :func:`kardemumma.prm.get_absolute_conc` (MultiIndex
    rows of ``qRePS`` / ``Peptide Sequence`` / ``Protein Name``, one column per replicate)
    plus the SDRF table, resolves the comparison group (default: the last
    ``factor value[...]`` column, matching SDRF's primary-variable convention -- same
    convention as :func:`kardemumma.prm.plot_peptide_concentration_by_group`), and runs a
    per-feature Welch's t-test via :meth:`differential_expression`, with :meth:`plot_volcano`
    for visualization -- the same way :class:`kardemumma.importer.MergeFiles` wraps
    skyline/SDRF merging.
    """

    def __init__(
        self,
        abs_df: pd.DataFrame,
        sdrf_data_file: pd.DataFrame,
        group_a: str,
        group_b: str,
        id_col: str = "Peptide Sequence",
        group_col: Optional[str] = None,
        log2_input: bool = False,
        equal_var: bool = False,
    ):
        """
        Args:
            abs_df: Wide output of :func:`kardemumma.prm.get_absolute_conc`.
            sdrf_data_file: SDRF table with ``'source name'`` and *group_col*.
            group_a: Group label (value in *group_col*) treated as the numerator (test) group.
            group_b: Group label treated as the denominator (reference) group.
            id_col: Feature column to test -- ``'Peptide Sequence'`` or ``'Protein Name'``.
            group_col: SDRF column to compare groups on. Defaults to the last
                ``'factor value[...]'`` column, matching SDRF's primary-variable convention.
            log2_input: If True, concentrations are already log2-scaled; log2fc is a mean difference.
                If False (default), values are log2-transformed internally (non-positive
                values are dropped with a warning before transforming).
            equal_var: Passed to scipy.stats.ttest_ind; False runs Welch's t-test.
        """
        self.abs_df = abs_df
        self.sdrf_data_file = sdrf_data_file
        self.group_a = group_a
        self.group_b = group_b
        self.id_col = id_col
        self.group_col = _resolve_group_col(sdrf_data_file, group_col)
        self.log2_input = log2_input
        self.equal_var = equal_var
        self._validate_groups()
        self.results = self.differential_expression()

    def _validate_groups(self) -> None:
        """
        Ensure ``group_col`` exists and has >=2 groups, and that ``group_a``/``group_b``
        are valid members of it -- ``group_col`` (e.g. the last SDRF ``factor value[...]``
        column) commonly has more than two groups (e.g. multiple sample-prep conditions
        plus a pool), so group_a/group_b must be given explicitly to pick which two to compare.
        """
        if "source name" not in self.sdrf_data_file.columns:
            raise KeyError("Column 'source name' not found in sdrf_data_file.")
        if self.group_col not in self.sdrf_data_file.columns:
            raise KeyError(f"Column '{self.group_col}' not found in sdrf_data_file.")

        available_groups = sorted(str(g) for g in self.sdrf_data_file[self.group_col].dropna().unique())
        if len(available_groups) < 2:
            msg = (
                f"SDRF column '{self.group_col}' has only {len(available_groups)} group(s) "
                f"({available_groups}) -- differential expression needs at least 2 groups to "
                "compare. Check the SDRF factor value assignments, or pass a different group_col."
            )
            logger.warning(msg)
            raise ValueError(msg)

        missing = [g for g in (self.group_a, self.group_b) if str(g) not in available_groups]
        if missing:
            raise ValueError(
                f"group_a/group_b must be among the {len(available_groups)} groups present in "
                f"'{self.group_col}': {available_groups}. Not found: {missing}."
            )

        print(
            f"Differential expression using SDRF column '{self.group_col}' "
            f"({len(available_groups)} group(s): {available_groups}) -- "
            f"comparing '{self.group_a}' vs '{self.group_b}'."
        )

    def differential_expression(self) -> pd.DataFrame:
        """
        Run the per-feature Welch's t-test comparing ``group_a`` vs ``group_b``.

        Returns:
            DataFrame with columns [id_col, 'mean_a', 'mean_b', 'log2fc', 'stat', 'p_value', 'padj'],
            sorted by 'padj' ascending. ``result.attrs`` carries 'group_col', 'group_a', 'group_b'.
        """
        value_col = "Protein conc [pmol]"
        long_df = _build_concentration_plot_df(self.abs_df, self.sdrf_data_file, self.group_col, self.group_col)

        annotate_with_protein = self.id_col != "Protein Name" and "Protein Name" in long_df.columns
        select_cols = [self.id_col, self.group_col, value_col] + (["Protein Name"] if annotate_with_protein else [])
        sub = long_df[long_df[self.group_col].isin([self.group_a, self.group_b])][select_cols].copy()
        sub[value_col] = pd.to_numeric(sub[value_col], errors="coerce")
        sub = sub.dropna(subset=[value_col])
        if sub.empty:
            raise ValueError("No rows left after filtering to group_a/group_b and dropping NaNs.")

        log_col = value_col
        if not self.log2_input:
            non_positive = int((sub[value_col] <= 0).sum())
            if non_positive:
                logger.warning(
                    "Dropping %d rows with non-positive %s prior to log2 transform.", non_positive, value_col
                )
                sub = sub[sub[value_col] > 0].copy()
            if sub.empty:
                raise ValueError("No positive values left after filtering; cannot log2-transform.")
            log_col = "_log2_value"
            sub[log_col] = np.log2(sub[value_col])

        results = []
        for feature_id, feat_df in sub.groupby(self.id_col):
            a = feat_df.loc[feat_df[self.group_col] == self.group_a, log_col].values
            b = feat_df.loc[feat_df[self.group_col] == self.group_b, log_col].values
            if len(a) < 2 or len(b) < 2:
                continue
            stat, p_value = stats.ttest_ind(a, b, equal_var=self.equal_var, nan_policy="omit")
            row = {
                self.id_col: feature_id,
                "mean_a": np.mean(a),
                "mean_b": np.mean(b),
                "log2fc": np.mean(a) - np.mean(b),
                "stat": stat,
                "p_value": p_value,
            }
            if annotate_with_protein:
                row["Protein Name"] = feat_df["Protein Name"].iloc[0]
            results.append(row)

        if not results:
            raise ValueError("No features had >=2 replicates in both groups; cannot test.")

        result_df = pd.DataFrame(results)
        result_df["padj"] = multipletests(result_df["p_value"], method="fdr_bh")[1]
        result_df = result_df.sort_values("padj").reset_index(drop=True)
        result_df.attrs.update(group_col=self.group_col, group_a=self.group_a, group_b=self.group_b)
        return result_df

    def plot_volcano(
        self,
        log2fc_col: str = "log2fc",
        pval_col: str = "padj",
        fc_threshold: float = 1.0,
        p_threshold: float = 0.05,
        label_top_n: int = 10,
    ):
        """
        Volcano plot of ``self.results``: log2 fold-change vs. -log10(p-value),
        highlighting significant hits.

        Args:
            log2fc_col: Column with log2 fold-change.
            pval_col: Column with (adjusted) p-value.
            fc_threshold: Absolute log2fc cutoff for significance.
            p_threshold: p-value cutoff for significance.
            label_top_n: Number of most-significant hits to label by name.

        Returns:
            matplotlib.figure.Figure
        """
        de_df = self.results
        for col in (self.id_col, log2fc_col, pval_col):
            if col not in de_df.columns:
                raise KeyError(f"Column '{col}' not found in DataFrame.")

        df = de_df.copy()
        df["neg_log10_p"] = -np.log10(df[pval_col].clip(lower=1e-300))
        sig = (
            (df[pval_col] < p_threshold) &
            (df[log2fc_col].abs() >= fc_threshold)
        )

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(
            df.loc[~sig, log2fc_col],
            df.loc[~sig, "neg_log10_p"],
            color="grey", alpha=0.5, s=20, label="Not significant",
        )
        ax.scatter(
            df.loc[sig, log2fc_col],
            df.loc[sig, "neg_log10_p"],
            color="crimson", alpha=0.7, s=20, label="Significant",
        )
        ax.axvline(
            fc_threshold, color="black", linestyle="--", linewidth=0.8
        )
        ax.axvline(
            -fc_threshold, color="black", linestyle="--", linewidth=0.8
        )
        ax.axhline(
            -np.log10(p_threshold), color="black",
            linestyle="--", linewidth=0.8
        )

        top_hits = df.loc[sig].sort_values(pval_col).head(label_top_n)
        # Use angles and offsets to reduce overlap
        angle_step = 360 / max(1, label_top_n)
        radius = 50
        for i, (_, row) in enumerate(top_hits.iterrows()):
            label = str(row[self.id_col])
            if "Protein Name" in df.columns:
                label = f"{label}\n{row['Protein Name']}"
            theta = np.deg2rad(i * angle_step)
            offset_x = int(np.cos(theta) * radius)
            offset_y = int(np.sin(theta) * radius) + 25
            # Stagger y offset a little further for even clearer separation
            offset_y += (i % 3) * 10
            ax.annotate(
                label,
                xy=(row[log2fc_col], row["neg_log10_p"]),
                xycoords='data',
                xytext=(offset_x, offset_y),
                textcoords='offset points',
                arrowprops=dict(arrowstyle="->", lw=0.7),
                fontsize=7,
                bbox=dict(
                    boxstyle="round,pad=0.2", fc="white",
                    alpha=0.7, lw=0
                ),
                ha='center'
            )

        ax.set_xlabel("log2 fold change")
        ax.set_ylabel(f"-log10({pval_col})")
        ax.set_title(
            f"Volcano plot: {self.group_a} vs "
            f"{self.group_b} ({self.group_col})"
        )
        ax.legend()
        fig.tight_layout()
        plt.close(fig)
        return fig



# ---------------------------------------------------------------------------
# Pathway / Gene-set Enrichment
# ---------------------------------------------------------------------------


def enrichment_analysis(
    hit_ids: Iterable[str],
    pathways: Dict[str, Sequence[str]],
    background_ids: Iterable[str],
    min_overlap: int = 2,
) -> pd.DataFrame:
    """
    Hypergeometric (Fisher's exact) enrichment test of hit_ids against pathway gene sets.

    Args:
        hit_ids: IDs of interest (e.g. significant proteins from DownStream.differential_expression).
        pathways: Mapping of pathway name -> iterable of member IDs.
        background_ids: Universe of IDs the test is drawn from (e.g. all quantified proteins).
        min_overlap: Minimum overlap between hit_ids and a pathway required to test it.

    Returns:
        DataFrame with columns ['pathway', 'overlap', 'pathway_size', 'hits_size',
        'background_size', 'p_value', 'padj'], sorted by 'padj' ascending.
    """
    hit_set = set(hit_ids)
    background_set = set(background_ids)
    if not hit_set.issubset(background_set):
        missing = hit_set - background_set
        raise ValueError(f"{len(missing)} hit_ids not found in background_ids, e.g. {list(missing)[:5]}")

    n_background = len(background_set)
    n_hits = len(hit_set)

    results = []
    for pathway_name, members in pathways.items():
        pathway_set = set(members) & background_set
        overlap = hit_set & pathway_set
        if len(overlap) < min_overlap or not pathway_set:
            continue
        table = [
            [len(overlap), len(pathway_set) - len(overlap)],
            [n_hits - len(overlap), n_background - len(pathway_set) - n_hits + len(overlap)],
        ]
        _, p_value = stats.fisher_exact(table, alternative="greater")
        results.append({
            "pathway": pathway_name,
            "overlap": len(overlap),
            "pathway_size": len(pathway_set),
            "hits_size": n_hits,
            "background_size": n_background,
            "p_value": p_value,
        })

    if not results:
        raise ValueError(f"No pathways met min_overlap={min_overlap} against hit_ids.")

    result_df = pd.DataFrame(results)
    result_df["padj"] = multipletests(result_df["p_value"], method="fdr_bh")[1]
    return result_df.sort_values("padj").reset_index(drop=True)


def plot_enrichment(
    enrichment_df: pd.DataFrame,
    top_n: int = 15,
    pval_col: str = "padj",
) -> None:
    """
    Horizontal bar plot of -log10(p-value) for the top enriched pathways.

    Args:
        enrichment_df: Output of :func:`enrichment_analysis`.
        top_n: Number of top pathways to display.
        pval_col: Column with (adjusted) p-value to rank/plot by.
    """
    if pval_col not in enrichment_df.columns or "pathway" not in enrichment_df.columns:
        raise KeyError(f"Expected 'pathway' and '{pval_col}' columns. Found: {list(enrichment_df.columns)}")

    plot_df = enrichment_df.sort_values(pval_col).head(top_n).copy()
    plot_df["neg_log10_p"] = -np.log10(plot_df[pval_col].clip(lower=1e-300))
    plot_df = plot_df.sort_values("neg_log10_p")

    plt.figure(figsize=(8, max(4, 0.4 * len(plot_df))))
    plt.barh(plot_df["pathway"], plot_df["neg_log10_p"], color="#7dc0a6")
    plt.xlabel(f"-log10({pval_col})")
    plt.title("Pathway enrichment")
    plt.tight_layout()
    plt.show()
