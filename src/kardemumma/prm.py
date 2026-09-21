"""
PRM (Parallel Reaction Monitoring) analysis and plotting functions for Skyline exports.

This module provides functions to assess quantification quality from Skyline
PRM reports. It expects a DataFrame produced by ``ImportFile.import_skyline_file``
and works with the following key columns:

- ``Precursor``               – precursor ion string (sequence + charge) 
Protein/Peptide/Precursors/Precursor
- ``Replicate``               – replicate / sample label
Replicates/Replicate
- ``File Name``               – raw data file name
- ``Peptide Sequence``        – stripped peptide sequence
- ``Peptide``                 – modified peptide string
Peptides/Peptide
- ``Protein Name``            – protein identifier
- ``Normalized Area``         – normalised peak area
- ``RatioLightToHeavy``       – light-to-heavy ratio (AQUA / SIS quantification)
- ``Ratio Dot Product``       – spectral dot product between light and heavy channels
- ``Library Dot Product``     – spectral dot product against the reference library
- ``Peptide Retention Time``  – observed retention time (minutes)
- ``Predicted Retention Time``– iRT-predicted retention time (minutes)
"""

import logging
from collections import Counter
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from sklearn.cluster import HDBSCAN
from statsmodels.stats.anova import anova_lm

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _to_python_scalar(x):
    """
    Reduce nested array/Series/list values to a plain Python scalar for patsy ``C()``.

    Patsy raises if any category cell is array-like with ndim > 1 or is non-scalar.
    """
    if x is None:
        return np.nan
    if isinstance(x, str):
        return x
    if isinstance(x, bytes):
        return x.decode()
    if isinstance(x, np.str_):
        return str(x)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return np.nan if pd.isna(x) else float(x)
    if isinstance(x, pd.Series):
        return np.nan if x.empty else _to_python_scalar(x.iloc[0])
    if isinstance(x, np.ndarray):
        if x.size == 0:
            return np.nan
        if x.ndim > 1:
            return _to_python_scalar(x.ravel()[0])
        if x.dtype == object:
            return _to_python_scalar(x.item() if x.shape == () else x.flat[0])
        out = x.flat[0]
        return _to_python_scalar(out) if isinstance(out, np.ndarray) else out.item()
    if isinstance(x, (list, tuple)):
        return np.nan if len(x) == 0 else _to_python_scalar(x[0])
    return x


def _as_1d_series(df: pd.DataFrame, col: str) -> pd.Series:
    """Return ``df[col]`` as a single Series (first column if duplicated)."""
    if col not in df.columns:
        raise KeyError(f"Missing column {col!r}. Found: {list(df.columns)}")
    obj = df[col]
    return obj.iloc[:, 0].copy() if isinstance(obj, pd.DataFrame) else obj.copy()


def _patsy_scalar_categorical(series: pd.Series) -> pd.Series:
    return series.map(_to_python_scalar)


def _resolve_group_col(sdrf_data_file: pd.DataFrame, group_col: Optional[str]) -> str:
    """Return group_col if provided; otherwise use the last column starting with 'factor value'."""
    if group_col is not None:
        return group_col
    factor_cols = [c for c in sdrf_data_file.columns if c.startswith("factor value")]
    if not factor_cols:
        raise KeyError("No column starting with 'factor value' found in sdrf_data_file.")
    return factor_cols[-1]


def _resolve_sdrf_file_key(sdrf_data_file: pd.DataFrame) -> str:
    """
    SDRF column that matches :func:`get_absolute_conc` wide-format column names.

    ``get_absolute_conc`` pivots on Skyline ``File Name`` (typically ``*.raw``), which
    aligns with ``comment[data file]`` in the SDRF, not ``source name``.
    """
    data_file_col = "comment[data file]"
    if data_file_col in sdrf_data_file.columns:
        return data_file_col
    if "source name" not in sdrf_data_file.columns:
        raise KeyError(
            "Expected 'comment[data file]' or 'source name' in sdrf_data_file."
        )
    return "source name"


_DISEASE_CATEGORY_PALETTE: dict[str, str] = {
    "Healthy": "#c9b28f",
    "Cardiovascular": "#ed936b",
    "Metabolic": "#e0c59a",
    "Cancer": "#919fc7",
    "Psychiatric": "#7dc0a6",
    "Neurologic": "#7dc0a6",
    "Autoimmune": "#da8ec0",
    "Infection": "#f9da56",
}


def _build_concentration_plot_df(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    group_col: str,
    color_col: str,
) -> pd.DataFrame:
    long_df = (
        abs_df.reset_index()
        .melt(
            id_vars=["qRePS", "Peptide Sequence", "Protein Name"],
            var_name="Replicate",
            value_name="Protein conc [pmol]",
        )
        .dropna(subset=["Protein conc [pmol]"])
    )
    merge_cols = list(dict.fromkeys(["source name", group_col, color_col]))
    return long_df.merge(
        sdrf_data_file[merge_cols].drop_duplicates(),
        left_on="Replicate", right_on="source name", how="left",
    )


def _make_group2color(color_groups, color_col: str) -> dict:
    if color_col == "characteristics[disease category]":
        return {g: _DISEASE_CATEGORY_PALETTE.get(g, "#cccccc") for g in color_groups}
    return dict(zip(color_groups, sns.color_palette(n_colors=len(color_groups))))


def _plot_protein_boxes(
    protein_df: pd.DataFrame,
    protein_name: str,
    group_col: str,
    color_col: str,
    group2color: dict,
) -> None:
    peptides = protein_df["Peptide Sequence"].unique()
    if len(peptides) == 0:
        return

    group_to_color_label = (
        protein_df.drop_duplicates(subset=[group_col])
        .set_index(group_col)[color_col]
    )
    if color_col == "characteristics[disease category]":
        cat_rank = {c: i for i, c in enumerate(_DISEASE_CATEGORY_PALETTE)}
        key_fn = lambda g: (cat_rank.get(group_to_color_label.get(g, ""), len(cat_rank)), str(g))
    else:
        key_fn = lambda g: (str(group_to_color_label.get(g, "")), str(g))
    group_order = sorted(protein_df[group_col].dropna().unique(), key=key_fn)

    fig, axes = plt.subplots(len(peptides), 1, figsize=(10, 3 * len(peptides)), sharex=True)
    if len(peptides) == 1:
        axes = [axes]

    for ax, pep in zip(axes, peptides):
        pep_df = protein_df[protein_df["Peptide Sequence"] == pep]
        sns.boxplot(
            data=pep_df, x=group_col, y="Protein conc [pmol]",
            hue=color_col, ax=ax, order=group_order, palette=group2color,
            dodge=False,
        )
        sns.stripplot(
            data=pep_df, x=group_col, y="Protein conc [pmol]",
            ax=ax, order=group_order, color="grey",
            dodge=False, jitter=True, alpha=0.5, size=3,
        )
        # for idx, group in enumerate(group_order):
        #     group_data = pep_df[pep_df[group_col] == group]["Protein conc [pmol]"]
        #     if not group_data.empty:
        #         mean_val = group_data.mean()
        #         ax.text(idx, mean_val, f"{mean_val:.2f}",
        #                 ha='center', va='center', fontsize=9, fontweight="bold", color='black',
        #                 bbox=dict(facecolor='white', edgecolor='none', pad=0.3, alpha=0.7))
        ax.set_title(f"{pep}|{protein_name}", fontsize=10, loc='left')
        ax.set_ylabel("Protein conc [pmol]", fontsize=5)
        ax.set_xlabel(group_col if ax == axes[-1] else "", fontsize=5)
        ax.tick_params(axis='x', rotation=90, labelsize=5)
        ax.tick_params(axis='y', labelsize=5)
        if ax.get_legend() is not None:
            ax.legend_.remove()

    handles, labels = axes[0].get_legend_handles_labels()
    if color_col == "characteristics[disease category]":
        cat_rank = {c: i for i, c in enumerate(_DISEASE_CATEGORY_PALETTE)}
        paired = sorted(
            zip(labels, handles),
            key=lambda x: cat_rank.get(x[0], len(cat_rank)),
        )
        labels, handles = zip(*paired) if paired else (labels, handles)
    fig.legend(
        handles, labels,
        loc="center right", bbox_to_anchor=(1.15, 0.5),
        fontsize=5, title=color_col, title_fontsize=5,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()


def _formula_clean_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Flatten MultiIndex columns, drop duplicates, reset index — so patsy sees a plain table."""
    out = df.copy()
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [
            "_".join(str(p) for p in tup if str(p) != "") if isinstance(tup, tuple) else tup
            for tup in out.columns
        ]
    out = out.loc[:, ~out.columns.duplicated(keep="first")]
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Quality Filtering
# ---------------------------------------------------------------------------


def filter_library_dot_product(
    df: pd.DataFrame,
    *,
    threshold: float = 0.8,
    col: str = "Library Dot Product",
) -> pd.DataFrame:
    """
    Keep rows where library dot product is **strictly greater** than *threshold*,
    then pivot to wide format (one column per isotope label type).

    Args:
        df: Skyline report DataFrame.
        threshold: Exclusive lower bound on Library Dot Product. Default ``0.8``.
        col: Column name for library dot product.

    Raises:
        KeyError: If *col* is missing from *df*.
    """
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame.")

    df = df.loc[df[col] > threshold].copy()
    df = df[df["Normalized Area"].notna()]
    
    # Print out Isotope Label Type count 
    print(df["Isotope Label Type"].value_counts())
    
    wide = df.pivot_table(
        index=["Replicate", "Protein Name", "Peptide"],
        columns="Isotope Label Type",
        values="Normalized Area",
        aggfunc="first",
    )
    wide = wide.reset_index()  # Moves "Isotope Label Type" from index to columns
    return wide


def filter_peptide_counts(
    peptide_counts_df: pd.DataFrame,
    light_cutoff: int = 700,
    heavy_cutoff: int = 700,
) -> pd.DataFrame:
    """Keep only rows where both light and heavy counts exceed their respective cutoffs."""
    mask = (
        (peptide_counts_df["light_count"] > light_cutoff) &
        (peptide_counts_df["heavy_count"] > heavy_cutoff)
    )
    return peptide_counts_df.loc[mask].reset_index(drop=True)


def cluster_abundant_peptides(
    peptide_counts: pd.DataFrame,
    min_cluster_size: int = 10,
    min_samples: Optional[int] = None,
) -> Tuple[pd.DataFrame, dict]:
    """
    Cluster peptides by their limiting detection count —
    ``min(heavy_count, light_count)`` — with HDBSCAN, and flag the
    most-abundant cluster as a data-driven alternative to a fixed count
    cutoff (see :func:`filter_peptide_counts`).

    A peptide's abundance is bottlenecked by whichever channel is
    worse-detected, so clustering on that single score (rather than raw
    ``(heavy_count, light_count)`` coordinates) correctly groups the whole
    "both channels well-detected" band together even when one channel's
    count is much more spread out than the other — a plain 2D HDBSCAN over
    the raw coordinates tends to fracture that band into small noise
    fragments instead of one cluster.

    The non-noise cluster with the highest mean score is taken as the
    "abundant" group. Points HDBSCAN labels as noise (-1) are never
    considered abundant.

    Args:
        peptide_counts: DataFrame with ``heavy_count`` and ``light_count`` columns.
        min_cluster_size: HDBSCAN ``min_cluster_size`` parameter.
        min_samples: HDBSCAN ``min_samples`` parameter (defaults to ``min_cluster_size``).

    Returns:
        ``(clustered_df, cutoff)`` where *clustered_df* is *peptide_counts* with
        added ``cluster`` and ``is_abundant`` columns, and *cutoff* is
        ``{"heavy_count": ..., "light_count": ...}`` — the minimum counts observed
        within the abundant cluster.
    """
    for col in ("heavy_count", "light_count"):
        if col not in peptide_counts.columns:
            raise KeyError(f"Column '{col}' not found in peptide_counts.")

    score = peptide_counts[["heavy_count", "light_count"]].min(axis=1).to_numpy(dtype=float).reshape(-1, 1)
    labels = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples).fit_predict(score)

    result = peptide_counts.copy()
    result["cluster"] = labels
    result["_score"] = score.ravel()

    non_noise = result[result["cluster"] != -1]
    if non_noise.empty:
        raise ValueError(
            "HDBSCAN found no clusters (all points labeled noise); try a smaller min_cluster_size."
        )

    cluster_scores = non_noise.groupby("cluster")["_score"].mean()
    abundant_cluster = cluster_scores.idxmax()
    result["is_abundant"] = result["cluster"] == abundant_cluster

    # min(heavy, light) >= c  <=>  heavy >= c AND light >= c, so the same
    # scalar threshold applies to both channels — using separate per-axis
    # minima here would be inconsistent with how the cluster was formed.
    c = int(result.loc[result["is_abundant"], "_score"].min())
    result = result.drop(columns="_score")
    cutoff = {"heavy_count": c, "light_count": c}

    logger.info(
        "HDBSCAN found %d cluster(s) (%d noise points); abundant cluster=%s, n=%d, cutoff=%s",
        len(cluster_scores), int((result["cluster"] == -1).sum()),
        abundant_cluster, int(result["is_abundant"].sum()), cutoff,
    )

    return result, cutoff


# ---------------------------------------------------------------------------
# Peptide Detection Summary
# ---------------------------------------------------------------------------


def summarise_peptide_counts(peptide_counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Count non-missing heavy and light measurements per peptide.

    Expects columns ``Protein Name``, ``Peptide``, and case-insensitive
    ``Heavy`` / ``Light`` intensity columns.

    Returns:
        DataFrame with columns ``['Protein Name', 'Peptide', 'heavy_count', 'light_count']``.
    """
    possible_heavy = [c for c in peptide_counts_df.columns if c.lower() == "heavy"]
    possible_light = [c for c in peptide_counts_df.columns if c.lower() == "light"]
    if not possible_heavy or not possible_light:
        raise KeyError("DataFrame must contain 'Heavy' and 'Light' columns (case-insensitive).")

    return (
        peptide_counts_df.groupby(["Protein Name", "Peptide"])
        .agg(
            heavy_count=(possible_heavy[0], lambda x: x.notna().sum()),
            light_count=(possible_light[0], lambda x: x.notna().sum()),
        )
        .reset_index()
    )


def report_peptide_protein_summary(
    peptide_counts: pd.DataFrame,
) -> Tuple[dict, list]:
    """
    Report unique peptide and protein counts.

    Args:
        peptide_counts: Output of :func:`summarise_peptide_counts`.

    Returns:
        ``(summary_dict, peptide_list)`` where *summary_dict* contains
        ``num_unique_peptides`` and ``num_unique_proteins``.
    """
    expected = ["Protein Name", "Peptide", "heavy_count", "light_count"]
    if not all(c in peptide_counts.columns for c in expected):
        raise ValueError(f"Expected columns {expected} not found in peptide_counts.")

    n_peptides = peptide_counts["Peptide"].nunique()
    n_proteins = peptide_counts["Protein Name"].nunique()
    peptide_list = list(pd.unique(peptide_counts["Peptide"]))

    logger.info("Unique peptides: %d", n_peptides)
    logger.info("Unique proteins: %d", n_proteins)
    logger.info("Selected peptides: %s", peptide_list)

    return {"num_unique_peptides": n_peptides, "num_unique_proteins": n_proteins}, peptide_list


# ---------------------------------------------------------------------------
# CV Analysis
# ---------------------------------------------------------------------------


def calculate_intra_plate_cv(
    pool_data: pd.DataFrame,
    col_name: str = "characteristics[plate]",
) -> pd.DataFrame:
    """
    Calculate intra-plate CV per peptide.

    Args:
        pool_data: DataFrame with at least ``[col_name, 'Peptide Sequence', 'RatioLightToHeavy']``.
        col_name: Column to group by. Default ``'characteristics[plate]'``.

    Returns:
        DataFrame with ``mean``, ``std``, and ``intra_plate_cv`` per plate/peptide.
    """
    if col_name not in pool_data.columns:
        raise KeyError(f"Column '{col_name}' not found in pool_data.")

    pool_data = pool_data[[col_name, "Peptide Sequence", "RatioLightToHeavy"]].drop_duplicates()
    stats = (
        pool_data.groupby([col_name, "Peptide Sequence"])["RatioLightToHeavy"]
        .agg(["mean", "std"])
        .reset_index()
    )
    stats["intra_plate_cv"] = stats["std"] / stats["mean"]
    return stats


def calculate_inter_plate_cv(peptide_plate_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate inter-plate CV from intra-plate stats.

    Args:
        peptide_plate_stats: Output of :func:`calculate_intra_plate_cv`.

    Returns:
        DataFrame with ``['Peptide Sequence', 'grand_mean', 'between_plate_sd', 'inter_plate_cv']``,
        sorted by ``inter_plate_cv`` ascending.
    """
    peptide_means = (
        peptide_plate_stats.groupby("Peptide Sequence")["mean"]
        .agg(["mean", "std"])
        .rename(columns={"mean": "grand_mean", "std": "between_plate_sd"})
        .reset_index()
    )
    peptide_means["inter_plate_cv"] = peptide_means["between_plate_sd"] / peptide_means["grand_mean"]
    return peptide_means.sort_values("inter_plate_cv").reset_index(drop=True)


def get_lowest_cv_peptides(
    interplate_cv_df: pd.DataFrame,
    cv_percentile: float,
) -> np.ndarray:
    """
    Return peptide sequences with inter-plate CV at or below *cv_percentile*.

    Args:
        interplate_cv_df: Output of :func:`calculate_inter_plate_cv`.
        cv_percentile: Threshold as a percentage (e.g. ``10`` for 10%).
    """
    cv_threshold = cv_percentile / 100.0
    return interplate_cv_df.loc[
        interplate_cv_df["inter_plate_cv"] <= cv_threshold, "Peptide Sequence"
    ].unique()


def extract_top_percentile(
    df: pd.DataFrame,
    column: str,
    percentile: float = 0.1,
    id_col: str = "Peptide Sequence",
    source_df: Optional[pd.DataFrame] = None,
    source_col: Optional[str] = None,
) -> Tuple[np.ndarray, Optional[pd.DataFrame]]:
    """
    Return IDs whose *column* value is at or below *percentile*, and optionally
    filter a source DataFrame to those IDs.

    Args:
        df: DataFrame with peptide statistics (e.g. output of :func:`calculate_inter_plate_cv`).
        column: Column to compute the quantile threshold from.
        percentile: Quantile threshold as a fraction in [0, 1] (e.g. ``0.1`` = lowest 10%).
        id_col: Column whose unique values to extract. Default ``'Peptide Sequence'``.
        source_df: If given, filter this DataFrame to matching IDs.
        source_col: Column of *source_df* to match against (defaults to *id_col*).

    Returns:
        ``(id_array, filtered_df)``; *filtered_df* is ``None`` if *source_df* is not provided.
    """
    threshold = df[column].quantile(percentile)
    id_list = df.loc[df[column] <= threshold, id_col].unique()
    if source_df is not None:
        key = source_col if source_col is not None else id_col
        return id_list, source_df[source_df[key].isin(id_list)].reset_index(drop=True)
    return id_list, None


# ---------------------------------------------------------------------------
# Plate Normalization
# ---------------------------------------------------------------------------


def plate_peptide_anova(
    selected_norm_peptides: pd.DataFrame,
    plate_col: str = "characteristics[plate]",
    log_transform: bool = False,
) -> Tuple:
    """
    Fit a two-way ANOVA: ``ratio ~ C(Plate) + C(Peptide)``.

    Args:
        selected_norm_peptides: DataFrame with ``RatioLightToHeavy``, plate, and peptide columns.
        plate_col: Column to use as the plate factor. Normalised internally to ``"Plate"``.
        log_transform: If ``True``, model ``log(RatioLightToHeavy)`` (positive rows only).

    Returns:
        ``(model, anova_res)`` — statsmodels OLS result and Type II ANOVA table.
    """
    df = _formula_clean_frame(selected_norm_peptides.copy())

    # Normalize plate column
    if plate_col in df.columns and plate_col != "Plate":
        df = df.drop(columns=["Plate"], errors="ignore").rename(columns={plate_col: "Plate"})
    elif "characteristics[plate]" in df.columns:
        df = df.drop(columns=["Plate"], errors="ignore").rename(
            columns={"characteristics[plate]": "Plate"}
        )

    # Normalize peptide column
    if "Peptide Sequence" in df.columns:
        df = df.drop(columns=["Peptide"], errors="ignore").rename(
            columns={"Peptide Sequence": "Peptide"}
        )

    df = _formula_clean_frame(df)

    required_cols = ["RatioLightToHeavy", "Plate", "Peptide"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise KeyError(f"plate_peptide_anova missing columns {missing!r}. Have: {list(df.columns)!r}")

    # Parse columns and drop NA
    df = df.loc[:, required_cols].dropna().reset_index(drop=True)
    df["Plate"] = df["Plate"].astype(str)
    df["Peptide"] = df["Peptide"].astype(str)
    df["RatioLightToHeavy"] = pd.to_numeric(df["RatioLightToHeavy"], errors="coerce")

    if df.empty:
        raise ValueError("plate_peptide_anova: no rows left after dropping missing values.")

    # Response column
    if log_transform:
        df = df[df["RatioLightToHeavy"] > 0].copy().reset_index(drop=True)
        if df.empty:
            raise ValueError("plate_peptide_anova: no positive RatioLightToHeavy rows for log_transform.")
        df["response"] = np.log(df["RatioLightToHeavy"])
    else:
        df["response"] = df["RatioLightToHeavy"]

    model = smf.ols("response ~ C(Plate) + C(Peptide)", data=df).fit()
    anova_res = anova_lm(model, typ=2)
    logger.info("OLS model summary:\n%s", model.summary())
    logger.info("Type II ANOVA table:\n%s", anova_res)

    def _get_p_value(keyword: str):
        # Try to get the p-value column (works for both PR(>F), PR(>F), or "P" columns)
        label_rows = [idx for idx in anova_res.index if keyword in str(idx).lower()]
        pr_col = next((c for c in anova_res.columns if c.upper().startswith("P")), None)
        return None if (not label_rows or not pr_col) else anova_res.loc[label_rows[0], pr_col]

    for label in ["Plate", "Peptide"]:
        p_val = _get_p_value(label.lower())
        if p_val is not None and pd.notna(p_val):
            logger.info("%s effect p-value (Type II ANOVA): %.4g", label, p_val)
            if label == "Plate":
                logger.info(
                    "Plate effect is %s.",
                    "significant — batch correction recommended" if p_val < 0.05 else "not significant",
                )
        else:
            logger.warning("%s effect p-value could not be read from the ANOVA table.", label)

    return model, anova_res

def get_plate_conversion_factors(
    df: pd.DataFrame,
    col_plate: str = "Plate",
    ratio_col: str = "RatioLightToHeavy",
    log_transform: bool = False,
) -> Tuple[pd.DataFrame, dict, object]:
    """
    Calculate per-plate multiplicative correction factors.

    Returns:
        (platemed_df, conversion_factors, model): plate medians and factors,
        {plate: factor} dict, and statsmodels OLS fit object.
    """
    cols = ["Peptide Sequence", "Replicate", col_plate, ratio_col]
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing required columns: {missing_cols}")

    d = df[cols].drop_duplicates().copy()
    d["ratio_fit"] = pd.to_numeric(d[ratio_col], errors="coerce")
    if log_transform:
        d["ratio_fit"] = np.log(d["ratio_fit"])

    model = smf.ols(
        f'ratio_fit ~ C(Q("{col_plate}")) + C(Q("Peptide Sequence"))', data=d
    ).fit()

    platemed = (
        d.groupby(col_plate, dropna=True, observed=True)["ratio_fit"]
        .median()
        .reset_index()
        .rename(columns={"ratio_fit": "plate_median"})
    )
    global_median = d["ratio_fit"].median()
    platemed["correction_factor"] = global_median / platemed["plate_median"]
    conversion_factors = dict(zip(platemed[col_plate].astype(str), platemed["correction_factor"]))
    
    # Print the conversion factors
    print(f"Conversion factors: {conversion_factors}")
    return platemed, conversion_factors, model


def plot_plate_conversion_factors(
    df: pd.DataFrame,
    col_plate: str = "Plate",
    ratio_col: str = "RatioLightToHeavy",
    log_transform: bool = False,
):
    """
    Plot peptide sequence vs ratio_fit colored by plate.
    """
    cols = ["Peptide Sequence", "Replicate", col_plate, ratio_col]
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing required columns: {missing_cols}")

    d = df[cols].drop_duplicates().copy()
    d["ratio_fit"] = pd.to_numeric(d[ratio_col], errors="coerce")
    if log_transform:
        d["ratio_fit"] = np.log(d["ratio_fit"])

    model = smf.ols(
        f'ratio_fit ~ C(Q("{col_plate}")) + C(Q("Peptide Sequence"))', data=d
    ).fit()

    plot_df = d[["Peptide Sequence", col_plate, "ratio_fit"]]
    peptide_means = plot_df.groupby("Peptide Sequence", observed=True)["ratio_fit"].mean().sort_values()
    ordered_peptides = peptide_means.index.tolist()
    # Use categorical dtype for peptide order
    plot_df["Peptide Sequence"] = pd.Categorical(
        plot_df["Peptide Sequence"], categories=ordered_peptides, ordered=True
    )
    plt.figure(figsize=(10, 6))
    ax = sns.boxplot(
        x="Peptide Sequence",
        y="ratio_fit",
        hue=col_plate,
        data=plot_df,
        dodge=True,
        showfliers=False,
    )
    plt.title("Peptide Sequence vs Ratio Fit (by Plate)")
    plt.xlabel("Peptide Sequence")
    plt.ylabel("Ratio Fit")
    plt.xticks(rotation=90)
    plt.tight_layout()
    # Place the legend at the bottom, horizontally
    plt.legend(title=col_plate, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=len(plot_df[col_plate].unique()), frameon=False)
    plt.show()
    return plt.gcf()


def adjust_ratio_by_plate(
    df: pd.DataFrame,
    conversion_factors: dict,
    col_match: str = "characteristics[plate]",
    ignore_nan_plates: bool = True,
) -> pd.DataFrame:
    """
    Adjust the RatioLightToHeavy values using provided plate conversion factors.

    Args:
        df: Input DataFrame (must contain "RatioLightToHeavy" and col_match columns).
        conversion_factors: Dict mapping plate values (col_match) to correction factors.
        col_match: Column whose values match the keys in conversion_factors.

    Returns:
        DataFrame with updated "RatioLightToHeavy" values normalized per plate.
    """
    df = df.copy()
    if "RatioLightToHeavy" not in df.columns:
        raise KeyError("DataFrame missing required column: 'RatioLightToHeavy'")
    if col_match not in df.columns:
        raise KeyError(f"DataFrame missing required column: '{col_match}'")
    if ignore_nan_plates:
        df = df[df[col_match].notna()].copy()
    plate_factors = df[col_match].astype(str).map(conversion_factors)
    if plate_factors.isnull().any():
        missing = df.loc[plate_factors.isnull(), col_match].unique()
        raise KeyError(
            f"Some plate values have no conversion factor: {missing}"
        )
    df["RatioLightToHeavy"] = df["RatioLightToHeavy"] / plate_factors.values
    return df



# ---------------------------------------------------------------------------
# Combine qRePS table with skyline_merge_adj
# ---------------------------------------------------------------------------

def get_absolute_conc(
    qreps_table: pd.DataFrame,
    skyline_df: pd.DataFrame,
    skyline_protein_col: str = 'Protein Name',
    ratio_cutoff_low: float = 10**-3,
    ratio_cutoff_high: float = 10**3,
) -> pd.DataFrame:
    """
    Combine qRePS table with skyline_merge_adj, using File Name (not Replicate) as the identifier for samples/files.
    Ratios below or above cutoffs will be set to np.nan.
    """
    # Check for necessary columns in both dataframes
    if skyline_protein_col not in skyline_df.columns:
        raise KeyError(f"Column '{skyline_protein_col}' not found in skyline_df.")
    # Check if qRePS is in qreps_table
    if 'qRePS' not in qreps_table.columns:
        raise KeyError("Column 'qRePS' not found in qreps_table.")
    # Check for File Name in skyline_df
    if 'File Name' not in skyline_df.columns:
        raise KeyError("Column 'File Name' not found in skyline_df.")

    # Extract qRePS id from skyline_df from the skyline_protein_col column
    skyline_df = skyline_df.copy()
    skyline_df['qRePS'] = skyline_df[skyline_protein_col].str.extract(r'(QR\d+)')

    # Merge skyline_df and qreps_table on 'qRePS'
    combined_df = pd.merge(skyline_df, qreps_table, on='qRePS', how='left')

    # Check for required columns after merge
    for col in ["RatioLightToHeavy", "Amount per well [pmol]", "Peptide Sequence", "Protein Name", "File Name"]:
        if col not in combined_df.columns:
            raise KeyError(f"Column '{col}' not found in the merged DataFrame.")

    # Filter RatioLightToHeavy by cutoffs (set out of bounds to np.nan)
    ratio = combined_df["RatioLightToHeavy"]
    out_of_bounds = (ratio < ratio_cutoff_low) | (ratio > ratio_cutoff_high)
    combined_df.loc[out_of_bounds, "RatioLightToHeavy"] = np.nan

    # print out of bounds ratio
    print(f"Number of out of bounds ratios: {out_of_bounds.sum()}")
    print(f"Ratio cutoff low: {ratio_cutoff_low}")
    print(f"Ratio cutoff high: {ratio_cutoff_high}")
    # Print numbers of data ratio that has been removed
    print(f"Number of data ratio that has been removed: {combined_df.loc[out_of_bounds, 'RatioLightToHeavy'].count()}")

    # Calculate absolute protein concentration [pmol]
    combined_df["Protein conc [pmol]"] = (
        combined_df["RatioLightToHeavy"] * combined_df["Amount per well [pmol]"]
    ).round(4)

    # Select columns to export (replace Replicate with File Name)
    export_cols = ["qRePS", "Peptide Sequence", "Protein Name", "File Name", "Protein conc [pmol]"]
    combined_df = combined_df[export_cols]

    # Pivot to wide format (columns=File Name)
    combined_df_wide = combined_df.pivot_table(
        index=["qRePS", "Peptide Sequence", "Protein Name"],
        columns="File Name",
        values="Protein conc [pmol]",
        aggfunc="first"
    )

    if combined_df_wide.empty:
        raise ValueError("Combined dataframe is empty. Please check the input data.")

    def _report_abs(df: pd.DataFrame) -> None:
        idx = df.index
        # idx is a MultiIndex with levels: 'qRePS', 'Peptide Sequence', 'Protein Name'
        n_qreps = len(idx.get_level_values('qRePS').unique())
        n_peptides = len(idx.get_level_values('Peptide Sequence').unique())
        n_proteins = len(idx.get_level_values('Protein Name').unique())
        n_files = len(df.columns.unique())
        print(f"Number of unique qRePS ids: {n_qreps}")
        print(f"Number of unique Peptide Sequences: {n_peptides}")
        print(f"Number of unique Protein Names: {n_proteins}")
        print(f"Number of unique File Names: {n_files}")

    _report_abs(combined_df_wide)

    return combined_df_wide

def compute_cv(
    df: pd.DataFrame,
    value_col: str = "RatioLightToHeavy",
    group_by: str = "Precursor",
) -> pd.DataFrame:
    """
    Compute %CV per group.

    Args:
        df: Skyline report DataFrame.
        value_col: Numeric column to compute CV on. Default ``'RatioLightToHeavy'``.
        group_by: Grouping column. Default ``'Precursor'``.

    Returns:
        DataFrame with ``[group_by, 'mean', 'std', 'cv_pct', 'n']``, sorted by ``cv_pct``.
    """
    for col in (value_col, group_by):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    result = df.groupby(group_by)[value_col].agg(["mean", "std", "count"]).reset_index()
    result = result.rename(columns={"count": "n"})
    result["cv_pct"] = (result["std"] / result["mean"]).abs() * 100
    result.loc[result["mean"].abs() < 1e-12, "cv_pct"] = float("nan")
    return result[[group_by, "mean", "std", "cv_pct", "n"]].sort_values(
        "cv_pct", na_position="last"
    ).reset_index(drop=True)


def flag_missing_values(
    df: pd.DataFrame,
    value_col: str = "RatioLightToHeavy",
    precursor_col: str = "Precursor",
    replicate_col: str = "Replicate",
) -> pd.DataFrame:
    """
    Identify precursors with missing quantification values across replicates.

    Returns:
        DataFrame with one row per precursor and columns
        ``n_detected``, ``n_total``, ``detection_rate``.
    """
    for col in (value_col, precursor_col, replicate_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    pivot = df.pivot_table(
        index=precursor_col, columns=replicate_col, values=value_col, aggfunc="first"
    )
    detected = pivot.notna() & (pivot != 0)
    n_total = len(pivot.columns)
    summary = pd.DataFrame({
        precursor_col: pivot.index,
        "n_detected": detected.sum(axis=1).values,
        "n_total": n_total,
    })
    summary["detection_rate"] = summary["n_detected"] / n_total
    return summary.sort_values("detection_rate").reset_index(drop=True)


def dot_product_summary(
    df: pd.DataFrame,
    lib_dot_col: str = "Library Dot Product",
    ratio_dot_col: str = "Ratio Dot Product",
    lib_threshold: float = 0.8,
    ratio_threshold: float = 0.8,
) -> pd.DataFrame:
    """
    Flag precursor–replicate pairs that fall below dot-product quality thresholds.

    Returns:
        DataFrame with ``lib_dot_pass``, ``ratio_dot_pass``, and ``both_pass`` columns.
    """
    for col in (lib_dot_col, ratio_dot_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")
    missing_cols = [c for c in ("Precursor", "Replicate") if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Required columns missing: {missing_cols}")

    result = df[["Precursor", "Replicate", lib_dot_col, ratio_dot_col]].copy()
    result["lib_dot_pass"] = result[lib_dot_col] >= lib_threshold
    result["ratio_dot_pass"] = result[ratio_dot_col] >= ratio_threshold
    result["both_pass"] = result["lib_dot_pass"] & result["ratio_dot_pass"]
    return result.reset_index(drop=True)


def retention_time_deviation(
    df: pd.DataFrame,
    observed_col: str = "Peptide Retention Time",
    predicted_col: str = "Predicted Retention Time",
) -> pd.DataFrame:
    """
    Compute observed vs. predicted RT deviation per precursor–replicate pair.

    Returns:
        DataFrame sorted by ``abs_rt_dev`` descending, with columns
        ``rt_observed``, ``rt_predicted``, ``rt_dev``, ``abs_rt_dev``.
    """
    for col in (observed_col, predicted_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    base_cols = [c for c in ("Precursor", "Replicate") if c in df.columns]
    result = df[base_cols + [observed_col, predicted_col]].copy().rename(
        columns={observed_col: "rt_observed", predicted_col: "rt_predicted"}
    )
    result["rt_dev"] = (result["rt_observed"] - result["rt_predicted"]).round(6)
    result["abs_rt_dev"] = result["rt_dev"].abs()
    return result.sort_values("abs_rt_dev", ascending=False).reset_index(drop=True)


def summarize_prm(
    df: pd.DataFrame,
    cv_col: str = "RatioLightToHeavy",
    lib_threshold: float = 0.8,
    ratio_threshold: float = 0.8,
    rt_dev_threshold: float = 2.0,
) -> dict:
    """
    Run all PRM quality checks and return a summary dictionary.

    Wraps :func:`compute_cv`, :func:`flag_missing_values`,
    :func:`dot_product_summary`, and :func:`retention_time_deviation`.

    Returns:
        dict with keys: ``cv``, ``missing``, ``dot_products``, ``rt_deviation``,
        ``n_precursors``, ``n_replicates``, ``pct_dot_pass``, ``pct_rt_within``,
        ``median_cv_pct``.
    """
    cv_df = compute_cv(df, value_col=cv_col)
    missing_df = flag_missing_values(df, value_col=cv_col)
    dot_df = dot_product_summary(df, lib_threshold=lib_threshold, ratio_threshold=ratio_threshold)
    rt_df = retention_time_deviation(df)

    pct_dot_pass = dot_df["both_pass"].mean() * 100 if len(dot_df) > 0 else float("nan")
    pct_rt_within = (rt_df["abs_rt_dev"] <= rt_dev_threshold).mean() * 100 if len(rt_df) > 0 else float("nan")
    median_cv = cv_df["cv_pct"].median() if len(cv_df) > 0 else float("nan")

    return {
        "cv": cv_df,
        "missing": missing_df,
        "dot_products": dot_df,
        "rt_deviation": rt_df,
        "n_precursors": df["Precursor"].nunique() if "Precursor" in df.columns else None,
        "n_replicates": df["Replicate"].nunique() if "Replicate" in df.columns else None,
        "pct_dot_pass": round(pct_dot_pass, 2) if pd.notnull(pct_dot_pass) else None,
        "pct_rt_within": round(pct_rt_within, 2) if pd.notnull(pct_rt_within) else None,
        "median_cv_pct": round(median_cv, 2) if pd.notnull(median_cv) else None,
    }

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def plot_library_dot_product_distribution(df: pd.DataFrame) -> None:
    """
    Histogram (+ KDE) of ``Library Dot Product`` colored by ``Isotope Label Type``.
    """
    col = "Library Dot Product"
    color_col = "Isotope Label Type"
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame.")
    if color_col not in df.columns:
        raise KeyError(f"Column '{color_col}' not found in DataFrame.")

    plt.figure(figsize=(10, 6))
    # Loop through each isotop label type and plot each as a separate histogram
    unique_types = df[color_col].dropna().unique()
    for i, iso_type in enumerate(unique_types):
        subset = df[df[color_col] == iso_type]
        sns.histplot(
            subset[col],
            bins=20, kde=True,
            label=str(iso_type),
            alpha=0.5,
            element="step"
        )
    plt.title("Distribution of Library Dot Product by Isotope Label Type")
    plt.xlabel("Library Dot Product")
    plt.legend(title=color_col)
    plt.show()


def plot_heavy_light_scatter(peptide_counts: pd.DataFrame) -> None:
    """
    Scatter plot of heavy vs. light peptide counts, colored by point density.

    Args:
        peptide_counts: DataFrame with ``heavy_count`` and ``light_count`` columns.
    """
    x = peptide_counts["heavy_count"].values
    y = peptide_counts["light_count"].values
    counts = Counter(zip(x, y))
    point_count = np.array([counts[(hx, ly)] for hx, ly in zip(x, y)])
    sort_idx = np.argsort(point_count)

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(
        x[sort_idx], y[sort_idx],
        c=point_count[sort_idx],
        cmap="viridis", alpha=0.7, s=60, edgecolors="k", linewidth=0.5,
    )
    plt.xlabel("Heavy Count")
    plt.ylabel("Light Count")
    plt.title("Heavy vs. Light Peptide Counts\n(colored by number of peptides at each point)")
    min_val, max_val = min(x.min(), y.min()), max(x.max(), y.max())
    plt.plot([min_val, max_val], [min_val, max_val], "r--", lw=1)
    plt.grid(True)
    plt.colorbar(sc, label="# Peptides at Point")
    plt.tight_layout()
    plt.show()


def plot_heavy_light_clusters(clustered_df: pd.DataFrame) -> None:
    """
    Scatter plot of heavy vs. light peptide counts colored by HDBSCAN cluster,
    with the most-abundant cluster highlighted and its data-driven cutoff
    lines drawn.

    Args:
        clustered_df: Output of :func:`cluster_abundant_peptides` (must contain
            ``heavy_count``, ``light_count``, ``cluster``, and ``is_abundant``).
    """
    for col in ("heavy_count", "light_count", "cluster", "is_abundant"):
        if col not in clustered_df.columns:
            raise KeyError(f"Column '{col}' not found; run cluster_abundant_peptides first.")

    plt.figure(figsize=(8, 6))

    noise = clustered_df[clustered_df["cluster"] == -1]
    if not noise.empty:
        plt.scatter(
            noise["heavy_count"], noise["light_count"],
            c="lightgray", s=30, alpha=0.6, label="noise",
        )

    other = clustered_df[(clustered_df["cluster"] != -1) & (~clustered_df["is_abundant"])]
    if not other.empty:
        plt.scatter(
            other["heavy_count"], other["light_count"],
            c=other["cluster"], cmap="tab20", s=40, alpha=0.7, label="other clusters",
        )

    abundant = clustered_df[clustered_df["is_abundant"]]
    plt.scatter(
        abundant["heavy_count"], abundant["light_count"],
        c="crimson", s=60, edgecolors="k", linewidth=0.5, label=f"abundant (n={len(abundant)})",
    )

    if not abundant.empty:
        # min(heavy, light) >= c defines the abundant cluster, so the same
        # scalar cutoff applies to both axes.
        c = abundant[["heavy_count", "light_count"]].min(axis=1).min()
        plt.axvline(c, color="crimson", ls="--", lw=1)
        plt.axhline(c, color="crimson", ls="--", lw=1)

    plt.xlabel("Heavy Count")
    plt.ylabel("Light Count")
    plt.title("Heavy vs. Light Peptide Counts\n(HDBSCAN clusters, abundant cluster highlighted)")
    plt.legend(loc="best", fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_peptide_counts(df: pd.DataFrame) -> None:
    """Histogram (+ KDE) of heavy peptide counts."""
    plt.figure(figsize=(10, 6))
    sns.histplot(df["heavy_count"], bins=20, kde=True)
    plt.title("Distribution of Heavy Peptide Counts")
    plt.xlabel("Heavy Peptide Counts")
    plt.show()


def plot_pool_boxplot(
    pool_df: pd.DataFrame,
    col_ratio: str = "RatioLightToHeavy",
    col_plate: str = "characteristics[plate]",
) -> None:
    """
    Boxplot of *col_ratio* (log scale) by Replicate, colored by Plate.

    Args:
        pool_df: DataFrame with ``Replicate``, *col_ratio*, and *col_plate* columns.
        col_ratio: Column for the light-to-heavy ratio.
        col_plate: Column for plate labels.
    """
    pool_df = pool_df[["Replicate", col_ratio, col_plate]].drop_duplicates()
    plt.figure(figsize=(12, 8))
    sns.boxplot(x="Replicate", y=col_ratio, data=pool_df, hue=col_plate)
    plt.title(f"Boxplot of {col_ratio} colored by {col_plate}")
    plt.xlabel(col_plate)
    plt.ylabel(col_ratio)
    plt.yscale("log")
    plt.xticks([], [])
    plt.legend(title=col_plate, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

def plot_pool_heatmap(pool_data: pd.DataFrame, aggfunc: str = "mean") -> None:
    """
    Heatmap of log(RatioLightToHeavy) for Pool samples.

    Rows (peptides) and columns (replicates) are ordered by their mean log-ratio.
    Duplicate peptide/replicate pairs are aggregated with *aggfunc*.

    Args:
        pool_data: DataFrame with ``Peptide Sequence``, ``Replicate``, ``RatioLightToHeavy``.
        aggfunc: Aggregation passed to ``pivot_table``. Default ``"mean"``.
    """
    required = ["Peptide Sequence", "Replicate", "RatioLightToHeavy"]
    missing = [c for c in required if c not in pool_data.columns]
    if missing:
        raise KeyError(f"pool_data missing columns {missing}. Found: {list(pool_data.columns)}")

    pivot = pool_data.pivot_table(
        index="Peptide Sequence", columns="Replicate", values="RatioLightToHeavy", aggfunc=aggfunc
    )
    heatmap_data = np.log(pivot)
    heatmap_data = heatmap_data.loc[
        heatmap_data.mean(axis=1).sort_values().index,
        heatmap_data.mean(axis=0).sort_values().index,
    ]
    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, cmap="viridis")
    plt.title("Heatmap of log(RatioLightToHeavy) for Pool")
    plt.xlabel("")
    plt.ylabel("Peptide Sequence")
    plt.xticks([], [])
    plt.show()


def plot_intra_plate_cv_stats(
    peptide_plate_stats: pd.DataFrame,
    col_name: str = "characteristics[plate]",
) -> None:
    """
    Boxplot of intra-plate CV per plate.

    Args:
        peptide_plate_stats: Output of :func:`calculate_intra_plate_cv`.
        col_name: Group column used in the stats DataFrame.
    """
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=col_name, y="intra_plate_cv", data=peptide_plate_stats)
    plt.title("Intra-Plate CV for Pool")
    plt.xlabel(col_name.replace("characteristics[", "").replace("]", "").replace("plate", "Plate").capitalize())
    plt.ylabel("Intra-Plate CV")
    plt.show()


def plot_inter_plate_cv_kde(interplate_cv: pd.DataFrame) -> None:
    """
    KDE of inter-plate CV across peptides, with a vertical line at the median.

    Args:
        interplate_cv: Output of :func:`calculate_inter_plate_cv`.
            Must have an ``inter_plate_cv`` column.
    """
    if "inter_plate_cv" not in interplate_cv.columns:
        raise KeyError(f"Expected column 'inter_plate_cv'. Found: {list(interplate_cv.columns)}")
    median_cv = interplate_cv["inter_plate_cv"].median()
    plt.figure(figsize=(12, 6))
    sns.kdeplot(interplate_cv["inter_plate_cv"].dropna(), fill=True)
    plt.axvline(median_cv, color="red", linestyle="--", label=f"Median = {median_cv:.2f}")
    plt.title("KDE of Inter-Plate CV Across Peptides")
    plt.xlabel("Inter-Plate CV")
    plt.ylabel("Density")
    plt.legend()
    plt.show()


def plot_cumulative_peptide_count_by_cv(peptide_means: pd.DataFrame) -> None:
    """
    Cumulative peptide count as a function of sorted inter-plate CV.

    Args:
        peptide_means: DataFrame with ``inter_plate_cv`` and ``Peptide Sequence`` columns.
    """
    cv_sorted = (
        peptide_means[["inter_plate_cv", "Peptide Sequence"]]
        .sort_values("inter_plate_cv")
        .reset_index(drop=True)
    )
    cv_sorted["cumulative_count"] = range(1, len(cv_sorted) + 1)

    print(f"Total peptides: {len(cv_sorted)}")
    for thresh in [0.10, 0.20]:
        n = (cv_sorted["inter_plate_cv"] < thresh).sum()
        print(f"  CV < {int(thresh * 100)}%: {n} peptides")

    plt.figure(figsize=(10, 6))
    plt.plot(cv_sorted["inter_plate_cv"], cv_sorted["cumulative_count"], marker="o", linestyle="-")
    plt.xlabel("Inter-Plate CV")
    plt.ylabel("Cumulative Number of Peptides")
    plt.title("Cumulative Peptide Count by Inter-Plate CV")

    for cv_mark in [0.1 * i for i in range(1, 11)]:
        mask = cv_sorted["inter_plate_cv"] >= cv_mark
        if mask.any():
            idx = mask.idxmax()
            x_val = cv_sorted.at[idx, "inter_plate_cv"]
            y_val = cv_sorted.at[idx, "cumulative_count"]
            plt.axvline(x_val, color="gray", linestyle="--", linewidth=0.8)
            plt.text(x_val, y_val, f"{int(y_val)} peptides\n{cv_mark:.1f} CV",
                     va="bottom", ha="left", fontsize=9, color="blue")

    plt.tight_layout()
    plt.show()


def plot_logratio_by_plate_boxplot(df: pd.DataFrame) -> None:
    """
    Boxplot of ``log_ratio`` by Replicate, colored by Plate.

    Args:
        df: DataFrame with columns ``['Replicate', 'log_ratio', 'Plate']``.
    """
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="Replicate", y="log_ratio", data=df.sort_values("Plate"), hue="Plate", dodge=False)
    plt.title("log(RatioLightToHeavy) by Replicate (colored by Plate)")
    plt.xlabel("")
    plt.ylabel("log(RatioLightToHeavy)")
    plt.xticks([], [])
    plt.legend(title="Plate", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Absolute Concentration Plots
# ---------------------------------------------------------------------------


def map_peptide_sequence(
    abs_df: pd.DataFrame,
    fasta_file: pd.DataFrame,
    protein_name: str,
    line_length: int = 30,
):
    """
    Plots for showing peptide-to-protein mapping.

    Produces two plots:
    1. Protein sequence as blocks, with highlighted detected peptide locations.
    2. Coverage bar indicating which sequence regions are covered by detected peptides.

    Args:
        abs_df: Peptide concentration data (wide or long format from get_absolute_conc).
        fasta_file: Fasta/sequence dataframe (output of fetch_fasta).
        protein_name: Protein name matching abs_df['Protein Name'].
        line_length: Amino acids per row.

    Returns:
        matplotlib.figure.Figure
    """
    if isinstance(abs_df.columns, pd.MultiIndex):
        abs_df_flat = abs_df.copy()
        abs_df_flat.columns = [
            '|'.join([str(x) for x in col if x != '' and x is not None])
            for col in abs_df_flat.columns.values
        ]
    else:
        abs_df_flat = abs_df

    def _get_col(df, colname):
        if colname in df.columns:
            return df[colname]
        elif df.index.names and colname in df.index.names:
            return df.index.get_level_values(colname)
        raise KeyError(f"'{colname}' not found in DataFrame columns or index.")

    try:
        protein_mask = _get_col(abs_df_flat, 'Protein Name') == protein_name
        if not protein_mask.any():
            raise ValueError(f"No rows found for Protein Name '{protein_name}'.")
        peptides = _get_col(abs_df_flat, 'Peptide Sequence')[protein_mask].unique()
    except Exception as e:
        raise RuntimeError(
            f"Could not extract peptide sequences for '{protein_name}': {e}"
        )

    matching = fasta_file[fasta_file['id'] == protein_name]
    if matching.empty:
        stripped = protein_name.split('|')[0]
        matching = fasta_file[fasta_file['id'].str.startswith(stripped)]
        if matching.empty:
            raise ValueError(f"Could not find fasta row for '{protein_name}'.")
    sequence = matching.iloc[0]['sequence']
    seq_len = len(sequence)
    n_lines = int(np.ceil(seq_len / line_length))

    pep_positions = np.zeros(seq_len, dtype=int)
    pep_indices: dict = {}
    for idx, pep in enumerate(peptides):
        start = sequence.find(pep)
        if start == -1:
            print(
                f"Warning: peptide '{pep}' not found in sequence for '{protein_name}'."
            )
            continue
        pep_indices[pep] = idx
        pep_positions[start : start + len(pep)] = idx + 1

    # Plot 1: Block diagram, as originally implemented
    fig_height = max(n_lines * 1.3 + 1.5, 2.8)
    fig_width = min(1.3 * line_length, 28)
    fig, axes = plt.subplots(
        n_lines, 1, figsize=(fig_width, fig_height), sharex=False, squeeze=False
    )
    fig.suptitle(str(protein_name), fontsize=18, y=1.0, weight='bold')
    colors = plt.cm.tab20.colors

    for line in range(n_lines):
        ax = axes[line, 0]
        start_idx = line * line_length
        end_idx = min((line + 1) * line_length, seq_len)
        n_aa = end_idx - start_idx
        seq_sub = sequence[start_idx:end_idx]

        for i in range(line_length):
            if i < n_aa:
                in_pep = pep_positions[start_idx + i]
                ax.add_patch(
                    plt.Rectangle(
                        (i, 0),
                        1,
                        1,
                        facecolor=colors[in_pep - 1] if in_pep else '#ffffff',
                        edgecolor='black',
                        lw=1.5 if in_pep else 0.7,
                        alpha=0.7 if in_pep else 1.0,
                        zorder=2 if in_pep else 1,
                        linewidth=1.3,
                    )
                )
                ax.text(
                    i + 0.5,
                    0.5,
                    seq_sub[i],
                    ha='center',
                    va='center',
                    color='black',
                    fontsize=16,
                    fontfamily='monospace',
                    weight='bold',
                    zorder=4,
                )
            else:
                ax.add_patch(
                    plt.Rectangle(
                        (i, 0),
                        1,
                        1,
                        facecolor='#f5f5f5',
                        edgecolor='black',
                        lw=0.6,
                        alpha=1.0,
                        zorder=1,
                        linewidth=1.0,
                    )
                )

        ax.set_xlim(0, line_length)
        ax.set_ylim(0, 1.23)
        ax.set_yticks([])
        residue_ticks = list(range(0, line_length, 10))
        xtick_labels = [
            str(start_idx + 1 + x) if x < n_aa else '' for x in residue_ticks
        ]
        ax.set_xticks([x + 0.5 for x in residue_ticks])
        ax.set_xticklabels(xtick_labels, fontsize=12)
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_visible(True)

    if pep_indices:
        legend_handles = [
            plt.Line2D([0], [0], color=colors[idx % len(colors)], lw=8)
            for idx in pep_indices.values()
        ]
        axes[-1, 0].legend(
            legend_handles,
            list(pep_indices.keys()),
            loc='upper center',
            bbox_to_anchor=(0.5, -0.30 + (0.3 / fig_height)),
            ncol=3,
            fontsize=11,
            frameon=False,
        )
    plt.tight_layout(h_pad=0.8, rect=[0, 0, 1, 0.98])

    # Plot 2: Coverage bar, simple 1d highlight
    cov_fig, cov_ax = plt.subplots(figsize=(max(7, min(seq_len // 5, 18)), 1.1))
    cov_ax.set_title('Peptide Coverage', fontsize=14, weight='bold', y=1.2)
    cov_ax.set_xlim(0, seq_len)
    cov_ax.set_ylim(0, 1)
    cov_ax.axis('off')

    # Draw full protein bar
    cov_ax.add_patch(
        plt.Rectangle((0, 0.4), seq_len, 0.2, facecolor='#E0E0E0', edgecolor='black', lw=1.1)
    )

    # Overlay covered regions by detected peptides
    cov_colors = [colors[idx % len(colors)] for idx in range(len(peptides))]
    for i, pep in enumerate(peptides):
        # Find all occurrences for peptides that may appear more than once
        starts = []
        s = 0
        while True:
            found = sequence.find(pep, s)
            if found == -1:
                break
            starts.append(found)
            s = found + 1
        for start in starts:
            cov_ax.add_patch(
                plt.Rectangle(
                    (start, 0.4),
                    len(pep),
                    0.2,
                    facecolor=cov_colors[i],
                    edgecolor='#202020',
                    lw=2,
                    zorder=2,
                    alpha=0.8,
                )
            )

    cov_ax.text(
        0,
        1.01,
        f"{protein_name} (length={seq_len})",
        va='bottom',
        ha='left',
        fontsize=12,
        fontweight='bold',
        fontfamily='monospace',
        color='#2C3140',
    )

    ticks = list(range(0, seq_len + 1, 50)) if seq_len > 150 else list(range(0, seq_len + 1, 20))
    cov_ax.set_xticks(ticks)
    cov_ax.set_xticklabels([str(t + 1) for t in ticks], fontsize=11)
    cov_ax.tick_params(axis='x', which='both', length=0)

    plt.tight_layout()

    return fig, cov_fig


def plot_peptide_concentration_by_group(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    protein_name: str,
    group_col: Optional[str] = None,
    color_col: Optional[str] = None,
) -> None:
    """
    Boxplot of absolute peptide concentration by group, one subplot per peptide.

    Args:
        abs_df: Wide output from get_absolute_conc (MultiIndex + replicate columns).
        sdrf_data_file: SDRF table with 'source name' and group_col.
        protein_name: Filter to this 'Protein Name' value.
        group_col: SDRF column to group x-axis by. Defaults to the last 'factor value[...]' column.
        color_col: SDRF column to color boxes by. Defaults to group_col.
    """
    group_col = _resolve_group_col(sdrf_data_file, group_col)
    color_col = color_col if color_col is not None else group_col
    for col in (group_col, color_col):
        if col not in sdrf_data_file.columns:
            raise KeyError(f"Column '{col}' not found in sdrf_data_file.")
    if "source name" not in sdrf_data_file.columns:
        raise KeyError("Column 'source name' not found in sdrf_data_file.")

    plot_df = _build_concentration_plot_df(abs_df, sdrf_data_file, group_col, color_col)

    file_key = _resolve_sdrf_file_key(sdrf_data_file)
    unmapped = set(plot_df["Replicate"]) - set(sdrf_data_file[file_key])
    if unmapped:
        print(f"Unmapped replicates ({len(unmapped)}): {sorted(unmapped)}")

    protein_df = plot_df[
        (plot_df["Protein Name"] == protein_name) & plot_df[color_col].notna()
    ].copy()
    if protein_df.empty:
        print(f"No data found for Protein Name: '{protein_name}'")
        return

    group2color = _make_group2color(protein_df[color_col].dropna().unique(), color_col)
    _plot_protein_boxes(protein_df, protein_name, group_col, color_col, group2color)


def plot_peptide_all(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    protein_name: str,
    group_col: Optional[str] = None,
    color_col: Optional[str] = None,
) -> None:
    """
    Combined view per protein: one boxplot+strip subplot per peptide, plus a
    median line plot (no error bars) of all peptides at the bottom.

    Args:
        abs_df: Wide output from get_absolute_conc.
        sdrf_data_file: SDRF table with 'source name' and group_col.
        protein_name: Filter to this 'Protein Name' value.
        group_col: SDRF column to group x-axis by.
        color_col: SDRF column to color boxes by. Defaults to group_col.
    """
    group_col = _resolve_group_col(sdrf_data_file, group_col)
    color_col = color_col if color_col is not None else group_col
    for col in (group_col, color_col):
        if col not in sdrf_data_file.columns:
            raise KeyError(f"Column '{col}' not found in sdrf_data_file.")
    if "source name" not in sdrf_data_file.columns:
        raise KeyError("Column 'source name' not found in sdrf_data_file.")

    plot_df = _build_concentration_plot_df(abs_df, sdrf_data_file, group_col, color_col)
    protein_df = plot_df[
        (plot_df["Protein Name"] == protein_name) & plot_df[color_col].notna()
    ].copy()
    if protein_df.empty:
        print(f"No data found for Protein Name: '{protein_name}'")
        return

    peptides = protein_df["Peptide Sequence"].unique()
    group2color = _make_group2color(protein_df[color_col].dropna().unique(), color_col)

    group_to_color_label = (
        protein_df.drop_duplicates(subset=[group_col])
        .set_index(group_col)[color_col]
    )
    if color_col == "characteristics[disease category]":
        cat_rank = {c: i for i, c in enumerate(_DISEASE_CATEGORY_PALETTE)}
        key_fn = lambda g: (cat_rank.get(group_to_color_label.get(g, ""), len(cat_rank)), str(g))
    else:
        key_fn = lambda g: (str(group_to_color_label.get(g, "")), str(g))
    group_order = sorted(protein_df[group_col].dropna().unique(), key=key_fn)

    agg = (
        protein_df.groupby(["Peptide Sequence", group_col])["Protein conc [pmol]"]
        .median()
        .reset_index()
    )
    x_positions = {g: i for i, g in enumerate(group_order)}

    n_rows = len(peptides)
    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3 * n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes]

    for ax, pep in zip(axes, peptides):
        pep_df = protein_df[protein_df["Peptide Sequence"] == pep]
        sns.boxplot(
            data=pep_df, x=group_col, y="Protein conc [pmol]",
            hue=color_col, ax=ax, order=group_order, palette=group2color,
            dodge=False,
        )
        sns.stripplot(
            data=pep_df, x=group_col, y="Protein conc [pmol]",
            ax=ax, order=group_order, color="grey",
            dodge=False, jitter=True, alpha=0.5, size=3,
        )
        pep_agg = agg[agg["Peptide Sequence"] == pep].set_index(group_col).reindex(group_order)
        ax.plot(
            [x_positions[g] for g in group_order],
            pep_agg["Protein conc [pmol]"].values,
            color="orange", linewidth=1.2, linestyle="--", marker="o", markersize=3, zorder=5,
        )
        ax.set_title(f"{pep}|{protein_name}", fontsize=10, loc='left')
        ax.set_ylabel("Protein conc [pmol]", fontsize=5)
        ax.set_xlabel(group_col if ax == axes[-1] else "", fontsize=5)
        ax.tick_params(axis='x', rotation=90, labelsize=5)
        ax.tick_params(axis='y', labelsize=5)
        if ax.get_legend() is not None:
            ax.legend_.remove()

    handles, labels = axes[0].get_legend_handles_labels()
    if color_col == "characteristics[disease category]":
        cat_rank = {c: i for i, c in enumerate(_DISEASE_CATEGORY_PALETTE)}
        paired = sorted(
            zip(labels, handles),
            key=lambda x: cat_rank.get(x[0], len(cat_rank)),
        )
        labels, handles = zip(*paired) if paired else (labels, handles)
    fig.legend(
        handles, labels,
        loc="center right", bbox_to_anchor=(1.15, 0.5),
        fontsize=5, title=color_col, title_fontsize=5,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()


def plot_median_peptide_concentration_by_group(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    protein_name: str,
    group_col: Optional[str] = None,
) -> None:
    """
    Line plot of median peptide concentration ± SEM by group for each peptide of a protein.

    Args:
        abs_df: Peptide concentration data (wide format from get_absolute_conc).
        sdrf_data_file: Sample metadata with 'source name' and group_col.
        protein_name: Protein Name value to filter on.
        group_col: Metadata column to group samples by. Defaults to the last 'factor value[...]' column.
    """
    group_col = _resolve_group_col(sdrf_data_file, group_col)
    if group_col not in sdrf_data_file.columns:
        raise KeyError(f"Column '{group_col}' not found in sdrf_data_file.")
    if "source name" not in sdrf_data_file.columns:
        raise KeyError("Column 'source name' not found in sdrf_data_file.")

    plot_df = _build_concentration_plot_df(abs_df, sdrf_data_file, group_col, group_col)

    protein_df = plot_df[plot_df["Protein Name"] == protein_name].copy()
    if protein_df.empty:
        print(f"No data found for Protein Name: '{protein_name}'")
        return

    peptides = protein_df["Peptide Sequence"].unique()
    peptide_colors = dict(zip(peptides, sns.color_palette("tab10", len(peptides))))

    agg_stats = (
        protein_df
        .groupby(["Peptide Sequence", group_col])["Protein conc [pmol]"]
        .agg(['median', 'std', 'count'])
        .reset_index()
    )
    agg_stats["sem"] = agg_stats["std"] / agg_stats["count"].pow(0.5)
    group_order = sorted(agg_stats[group_col].dropna().unique(), key=str)

    plt.figure(figsize=(12, 8))
    for pep in peptides:
        pep_stats = agg_stats[agg_stats["Peptide Sequence"] == pep].set_index(group_col).reindex(group_order)
        plt.errorbar(
            group_order,
            pep_stats["median"].values,
            yerr=pep_stats["sem"].values,
            marker="o", label=pep, color=peptide_colors[pep],
            capsize=4, linestyle='-',
        )

    plt.title(f"Median peptide concentration by group for {protein_name}")
    plt.xlabel(group_col)
    plt.ylabel("Median Protein/Peptide conc [pmol]")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title="Peptide Sequence")
    plt.tight_layout()
    plt.show()


def plot_all_median_peptide_concentration_by_group(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    group_col: Optional[str] = None,
    pdf_path: str = "median_peptide_concentration_by_group.pdf",
) -> None:
    """
    Save one median-concentration line plot per protein to a multi-page PDF (A4 landscape).

    Args:
        abs_df: Wide output from get_absolute_conc.
        sdrf_data_file: SDRF metadata table.
        group_col: Column to group samples by. Defaults to the last 'factor value[...]' column.
        pdf_path: Output PDF path.
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import warnings

    group_col = _resolve_group_col(sdrf_data_file, group_col)
    plot_df = _build_concentration_plot_df(abs_df, sdrf_data_file, group_col, group_col)
    protein_names = plot_df["Protein Name"].dropna().unique()
    total = len(protein_names)

    A4_WIDTH, A4_HEIGHT = 11.69, 8.27
    original_show = plt.show
    plt.show = lambda *a, **kw: None
    skipped = []
    saved = 0
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
            with PdfPages(pdf_path) as pdf:
                for i, pname in enumerate(protein_names, 1):
                    print(f"Generating plot {i}/{total}: {pname}")
                    protein_df = plot_df[plot_df["Protein Name"] == pname].copy()
                    if protein_df.empty:
                        skipped.append(pname)
                        continue
                    try:
                        peptides = protein_df["Peptide Sequence"].unique()
                        peptide_colors = dict(zip(peptides, sns.color_palette("tab10", len(peptides))))
                        agg_stats = (
                            protein_df
                            .groupby(["Peptide Sequence", group_col])["Protein conc [pmol]"]
                            .agg(['median', 'std', 'count'])
                            .reset_index()
                        )
                        agg_stats["sem"] = agg_stats["std"] / agg_stats["count"].pow(0.5)
                        group_order = sorted(agg_stats[group_col].dropna().unique(), key=str)
                        plt.figure(figsize=(A4_WIDTH, A4_HEIGHT))
                        for pep in peptides:
                            pep_stats = agg_stats[agg_stats["Peptide Sequence"] == pep].set_index(group_col).reindex(group_order)
                            plt.errorbar(
                                group_order, pep_stats["median"].values, yerr=pep_stats["sem"].values,
                                marker="o", label=pep, color=peptide_colors[pep], capsize=4, linestyle='-',
                            )
                        plt.title(f"Median peptide concentration by group for {pname}")
                        plt.xlabel(group_col)
                        plt.ylabel("Median Protein/Peptide conc [pmol]")
                        plt.xticks(rotation=45, ha='right')
                        plt.legend(title="Peptide Sequence")
                        plt.tight_layout()
                        fig = plt.gcf()
                        pdf.savefig(fig)
                        plt.close(fig)
                        saved += 1
                    except Exception as e:
                        print(f"  Skipping {pname}: {e}")
                        skipped.append(pname)
    finally:
        plt.show = original_show

    print(f"Saved {saved} plots to {pdf_path}")
    if skipped:
        print(f"Skipped {len(skipped)} protein(s): {skipped}")


def plot_all_peptide_concentration_by_group(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    group_col: Optional[str] = None,
    color_col: Optional[str] = None,
    pdf_path: str = "peptide_concentration_by_group.pdf",
) -> str:
    """
    Save one boxplot per protein to a multi-page PDF (A4 portrait).

    Args:
        abs_df: Wide output from get_absolute_conc.
        sdrf_data_file: SDRF metadata table.
        group_col: Column to group samples by. Defaults to the last 'factor value[...]' column.
        color_col: Column to color boxes by. Defaults to group_col.
        pdf_path: Output PDF path.

    Returns:
        pdf_path
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import warnings

    group_col = _resolve_group_col(sdrf_data_file, group_col)
    color_col = color_col if color_col is not None else group_col

    plot_df = _build_concentration_plot_df(abs_df, sdrf_data_file, group_col, color_col)
    plot_df = plot_df[plot_df[color_col].notna()]

    file_key = _resolve_sdrf_file_key(sdrf_data_file)
    unmapped = set(plot_df["Replicate"]) - set(sdrf_data_file[file_key])
    if unmapped:
        print(f"Unmapped replicates ({len(unmapped)}): {sorted(unmapped)}")

    group2color = _make_group2color(plot_df[color_col].dropna().unique(), color_col)
    protein_names = plot_df["Protein Name"].dropna().unique()
    total = len(protein_names)

    A4_WIDTH, A4_HEIGHT = 8.27, 11.69
    original_show = plt.show
    plt.show = lambda *a, **kw: None
    saved = 0
    skipped = []
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
            with PdfPages(pdf_path) as pdf:
                for i, pname in enumerate(protein_names, 1):
                    print(f"Generating plot {i}/{total}: {pname}")
                    protein_df = plot_df[plot_df["Protein Name"] == pname].copy()
                    if protein_df.empty:
                        skipped.append(pname)
                        continue
                    try:
                        plt.figure(figsize=(A4_WIDTH, A4_HEIGHT))
                        _plot_protein_boxes(protein_df, pname, group_col, color_col, group2color)
                        fig = plt.gcf()
                        pdf.savefig(fig)
                        plt.close(fig)
                        saved += 1
                    except Exception as e:
                        print(f"  Skipping {pname}: {e}")
                        skipped.append(pname)
    finally:
        plt.show = original_show

    print(f"Saved {saved} plots to {pdf_path}")
    if skipped:
        print(f"Skipped {len(skipped)} protein(s): {skipped}")
    return pdf_path


def plot_pca(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    color_col: Optional[str] = None,
    feature_col: str = "Protein Name",
    figsize: Tuple[int, int] = (8, 6),
    title: str = "PCA",
) -> None:
    """
    PCA scatter plot of samples colored by a metadata column.

    Args:
        abs_df: Wide output from get_absolute_conc.
        sdrf_data_file: SDRF metadata table with 'source name' and color_col.
        color_col: Column in sdrf_data_file to color points by. Defaults to last 'factor value[...]' column.
        feature_col: Feature axis for the matrix ('Protein Name' or 'Peptide Sequence').
        figsize: Figure size.
        title: Plot title.
    """
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        raise ImportError("scikit-learn is required: pip install scikit-learn")

    color_col = _resolve_group_col(sdrf_data_file, color_col)
    if color_col not in sdrf_data_file.columns:
        raise KeyError(f"Column '{color_col}' not found in sdrf_data_file.")
    if "source name" not in sdrf_data_file.columns:
        raise KeyError("Column 'source name' not found in sdrf_data_file.")

    long_df = (
        abs_df.reset_index()
        .melt(id_vars=["qRePS", "Peptide Sequence", "Protein Name"],
              var_name="Replicate", value_name="Protein conc [pmol]")
        .dropna(subset=["Protein conc [pmol]"])
    )
    pivot = long_df.pivot_table(
        index="Replicate", columns=feature_col, values="Protein conc [pmol]", aggfunc="mean"
    )
    pivot = pivot.dropna(axis=1, thresh=max(1, int(0.5 * len(pivot))))
    pivot = pivot.fillna(pivot.median())

    if pivot.shape[0] < 2:
        raise ValueError("Not enough samples for PCA after dropping missing values.")

    X = StandardScaler().fit_transform(pivot)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)
    var_explained = pca.explained_variance_ratio_

    file_key = _resolve_sdrf_file_key(sdrf_data_file)
    rep_meta = sdrf_data_file[[file_key, color_col]].drop_duplicates().set_index(file_key)
    pca_df = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1]}, index=pivot.index)
    pca_df[color_col] = rep_meta.reindex(pca_df.index)[color_col].values

    group2color = _make_group2color(pca_df[color_col].dropna().unique(), color_col)

    plt.figure(figsize=figsize)
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue=color_col, palette=group2color, s=80, edgecolor="k")
    plt.xlabel(f"PC1 ({100 * var_explained[0]:.1f}%)")
    plt.ylabel(f"PC2 ({100 * var_explained[1]:.1f}%)")
    plt.title(title)
    plt.legend(title=color_col, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


def plot_umap(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    color_col: Optional[str] = None,
    feature_col: str = "Protein Name",
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    metric: str = "euclidean",
    random_state: int = 42,
    figsize: Tuple[int, int] = (8, 6),
    title: str = "UMAP",
) -> None:
    """
    UMAP scatter plot of samples colored by a metadata column.

    Args:
        abs_df: Wide output from get_absolute_conc.
        sdrf_data_file: SDRF metadata table with 'source name' and color_col.
        color_col: Column in sdrf_data_file to color points by. Defaults to last 'factor value[...]' column.
        feature_col: Feature axis for the matrix ('Protein Name' or 'Peptide Sequence').
        n_neighbors: UMAP neighborhood size.
        min_dist: UMAP minimum distance between embedded points.
        metric: Distance metric for UMAP.
        random_state: Random seed for reproducibility.
        figsize: Figure size.
        title: Plot title.
    """
    try:
        import umap as umap_lib
    except ImportError:
        raise ImportError("umap-learn is required: pip install umap-learn")
    try:
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        raise ImportError("scikit-learn is required: pip install scikit-learn")

    color_col = _resolve_group_col(sdrf_data_file, color_col)
    if color_col not in sdrf_data_file.columns:
        raise KeyError(f"Column '{color_col}' not found in sdrf_data_file.")
    if "source name" not in sdrf_data_file.columns:
        raise KeyError("Column 'source name' not found in sdrf_data_file.")

    long_df = (
        abs_df.reset_index()
        .melt(id_vars=["qRePS", "Peptide Sequence", "Protein Name"],
              var_name="Replicate", value_name="Protein conc [pmol]")
        .dropna(subset=["Protein conc [pmol]"])
    )
    pivot = long_df.pivot_table(
        index="Replicate", columns=feature_col, values="Protein conc [pmol]", aggfunc="mean"
    )
    pivot = pivot.dropna(axis=1, thresh=max(1, int(0.5 * len(pivot))))
    pivot = pivot.fillna(pivot.median())

    if pivot.shape[0] < 2:
        raise ValueError("Not enough samples for UMAP after dropping missing values.")

    X = StandardScaler().fit_transform(pivot)
    embedding = umap_lib.UMAP(
        n_neighbors=n_neighbors, min_dist=min_dist, n_components=2,
        metric=metric, random_state=random_state,
    ).fit_transform(X)

    file_key = _resolve_sdrf_file_key(sdrf_data_file)
    rep_meta = sdrf_data_file[[file_key, color_col]].drop_duplicates().set_index(file_key)
    umap_df = pd.DataFrame({"UMAP1": embedding[:, 0], "UMAP2": embedding[:, 1]}, index=pivot.index)
    umap_df[color_col] = rep_meta.reindex(umap_df.index)[color_col].values

    group2color = _make_group2color(umap_df[color_col].dropna().unique(), color_col)

    plt.figure(figsize=figsize)
    sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue=color_col, palette=group2color, s=80, edgecolor="k")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(title)
    plt.legend(title=color_col, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


def plot_all_all(
    abs_df: pd.DataFrame,
    sdrf_data_file: pd.DataFrame,
    group_col: Optional[str] = None,
    color_col: Optional[str] = None,
    pdf_path: str = "peptide_all_by_group.pdf",
) -> str:
    """
    Save one combined boxplot+strip+median-line figure per protein to a multi-page PDF.

    Args:
        abs_df: Wide output from get_absolute_conc.
        sdrf_data_file: SDRF metadata table.
        group_col: Column to group samples by. Defaults to the last 'factor value[...]' column.
        color_col: Column to color boxes by. Defaults to group_col.
        pdf_path: Output PDF path.

    Returns:
        pdf_path
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import warnings

    group_col = _resolve_group_col(sdrf_data_file, group_col)
    color_col = color_col if color_col is not None else group_col

    plot_df = _build_concentration_plot_df(abs_df, sdrf_data_file, group_col, color_col)
    plot_df = plot_df[plot_df[color_col].notna()]

    file_key = _resolve_sdrf_file_key(sdrf_data_file)
    unmapped = set(plot_df["Replicate"]) - set(sdrf_data_file[file_key])
    if unmapped:
        print(f"Unmapped replicates ({len(unmapped)}): {sorted(unmapped)}")

    group2color = _make_group2color(plot_df[color_col].dropna().unique(), color_col)

    group_to_color_label = (
        plot_df.drop_duplicates(subset=[group_col])
        .set_index(group_col)[color_col]
    )
    if color_col == "characteristics[disease category]":
        cat_rank = {c: i for i, c in enumerate(_DISEASE_CATEGORY_PALETTE)}
        key_fn = lambda g: (cat_rank.get(group_to_color_label.get(g, ""), len(cat_rank)), str(g))
    else:
        key_fn = lambda g: (str(group_to_color_label.get(g, "")), str(g))
    group_order = sorted(plot_df[group_col].dropna().unique(), key=key_fn)
    x_positions = {g: i for i, g in enumerate(group_order)}

    agg_all = (
        plot_df.groupby(["Protein Name", "Peptide Sequence", group_col])["Protein conc [pmol]"]
        .median()
        .reset_index()
    )

    protein_names = plot_df["Protein Name"].dropna().unique()
    total = len(protein_names)

    A4_WIDTH, A4_HEIGHT = 8.27, 11.69
    original_show = plt.show
    plt.show = lambda *a, **kw: None
    saved = 0
    skipped = []
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
            with PdfPages(pdf_path) as pdf:
                for i, pname in enumerate(protein_names, 1):
                    print(f"Generating plot {i}/{total}: {pname}")
                    protein_df = plot_df[plot_df["Protein Name"] == pname].copy()
                    if protein_df.empty:
                        skipped.append(pname)
                        continue
                    try:
                        peptides = protein_df["Peptide Sequence"].unique()
                        agg = agg_all[agg_all["Protein Name"] == pname]
                        pep_group_order = [g for g in group_order if g in protein_df[group_col].values]

                        n_rows = len(peptides)
                        fig, axes = plt.subplots(n_rows, 1, figsize=(A4_WIDTH, max(A4_HEIGHT, 3 * n_rows)), sharex=True)
                        if n_rows == 1:
                            axes = [axes]

                        for ax, pep in zip(axes, peptides):
                            pep_df = protein_df[protein_df["Peptide Sequence"] == pep]
                            sns.boxplot(
                                data=pep_df, x=group_col, y="Protein conc [pmol]",
                                hue=color_col, ax=ax, order=pep_group_order, palette=group2color,
                                dodge=False,
                            )
                            sns.stripplot(
                                data=pep_df, x=group_col, y="Protein conc [pmol]",
                                ax=ax, order=pep_group_order, color="grey",
                                dodge=False, jitter=True, alpha=0.5, size=3,
                            )
                            pep_agg = agg[agg["Peptide Sequence"] == pep].set_index(group_col).reindex(pep_group_order)
                            ax.plot(
                                [x_positions[g] for g in pep_group_order],
                                pep_agg["Protein conc [pmol]"].values,
                                color="orange", linewidth=1.2, linestyle="--", marker="o", markersize=3, zorder=5,
                            )
                            ax.set_title(f"{pep}|{pname}", fontsize=10, loc='left')
                            ax.set_ylabel("Protein conc [pmol]", fontsize=5)
                            ax.set_xlabel(group_col if ax == axes[-1] else "", fontsize=5)
                            ax.tick_params(axis='x', rotation=90, labelsize=5)
                            ax.tick_params(axis='y', labelsize=5)
                            if ax.get_legend() is not None:
                                ax.legend_.remove()

                        handles, labels = axes[0].get_legend_handles_labels()
                        if color_col == "characteristics[disease category]":
                            cat_rank = {c: i for i, c in enumerate(_DISEASE_CATEGORY_PALETTE)}
                            paired = sorted(
                                zip(labels, handles),
                                key=lambda x: cat_rank.get(x[0], len(cat_rank)),
                            )
                            labels, handles = zip(*paired) if paired else (labels, handles)
                        fig.legend(
                            handles, labels,
                            loc="center right", bbox_to_anchor=(1.15, 0.5),
                            fontsize=5, title=color_col, title_fontsize=5,
                        )
                        plt.tight_layout(rect=[0, 0, 1, 0.98])
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)
                        saved += 1
                    except Exception as e:
                        print(f"  Skipping {pname}: {e}")
                        skipped.append(pname)
    finally:
        plt.show = original_show

    print(f"Saved {saved} plots to {pdf_path}")
    if skipped:
        print(f"Skipped {len(skipped)} protein(s): {skipped}")
    return pdf_path
