"""
PRM (Parallel Reaction Monitoring) analysis and plotting functions for Skyline exports.

This module provides functions to assess quantification quality from Skyline
PRM reports. It expects a DataFrame produced by ``ImportFile.import_skyline_file``
and works with the following key columns:

- ``Precursor``               – precursor ion string (sequence + charge)
- ``Replicate``               – replicate / sample label
- ``File Name``               – raw data file name
- ``Peptide Sequence``        – stripped peptide sequence
- ``Peptide``                 – modified peptide string
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
    return df.pivot_table(
        index=["Replicate", "Protein Name", "Peptide"],
        columns="Isotope Label Type",
        values="Normalized Area",
        aggfunc="first",
    ).reset_index()


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
    col_name: str = "characteristics[Plate]",
) -> pd.DataFrame:
    """
    Calculate intra-plate CV per peptide.

    Args:
        pool_data: DataFrame with at least ``[col_name, 'Peptide Sequence', 'RatioLightToHeavy']``.
        col_name: Column to group by. Default ``'characteristics[Plate]'``.

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
    plate_col: str = "characteristics[Plate]",
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
    elif "characteristics[Plate]" in df.columns:
        df = df.drop(columns=["Plate"], errors="ignore").rename(
            columns={"characteristics[Plate]": "Plate"}
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
    sns.boxplot(
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
    plt.show()
    return plt.gcf()


def adjust_ratio_by_plate(
    df: pd.DataFrame,
    conversion_factors: dict,
    col_match: str = 'characteristics[Plate]'
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
    plate_factors = df[col_match].astype(str).map(conversion_factors)
    if plate_factors.isnull().any():
        missing = df.loc[plate_factors.isnull(), col_match].unique()
        raise KeyError(f"Some plate values have no conversion factor: {missing}")
    df["RatioLightToHeavy"] = df["RatioLightToHeavy"] / plate_factors.values
    return df

# ---------------------------------------------------------------------------
# General QC
# ---------------------------------------------------------------------------


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
    """Histogram (+ KDE) of ``Library Dot Product``."""
    col = "Library Dot Product"
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame.")
    plt.figure(figsize=(10, 6))
    sns.histplot(df[col], bins=20, kde=True)
    plt.title("Distribution of Library Dot Product")
    plt.xlabel("Library Dot Product")
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
    col_plate: str = "characteristics[Plate]",
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
    col_name: str = "characteristics[Plate]",
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
    plt.xlabel(col_name.replace("characteristics[", "").replace("]", "").capitalize())
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
