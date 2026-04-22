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

from __future__ import annotations

import logging
from collections import Counter

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
    Patsy raises if any category cell is array-like with ``ndim > 1`` or remains non-scalar.
    """
    if x is None:
        return np.nan
    if isinstance(x, (float, np.floating)) and pd.isna(x):
        return np.nan
    if isinstance(x, str):
        return x
    if isinstance(x, (bytes, np.str_)):
        return x.decode() if isinstance(x, bytes) else str(x)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x)
    if isinstance(x, pd.Series):
        return np.nan if x.empty else _to_python_scalar(x.iloc[0])
    if isinstance(x, np.ndarray):
        if x.size == 0:
            return np.nan
        if x.ndim > 1:
            return _to_python_scalar(x.ravel()[0])
        if x.dtype == object:
            return _to_python_scalar(x.item()) if x.shape == () else _to_python_scalar(x.flat[0])
        out = x.flat[0]
        return _to_python_scalar(out) if isinstance(out, np.ndarray) else out.item()
    if isinstance(x, (list, tuple)):
        return np.nan if len(x) == 0 else _to_python_scalar(x[0])
    return x


def _as_1d_series(df: pd.DataFrame, col: str) -> pd.Series:
    """Return ``df[col]`` as a single Series (first duplicate column if needed)."""
    if col not in df.columns:
        raise KeyError(f"Missing column {col!r}. Found: {list(df.columns)}")
    obj = df[col]
    if isinstance(obj, pd.DataFrame):
        return obj.iloc[:, 0].copy()
    return obj.copy()


def _patsy_scalar_categorical(series: pd.Series) -> pd.Series:
    return series.map(_to_python_scalar)


def _formula_clean_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Drop duplicate / MultiIndex column names and reset index so patsy sees a flat table."""
    out = df.copy()
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [
            "_".join(str(p) for p in tup if str(p) != "")
            if isinstance(tup, tuple)
            else tup
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
    Keep rows where library dot product is **strictly greater** than *threshold*.

    Args:
        df: Skyline report DataFrame.
        threshold: Minimum Library Dot Product (exclusive below; rows with
            value equal to *threshold* are removed). Default ``0.8``.
        col: Column name for library dot product.

    Raises:
        KeyError: If *col* is missing from *df*.
    """
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame.")

    df = df.loc[df[col] > threshold].copy()
    df = df[df['Normalized Area'].notna()]

    pivot_df = df.pivot_table(
        index=['Replicate', 'Protein Name', 'Peptide'],
        columns='Isotope Label Type',
        values='Normalized Area',
        aggfunc='first').reset_index()
    return pivot_df


def filter_peptide_counts(
    peptide_counts_df: pd.DataFrame,
    light_cutoff: int = 700,
    heavy_cutoff: int = 700,
) -> pd.DataFrame:
    """Filter the peptide counts DataFrame by light and heavy count thresholds."""
    peptide_counts_df = peptide_counts_df.loc[
        (peptide_counts_df['light_count'] > light_cutoff) &
        (peptide_counts_df['heavy_count'] > heavy_cutoff)
    ]
    return peptide_counts_df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Peptide Detection Summary
# ---------------------------------------------------------------------------


def summarise_peptide_counts(peptide_counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise the peptide counts in the DataFrame.

    Counts the number of non-missing heavy/light measurements per peptide.

    Expects input DataFrame to have columns:
      - 'Protein Name'
      - 'Peptide'
      - one column for heavy intensity (e.g., 'Heavy' or 'heavy'), and one for light (e.g., 'Light' or 'light').

    Returns:
        DataFrame with ['Protein Name', 'Peptide', 'heavy_count', 'light_count']
    """
    possible_heavy = [col for col in peptide_counts_df.columns if col.lower() == 'heavy']
    possible_light = [col for col in peptide_counts_df.columns if col.lower() == 'light']

    if not possible_heavy or not possible_light:
        raise KeyError(
            "Input DataFrame must contain columns 'Heavy' and 'Light' (case-insensitive)"
        )
    heavy_col = possible_heavy[0]
    light_col = possible_light[0]

    pept_sum = (
        peptide_counts_df.groupby(['Protein Name', 'Peptide'])
        .agg(
            heavy_count=(heavy_col, lambda x: x.notna().sum()),
            light_count=(light_col, lambda x: x.notna().sum()),
        )
        .reset_index()
    )
    return pept_sum


def report_peptide_protein_summary(peptide_counts: pd.DataFrame):
    """
    Report summary statistics on peptide and protein detection.

    Args:
        peptide_counts: DataFrame with summarised peptide counts, must include
            'Peptide', 'Protein Name', 'heavy_count', and 'light_count'.

    Returns:
        tuple: ``(summary_dict, peptide_list)``
    """
    expected_columns = ['Protein Name', 'Peptide', 'heavy_count', 'light_count']
    if not all(col in peptide_counts.columns for col in expected_columns):
        raise ValueError(f"Expected columns {expected_columns} not found in peptide_counts")

    num_unique_peptides = peptide_counts['Peptide'].nunique()
    num_unique_proteins = peptide_counts['Protein Name'].nunique()

    summary_dict = {
        'num_unique_peptides': num_unique_peptides,
        'num_unique_proteins': num_unique_proteins,
    }

    peptide_list = list(pd.unique(peptide_counts['Peptide']))

    logger.info("Number of unique peptides: %d", num_unique_peptides)
    logger.info("Number of unique proteins: %d", num_unique_proteins)
    logger.info("Selected peptides: %s", peptide_list)

    return summary_dict, peptide_list


# ---------------------------------------------------------------------------
# CV Analysis
# ---------------------------------------------------------------------------


def calculate_intra_plate_cv(
    pool_data: pd.DataFrame,
    col_name: str = 'characteristics[Plate]',
) -> pd.DataFrame:
    """
    Calculate intra-plate CV per peptide grouped by the given column.

    Args:
        pool_data: DataFrame containing at least [col_name, 'Peptide Sequence', 'RatioLightToHeavy'].
        col_name: Column to group by (default: 'characteristics[Plate]').

    Returns:
        DataFrame summarizing mean, std, and intra_plate_cv per group/peptide.
    """
    if col_name not in pool_data.columns:
        raise KeyError(f"Column '{col_name}' not found in pool_data")

    pool_data = pool_data[[col_name, 'Peptide Sequence', 'RatioLightToHeavy']].drop_duplicates()

    peptide_plate_stats = (
        pool_data.groupby([col_name, 'Peptide Sequence'])['RatioLightToHeavy']
        .agg(['mean', 'std'])
        .reset_index()
    )
    peptide_plate_stats['intra_plate_cv'] = peptide_plate_stats['std'] / peptide_plate_stats['mean']
    return peptide_plate_stats


def calculate_inter_plate_cv(peptide_plate_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate inter-plate CV from intra-plate stats.

    Args:
        peptide_plate_stats: Output of :func:`calculate_intra_plate_cv`.

    Returns:
        DataFrame with ['Peptide Sequence', 'grand_mean', 'between_plate_sd', 'inter_plate_cv'],
        sorted by inter_plate_cv ascending.
    """
    peptide_means = (
        peptide_plate_stats.groupby('Peptide Sequence')['mean']
        .agg(['mean', 'std'])
        .rename(columns={'mean': 'grand_mean', 'std': 'between_plate_sd'})
        .reset_index()
    )
    peptide_means['inter_plate_cv'] = peptide_means['between_plate_sd'] / peptide_means['grand_mean']
    return peptide_means.sort_values('inter_plate_cv').reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plate Normalization
# ---------------------------------------------------------------------------


def get_lowest_cv_peptides(interplate_cv_df: pd.DataFrame, cv_percentile: float):
    """
    Return peptide sequences with inter-plate CV at or below the specified percentile.

    Args:
        interplate_cv_df: Output of :func:`calculate_inter_plate_cv`.
        cv_percentile: Percentile threshold as a percentage (e.g. ``10`` for 10%).
    """
    cv_threshold = cv_percentile / 100.0
    selected_peptides = interplate_cv_df.loc[
        interplate_cv_df['inter_plate_cv'] <= cv_threshold,
        'Peptide Sequence'
    ].unique()
    return selected_peptides


def extract_top_percentile(
    df: pd.DataFrame,
    column: str,
    percentile: float = 0.1,
    id_col: str = 'Peptide Sequence',
    source_df: pd.DataFrame = None,
    source_col: str = None,
):
    """
    Return IDs from ``id_col`` whose ``column`` value is at or below the given
    quantile, and optionally filter a source DataFrame to those IDs.

    Args:
        df: DataFrame with peptide statistics (e.g. output of :func:`calculate_inter_plate_cv`).
        column: Column to compute the quantile threshold from (e.g. ``'inter_plate_cv'``).
        percentile: Quantile threshold as a fraction in [0, 1] (e.g. ``0.1`` for the lowest 10%).
        id_col: Column whose unique values to extract (default ``'Peptide Sequence'``).
        source_df: If given, filter this DataFrame to rows whose ``source_col`` is in the returned ID list.
        source_col: Column of ``source_df`` to match against IDs (defaults to ``id_col``).

    Returns:
        tuple: ``(id_array, filtered_df)``; ``filtered_df`` is ``None`` if ``source_df`` is not provided.
    """
    threshold = df[column].quantile(percentile)
    id_list = df.loc[df[column] <= threshold, id_col].unique()
    if source_df is not None:
        if source_col is None:
            source_col = id_col
        filtered = source_df[source_df[source_col].isin(id_list)].reset_index(drop=True)
        return id_list, filtered
    return id_list, None


def plate_peptide_anova(
    selected_norm_peptides: pd.DataFrame,
    plate_col: str = "characteristics[Plate]",
    log_transform: bool = False,
):
    """
    Fit a two-way ANOVA on ratio ~ plate + peptide.

    Args:
        selected_norm_peptides: DataFrame with ``RatioLightToHeavy`` and peptide / plate columns.
        plate_col: Source column for plate (e.g. ``'characteristics[Plate]'`` or ``'Plate'``).
            Renamed internally to ``Plate`` for the formula ``C(Plate)``.
        log_transform: If True, model ``log(RatioLightToHeavy)`` (only rows with ratio > 0).

    Returns:
        ``model`` (statsmodels OLS result), ``anova_res`` (Type II ANOVA table).

    See :func:`get_plate_conversion_factors` for plate correction factors.
    """
    df = _formula_clean_frame(selected_norm_peptides.copy())

    if plate_col in df.columns and plate_col != "Plate":
        if "Plate" in df.columns:
            df = df.drop(columns=["Plate"])
        df = df.rename(columns={plate_col: "Plate"})
    elif "characteristics[Plate]" in df.columns and "Plate" not in df.columns:
        df = df.rename(columns={"characteristics[Plate]": "Plate"})
    elif "characteristics[Plate]" in df.columns and "Plate" in df.columns:
        df = df.drop(columns=["characteristics[Plate]"])

    if "Peptide Sequence" in df.columns and "Peptide" not in df.columns:
        df = df.rename(columns={"Peptide Sequence": "Peptide"})
    elif "Peptide Sequence" in df.columns and "Peptide" in df.columns:
        df = df.drop(columns=["Peptide Sequence"])

    df = _formula_clean_frame(df)

    missing = [c for c in ("RatioLightToHeavy", "Plate", "Peptide") if c not in df.columns]
    if missing:
        raise KeyError(
            f"plate_peptide_anova missing columns {missing!r}. Have: {list(df.columns)!r}"
        )

    ratio = pd.to_numeric(
        _patsy_scalar_categorical(_as_1d_series(df, "RatioLightToHeavy")),
        errors="coerce",
    )
    plate = _patsy_scalar_categorical(_as_1d_series(df, "Plate"))
    peptide = _patsy_scalar_categorical(_as_1d_series(df, "Peptide"))
    df = pd.DataFrame(
        {"RatioLightToHeavy": ratio, "Plate": plate, "Peptide": peptide}
    ).dropna(subset=["RatioLightToHeavy", "Plate", "Peptide"]).reset_index(drop=True)
    df["Plate"] = df["Plate"].astype(str)
    df["Peptide"] = df["Peptide"].astype(str)

    if df.empty:
        raise ValueError(
            "plate_peptide_anova: no rows left after dropping missing Plate/Peptide/ratio."
        )

    if log_transform:
        df = df.loc[df["RatioLightToHeavy"] > 0].copy().reset_index(drop=True)
        if df.empty:
            raise ValueError(
                "plate_peptide_anova: no positive RatioLightToHeavy rows for log_transform."
            )
        df["log_ratio"] = np.log(df["RatioLightToHeavy"])
        df["response"] = df["log_ratio"]
    else:
        df["response"] = df["RatioLightToHeavy"]

    model = smf.ols("response ~ C(Plate) + C(Peptide)", data=df).fit()
    anova_res = anova_lm(model, typ=2)
    logger.info("OLS model summary:\n%s", model.summary())
    logger.info("Type II ANOVA table:\n%s", anova_res)

    def _p_from_anova(label_substr: str):
        rows = [i for i in anova_res.index if label_substr in str(i).lower()]
        if not rows:
            return None
        pr_col = next(
            (c for c in anova_res.columns if "PR" in c or c.startswith("P")),
            None,
        )
        return anova_res.loc[rows[0], pr_col] if pr_col else None

    plate_p = _p_from_anova("plate")
    peptide_p = _p_from_anova("peptide")

    if plate_p is not None and pd.notna(plate_p):
        logger.info("Plate effect p-value (Type II ANOVA): %.4g", plate_p)
        if plate_p < 0.05:
            logger.info("Plate effect is significant — batch correction recommended.")
        else:
            logger.info("Plate effect is not significant — no batch correction needed.")
    else:
        logger.warning("Plate effect p-value could not be read from the ANOVA table.")

    if peptide_p is not None and pd.notna(peptide_p):
        logger.info("Peptide effect p-value (Type II ANOVA): %.4g", peptide_p)
        if peptide_p < 0.05:
            logger.info("Peptide effect is significant.")
        else:
            logger.info("Peptide effect is not significant.")
    else:
        logger.warning("Peptide effect p-value could not be read from the ANOVA table.")

    return model, anova_res


def get_plate_conversion_factors(
    df: pd.DataFrame,
    col_plate: str = "Plate",
    ratio_col: str = "RatioLightToHeavy",
    log_transform: bool = False,
):
    """
    Calculate plate conversion factors for normalisation.

    Returns:
        conv_factors_df (pd.DataFrame): DataFrame of plate median and conversion factor.
        conversion_factors (dict): Mapping plate -> multiplicative correction factor.
        model: Fitted OLS model from statsmodels.
    """
    cols = ['Peptide Sequence', 'Replicate', col_plate, ratio_col]
    missing_cols = [col for col in cols if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing required columns in df: {missing_cols}")

    d = df[cols].drop_duplicates().copy()
    d = d.assign(ratio_fit=pd.to_numeric(d[ratio_col], errors='coerce'))
    if log_transform:
        d['ratio_fit'] = np.log(d['ratio_fit'])

    model = smf.ols(
        f'ratio_fit ~ C(Q("{col_plate}")) + C(Q("Peptide Sequence"))',
        data=d,
    ).fit()

    platemed = (
        d.groupby(col_plate, dropna=True, observed=True)['ratio_fit']
        .median()
        .reset_index()
        .rename(columns={'ratio_fit': 'plate_median'})
    )
    global_median = d['ratio_fit'].median()
    platemed['correction_factor'] = global_median / platemed['plate_median']

    conversion_factors = {
        str(row[col_plate]): row['correction_factor']
        for _, row in platemed.iterrows()
    }
    return platemed, conversion_factors, model


def adjust_ratio_by_plate(df: pd.DataFrame, conversion_factors: dict) -> pd.DataFrame:
    """
    Normalise RatioLightToHeavy by the plate-specific conversion factor.

    Args:
        df: DataFrame with columns 'RatioLightToHeavy' and 'Plate'.
        conversion_factors: Dict mapping plate label (str) to multiplicative factor.

    Returns:
        Copy of *df* with an added 'RatioLightToHeavy_adj' column.
    """
    def get_factor(plate):
        key = str(plate)
        if key not in conversion_factors:
            try:
                key = str(int(float(plate)))
            except (ValueError, TypeError):
                pass
        if key not in conversion_factors:
            raise KeyError(
                f"Plate {plate!r} not found in conversion_factors. "
                f"Available plates: {list(conversion_factors.keys())}"
            )
        return conversion_factors[key]

    df = df.copy()
    df['RatioLightToHeavy_adj'] = [
        r * get_factor(p) for r, p in zip(df['RatioLightToHeavy'], df['Plate'])
    ]
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
    Compute the percent coefficient of variation (%CV) for a numeric column,
    grouped by a categorical column (default: per precursor across replicates).

    Args:
        df: Skyline report DataFrame.
        value_col: Column to compute CV on. Defaults to ``'RatioLightToHeavy'``.
        group_by: Column used to form groups. Defaults to ``'Precursor'``.

    Returns:
        DataFrame with columns ``[group_by, 'mean', 'std', 'cv_pct', 'n']``,
        one row per group, sorted by ``cv_pct`` ascending.
    """
    for col in (value_col, group_by):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    grouped = df.groupby(group_by)[value_col]
    result = grouped.agg(['mean', 'std', 'count']).reset_index()
    result = result.rename(columns={'count': 'n'})
    result["cv_pct"] = (result["std"] / result["mean"]).abs() * 100
    result.loc[result["mean"].abs() < 1e-12, "cv_pct"] = float("nan")
    result = result[[group_by, 'mean', 'std', 'cv_pct', 'n']]
    return result.sort_values("cv_pct", na_position='last').reset_index(drop=True)


def flag_missing_values(
    df: pd.DataFrame,
    value_col: str = "RatioLightToHeavy",
    precursor_col: str = "Precursor",
    replicate_col: str = "Replicate",
) -> pd.DataFrame:
    """
    Identify precursors with missing quantification values across replicates.

    Args:
        df: Skyline report DataFrame.
        value_col: Column holding the quantitative values.
        precursor_col: Column identifying precursors.
        replicate_col: Column identifying replicates.

    Returns:
        DataFrame with one row per precursor, plus columns
        ``n_detected``, ``n_total``, ``detection_rate``.
    """
    for col in (value_col, precursor_col, replicate_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    pivot = df.pivot_table(
        index=precursor_col,
        columns=replicate_col,
        values=value_col,
        aggfunc="first",
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
    Summarise library and ratio dot-product scores per precursor and flag
    transitions that fall below quality thresholds.

    Args:
        df: Skyline report DataFrame.
        lib_dot_col: Column for library dot product.
        ratio_dot_col: Column for ratio dot product.
        lib_threshold: Minimum acceptable library dot product.
        ratio_threshold: Minimum acceptable ratio dot product.

    Returns:
        DataFrame with ``lib_dot_pass``, ``ratio_dot_pass``, ``both_pass`` columns.
    """
    for col in (lib_dot_col, ratio_dot_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    missing_cols = [c for c in ("Precursor", "Replicate") if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Required columns for dot_product_summary missing: {missing_cols}")

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
    Calculate the absolute and relative deviation between observed and
    predicted retention times for each precursor–replicate pair.

    Args:
        df: Skyline report DataFrame.
        observed_col: Column with observed retention time.
        predicted_col: Column with predicted (iRT-based) retention time.

    Returns:
        DataFrame sorted by ``abs_rt_dev`` descending, with columns
        ``rt_observed``, ``rt_predicted``, ``rt_dev``, ``abs_rt_dev``.
    """
    for col in (observed_col, predicted_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    base_cols = [c for c in ("Precursor", "Replicate") if c in df.columns]
    base_cols += [observed_col, predicted_col]

    result = df[base_cols].copy().rename(
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

    Convenience wrapper around :func:`compute_cv`, :func:`flag_missing_values`,
    :func:`dot_product_summary`, and :func:`retention_time_deviation`.

    Args:
        df: Skyline report DataFrame.
        cv_col: Column used for CV calculation.
        lib_threshold: Minimum acceptable library dot product.
        ratio_threshold: Minimum acceptable ratio dot product.
        rt_dev_threshold: Maximum acceptable absolute RT deviation in minutes.

    Returns:
        dict with keys: ``cv``, ``missing``, ``dot_products``, ``rt_deviation``,
        ``n_precursors``, ``n_replicates``, ``pct_dot_pass``, ``pct_rt_within``, ``median_cv_pct``.
    """
    cv_df = compute_cv(df, value_col=cv_col)
    missing_df = flag_missing_values(df, value_col=cv_col)
    dot_df = dot_product_summary(df, lib_threshold=lib_threshold, ratio_threshold=ratio_threshold)
    rt_df = retention_time_deviation(df)

    n_precursors = df["Precursor"].nunique() if "Precursor" in df.columns else None
    n_replicates = df["Replicate"].nunique() if "Replicate" in df.columns else None
    pct_dot_pass = dot_df["both_pass"].mean() * 100 if len(dot_df) > 0 else float("nan")
    pct_rt_within = (rt_df["abs_rt_dev"] <= rt_dev_threshold).mean() * 100 if len(rt_df) > 0 else float("nan")
    median_cv = cv_df["cv_pct"].median() if len(cv_df) > 0 else float("nan")

    return {
        "cv": cv_df,
        "missing": missing_df,
        "dot_products": dot_df,
        "rt_deviation": rt_df,
        "n_precursors": n_precursors,
        "n_replicates": n_replicates,
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
        peptide_counts: DataFrame with columns 'heavy_count' and 'light_count'.
    """
    x = peptide_counts['heavy_count'].values
    y = peptide_counts['light_count'].values

    xy = list(zip(x, y))
    counts = Counter(xy)
    point_count = np.array([counts[(hx, ly)] for hx, ly in xy])
    sort_idx = np.argsort(point_count)

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(
        x[sort_idx], y[sort_idx],
        c=point_count[sort_idx],
        cmap='viridis',
        alpha=0.7,
        s=60,
        edgecolors='k',
        linewidth=0.5,
    )
    plt.xlabel("Heavy Count")
    plt.ylabel("Light Count")
    plt.title("Scatter plot of Heavy vs. Light Peptide Counts\n(colored by number of peptides at each point)")
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=1)
    plt.grid(True)
    plt.colorbar(sc, label='# Peptides at Point')
    plt.tight_layout()
    plt.show()


def plot_peptide_counts(df: pd.DataFrame) -> None:
    """Histogram (+ KDE) of heavy peptide counts."""
    plt.figure(figsize=(10, 6))
    sns.histplot(df['heavy_count'], bins=20, kde=True)
    plt.title("Distribution of Heavy Peptide Counts")
    plt.xlabel("Heavy Peptide Counts")
    plt.show()


def plot_pool_boxplot(
    pool_df: pd.DataFrame,
    col_ratio: str = 'RatioLightToHeavy',
    col_plate: str = 'characteristics[Plate]',
) -> None:
    """
    Boxplot of log(RatioLightToHeavy) by Replicate, colored by Plate.

    Args:
        pool_df: DataFrame with columns 'Replicate', col_ratio, and col_plate.
        col_ratio: Column for the light-to-heavy ratio.
        col_plate: Column for plate labels.
    """
    pool_df = pool_df[['Replicate', col_ratio, col_plate]].drop_duplicates()

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Replicate', y=col_ratio, data=pool_df, hue=col_plate)
    plt.title(f'Boxplot of {col_ratio} by {col_plate} (colored by {col_plate})')
    plt.xlabel(col_plate)
    plt.ylabel(col_ratio)
    plt.yscale('log')
    plt.xticks([], [])
    plt.legend(title=col_plate, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()


def plot_pool_heatmap(pool_data: pd.DataFrame, aggfunc: str = "mean") -> None:
    """
    Heatmap of log(RatioLightToHeavy) for Pool samples.

    Peptides (rows) and Replicates (columns) are both ordered by their mean log-ratio.
    Duplicate peptide/replicate combinations are aggregated with *aggfunc*.

    Args:
        pool_data: DataFrame with columns ``Peptide Sequence``, ``Replicate``, ``RatioLightToHeavy``.
        aggfunc: Aggregation function passed to ``pivot_table`` (e.g. ``"mean"``, ``"median"``).
    """
    required = ["Peptide Sequence", "Replicate", "RatioLightToHeavy"]
    missing = [c for c in required if c not in pool_data.columns]
    if missing:
        raise KeyError(f"pool_data missing columns {missing}. Found: {list(pool_data.columns)}")

    pivot = pool_data.pivot_table(
        index="Peptide Sequence",
        columns="Replicate",
        values="RatioLightToHeavy",
        aggfunc=aggfunc,
    )
    heatmap_data = np.log(pivot)
    col_order = heatmap_data.mean(axis=0).sort_values().index
    row_order = heatmap_data.mean(axis=1).sort_values().index
    heatmap_data = heatmap_data.loc[row_order, col_order]

    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, cmap='viridis')
    plt.title('Heatmap of log(RatioLightToHeavy) for Pool')
    plt.xlabel('')
    plt.ylabel('Peptide Sequence')
    plt.xticks([], [])
    plt.show()


def plot_intra_plate_cv_stats(
    peptide_plate_stats: pd.DataFrame,
    col_name: str = 'characteristics[Plate]',
) -> None:
    """
    Boxplot of intra-plate CV per group/peptide.

    Args:
        peptide_plate_stats: Output of :func:`calculate_intra_plate_cv`.
        col_name: The group column used in the stats DataFrame.
    """
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=col_name, y='intra_plate_cv', data=peptide_plate_stats)
    plt.title('Boxplot of Intra Plate CV for Pool')
    plt.xlabel(col_name.replace('characteristics[', '').replace(']', '').capitalize())
    plt.ylabel('Intra Plate CV')
    plt.show()


def plot_inter_plate_cv_kde(interplate_cv: pd.DataFrame) -> plt.Figure:
    """
    KDE of inter-plate CV across peptides, with a vertical line at the median.

    Args:
        interplate_cv: Output of :func:`calculate_inter_plate_cv`.
            Must have an ``inter_plate_cv`` column.

    Returns:
        The matplotlib Figure object.
    """
    if 'inter_plate_cv' not in interplate_cv.columns:
        raise KeyError(
            f"Expected column 'inter_plate_cv'. Found: {list(interplate_cv.columns)}"
        )
    fig = plt.figure(figsize=(12, 6))
    sns.kdeplot(interplate_cv['inter_plate_cv'].dropna(), fill=True)
    median_cv = interplate_cv['inter_plate_cv'].median()
    plt.axvline(median_cv, color='red', linestyle='--', label=f'Median = {median_cv:.2f}')
    plt.title('KDE Plot of Inter-Plate CV Across Peptides')
    plt.xlabel('Inter-Plate CV')
    plt.ylabel('Density')
    plt.legend()
    plt.show()
    return fig


def plot_cumulative_peptide_count_by_cv(peptide_means: pd.DataFrame) -> None:
    """
    Cumulative number of peptides as a function of sorted inter-plate CV.

    Args:
        peptide_means: DataFrame with columns ``inter_plate_cv`` and ``Peptide Sequence``.
    """
    cv_sorted = (
        peptide_means[['inter_plate_cv', 'Peptide Sequence']]
        .sort_values('inter_plate_cv')
        .reset_index(drop=True)
    )
    cv_sorted['cumulative_count'] = range(1, len(cv_sorted) + 1)

    print(f"Total number of peptides: {len(cv_sorted)}")
    for thresh in [0.10, 0.20]:
        count_below = (cv_sorted['inter_plate_cv'] < thresh).sum()
        print(f"Number of peptides with inter-plate CV < {int(thresh*100)}%: {count_below}")

    plt.figure(figsize=(10, 6))
    plt.plot(cv_sorted['inter_plate_cv'], cv_sorted['cumulative_count'], marker='o', linestyle='-')
    plt.xlabel('Inter-Plate CV')
    plt.ylabel('Cumulative Number of Peptides')
    plt.title('Cumulative Peptide Count by Inter-Plate CV')

    for cv_mark in [0.1 * i for i in range(1, 11)]:
        mask = cv_sorted['inter_plate_cv'] >= cv_mark
        if mask.any():
            idx = mask.idxmax()
            x = cv_sorted.at[idx, 'inter_plate_cv']
            y = cv_sorted.at[idx, 'cumulative_count']
            plt.axvline(x, color='gray', linestyle='--', linewidth=0.8)
            plt.text(x, y, f"{int(y)} peptides\n{cv_mark:.1f} CV", va='bottom', ha='left', fontsize=9, color='blue')

    plt.tight_layout()
    plt.show()


def plot_logratio_by_plate_boxplot(df: pd.DataFrame) -> None:
    """
    Boxplot of log(RatioLightToHeavy) by Replicate, colored by Plate.

    Args:
        df: DataFrame with columns ['Replicate', 'log_ratio', 'Plate'].
    """
    df_sorted = df.sort_values('Plate')
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Replicate', y='log_ratio', data=df_sorted, hue='Plate', dodge=False)
    plt.title('Boxplot of log(RatioLightToHeavy) by Replicate (colored by Plate)')
    plt.xlabel('')
    plt.ylabel('log(RatioLightToHeavy)')
    plt.xticks([], [])
    plt.legend(title='Plate', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
