"""
PRM (Parallel Reaction Monitoring) analysis functions for Skyline exports.

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

Plotting functions are in :mod:`skyline_qc.prm_plots`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm


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


def filter_peptide_counts(peptide_counts_df: pd.DataFrame, light_cutoff: int = 700, heavy_cutoff: int = 700) -> pd.DataFrame:
    """
    Filter the peptide counts in the DataFrame.
    """
    peptide_counts_df = peptide_counts_df.loc[(peptide_counts_df['light_count'] > light_cutoff) & (peptide_counts_df['heavy_count'] > heavy_cutoff)]
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


def report_peptide_protein_summary(peptide_counts):
    """
    Report summary statistics on peptide and protein detection.

    Args:
        peptide_counts: DataFrame containing summarised peptide counts, must include 'Peptide', 'heavy_count', and 'light_count'.

    Returns:
        summary_dict: Dictionary with summary statistics.
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

    print(f"Number of unique peptides: {num_unique_peptides}")
    print(f"Number of unique proteins: {num_unique_proteins}")

    return summary_dict


# ---------------------------------------------------------------------------
# CV Analysis
# ---------------------------------------------------------------------------


def calculate_intra_plate_cv(pool_data: pd.DataFrame, col_name: str = 'characteristics[Plate]') -> pd.DataFrame:
    """
    Calculate intra-plate CV per peptide grouped by the given column.

    Args:
        pool_data (pd.DataFrame): DataFrame containing at least [col_name, 'Peptide Sequence', 'RatioLightToHeavy'] columns.
        col_name (str): The column name to group by (default: 'characteristics[Plate]').

    Returns:
        pd.DataFrame: DataFrame summarizing mean, std, and intra_plate_cv per group/peptide.
    """
    peptide_plate_stats = (
        pool_data.groupby([col_name, 'Peptide Sequence'])['RatioLightToHeavy']
        .agg(['mean', 'std'])
        .reset_index()
    )
    peptide_plate_stats['intra_plate_cv'] = peptide_plate_stats['std'] / peptide_plate_stats['mean']
    return peptide_plate_stats


def get_peptide_means(peptide_plate_stats):
    """
    Given a DataFrame of peptide_plate_stats (output of calculate_intra_plate_cv),
    return a DataFrame of peptide_means with columns:
    ['Peptide Sequence', 'grand_mean', 'between_plate_sd', 'inter_plate_cv']
    """
    peptide_means = (
        peptide_plate_stats.groupby('Peptide Sequence')['mean']
        .agg(['mean', 'std'])
        .rename(columns={'mean': 'grand_mean', 'std': 'between_plate_sd'})
        .reset_index()
    )
    peptide_means['inter_plate_cv'] = peptide_means['between_plate_sd'] / peptide_means['grand_mean']

    return peptide_means.sort_values('inter_plate_cv').reset_index(drop=True)


def get_peptides_below_cv_percentile(peptide_means, percentile):
    """
    Return a list of peptide sequences with inter-plate CV
    at or below the specified percentile.

    Args:
        peptide_means (pd.DataFrame): DataFrame with 'inter_plate_cv' and 'Peptide Sequence' columns.
        percentile (float): Percentile threshold (between 0 and 100).

    Returns:
        np.ndarray: Array of peptide sequences meeting the criterion.
    """
    cv_threshold = peptide_means['inter_plate_cv'].quantile(percentile / 100.0)
    selected_peptides = peptide_means.loc[
        peptide_means['inter_plate_cv'] <= cv_threshold,
        'Peptide Sequence'
    ].unique()
    return selected_peptides


# ---------------------------------------------------------------------------
# Plate Normalization
# ---------------------------------------------------------------------------


def plate_peptide_anova(selected_norm_peptides):
    """
    Fits a two-way ANOVA linear model to decompose variation in ratio measurements
    into contributions from plate and peptide effects.

    Args:
        selected_norm_peptides (pd.DataFrame): DataFrame with
            'characteristics[Plate]', 'Peptide Sequence', and 'RatioLightToHeavy' columns.

    Returns:
        model: The fitted statsmodels OLS linear model.
        anova_res: The ANOVA decomposition results (as a DataFrame).
    """
    df = selected_norm_peptides.rename(
        columns={
            'characteristics[Plate]': 'Plate',
            'Peptide Sequence': 'Peptide',
        }
    ).copy()

    df['RatioLightToHeavy'] = pd.to_numeric(df['RatioLightToHeavy'], errors='coerce')
    df = df.dropna(subset=['RatioLightToHeavy', 'Plate', 'Peptide'])

    model = smf.ols('RatioLightToHeavy ~ C(Plate) + C(Peptide)', data=df).fit()
    print(model.summary())

    anova_res = anova_lm(model, typ=2)
    print(anova_res)

    return model, anova_res


def get_plate_conversion_factors(df):
    """
    Fit a model: log_ratio ~ C(Plate) to extract plate conversion factors.
    This gives the correction (as multiplicative factor) for each plate to equalize across plates.

    Args:
        df (pd.DataFrame): DataFrame containing at least 'log_ratio' and 'Plate' columns.

    Returns:
        conv_factors_df (pd.DataFrame): DataFrame with columns ['Plate', 'conversion_factor'].
        conversion_factors (dict): Dictionary mapping Plate value to its conversion factor.
        model (statsmodels.regression.linear_model.RegressionResultsWrapper): The fitted model.
    """
    m = smf.ols('log_ratio ~ C(Plate)', data=df).fit()

    plate_effects = m.params.filter(like='C(Plate)')
    reference_plate = df['Plate'].unique()
    if isinstance(df['Plate'].iloc[0], str):
        reference_plate = sorted(reference_plate, key=lambda x: str(x))
    else:
        reference_plate = sorted(reference_plate)
    reference_plate = reference_plate[0]

    conversion_factors = {}
    conversion_factors[reference_plate] = 1.0
    for term, coef in plate_effects.items():
        plate_number = term.replace('C(Plate)[T.', '').replace(']', '')
        conversion_factors[plate_number] = np.exp(coef)

    conv_factors_df = pd.DataFrame(
        list(conversion_factors.items()),
        columns=['Plate', 'conversion_factor']
    )
    print("Conversion factors for each plate (to equalize them):")
    print(conv_factors_df)

    return conv_factors_df, conversion_factors, m


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
        df: Skyline report DataFrame (from ``ImportFile.import_skyline_file``).
        value_col: Column to compute CV on. Defaults to ``'RatioLightToHeavy'``.
        group_by: Column used to form groups. Defaults to ``'Precursor'``.

    Returns:
        pd.DataFrame with columns ``[group_by, 'mean', 'std', 'cv_pct', 'n']``,
        one row per group, sorted by ``cv_pct`` ascending.

    Raises:
        KeyError: If ``value_col`` or ``group_by`` are not found in *df*.
    """
    for col in (value_col, group_by):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    result = (
        df.groupby(group_by)[value_col]
        .agg(
            mean="mean",
            std="std",
            n="count",
        )
        .reset_index()
    )
    result["cv_pct"] = (result["std"] / result["mean"]).abs() * 100
    return result.sort_values("cv_pct").reset_index(drop=True)


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
        value_col: Column holding the quantitative values. Defaults to ``'RatioLightToHeavy'``.
        precursor_col: Column identifying precursors. Defaults to ``'Precursor'``.
        replicate_col: Column identifying replicates. Defaults to ``'Replicate'``.

    Returns:
        pd.DataFrame with one row per precursor, plus columns:
        ``n_detected``, ``n_total``, ``detection_rate``.

    Raises:
        KeyError: If any of the required columns are absent from *df*.
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
    n_total = pivot.shape[1]

    pivot["n_detected"] = detected.sum(axis=1)
    pivot["n_total"] = n_total
    pivot["detection_rate"] = pivot["n_detected"] / n_total

    return pivot.reset_index().sort_values("detection_rate").reset_index(drop=True)


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
        lib_dot_col: Column for library dot product. Defaults to ``'Library Dot Product'``.
        ratio_dot_col: Column for ratio dot product. Defaults to ``'Ratio Dot Product'``.
        lib_threshold: Minimum acceptable library dot product. Defaults to ``0.8``.
        ratio_threshold: Minimum acceptable ratio dot product. Defaults to ``0.8``.

    Returns:
        pd.DataFrame with ``lib_dot_pass``, ``ratio_dot_pass``, ``both_pass`` columns.

    Raises:
        KeyError: If *lib_dot_col* or *ratio_dot_col* are not found in *df*.
    """
    for col in (lib_dot_col, ratio_dot_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    cols = ["Precursor", "Replicate", lib_dot_col, ratio_dot_col]
    result = df[cols].copy()
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
        observed_col: Column with observed retention time. Defaults to ``'Peptide Retention Time'``.
        predicted_col: Column with predicted (iRT-based) retention time. Defaults to ``'Predicted Retention Time'``.

    Returns:
        pd.DataFrame sorted by ``abs_rt_dev`` descending, with columns
        ``rt_observed``, ``rt_predicted``, ``rt_dev``, ``abs_rt_dev``.

    Raises:
        KeyError: If *observed_col* or *predicted_col* are absent from *df*.
    """
    for col in (observed_col, predicted_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

    result = df[["Precursor", "Replicate", observed_col, predicted_col]].copy()
    result = result.rename(
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
        cv_col: Column used for CV calculation. Defaults to ``'RatioLightToHeavy'``.
        lib_threshold: Minimum acceptable library dot product. Defaults to ``0.8``.
        ratio_threshold: Minimum acceptable ratio dot product. Defaults to ``0.8``.
        rt_dev_threshold: Maximum acceptable absolute RT deviation in minutes. Defaults to ``2.0``.

    Returns:
        dict with keys: ``cv``, ``missing``, ``dot_products``, ``rt_deviation``,
        ``n_precursors``, ``n_replicates``, ``pct_dot_pass``, ``pct_rt_within``, ``median_cv_pct``.
    """
    cv_df = compute_cv(df, value_col=cv_col)
    missing_df = flag_missing_values(df)
    dot_df = dot_product_summary(
        df, lib_threshold=lib_threshold, ratio_threshold=ratio_threshold
    )
    rt_df = retention_time_deviation(df)

    n_precursors = df["Precursor"].nunique()
    n_replicates = df["Replicate"].nunique()
    pct_dot_pass = dot_df["both_pass"].mean() * 100
    pct_rt_within = (rt_df["abs_rt_dev"] <= rt_dev_threshold).mean() * 100
    median_cv = cv_df["cv_pct"].median()

    return {
        "cv": cv_df,
        "missing": missing_df,
        "dot_products": dot_df,
        "rt_deviation": rt_df,
        "n_precursors": n_precursors,
        "n_replicates": n_replicates,
        "pct_dot_pass": round(pct_dot_pass, 2),
        "pct_rt_within": round(pct_rt_within, 2),
        "median_cv_pct": round(median_cv, 2),
    }
