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
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# ---------------------------------------------------------------------------
# Coefficient of Variation
# ---------------------------------------------------------------------------


def compute_cv(
    df: pd.DataFrame,
    value_col: str = "RatioLightToHeavy",
    group_by: str = "Precursor",
) -> pd.DataFrame:
    """
    Compute the percent coefficient of variation (%CV) for a numeric column,
    grouped by a categorical column (default: per precursor across replicates).

    CV is calculated as ``(std / mean) * 100``.  Groups with fewer than two
    non-null observations will have ``NaN`` in the CV column.

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


# ---------------------------------------------------------------------------
# Missing-value / detection summary
# ---------------------------------------------------------------------------


def flag_missing_values(
    df: pd.DataFrame,
    value_col: str = "RatioLightToHeavy",
    precursor_col: str = "Precursor",
    replicate_col: str = "Replicate",
) -> pd.DataFrame:
    """
    Identify precursors with missing quantification values across replicates.

    A value is considered missing when it is ``NaN`` or zero.  The function
    pivots the data into a precursor × replicate matrix and annotates each
    precursor with its detection rate.

    Args:
        df: Skyline report DataFrame.
        value_col: Column holding the quantitative values. Defaults to
            ``'RatioLightToHeavy'``.
        precursor_col: Column identifying precursors. Defaults to
            ``'Precursor'``.
        replicate_col: Column identifying replicates. Defaults to
            ``'Replicate'``.

    Returns:
        pd.DataFrame with one row per precursor.  Columns include the
        replicate-level values (wide format), plus:

        - ``n_detected``     – number of replicates with a non-missing value
        - ``n_total``        – total number of replicates
        - ``detection_rate`` – fraction of replicates detected (0.0 – 1.0)

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


# ---------------------------------------------------------------------------
# Dot-product quality summary
# ---------------------------------------------------------------------------


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

    Dot products range from 0 to 1; values below the threshold indicate poor
    spectral similarity to the reference library or between light and heavy
    channels.

    Args:
        df: Skyline report DataFrame.
        lib_dot_col: Column for library dot product. Defaults to
            ``'Library Dot Product'``.
        ratio_dot_col: Column for ratio dot product. Defaults to
            ``'Ratio Dot Product'``.
        lib_threshold: Minimum acceptable library dot product. Defaults to
            ``0.8``.
        ratio_threshold: Minimum acceptable ratio dot product. Defaults to
            ``0.8``.

    Returns:
        pd.DataFrame with one row per (Precursor, Replicate) pair containing:

        - ``Library Dot Product``     – original value
        - ``Ratio Dot Product``       – original value
        - ``lib_dot_pass``            – True if ≥ *lib_threshold*
        - ``ratio_dot_pass``          – True if ≥ *ratio_threshold*
        - ``both_pass``               – True if both dot products pass

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


# ---------------------------------------------------------------------------
# Retention-time deviation
# ---------------------------------------------------------------------------


def retention_time_deviation(
    df: pd.DataFrame,
    observed_col: str = "Peptide Retention Time",
    predicted_col: str = "Predicted Retention Time",
) -> pd.DataFrame:
    """
    Calculate the absolute and relative deviation between observed and
    predicted retention times for each precursor–replicate pair.

    Large deviations may indicate co-elution interference, incorrect
    peak-picking, or iRT calibration issues.

    Args:
        df: Skyline report DataFrame.
        observed_col: Column with observed retention time. Defaults to
            ``'Peptide Retention Time'``.
        predicted_col: Column with predicted (iRT-based) retention time.
            Defaults to ``'Predicted Retention Time'``.

    Returns:
        pd.DataFrame with one row per (Precursor, Replicate) pair, sorted by
        ``abs_rt_dev`` descending, with columns:

        - ``Precursor``    – precursor string
        - ``Replicate``    – replicate label
        - ``rt_observed``  – observed retention time (minutes)
        - ``rt_predicted`` – predicted retention time (minutes)
        - ``rt_dev``       – signed deviation (observed − predicted, minutes)
        - ``abs_rt_dev``   – absolute deviation (minutes)

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


# ---------------------------------------------------------------------------
# High-level summary
# ---------------------------------------------------------------------------


def summarize_prm(
    df: pd.DataFrame,
    cv_col: str = "RatioLightToHeavy",
    lib_threshold: float = 0.8,
    ratio_threshold: float = 0.8,
    rt_dev_threshold: float = 2.0,
) -> dict:
    """
    Run all PRM quality checks and return a summary dictionary.

    This is a convenience wrapper around :func:`compute_cv`,
    :func:`flag_missing_values`, :func:`dot_product_summary`, and
    :func:`retention_time_deviation`.

    Args:
        df: Skyline report DataFrame (from ``ImportFile.import_skyline_file``).
        cv_col: Column used for CV calculation. Defaults to
            ``'RatioLightToHeavy'``.
        lib_threshold: Minimum acceptable library dot product. Defaults to
            ``0.8``.
        ratio_threshold: Minimum acceptable ratio dot product. Defaults to
            ``0.8``.
        rt_dev_threshold: Maximum acceptable absolute RT deviation in minutes.
            Defaults to ``2.0``.

    Returns:
        dict with the following keys:

        - ``"cv"``               – DataFrame from :func:`compute_cv`
        - ``"missing"``          – DataFrame from :func:`flag_missing_values`
        - ``"dot_products"``     – DataFrame from :func:`dot_product_summary`
        - ``"rt_deviation"``     – DataFrame from :func:`retention_time_deviation`
        - ``"n_precursors"``     – total unique precursors
        - ``"n_replicates"``     – total unique replicates
        - ``"pct_dot_pass"``     – percentage of rows passing both dot-product thresholds
        - ``"pct_rt_within"``    – percentage of rows within *rt_dev_threshold*
        - ``"median_cv_pct"``    – median %CV across all precursors

    Example::

        from skyline_qc import ImportFile
        from skyline_qc.prm import summarize_prm

        df = ImportFile("results.csv").import_skyline_file()
        report = summarize_prm(df)
        print(report["median_cv_pct"])
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


# ---------------------------------------------------------------------------
# Filtering
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

    # Remove rows where 'Normalized Area' is NaN, matching ratio_picking.ipynb data cleaning
    df = df[df['Normalized Area'].notna()]
    # prepare pivot table 

    index_cols = [col for col in df.columns if col not in ['Isotope Label Type', 'Intensity']]  
    pivot_df = df.pivot_table(
        index=['Replicate', 'Protein Name', 'Peptide'],   
        columns='Isotope Label Type',
        values='Normalized Area',
        aggfunc='first').reset_index()
    return pivot_df 

def summarise_peptide_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise the peptide counts in the DataFrame.
    """
    pept_sum = df.groupby(['Protein Name','Peptide']).agg(
        heavy_count=pd.NamedAgg(column='heavy', aggfunc=lambda x: x.notna().sum()),
        light_count=pd.NamedAgg(column='light', aggfunc=lambda x: x.notna().sum())
    ).reset_index() 

    # Order by sum of heavy and light counts
    pept_sum = pept_sum.sort_values(by=['heavy_count', 'light_count'], ascending=False)
    return pept_sum.reset_index(drop=True) 

def report_peptide_protein_summary(peptide_counts):
    """
    Report summary statistics on peptide and protein detection.

    Args:
        skyline_data: DataFrame containing raw skyline data.
        peptide_counts: DataFrame containing summarised peptide counts, must include 'Peptide', 'heavy_count', and 'light_count'.

    Returns:
        summary_dict: Dictionary with summary statistics.
    """
    # Count the number of unique peptides and proteins in the data
    expected_columns = ['Protein Name', 'Peptide', 'heavy_count', 'light_count']
    if not all(col in peptide_counts.columns for col in expected_columns):
        raise ValueError(f"Expected columns {expected_columns} not found in peptide_counts")

    # For peptide_counts, Protein Name may be missing, so map peptide->protein using skyline_data
    num_unique_peptides = peptide_counts['Peptide'].nunique()
    num_unique_proteins = peptide_counts['Protein Name'].nunique()

    summary_dict = {
        'num_unique_peptides': num_unique_peptides,
        'num_unique_proteins': num_unique_proteins,
    }

    print(f"Number of unique peptides: {num_unique_peptides}")
    print(f"Number of unique proteins: {num_unique_proteins}")

    return summary_dict

def plot_heavy_light_scatter(peptide_counts):
    """
    Scatter plot of heavy vs. light peptide counts for each peptide.

    Args:
        peptide_counts (pd.DataFrame): DataFrame with columns 'Peptide', 'heavy_count', 'light_count'
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(peptide_counts['heavy_count'], peptide_counts['light_count'], alpha=0.6)
    plt.xlabel("Heavy Count")
    plt.ylabel("Light Count")
    plt.title("Scatter plot of Heavy vs. Light Peptide Counts")
    plt.grid(True)
    # Optionally draw y=x reference line
    min_val = min(peptide_counts['heavy_count'].min(), peptide_counts['light_count'].min())
    max_val = max(peptide_counts['heavy_count'].max(), peptide_counts['light_count'].max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=1)
    plt.tight_layout()
    plt.show()

# ---------------------------------------------------------------------------
# Plotting
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

def plot_peptide_counts(df: pd.DataFrame) -> None:
    """
    Plot the peptide counts in the DataFrame.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(df['heavy_count'], bins=20, kde=True)
    plt.title("Distribution of Heavy Peptide Counts")
    plt.xlabel("Heavy Peptide Counts")
    plt.show()