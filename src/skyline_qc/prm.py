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
    """
    Drop duplicate / MultiIndex column names and reset index so patsy sees a flat table.
    """
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
    Report summary statistics on peptide and protein detection. Then export list of 
    peptides. 
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

    peptide_list = list(pd.unique(peptide_counts['Peptide']))

    print(f"Number of unique peptides: {num_unique_peptides}")
    print(f"Number of unique proteins: {num_unique_proteins}")
    print(f"List of selected peptides: {peptide_list}")

    return summary_dict, peptide_list


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


def calculate_inter_plate_cv(peptide_plate_stats):
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



def extract_top_percentile(df, column, percentile=0.1, id_col='Peptide Sequence', source_df=None, source_col=None):
    """
    Extract unique IDs from `id_col` where values in `column` are at or below the given percentile,
    and return both the ID list and filtered DataFrame from `source_df` (if provided).

    Args:
        df (pd.DataFrame): DataFrame containing summary/statistics (e.g. interplate_cv).
        column (str): Name of column to compute percentile threshold over (e.g. 'inter_plate_cv').
        percentile (float): Fraction for percentile threshold (e.g. 0.1 for 10% lowest values).
        id_col (str): Column in `df` whose unique values to extract (peptide identifier).
        source_df (pd.DataFrame, optional): DataFrame to filter based on the returned ID list.
        source_col (str, optional): Column of `source_df` to match IDs (default: id_col).

    Returns:
        tuple: (ID list, filtered DataFrame [if source_df given, else None])
    """
    threshold = df[column].quantile(percentile)
    id_list = df[df[column] <= threshold][id_col].unique()
    if source_df is not None:
        if source_col is None:
            source_col = id_col
        filtered_df = source_df[source_df[source_col].isin(id_list)]
        return id_list, filtered_df.reset_index(drop=True)
    else:
        return id_list, None
        

def plate_peptide_anova(
    selected_norm_peptides,
    log_transform: bool = False,
):
    """
    Fit a two-way ANOVA to assess plate and peptide effects on RatioLightToHeavy, or log-transformed ratio.

    Args:
        selected_norm_peptides (pd.DataFrame): Must contain columns for
            'characteristics[Plate]' or 'Plate',
            'Peptide Sequence' or 'Peptide',
            and 'RatioLightToHeavy'.
        log_transform (bool): If True, analyze log(RatioLightToHeavy) instead of RatioLightToHeavy.

    Returns:
        model: statsmodels OLS result
        anova_res: ANOVA table (pd.DataFrame)
        df: Cleaned DataFrame used in the model (with column 'response')
    """
    # Flexible column handling
    df = selected_norm_peptides.copy()
    if "characteristics[Plate]" in df.columns:
        df = df.rename(columns={"characteristics[Plate]": "Plate"})
    if "Peptide Sequence" in df.columns:
        df = df.rename(columns={"Peptide Sequence": "Peptide"})

    # Only keep required cols
    for c in ("Plate", "Peptide", "RatioLightToHeavy"):
        if c not in df.columns:
            raise KeyError(f"Missing column: {c}")
    df = df[["Plate", "Peptide", "RatioLightToHeavy"]].copy()

    # Drop NaNs & ensure all values are strings for categoricals
    df = df.dropna(subset=["Plate", "Peptide", "RatioLightToHeavy"])
    df["Plate"] = df["Plate"].astype(str)
    df["Peptide"] = df["Peptide"].astype(str)
    df["RatioLightToHeavy"] = pd.to_numeric(df["RatioLightToHeavy"], errors="coerce")

    # Set response
    if log_transform:
        df = df[df["RatioLightToHeavy"] > 0]  # drop nonpositive ratios for log
        df["response"] = np.log(df["RatioLightToHeavy"])
        model_formula = "response ~ C(Plate) + C(Peptide)"
    else:
        df["response"] = df["RatioLightToHeavy"]
        model_formula = "response ~ C(Plate) + C(Peptide)"

    # Fit model
    model = smf.ols(model_formula, data=df).fit()
    anova_res = anova_lm(model, typ=2)
    print(model.summary())
    print(anova_res)

    # Plate p-value printout (robust logic)
    f_col = next((c for c in anova_res.columns if "PR" in c or "p" in c), None)
    plate_row = next((i for i in anova_res.index if "plate" in str(i).lower()), None)
    if f_col and plate_row:
        plate_pval = anova_res.loc[plate_row, f_col]
        print(f"Plate effect p-value: {plate_pval}")
        if pd.isnull(plate_pval):
            print("Warning: Plate effect p-value is NaN.")
        elif plate_pval < 0.05:
            print("Plate effect is significant. There is a significant effect of plate on the response.")
        else:
            print("Plate effect is not significant. There is no significant effect of plate on the response.")
    else:
        print("Plate effect p-value could not be determined.")

    return model, anova_res, df

def fit_plate_logratio_model(selected_norm_peptides):
    """
    Fits a linear model log(RatioLightToHeavy) ~ Plate (Plate as categorical).
    Returns the fitted model and the cleaned DataFrame (with Plate and log_ratio columns).

    Args:
        selected_norm_peptides (pd.DataFrame): DataFrame containing at least 'characteristics[Plate]', 'Replicate', and 'RatioLightToHeavy'
    
    Returns:
        model: statsmodels OLS fitted model
        df: DataFrame with columns ['Plate', 'Replicate', 'RatioLightToHeavy', 'log_ratio', ...]
    """
    # Ensure required columns exist
    required_cols = ['characteristics[Plate]', 'RatioLightToHeavy', 'Replicate']
    for col in required_cols:
        if col not in selected_norm_peptides.columns:
            raise ValueError(f"Missing required column: '{col}'")

    df = selected_norm_peptides.rename(
        columns={'characteristics[Plate]': 'Plate'}
    ).copy()

    # Ensure numeric RatioLightToHeavy
    df['RatioLightToHeavy'] = pd.to_numeric(df['RatioLightToHeavy'], errors='coerce')
    df = df.dropna(subset=['RatioLightToHeavy', 'Plate'])

    # Create log_ratio column
    df['log_ratio'] = np.log(df['RatioLightToHeavy'])

    # Model: log(ratio) ~ Plate (as categorical)
    m = smf.ols('log_ratio ~ C(Plate)', data=df).fit()

    # print the verdict from the model
    print(m.summary())
    print(f"Plate effect p-value: {m.pvalues['C(Plate)']}")
    if m.pvalues['C(Plate)'] < 0.05:
        print("Plate effect is significant. There is a significant effect of plate on the ratio.")
    else:
        print("Plate effect is not significant. There is no significant effect of plate on the ratio.")

    return m, df



def plot_logratio_by_plate_boxplot(df):
    """
    Plots boxplots of log(RatioLightToHeavy) by Replicate, colored by Plate.

    Args:
        df (pd.DataFrame): DataFrame with columns ['Replicate', 'log_ratio', 'Plate']
    """
    # Sort dataframe by Plate for plotting (optional)
    df_sorted = df.sort_values('Plate')
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Replicate', y='log_ratio', data=df_sorted, hue='Plate', dodge=False)
    plt.title('Boxplot of log(RatioLightToHeavy) by Replicate (colored by Plate)')
    plt.xlabel('')
    plt.ylabel('log(RatioLightToHeavy)')
    plt.xticks([], [])  # Remove x tick labels and marks
    plt.legend(title='Plate', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
    return df_sorted

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
