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

Plotting functions are in :mod:`kardemumma.prm_plots`.
"""

import logging

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

logger = logging.getLogger(__name__)


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

    logger.info("Number of unique peptides: %d", num_unique_peptides)
    logger.info("Number of unique proteins: %d", num_unique_proteins)
    logger.info("Selected peptides: %s", peptide_list)

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

    if col_name not in pool_data.columns:
        raise KeyError(f"Column '{col_name}' not found in pool_data")

    # Select only the columns we need
    pool_data = pool_data[[col_name, 'Peptide Sequence', 'RatioLightToHeavy']]
    pool_data = pool_data.drop_duplicates()

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


# ---------------------------------------------------------------------------
# Plate Normalization
# ---------------------------------------------------------------------------


def get_lowest_cv_peptides(interplate_cv_df, cv_percentile: float):
    """
    Return a list of peptide sequences with inter-plate CV
    at or below the specified percentile.
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
        percentile: Quantile threshold as a fraction in [0, 1] (e.g. ``0.1`` for the
            lowest 10% of values).
        id_col: Column whose unique values to extract (default ``'Peptide Sequence'``).
        source_df: If given, filter this DataFrame to rows whose ``source_col`` is in
            the returned ID list.
        source_col: Column of ``source_df`` to match against IDs (defaults to ``id_col``).

    Returns:
        tuple: ``(id_array, filtered_df)``; ``filtered_df`` is ``None`` if ``source_df``
        is not provided.
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
    selected_norm_peptides,
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

    # Unify Plate column name
    if plate_col in df.columns and plate_col != "Plate":
        if "Plate" in df.columns:
            df = df.drop(columns=["Plate"])
        df = df.rename(columns={plate_col: "Plate"})
    elif "characteristics[Plate]" in df.columns and "Plate" not in df.columns:
        df = df.rename(columns={"characteristics[Plate]": "Plate"})
    elif "characteristics[Plate]" in df.columns and "Plate" in df.columns:
        df = df.drop(columns=["characteristics[Plate]"])

    # Unify Peptide column name
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

    # Clean and recode columns
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

    if df.empty or df.shape[0] == 0:
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

    model_formula = "response ~ C(Plate) + C(Peptide)"
    model = smf.ols(model_formula, data=df).fit()
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
        if pr_col is None:
            return None
        return anova_res.loc[rows[0], pr_col]

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
    df,
    col_plate: str = "Plate",
    ratio_col: str = "RatioLightToHeavy",
    log_transform: bool = False,
):
    """
    Calculate plate conversion factors for the given DataFrame.

    Returns:
        conv_factors_df (pd.DataFrame): DataFrame of plate median and conversion factor.
        conversion_factors (dict): Mapping plate -> multiplicative correction factor.
        model: Fitted OLS model object from statsmodels.
    """
    # Defensive: required columns
    cols = ['Peptide Sequence', 'Replicate', col_plate, ratio_col]
    missing_cols = [col for col in cols if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing required columns in df: {missing_cols}")

    # Work on a copy, drop duplicate rows
    d = df[cols].drop_duplicates().copy()

    # Compute ratio_fit, handle coercion, log if requested
    d = d.assign(
        ratio_fit=pd.to_numeric(d[ratio_col], errors='coerce')
    )
    if log_transform:
        d['ratio_fit'] = np.log(d['ratio_fit'])

    # Construct formula: always reference everything with Q()
    model_formula = f'ratio_fit ~ C(Q("{col_plate}")) + C(Q("Peptide Sequence"))'
    model = smf.ols(model_formula, data=d).fit()

    # Plate medians: median ratio_fit for each plate
    platemed = (
        d.groupby(col_plate, dropna=True, observed=True)['ratio_fit']
        .median()
        .reset_index()
        .rename(columns={'ratio_fit': 'plate_median'})
    )

    # Global peptide median (used as normalization target)
    global_median = d['ratio_fit'].median()

    # Correction factor: for each plate = global median / plate median
    platemed['correction_factor'] = global_median / platemed['plate_median'] 

    # Return as dict: keys must match how Plate appears in the DataFrame
    conversion_factors = {
        str(row[col_plate]): row['correction_factor']
        for _, row in platemed.iterrows()
    }

    # For reproducibility, also output the DataFrame form
    return platemed, conversion_factors, model

# ---------------------------------------------------------------------------
# Plate Normalization
# ---------------------------------------------------------------------------


def adjust_ratio_by_plate(df, conversion_factors):
    """
    Normalize RatioLightToHeavy by dividing by the plate-specific conversion factor.

    Args:
        df (pd.DataFrame): DataFrame with columns 'RatioLightToHeavy' and 'Plate'.
        conversion_factors (dict): Dict mapping plate (as str or int) to factor.

    Returns:
        pd.DataFrame: Input DataFrame with added 'RatioLightToHeavy_adj' column.
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

    # Factor is median_global/median_plate, so we multiply by the factor to get the adjusted ratio
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

    grouped = df.groupby(group_by)[value_col]
    result = grouped.agg(['mean', 'std', 'count']).reset_index()
    result = result.rename(columns={'count': 'n'})
    # Handle division by zero for mean == 0
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

    # Detected means not NA and not exactly zero (to ignore missing and zero/failed quant.)
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

    req_cols = ["Precursor", "Replicate"]
    missing_cols = [c for c in req_cols if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Required columns for dot_product_summary missing: {missing_cols}")

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

    base_cols = []
    for base in ["Precursor", "Replicate"]:
        if base in df.columns:
            base_cols.append(base)
    base_cols += [observed_col, predicted_col]

    result = df[base_cols].copy()
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
    missing_df = flag_missing_values(df, value_col=cv_col)
    dot_df = dot_product_summary(
        df, lib_threshold=lib_threshold, ratio_threshold=ratio_threshold
    )
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
