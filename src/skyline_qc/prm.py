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
    DEPRECATED/MISLEADING: Use get_peptides_below_cv_percentile instead.

    This function is not the correct approach for extracting peptides by percentile,
    and its semantics (percentile as a fraction, not value) may be confusing.

    Please use get_peptides_below_cv_percentile(peptide_means, percentile) for correct peptide selection.

    Args:
        df (pd.DataFrame): DataFrame with peptide statistics.
        column (str): Column to compute the percentile from (e.g. 'inter_plate_cv').
        percentile (float): Percentile (as 0-100); will be interpreted as *fraction* here, which is confusing.
        id_col (str): Column to extract IDs from (e.g. 'Peptide Sequence').
        source_df (pd.DataFrame, optional): DataFrame to filter, matching on provided peptides.
        source_col (str, optional): Column in source_df to match IDs (defaults to id_col).

    Returns:
        tuple: (array of peptide IDs, filtered DataFrame if source_df is given else None)
    """
    # Warn about function misuse
    import warnings
    warnings.warn(
        "extract_top_percentile is deprecated/wrong. "
        "Use get_peptides_below_cv_percentile instead.", 
        DeprecationWarning
    )
    # This is intentionally inconsistent with percentile semantics elsewhere.
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
    return_conversion: bool = False,
):
    """
    Fit a two-way ANOVA on ratio ~ plate + peptide.

    Args:
        selected_norm_peptides: DataFrame with ``RatioLightToHeavy`` and peptide / plate columns.
        plate_col: Source column for plate (e.g. ``'characteristics[Plate]'`` or ``'Plate'``).
            Renamed internally to ``Plate`` for the formula ``C(Plate)``.
        log_transform: If True, model ``log(RatioLightToHeavy)`` (only rows with ratio > 0).

    Returns:
        By default: ``model`` (statsmodels OLS result), ``anova_res`` (Type II ANOVA table).
        If ``return_conversion=True``: also returns
        ``conv_factors_df`` and ``conversion_factors``.

    See :func:`get_plate_conversion_factors` for plate correction factors (e.g. on the
    dataframe returned by :func:`fit_plate_logratio_model`).
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

    # Calculate conversion factors
    conv_factors_df, conversion_factors, _ = get_plate_conversion_factors(df, log_transform=log_transform)

    if return_conversion:
        return model, anova_res, conv_factors_df, conversion_factors
    return model, anova_res


def get_plate_conversion_factors(df, col_plate: str = "Plate", log_transform: bool = False):
    """
    Fit ``log_ratio ~ C(Plate)`` and return multiplicative correction factors per plate using the *global median* as the reference.

    Each plate's conversion factor is: median_of_all / median_of_plate, effectively scaling each plate's values to the median of all measurements.

    Args:
        df: DataFrame with ``Plate`` and ``log_ratio``, or ``Plate`` and
            ``RatioLightToHeavy`` (positive values used; log is taken for the fit).
        log_transform: Reserved for API compatibility; correction is always defined
            in log space (equivalent to multiplicative factors on the raw ratio).

    Returns:
        Tuple: (conv_factors_df, conversion_factors dict, None)
    """
    _ = log_transform
    work = _formula_clean_frame(df.copy())
    if col_plate in work.columns and col_plate != "Plate":
        work = work.rename(columns={col_plate: "Plate"})
    elif col_plate in work.columns and "Plate" in work.columns:
        work = work.drop(columns=[col_plate])
    work = _formula_clean_frame(work)

    if "Plate" not in work.columns:
        raise ValueError(
            f"get_plate_conversion_factors needs 'Plate' (or {col_plate}). "
            f"Have: {list(work.columns)!r}"
        )

    # Plate as string categorical (for consistent grouping)
    plate = _patsy_scalar_categorical(_as_1d_series(work, "Plate")).astype(str)

    if "log_ratio" in work.columns:
        lr = pd.to_numeric(_as_1d_series(work, "log_ratio"), errors="coerce")
        conv_df = pd.DataFrame({"Plate": plate, "log_ratio": lr})
    elif "RatioLightToHeavy" in work.columns:
        ratio = pd.to_numeric(
            _patsy_scalar_categorical(_as_1d_series(work, "RatioLightToHeavy")),
            errors="coerce",
        )
        conv_df = pd.DataFrame({"Plate": plate, "RatioLightToHeavy": ratio})
        conv_df = conv_df.dropna(subset=["Plate", "RatioLightToHeavy"])
        conv_df = conv_df.loc[conv_df["RatioLightToHeavy"] > 0].copy()
        conv_df["log_ratio"] = np.log(conv_df["RatioLightToHeavy"])
        conv_df = conv_df[["Plate", "log_ratio"]]
    else:
        raise ValueError(
            "get_plate_conversion_factors needs 'log_ratio' or 'RatioLightToHeavy'. "
            f"Have: {list(work.columns)!r}"
        )

    conv_df = conv_df.dropna(subset=["log_ratio", "Plate"]).reset_index(drop=True)
    if conv_df.empty:
        raise ValueError("get_plate_conversion_factors: no rows left after cleaning.")

    # Compute global median and per-plate medians (work in *raw* ratio space for factors)
    # Undo log-transform for medians
    conv_df['raw_ratio'] = np.exp(conv_df['log_ratio'])
    global_median = conv_df['raw_ratio'].median()
    plate_medians = conv_df.groupby('Plate')['raw_ratio'].median()

    # Conversion factors scale each plate's median to global median
    conversion_factors: dict[str, float] = {
        plate: (global_median / median_val) if median_val > 0 else 1.0
        for plate, median_val in plate_medians.items()
    }
    conv_factors_df = pd.DataFrame(
        list(conversion_factors.items()),
        columns=["Plate", "conversion_factor"],
    )
    logger.info("Conversion factors per plate:\n%s", conv_factors_df.to_string(index=False))

    # No statsmodels model is used here (None for API compatibility)
    return conv_factors_df, conversion_factors, None


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
    df['RatioLightToHeavy_adj'] = [
        r / get_factor(p) for r, p in zip(df['RatioLightToHeavy'], df['Plate'])
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
