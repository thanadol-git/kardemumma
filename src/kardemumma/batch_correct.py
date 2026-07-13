"""
Batch correction visualisation for pool QC samples.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_pool_pca(
    pool_df: pd.DataFrame,
    col_ratio: str = "RatioLightToHeavy",
    col_plate: str = "characteristics[plate]",
) -> None:
    """
    Combined PCA plot: Scores scatter + top peptide loadings for PC1/PC2 as lollipop.

    Replicates are projected onto the first two principal components computed
    from log-transformed peptide ratios. Each point is one replicate.
    Alongside, the *top* peptide loadings for each PC are shown in side-by-side lollipop plots.

    Args:
        pool_df: DataFrame with ``Replicate``, ``Peptide Sequence``, *col_ratio*, and *col_plate* columns.
        col_ratio: Column for the light-to-heavy ratio.
        col_plate: Column for plate labels.
    """
    import matplotlib.gridspec as gridspec

    required = ["Replicate", "Peptide Sequence", col_ratio, col_plate]
    missing = [c for c in required if c not in pool_df.columns]
    if missing:
        raise KeyError(f"pool_df missing columns {missing}. Found: {list(pool_df.columns)}")

    pivot = pool_df.pivot_table(
        index="Replicate", columns="Peptide Sequence", values=col_ratio, aggfunc="mean"
    )
    log_pivot = np.log(pivot.clip(lower=1e-12)).dropna(axis=1).dropna(axis=0)

    if log_pivot.shape[0] < 2:
        raise ValueError("Not enough replicates for PCA after dropping missing values.")

    X = log_pivot.values
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-12)

    # SVD: X = U S Vh
    U, S, Vh = np.linalg.svd(X, full_matrices=False)
    coords = U * S
    var_explained = S ** 2 / np.sum(S ** 2) * 100

    # Principal component loadings (columns: peptides, rows: PC)
    # The loading for peptide j on PC k is Vh[k, j]
    peptides = log_pivot.columns
    pc1_loadings = pd.Series(Vh[0, :], index=peptides, name="PC1_loading")
    pc2_loadings = pd.Series(Vh[1, :], index=peptides, name="PC2_loading")
    loadings_df = pd.DataFrame({"PC1_loading": pc1_loadings, "PC2_loading": pc2_loadings})

    # Map replicate to plate
    rep_plate = (
        pool_df[["Replicate", col_plate]]
        .drop_duplicates()
        .set_index("Replicate")[col_plate]
    )
    plot_df = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1]}, index=log_pivot.index)
    plot_df[col_plate] = rep_plate.reindex(log_pivot.index).values

    # Params for subplot arrangement
    top_n = min(12, len(loadings_df))
    fig = plt.figure(figsize=(15, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[2.3, 1.1, 1.1])

    # ---- PCA SCORE (SCATTER) ----
    ax0 = fig.add_subplot(gs[0])
    for plate, grp in plot_df.groupby(col_plate):
        ax0.scatter(grp["PC1"], grp["PC2"], label=plate, s=60, alpha=0.8)
    ax0.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    ax0.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
    ax0.set_title(f"PCA of pool replicates\n(colored by {col_plate})")
    ax0.legend(title=col_plate, bbox_to_anchor=(1.03, 1), loc="upper left", borderaxespad=0.)

    # ---- PC1 LOADINGS ----
    ax1 = fig.add_subplot(gs[1])
    load = loadings_df["PC1_loading"]
    load_abs = load.abs().sort_values(ascending=False).head(top_n)
    load_top = load.loc[load_abs.index]  # Keep sign
    color = "#377eb8"
    ax1.hlines(
        y=range(top_n), xmin=0, xmax=load_top.values,
        color=color, alpha=0.6, linewidth=2,
    )
    ax1.plot(
        load_top.values, range(top_n),
        "o", color=color, markersize=8,
    )
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels(load_top.index)
    ax1.set_xlabel("PC1 Loading Value")
    ax1.set_title(f"Top {top_n} peptides\nPC1")
    ax1.axvline(0, color="grey", linestyle="--", lw=1)

    # ---- PC2 LOADINGS ----
    ax2 = fig.add_subplot(gs[2])
    load = loadings_df["PC2_loading"]
    load_abs = load.abs().sort_values(ascending=False).head(top_n)
    load_top = load.loc[load_abs.index]
    color = "#e41a1c"
    ax2.hlines(
        y=range(top_n), xmin=0, xmax=load_top.values,
        color=color, alpha=0.6, linewidth=2,
    )
    ax2.plot(
        load_top.values, range(top_n),
        "o", color=color, markersize=8,
    )
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels(load_top.index)
    ax2.set_xlabel("PC2 Loading Value")
    ax2.set_title(f"Top {top_n} peptides\nPC2")
    ax2.axvline(0, color="grey", linestyle="--", lw=1)

    # Adjust layout
    plt.suptitle("PCA of pool replicates and top peptide loadings", fontsize=15, y=1.03)
    plt.tight_layout(rect=[0, 0, 1, 0.96], w_pad=2)
    plt.show()

    # If users want loadings, return as well (comment/uncomment as needed)
    # return loadings_df


def _permanova_f(dist_sq: np.ndarray, labels: np.ndarray) -> float:
    """Compute the PERMANOVA pseudo-F statistic from a squared distance matrix."""
    n = len(labels)
    groups = np.unique(labels)
    k = len(groups)

    # Total sum of squares: SS_T = sum of all squared distances / n
    ss_total = np.sum(dist_sq) / n

    # Within-group sum of squares
    ss_within = 0.0
    for g in groups:
        idx = np.where(labels == g)[0]
        ng = len(idx)
        if ng > 1:
            ss_within += np.sum(dist_sq[np.ix_(idx, idx)]) / ng

    ss_between = ss_total - ss_within
    if ss_within == 0:
        return np.inf
    return (ss_between / (k - 1)) / (ss_within / (n - k))


def permanova_batch_effects(
    df: pd.DataFrame,
    col_ratio: str = "RatioLightToHeavy",
    n_permutations: int = 999,
    log_transform: bool = True,
) -> pd.DataFrame:
    """
    Test batch effects for all metadata columns starting with 'characteristics'
    using PERMANOVA on log-transformed peptide ratios.

    For each characteristics column, the function builds a replicate × peptide
    matrix, computes a Euclidean distance matrix, and runs PERMANOVA to test
    whether replicates cluster by that metadata variable.

    Args:
        df: DataFrame with ``Replicate``, ``Peptide Sequence``, *col_ratio*, and
            one or more ``characteristics[*]`` columns.
        col_ratio: Column containing the light-to-heavy ratio.
        n_permutations: Number of label permutations for the p-value (default 999).
        log_transform: Log-transform ratios before computing distances (default True).

    Returns:
        DataFrame with one row per characteristics column and columns:
        ``variable``, ``F_statistic``, ``p_value``, ``R2``, ``n_groups``,
        ``n_samples``, ``n_permutations``.
    """
    char_cols = [c for c in df.columns if c.startswith("characteristics")]
    if not char_cols:
        raise ValueError("No columns starting with 'characteristics' found in df.")

    required = ["Replicate", "Peptide Sequence", col_ratio]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # Build replicate × peptide pivot (one value per cell)
    pivot = df.pivot_table(
        index="Replicate", columns="Peptide Sequence", values=col_ratio, aggfunc="mean"
    ).dropna(axis=1).dropna(axis=0)

    if log_transform:
        pivot = np.log(pivot.clip(lower=1e-12))

    X = pivot.values
    replicates = pivot.index.tolist()

    # Euclidean distance matrix → squared distances
    from sklearn.metrics import pairwise_distances
    dist = pairwise_distances(X, metric="euclidean")
    dist_sq = dist ** 2

    rng = np.random.default_rng(42)
    results = []

    for col in char_cols:
        # Map each replicate to its group label (take first occurrence)
        rep_meta = (
            df[["Replicate", col]]
            .drop_duplicates(subset="Replicate")
            .set_index("Replicate")[col]
        )
        labels = np.array([str(rep_meta.get(r, np.nan)) for r in replicates])

        # Drop replicates with missing metadata
        valid = labels != "nan"
        if valid.sum() < 3 or len(np.unique(labels[valid])) < 2:
            results.append({
                "variable": col,
                "F_statistic": np.nan,
                "p_value": np.nan,
                "R2": np.nan,
                "n_groups": int(np.unique(labels[valid]).size),
                "n_samples": int(valid.sum()),
                "n_permutations": n_permutations,
            })
            continue

        labels_v = labels[valid]
        dist_sq_v = dist_sq[np.ix_(valid, valid)]
        n = len(labels_v)
        k = len(np.unique(labels_v))

        f_obs = _permanova_f(dist_sq_v, labels_v)
        ss_total = np.sum(dist_sq_v) / n
        ss_within = ss_total - (f_obs * ss_total * (k - 1)) / (f_obs * (k - 1) + (n - k))
        r2 = 1.0 - ss_within / ss_total if ss_total > 0 else np.nan

        # Permutation test
        count_ge = 0
        for _ in range(n_permutations):
            perm_labels = rng.permutation(labels_v)
            if _permanova_f(dist_sq_v, perm_labels) >= f_obs:
                count_ge += 1
        p_value = (count_ge + 1) / (n_permutations + 1)

        results.append({
            "variable": col,
            "F_statistic": round(f_obs, 4),
            "p_value": round(p_value, 4),
            "R2": round(r2, 4),
            "n_groups": k,
            "n_samples": n,
            "n_permutations": n_permutations,
        })

    result_df = pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)

    _print_permanova_summary(result_df, n_permutations)

    return result_df


def _r2_label(r2: float) -> str:
    if r2 >= 0.25:
        return "large"
    if r2 >= 0.06:
        return "medium"
    return "small"


def _print_permanova_summary(df: pd.DataFrame, n_permutations: int) -> None:
    sig = df[df["p_value"] < 0.05].copy()
    ns  = df[df["p_value"] >= 0.05].copy()
    skipped = df[df["p_value"].isna()].copy()

    print(f"\nPERMANOVA batch effect summary  ({n_permutations} permutations)")
    print("=" * 60)

    if sig.empty:
        print("No significant batch effects detected (all p >= 0.05).")
    else:
        print(f"Significant factors (p < 0.05):  {len(sig)}")
        for _, row in sig.iterrows():
            effect = _r2_label(row["R2"])
            stars = "***" if row["p_value"] < 0.001 else ("**" if row["p_value"] < 0.01 else "*")
            print(
                f"  {stars}  {row['variable']}\n"
                f"       p = {row['p_value']:.4f}  |  F = {row['F_statistic']:.3f}"
                f"  |  R² = {row['R2']:.3f} ({effect} effect)"
                f"  |  {row['n_groups']} groups, {row['n_samples']} samples"
            )

    if not ns.empty:
        print(f"\nNon-significant factors (p >= 0.05):  {len(ns)}")
        for _, row in ns.iterrows():
            print(
                f"       {row['variable']}\n"
                f"       p = {row['p_value']:.4f}  |  F = {row['F_statistic']:.3f}"
                f"  |  R² = {row['R2']:.3f}"
                f"  |  {row['n_groups']} groups, {row['n_samples']} samples"
            )

    if not skipped.empty:
        print(f"\nSkipped (insufficient data):")
        for _, row in skipped.iterrows():
            print(f"       {row['variable']}  ({row['n_groups']} groups, {row['n_samples']} samples)")

    print("=" * 60)


def correct_ratio_by_factors(
    df: pd.DataFrame,
    factors: list[str],
    col_ratio: str = "RatioLightToHeavy",
) -> pd.DataFrame:
    """
    Correct peptide ratios for batch factors using per-peptide median centering.

    Factors are applied sequentially in the order given (most significant first
    is recommended). For each factor, the correction divides each ratio by its
    peptide-level group median relative to the peptide-level overall median,
    effectively removing the group-level shift while preserving within-group
    biological variation.

    Args:
        df: DataFrame with ``Replicate``, ``Peptide Sequence``, *col_ratio*,
            and all columns listed in *factors*.
        factors: Ordered list of characteristics columns to correct for
            (e.g. from ``permanova_batch_effects``, most significant first).
        col_ratio: Column containing the light-to-heavy ratio.

    Returns:
        Copy of *df* with an additional column ``{col_ratio}_corrected``.
    """
    missing = [c for c in ["Peptide Sequence", col_ratio] + factors if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    out = df.copy()
    corrected_col = f"{col_ratio}_corrected"
    out[corrected_col] = pd.to_numeric(out[col_ratio], errors="coerce")

    for factor in factors:
        # Per-peptide overall median
        peptide_overall = (
            out.groupby("Peptide Sequence")[corrected_col]
            .median()
            .rename("_overall_median")
        )
        # Per-peptide per-group median
        peptide_group = (
            out.groupby(["Peptide Sequence", factor])[corrected_col]
            .median()
            .rename("_group_median")
            .reset_index()
        )
        peptide_group = peptide_group.merge(
            peptide_overall.reset_index(), on="Peptide Sequence"
        )
        # Correction factor: group_median / overall_median
        peptide_group["_factor"] = peptide_group["_group_median"] / peptide_group["_overall_median"]

        out = out.merge(
            peptide_group[["Peptide Sequence", factor, "_factor"]],
            on=["Peptide Sequence", factor],
            how="left",
        )
        valid = out["_factor"] > 0
        out.loc[valid, corrected_col] = out.loc[valid, corrected_col] / out.loc[valid, "_factor"]
        out.drop(columns=["_factor"], inplace=True)

        print(f"  Corrected for: {factor}")

    print(f"\nCorrected ratios written to '{corrected_col}'.")
    return out


def correct_ratio_by_irt(
    df: pd.DataFrame,
    irt_peptides: list[str],
    col_ratio: str = "RatioLightToHeavy",
) -> pd.DataFrame:
    """
    Correct peptide ratios using iRT peptides as negative controls (RUV-style).

    iRT peptides are synthetic standards unaffected by biology. Their per-replicate
    median log-ratio deviation from the overall iRT median is used as an estimate
    of technical noise, which is then subtracted from all peptides in that replicate.

    Args:
        df: DataFrame with ``Replicate``, ``Peptide Sequence``, and *col_ratio*.
        irt_peptides: List of iRT peptide sequences from ``get_irt_peptides()``.
        col_ratio: Column containing the light-to-heavy ratio.

    Returns:
        Copy of *df* with an additional column ``{col_ratio}_corrected``.
    """
    missing = [c for c in ["Replicate", "Peptide Sequence", col_ratio] if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    irt_df = df[df["Peptide Sequence"].isin(irt_peptides)].copy()
    if irt_df.empty:
        raise ValueError("None of the provided iRT peptides were found in df.")

    n_found = irt_df["Peptide Sequence"].nunique()
    print(f"Using {n_found} iRT peptides as negative controls.")

    irt_df["_log_ratio"] = np.log(pd.to_numeric(irt_df[col_ratio], errors="coerce").clip(lower=1e-12))

    # Per-replicate median log-ratio of iRT peptides
    rep_irt_median = irt_df.groupby("Replicate")["_log_ratio"].median().rename("_rep_offset")

    # Overall median log-ratio of iRT peptides across all replicates
    global_irt_median = irt_df["_log_ratio"].median()

    # Offset to subtract per replicate: how much each replicate deviates from global
    rep_offset = (rep_irt_median - global_irt_median).rename("_offset")

    out = df.copy()
    corrected_col = f"{col_ratio}_corrected"
    out["_log_ratio"] = np.log(pd.to_numeric(out[col_ratio], errors="coerce").clip(lower=1e-12))
    out = out.join(rep_offset, on="Replicate")
    out[corrected_col] = np.exp(out["_log_ratio"] - out["_offset"])
    out.drop(columns=["_log_ratio", "_offset"], inplace=True)

    print(f"Corrected ratios written to '{corrected_col}'.")
    return out
