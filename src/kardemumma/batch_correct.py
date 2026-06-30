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
