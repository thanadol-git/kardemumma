"""
Plotting functions for PRM (Parallel Reaction Monitoring) analysis.

All analysis functions are in :mod:`skyline_qc.prm`.
"""

from __future__ import annotations

from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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


def plot_heavy_light_scatter(peptide_counts):
    """
    Scatter plot of heavy vs. light peptide counts for each peptide,
    colored by the density of peptides at each point. Highest density points are plotted on top.

    Args:
        peptide_counts (pd.DataFrame): DataFrame with columns 'Peptide', 'heavy_count', 'light_count'
    """
    x = peptide_counts['heavy_count'].values
    y = peptide_counts['light_count'].values

    xy = list(zip(x, y))
    counts = Counter(xy)
    point_count = np.array([counts[(hx, ly)] for hx, ly in xy])
    sort_idx = np.argsort(point_count)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]
    point_count_sorted = point_count[sort_idx]

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(
        x_sorted, y_sorted,
        c=point_count_sorted,
        cmap='viridis',
        alpha=0.7,
        s=60,
        edgecolors='k',
        linewidth=0.5
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
    """
    Plot the distribution of heavy peptide counts.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(df['heavy_count'], bins=20, kde=True)
    plt.title("Distribution of Heavy Peptide Counts")
    plt.xlabel("Heavy Peptide Counts")
    plt.show()


def plot_pool_boxplot(pool_df):
    """
    Plot a boxplot of log(RatioLightToHeavy) by Replicate, colored by Plate.

    Args:
        pool_df (pd.DataFrame): DataFrame filtered to include only pool samples.
                               Must have columns 'Replicate', 'RatioLightToHeavy', and 'characteristics[Plate]'.
    """
    sns.boxplot(
        x='Replicate',
        y='RatioLightToHeavy',
        data=pool_df,
        hue='characteristics[Plate]'
    )
    plt.title('Boxplot of log(RatioLightToHeavy) by Replicate (colored by Plate)')
    plt.xlabel('')
    plt.ylabel('log(RatioLightToHeavy)')
    plt.yscale('log')
    plt.xticks([], [])
    plt.legend(title='Plate', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()


def plot_pool_heatmap(pool_data, aggfunc: str = "mean"):
    """
    Plot a heatmap of log(RatioLightToHeavy) for Pool samples.
    Peptides (rows) and Replicates (columns) are both ordered by their mean log-ratio.

    If several rows share the same peptide and replicate (e.g. multiple precursors),
    values are aggregated with *aggfunc* (default ``mean``) so the table is unique
    for ``pivot_table``.

    Args:
        pool_data (pd.DataFrame): DataFrame filtered for Pool samples, must have columns
            ``Peptide Sequence``, ``Replicate``, ``RatioLightToHeavy``.
        aggfunc: Passed to :meth:`pandas.DataFrame.pivot_table` (e.g. ``\"mean\"``, ``\"median\"``, ``\"first\"``).
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


def plot_intra_plate_cv_stats(peptide_plate_stats: pd.DataFrame, col_name: str = 'characteristics[Plate]'):
    """
    Plot a boxplot of intra-plate CV per group/peptide from the given stats DataFrame.

    Args:
        peptide_plate_stats (pd.DataFrame): Output from calculate_intra_plate_cv.
        col_name (str): The group column used in the stats DataFrame.
    """
    display(peptide_plate_stats.head())

    plt.figure(figsize=(10, 6))
    sns.boxplot(x=col_name, y='intra_plate_cv', data=peptide_plate_stats)
    plt.title('Boxplot of Intra Plate CV for Pool')
    plt.xlabel(col_name.replace('characteristics[', '').replace(']', '').capitalize())
    plt.ylabel('Intra Plate CV')
    plt.show()


def plot_inter_plate_cv_kde(peptide_plate_stats):
    """
    Plot a KDE of inter-plate CV for each peptide, with a vertical line at the median.

    Args:
        peptide_plate_stats (pd.DataFrame): Output from calculate_intra_plate_cv.
            Must have columns ['Peptide Sequence', 'mean'] at a minimum.

    Returns:
        matplotlib.figure.Figure: The figure object containing the plot.
    """
    peptide_means = (
        peptide_plate_stats.groupby('Peptide Sequence')['mean']
        .agg(['mean', 'std'])
        .rename(columns={'mean': 'grand_mean', 'std': 'between_plate_sd'})
        .reset_index()
    )
    peptide_means['inter_plate_cv'] = peptide_means['between_plate_sd'] / peptide_means['grand_mean']

    fig = plt.figure(figsize=(12, 6))
    sns.kdeplot(peptide_means['inter_plate_cv'].dropna(), fill=True)
    median_cv = peptide_means['inter_plate_cv'].median()
    plt.axvline(median_cv, color='red', linestyle='--', label=f'Median = {median_cv:.2f}')
    plt.title('KDE Plot of Inter-Plate CV Across Peptides')
    plt.xlabel('Inter-Plate CV')
    plt.ylabel('Density')
    plt.legend()
    plt.show()
    return fig



def plot_cumulative_peptide_count_by_cv(peptide_means):
    """
    Plot cumulative number of peptides as a function of sorted inter-plate CV.

    Args:
        peptide_means (pd.DataFrame): DataFrame with at least ['inter_plate_cv', 'Peptide Sequence'] columns.
    """
    cv_sorted = peptide_means[['inter_plate_cv', 'Peptide Sequence']].sort_values('inter_plate_cv').reset_index(drop=True)
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

