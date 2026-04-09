"""

OpenSWATH analysis functions for Skyline exports.

This module provides functions to assess quantification quality from OpenSWATH results.
It expects a DataFrame produced by ``ImportFile.import_skyline_file`` and works with the following key columns:

- ``run_id``                       – unique identifier for each run
- ``filename``                     – original raw data file name
- ``peptide_id``                   – unique peptide identifier
- ``transition_group_id``          – identifier for transition group
- ``decoy``                        – boolean or indicator for decoy entries
- ``RT``                           – observed retention time (minutes)
- ``assay_rt``                     – assay/reference retention time (minutes)
- ``delta_rt``                     – RT difference between observed and assay
- ``assay_RT``                     – alternative reference RT (if present)
- ``delta_RT``                     – alternative RT difference (if present)
- ``feature_id``                   – internal feature identifier
- ``ProteinId``                    – protein identifier
- ``Sequence``                     – stripped peptide sequence
- ``FullPeptideName``              – unmodified or modified peptide sequence
- ``Charge``                       – peptide charge state
- ``mz``                           – mass/charge of precursor
- ``VAR_BSERIES_SCORE``            – b-ion series score
- ``VAR_DOTPROD_SCORE``            – observed dot product score
- ``VAR_INTENSITY_SCORE``          – intensity fit score
- ``VAR_ISOTOPE_CORRELATION_SCORE``– isotopic pattern correlation score
- ``VAR_ISOTOPE_OVERLAP_SCORE``    – isotopic overlap score
- ``VAR_LIBRARY_CORR``             – library correlation score
- ``VAR_LIBRARY_DOTPROD``          – spectral library dot product
- ``VAR_LIBRARY_MANHATTAN``        – Manhattan distance to library
- ``VAR_LIBRARY_RMSD``             – RMSD to library spectrum
- ``VAR_LIBRARY_ROOTMEANSQUARE``   – root-mean-square error to library
- ``VAR_LIBRARY_SANGLE``           – spectral angle to library
- ``VAR_LOG_SN_SCORE``             – log signal-to-noise score
- ``VAR_MANHATTAN_SCORE``          – Manhattan distance score
- ``VAR_MASSDEV_SCORE``            – mass deviation score
- ``VAR_MASSDEV_SCORE_WEIGHTED``   – weighted mass deviation
- ``VAR_MI_SCORE``                 – mutual information score
- ``VAR_MI_WEIGHTED_SCORE``        – weighted mutual information score
- ``VAR_NORM_RT_SCORE``            – normalized RT score
- ``VAR_XCORR_COELUTION``          – cross-correlation (co-elution)
- ``VAR_XCORR_COELUTION_WEIGHTED`` – weighted co-elution score
- ``VAR_XCORR_SHAPE``              – cross-correlation shape score
- ``VAR_XCORR_SHAPE_WEIGHTED``     – weighted shape score
- ``VAR_YSERIES_SCORE``            – y-ion series score
- ``d_score``                      – discriminant score
- ``peak_group_rank``              – rank of the peak group by score
- ``p_value``                      – statistical p-value
- ``q_value``                      – estimated false discovery rate (FDR)
- ``pep``                          – posterior error probability

Plotting functions are in :mod:`skyline_qc.openswath_plots`.
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def import_openswath_file(file_path: str, remove_file_path: bool = False) -> pd.DataFrame:
    """
    Import an OpenSWATH results file.
    Args:
        file_path: Path to the OpenSWATH results file.
    Returns:
        pd.DataFrame: The OpenSWATH results DataFrame.
    """


    # 1. Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    # 2. Check file type
    if not file_path.lower().endswith('.tsv'):
        raise ValueError("Provided file is not a TSV.")

    # 3. Read the file
    df = pd.read_csv(file_path, sep='\t')

    # 4. Check if the file is empty
    if df.empty:
        raise ValueError(f"The file {file_path} is empty.")

    # 5. Check if the file has the expected columns
    
    expected_columns = [
        'run_id',
        'transition_group_id',
        'feature_id',
        'decoy',
        'RT',
        'assay_rt',
        'delta_rt',
        'assay_RT',
        'delta_RT',
        'peak_group_rank',
        'p_value',  
        'VAR_LIBRARY_DOTPROD',
        'VAR_LIBRARY_CORR',
        'VAR_LIBRARY_MANHATTAN',
        'VAR_LIBRARY_RMSD',
        'VAR_LIBRARY_ROOTMEANSQUARE',
        'VAR_LIBRARY_SANGLE',
        'VAR_LIBRARY_DOTPROD',
        'VAR_LIBRARY_CORR',
        'VAR_LIBRARY_MANHATTAN',
        'VAR_LIBRARY_RMSD',
    ]   
    missing_columns = [col for col in expected_columns if col not in df.columns]
    if missing_columns: 
        raise ValueError(f"The file {file_path} is missing expected columns: {missing_columns}. "
                         f"Actual columns: {list(df.columns)}")

    if remove_file_path:
        # In column 'filename', extract string after the last '/', only if it ends with .mzML or .raw
        df['filename'] = df['filename'].str.extract(r'([^/]+\.mzML|[^/]+\.raw)$', expand=False)
 
    # Sort by column 'filename'
    sorting_column = ['filename', 'FullPeptideName', 'feature_id']
    df = df.sort_values(sorting_column)

    return df.reset_index(drop=True)


def filter_best_peak_group(df: pd.DataFrame, threshold: float = 0.95) -> pd.DataFrame:
    """
    1. filter tsv for `decoy = 0`
    2. filter tsv for `VAR_LIBRARY_DOTPROD >= 0.95` (or whatever threshold you use in Skyline)
    3. group by `run_id`, and `transition_group_id` to filter for the best peak_group per precursor, run.

    Args:
        df: The OpenSWATH results DataFrame.
    Returns:
        pd.DataFrame: The filtered OpenSWATH results DataFrame.
    """

    # If FullPeptideName contains `K(UniMod:259)` or `R(UniMod:267)` add new column `Isotope Label Type` with value `heavy` else 'light'
    df['Isotope Label Type'] = df['FullPeptideName'].apply(lambda x: 'heavy' if 'K(UniMod:259)' in x or 'R(UniMod:267)' in x else 'light')


    # Remove decoys
    print(f"** Removing decoys **")
    df = df[df['decoy'].astype(str) == '0']

    # # Remove rows where VAR_LIBRARY_DOTPROD is less than the threshold
    # df = df[df['VAR_LIBRARY_DOTPROD'] >= threshold]

    # Filter Peak group rank = 1
    print(f"** Filtering Peak group rank = 1 **")
    df = df[df['peak_group_rank'] == 1]

    # Filter dotprod > threshold
    print(f"** Filtering dotprod > threshold: {threshold} **")
    df = df[df['VAR_LIBRARY_DOTPROD'] >= threshold]

    # # Group by filename, Sequence, transition_group_id
    # df = df.groupby(['filename', 'Sequence', 'transition_group_id']).agg({
    #     'VAR_LIBRARY_DOTPROD': 'mean',
    #     'Isotope Label Type': 'first'
    # }).reset_index()

    # Arrange by filename, Sequene, transition_group_id
    # df = df.sort_values(['filename', 'Sequence', 'transition_group_id'])
    df = df.sort_values(['filename', 'Sequence'])

    return df.reset_index(drop=True)

def remove_precursor(df: pd.DataFrame, charges = list[int]):
    """
    Remove precursors with charges not in the list
    """
    total_precursors = df.shape[0]
    print(f"** Total precursors: {total_precursors} **")
    print(f"** Removing precursors with charges not in the list: {charges} **")
    
    # Removing
    df = df[df['Charge'].isin(charges)]

    after_total_precursors = df.shape[0]
    print(f"** After removing precursors with charges not in the list: {after_total_precursors} **")
    return df.reset_index(drop=True)    


def plot_dotprod_kde(df: pd.DataFrame) -> None:
    """
    Plot the KDE of VAR_LIBRARY_DOTPROD for each Isotope Label Type
    """
    
    plt.figure(figsize=(10, 6))
    # Use long-form DataFrame with explicit x and hue
    sns.kdeplot(
        data=df,
        x='VAR_LIBRARY_DOTPROD',
        hue='Isotope Label Type',
        fill=True,
        common_norm=False,
        alpha=0.6
    )
    plt.title('KDE of VAR_LIBRARY_DOTPROD by Isotope Label Type')
    plt.xlabel('VAR_LIBRARY_DOTPROD')
    plt.legend()
    plt.show()

# Extract RT (later)


def pair_ions_matching(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate heavy/light intensity ratio per peptide in each file.
    """

    # Some inputs can contain duplicated column names (e.g. two "Isotope Label Type"
    # columns), which makes pandas groupby/pivot fail with "not 1-dimensional".
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()].copy()

    required_columns = ['filename', 'ProteinId', 'Sequence', 'Charge', 'Isotope Label Type', 'Intensity']
    
    # Check for missing columns
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns for ratio calculation: {missing_columns}")

    # Use peptide identity as index and pivot isotope label to columns.
    group_col = ['filename', 'ProteinId', 'Sequence', 'Charge']
    df = (
        df.pivot_table(
            index=group_col,
            columns='Isotope Label Type',
            values='Intensity',
            aggfunc='sum',
        )
        .reset_index()
    )

    # Ensure expected isotope columns exist even if one class is absent.
    if 'heavy' not in df.columns:
        df['heavy'] = np.nan
    if 'light' not in df.columns:
        df['light'] = np.nan

    df['ratio_heavy_light'] = df['heavy'] / df['light']

    # Sort by filename, Sequence
    df = df.sort_values(['filename', 'Sequence'])
    return df

# ================================
# Counting ions channel
# ================================

def _select_ions_channel(df: pd.DataFrame) -> pd.DataFrame:
    """
    Count the number of ions in each channel per peptide in each file.
    Input should be in the same format that get from OpenSWATH analysis. Row filtering is applied in the function.
    """

    required_columns = ['filename', 'ProteinId', 'Sequence', 'Charge', 'Isotope Label Type', 'Intensity']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns for count calculation: {missing_columns}")

    # Summarise count of Isotope Label Type for each peptide in each file and charge
    df = (
        df.groupby(['filename', 'ProteinId', 'Sequence', 'Charge', 'Isotope Label Type'], as_index=False)
        .agg(n_ions=('Intensity', 'sum'))
    )

    # Sort by filename, Sequence, Charge
    df = df.sort_values(['filename', 'Sequence', 'Charge', 'Isotope Label Type'])
    # Check for missing columns
    return df.reset_index(drop=True)

def _summarise_ions_channel(df: pd.DataFrame) -> pd.DataFrame:   
    """
    Count the number of ions in each channel per peptide in each file.
    """

    df = _select_ions_channel(df)

    group_col = ['ProteinId', 'Sequence', 'Charge', 'Isotope Label Type']
    # Count the number of rows (ions) in each channel per peptide in each file
    df_count = (
        df.groupby(group_col)
        .size() # count the number of rows (ions) in each channel per peptide in each file
        .reset_index(name='n_ions')
    )

    return df_count.reset_index(drop=True)

def count_ions_channel(df: pd.DataFrame) -> pd.DataFrame:

    df_count = _summarise_ions_channel(df)
    # Pivot Isotope Label Type to columns and values are n_ions
    df_count = df_count.pivot_table(index=['ProteinId', 'Sequence', 'Charge'], columns='Isotope Label Type', values='n_ions', fill_value=0)
    df_count = df_count.reset_index()

    return df_count

def _summarise_ions_channel_count(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise the number of ions in each channel per peptide in each file.
    """
      # Count the number of ions in each channel per peptide in each file
    df_count = count_ions_channel(df)
    
    group_col = ['heavy', 'light']
    df_count = df_count.groupby(group_col).size().reset_index(name='count')
    # Remove column index
    df_count = df_count.reset_index(drop=True)

    return df_count

def plot_ions_channel(df: pd.DataFrame) -> None:
    """
    Plot the number of ions in each channel per peptide in each file.
    """

    df = _summarise_ions_channel_count(df)

    # Plot dot plot where heavy is on x axis and light is on y axis and color by count
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='heavy', y='light', data=df, hue='count')
    plt.title('Dot plot of heavy vs. light ions')
    plt.xlabel('Heavy')
    plt.ylabel('Light')
    plt.show()
    return df

def filter_ions_channel(df: pd.DataFrame, light_cutoff: int , heavy_cutoff: int ) -> pd.DataFrame:

    df = count_ions_channel(df)

    df = df[df['light'] >= light_cutoff]
    df = df[df['heavy'] >= heavy_cutoff]

    # Select Sequence and Charge
    df = df[['Sequence', 'Charge']]
    df = df.drop_duplicates()

    return df
# ========================================================
# Extracting ratio
# ========================================================
def filter_from_peptides(df: pd.DataFrame, peptides: pd.DataFrame) -> pd.DataFrame: 
    """
    Inner join df with peptides on Sequence and Charge
    """
    df = df.merge(peptides, on=['Sequence', 'Charge'], how='inner')
    return df.reset_index(drop=True)

def get_ratio(df: pd.DataFrame, level: str = 'peptide') -> pd.DataFrame:
    """
    Calculate the ratio of heavy to light ions
    """

    # check if level is valid
    valid_levels = ['peptide', 'precursor']
    if level not in valid_levels:
        raise ValueError(f"Invalid level: {level}. Valid levels are: {valid_levels}")

    # Summarise signal quant by level of grouping
    if level == 'peptide':
        group_col = ['filename', 'ProteinId', 'Sequence', 'Isotope Label Type']
    elif level == 'precursor':
        group_col = ['filename', 'ProteinId', 'Sequence', 'Charge', 'Isotope Label Type']
    else:
        raise ValueError(f"Invalid level: {level}")

    df = df.groupby(group_col).agg({'Intensity': 'sum'}).reset_index()

    # Spread Isotope Label Type to columns
    # Remove Isotope Label Type from group_col
    spread_col = list(group_col)
    spread_col.remove('Isotope Label Type')
    
    df = df.pivot_table(index=spread_col, columns='Isotope Label Type', values='Intensity')
    df = df.reset_index()

    # Remove row with NA in heavy or light
    df = df.dropna(subset=['heavy', 'light'])
    
    # Calculate ratio of heavy to light
    df['ratio_heavy_light'] = df['heavy'] / df['light']
    

    # Sort by filename, Sequence
    df = df.sort_values(['filename', 'Sequence'])
    return df

