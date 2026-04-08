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
 

    return df


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
    df = df[df['decoy'] == 0]
    df = df[df['VAR_LIBRARY_DOTPROD'] >= 0.95]
    df = df.groupby(['run_id', 'transition_group_id']).apply(lambda x: x.loc[x['peak_group_rank'] == 1]).reset_index(drop=True)   
    return df