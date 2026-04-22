"""
SDRF (Sample and Data Relationship Format) validation using sdrf-pipelines.
Requires: pip install sdrf-pipelines[ontology]
"""
import os
import subprocess
import sys
from typing import Tuple

import numpy as np
import pandas as pd

def validate_sdrf(sdrf_file: str) -> Tuple[bool, str]:
    """
    Validate an SDRF file using sdrf-pipelines (parse_sdrf validate-sdrf).

    Args:
        sdrf_file: Path to the SDRF file (e.g. .sdrf.tsv).

    Returns:
        Tuple of (success: bool, message: str). success is True if validation passed.

    Raises:
        FileNotFoundError: If sdrf_file does not exist.
        FileNotFoundError: If parse_sdrf CLI is not installed (install with: pip install sdrf-pipelines).
    """
    if not os.path.exists(sdrf_file):
        raise FileNotFoundError(f"SDRF file not found: {sdrf_file}")

    # Prefer the parse_sdrf that belongs to the current Python environment
    scripts_dir = os.path.dirname(sys.executable)
    parse_sdrf_path = os.path.join(scripts_dir, "parse_sdrf")
    cmd = [parse_sdrf_path, "validate-sdrf", "--sdrf_file", os.path.abspath(sdrf_file)]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "Could not find 'parse_sdrf' in the current environment. "
            "Make sure sdrf-pipelines is installed in this environment "
            "and that you are using the correct Python/venv."
        ) from exc

    # Check the subprocess result for SDRF validation
    if result.returncode == 0:
        msg = result.stdout.strip() or "SDRF validation passed."
        return True, msg
    else:
        # Prefer stderr for error messages, fallback to stdout, then make a generic message
        msg = (
            result.stderr.strip()
            or result.stdout.strip()
            or f"Validation failed with exit code {result.returncode}."
        )
        return False, msg

def readout_ms_type(sdrf_file: str) -> list:
    """
    Return a list of unique values from the
    'comment[proteomics data acquisition method]' column in the SDRF file,
    preserving their original row order.

    Args:
        sdrf_file: Path to the SDRF file (e.g. .sdrf.tsv).

    Returns:
        list: Unique MS types, in order of appearance.
    """
    if not os.path.exists(sdrf_file):
        raise FileNotFoundError(f"SDRF file not found: {sdrf_file}")

    df = pd.read_csv(sdrf_file, sep="\t")
    col_name = "comment[proteomics data acquisition method]"
    if col_name not in df.columns:
        raise ValueError(
            f"Column '{col_name}' not found in SDRF file: {sdrf_file}"
        )

    return df[col_name].drop_duplicates().tolist()

def csv_to_tsv(csv_file: str, tsv_file: str) -> None:
    """
    Convert a CSV file to a TSV file.

    Example:
        >>> csv_to_tsv("sample.csv", "sample.tsv")
        # This will read 'sample.csv' (comma-separated) and
        # output it as 'sample.tsv' (tab-separated), preserving columns and data.
    """
    df = pd.read_csv(csv_file, sep=",")
    df.to_csv(tsv_file, sep="\t", index=False)

def remove_whitespace(sdrf_file: str, out_file: str | None = None) -> None:
    """
    Remove leading and trailing whitespace from all string cells in the SDRF file.

    Args:
        sdrf_file: Path to the input SDRF file.
        out_file: Path to write the cleaned file. If None, overwrites the input file.
    """
    df = pd.read_csv(sdrf_file, sep="\t")
    mask = df.apply(
        lambda col: col.map(lambda x: isinstance(x, str) and (x != x.strip()))
    )
    for col in df.columns:
        str_mask = mask[col]
        if str_mask.any():
            df.loc[str_mask, col] = df.loc[str_mask, col].map(
                lambda x: x.strip() if isinstance(x, str) else x
            )
    output_path = out_file if out_file is not None else sdrf_file
    df.to_csv(output_path, sep="\t", index=False)

def detect_trailing_whitespace(sdrf_file: str) -> bool:
    """
    Detect trailing whitespace in the SDRF file. Print and highlight locations and columns in the dataframe if found.

    Returns True if any leading/trailing whitespace is detected, False if the file is clean.
    """
    df = pd.read_csv(sdrf_file, sep="\t")
    has_whitespace = False
    # pandas >= 3.0 removed DataFrame.applymap; build a boolean mask column-wise
    mask = df.apply(
        lambda col: col.map(lambda x: isinstance(x, str) and (x != x.strip()))
    )
    problem_columns = set()

    if mask.values.any():
        has_whitespace = True
        print("Trailing or leading whitespace detected at the following locations:")
        for (row, col), flagged in np.ndenumerate(mask.values):
            if flagged:
                row_label = row + 2  # +2 for header and zero indexing
                col_name = df.columns[col]
                problem_columns.add(col_name)
                original = df.iloc[row, col]
                print(
                    f"  Row {row_label}, Column '{col_name}': ->{original!r}<-"
                )
        if problem_columns:
            print(
                "\nColumns with one or more cells containing leading/trailing whitespace:"
            )
            for col_name in problem_columns:
                print(f"  - {col_name}")
    else:
        print("No leading or trailing whitespace detected in the SDRF file.")
    return has_whitespace