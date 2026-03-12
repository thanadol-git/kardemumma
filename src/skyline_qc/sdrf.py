"""
SDRF (Sample and Data Relationship Format) validation using sdrf-pipelines.
Requires: pip install sdrf-pipelines[ontology]
"""
import os
import subprocess
import sys
from typing import Tuple

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
        raise ValueError(f"Column '{col_name}' not found in SDRF file: {sdrf_file}")

    return df[col_name].drop_duplicates().tolist()