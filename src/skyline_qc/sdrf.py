"""
SDRF (Sample and Data Relationship Format) validation using sdrf-pipelines.
Requires: pip install sdrf-pipelines
"""
import os
import subprocess
from typing import Tuple


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

    result = subprocess.run(
        ["parse_sdrf", "validate-sdrf", "--sdrf_file", os.path.abspath(sdrf_file)],
        capture_output=True,
        text=True,
        timeout=120,
    )

    # Check the subprocess result for SDRF validation
    if result.returncode == 0:
        msg = result.stdout.strip() or "SDRF validation passed."
        return True, msg
    else:
        # Prefer stderr for error messages, fallback to stdout, then make a generic message
        msg = result.stderr.strip() or result.stdout.strip() or f"Validation failed with exit code {result.returncode}."
        return False, msg

