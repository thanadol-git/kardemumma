import argparse
import os
import shutil
import tempfile

import pandas as pd

from skyline_qc.sdrf import remove_whitespace, detect_trailing_whitespace


def run_remove_whitespace(sdrf_path: str) -> None:
    """
    Test the remove_whitespace function on a copy of the given SDRF file,
    and use detect_trailing_whitespace before and after.
    """
    src_abs = os.path.abspath(sdrf_path)
    print(f"Testing remove_whitespace on: {src_abs}")

    if not os.path.exists(src_abs):
        print(f"File error: SDRF file not found: {src_abs}")
        return

    # Work on a temporary copy so the original file is untouched
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_sdrf = os.path.join(tmpdir, os.path.basename(src_abs))
        shutil.copy2(src_abs, tmp_sdrf)

        print(f"  Working on temporary copy: {tmp_sdrf}")

        # Show a small sample *before* cleaning
        df_before = pd.read_csv(tmp_sdrf, sep="\t")
        print("\nBefore cleaning (first 3 rows):")
        print(df_before.head(3))

        print("\nDetecting whitespace BEFORE cleaning:")
        detect_trailing_whitespace(tmp_sdrf)

        # Run the function under test
        remove_whitespace(tmp_sdrf)

        # Show a small sample *after* cleaning
        df_after = pd.read_csv(tmp_sdrf, sep="\t")
        print("\nAfter cleaning (first 3 rows):")
        print(df_after.head(3))

        print("\nDetecting whitespace AFTER cleaning:")
        detect_trailing_whitespace(tmp_sdrf)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test the skyline_qc.sdrf.remove_whitespace function."
    )
    parser.add_argument(
        "sdrf_file",
        nargs="?",
        default="sdrf/20260227_Project ALS _AD_combined.sdrf.tsv",
        help=(
            "Path to an SDRF .sdrf.tsv file to test on "
            "(default: sdrf/20260227_Project ALS _AD_combined.sdrf.tsv)"
        ),
    )
    args = parser.parse_args()

    run_remove_whitespace(args.sdrf_file)


if __name__ == "__main__":
    main()

