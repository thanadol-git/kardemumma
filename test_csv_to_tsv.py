import argparse
import os

import pandas as pd

from skyline_qc.sdrf import csv_to_tsv


def run_csv_to_tsv(csv_path: str, tsv_path: str | None = None) -> None:
    """
    Call csv_to_tsv on the given file and print a short summary.
    """
    csv_abs = os.path.abspath(csv_path)

    if tsv_path is None:
        base, _ = os.path.splitext(csv_path)
        tsv_path = base + ".tsv"
    tsv_abs = os.path.abspath(tsv_path)

    print(f"Converting CSV to TSV:")
    print(f"  input : {csv_abs}")
    print(f"  output: {tsv_abs}")

    try:
        csv_to_tsv(csv_abs, tsv_abs)
    except FileNotFoundError as e:
        print(f"File error: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during conversion: {e}")
        return

    # Quick check: read back the TSV and show shape and first few rows
    try:
        df = pd.read_csv(tsv_abs, sep="\t")
        print(f"\nTSV created successfully. Shape: {df.shape}")
        print("\nFirst 5 rows:")
        print(df.head())
    except Exception as e:
        print(f"\nTSV was written but could not be read back: {e}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test the skyline_qc.sdrf.csv_to_tsv function."
    )
    parser.add_argument(
        "csv_file",
        nargs="?",
        default="ratio/DE17501_ratio.csv",
        help="Path to an input CSV file (default: ratio/DE17501_ratio.csv)",
    )
    parser.add_argument(
        "--out",
        dest="tsv_file",
        default=None,
        help="Optional output TSV path (default: same as CSV but with .tsv extension)",
    )
    args = parser.parse_args()

    run_csv_to_tsv(args.csv_file, args.tsv_file)


if __name__ == "__main__":
    main()

