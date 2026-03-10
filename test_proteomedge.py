#!/usr/bin/env python
"""
Test fetch_qreps_table from skyline_qc.proteomedge.

Usage (from repo root):
  python test_proteomedge.py LOT_NUMBER
  python test_proteomedge.py https://proteomedge.com/lotdata/23002/
  python test_proteomedge.py        # prompts for input
"""

import sys

from skyline_qc.proteomedge import fetch_qreps_table


def main() -> None:
    if len(sys.argv) >= 2:
        link_or_lot = sys.argv[1]
    else:
        link_or_lot = input("Enter ProteomEdge lot number or full URL: ").strip()

    if not link_or_lot:
        print("No lot or URL provided. Exiting.")
        sys.exit(1)

    print(f"Fetching qRePS table for: {link_or_lot}")
    try:
        df = fetch_qreps_table(link_or_lot)
    except Exception as exc:  # requests errors, ValueError, etc.
        print(f"Error fetching qRePS table: {exc}")
        sys.exit(1)

    # Save the DataFrame to a local CSV file instead of displaying it or showing HTML
    outfilename = "qreps_ratio_table.csv"
    df.to_csv(outfilename, index=False)
    print(f"OK: loaded {len(df)} rows, {len(df.columns)} columns")
    print("Columns:", list(df.columns))
    print(f"Saved table as {outfilename}")


if __name__ == "__main__":
    main()

