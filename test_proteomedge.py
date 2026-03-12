#!/usr/bin/env python
"""
Test fetch_qreps_table from skyline_qc.proteomedge.

Usage (from repo root):
  python test_proteomedge.py LOT_NUMBER
  python test_proteomedge.py https://proteomedge.com/lotdata/23002/
  python test_proteomedge.py        # prompts for input
"""

import sys
import os
import re

from skyline_qc.proteomedge import fetch_qreps_table


from skyline_qc.proteomedge import fetch_qreps_table, extract_lot_number, load_qRePs

def test_load_qRePs():
    # Check if the user provided a lot number or URL
    if len(sys.argv) >= 2:
        link_or_lot = sys.argv[1]
    else:
        link_or_lot = input("Enter ProteomEdge lot number or full URL: ").strip()

    if not link_or_lot:
        print("No lot or URL provided. Exiting.")
        sys.exit(1)

    df, out_file = load_qRePs(link_or_lot)

    try:
        print(f"OK: loaded {len(df)} rows, {len(df.columns)} columns")
        print(df.head())
        print(f"Table was saved to CSV: {out_file}")
    except Exception as exc:
        print(f"Could not display output: {exc}")
        print("Displaying first rows only:")
        print(df.head())


if __name__ == "__main__":
    main()

