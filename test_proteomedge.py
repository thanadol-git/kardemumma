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


def extract_lot_number(link_or_lot: str) -> str:
    # Match ProteomEdge lot URLs or return the input if it looks like a lot number
    url_pattern = r'/lotdata/(\w+)/'
    match = re.search(url_pattern, link_or_lot)
    if match:
        return match.group(1)
    # If only the lot number is passed, ensure only valid filename chars
    lot = re.sub(r'\W+', '', link_or_lot.strip())
    return lot


def main() -> None:
    # Check if the user provided a lot number or URL
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
    except Exception as exc:
        print(f"Error fetching qRePS table: {exc}")
        sys.exit(1)

    from datetime import datetime

    lot_number = extract_lot_number(link_or_lot)
    if not lot_number:
        print("Could not determine lot number for filename.")
        sys.exit(1)
    today_str = datetime.now().strftime("%Y%m%d")
    out_file = f"{today_str}_{lot_number}_qRePs.csv"

    try:
        df.to_csv(out_file, index=False)
        print(f"OK: loaded {len(df)} rows, {len(df.columns)} columns")
        print(df.head())
        print(f"Table was saved to CSV: {out_file}")
    except Exception as exc:
        print(f"Could not save CSV file: {exc}")
        print("Displaying first rows only:")
        print(df.head())

if __name__ == "__main__":
    main()

