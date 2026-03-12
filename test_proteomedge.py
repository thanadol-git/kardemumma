#!/usr/bin/env python

import sys
from skyline_qc.proteomedge import load_qRePs

def main() -> None:
    if len(sys.argv) >= 2:
        link_or_lot = sys.argv[1]
    else:
        link_or_lot = input("Enter ProteomEdge lot number or full URL: ").strip()

    if not link_or_lot:
        print("No lot or URL provided. Exiting.")
        sys.exit(1)

    print(f"Fetching and saving qRePS table for: {link_or_lot}")
    try:
        df, out_file = load_qRePs(link_or_lot)
    except Exception as exc:
        print(f"Error loading qRePS table: {exc}")
        sys.exit(1)

    print(f"OK: loaded {len(df)} rows, {len(df.columns)} columns")
    print("Columns:", list(df.columns))
    print(df.head())
    print(f"Saved table as {out_file}")

if __name__ == "__main__":
    main()