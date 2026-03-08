#!/usr/bin/env python
"""
Test skyline_qc with your Skyline export CSV.
Usage:
  python test_import_file.py path/to/your_skyline_export.csv
  python test_import_file.py   # prompts for path
"""
import sys
from skyline_qc import ImportFile


def main():
    if len(sys.argv) >= 2:
        file_path = sys.argv[1]
    else:
        file_path = input("Path to Skyline CSV: ").strip()

    if not file_path:
        print("No file path given. Exiting.")
        sys.exit(1)

    imp = ImportFile(file_path)
    try:
        df = imp.import_skyline_file()
        print(f"OK: loaded {len(df)} rows, {len(df.columns)} columns")
        print("Columns:", list(df.columns))
        print(df.head())
    except FileNotFoundError as e:
        print(f"File error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Validation error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
