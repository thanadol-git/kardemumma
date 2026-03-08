#!/usr/bin/env python
"""
Test skyline_qc with your Skyline export CSV.
Usage:
  python test_import_file.py path/to/your_skyline_export.csv
  python test_import_file.py   # prompts for path
"""
import sys
from skyline_qc import ImportFile, validate_sdrf


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

def test_import_sdrf_file(sdrf_path: str = "sdrf/sdrf_MARTHA_combined.tsv"):
    imp = ImportFile(sdrf_path)
    df = imp.import_sdrf_file()
    print(f"OK: loaded {len(df)} rows, {len(df.columns)} columns")
    return df


def test_validate_sdrf(sdrf_path: str = "sdrf/sdrf_MARTHA_combined.tsv"):
    """Call validate_sdrf() and print result."""
    ok, msg = validate_sdrf(sdrf_path)
    print("Success:", ok)
    print("Message:", msg)
    return ok


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1].lower() == "sdrf":
        path = sys.argv[2] if len(sys.argv) >= 3 else "sdrf/sdrf_MARTHA_combined.tsv"
        try:
            test_import_sdrf_file(path)
        except FileNotFoundError as e:
            print(f"File error: {e}")
            sys.exit(1)
        except ValueError as e:
            print(f"Validation error: {e}")
            sys.exit(1)
    elif len(sys.argv) >= 2 and sys.argv[1].lower() == "validate-sdrf":
        path = sys.argv[2] if len(sys.argv) >= 3 else "sdrf/sdrf_MARTHA_combined.tsv"
        try:
            ok = test_validate_sdrf(path)
            sys.exit(0 if ok else 1)
        except FileNotFoundError as e:
            print(f"File error: {e}")
            sys.exit(1)
    else:
        main()
