import argparse
import os

from skyline_qc.sdrf import validate_sdrf


def run_validate_sdrf(sdrf_path: str) -> None:
    """
    Call validate_sdrf on the given file and print the result.
    """
    abs_path = os.path.abspath(sdrf_path)
    print(f"Validating SDRF file: {abs_path}")

    try:
        ok, message = validate_sdrf(sdrf_path)
    except FileNotFoundError as e:
        print(f"File error: {e}")
        return
    except Exception as e:
        print(f"Unexpected error while validating: {e}")
        return

    status = "SUCCESS" if ok else "FAIL"
    print(f"\nValidation status: {status}")
    print("Message from validator:")
    print(message)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test the skyline_qc.sdrf.validate_sdrf function."
    )
    parser.add_argument(
        "sdrf_file",
        nargs="?",
        default="sdrf/x.sdrf.tsv",
        help="Path to an SDRF .sdrf.tsv file (default: sdrf/x.sdrf.tsv)",
    )
    args = parser.parse_args()

    run_validate_sdrf(args.sdrf_file)


if __name__ == "__main__":
    main()

