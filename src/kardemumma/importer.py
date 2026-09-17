import os
import re
from typing import Iterable, List, Optional

import pandas as pd

__all__ = [
    "ImportSkylineFile",
    "ImportSDRFFile",
    "MergeFiles",
    "get_irt_peptides",
    "normalize_data_filename",
    "remove_qc_samples",
]

_ACQSUFFIX_RE = re.compile(r"_\d{8,14}$")
_HEAVY_RE = re.compile(r"heavy", re.IGNORECASE)

_SKYLINE_EXPECTED_COLS = [
    "Precursor",
    "Replicate",
    "File Name",
    "Peptide Retention Time",
    "Predicted Retention Time",
    "Precursor Charge",
    "Peptide Sequence",
    "Peptide",
    "Normalized Area",
    "RatioLightToHeavy",
    "Ratio Dot Product",
    "Library Dot Product",
    "Protein Name",
]

_BIOGNOSYS_IRT_PEPTIDES = [
    "LGGNEQVTR",
    "GAGSSEPVTGLDAK",
    "VEATFGVDESNAK",
    "YILAGVENSK",
    "TPVISGGPYEYR",
    "TPVITGAPYEYR",
    "DGLDAASYYAPVR",
    "ADVTPADFSEWSK",
    "GTFIIDPGGVIR",
    "GTFIIDPAAVIR",
    "LFLQFGAQGSPFLK",
]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_path(path: str, ext: str) -> None:
    """Raise if *path* does not exist or does not end with *ext*."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    if not path.lower().endswith(ext):
        raise ValueError(f"Expected a {ext.upper()} file, got: {path}")


def _suggest_qc_replicate_names(names: Iterable[object]) -> List[str]:
    """Return sorted list of names that look like QC samples (contain 'qc' or 'quality')."""
    if names is None:
        return []
    unique_names = {str(x) for x in names if pd.notna(x)}
    targets = ("qc", "quality control", "quality")
    qc_samples = [n for n in unique_names if any(t in n.lower() for t in targets)]
    print(f"Number of possible QC samples: {len(qc_samples)}")
    print(f"Possible QC samples: {qc_samples}")
    print("Consider removing these QC samples before further analysis.")
    return sorted(qc_samples)


def _match_biognosys(peptide_list: List[str]) -> None:
    """Print which peptides in *peptide_list* match the Biognosys iRT kit."""
    unique = list(set(peptide_list))
    matched = [p for p in unique if p in _BIOGNOSYS_IRT_PEPTIDES]
    not_matched = [p for p in unique if p not in _BIOGNOSYS_IRT_PEPTIDES]
    print("Checking if the peptides are from Biognosys...")
    if matched:
        print(f"Matched Biognosys iRT peptides: {', '.join(matched)}")
    if not_matched:
        print(f"Not matched: {', '.join(not_matched)}")
    print(f"Number of peptides not matching Biognosys: {len(not_matched)}")


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def get_irt_peptides(
    df: pd.DataFrame,
    *,
    protein_col: str = "Protein Name",
    peptide_col: str = "Peptide Sequence",
    tag_substring: str = "iRT_Tag",
    match_biognosys: bool = True,
) -> List[str]:
    """
    Return unique peptide sequences for rows whose protein name contains *tag_substring*.

    Args:
        df: Skyline report table.
        protein_col: Column with protein / group name.
        peptide_col: Column with peptide sequence.
        tag_substring: Literal substring to match in ``protein_col`` (not a regex).
        match_biognosys: If True, check matched peptides against the Biognosys iRT kit.

    Returns:
        List of distinct peptide sequences, in first-seen order.
    """
    if protein_col not in df.columns:
        raise ValueError(f"Protein column '{protein_col}' not found in dataframe.")
    if peptide_col not in df.columns:
        raise ValueError(f"Peptide column '{peptide_col}' not found in dataframe.")
    if tag_substring not in df[protein_col].astype(str).unique():
        raise ValueError(f"Tag substring '{tag_substring}' not found in '{protein_col}'.")

    peptide_list = (
        df.loc[
            df[protein_col].astype(str).str.contains(tag_substring, na=False, regex=False),
            peptide_col,
        ]
        .drop_duplicates()
        .tolist()
    )

    if match_biognosys:
        _match_biognosys(peptide_list)

    print(f"\nNumber of iRT peptides: {len(peptide_list)}")
    print(f"iRT peptides: {peptide_list}\n")
    return peptide_list


def normalize_data_filename(
    name: str,
    *,
    strip_acquisition_suffix: bool = True,
) -> str:
    """
    Return the canonical basename for comparing raw filenames between Skyline and SDRF.

    For ``.raw`` files, optionally removes trailing ``_`` + 8–14 digit acquisition
    suffixes, e.g. ``..._Plate_5_C6_20251211114043.raw`` → ``..._Plate_5_C6.raw``.
    """
    s = os.path.basename(str(name).strip().replace("\\", "/"))
    if not strip_acquisition_suffix:
        return s
    base, ext = os.path.splitext(s)
    if ext.lower() != ".raw":
        return s
    while _ACQSUFFIX_RE.search(base):
        base = _ACQSUFFIX_RE.sub("", base, count=1)
    return f"{base}{ext}"


def remove_qc_samples(
    df: pd.DataFrame,
    qc_samples: List[str],
    file_col: str = "File Name",
    ignore_filetype: bool = False,
    print_summary: bool = False,
) -> pd.DataFrame:
    """
    Remove QC samples from the Skyline DataFrame.

    Args:
        df: Input DataFrame containing Skyline data.
        qc_samples: List of sample names to remove (matched against *file_col*).
        file_col: Column containing sample names or filenames. Default: ``"File Name"``.
        ignore_filetype: Strip file extensions before matching. Default: ``False``.
        print_summary: Print a removal summary. Default: ``False``.

    Returns:
        DataFrame with QC samples removed.
    """
    if file_col not in df.columns:
        raise ValueError(f"Column '{file_col}' not found. Columns: {list(df.columns)}")

    working_df = df.copy()
    if ignore_filetype:
        working_df[file_col] = working_df[file_col].astype(str).str.split(".").str[0]
        qc_samples_stripped = [str(n).split(".")[0] for n in qc_samples]
    else:
        qc_samples_stripped = qc_samples

    filtered_df = working_df[~working_df[file_col].isin(qc_samples_stripped)]

    if print_summary:
        n_removed = len(working_df) - len(filtered_df)
        print(f"Rows removed: {n_removed}")
        print(f"QC samples excluded: {sorted(set(qc_samples_stripped))}")

    return filtered_df


# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------


class ImportSkylineFile:
    def __init__(self, file_path: str):
        """
        Args:
            file_path: Path to the Skyline CSV export file.
        """
        self.file_path = file_path
        print(f"Importing Skyline file: {file_path}")

    def import_skyline_file(self) -> pd.DataFrame:
        """Read and validate the Skyline CSV. Adds an ``Isotope Label Type`` column."""
        _validate_path(self.file_path, ".csv")
        df = pd.read_csv(self.file_path)
        if df.empty:
            raise ValueError(f"File is empty: {self.file_path}")
        missing = [c for c in _SKYLINE_EXPECTED_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"Missing expected columns: {missing}. Actual: {list(df.columns)}")
        df["Isotope Label Type"] = (
            df["Precursor"].astype(str)
            .str.contains(_HEAVY_RE.pattern, case=False, regex=True)
            .map({True: "heavy", False: "light"})
        )
        print(f"File {self.file_path} is valid.")
        return df
    
    def irt_signals(self) -> pd.DataFrame:
        """
        Return a DataFrame with the iRT signals from function get_irt_peptides().
        """
        main_df = self.import_skyline_file()
        irt_peptides = get_irt_peptides(main_df)
        
        df_irt = main_df[main_df["Peptide Sequence"].isin(irt_peptides)]
        
        return df_irt

    def suggest_qc_samples(self, df: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Suggest QC replicate names (those containing ``qc`` or ``quality``).

        Args:
            df: Optional already-loaded Skyline table. If omitted, reads from ``self.file_path``.

        Returns:
            Sorted list of suspected QC replicate names.
        """
        if df is None:
            _validate_path(self.file_path, ".csv")
            df = pd.read_csv(self.file_path)
        if "Replicate" not in df.columns:
            raise ValueError(f"Column 'Replicate' not found. Columns: {list(df.columns)}")
        return _suggest_qc_replicate_names(df["Replicate"].unique())


class ImportSDRFFile:
    def __init__(self, file_path: str):
        """
        Args:
            file_path: Path to the SDRF ``.tsv`` file.
        """
        self.file_path = file_path
        _validate_path(file_path, ".tsv")
        print(f"Importing SDRF file: {file_path}")

    def import_sdrf_file(self) -> pd.DataFrame:
        """Read and return the SDRF file as a DataFrame."""
        _validate_path(self.file_path, ".tsv")
        
        
        return pd.read_csv(self.file_path, sep="\t")

    def info(self) -> None:
        """Print a brief summary of the SDRF file (shape and column types)."""
        df = self.import_sdrf_file()
        print(f"Summary for: {self.file_path}")
        print(f"Rows: {df.shape[0]}  Columns: {df.shape[1]}")
        print(df.dtypes.to_string())

    def suggest_qc_samples(self) -> List[str]:
        """Return suspected QC sample names from ``comment[data file]``."""
        df = self.import_sdrf_file()
        return _suggest_qc_replicate_names(df["comment[data file]"].unique())

    def extract_qreps_lot_number(self) -> str:
        """
        Extract the qREPs lot number from ``comment[ProteomEdge]``.

        Returns:
            The unique lot number string.

        Raises:
            ValueError: If the column is missing or contains more than one distinct value.
        """
        df = self.import_sdrf_file()
        col = "comment[ProteomEdge]"
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in SDRF file.")
        lot_numbers = df[col].dropna().unique()
        if len(lot_numbers) != 1:
            raise ValueError(
                f"Expected exactly one lot number in '{col}', "
                f"found {len(lot_numbers)}: {lot_numbers}"
            )
        lot_number = lot_numbers[0]
        if hasattr(lot_number, "item"):
            lot_number = lot_number.item()
        if isinstance(lot_number, float) and lot_number.is_integer():
            lot_number = int(lot_number)
        lot_number = str(lot_number).strip()
        print(f"qREPs lot number: {lot_number}")
        return lot_number

class MergeFiles:
    @staticmethod
    def compare_file_names(
        skyline_df: pd.DataFrame,
        sdrf_df: pd.DataFrame
    ) -> None:
        """
        Print overlaps and differences of sample file names between Skyline and SDRF DataFrames.
        """
        skyline_files = set(skyline_df["File Name"])
        sdrf_files = set(sdrf_df["comment[data file]"])
        not_in_sdrf = skyline_files - sdrf_files
        not_in_skyline = sdrf_files - skyline_files

        print(f"Sample names in Skyline but not in SDRF ({len(not_in_sdrf)}).")
        print("-" * 40)
        print(f"Sample names in SDRF but not in Skyline ({len(not_in_skyline)}).")
   

    def __init__(
        self,
        skyline_df: pd.DataFrame,
        sdrf_df: pd.DataFrame,
        selected_peptides: Optional[List[str]] = None,
    ):
        """
        Args:
            skyline_df: Skyline long-format DataFrame (must include ``Replicate``, ``Peptide``).
            sdrf_df: SDRF DataFrame (must include ``source name``, ``characteristics[Sample]``, characteristics[plate]).
            selected_peptides: Optional; restrict ``skyline_df`` to these ``Peptide`` values.
        """
        self.compare_file_names(skyline_df, sdrf_df)
        self.skyline_df = skyline_df
        self.sdrf_df = sdrf_df
        self.selected_peptides = selected_peptides

    def merge_files(self) -> pd.DataFrame:
        """
        Merge Skyline with SDRF on Replicate and source name (left join). Raise if columns are missing.
        Also merges any SDRF columns that start with 'factor value'.
        """
        # Check if selected_peptides is provided and subset if so
        if self.selected_peptides is not None:
            skyline_df = self.skyline_df[self.skyline_df["Peptide"].isin(self.selected_peptides)].copy()
        else:
            skyline_df = self.skyline_df

        # Columns required for merging
        required_skyline_cols = ["File Name"]
        required_sdrf_cols = ["comment[data file]"] + [
            col for col in self.sdrf_df.columns if col.startswith("characteristics[")
        ]
   

        # Add all columns in sdrf_df that start with 'factor value'
        factor_value_cols = [col for col in self.sdrf_df.columns if col.startswith("factor value")]
        all_sdrf_cols = required_sdrf_cols + factor_value_cols

        # Check for required columns in Skyline dataframe
        missing_skyline = [c for c in required_skyline_cols if c not in skyline_df.columns]
        if missing_skyline:
            raise ValueError(f"Missing column(s) {missing_skyline} in skyline_df")

        # Check for required columns in SDRF dataframe
        missing_sdrf = [c for c in required_sdrf_cols if c not in self.sdrf_df.columns]
        if missing_sdrf:
            raise ValueError(f"Missing column(s) {missing_sdrf} in sdrf_df")

        # Perform left merge of Skyline with selected columns from SDRF
        merged_df = pd.merge(
            skyline_df,
            self.sdrf_df[all_sdrf_cols],
            left_on="File Name",
            right_on="comment[data file]",
            how="left"
        )
        return merged_df

    def select_pool_data(
        self,
        label_col: str = "factor value[Sample]",
        pool_label: str = "Pool",
    ) -> pd.DataFrame:
        """
        Filter merged DataFrame to rows where *label_col* equals *pool_label*.

        Prints a summary of unique sample counts in each plate for the filtered data.
        """
        def _summarize_samples_per_plate(
            input_df: pd.DataFrame,
            plate_col: str = "characteristics[plate]",
            source_col: str = "comment[data file]",
        ) -> None:
            if plate_col not in input_df.columns:
                print(f"  Column '{plate_col}' not in DataFrame; cannot summarize by plate.")
                return
            counts = input_df.groupby(plate_col)[source_col].nunique()
            if counts.empty:
                print("  No plates found.")
            else:
                for plate, count in counts.items():
                    print(f"Plate {plate}: {count} unique comment[data file] names")

        merged_df = self.merge_files()
        if label_col not in merged_df.columns:
            raise ValueError(f"Column '{label_col}' not found in merged DataFrame.")
        filtered = (
            merged_df[merged_df[label_col] == pool_label]
            .sort_values(["Replicate", "Peptide", "Isotope Label Type"])
            .reset_index(drop=True)
        )
        print("Summary of samples per plate:")
        _summarize_samples_per_plate(filtered)
        return filtered

    def mistmatched_samples(self) -> pd.DataFrame:
        """
        Report a DataFrame with file names that are present in only one of:
        Skyline (File Name) or SDRF (comment[data file]). Excludes files that
        are present in both.

        Columns: Files, Skyline, SDRF, suggested date (parsed as 8 consecutive digits from Files).
        """
        import re

        skyline_names = set(self.skyline_df["File Name"].unique())
        sdrf_names = set(self.sdrf_df["comment[data file]"].unique())
        # Only keep files present in exactly one set
        mismatched_files = sorted(skyline_names ^ sdrf_names)
        # Extract the first occurrence of 8 consecutive digits from the filename, if present
        suggested_dates = [
            (re.search(r"\d{8}", f).group(0) if re.search(r"\d{8}", f) else "")
            for f in mismatched_files
        ]
        result_df = pd.DataFrame({
            "Files": mismatched_files,
            "Skyline": [1 if f in skyline_names else 0 for f in mismatched_files],
            "SDRF": [1 if f in sdrf_names else 0 for f in mismatched_files],
            "suggested date": suggested_dates,
        })
        # Sort by Skyline column so all Skyline-only appear together
        result_df = result_df.sort_values("Skyline", ascending=False).reset_index(drop=True)
        print("Mismatch summary (Files present in ONLY Skyline or ONLY SDRF):")
        print(result_df)
        return result_df
