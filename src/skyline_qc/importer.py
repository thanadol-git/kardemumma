import os
import re
from typing import List, Literal, Optional, Set

import pandas as pd
from .sdrf import validate_sdrf as validate_sdrf_file

__all__ = [
    "CheckSkylineFile",
    "ImportFile",
    "MergeFiles",
    "cross_check_skyline_sdrf",
    "get_irt_peptides",
    "normalize_data_filename",
]


def get_irt_peptides(
    df: pd.DataFrame,
    *,
    protein_col: str = "Protein Name",
    peptide_col: str = "Peptide Sequence",
    tag_substring: str = "iRT_Tag",
) -> List[str]:
    """
    Unique peptide sequences for rows whose protein name contains the iRT tag
    marker (default substring ``iRT_Tag``), in first-seen order.

    Args:
        df: Skyline report table (e.g. from ``ImportFile.import_skyline_file``).
        protein_col: Column with protein / group name.
        peptide_col: Column with peptide sequence.
        tag_substring: Literal substring to match in ``protein_col`` (not a regex).

    Returns:
        List of distinct peptide sequences.
    """
    missing = [c for c in (protein_col, peptide_col) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing column(s) {missing}. Found: {list(df.columns)}"
        )
    mask = (
        df[protein_col]
        .astype(str)
        .str.contains(tag_substring, na=False, regex=False)
    )
    return df.loc[mask, peptide_col].drop_duplicates().tolist()


def _suggest_qc_replicate_names(df: pd.DataFrame) -> List[str]:
    """
    Distinct Replicate labels that look like QC (substring 'qc', case-insensitive).
    """
    if "Replicate" not in df.columns:
        raise ValueError(
            f"Expected column 'Replicate'. Found: {list(df.columns)}"
        )
    reps = df["Replicate"].dropna().astype(str)
    return sorted({r for r in reps.unique() if "qc" in r.lower()})


def normalize_data_filename(
    name: str,
    *,
    strip_acquisition_suffix: bool = True,
) -> str:
    """
    Canonical form for comparing raw names between Skyline and SDRF.

    - Uses basename only (strips any path / URI path segments).
    - For ``.raw`` / ``.RAW``, optionally removes one or more trailing
      ``_`` + 8–14 digit blocks (common acquisition / replicate suffixes),
      e.g. ``..._Plate_5_C6_20251211114043.raw`` → ``..._Plate_5_C6.raw``.
    """
    s = str(name).strip().replace("\\", "/")
    s = os.path.basename(s)
    if not strip_acquisition_suffix:
        return s
    base, ext = os.path.splitext(s)
    if ext.lower() != ".raw":
        return s
    pat = re.compile(r"_\d{8,14}$")
    while pat.search(base):
        base = pat.sub("", base, count=1)
    return f"{base}{ext}"


def _normalized_file_set(
    series: pd.Series,
    *,
    strip_acquisition_suffix: bool,
) -> Set[str]:
    return {
        normalize_data_filename(v, strip_acquisition_suffix=strip_acquisition_suffix)
        for v in series.dropna().astype(str).unique()
    }


def cross_check_skyline_sdrf(
    skyline_df: pd.DataFrame,
    sdrf_df: pd.DataFrame,
    *,
    skyline_file_col: str = "File Name",
    sdrf_file_col: str = "comment[data file]",
    strip_acquisition_suffix: bool = True,
    match_mode: Literal["exact", "sdrf_in_skyline"] = "sdrf_in_skyline",
) -> bool:
    """
    Cross-check Skyline ``File Name`` values vs SDRF ``comment[data file]``.

    By default, names are normalized (basename + strip trailing ``_########`` before
    ``.raw``) and the check is **``sdrf_in_skyline``**: every distinct SDRF file must
    appear in Skyline; extra Skyline files (e.g. repeated QC acquisitions) are allowed.

    Use ``match_mode="exact"`` for strict set equality after normalization.
    Use ``strip_acquisition_suffix=False`` to compare strings verbatim.

    Raises:
        KeyError: If required columns are missing.
        ValueError: If the check fails (see message for only-Skyline / only-SDRF lists).

    Returns:
        ``True`` if the check passes.
    """
    for name, df, col in (
        ("skyline", skyline_df, skyline_file_col),
        ("sdrf", sdrf_df, sdrf_file_col),
    ):
        if col not in df.columns:
            raise KeyError(f"{name}_df: missing column {col!r}. Found: {list(df.columns)}")

    skyline_set = _normalized_file_set(
        skyline_df[skyline_file_col],
        strip_acquisition_suffix=strip_acquisition_suffix,
    )
    sdrf_set = _normalized_file_set(
        sdrf_df[sdrf_file_col],
        strip_acquisition_suffix=strip_acquisition_suffix,
    )

    if match_mode == "exact":
        if skyline_set != sdrf_set:
            only_skyline = sorted(skyline_set - sdrf_set)
            only_sdrf = sorted(sdrf_set - skyline_set)
            print("❌ File sets do not match after normalization.")
            raise ValueError(
                "Skyline vs SDRF raw file names do not match (exact, after normalization). "
                f"Only in Skyline: {only_skyline}. Only in SDRF: {only_sdrf}."
            )
        print("✅ File sets match exactly after normalization! The file is correct!")
        return True

    if match_mode == "sdrf_in_skyline":
        missing_in_skyline = sorted(sdrf_set - skyline_set)
        if missing_in_skyline:
            extra_skyline = sorted(skyline_set - sdrf_set)
            print("❌ Some SDRF files are missing from Skyline.")
            raise ValueError(
                "After normalization, some SDRF raw files are missing from Skyline "
                f"File Name: {missing_in_skyline}. "
                f"(Skyline-only extras, allowed in this mode: {extra_skyline})"
            )
        print("✅ All SDRF files are present in Skyline after normalization! The file is correct!")
       
        return True

    raise ValueError(f"Unknown match_mode: {match_mode!r}")


class ImportFile():
    def __init__(self, file_path: str):
        """
        Args:
            file_path: Path to the file to import.
        """
        self.file_path = file_path
    ### Skyline file ###
    def import_skyline_file(self):
        
        # 1. Check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"The file {self.file_path} does not exist.")
        
        # 2. Check file type
        if not self.file_path.lower().endswith('.csv'):
            raise ValueError("Provided file is not a CSV.")

        # 3. Read the file
        df = pd.read_csv(self.file_path)

        # 4. Check if the file is empty
        if df.empty:
            raise ValueError(f"The file {self.file_path} is empty.")

        # 5.  Check if the file has the expected columns
        expected_columns = [
            "Precursor",
            "Replicate",
            "File Name",
            "Peptide Retention Time",
            "Retention Time Calculator Score",
            "Predicted Retention Time",
            "Peptide Sequence",
            "Peptide",
            "Normalized Area",
            "RatioLightToHeavy",
            "Ratio Dot Product",
            "Library Dot Product",
            "Protein Name"
        ]
        
        # 6. Check if the file has the expected columns
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(
                f"The file {self.file_path} is missing expected columns: {missing_columns}. "
                f"Actual columns: {list(df.columns)}"
            )

        df = df.copy()
        df["Isotope Label Type"] = df["Precursor"].apply(
            lambda x: "heavy"
            if re.search(r"heavy", str(x), re.IGNORECASE)
            else "light"
        )

        # Print if the file is valid
        print(f"The file {self.file_path} is valid.")
        return df

    def suggest_qc_samples(self, df: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Suggest QC samples from distinct ``Replicate`` values whose name contains
        ``qc`` (case-insensitive), e.g. ``Pool_QC``, ``qc-pool``.

        Args:
            df: Optional already-loaded Skyline table. If omitted, the CSV at
                ``self.file_path`` is read (must exist and be ``.csv``).

        Returns:
            Sorted list of replicate names.
        """
        if df is None:
            if not os.path.exists(self.file_path):
                raise FileNotFoundError(f"The file {self.file_path} does not exist.")
            if not self.file_path.lower().endswith(".csv"):
                raise ValueError("Provided file is not a CSV.")
            df = pd.read_csv(self.file_path)
        return _suggest_qc_replicate_names(df)
    ### qREPs file ###
    def import_qreps_file(self):
        """
        Import the qREPs file.
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"The file {self.file_path} does not exist.")
        if not self.file_path.lower().endswith(".csv"):
            raise ValueError("Provided file is not a CSV.")
        df = pd.read_csv(self.file_path)
        return df

    ### SDRF file ###
    def import_sdrf_file(self):
        # 1. Check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"The file {self.file_path} does not exist.")
        
        # 2. Check file type
        if not self.file_path.lower().endswith('.tsv'):
            raise ValueError("Provided file is not a TSV.")
        
        # 3. Read the file
        df = pd.read_csv(self.file_path, sep='\t')

        # 4. Validate the SDRF file (via sdrf-pipelines CLI)
        ok, msg = validate_sdrf_file(self.file_path)
        if not ok:
            raise ValueError(f"SDRF validation failed: {msg}")
        return df

class CheckSkylineFile():
    def __init__(self, file_path: str):
        """
        Args:
            file_path: Path to the file to check.
        """
        self.file_path = file_path
    
    def check_skyline_file(self):
        """
        Check if the Skyline file is valid and return a DataFrame.
        """
        # 1. Check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"The file {self.file_path} does not exist.")
        
        # 2. Check file type
        if not self.file_path.lower().endswith('.csv'):
            raise ValueError("Provided file is not a CSV.")
        
        # 3. Read the file
        df = pd.read_csv(self.file_path)

        # Sort by Replicate and Peptide
        df = df.sort_values(by=['Replicate', 'Peptide'])
        return df

    def suggest_qc_samples(self) -> List[str]:
        """
        Suggest QC samples: distinct ``Replicate`` values whose name contains ``qc``
        (case-insensitive). Reads the CSV at ``self.file_path``.
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"The file {self.file_path} does not exist.")
        if not self.file_path.lower().endswith(".csv"):
            raise ValueError("Provided file is not a CSV.")
        df = pd.read_csv(self.file_path)
        return _suggest_qc_replicate_names(df)

    def get_qc_data(self, qc_samples: List[str], df: pd.DataFrame) -> pd.DataFrame:
        """
        Get QC data from the Skyline data object.
        """
        return df[df['Replicate'].isin(qc_samples)]

    def get_irt_peptides(self) -> List[str]:
        """
        Load this CSV and return iRT peptide sequences (same rules as
        :func:`get_irt_peptides`).
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"The file {self.file_path} does not exist.")
        if not self.file_path.lower().endswith(".csv"):
            raise ValueError("Provided file is not a CSV.")
        df = pd.read_csv(self.file_path)
        return get_irt_peptides(df)

    def get_test_samples(self, qc_samples: List[str], df: pd.DataFrame) -> List[str]:
        """
        Return a list of test sample names by excluding QC samples from the distinct Replicate names in the given DataFrame.
        
        Args:
            qc_samples: List of QC sample names to exclude.
            df: DataFrame containing Skyline data, must include a 'Replicate' column.
        
        Returns:
            List of unique test sample names.
        """
        all_samples = set(df["Replicate"].unique())
        test_samples = sorted(list(all_samples - set(qc_samples)))
        return test_samples

    def get_test_data(self, test_samples: List[str], df: pd.DataFrame) -> pd.DataFrame:
        """
        Get test data from the Skyline data object and list of test samples. It basically remove test samples from the skyline data object.
        """
        df = df[df['Replicate'].isin(test_samples)]
        # Sort by Replicate and Peptide
        df = df.sort_values(by=['Replicate', 'Peptide'])
        return df

class MergeFiles:
    def __init__(
        self,
        skyline_df: pd.DataFrame,
        sdrf_df: pd.DataFrame,
        selected_peptides: Optional[List[str]] = None,
    ):
        """
        Args:
            skyline_df: Skyline long-format report (must include ``Replicate``, ``Peptide``).
            sdrf_df: SDRF table (must include ``source name``, ``characteristics[Sample]``).
            selected_peptides: If given, restrict ``skyline_df`` to these ``Peptide`` values
                before merging. If ``None``, all peptides are kept.
        """
        self.skyline_df = skyline_df
        self.sdrf_df = sdrf_df
        self.selected_peptides = selected_peptides

    def merge_files(self) -> pd.DataFrame:
        """
        Merge Skyline with SDRF on replicate / source name and add ``characteristics[Plate]``.
        """
        if self.selected_peptides is not None:
            self.skyline_df = self.skyline_df[
                self.skyline_df["Peptide"].isin(self.selected_peptides)
            ].copy()

        # 1. Check that both DataFrames contain the required columns before merging
        required_skyline_cols = ['Replicate']
        required_sdrf_cols = ['source name', 'characteristics[Sample]']
        for col in required_skyline_cols:
            if col not in self.skyline_df.columns:
                raise ValueError(f"Missing column '{col}' in skyline_df")
        for col in required_sdrf_cols:
            if col not in self.sdrf_df.columns:
                raise ValueError(f"Missing column '{col}' in sdrf_df")

        # 2. Merge skyline_df with sdrf_df using 'Replicate' from skyline and 'comment[data file]' from sdrf_df
        skyline_merge = pd.merge(self.skyline_df, self.sdrf_df[['source name', 'characteristics[Sample]']], left_on='Replicate', right_on='source name', how='left')

        # 3. Plate id from Replicate (digits after "Plate_"). Use expand=False so this is
        #    a scalar Series; str.extract(..., expand=True) returns a DataFrame and can
        #    break statsmodels/patsy when used as a categorical column.
        plate = skyline_merge["Replicate"].str.extract(r"Plate_(\d+)", expand=False)
        if isinstance(plate, pd.DataFrame):
            plate = plate.squeeze(axis=1)
        skyline_merge["characteristics[Plate]"] = plate
        return skyline_merge