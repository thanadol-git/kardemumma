import os
import re
from typing import List, Literal, Optional, Set

import pandas as pd

__all__ = [
    "CheckSkylineFile",
    "ImportFile",
    "MergeFiles",
    "cross_check_skyline_sdrf",
    "get_irt_peptides",
    "normalize_data_filename",
    "import_sdrf_file",
]

# Compiled once at import time
_ACQSUFFIX_RE = re.compile(r"_\d{8,14}$")
_HEAVY_RE = re.compile(r"heavy", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_path(path: str, ext: str) -> None:
    """Raise if *path* does not exist or does not end with *ext*."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    if not path.lower().endswith(ext):
        raise ValueError(f"Expected a {ext.upper()} file, got: {path}")


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# iRT peptide helpers
# ---------------------------------------------------------------------------
def _match_bionosys(peptide_list: List[str], ):
    """
    Print out suggestion that the iRT peptides are from Bionosys.
    Args:
        protein_list: List of protein names.
    Returns:
        None
    """ 
    
    # Make peptide list unique
    peptide_list = list(set(peptide_list))

    # Bionosys iRT peptides
    bionosys_list = ["LGGNEQVTR", 
                     "GAGSSEPVTGLDAK",
                     "VEATFGVDESNAK",
                     "YILAGVENSK",
                     "TPVISGGPYEYR",
                     "TPVITGAPYEYR",
                     "DGLDAASYYAPVR",
                     "ADVTPADFSEWSK",
                     "GTFIIDPGGVIR",
                     "GTFIIDPAAVIR",
                     "LFLQFGAQGSPFLK"]

    matched = []
    not_matched = []
    
    # Check if the peptide is in the Bionosys list
    for peptide in peptide_list:
        if peptide in bionosys_list:
            matched.append(peptide)
        else:
            not_matched.append(peptide)
    if matched:
        print(f"The following peptides are from Bionosys: {', '.join(matched)}.")
    if not_matched:
        print(f"The following peptides are not from Bionosys: {', '.join(not_matched)}.")
       
    
    # Print out the number of peptides that matched with Bionosys
    print(f"Number of peptides that did not match with Bionosys: {len(not_matched)}")

def get_irt_peptides(
    df: pd.DataFrame,
    *,
    protein_col: str = "Protein Name",
    peptide_col: str = "Peptide Sequence",
    tag_substring: str = "iRT_Tag",
    match_bionosys: bool = True,
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
    
    # Check if Protein and Peptide columns are in the dataframe
    if protein_col not in df.columns:
        raise ValueError(f"Protein column {protein_col} not found in dataframe.")
    if peptide_col not in df.columns:
        raise ValueError(f"Peptide column {peptide_col} not found in dataframe.")
    
    # Check if tag substring is in the protein column
    if tag_substring not in df[protein_col].astype(str).unique():
        raise ValueError(f"Tag substring {tag_substring} not found in protein column.")
    
    # Extract peptide list that contains the tag substring in protein column
    peptide_list = df.loc[df[protein_col].astype(str).str.contains(tag_substring, na=False, regex=False), peptide_col].drop_duplicates().tolist()
    
    # Match iRT peptides with Bionosys
    if match_bionosys:
        _match_bionosys(peptide_list)
    
    print("\n")
    print(f"Number of iRT peptides: {len(peptide_list)}")
    print(f"These are the iRT peptides: {peptide_list}")
    print("\n")
    print("iRT peptide analysis completed.")

    return peptide_list

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
    s = os.path.basename(str(name).strip().replace("\\", "/"))
    if not strip_acquisition_suffix:
        return s
    base, ext = os.path.splitext(s)
    if ext.lower() != ".raw":
        return s
    while _ACQSUFFIX_RE.search(base):
        base = _ACQSUFFIX_RE.sub("", base, count=1)
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

    Returns:
        ``True`` after printing a comparison report.
    """
    for label, df, col in (
        ("skyline", skyline_df, skyline_file_col),
        ("sdrf", sdrf_df, sdrf_file_col),
    ):
        if col not in df.columns:
            print(f"Missing required column {col!r} in {label}_df. Found: {list(df.columns)}")
            return True

    skyline_set = _normalized_file_set(
        skyline_df[skyline_file_col], strip_acquisition_suffix=strip_acquisition_suffix
    )
    sdrf_set = _normalized_file_set(
        sdrf_df[sdrf_file_col], strip_acquisition_suffix=strip_acquisition_suffix
    )

    if match_mode not in {"exact", "sdrf_in_skyline"}:
        print(f"Unknown match_mode {match_mode!r}; continuing with a generic set comparison report.")

    print("Skyline files:", len(skyline_set))
    print("SDRF files:", len(sdrf_set))
    print("Overlapping files:", len(skyline_set & sdrf_set))
    print("Only in Skyline:", sorted(skyline_set - sdrf_set))
    print("Only in SDRF:", sorted(sdrf_set - skyline_set))
    return True


def import_sdrf_file(file_path: str) -> pd.DataFrame:
    """Backward-compatible module-level wrapper for SDRF import."""
    return ImportFile(file_path).import_sdrf_file()


# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------

_SKYLINE_EXPECTED_COLS = [
    "Precursor",
    "Replicate",
    "File Name",
    "Peptide Retention Time",
    "Retention Time Calculator Score",
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


class ImportFile:
    def __init__(self, file_path: str):
        """
        Args:
            file_path: Path to the file to import.
        """
        self.file_path = file_path

    def import_skyline_file(self) -> pd.DataFrame:
        _validate_path(self.file_path, ".csv")
        df = pd.read_csv(self.file_path)
        if df.empty:
            raise ValueError(f"File is empty: {self.file_path}")
        missing = [c for c in _SKYLINE_EXPECTED_COLS if c not in df.columns]
        if missing:
            raise ValueError(
                f"Missing expected columns: {missing}. Actual: {list(df.columns)}"
            )
        # Vectorized heavy/light label — faster than row-wise apply + re.search
        df["Isotope Label Type"] = (
            df["Precursor"].astype(str).str.contains(_HEAVY_RE.pattern, case=False, regex=True)
            .map({True: "heavy", False: "light"})
        )
        print(f"The file {self.file_path} is valid.")
        return df

    def suggest_qc_samples(self, df: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Suggest QC samples from distinct ``Replicate`` values whose name contains
        ``qc`` (case-insensitive), e.g. ``Pool_QC``, ``qc-pool``.

        Args:
            df: Optional already-loaded Skyline table. If omitted, the CSV at
                ``self.file_path`` is read.

        Returns:
            Sorted list of replicate names.
        """
        if df is None:
            _validate_path(self.file_path, ".csv")
            df = pd.read_csv(self.file_path)
        suspicious = _suggest_qc_replicate_names(df)
        if suspicious:
            print(f"Possible QC samples: {suspicious}")
            print("Consider removing these before further analysis.")
        else:
            print("No suspicious replicates found.")
        return suspicious

    def import_qreps_file(self) -> pd.DataFrame:
        """Import the qREPs file."""
        _validate_path(self.file_path, ".csv")
        return pd.read_csv(self.file_path)

    def import_sdrf_file(self) -> pd.DataFrame:
        _validate_path(self.file_path, ".tsv")
        df = pd.read_csv(self.file_path, sep="\t")
        qc_samples = _suggest_qc_replicate_names(df, col_files="comment[data file]")
        print(f"QC samples from SDRF: {qc_samples}")
        print("Consider removing these before further analysis.")
        return df


class CheckSkylineFile:
    def __init__(self, file_path: str):
        """
        Args:
            file_path: Path to the file to check.
        """
        self.file_path = file_path

    def _load_csv(self) -> pd.DataFrame:
        _validate_path(self.file_path, ".csv")
        return pd.read_csv(self.file_path)

    def check_skyline_file(self) -> pd.DataFrame:
        """Check if the Skyline file is valid and return a sorted DataFrame."""
        return self._load_csv().sort_values(["Replicate", "Peptide"])

    def suggest_qc_samples(self, df: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Suggest QC samples: distinct ``Replicate`` values whose name contains ``qc``
        (case-insensitive), either from a given DataFrame or by reading from file.
        """
        if df is None:
            df = self._load_csv()
        return _suggest_qc_replicate_names(df)

    def get_irt_peptides(self) -> List[str]:
        """Load this CSV and return iRT peptide sequences."""
        return get_irt_peptides(self._load_csv())

    def get_test_samples(self, removed_samples: List[str], df: pd.DataFrame) -> List[str]:
        """
        Return sample names from ``df['Replicate']`` excluding ``removed_samples``.
        """
        return sorted(set(df["Replicate"].unique()) - set(removed_samples))

    def get_test_data(self, test_samples: List[str], df: pd.DataFrame) -> pd.DataFrame:
        """Filter ``df`` to ``test_samples`` and sort by Replicate and Peptide."""
        return (
            df[df["Replicate"].isin(test_samples)]
            .sort_values(["Replicate", "Peptide"])
        )


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
        """Merge Skyline with SDRF on replicate / source name and add ``characteristics[Plate]``."""
        skyline_df = (
            self.skyline_df[self.skyline_df["Peptide"].isin(self.selected_peptides)].copy()
            if self.selected_peptides is not None
            else self.skyline_df
        )

        required = {"skyline_df": (skyline_df, ["Replicate"]),
                    "sdrf_df": (self.sdrf_df, ["source name", "characteristics[Sample]"])}
        for label, (frame, cols) in required.items():
            missing = [c for c in cols if c not in frame.columns]
            if missing:
                raise ValueError(f"Missing column(s) {missing} in {label}")

        merged = pd.merge(
            skyline_df,
            self.sdrf_df[["source name", "characteristics[Sample]"]],
            left_on="Replicate",
            right_on="source name",
            how="left",
        )
        # expand=False guarantees a Series, no isinstance guard needed
        merged["characteristics[Plate]"] = merged["Replicate"].str.extract(
            r"Plate_(\d+)", expand=False
        )
        return merged

    def select_pool_data(
        self,
        df: pd.DataFrame,
        col_sample: str = "characteristics[Sample]",
        sample_value: str = "Pool",
    ) -> pd.DataFrame:
        """Select pool data from the DataFrame."""
        return (
            df[df[col_sample] == sample_value]
            .sort_values(["Replicate", "Peptide", "Isotope Label Type"])
            .reset_index(drop=True)
        )
