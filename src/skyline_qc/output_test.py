import os
import re
from typing import List, Optional

import pandas as pd
from .sdrf import validate_sdrf as validate_sdrf_file


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