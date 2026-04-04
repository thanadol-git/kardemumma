import os
import re
from typing import List, Optional

import pandas as pd
from .sdrf import validate_sdrf as validate_sdrf_file


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
