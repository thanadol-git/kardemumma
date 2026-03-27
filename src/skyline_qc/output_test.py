import os
import pandas as pd
from .sdrf import validate_sdrf as validate_sdrf_file

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
        return df

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