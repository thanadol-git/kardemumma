import os
import unittest
import pandas as pd
from .sdrf import validate_sdrf as validate_sdrf_file

class ImportFile():
    def __init__(self, file_path: str):
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
        if not all(col in df.columns for col in expected_columns):
            print(f"The file {self.file_path} does not have the expected columns. Expected:")
            for col in expected_columns:
                print(col)
            raise ValueError(f"The file {self.file_path} does not have the expected columns.")
        
        # 7. Cross-check expected columns with actual column names from the file.
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            print(f"The following expected columns are missing from the file: {missing_columns}")
            print(f"Actual columns in the file: {list(df.columns)}")
            raise ValueError(f"The file {self.file_path} does not have the expected columns: {expected_columns}")
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