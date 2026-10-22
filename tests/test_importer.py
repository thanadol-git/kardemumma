"""
Tests for kardemumma.importer

Run with:  pytest tests/test_importer.py -v
"""

import pandas as pd
import pytest

from kardemumma.importer import ImportSkylineFile

_REQUIRED_COLS = [
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


def _write_skyline_csv(path, precursors, library_dot_products):
    """Write a minimal, valid Skyline-shaped CSV. One row per entry."""
    n = len(precursors)
    df = pd.DataFrame(
        {
            "Precursor": precursors,
            "Replicate": [f"R{i}" for i in range(n)],
            "File Name": [f"f{i}.raw" for i in range(n)],
            "Peptide Retention Time": [10.0] * n,
            "Predicted Retention Time": [10.0] * n,
            "Precursor Charge": [2] * n,
            "Peptide Sequence": ["PEPTIDE"] * n,
            "Peptide": ["PEPTIDE"] * n,
            "Normalized Area": [1000.0] * n,
            "RatioLightToHeavy": [1.0] * n,
            "Ratio Dot Product": [0.9] * n,
            "Library Dot Product": library_dot_products,
            "Protein Name": ["ProtA"] * n,
        }
    )
    df.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# ImportSkylineFile
# ---------------------------------------------------------------------------

class TestImportSkylineFile:
    def test_derives_isotope_label_type(self, tmp_path):
        path = tmp_path / "skyline.csv"
        _write_skyline_csv(
            path,
            precursors=["PEP_2 (heavy)", "PEP_2"],
            library_dot_products=[0.9, 0.85],
        )
        df = ImportSkylineFile(str(path)).import_skyline_file()
        assert list(df["Isotope Label Type"]) == ["heavy", "light"]

    def test_warns_when_one_label_type_has_no_dot_product_data(self, tmp_path, caplog):
        """Regression test for a real export we hit: Kardemumma_v4_oct05.csv
        had 'Library Dot Product' entirely empty for every light-labeled row
        (heavy was fully populated). The data imported without error and the
        problem only surfaced later as a silently-empty plot. This pins down
        that import_skyline_file now warns immediately, at import time.
        """
        path = tmp_path / "skyline.csv"
        _write_skyline_csv(
            path,
            precursors=["PEP_2 (heavy)", "PEP_2 (heavy)", "PEP_2", "PEP_2"],
            library_dot_products=[0.9, 0.92, None, None],
        )
        with caplog.at_level("WARNING", logger="kardemumma.importer"):
            ImportSkylineFile(str(path)).import_skyline_file()
        assert any(
            "light" in record.message and "Library Dot Product" in record.message
            for record in caplog.records
        )

    def test_no_warning_when_both_label_types_have_data(self, tmp_path, caplog):
        path = tmp_path / "skyline.csv"
        _write_skyline_csv(
            path,
            precursors=["PEP_2 (heavy)", "PEP_2 (heavy)", "PEP_2", "PEP_2"],
            library_dot_products=[0.9, 0.92, 0.7, 0.75],
        )
        with caplog.at_level("WARNING", logger="kardemumma.importer"):
            ImportSkylineFile(str(path)).import_skyline_file()
        assert caplog.records == []

    def test_missing_expected_column_raises(self, tmp_path):
        path = tmp_path / "skyline.csv"
        pd.DataFrame({"Precursor": ["PEP_2"]}).to_csv(path, index=False)
        with pytest.raises(ValueError):
            ImportSkylineFile(str(path)).import_skyline_file()

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            ImportSkylineFile(str(tmp_path / "does_not_exist.csv")).import_skyline_file()

    def test_wrong_extension_raises(self, tmp_path):
        path = tmp_path / "skyline.tsv"
        _write_skyline_csv(path, precursors=["PEP_2"], library_dot_products=[0.9])
        with pytest.raises(ValueError):
            ImportSkylineFile(str(path)).import_skyline_file()
