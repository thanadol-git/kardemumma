"""
Tests for kardemumma.teapot

Run with:  pytest tests/test_teapot.py -v
"""

import pandas as pd
import pytest

from kardemumma.importer import MergeFiles
from kardemumma.prm import summarise_peptide_counts
from kardemumma.teapot import ImportTeapotFile, find_teapot_confident_table, pivot_confident_pairs

_PROTEIN = {"PEPAK": "QR0000001_ABC|P11111", "PEPCR": "QR0000002_DEF|P22222"}


def _pair(sample, peptide, charge, light=100.0, heavy=50.0, peptidoform=None):
    ratio = light / heavy
    return {
        "sample_id": sample, "protein": _PROTEIN[peptide], "peptide": peptide,
        "peptidoform": peptidoform or peptide, "charge": charge,
        "light_abundance": light, "heavy_abundance": heavy,
        "ratio_L_H": ratio, "ratio_H_L": 1 / ratio,
        "heavy_pmol_known": 0.5, "abs_pmol_light": ratio * 0.5,
        "light_q_value": 0.001, "heavy_q_value": 0.002,
        "light_rt_seconds": 600.0, "heavy_rt_seconds": 603.0,
    }


@pytest.fixture()
def confident_table() -> pd.DataFrame:
    return pd.DataFrame([
        _pair("Plate_1_A1", "PEPAK", 2),
        _pair("Plate_1_A1", "PEPCR", 3, light=30.0, heavy=60.0, peptidoform="PEPC[+57.0]R"),
        _pair("Plate_2_A1", "PEPAK", 2, light=120.0),
    ])


@pytest.fixture()
def teapot_outdir(tmp_path, confident_table):
    """A teapot --outdir with one route holding the confident table."""
    route = tmp_path / "results" / "quant" / "enc_dia"
    route.mkdir(parents=True)
    confident_table.to_csv(route / "peptide.abs_quant.confident.tsv", sep="\t", index=False)
    return tmp_path / "results"


@pytest.fixture()
def sdrf_df() -> pd.DataFrame:
    return pd.DataFrame({
        "source name": ["sample 1", "sample 2"],
        "characteristics[plate]": ["Plate_1", "Plate_2"],
        "characteristics[Sample]": ["Pool", "Pool"],
        "comment[data file]": ["Plate_1_A1.raw", "Plate_2_A1.raw"],
        "factor value[Sample]": ["Pool", "Pool"],
    })


class TestFindTeapotConfidentTable:
    def test_outdir_with_one_route(self, teapot_outdir):
        path = find_teapot_confident_table(str(teapot_outdir))
        assert path == teapot_outdir / "quant" / "enc_dia" / "peptide.abs_quant.confident.tsv"

    def test_file_itself(self, teapot_outdir):
        f = teapot_outdir / "quant" / "enc_dia" / "peptide.abs_quant.confident.tsv"
        assert find_teapot_confident_table(str(f)) == f

    def test_several_routes_need_a_route(self, teapot_outdir, confident_table):
        other = teapot_outdir / "quant" / "enc_prm"
        other.mkdir()
        confident_table.to_csv(other / "peptide.abs_quant.confident.tsv", sep="\t", index=False)
        with pytest.raises(FileNotFoundError, match="choose one with route"):
            find_teapot_confident_table(str(teapot_outdir))
        assert find_teapot_confident_table(str(teapot_outdir), route="enc_prm").parent.name == "enc_prm"

    def test_unknown_route_raises(self, teapot_outdir):
        with pytest.raises(FileNotFoundError, match="enc_prm"):
            find_teapot_confident_table(str(teapot_outdir), route="enc_prm")


class TestImportTeapotFile:
    def test_one_light_and_one_heavy_row_per_pair(self, teapot_outdir):
        importer = ImportTeapotFile(str(teapot_outdir))
        df = importer.import_teapot_file()
        assert importer.route == "enc_dia"
        assert (df["Isotope Label Type"].value_counts() == 3).all()
        heavy = df[df["Isotope Label Type"] == "heavy"]
        assert heavy["RatioLightToHeavy"].isna().all()
        assert set(heavy["Precursor"]) == {"PEPAK++ (heavy)", "PEPC[+57.0]R+++ (heavy)"}

    def test_columns_kardemumma_reads(self, teapot_outdir):
        df = ImportTeapotFile(str(teapot_outdir)).import_teapot_file()
        light = df[df["Isotope Label Type"] == "light"].set_index(["Replicate", "Peptide"])
        row = light.loc[("Plate_1_A1", "PEPC[+57.0]R")]
        assert row["Peptide Sequence"] == "PEPCR"
        assert row["RatioLightToHeavy"] == pytest.approx(0.5)
        assert row["Normalized Area"] == pytest.approx(30.0)
        assert row["Peptide Retention Time"] == pytest.approx(10.0)
        assert row["Q Value"] == pytest.approx(0.001)
        assert row["Heavy Amount [pmol]"] == pytest.approx(0.5)
        assert light.loc[("Plate_2_A1", "PEPAK"), "RatioLightToHeavy"] == pytest.approx(2.4)
        heavy = df[df["Isotope Label Type"] == "heavy"]
        assert heavy["Peptide Retention Time"].iloc[0] == pytest.approx(10.05)
        assert heavy["Q Value"].iloc[0] == pytest.approx(0.002)

    def test_not_the_confident_table_raises(self, tmp_path):
        f = tmp_path / "peptide.abs_quant.confident.tsv"
        f.write_text("sample_id\tprotein\nS1\tP\n")
        with pytest.raises(ValueError, match="not teapot's"):
            ImportTeapotFile(str(f)).import_teapot_file()

    def test_without_sdrf_file_name_is_sample_id(self, teapot_outdir):
        df = ImportTeapotFile(str(teapot_outdir)).import_teapot_file()
        assert set(df["File Name"]) == {"Plate_1_A1", "Plate_2_A1"}

    def test_sdrf_gives_data_files(self, teapot_outdir, sdrf_df):
        df = ImportTeapotFile(str(teapot_outdir)).import_teapot_file(sdrf_df=sdrf_df)
        assert set(df["File Name"]) == {"Plate_1_A1.raw", "Plate_2_A1.raw"}
        merged = MergeFiles(df, sdrf_df).merge_files()
        assert merged["characteristics[plate]"].notna().all()

    def test_runs_the_sdrf_does_not_list_are_left_out(self, teapot_outdir, sdrf_df, capsys):
        df = ImportTeapotFile(str(teapot_outdir)).import_teapot_file(sdrf_df=sdrf_df.head(1))
        assert set(df["Replicate"]) == {"Plate_1_A1"}
        assert "Plate_2_A1" in capsys.readouterr().out

    def test_no_run_in_the_sdrf_raises(self, teapot_outdir, sdrf_df):
        sdrf = sdrf_df.assign(**{"comment[data file]": ["x.raw", "y.raw"]})
        with pytest.raises(ValueError, match="No teapot sample_id is an SDRF data file"):
            ImportTeapotFile(str(teapot_outdir)).import_teapot_file(sdrf_df=sdrf)


class TestPivotConfidentPairs:
    def test_counts_runs_with_a_confident_pair(self, teapot_outdir):
        df = ImportTeapotFile(str(teapot_outdir)).import_teapot_file()
        wide = pivot_confident_pairs(df)
        assert list(wide.columns) == ["Replicate", "Protein Name", "Peptide", "heavy", "light"]
        counts = summarise_peptide_counts(wide).set_index("Peptide")
        assert counts.loc["PEPAK", ["heavy_count", "light_count"]].tolist() == [2, 2]
        assert counts.loc["PEPC[+57.0]R", ["heavy_count", "light_count"]].tolist() == [1, 1]
