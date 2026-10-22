"""
Tests for skyline_qc.prm

Run with:  pytest tests/test_prm.py -v
"""

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest

from kardemumma.prm import (
    adjust_ratio_by_plate,
    compute_cv,
    dot_product_summary,
    filter_library_dot_product,
    flag_missing_values,
    get_absolute_conc,
    get_plate_conversion_factors,
    plot_all_all,
    retention_time_deviation,
    summarise_peptide_counts,
    summarize_prm,
)

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

@pytest.fixture()
def sample_df() -> pd.DataFrame:
    """Minimal Skyline-style DataFrame with three precursors and three replicates."""
    return pd.DataFrame(
        {
            "Precursor": ["PEP1_2", "PEP1_2", "PEP1_2",
                          "PEP2_2", "PEP2_2", "PEP2_2",
                          "PEP3_3", "PEP3_3", "PEP3_3"],
            "Replicate": ["R1", "R2", "R3",
                          "R1", "R2", "R3",
                          "R1", "R2", "R3"],
            "File Name": ["f1.raw", "f2.raw", "f3.raw"] * 3,
            "Peptide Sequence": ["PEPTIDE", "PEPTIDE", "PEPTIDE",
                                 "SEQUENCE", "SEQUENCE", "SEQUENCE",
                                 "ANOTHER", "ANOTHER", "ANOTHER"],
            "Peptide": ["PEPTIDE", "PEPTIDE", "PEPTIDE",
                        "SEQUENCE", "SEQUENCE", "SEQUENCE",
                        "ANOTHER", "ANOTHER", "ANOTHER"],
            "Protein Name": ["ProtA"] * 3 + ["ProtB"] * 3 + ["ProtC"] * 3,
            "Normalized Area": [1000.0, 1100.0, 950.0,
                                500.0, 490.0, 510.0,
                                800.0, 820.0, 780.0],
            # PEP1: low CV; PEP2: high CV; PEP3: perfect (constant)
            "RatioLightToHeavy": [1.0, 1.1, 0.9,
                                  1.0, 2.0, 3.0,
                                  2.0, 2.0, 2.0],
            "Library Dot Product": [0.95, 0.90, 0.85,
                                    0.75, 0.70, 0.65,
                                    0.92, 0.88, 0.91],
            "Ratio Dot Product":   [0.93, 0.91, 0.89,
                                    0.60, 0.55, 0.50,
                                    0.94, 0.96, 0.95],
            "Peptide Retention Time":  [10.0, 10.2, 9.8,
                                        20.0, 20.1, 19.9,
                                        15.0, 15.5, 14.5],
            "Predicted Retention Time": [10.0, 10.0, 10.0,
                                         20.0, 20.0, 20.0,
                                         15.0, 15.0, 15.0],
            "Retention Time Calculator Score": [0.99] * 9,
        }
    )


@pytest.fixture()
def sample_long_df(sample_df) -> pd.DataFrame:
    """Long-format PRM table with light/heavy rows for filter tests."""
    heavy = sample_df.copy()
    heavy["Precursor"] = heavy["Precursor"] + " (heavy)"
    heavy["Isotope Label Type"] = "heavy"
    light = sample_df.copy()
    light["Isotope Label Type"] = "light"
    return pd.concat([heavy, light], ignore_index=True)


# ---------------------------------------------------------------------------
# filter_library_dot_product
# ---------------------------------------------------------------------------

class TestFilterLibraryDotProduct:
    def test_pivots_to_heavy_light_columns(self, sample_long_df):
        result = filter_library_dot_product(sample_long_df, threshold=0.6)
        assert list(result.columns) == [
            "Replicate", "Protein Name", "Peptide", "heavy", "light",
        ]
        assert len(result) == 9

    def test_coerces_string_normalized_area(self, sample_long_df):
        df = sample_long_df.copy()
        df["Normalized Area"] = df["Normalized Area"].astype(str)
        result = filter_library_dot_product(df, threshold=0.6)
        counts = summarise_peptide_counts(result)
        assert not counts.empty

    def test_single_channel_input_has_no_other_channel_column(self, sample_long_df):
        light_only = sample_long_df[sample_long_df["Isotope Label Type"] == "light"]
        result = filter_library_dot_product(light_only, threshold=0.6)
        assert list(result.columns) == ["Replicate", "Protein Name", "Peptide", "light"]
        assert "heavy" not in result.columns

    def test_empty_result_when_no_rows_pass(self, sample_long_df):
        result = filter_library_dot_product(sample_long_df, threshold=0.99)
        assert result.empty

    def test_rows_with_missing_normalized_area_are_dropped(self):
        df = pd.DataFrame(
            {
                "Precursor": ["PEP_2", "PEP_2 (heavy)"],
                "Replicate": ["R1", "R1"],
                "Protein Name": ["ProtA", "ProtA"],
                "Peptide": ["PEP", "PEP"],
                "Isotope Label Type": ["light", "heavy"],
                "Library Dot Product": [0.95, 0.95],
                "Normalized Area": [pd.NA, pd.NA],
                "Total Area": [1000.0, 2000.0],
            }
        )
        result = filter_library_dot_product(df, threshold=0.6)
        assert result.empty


# ---------------------------------------------------------------------------
# compute_cv
# ---------------------------------------------------------------------------

class TestComputeCV:
    def test_returns_one_row_per_precursor(self, sample_df):
        result = compute_cv(sample_df)
        assert len(result) == 3

    def test_columns_present(self, sample_df):
        result = compute_cv(sample_df)
        assert {"Precursor", "mean", "std", "cv_pct", "n"}.issubset(result.columns)

    def test_zero_cv_for_constant_values(self, sample_df):
        result = compute_cv(sample_df)
        pep3_row = result[result["Precursor"] == "PEP3_3"].iloc[0]
        assert pep3_row["cv_pct"] == pytest.approx(0.0, abs=1e-6)

    def test_sorted_ascending_by_cv(self, sample_df):
        result = compute_cv(sample_df)
        assert list(result["cv_pct"]) == sorted(result["cv_pct"])

    def test_high_cv_for_variable_precursor(self, sample_df):
        result = compute_cv(sample_df)
        pep2_cv = result[result["Precursor"] == "PEP2_2"]["cv_pct"].iloc[0]
        assert pep2_cv >= 50

    def test_custom_value_col(self, sample_df):
        result = compute_cv(sample_df, value_col="Normalized Area")
        assert "cv_pct" in result.columns

    def test_missing_column_raises(self, sample_df):
        with pytest.raises(KeyError):
            compute_cv(sample_df, value_col="NonExistentColumn")

    def test_missing_group_col_raises(self, sample_df):
        with pytest.raises(KeyError):
            compute_cv(sample_df, group_by="NonExistentColumn")


# ---------------------------------------------------------------------------
# flag_missing_values
# ---------------------------------------------------------------------------

class TestFlagMissingValues:
    def test_full_detection_rate(self, sample_df):
        result = flag_missing_values(sample_df)
        assert (result["detection_rate"] == 1.0).all()

    def test_detection_rate_with_nans(self, sample_df):
        df = sample_df.copy()
        # Make PEP2_2 / R1 missing
        df.loc[(df["Precursor"] == "PEP2_2") & (df["Replicate"] == "R1"),
               "RatioLightToHeavy"] = float("nan")
        result = flag_missing_values(df)
        pep2_rate = result[result["Precursor"] == "PEP2_2"]["detection_rate"].iloc[0]
        assert pep2_rate == pytest.approx(2 / 3)

    def test_zero_treated_as_missing(self, sample_df):
        df = sample_df.copy()
        df.loc[(df["Precursor"] == "PEP1_2") & (df["Replicate"] == "R2"),
               "RatioLightToHeavy"] = 0
        result = flag_missing_values(df)
        pep1_rate = result[result["Precursor"] == "PEP1_2"]["detection_rate"].iloc[0]
        assert pep1_rate == pytest.approx(2 / 3)

    def test_sorted_ascending_by_detection_rate(self, sample_df):
        df = sample_df.copy()
        df.loc[(df["Precursor"] == "PEP2_2") & (df["Replicate"] == "R1"),
               "RatioLightToHeavy"] = float("nan")
        result = flag_missing_values(df)
        rates = result["detection_rate"].tolist()
        assert rates == sorted(rates)

    def test_missing_column_raises(self, sample_df):
        with pytest.raises(KeyError):
            flag_missing_values(sample_df, value_col="NoSuchCol")


# ---------------------------------------------------------------------------
# dot_product_summary
# ---------------------------------------------------------------------------

class TestDotProductSummary:
    def test_all_pass_above_threshold(self, sample_df):
        # PEP3 rows all have lib≥0.88 and ratio≥0.94 — should all pass at 0.8
        result = dot_product_summary(sample_df)
        pep3 = result[result["Precursor"] == "PEP3_3"]
        assert pep3["both_pass"].all()

    def test_fails_below_threshold(self, sample_df):
        # PEP2 lib dot products are 0.75/0.70/0.65 — fail at default 0.8
        result = dot_product_summary(sample_df)
        pep2 = result[result["Precursor"] == "PEP2_2"]
        assert not pep2["lib_dot_pass"].any()

    def test_custom_threshold(self, sample_df):
        result = dot_product_summary(sample_df, lib_threshold=0.5, ratio_threshold=0.5)
        assert result["both_pass"].all()

    def test_columns_present(self, sample_df):
        result = dot_product_summary(sample_df)
        assert {"lib_dot_pass", "ratio_dot_pass", "both_pass"}.issubset(result.columns)

    def test_missing_column_raises(self, sample_df):
        with pytest.raises(KeyError):
            dot_product_summary(sample_df, lib_dot_col="NoSuchCol")


# ---------------------------------------------------------------------------
# retention_time_deviation
# ---------------------------------------------------------------------------

class TestRetentionTimeDeviation:
    def test_zero_deviation_when_observed_equals_predicted(self, sample_df):
        df = sample_df.copy()
        df["Peptide Retention Time"] = df["Predicted Retention Time"]
        result = retention_time_deviation(df)
        assert (result["abs_rt_dev"] == 0.0).all()

    def test_deviation_values(self, sample_df):
        result = retention_time_deviation(sample_df)
        # PEP3_3 R2: observed=15.5, predicted=15.0 → dev=0.5
        row = result[
            (result["Precursor"] == "PEP3_3") & (result["Replicate"] == "R2")
        ].iloc[0]
        assert row["rt_dev"] == pytest.approx(0.5)
        assert row["abs_rt_dev"] == pytest.approx(0.5)

    def test_sorted_descending_by_abs_rt_dev(self, sample_df):
        result = retention_time_deviation(sample_df)
        devs = result["abs_rt_dev"].tolist()
        assert devs == sorted(devs, reverse=True)

    def test_columns_present(self, sample_df):
        result = retention_time_deviation(sample_df)
        assert {"rt_observed", "rt_predicted", "rt_dev", "abs_rt_dev"}.issubset(
            result.columns
        )

    def test_missing_column_raises(self, sample_df):
        with pytest.raises(KeyError):
            retention_time_deviation(sample_df, observed_col="NoSuchCol")


# ---------------------------------------------------------------------------
# summarize_prm
# ---------------------------------------------------------------------------

class TestSummarizePRM:
    def test_returns_all_keys(self, sample_df):
        report = summarize_prm(sample_df)
        expected_keys = {
            "cv", "missing", "dot_products", "rt_deviation",
            "n_precursors", "n_replicates",
            "pct_dot_pass", "pct_rt_within", "median_cv_pct",
        }
        assert expected_keys == set(report.keys())

    def test_counts(self, sample_df):
        report = summarize_prm(sample_df)
        assert report["n_precursors"] == 3
        assert report["n_replicates"] == 3

    def test_pct_rt_within_all_small_deviations(self, sample_df):
        # All |dev| ≤ 0.5 which is within default threshold of 2.0
        report = summarize_prm(sample_df)
        assert report["pct_rt_within"] == pytest.approx(100.0)

    def test_pct_rt_within_strict_threshold(self, sample_df):
        # |RT dev| ≤ 0.1 after stable rounding (avoids 20.1−20.0 float noise).
        report = summarize_prm(sample_df, rt_dev_threshold=0.1)
        assert report["pct_rt_within"] == pytest.approx(5 / 9 * 100, rel=1e-3)

    def test_median_cv_is_finite(self, sample_df):
        report = summarize_prm(sample_df)
        assert math.isfinite(report["median_cv_pct"])

    def test_sub_dataframes_are_dataframes(self, sample_df):
        report = summarize_prm(sample_df)
        for key in ("cv", "missing", "dot_products", "rt_deviation"):
            assert isinstance(report[key], pd.DataFrame)


# ---------------------------------------------------------------------------
# Plate normalization
# ---------------------------------------------------------------------------

PLATE = "characteristics[plate]"


def _pools(plate_level, missing=(), runs_per_plate=3):
    """Pool rows: each peptide's ratio is its own level times its plate's."""
    peptides = {"PEPA": 0.5, "PEPB": 2.0, "PEPC": 8.0}
    rows = []
    for plate, level in plate_level.items():
        for r in range(runs_per_plate):
            run = f"{plate}_R{r}"
            for pep, base in peptides.items():
                if (plate, pep) not in missing:
                    rows.append({"Replicate": run, "File Name": f"{run}.raw", PLATE: plate,
                                 "Peptide Sequence": pep, "RatioLightToHeavy": base * level})
    return pd.DataFrame(rows)


def _spread_after(pools, factors):
    """Largest max/min ratio of a peptide across runs after correction (1 = plates agree)."""
    adjusted = adjust_ratio_by_plate(pools, factors, col_match=PLATE)
    per_peptide = adjusted.groupby("Peptide Sequence")["RatioLightToHeavy"]
    return (per_peptide.max() / per_peptide.min()).max()


class TestPlateConversionFactors:
    @pytest.mark.parametrize("log_transform", [True, False])
    def test_plate_reading_high_is_lowered(self, log_transform):
        pools = _pools({"P1": 1.0, "P2": 1.0, "P3": 1.5})
        table, factors, _ = get_plate_conversion_factors(pools, col_plate=PLATE, log_transform=log_transform)
        assert factors["P1"] == pytest.approx(1.0)
        assert factors["P3"] == pytest.approx(1 / 1.5)
        assert table.set_index(PLATE).loc["P3", "plate_median"] == pytest.approx(1.5)
        assert _spread_after(pools, factors) == pytest.approx(1.0)

    @pytest.mark.parametrize("log_transform", [True, False])
    def test_plates_end_up_level(self, log_transform):
        pools = _pools({"P1": 1.0, "P2": 2.0})
        _, factors, _ = get_plate_conversion_factors(pools, col_plate=PLATE, log_transform=log_transform)
        assert factors["P1"] > 1 > factors["P2"]
        assert _spread_after(pools, factors) == pytest.approx(1.0)

    def test_peptides_missing_from_a_plate_do_not_shift_it(self):
        # P2 lacks the high-ratio peptide; the peptides it has read as on P1
        pools = _pools({"P1": 1.0, "P2": 1.0}, missing={("P2", "PEPC")})
        _, factors, _ = get_plate_conversion_factors(pools, col_plate=PLATE, log_transform=True)
        assert factors["P2"] == pytest.approx(1.0)

    def test_non_positive_and_missing_ratios_are_left_out(self):
        pools = _pools({"P1": 1.0, "P2": 2.0})
        bad = pools.head(2).assign(RatioLightToHeavy=[0.0, float("nan")])
        _, factors, _ = get_plate_conversion_factors(pd.concat([pools, bad]), col_plate=PLATE,
                                                     log_transform=True)
        assert np.isfinite(list(factors.values())).all()
        assert _spread_after(pools, factors) == pytest.approx(1.0)


class TestAdjustRatioByPlate:
    def test_multiplies_by_the_factor(self):
        df = pd.DataFrame({PLATE: ["P1", "P2"], "RatioLightToHeavy": [1.0, 3.0]})
        out = adjust_ratio_by_plate(df, {"P1": 2.0, "P2": 0.5})
        assert out["RatioLightToHeavy"].tolist() == [2.0, 1.5]

    def test_rows_without_a_plate_are_left_out(self):
        df = pd.DataFrame({"Replicate": ["R1", "R2"], PLATE: ["P1", None],
                           "RatioLightToHeavy": [1.0, 3.0]})
        out = adjust_ratio_by_plate(df, {"P1": 2.0})
        assert out["Replicate"].tolist() == ["R1"]

    def test_plate_without_a_factor_raises(self):
        df = pd.DataFrame({PLATE: ["P1", "P9"], "RatioLightToHeavy": [1.0, 3.0]})
        with pytest.raises(KeyError, match="P9"):
            adjust_ratio_by_plate(df, {"P1": 2.0})


# ---------------------------------------------------------------------------
# Absolute concentrations
# ---------------------------------------------------------------------------

class TestGetAbsoluteConc:
    KEY_A = ("QR0000001", "PEPAK", "QR0000001_ABC|P11111")
    KEY_B = ("QR0000002", "PEPBR", "QR0000002_DEF|P22222")

    @pytest.fixture()
    def light_rows(self):
        return pd.DataFrame({
            "Protein Name": [self.KEY_A[2]] * 3 + [self.KEY_B[2]],
            "Peptide Sequence": ["PEPAK", "PEPAK", "PEPAK", "PEPBR"],
            "Precursor Charge": [2, 3, 2, 2],
            "File Name": ["f1.raw", "f1.raw", "f2.raw", "f1.raw"],
            "RatioLightToHeavy": [0.0123, 0.0125, 2.0, 1.0],
            "Heavy Amount [pmol]": [0.0123, 0.0123, 0.0123, 0.5],
        })

    def test_charge_states_are_averaged_and_not_rounded(self, light_rows):
        qreps = pd.DataFrame({"qRePS": ["QR0000001", "QR0000002"], "Amount per well [pmol]": [0.0123, 0.5]})
        wide = get_absolute_conc(qreps, light_rows)
        assert wide.loc[self.KEY_A, "f1.raw"] == pytest.approx(0.0124 * 0.0123, rel=1e-9)
        assert wide.loc[self.KEY_B, "f1.raw"] == pytest.approx(0.5)

    def test_amount_of_each_row(self, light_rows):
        wide = get_absolute_conc(None, light_rows, amount_col="Heavy Amount [pmol]")
        assert wide.loc[self.KEY_A, "f2.raw"] == pytest.approx(2.0 * 0.0123)
        assert wide.loc[self.KEY_B, "f1.raw"] == pytest.approx(0.5)

    def test_amount_of_each_row_needs_no_qreps_id(self, light_rows):
        rows = light_rows.assign(**{"Protein Name": light_rows["Protein Name"].str.replace("QR", "XX")})
        wide = get_absolute_conc(None, rows, amount_col="Heavy Amount [pmol]")
        assert wide.loc[("XX0000002_DEF|P22222", "PEPBR", "XX0000002_DEF|P22222"), "f1.raw"] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Report plots
# ---------------------------------------------------------------------------

class TestPlotAllAll:
    @pytest.fixture()
    def abs_df(self):
        index = pd.MultiIndex.from_tuples(
            [("QR0000001", "PEPAK", "QR0000001_ABC|P11111"),
             ("QR0000002", "PEPBR", "QR0000002_DEF|P22222")],
            names=["qRePS", "Peptide Sequence", "Protein Name"],
        )
        # the second protein has values in group B only
        return pd.DataFrame({"a1.raw": [1.0, np.nan], "a2.raw": [1.2, np.nan],
                             "b1.raw": [2.0, 3.0], "b2.raw": [2.2, 3.3]},
                            index=index).rename_axis(columns="File Name")

    @pytest.fixture()
    def sdrf(self):
        return pd.DataFrame({
            "source name": ["a1", "a2", "b1", "b2"],
            "comment[data file]": ["a1.raw", "a2.raw", "b1.raw", "b2.raw"],
            "characteristics[plate]": ["P1", "P2", "P1", "P2"],
            "factor value[group]": ["A", "A", "B", "B"],
        })

    def test_colour_defaults_to_the_group_column(self, abs_df, sdrf, tmp_path, capsys):
        pdf = plot_all_all(abs_df, sdrf, pdf_path=str(tmp_path / "all.pdf"))
        assert Path(pdf).stat().st_size > 0
        assert "Saved 2 plots" in capsys.readouterr().out

    def test_other_colour_column(self, abs_df, sdrf, tmp_path, capsys):
        plot_all_all(abs_df, sdrf, color_col="characteristics[plate]", pdf_path=str(tmp_path / "all.pdf"))
        assert "Saved 2 plots" in capsys.readouterr().out
