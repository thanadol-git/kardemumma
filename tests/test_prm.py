"""
Tests for skyline_qc.prm

Run with:  pytest tests/test_prm.py -v
"""

import math

import pandas as pd
import pytest

from kardemumma.prm import (
    compute_cv,
    dot_product_summary,
    filter_library_dot_product,
    flag_missing_values,
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
