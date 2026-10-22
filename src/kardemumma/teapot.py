"""
Import of teapot absolute-quantification results.
"""

import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

__all__ = ["ImportTeapotFile", "find_teapot_confident_table", "pivot_confident_pairs"]

_CONFIDENT_FILE = "peptide.abs_quant.confident.tsv"

_TEAPOT_COLS = [
    "sample_id", "protein", "peptide", "peptidoform", "charge",
    "light_abundance", "heavy_abundance", "ratio_L_H", "ratio_H_L",
    "heavy_pmol_known", "abs_pmol_light",
    "light_q_value", "heavy_q_value", "light_rt_seconds", "heavy_rt_seconds",
]


def find_teapot_confident_table(path: str, route: Optional[str] = None) -> Path:
    p = Path(path)
    if p.is_file():
        return p
    tables = {t.parent.name: t for t in sorted(p.glob(f"quant/*/{_CONFIDENT_FILE}"))}
    if route is not None:
        if route not in tables:
            raise FileNotFoundError(
                f"No {_CONFIDENT_FILE} for route '{route}' under {p / 'quant'}; "
                f"routes: {sorted(tables)}"
            )
        return tables[route]
    if not tables:
        raise FileNotFoundError(f"No {_CONFIDENT_FILE} under {p / 'quant'}")
    if len(tables) > 1:
        raise FileNotFoundError(f"Routes {sorted(tables)} under {p / 'quant'}: choose one with route.")
    return next(iter(tables.values()))


def pivot_confident_pairs(df: pd.DataFrame) -> pd.DataFrame:
    print(df["Isotope Label Type"].value_counts())
    wide = df.pivot_table(
        index=["Replicate", "Protein Name", "Peptide"],
        columns="Isotope Label Type",
        values="Normalized Area",
        aggfunc="first",
    )
    return wide.reset_index()


class ImportTeapotFile:
    def __init__(self, file_path: str, route: Optional[str] = None):
        """
        Args:
            file_path: teapot's ``--outdir``, or ``peptide.abs_quant.confident.tsv`` itself.
            route: Route under ``<outdir>/quant/``, when the outdir holds more than one.
        """
        self.file_path = find_teapot_confident_table(file_path, route)
        self.route = self.file_path.parent.name
        print(f"Importing teapot file: {self.file_path}")

    def import_teapot_file(self, sdrf_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        df = pd.read_csv(self.file_path, sep="\t", dtype={"sample_id": str})
        missing = [c for c in _TEAPOT_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"{self.file_path} is not teapot's {_CONFIDENT_FILE}: no {missing}")
        if df.empty:
            raise ValueError(f"No confident pair in {self.file_path}")

        df["File Name"] = df["sample_id"]
        if sdrf_df is not None:
            data_files = sdrf_df["comment[data file]"].dropna().astype(str)
            file_of = dict(zip(data_files.map(lambda f: os.path.splitext(f)[0]), data_files))
            listed = df["sample_id"].isin(file_of)
            if not listed.any():
                raise ValueError(
                    "No teapot sample_id is an SDRF data file (comment[data file] without its "
                    f"extension). sample_ids: {sorted(df['sample_id'].unique())[:5]} ..."
                )
            if not listed.all():
                unlisted = sorted(df.loc[~listed, "sample_id"].unique())
                print(f"Runs the SDRF does not list, left out ({len(unlisted)}): {unlisted}")
            df = df[listed].copy()
            df["File Name"] = df["sample_id"].map(file_of)

        channels = []
        for label in ("light", "heavy"):
            heavy = label == "heavy"
            channels.append(pd.DataFrame({
                "Precursor": df["peptidoform"] + df["charge"].map(lambda z: "+" * z)
                             + (" (heavy)" if heavy else ""),
                "Replicate": df["sample_id"],
                "File Name": df["File Name"],
                "Peptide Retention Time": df[f"{label}_rt_seconds"] / 60,
                "Precursor Charge": df["charge"],
                "Peptide Sequence": df["peptide"],
                "Peptide": df["peptidoform"],
                "Normalized Area": df[f"{label}_abundance"],
                "RatioLightToHeavy": np.nan if heavy else df["ratio_L_H"],
                "Protein Name": df["protein"],
                "Isotope Label Type": label,
                "Q Value": df[f"{label}_q_value"],
                "Heavy Amount [pmol]": df["heavy_pmol_known"],
            }))
        long_df = (
            pd.concat(channels, ignore_index=True)
            .sort_values(["Protein Name", "Peptide", "Precursor Charge", "Isotope Label Type", "Replicate"])
            .reset_index(drop=True)
        )

        print(f"Route: {self.route}")
        print(f"Runs: {df['sample_id'].nunique()}  Proteins: {df['protein'].nunique()}  "
              f"Peptidoforms: {df['peptidoform'].nunique()}")
        print(f"Confident light/heavy pairs: {len(df)}")
        return long_df
