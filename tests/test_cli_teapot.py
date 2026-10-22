"""
kardemumma-preview, -cutoff, -report and -ipynb on a synthetic teapot table with a known plate effect.

Run with:  pytest tests/test_cli_teapot.py -v
"""

import sys

import matplotlib
matplotlib.use("Agg")
import nbformat
import numpy as np
import pandas as pd
import pytest

import kardemumma as kdm
from kardemumma.cli import cutoff, ipynb, preview, report

# Every ratio on a plate reads this much off
PLATE_LEVEL = {"Plate_1": 1.0, "Plate_2": 1.25, "Plate_3": 0.8}
# peptide: (protein, light/heavy ratio, known heavy amount in pmol)
PEPTIDES = {
    "PEPAAK": ("QR0000001_ABC|P11111", 0.5, 0.2),
    "PEPBBK": ("QR0000001_ABC|P11111", 0.6, 0.2),
    "PEPCCR": ("QR0000002_DEF|P22222", 2.0, 1.0),
    "PEPDDR": ("QR0000002_DEF|P22222", 1.8, 1.0),
    "PEPEEK": ("QR0000003_GHI|P33333", 8.0, 0.05),
    "PEPFFK": ("QR0000003_GHI|P33333", 7.0, 0.05),
}


def _runs():
    """(sample_id, plate, characteristics[Sample], group): 3 pools and 4 patients per plate."""
    runs = []
    for p, plate in enumerate(PLATE_LEVEL, 1):
        runs += [(f"Study_{plate}_pool{i}", plate, "Pool", "Pool") for i in range(3)]
        runs += [(f"Study_{plate}_s{i}", plate, "Patient", "A" if i < 2 else "B") for i in range(4)]
    return runs


@pytest.fixture()
def inputs(tmp_path):
    """teapot's confident pairs (and a run the SDRF does not list), the runs' SDRF."""
    rng = np.random.default_rng(0)
    rows = []
    for sid, plate, _, group in _runs() + [("Study_QC_injection_1", "Plate_1", None, None)]:
        for pep, (protein, ratio, amount) in PEPTIDES.items():
            ratio *= PLATE_LEVEL[plate] * (2.0 if group == "B" and "ABC" in protein else 1.0)
            ratio *= float(np.exp(rng.normal(0, 0.02)))
            rows.append({
                "sample_id": sid, "protein": protein, "peptide": pep, "peptidoform": pep, "charge": 2,
                "light_abundance": ratio * 1e6, "heavy_abundance": 1e6,
                "ratio_L_H": ratio, "ratio_H_L": 1 / ratio, "heavy_pmol_known": amount,
                "abs_pmol_light": ratio * amount, "light_q_value": 1e-4, "heavy_q_value": 1e-4,
                "light_rt_seconds": 600.0, "heavy_rt_seconds": 600.5,
            })
    route = tmp_path / "results" / "quant" / "enc_dia"
    route.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(route / "peptide.abs_quant.confident.tsv", sep="\t", index=False)

    sdrf = pd.DataFrame([{
        "source name": sid,
        "characteristics[plate]": plate,
        "characteristics[Sample]": sample,
        "assay name": f"run {i}",
        "comment[data file]": f"{sid}.raw",
        "comment[proteomics data acquisition method]": "NT=Parallel Reaction Monitoring;AC=MS:1002956",
        "comment[ProteomEdge]": 24001,
        "factor value[group]": group,
    } for i, (sid, plate, sample, group) in enumerate(_runs(), 1)])
    sdrf_path = tmp_path / "experiment.sdrf.tsv"
    sdrf.to_csv(sdrf_path, sep="\t", index=False)
    return tmp_path / "results", sdrf_path


def _run(module, monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", [module.__name__, *map(str, args)])
    module.main()


def test_preview_cutoff_report(inputs, tmp_path, monkeypatch):
    results, sdrf = inputs
    out = tmp_path / "out"
    monkeypatch.setattr(kdm, "plot_umap", lambda *a, **kw: None)  # UMAP compiles slowly
    monkeypatch.setattr(kdm, "fetch_fasta", lambda lot: pd.DataFrame({"id": [], "sequence": []}))

    _run(preview, monkeypatch, "--teapot", results, "--sdrf", sdrf, "--output", out)
    assert (out / "teapot_raw.csv").exists() and not (out / "dot_product_distribution.pdf").exists()

    _run(cutoff, monkeypatch, "--teapot", results, "--sdrf", sdrf, "--output", out,
         "--light-cutoff", 10, "--heavy-cutoff", 10, "--cv-percentile", 50, "--pool-value", "Pool")
    factors = pd.read_csv(out / "conversion_factors.csv").set_index("plate")["conversion_factor"]
    assert factors["Plate_2"] / factors["Plate_1"] == pytest.approx(1 / 1.25, rel=0.03)
    assert factors["Plate_3"] / factors["Plate_1"] == pytest.approx(1 / 0.8, rel=0.03)
    assert not list(out.glob("skyline_*"))
    adj = pd.read_csv(out / "teapot_merge_adj.csv")
    assert "Study_QC_injection_1" not in set(adj["Replicate"])
    pools = adj[(adj["characteristics[Sample]"] == "Pool") & adj["RatioLightToHeavy"].notna()]
    plate_medians = (
        pools.assign(log_ratio=np.log(pools["RatioLightToHeavy"]))
        .groupby(["Peptide Sequence", "characteristics[plate]"])["log_ratio"].median().unstack()
    )
    assert (plate_medians.max(axis=1) - plate_medians.min(axis=1)).max() < 0.1

    _run(report, monkeypatch, "--output", out, "--sdrf", sdrf, "--group-col", "characteristics[plate]")
    abs_df = pd.read_csv(out / "abs_df.csv")
    assert set(abs_df["qRePS"]) == {"QR0000001", "QR0000002", "QR0000003"}
    assert not (out / "qreps_table.csv").exists()
    for pdf in ("protein_concentrations.pdf", "protein_median_concentrations.pdf", "pca.pdf"):
        assert (out / pdf).stat().st_size > 0


def test_cutoff_stops_when_no_peptide_passes(inputs, tmp_path, monkeypatch, capsys):
    results, sdrf = inputs
    with pytest.raises(SystemExit):
        _run(cutoff, monkeypatch, "--teapot", results, "--sdrf", sdrf, "--output", tmp_path / "out",
             "--light-cutoff", 600, "--heavy-cutoff", 600)
    assert "heavy_light_scatter.pdf" in capsys.readouterr().err


def test_ipynb_for_teapot(inputs, tmp_path, monkeypatch):
    results, sdrf = inputs
    _run(ipynb, monkeypatch, "--teapot", results, "--sdrf", sdrf, "--output", tmp_path,
         "--light-cutoff", 10, "--heavy-cutoff", 10)
    nb = nbformat.read(tmp_path / "ratio_analysis.ipynb", as_version=4)
    sources = [c.source for c in nb.cells if c.cell_type == "code"]
    for src in sources:
        compile(src, "<cell>", "exec")
    text = "\n".join(sources)
    assert "kdm.ImportTeapotFile(" in text and "pivot_confident_pairs" in text
    assert "skyline" not in text and "qc_samples" not in text
    assert "%pip" not in text and "DownStream" not in text


def test_ipynb_needs_both_groups(inputs, tmp_path, monkeypatch):
    results, sdrf = inputs
    with pytest.raises(SystemExit):
        _run(ipynb, monkeypatch, "--teapot", results, "--sdrf", sdrf, "--output", tmp_path,
             "--group-a", "A")
