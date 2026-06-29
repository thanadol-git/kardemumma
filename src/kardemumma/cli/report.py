import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.show = lambda *a, **kw: None

import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        prog="kardemumma-report",
        description=(
            "Step 3: Fetch qRePS data, calculate absolute protein concentrations, "
            "and generate all report plots. Reads skyline_merge_adj.csv written by "
            "kardemumma-cutoff."
        ),
    )
    parser.add_argument("--output", required=True, metavar="DIR",
                        help="Output directory (must contain skyline_merge_adj.csv "
                             "from kardemumma-cutoff)")
    parser.add_argument("--sdrf", required=True, metavar="TSV",
                        help="SDRF .tsv file (provides sample metadata and qRePS lot number)")
    parser.add_argument("--color-col", default=None, metavar="COL",
                        help="SDRF column to color concentration plots by "
                             "(default: last 'factor value[...]' column)")
    args = parser.parse_args()

    import kardemumma as kdm

    out = Path(args.output)
    adj_path = out / "skyline_merge_adj.csv"
    if not adj_path.exists():
        print(f"Error: {adj_path} not found. Run kardemumma-cutoff first.", file=sys.stderr)
        sys.exit(1)

    print("=== kardemumma-report ===")

    # --- Load inputs ---
    skyline_merge_adj = pd.read_csv(adj_path)
    sdrf_obj = kdm.ImportSDRFFile(args.sdrf)
    sdrf_data = sdrf_obj.import_sdrf_file()

    # --- Fetch qRePS data from ProteomEdge (no authentication required) ---
    qreps_lot = sdrf_obj.extract_qreps_lot_number()
    print(f"Fetching qRePS table for lot {qreps_lot} ...")
    qreps_table = kdm.fetch_qreps_table(qreps_lot)
    print(f"Fetching FASTA for lot {qreps_lot} ...")
    fasta_df = kdm.fetch_fasta(qreps_lot)

    # --- Absolute quantification ---
    print("Calculating absolute concentrations ...")
    abs_df = kdm.get_absolute_conc(qreps_table, skyline_merge_adj)

    # --- Save tables ---
    qreps_table.to_csv(out / "qreps_table.csv", index=False)
    fasta_df.to_csv(out / "fasta.csv", index=False)
    abs_df.reset_index().to_csv(out / "abs_df.csv", index=False)
    print("  Saved: qreps_table.csv, fasta.csv, abs_df.csv")

    # --- Per-protein concentration plots (boxplot + strip + median line) ---
    print("Generating protein concentration PDF ...")
    kdm.plot_all_all(
        abs_df,
        sdrf_data,
        color_col=args.color_col,
        pdf_path=str(out / "protein_concentrations.pdf"),
    )

    # --- Per-protein median concentration plots ---
    print("Generating median concentration PDF ...")
    kdm.plot_all_median_peptide_concentration_by_group(
        abs_df,
        sdrf_data,
        group_col=args.color_col,
        pdf_path=str(out / "protein_median_concentrations.pdf"),
    )

    # --- PCA ---
    print("Generating PCA ...")
    try:
        kdm.plot_pca(abs_df, sdrf_data, color_col=args.color_col)
        plt.savefig(out / "pca.pdf", bbox_inches="tight")
        plt.close("all")
    except Exception as exc:
        print(f"  PCA skipped: {exc}")
        plt.close("all")

    # --- UMAP ---
    print("Generating UMAP ...")
    try:
        kdm.plot_umap(abs_df, sdrf_data, color_col=args.color_col)
        plt.savefig(out / "umap.pdf", bbox_inches="tight")
        plt.close("all")
    except Exception as exc:
        print(f"  UMAP skipped: {exc}")
        plt.close("all")

    print(f"\nReport complete → {out}")
    print("  Tables : qreps_table.csv, fasta.csv, abs_df.csv")
    print("  Plots  : protein_concentrations.pdf, protein_median_concentrations.pdf,")
    print("           pca.pdf, umap.pdf")
