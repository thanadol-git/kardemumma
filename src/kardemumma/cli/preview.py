import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.show = lambda *a, **kw: None

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        prog="kardemumma-preview",
        description=(
            "Step 1: Import Skyline CSV + SDRF and generate an initial QC summary "
            "to inform cutoff choices for kardemumma-cutoff."
        ),
    )
    parser.add_argument("--skyline", required=True, metavar="CSV",
                        help="Skyline CSV export file")
    parser.add_argument("--sdrf", required=True, metavar="TSV",
                        help="SDRF .tsv file")
    parser.add_argument("--output", required=True, metavar="DIR",
                        help="Directory to write results into")
    parser.add_argument("--dotp", type=float, default=0.6, metavar="THRESHOLD",
                        help="Library dot product threshold used for the heavy/light scatter "
                             "(default: 0.6). The full distribution is always shown unfiltered.")
    args = parser.parse_args()

    import kardemumma as kdm

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("=== kardemumma-preview ===")

    # Import
    skyline_importer = kdm.ImportSkylineFile(args.skyline)
    skyline_data = skyline_importer.import_skyline_file()

    sdrf_obj = kdm.ImportSDRFFile(args.sdrf)
    sdrf_data = sdrf_obj.import_sdrf_file()
    kdm.readout_ms_type(args.sdrf)

    # Collect stats from raw data before any removal
    total_files = skyline_data["File Name"].nunique()
    total_replicates = skyline_data["Replicate"].nunique()
    total_proteins_raw = skyline_data["Protein Name"].nunique()
    total_peptides_raw = skyline_data["Peptide Sequence"].nunique()

    # SDRF stats
    sdrf_samples = sdrf_data["source name"].nunique() if "source name" in sdrf_data.columns else "n/a"
    sdrf_files = sdrf_data["comment[data file]"].nunique() if "comment[data file]" in sdrf_data.columns else "n/a"

    # Suggest and remove QC samples (matches notebook: remove before any analysis)
    qc_samples = skyline_importer.suggest_qc_samples(skyline_data)
    (out / "qc_samples.txt").write_text("\n".join(qc_samples))
    skyline_data = kdm.remove_qc_samples(skyline_data, qc_samples)
    study_files = skyline_data["File Name"].nunique()

    # Peptide counts at the chosen dotp threshold (matches notebook: filter then count)
    skyline_pivoted = kdm.filter_library_dot_product(skyline_data, threshold=args.dotp)
    peptide_counts = kdm.summarise_peptide_counts(skyline_pivoted)
    report_summary, _ = kdm.report_peptide_protein_summary(peptide_counts)

    # Save tables
    skyline_data.to_csv(out / "skyline_raw.csv", index=False)
    sdrf_data.to_csv(out / "sdrf.csv", index=False)
    peptide_counts.to_csv(out / "peptide_counts_raw.csv", index=False)

    # Write descriptive text report
    import datetime
    lines = [
        "=== kardemumma preview report ===",
        f"Date : {datetime.date.today()}",
        "",
        "Input files",
        "-----------",
        f"Skyline export : {Path(args.skyline).name}",
        f"SDRF file      : {Path(args.sdrf).name}",
        "",
        "Skyline — sample overview",
        "--------------------------",
        f"Total raw files (File Name)  : {total_files}",
        f"Total replicates             : {total_replicates}",
        f"QC samples detected          : {len(qc_samples)}",
        f"  " + (", ".join(qc_samples) if qc_samples else "(none)"),
        f"Study samples (after QC rm)  : {study_files}",
        "",
        "Skyline — content (before QC removal, no dotp filter)",
        "------------------------------------------------------",
        f"Unique proteins  : {total_proteins_raw}",
        f"Unique peptides  : {total_peptides_raw}",
        "",
        f"After QC removal + dotp > {args.dotp}",
        "-" * 38,
        f"Unique proteins  : {report_summary['num_unique_proteins']}",
        f"Unique peptides  : {report_summary['num_unique_peptides']}",
        "",
        "SDRF — sample overview",
        "-----------------------",
        f"Unique source names  : {sdrf_samples}",
        f"Unique data files    : {sdrf_files}",
        "",
        "Notes",
        "-----",
        f"dotp threshold used for scatter plot : {args.dotp}",
        "Run kardemumma-cutoff to apply final cutoffs and batch correction.",
    ]
    (out / "report_summary.txt").write_text("\n".join(lines))

    # Dot product distribution — full unfiltered view with threshold marker
    kdm.plot_library_dot_product_distribution(skyline_data)
    plt.axvline(args.dotp, color="red", linestyle="--", linewidth=1.5,
                label=f"cutoff = {args.dotp}")
    plt.legend()
    plt.savefig(out / "dot_product_distribution.pdf", bbox_inches="tight")
    plt.close("all")

    # Heavy/light scatter — at the chosen dotp threshold
    kdm.plot_heavy_light_scatter(peptide_counts)
    plt.savefig(out / "heavy_light_scatter.pdf", bbox_inches="tight")
    plt.close("all")

    print(f"\nPreview complete → {out}")
    print("  Report : report_summary.txt")
    print("  Tables : skyline_raw.csv, sdrf.csv, peptide_counts_raw.csv, qc_samples.txt")
    print("  Plots  : dot_product_distribution.pdf, heavy_light_scatter.pdf")
    print()
    print("Review the plots, then run kardemumma-cutoff with your chosen thresholds.")
