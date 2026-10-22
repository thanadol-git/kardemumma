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
            "Step 1: Import a Skyline CSV or teapot results + SDRF and generate an "
            "initial QC summary to inform cutoff choices for kardemumma-cutoff."
        ),
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--skyline", metavar="CSV",
                        help="Skyline CSV export file")
    source.add_argument("--teapot", metavar="PATH",
                        help="teapot's --outdir, or its peptide.abs_quant.confident.tsv. "
                             "A run's sample_id is its SDRF data file without the extension")
    parser.add_argument("--route", metavar="ROUTE",
                        help="teapot route (e.g. enc_dia), when the outdir holds more than one")
    parser.add_argument("--sdrf", required=True, metavar="TSV",
                        help="SDRF .tsv file")
    parser.add_argument("--output", required=True, metavar="DIR",
                        help="Directory to write results into")
    parser.add_argument("--dotp", type=float, default=0.6, metavar="THRESHOLD",
                        help="Skyline: library dot product threshold used for the heavy/light "
                             "scatter (default: 0.6). The full distribution is always shown "
                             "unfiltered. teapot: every confident pair counts, in both channels.")
    args = parser.parse_args()

    import kardemumma as kdm

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("=== kardemumma-preview ===")

    # Import
    sdrf_obj = kdm.ImportSDRFFile(args.sdrf)
    sdrf_data = sdrf_obj.import_sdrf_file()
    kdm.readout_ms_type(args.sdrf)

    if args.teapot:
        skyline_importer = kdm.ImportTeapotFile(args.teapot, route=args.route)
        skyline_data = skyline_importer.import_teapot_file(sdrf_df=sdrf_data)
    else:
        skyline_importer = kdm.ImportSkylineFile(args.skyline)
        skyline_data = skyline_importer.import_skyline_file()

    # Collect stats from raw data before any removal
    total_files = skyline_data["File Name"].nunique()
    total_replicates = skyline_data["Replicate"].nunique()
    total_proteins_raw = skyline_data["Protein Name"].nunique()
    total_peptides_raw = skyline_data["Peptide Sequence"].nunique()

    # SDRF stats
    sdrf_samples = sdrf_data["source name"].nunique() if "source name" in sdrf_data.columns else "n/a"
    sdrf_files = sdrf_data["comment[data file]"].nunique() if "comment[data file]" in sdrf_data.columns else "n/a"

    # Suggest and remove QC samples (matches notebook: remove before any analysis);
    # teapot's runs are the SDRF's
    qc_samples = [] if args.teapot else skyline_importer.suggest_qc_samples(skyline_data)
    (out / "qc_samples.txt").write_text("\n".join(qc_samples))
    skyline_data = kdm.remove_qc_samples(skyline_data, qc_samples)
    study_files = skyline_data["File Name"].nunique()

    # Peptide counts at the chosen dotp threshold (matches notebook: filter then count);
    # teapot reports no dot product: its confident pairs count, in both channels
    if args.teapot:
        skyline_pivoted = kdm.pivot_confident_pairs(skyline_data)
        filter_desc = "teapot's confident pairs"
    else:
        skyline_pivoted = kdm.filter_library_dot_product(skyline_data, threshold=args.dotp)
        filter_desc = f"dotp > {args.dotp}"
    peptide_counts = kdm.summarise_peptide_counts(skyline_pivoted)
    report_summary, _ = kdm.report_peptide_protein_summary(peptide_counts)

    # Save tables (named after the source: skyline_raw.csv or teapot_raw.csv)
    prefix = "teapot" if args.teapot else "skyline"
    skyline_data.to_csv(out / f"{prefix}_raw.csv", index=False)
    sdrf_data.to_csv(out / "sdrf.csv", index=False)
    peptide_counts.to_csv(out / "peptide_counts_raw.csv", index=False)

    # Write descriptive text report
    import datetime
    source = "teapot" if args.teapot else "Skyline"
    source_file = f"teapot results : {skyline_importer.file_path}" if args.teapot \
        else f"Skyline export : {Path(args.skyline).name}"
    lines = [
        "=== kardemumma preview report ===",
        f"Date : {datetime.date.today()}",
        "",
        "Input files",
        "-----------",
        source_file,
        f"SDRF file      : {Path(args.sdrf).name}",
        "",
        f"{source} — sample overview",
        "--------------------------",
        f"Total raw files (File Name)  : {total_files}",
        f"Total replicates             : {total_replicates}",
        f"QC samples detected          : {len(qc_samples)}",
        f"  " + (", ".join(qc_samples) if qc_samples else "(none)"),
        f"Study samples (after QC rm)  : {study_files}",
        "",
        f"{source} — content (before QC removal, no dotp filter)",
        "------------------------------------------------------",
        f"Unique proteins  : {total_proteins_raw}",
        f"Unique peptides  : {total_peptides_raw}",
        "",
        f"After QC removal + {filter_desc}",
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
        f"Filter used for scatter plot : {filter_desc}",
        "Run kardemumma-cutoff to apply final cutoffs and batch correction.",
    ]
    (out / "report_summary.txt").write_text("\n".join(lines))

    # Dot product distribution — full unfiltered view with threshold marker
    if not args.teapot:
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
    print(f"  Tables : {prefix}_raw.csv, sdrf.csv, peptide_counts_raw.csv, qc_samples.txt")
    print("  Plots  : " + ("" if args.teapot else "dot_product_distribution.pdf, ")
          + "heavy_light_scatter.pdf")
    print()
    print("Review the plots, then run kardemumma-cutoff with your chosen thresholds.")
