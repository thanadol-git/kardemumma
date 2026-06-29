import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.show = lambda *a, **kw: None

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        prog="kardemumma-cutoff",
        description=(
            "Step 2: Apply QC filters, compute pool CV, select normalization peptides, "
            "and batch-correct ratios. Reads the Skyline CSV + SDRF directly."
        ),
    )
    parser.add_argument("--skyline", required=True, metavar="CSV",
                        help="Skyline CSV export file")
    parser.add_argument("--sdrf", required=True, metavar="TSV",
                        help="SDRF .tsv file")
    parser.add_argument("--output", required=True, metavar="DIR",
                        help="Directory to write results into")
    parser.add_argument("--dotp", type=float, default=0.6, metavar="THRESHOLD",
                        help="Library dot product threshold (default: 0.6)")
    parser.add_argument("--light-cutoff", type=int, default=600, metavar="N",
                        help="Min light sample count per peptide (default: 600)")
    parser.add_argument("--heavy-cutoff", type=int, default=600, metavar="N",
                        help="Min heavy sample count per peptide (default: 600)")
    parser.add_argument("--pool-value", default="PlasmaPool", metavar="STR",
                        help="Value in characteristics[Sample] that marks pool samples "
                             "(default: PlasmaPool)")
    parser.add_argument("--bad-pools", nargs="+", default=[], metavar="FILE",
                        help="Raw file names to exclude from pool normalization")
    parser.add_argument("--bad-peptides", nargs="+", default=[], metavar="SEQ",
                        help="Peptide sequences to exclude from normalization "
                             "(matched against Peptide Sequence column)")
    parser.add_argument("--exclude-plates", nargs="+", default=[], metavar="PLATE",
                        help="Plate names to drop from the final adjusted dataset")
    parser.add_argument("--cv-percentile", type=float, default=10.0, metavar="PCT",
                        help="Inter-plate CV percentile threshold for normalization "
                             "peptide selection (default: 10)")
    parser.add_argument("--no-log-transform", dest="log_transform", action="store_false",
                        help="Disable log-transform in plate correction (default: enabled)")
    parser.set_defaults(log_transform=True)
    args = parser.parse_args()

    import kardemumma as kdm

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("=== kardemumma-cutoff ===")

    # --- Import ---
    skyline_importer = kdm.ImportSkylineFile(args.skyline)
    skyline_data = skyline_importer.import_skyline_file()

    sdrf_obj = kdm.ImportSDRFFile(args.sdrf)
    sdrf_data = sdrf_obj.import_sdrf_file()

    # --- Remove QC samples (by File Name) ---
    qc_samples = skyline_importer.suggest_qc_samples(skyline_data)
    skyline_data = kdm.remove_qc_samples(skyline_data, qc_samples)

    # --- Library dot product filter (also pivots to heavy/light wide format) ---
    skyline_clean = kdm.filter_library_dot_product(skyline_data, threshold=args.dotp)

    # --- Peptide detection summary (on pivoted data) ---
    peptide_counts = kdm.summarise_peptide_counts(skyline_clean)
    kdm.report_peptide_protein_summary(peptide_counts)

    # --- Count cutoff → selected_peptides are Peptide column values for merge filtering ---
    filtered_peptide_counts = kdm.filter_peptide_counts(
        peptide_counts, light_cutoff=args.light_cutoff, heavy_cutoff=args.heavy_cutoff
    )
    selected_report, selected_peptides = kdm.report_peptide_protein_summary(filtered_peptide_counts)

    # --- Merge with SDRF (uses raw long-format skyline_data, filtered to selected_peptides) ---
    merge_obj = kdm.MergeFiles(skyline_data, sdrf_data, selected_peptides)
    skyline_merge = merge_obj.merge_files()

    # --- Pool data ---
    skyline_pool = (
        skyline_merge[skyline_merge["characteristics[Sample]"] == args.pool_value]
        .sort_values(["Replicate", "Peptide", "Isotope Label Type"])
        .reset_index(drop=True)
    )
    print(f"Pool samples: {skyline_pool['Replicate'].nunique()} unique replicates")

    # --- Pool QC plot ---
    kdm.plot_pool_boxplot(skyline_pool)
    plt.savefig(out / "pool_boxplot.pdf", bbox_inches="tight")
    plt.close("all")

    # --- Intra-plate CV ---
    peptide_plate_stats = kdm.calculate_intra_plate_cv(
        skyline_pool, col_name="characteristics[plate]"
    )
    kdm.plot_intra_plate_cv_stats(peptide_plate_stats, col_name="characteristics[plate]")
    plt.savefig(out / "intra_plate_cv.pdf", bbox_inches="tight")
    plt.close("all")

    # --- Inter-plate CV ---
    interplate_cv = kdm.calculate_inter_plate_cv(peptide_plate_stats)
    kdm.plot_inter_plate_cv_kde(interplate_cv)
    plt.savefig(out / "inter_plate_cv_kde.pdf", bbox_inches="tight")
    plt.close("all")

    kdm.plot_cumulative_peptide_count_by_cv(interplate_cv)
    plt.savefig(out / "cumulative_cv.pdf", bbox_inches="tight")
    plt.close("all")

    # --- Normalization peptides (lowest inter-plate CV) ---
    selected_norm_peptides = kdm.get_lowest_cv_peptides(interplate_cv, args.cv_percentile)
    print(f"Normalization peptides selected: {len(selected_norm_peptides)}")

    pool_selected = skyline_pool[
        skyline_pool["Peptide Sequence"].isin(selected_norm_peptides)
    ].copy()
    if args.bad_pools:
        pool_selected = pool_selected[~pool_selected["File Name"].isin(args.bad_pools)]
        print(f"Excluded {len(args.bad_pools)} bad pool(s)")
    if args.bad_peptides:
        pool_selected = pool_selected[~pool_selected["Peptide Sequence"].isin(args.bad_peptides)]
        print(f"Excluded {len(args.bad_peptides)} bad peptide(s)")

    kdm.plot_pool_boxplot(pool_selected)
    plt.savefig(out / "norm_peptides_boxplot.pdf", bbox_inches="tight")
    plt.close("all")

    # --- Plate conversion factors ---
    plate_factor_table, conversion_factors, _model = kdm.get_plate_conversion_factors(
        pool_selected,
        col_plate="characteristics[plate]",
        log_transform=args.log_transform,
    )
    kdm.plot_plate_conversion_factors(
        pool_selected,
        col_plate="characteristics[plate]",
        log_transform=args.log_transform,
    )
    plt.savefig(out / "plate_conversion_factors.pdf", bbox_inches="tight")
    plt.close("all")

    # --- Build batch-corrected dataset ---
    # Drop QC samples by Replicate name (belt-and-suspenders alongside the File Name removal above)
    skyline_merge_adj = skyline_merge[
        ~skyline_merge["Replicate"].isin(qc_samples)
    ].copy()
    for plate in args.exclude_plates:
        skyline_merge_adj = skyline_merge_adj[
            skyline_merge_adj["characteristics[plate]"] != plate
        ]
    skyline_merge_adj = kdm.adjust_ratio_by_plate(skyline_merge_adj, conversion_factors)

    # --- Save tables ---
    peptide_counts.to_csv(out / "peptide_counts.csv", index=False)
    filtered_peptide_counts.to_csv(out / "filtered_peptide_counts.csv", index=False)
    pd.DataFrame([selected_report]).to_csv(out / "selected_peptides_report.csv", index=False)
    pd.Series(list(selected_norm_peptides), name="Peptide Sequence").to_csv(
        out / "norm_peptides.csv", index=False
    )
    skyline_merge.to_csv(out / "skyline_merge.csv", index=False)
    skyline_pool.to_csv(out / "pool_data.csv", index=False)
    peptide_plate_stats.to_csv(out / "peptide_plate_stats.csv", index=False)
    interplate_cv.to_csv(out / "interplate_cv.csv", index=False)
    plate_factor_table.to_csv(out / "plate_factor_table.csv", index=False)
    pd.DataFrame(
        list(conversion_factors.items()), columns=["plate", "conversion_factor"]
    ).to_csv(out / "conversion_factors.csv", index=False)
    skyline_merge_adj.to_csv(out / "skyline_merge_adj.csv", index=False)

    # Save parameters for reproducibility
    params = {
        "skyline": str(Path(args.skyline).resolve()),
        "sdrf": str(Path(args.sdrf).resolve()),
        "dotp_threshold": args.dotp,
        "light_cutoff": args.light_cutoff,
        "heavy_cutoff": args.heavy_cutoff,
        "pool_value": args.pool_value,
        "bad_pools": args.bad_pools,
        "bad_peptides": args.bad_peptides,
        "exclude_plates": args.exclude_plates,
        "cv_percentile": args.cv_percentile,
        "log_transform": args.log_transform,
    }
    (out / "cutoff_params.json").write_text(json.dumps(params, indent=2))

    print(f"\nCutoff results → {out}")
    print("  Tables : peptide_counts.csv, filtered_peptide_counts.csv,")
    print("           skyline_merge.csv, pool_data.csv, peptide_plate_stats.csv,")
    print("           interplate_cv.csv, conversion_factors.csv, skyline_merge_adj.csv")
    print("  Plots  : pool_boxplot.pdf, intra_plate_cv.pdf, inter_plate_cv_kde.pdf,")
    print("           cumulative_cv.pdf, norm_peptides_boxplot.pdf, plate_conversion_factors.pdf")
    print("  Config : cutoff_params.json")
    print()
    print(f"Next: kardemumma-report --output {args.output} --sdrf {args.sdrf}")
