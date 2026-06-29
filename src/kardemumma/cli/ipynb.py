import argparse
import sys
from pathlib import Path


def _make_notebook(skyline_path, sdrf_path, dotp, light_cutoff, heavy_cutoff,
                   pool_value, cv_percentile):
    import nbformat

    nb = nbformat.v4.new_notebook()
    code = nbformat.v4.new_code_cell
    md = nbformat.v4.new_markdown_cell

    cells = [
        # --- Environment check ---
        code(
            "import sys, kardemumma\n"
            "print(sys.executable)\n"
            "print(kardemumma.__file__)"
        ),

        # --- Imports ---
        code(
            "from pathlib import Path\n"
            "import os, re, sys, importlib\n"
            "import matplotlib.pyplot as plt\n"
            "import numpy as np\n"
            "import pandas as pd\n"
            "import seaborn as sns\n"
            "import kardemumma as kdm"
        ),

        # --- Skyline import ---
        code(
            f'skyline_path = "{skyline_path}"\n\n'
            "skyline_importer = kdm.ImportSkylineFile(skyline_path)\n"
            "skyline_data = skyline_importer.import_skyline_file()"
        ),

        # --- SDRF readout ---
        code(
            f'sdrf_path = "{sdrf_path}"\n\n'
            "kdm.readout_ms_type(sdrf_path)"
        ),

        # --- iRT peptides ---
        code("kdm.get_irt_peptides(skyline_data)"),

        # --- SDRF import + file comparison ---
        code(
            "sdrf_data = kdm.ImportSDRFFile(sdrf_path)\n"
            "sdrf_data_file = sdrf_data.import_sdrf_file()\n\n"
            "not_in_sdrf = set(skyline_data['File Name']) - set(sdrf_data_file['comment[data file]'])\n"
            "not_in_skyline = set(sdrf_data_file['comment[data file]']) - set(skyline_data['File Name'])\n"
            "print('In Skyline but not SDRF:', not_in_sdrf)\n"
            "print('In SDRF but not Skyline:', not_in_skyline)"
        ),

        # --- QC samples ---
        code(
            "qc_samples = skyline_importer.suggest_qc_samples(skyline_data)\n"
            "print(qc_samples)"
        ),

        # --- Remove QC ---
        code("skyline_data = kdm.remove_qc_samples(skyline_data, qc_samples)"),

        # --- Dot product distribution ---
        code("kdm.plot_library_dot_product_distribution(skyline_data)"),

        # --- Dot product filter ---
        code(f"skyline_clean = kdm.filter_library_dot_product(skyline_data, threshold={dotp})"),

        # --- Peptide counts ---
        code("peptide_counts = kdm.summarise_peptide_counts(skyline_clean)"),

        # --- Report summary ---
        code("report_summary, peptide_list = kdm.report_peptide_protein_summary(peptide_counts)"),

        # --- Heavy/light scatter ---
        code("kdm.plot_heavy_light_scatter(peptide_counts)"),

        # --- Count cutoff ---
        code(
            f"filtered_peptide_counts = kdm.filter_peptide_counts(\n"
            f"    peptide_counts, light_cutoff={light_cutoff}, heavy_cutoff={heavy_cutoff}\n"
            f")\n"
            "filtered_peptide_counts.head()"
        ),

        # --- Selected peptides ---
        code(
            "selected_peptides_report, selected_peptides = "
            "kdm.report_peptide_protein_summary(filtered_peptide_counts)"
        ),

        # --- Merge with SDRF ---
        code(
            "skyline_merge_obj = kdm.MergeFiles(skyline_data, sdrf_data_file, selected_peptides)\n"
            "skyline_merge = skyline_merge_obj.merge_files()"
        ),

        # --- Pool data ---
        code(
            f"skyline_pool = skyline_merge_obj.select_pool_data(\n"
            f"    col_sample='characteristics[Sample]', sample_value='{pool_value}'\n"
            f")"
        ),

        # --- Pool boxplot ---
        code("kdm.plot_pool_boxplot(skyline_pool)"),

        # --- Intra-plate CV ---
        code(
            "peptide_plate_stats = kdm.calculate_intra_plate_cv(\n"
            "    skyline_pool, col_name='characteristics[plate]'\n"
            ")\n"
            "kdm.plot_intra_plate_cv_stats(peptide_plate_stats, col_name='characteristics[plate]')"
        ),

        # --- Inter-plate CV ---
        code(
            "interplate_cv = kdm.calculate_inter_plate_cv(peptide_plate_stats)\n"
            "kdm.plot_inter_plate_cv_kde(interplate_cv)"
        ),

        # --- Inter-plate CV table ---
        code(
            "interplate_cv.head()\n"
            "print('Inter-plate CV (low to high):')\n"
            "print(interplate_cv[interplate_cv['inter_plate_cv'] < 0.1])"
        ),

        # --- Cumulative CV ---
        code("kdm.plot_cumulative_peptide_count_by_cv(interplate_cv)"),

        # --- Normalization peptides ---
        code(
            f"selected_norm_peptides = kdm.get_lowest_cv_peptides(interplate_cv, {cv_percentile})\n"
            "print(selected_norm_peptides)"
        ),

        # --- Pool selected df ---
        code(
            "pool_selected_df = skyline_pool[\n"
            "    skyline_pool['Peptide Sequence'].isin(selected_norm_peptides)\n"
            "]\n"
            "pool_selected_df"
        ),

        # --- Pool PCA (full merge) ---
        code("kdm.plot_pool_pca(skyline_merge)"),

        # --- Pool PCA (norm peptides only) ---
        code("kdm.plot_pool_pca(pool_selected_df)"),

        # --- Pool boxplot (norm peptides) ---
        code("kdm.plot_pool_boxplot(pool_selected_df)"),

        # --- Plate conversion factors ---
        code(
            "plate_factor_table, conversion_factors, model = kdm.get_plate_conversion_factors(\n"
            "    pool_selected_df, col_plate='characteristics[plate]', log_transform=True\n"
            ")"
        ),

        # --- Batch effect test ---
        code(
            "batch_dict = kdm.detect_batch_effect(\n"
            "    pool_selected_df, col_plate='characteristics[plate]',\n"
            "    ratio_col='RatioLightToHeavy', log_transform=True\n"
            ")"
        ),

        # --- Show factor tables ---
        code("plate_factor_table"),
        code("conversion_factors"),

        # --- Plot conversion factors ---
        code(
            "kdm.plot_plate_conversion_factors(\n"
            "    pool_selected_df, col_plate='characteristics[plate]', log_transform=True\n"
            ")"
        ),

        # --- Build adjusted dataset ---
        code(
            "skyline_merge_adj = skyline_merge.copy()\n"
            "skyline_merge_adj = skyline_merge_adj[\n"
            "    ~skyline_merge_adj['Replicate'].isin(qc_samples)\n"
            "]"
        ),

        # --- Apply batch correction ---
        code(
            "skyline_merge_adj = kdm.adjust_ratio_by_plate(skyline_merge_adj, conversion_factors)\n"
            "skyline_merge_adj.head()"
        ),

        # --- Pool adjusted data ---
        code(
            f"pool_data_adj = skyline_merge_adj[\n"
            f"    skyline_merge_adj['characteristics[Sample]'] == '{pool_value}'\n"
            f"].copy()\n"
            "pool_data_adj = pool_data_adj.sort_values(by='characteristics[plate]')\n"
            "pool_data_adj.head()"
        ),

        # --- Pool PCA (adjusted) ---
        code("kdm.plot_pool_pca(pool_data_adj)"),

        # --- Pool boxplot (adjusted) ---
        code("kdm.plot_pool_boxplot(pool_data_adj)"),

        # --- Verify normalization peptides after correction ---
        code(
            "pool_selected_adj = pool_data_adj[\n"
            "    pool_data_adj['Peptide Sequence'].isin(selected_norm_peptides)\n"
            "]\n"
            "kdm.plot_plate_conversion_factors(\n"
            "    pool_selected_adj, col_plate='characteristics[plate]', log_transform=True\n"
            ")"
        ),

        # --- Absolute quantification section ---
        md("# Calculate absolute quantification"),

        # --- Fetch qRePS ---
        code(
            "qreps_lot = sdrf_data.extract_qreps_lot_number()\n"
            "qreps_table = kdm.fetch_qreps_table(qreps_lot)\n"
            "qreps_table.head()"
        ),

        # --- Inspect adjusted data ---
        code("skyline_merge_adj.head()"),

        # --- Pool PCA on adjusted merge ---
        code("kdm.plot_pool_pca(skyline_merge_adj)"),

        # --- Absolute concentrations ---
        code(
            "abs_df = kdm.get_absolute_conc(qreps_table, skyline_merge_adj)\n"
            "abs_df.head()"
        ),
    ]

    nb.cells = cells
    return nb


def main():
    parser = argparse.ArgumentParser(
        prog="kardemumma-ipynb",
        description=(
            "Generate a ratio analysis Jupyter notebook with pre-filled file paths "
            "and parameters. The output notebook mirrors the ratio_MORPHEUS analysis flow."
        ),
    )
    parser.add_argument("--skyline", required=True, metavar="CSV",
                        help="Path to Skyline CSV export")
    parser.add_argument("--sdrf", required=True, metavar="TSV",
                        help="Path to SDRF .tsv file")
    parser.add_argument("--output", required=True, metavar="DIR",
                        help="Directory to write the notebook into")
    parser.add_argument("--name", default="ratio_analysis.ipynb", metavar="FILENAME",
                        help="Notebook filename (default: ratio_analysis.ipynb)")
    parser.add_argument("--dotp", type=float, default=0.6, metavar="THRESHOLD",
                        help="Library dot product threshold (default: 0.6)")
    parser.add_argument("--light-cutoff", type=int, default=700, metavar="N",
                        help="Min light count per peptide (default: 700)")
    parser.add_argument("--heavy-cutoff", type=int, default=700, metavar="N",
                        help="Min heavy count per peptide (default: 700)")
    parser.add_argument("--pool-value", default="Pool", metavar="STR",
                        help="Value in characteristics[Sample] marking pool samples "
                             "(default: Pool)")
    parser.add_argument("--cv-percentile", type=float, default=10.0, metavar="PCT",
                        help="Inter-plate CV percentile for normalization peptides (default: 10)")
    args = parser.parse_args()

    try:
        import nbformat
    except ImportError:
        print("Error: nbformat is required. Install with: pip install nbformat",
              file=sys.stderr)
        sys.exit(1)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    nb = _make_notebook(
        skyline_path=args.skyline,
        sdrf_path=args.sdrf,
        dotp=args.dotp,
        light_cutoff=args.light_cutoff,
        heavy_cutoff=args.heavy_cutoff,
        pool_value=args.pool_value,
        cv_percentile=args.cv_percentile,
    )

    nb_path = out / args.name
    with open(nb_path, "w") as fh:
        nbformat.write(nb, fh)

    print(f"Notebook written → {nb_path}")
    print()
    print("Parameters embedded:")
    print(f"  Skyline       : {args.skyline}")
    print(f"  SDRF          : {args.sdrf}")
    print(f"  dotp          : {args.dotp}")
    print(f"  light/heavy   : {args.light_cutoff} / {args.heavy_cutoff}")
    print(f"  pool value    : {args.pool_value}")
    print(f"  cv percentile : {args.cv_percentile}")
    print()
    print("Open with:  jupyter lab " + str(nb_path))
