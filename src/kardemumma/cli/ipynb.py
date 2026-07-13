import argparse
import sys
from pathlib import Path


def _make_notebook(skyline_path, sdrf_path, dotp, light_cutoff, heavy_cutoff,
                   pool_value, group_a, group_b, group_col, id_col):
    import nbformat

    nb = nbformat.v4.new_notebook()
    code = nbformat.v4.new_code_cell
    md = nbformat.v4.new_markdown_cell

    group_col_kwarg = f",\n    group_col={group_col!r}" if group_col else ""

    cells = [
        # --- Install ---
        code("%pip install kardemumma"),

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
            f"skyline_path = {skyline_path!r}\n\n"
            "skyline_importer = kdm.ImportSkylineFile(skyline_path)\n"
            "skyline_data = skyline_importer.import_skyline_file()"
        ),

        # --- SDRF readout ---
        code(
            f"sdrf_path = {sdrf_path!r}\n\n"
            "kdm.readout_ms_type(sdrf_path)"
        ),

        # --- SDRF import + file comparison ---
        code(
            "sdrf_data = kdm.ImportSDRFFile(sdrf_path)\n"
            "sdrf_data_file = sdrf_data.import_sdrf_file()\n\n"
            "not_in_sdrf = set(skyline_data['File Name']) - set(sdrf_data_file['comment[data file]'])\n"
            "not_in_skyline = set(sdrf_data_file['comment[data file]']) - set(skyline_data['File Name'])\n"
            "print('Sample names in the Skyline file but not in the SDRF file:', not_in_sdrf)\n"
            "print('Sample names in the SDRF file but not in the Skyline file:', not_in_skyline)"
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
            f"    label_col='characteristics[Sample]', pool_label={pool_value!r}\n"
            f")"
        ),

        # --- Pool boxplot ---
        code("kdm.plot_pool_boxplot(skyline_pool)"),

        # --- Batch correction section ---
        md(
            "## Batch correction: PERMANOVA-driven per-peptide correction (batch_correct.py)\n\n"
            "Test every `characteristics[*]` column for batch effects with PERMANOVA, then use "
            "`correct_ratio_by_factors` to median-center ratios per peptide for each significant "
            "factor (sequentially)."
        ),

        # --- PERMANOVA test ---
        code(
            "# PERMANOVA test for batch effects across all characteristics columns\n"
            "permanova_results = kdm.permanova_batch_effects(skyline_pool)\n\n"
            "# Suggested factors to correct for (p < 0.05)\n"
            'sig_factors = permanova_results[permanova_results["p_value"] < 0.05]["variable"].tolist()\n'
            'print("Suggested factors to correct for batch effect (p < 0.05):", sig_factors)\n\n'
            "permanova_results"
        ),

        # --- Apply correction ---
        code(
            "# Correct ratios by removing identified batch factors for pool data\n"
            "skyline_pool_corrected = kdm.correct_ratio_by_factors(\n"
            "    skyline_pool,\n"
            "    factors=sig_factors,\n"
            '    col_ratio="RatioLightToHeavy",\n'
            ")\n\n"
            "# Correct ratio for all data\n"
            "skyline_corrected = kdm.correct_ratio_by_factors(\n"
            "    skyline_merge,\n"
            "    factors=sig_factors,\n"
            '    col_ratio="RatioLightToHeavy",\n'
            ")"
        ),

        # --- Corrected PCA ---
        code('kdm.plot_pool_pca(skyline_corrected, col_ratio="RatioLightToHeavy_corrected")'),

        # --- Verify correction ---
        code(
            "# Verify batch correction — PERMANOVA on corrected ratios\n"
            "permanova_corrected = kdm.permanova_batch_effects(\n"
            "    skyline_pool_corrected,\n"
            '    col_ratio="RatioLightToHeavy_corrected",\n'
            ")\n\n"
            "# Visual check\n"
            'kdm.plot_pool_pca(skyline_pool_corrected, col_ratio="RatioLightToHeavy_corrected")'
        ),

        # --- Build adjusted dataset from the PERMANOVA-corrected ratios ---
        code(
            "# get_absolute_conc reads from 'RatioLightToHeavy', so carry the\n"
            "# PERMANOVA-corrected ratio into that column for the downstream steps.\n"
            "skyline_merge_adj = skyline_corrected.copy()\n"
            "# Some QC replicates can still be present in skyline_merge even after the\n"
            "# earlier remove_qc_samples() call -- scrub them again before absolute quant.\n"
            "skyline_merge_adj = skyline_merge_adj[~skyline_merge_adj['Replicate'].isin(qc_samples)]\n"
            'skyline_merge_adj["RatioLightToHeavy"] = skyline_merge_adj["RatioLightToHeavy_corrected"]\n'
            "skyline_merge_adj.head()"
        ),

        # --- Absolute quantification section ---
        md("# Calculate absolute quantification"),

        # --- Fetch qRePS ---
        code(
            "qreps_lot = sdrf_data.extract_qreps_lot_number()\n"
            "qreps_table = kdm.fetch_qreps_table(qreps_lot)\n"
            "qreps_table.head()"
        ),

        # --- Absolute concentrations ---
        code(
            "abs_df = kdm.get_absolute_conc(qreps_table, skyline_merge_adj)\n"
            "abs_df.head()"
        ),

        # --- Downstream analysis section ---
        md(
            "# Downstream analysis\n\n"
            "Differential expression and pathway/gene-set enrichment on the absolute peptide "
            "concentrations computed above (`abs_df`)."
        ),
        md(
            "## Differential expression\n\n"
            f"Compare peptide-level absolute concentrations between two SDRF groups "
            f"(`{group_a}` vs `{group_b}`) with `DownStream` "
            "(Welch's t-test per peptide, Benjamini-Hochberg adjusted)."
        ),

        # --- DE ---
        code(
            "# Differential expression between two groups of the SDRF primary variable column\n"
            "# (defaults to the last 'factor value[...]' column). abs_df is the wide output\n"
            "# of get_absolute_conc above -- DownStream handles melting/merging with sdrf internally.\n"
            "de = kdm.DownStream(\n"
            "    abs_df, sdrf_data_file,\n"
            f"    group_a={group_a!r}, group_b={group_b!r},\n"
            f"    id_col={id_col!r}{group_col_kwarg}\n"
            ")\n"
            "de.results.head()"
        ),

        # --- Volcano plot ---
        code("de.plot_volcano()"),

        # --- Enrichment section ---
        md(
            "## Pathway / gene-set enrichment\n\n"
            "Test the significant peptide hits from the differential expression step against "
            "pathway/gene-set definitions with a hypergeometric (Fisher's exact) test. Replace "
            "the `pathways` placeholder below with a real gene-set database (e.g. Reactome, "
            "KEGG, GO) keyed by peptide/protein ID."
        ),

        # --- Enrichment test ---
        code(
            f'sig_hits = de.results.loc[de.results["padj"] < 0.05, {id_col!r}].tolist()\n'
            f'background = de.results[{id_col!r}].tolist()\n'
            'print(f"{len(sig_hits)} significant peptides (padj < 0.05) out of {len(background)} tested")\n\n'
            "# TODO: replace with a real pathway / gene-set database mapping\n"
            "# pathway name -> list of member peptide/protein IDs drawn from `background`.\n"
            "pathways = {\n"
            '    "Example pathway A": background[:5],\n'
            '    "Example pathway B": background[5:10],\n'
            "}\n\n"
            "if sig_hits:\n"
            "    enrichment_results = kdm.enrichment_analysis(\n"
            "        hit_ids=sig_hits,\n"
            "        pathways=pathways,\n"
            "        background_ids=background,\n"
            "        min_overlap=1,\n"
            "    )\n"
            "    display(enrichment_results)\n"
            "else:\n"
            "    enrichment_results = None\n"
            '    print("No significant hits at padj < 0.05 -- skip enrichment test.")'
        ),

        # --- Enrichment plot ---
        code(
            "if enrichment_results is not None and not enrichment_results.empty:\n"
            "    kdm.plot_enrichment(enrichment_results)"
        ),
    ]

    nb.cells = cells
    return nb


def main():
    parser = argparse.ArgumentParser(
        prog="kardemumma-ipynb",
        description=(
            "Generate a ratio analysis Jupyter notebook with pre-filled file paths "
            "and parameters. The output notebook mirrors the ratio_MORPHEUS analysis "
            "flow, using PERMANOVA-driven per-peptide batch correction "
            "(kardemumma.batch_correct) followed by differential expression / "
            "pathway enrichment on the absolute quantification."
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
    parser.add_argument("--light-cutoff", type=int, default=0, metavar="N",
                        help="Min light count per peptide (default: 0)")
    parser.add_argument("--heavy-cutoff", type=int, default=0, metavar="N",
                        help="Min heavy count per peptide (default: 0)")
    parser.add_argument("--pool-value", default="Pool", metavar="STR",
                        help="Value in characteristics[Sample] marking pool samples "
                             "(default: Pool)")
    parser.add_argument("--group-a", required=True, metavar="LABEL",
                        help="First group label to compare for differential expression "
                             "(value in the SDRF group column)")
    parser.add_argument("--group-b", required=True, metavar="LABEL",
                        help="Second group label to compare for differential expression")
    parser.add_argument("--group-col", default=None, metavar="COLUMN",
                        help="SDRF column to compare group_a/group_b on. Defaults to the "
                             "last 'factor value[...]' column")
    parser.add_argument("--id-col", default="Peptide Sequence", metavar="COLUMN",
                        help="Feature column to test in differential expression: "
                             "'Peptide Sequence' or 'Protein Name' (default: Peptide Sequence)")
    args = parser.parse_args()

    try:
        import nbformat
    except ImportError:
        print("Error: nbformat is required. Install with: pip install nbformat",
              file=sys.stderr)
        sys.exit(1)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    skyline_path = str(Path(args.skyline).resolve())
    sdrf_path = str(Path(args.sdrf).resolve())

    nb = _make_notebook(
        skyline_path=skyline_path,
        sdrf_path=sdrf_path,
        dotp=args.dotp,
        light_cutoff=args.light_cutoff,
        heavy_cutoff=args.heavy_cutoff,
        pool_value=args.pool_value,
        group_a=args.group_a,
        group_b=args.group_b,
        group_col=args.group_col,
        id_col=args.id_col,
    )

    nb_path = out / args.name
    with open(nb_path, "w") as fh:
        nbformat.write(nb, fh)

    print(f"Notebook written → {nb_path}")
    print()
    print("Parameters embedded:")
    print(f"  Skyline       : {skyline_path}")
    print(f"  SDRF          : {sdrf_path}")
    print(f"  dotp          : {args.dotp}")
    print(f"  light/heavy   : {args.light_cutoff} / {args.heavy_cutoff}")
    print(f"  pool value    : {args.pool_value}")
    print(f"  group a / b   : {args.group_a} / {args.group_b}")
    print(f"  group col     : {args.group_col or '(auto: last factor value[...] column)'}")
    print(f"  id col        : {args.id_col}")
    print()
    print("Open with:  jupyter lab " + str(nb_path))
