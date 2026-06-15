"""
KARDEMUMMA — Key Analysis of Reproducible Data for Efficient Monitoring
in Unified Mass Spectrometry Methods and Assays.

Python package for processing and quality-checking targeted mass spectrometry
outputs (Skyline / OpenSWATH-style exports).
"""

from .client import SkylineClient
from .importer import (
    ImportSkylineFile,
    ImportSDRFFile,
    MergeFiles,
    get_irt_peptides,
    normalize_data_filename,
    remove_qc_samples,
)
from .sdrf import (
    validate_sdrf,
    readout_ms_type,
    csv_to_tsv,
    remove_whitespace,
    detect_trailing_whitespace,
)
from .prm import (
    # Quality Filtering
    filter_library_dot_product,
    filter_peptide_counts,
    # Peptide Detection Summary
    summarise_peptide_counts,
    report_peptide_protein_summary,
    # CV Analysis
    calculate_intra_plate_cv,
    calculate_inter_plate_cv,
    get_lowest_cv_peptides,
    extract_top_percentile,
    # Plate Normalization
    plate_peptide_anova,
    get_plate_conversion_factors,
    detect_batch_effect,
    plot_plate_conversion_factors,
    adjust_ratio_by_plate,
    get_absolute_conc,
    # General QC
    compute_cv,
    flag_missing_values,
    dot_product_summary,
    retention_time_deviation,
    summarize_prm,
    # Plots
    plot_library_dot_product_distribution,
    plot_heavy_light_scatter,
    plot_peptide_counts,
    plot_pool_boxplot,
    plot_pool_pca,
    plot_pool_heatmap,
    plot_intra_plate_cv_stats,
    plot_inter_plate_cv_kde,
    plot_cumulative_peptide_count_by_cv,
    plot_logratio_by_plate_boxplot,
    # Absolute Concentration Plots
    map_peptide_sequence,
    plot_peptide_concentration_by_group,
    plot_peptide_all,
    plot_median_peptide_concentration_by_group,
    plot_all_median_peptide_concentration_by_group,
    plot_all_peptide_concentration_by_group,
    plot_all_all,
)
from .openswath import (
    import_openswath_file,
    filter_best_peak_group,
    plot_dotprod_kde,
    pair_ions_matching,
    remove_precursor,
    count_ions_channel,
    plot_ions_channel,
    filter_ions_channel,
    get_ratio,
)
from .proteomedge import (
    fetch_qreps_table,
    extract_lot_number,
    summarise_qRePs,
    load_qRePs,
    load_qRePs_to_csv,
    fetch_fasta,
    save_fasta,
)
from .uniprot import (
    get_protein_sequences_batch,
    get_swissprot_sequences_batch,
    query_human_proteome,
    write_fasta,
    create_human_proteome_fasta,
    validate_fasta,
)

__all__ = [
    # Client
    "SkylineClient",
    # Import
    "ImportSkylineFile",
    "ImportSDRFFile",
    "MergeFiles",
    "get_irt_peptides",
    "normalize_data_filename",
    "remove_qc_samples",
    # SDRF
    "validate_sdrf",
    "readout_ms_type",
    "csv_to_tsv",
    "remove_whitespace",
    "detect_trailing_whitespace",
    # PRM — quality filtering
    "filter_library_dot_product",
    "filter_peptide_counts",
    # PRM — peptide detection
    "summarise_peptide_counts",
    "report_peptide_protein_summary",
    # PRM — CV analysis
    "calculate_intra_plate_cv",
    "calculate_inter_plate_cv",
    "get_lowest_cv_peptides",
    "extract_top_percentile",
    # PRM — plate normalisation
    "plate_peptide_anova",
    "get_plate_conversion_factors",
    "detect_batch_effect",
    "plot_plate_conversion_factors",
    "adjust_ratio_by_plate",
    "get_absolute_conc",
    # PRM — general QC
    "compute_cv",
    "flag_missing_values",
    "dot_product_summary",
    "retention_time_deviation",
    "summarize_prm",
    # PRM — plots
    "plot_library_dot_product_distribution",
    "plot_heavy_light_scatter",
    "plot_peptide_counts",
    "plot_pool_boxplot",
    "plot_pool_pca",
    "plot_pool_heatmap",
    "plot_intra_plate_cv_stats",
    "plot_inter_plate_cv_kde",
    "plot_cumulative_peptide_count_by_cv",
    "plot_logratio_by_plate_boxplot",
    # PRM — absolute concentration plots
    "map_peptide_sequence",
    "plot_peptide_concentration_by_group",
    "plot_peptide_all",
    "plot_median_peptide_concentration_by_group",
    "plot_all_median_peptide_concentration_by_group",
    "plot_all_peptide_concentration_by_group",
    "plot_all_all",
    # OpenSWATH
    "import_openswath_file",
    "filter_best_peak_group",
    "plot_dotprod_kde",
    "pair_ions_matching",
    "remove_precursor",
    "count_ions_channel",
    "plot_ions_channel",
    "filter_ions_channel",
    "get_ratio",
    # ProteomEdge
    "fetch_qreps_table",
    "extract_lot_number",
    "summarise_qRePs",
    "load_qRePs",
    "load_qRePs_to_csv",
    "fetch_fasta",
    "save_fasta",
    # UniProt
    "get_protein_sequences_batch",
    "get_swissprot_sequences_batch",
    "query_human_proteome",
    "write_fasta",
    "create_human_proteome_fasta",
    "validate_fasta",
]
