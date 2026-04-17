"""
KARDEMUMMA — Key Analysis of Reproducible Data for Efficient Monitoring
in Unified Mass Spectrometry Methods and Assays.

Python package for processing and quality-checking targeted mass spectrometry
outputs (Skyline / OpenSWATH-style exports).
"""

from .client import SkylineClient
from .importer import (
    ImportFile,
    CheckSkylineFile,
    MergeFiles,
    cross_check_skyline_sdrf,
    get_irt_peptides,
    normalize_data_filename,
    import_sdrf_file,
    remove_qc_samples,  # Added here
)
from .sdrf import validate_sdrf
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
    adjust_ratio_by_plate,
    # General QC
    compute_cv,
    flag_missing_values,
    dot_product_summary,
    retention_time_deviation,
    summarize_prm,
)
from .prm_plots import (
    plot_library_dot_product_distribution,
    plot_heavy_light_scatter,
    plot_peptide_counts,
    plot_pool_boxplot,
    plot_pool_heatmap,
    plot_intra_plate_cv_stats,
    plot_inter_plate_cv_kde,
    plot_cumulative_peptide_count_by_cv,
    plot_logratio_by_plate_boxplot,
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
    load_qRePs,
    load_qRePs_to_csv,
    fetch_fasta,
    save_fasta,
)

__all__ = [
    # Client
    "SkylineClient",
    # Import / merge
    "ImportFile",
    "CheckSkylineFile",
    "MergeFiles",
    "cross_check_skyline_sdrf",
    "get_irt_peptides",
    "normalize_data_filename",
    "import_sdrf_file",
    "remove_qc_samples",  # Added here
    # SDRF
    "validate_sdrf",
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
    "adjust_ratio_by_plate",
    # PRM — general QC
    "compute_cv",
    "flag_missing_values",
    "dot_product_summary",
    "retention_time_deviation",
    "summarize_prm",
    # PRM plots
    "plot_library_dot_product_distribution",
    "plot_heavy_light_scatter",
    "plot_peptide_counts",
    "plot_pool_boxplot",
    "plot_pool_heatmap",
    "plot_intra_plate_cv_stats",
    "plot_inter_plate_cv_kde",
    "plot_cumulative_peptide_count_by_cv",
    "plot_logratio_by_plate_boxplot",
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
    "load_qRePs",
    "load_qRePs_to_csv",
    "fetch_fasta",
    "save_fasta",
]
