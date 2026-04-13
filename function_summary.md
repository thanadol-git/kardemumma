# Function Summary

## `skyline_qc.client`
- `SkylineClient._headers`: Build request headers and include bearer token when an API key is provided.
- `SkylineClient.get`: Execute a GET request against the Skyline API and return parsed JSON.
- `SkylineClient.list_documents`: List available Skyline documents from `/api/documents`.
- `SkylineClient.get_document`: Fetch one document by document id.
- `SkylineClient.fetch_qc_metrics`: Fetch QC metrics for a document.
- `SkylineClient.check_col_names`: Fetch column-name metadata for a document.

## `skyline_qc.importer`
- `get_irt_peptides`: Extract distinct iRT peptide sequences from Skyline-style data.
- `_suggest_qc_replicate_names`: Detect replicate names likely to be QC based on `"qc"` substring matching.
- `normalize_data_filename`: Canonicalize raw file names for Skyline/SDRF cross-comparison.
- `_normalized_file_set`: Convert a filename series into a normalized unique-name set.
- `cross_check_skyline_sdrf`: Compare Skyline and SDRF file-name sets and print overlap/difference report.
- `import_sdrf_file`: Backward-compatible module wrapper that delegates to `ImportFile.import_sdrf_file`.
- `ImportFile.import_skyline_file`: Validate and read Skyline CSV, then annotate isotope label type.
- `ImportFile.suggest_qc_samples`: Suggest suspicious QC-like replicate names from Skyline data.
- `ImportFile.import_qreps_file`: Read qRePS CSV file.
- `ImportFile.import_sdrf_file`: Read SDRF TSV file and print detected QC-like data-file entries.
- `CheckSkylineFile.check_skyline_file`: Validate and read Skyline CSV, sorted by replicate and peptide.
- `CheckSkylineFile.suggest_qc_samples`: Suggest QC-like replicate names from file or provided DataFrame.
- `CheckSkylineFile.get_irt_peptides`: Load Skyline CSV and return iRT peptides.
- `CheckSkylineFile.get_test_samples`: Return non-QC sample names after removing selected samples.
- `CheckSkylineFile.get_test_data`: Filter and sort Skyline data to selected test samples.
- `MergeFiles.merge_files`: Merge Skyline and SDRF tables and derive `characteristics[Plate]` from replicate names.
- `MergeFiles.select_pool_data`: Filter merged data to pool samples and return sorted subset.

## `skyline_qc.sdrf`
- `validate_sdrf`: Run `parse_sdrf validate-sdrf` and return pass/fail with message.
- `readout_ms_type`: Extract unique acquisition-method values from SDRF in first-seen order.
- `csv_to_tsv`: Convert a CSV file to a TSV file.
- `remove_whitespace`: Trim leading/trailing whitespace in all string cells of an SDRF table.
- `detect_trailing_whitespace`: Detect and print cells/columns containing leading or trailing whitespace.

## `skyline_qc.proteomedge`
- `fetch_qreps_table`: Download and parse the qRePS table from ProteomEdge lot pages.
- `extract_lot_number`: Parse a lot number from a lot URL or return the provided lot string.
- `load_qRePs`: Fetch qRePS and save it as a date-prefixed CSV filename.
- `load_qRePs_to_csv`: Alias behavior of `load_qRePs` for explicit CSV-saving usage.

## `skyline_qc.prm`
- `_to_python_scalar`: Coerce nested/non-scalar values into scalars for formula modeling.
- `_as_1d_series`: Return a single Series for a column, even if duplicate columns exist.
- `_patsy_scalar_categorical`: Apply scalar coercion across categorical series values.
- `_formula_clean_frame`: Flatten/drop duplicate columns for statsmodels/patsy compatibility.
- `filter_library_dot_product`: Keep rows above a library dot-product threshold and pivot heavy/light areas.
- `filter_peptide_counts`: Filter peptide summaries by minimum heavy and light count cutoffs.
- `summarise_peptide_counts`: Count non-missing heavy/light measurements per protein-peptide.
- `report_peptide_protein_summary`: Report unique peptide/protein counts and return selected peptide list.
- `calculate_intra_plate_cv`: Compute per-peptide intra-plate CV from pool data.
- `calculate_inter_plate_cv`: Compute inter-plate CV across peptide plate-level means.
- `get_lowest_cv_peptides`: Select peptides at or below a user-defined CV percentile threshold.
- `extract_top_percentile`: Extract low-percentile IDs by a metric and optionally filter a source DataFrame.
- `plate_peptide_anova`: Fit two-way ANOVA (`Plate` + `Peptide`) on ratio (or log-ratio) response.
- `get_plate_conversion_factors`: Estimate plate-level correction factors from median-centered ratios.
- `adjust_ratio_by_plate`: Apply plate-specific correction factors to `RatioLightToHeavy`.
- `compute_cv`: Compute grouped percent CV summaries for a chosen value column.
- `flag_missing_values`: Summarize per-precursor missing/detected rates across replicates.
- `dot_product_summary`: Flag transitions passing/failing library and ratio dot-product thresholds.
- `retention_time_deviation`: Compute observed-vs-predicted RT deviation metrics.
- `summarize_prm`: Run major PRM QC summaries and return combined outputs/metrics.

## `skyline_qc.prm_plots`
- `plot_library_dot_product_distribution`: Plot histogram and KDE for `Library Dot Product`.
- `plot_heavy_light_scatter`: Plot heavy-vs-light count scatter colored by point density.
- `plot_peptide_counts`: Plot heavy peptide count distribution.
- `plot_pool_boxplot`: Plot log-scale ratio boxplots by replicate and plate.
- `plot_pool_heatmap`: Plot heatmap of log ratios for pool samples by peptide and replicate.
- `plot_intra_plate_cv_stats`: Plot intra-plate CV distributions by grouping column.
- `plot_inter_plate_cv_kde`: Plot KDE of inter-plate CV with median marker.
- `plot_cumulative_peptide_count_by_cv`: Plot cumulative peptide count as CV threshold increases.
- `plot_logratio_by_plate_boxplot`: Plot replicate-level log-ratio boxplots colored by plate.

## `openms_qc.openswath`
- `_extract_modification_location`: Identify sequence positions of a given UniMod annotation.
- `_n_term_acetylation_annotation`: Annotate whether peptide carries N-terminal acetylation marker.
- `import_openswath_file`: Validate/read OpenSWATH TSV and annotate isotope labels.
- `filter_best_peak_group`: Filter decoys, best-rank peak groups, and dot-product threshold.
- `remove_precursor`: Keep only precursor rows with selected charge states.
- `plot_dotprod_kde`: Plot KDE of OpenSWATH library dot product by isotope channel.
- `pair_ions_matching`: Pivot heavy/light intensities and compute heavy-to-light ratios.
- `_select_ions_channel`: Select and validate columns needed for ion-channel counting.
- `_summarise_ions_channel`: Count ion-channel observations per peptide/isotope group.
- `count_ions_channel`: Produce wide heavy/light ion counts per peptide/charge.
- `_summarise_ions_channel_count`: Summarize frequency of heavy/light ion-count combinations.
- `plot_ions_channel`: Plot heavy-vs-light ion-count scatter colored by frequency.
- `filter_ions_channel`: Filter rows by minimum heavy and light ion-channel counts.
- `get_ratio`: Aggregate intensity and compute light-to-heavy ratios at peptide or precursor level.

## Re-export Modules
- `skyline_qc.__init__`: Re-exports the main public API (client, import, PRM, and plotting helpers).
- `openms_qc.__init__`: Re-exports the OpenSWATH analysis helpers for package-level imports.
- `skyline_qc.output_test`: Backward-compatible re-export module for legacy notebooks/scripts.
