# Package index

## All functions

- [`annotateISF(`*`<PeakAnnotationISF>`*`)`](https://ruibings.github.io/MS1Match_Rpackage/reference/annotateISF-PeakAnnotationISF-method.md)
  : Annotate features as MS2 fragment ions (Tier 1c / 1c_exp)
- [`annotateLibMatch(`*`<PeakAnnotationFromLibraryMatch>`*`)`](https://ruibings.github.io/MS1Match_Rpackage/reference/annotateLibMatch-PeakAnnotationFromLibraryMatch-method.md)
  : Annotate features from spectral library match results
- [`annotateLibMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/annotateLibMatch.md)
  : Annotate peaks
- [`annotations_df()`](https://ruibings.github.io/MS1Match_Rpackage/reference/annotations_df.md)
  : Flatten annotate result (PeakAnnotation) to a data.frame
- [`assign_feature_names()`](https://ruibings.github.io/MS1Match_Rpackage/reference/assign_feature_names.md)
  : Give a grouping result the feature ids of the MS1Match feature table
- [`buildFeatureTable()`](https://ruibings.github.io/MS1Match_Rpackage/reference/buildFeatureTable.md)
  : Build simple feature table
- [`build_anchors()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md)
  : Define the M+H+ anchor feature for each compound instance
- [`build_annot_lookup()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_annot_lookup.md)
  : Build a per-feature annotation lookup table
- [`build_anno_accuracy()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md)
  : Assess Tier 1a adduct annotation accuracy
- [`build_benchmark_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table.md)
  : Build the tier-stratified benchmark summary table
- [`build_benchmark_table_from_full()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table_from_full.md)
  : Recompute the tier-stratified benchmark table from
  build_full_feature_table()
- [`build_compound_eval()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_eval.md)
  : Build the per-compound tier-wise grouping recall table
- [`build_compound_instance_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_instance_table.md)
  : Build a flat, all-compounds feature table with CAMERA grouping
  context
- [`build_eval_full()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md)
  : Build the feature × compound evaluation table
- [`build_full_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_full_feature_table.md)
  : Build the complete CAMERA feature table with reference-relative
  status
- [`build_log10_yaxis()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_log10_yaxis.md)
  : Build y-axis breaks and labels for the log10(x+1) intensity scale
- [`build_mirror_df()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_mirror_df.md)
  : Build the mirror plot data frame for one spectrum pair
- [`build_pcgroup_purity()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_pcgroup_purity.md)
  : Assess CAMERA pcgroup purity (over-grouping)
- [`build_plot2_data()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_plot2_data.md)
  : Build per-pcgroup purity summary for the over-grouping scatter plot
- [`build_query_label()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_query_label.md)
  : Build the label string for a single query peak
- [`calculate_exact_mass()`](https://ruibings.github.io/MS1Match_Rpackage/reference/calculate_exact_mass.md)
  : Calculate monoisotopic exact mass from molecular formulae
- [`camera_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/camera_to_grouping.md)
  : Grouping table from a CAMERA peak list
- [`check_feature_match()`](https://ruibings.github.io/MS1Match_Rpackage/reference/check_feature_match.md)
  : Check that MS1Match and a grouping tool used the same feature table
- [`classify_ref_status()`](https://ruibings.github.io/MS1Match_Rpackage/reference/classify_ref_status.md)
  : Classify features by reference-relative evaluation status
- [`compare_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/compare_grouping.md)
  : Compare several grouping evaluations side by side
- [`compounds()`](https://ruibings.github.io/MS1Match_Rpackage/reference/compounds.md)
  : Access library compounds
- [`compute_xlim()`](https://ruibings.github.io/MS1Match_Rpackage/reference/compute_xlim.md)
  : Compute a shared x-axis range for the MS1 + MS2 mirror plot pair
- [`convert_to_nist_msp()`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md)
  : Convert a parsed library to standard NIST MSP format
- [`createDirectLibMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createDirectLibMatch.md)
  : Create a DirectLibMatch object
- [`createDirectTargetMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createDirectTargetMatch.md)
  : Create a DirectTargetMatch object
- [`createPeakAnnotationFromLibraryMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationFromLibraryMatch.md)
  : Create a PeakAnnotationFromLibraryMatch object
- [`createPeakAnnotationISF()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationISF.md)
  : Create a PeakAnnotationISF object
- [`DirectLibMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/DirectLibMatch-class.md)
  : Class for Direct Library Matching
- [`DirectTargetMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/DirectTargetMatch-class.md)
  : Class for Direct Target Matching
- [`.strip_one_adduct()`](https://ruibings.github.io/MS1Match_Rpackage/reference/dot-strip_one_adduct.md)
  : Strip trailing mass value from a single adduct token
- [`evaluate_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md)
  : Evaluate a feature-grouping tool against the MS1Match reference set
- [`exportFeatureTable()`](https://ruibings.github.io/MS1Match_Rpackage/reference/exportFeatureTable.md)
  : Export feature table to CSV
- [`export_mirror_plots()`](https://ruibings.github.io/MS1Match_Rpackage/reference/export_mirror_plots.md)
  : Batch-export mirror plots for all compounds in a direct_match list
- [`export_msp()`](https://ruibings.github.io/MS1Match_Rpackage/reference/export_msp.md)
  : Export a compound library list to NIST MSP format
- [`extract_nl_label_from_adduct()`](https://ruibings.github.io/MS1Match_Rpackage/reference/extract_nl_label_from_adduct.md)
  : Extract a clean neutral-loss label from a pipeline adduct_match
  string
- [`extract_precursor_mz()`](https://ruibings.github.io/MS1Match_Rpackage/reference/extract_precursor_mz.md)
  : Extract precursor m/z from a direct_match entry
- [`features()`](https://ruibings.github.io/MS1Match_Rpackage/reference/features.md)
  : Access features data
- [`find_by_name()`](https://ruibings.github.io/MS1Match_Rpackage/reference/find_by_name.md)
  : Find the first entry in a list matching a compound name
- [`flatten()`](https://ruibings.github.io/MS1Match_Rpackage/reference/flatten.md)
  : Flatten matched PI result list into a data.frame
- [`getChromPeaks()`](https://ruibings.github.io/MS1Match_Rpackage/reference/getChromPeaks.md)
  : Extract chromPeaks as a plain data.frame
- [`get_ms2_peaks()`](https://ruibings.github.io/MS1Match_Rpackage/reference/get_ms2_peaks.md)
  : Get MS2 peaks with feature_name assigned
- [`hello()`](https://ruibings.github.io/MS1Match_Rpackage/reference/hello.md)
  : Hello, World!
- [`HMDB_target_list_import()`](https://ruibings.github.io/MS1Match_Rpackage/reference/HMDB_target_list_import.md)
  : Import HUMD Metabolite Data and Calculate Adducts
- [`import_msp()`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md)
  : Import a standardized NIST MSP spectral library
- [`inspect_compound()`](https://ruibings.github.io/MS1Match_Rpackage/reference/inspect_compound.md)
  : Inspect a compound's CAMERA grouping performance
- [`is_primary_ion()`](https://ruibings.github.io/MS1Match_Rpackage/reference/is_primary_ion.md)
  : Classify whether a feature qualifies as a primary ion
- [`is_singly_charged()`](https://ruibings.github.io/MS1Match_Rpackage/reference/is_singly_charged.md)
  : Test whether a single adduct token is singly charged
- [`join_annotation_results()`](https://ruibings.github.io/MS1Match_Rpackage/reference/join_annotation_results.md)
  : Join lists of annotation results to library match output
- [`libAnnotate()`](https://ruibings.github.io/MS1Match_Rpackage/reference/libAnnotate.md)
  : Annotate features using library matching and cosine similarity
- [`linkCPtoFeatures()`](https://ruibings.github.io/MS1Match_Rpackage/reference/linkCPtoFeatures.md)
  : Link chromPeaks to feature IDs
- [`long_to_wide.fun()`](https://ruibings.github.io/MS1Match_Rpackage/reference/long_to_wide.fun.md)
  : Convert long annotations to wide or collapsed format
- [`mzmine_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/mzmine_to_grouping.md)
  : Grouping table from an MZmine feature-table export
- [`normalize_query()`](https://ruibings.github.io/MS1Match_Rpackage/reference/normalize_query.md)
  : Normalize raw query feature intensities to 0-100 scale
- [`normalize_target_list()`](https://ruibings.github.io/MS1Match_Rpackage/reference/normalize_target_list.md)
  : Normalize a target list data frame to a standard schema
- [`PAMD_target_list_import()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PAMD_target_list_import.md)
  : Import PAMD Metabolite Data and Calculate Adducts
- [`PeakAnnotationFromLibraryMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch-class.md)
  : PeakAnnotationFromLibraryMatch S4 class
- [`PeakAnnotationFromLibraryMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch.md)
  : Constructor for PeakAnnotationFromLibraryMatch
- [`PeakAnnotationISF-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF-class.md)
  : PeakAnnotationISF S4 class
- [`PeakAnnotationISF()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF.md)
  : Constructor for PeakAnnotationISF
- [`PImatch_df()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PImatch_df.md)
  : Flatten PI target match result list to a data.frame
- [`PImatch_target()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PImatch_target.md)
  : Perform primary-ion matching against a direct target list
- [`plotMirrorSpectra()`](https://ruibings.github.io/MS1Match_Rpackage/reference/plotMirrorSpectra.md)
  : Plot mirror plots for a list of spectral matches
- [`plotMirrorSpectraSingle()`](https://ruibings.github.io/MS1Match_Rpackage/reference/plotMirrorSpectraSingle.md)
  : Plot a mirror plot between a query feature and its reference
  compound
- [`plot_mirror()`](https://ruibings.github.io/MS1Match_Rpackage/reference/plot_mirror.md)
  : Draw one mirror-plot panel (query up, reference down)
- [`plot_one_compound()`](https://ruibings.github.io/MS1Match_Rpackage/reference/plot_one_compound.md)
  : Plot MS1 and MS2 mirror panels stacked for one compound
- [`prepare_grouping_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/prepare_grouping_table.md)
  : Standardize a grouping tool's output for evaluation
- [`processFeatureTable()`](https://ruibings.github.io/MS1Match_Rpackage/reference/processFeatureTable.md)
  : Post-process feature table
- [`ramclustr_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/ramclustr_to_grouping.md)
  : Grouping table from a RAMClustR result
- [`read_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_feature_table.md)
  : Read a pre-processed feature table for MS1Match
- [`read_library()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_library.md)
  : Read any MSP-like spectral library file
- [`read_ms2_consensus()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_ms2_consensus.md)
  : Read MS2 consensus data into a named list
- [`read_mzmine_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md)
  : Read and Filter an MZmine Full Feature Table
- [`results()`](https://ruibings.github.io/MS1Match_Rpackage/reference/results.md)
  : Access results from DirectLibMatch object
- [`results_df()`](https://ruibings.github.io/MS1Match_Rpackage/reference/results_df.md)
  : Flatten DirectLibMatch results to a data.frame
- [`rt_thr_coelution()`](https://ruibings.github.io/MS1Match_Rpackage/reference/rt_thr_coelution.md)
  : Access strict RT co-elution threshold
- [`runxcmsFT()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runxcmsFT.md)
  : Pipeline function feature table pipeline
- [`runXCMSProcessing()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runXCMSProcessing.md)
  : Run XCMS preprocessing
- [`run_isf_annotation()`](https://ruibings.github.io/MS1Match_Rpackage/reference/run_isf_annotation.md)
  : Run MS2 fragment annotation (Phase 1 + Phase 2)
- [`sanitize_filename()`](https://ruibings.github.io/MS1Match_Rpackage/reference/sanitize_filename.md)
  : Extract a string for use as a filename
- [`scale_intensity()`](https://ruibings.github.io/MS1Match_Rpackage/reference/scale_intensity.md)
  : Apply intensity scaling transformation
- [`strip_adduct_mass()`](https://ruibings.github.io/MS1Match_Rpackage/reference/strip_adduct_mass.md)
  : Strip trailing mass values from an adduct annotation string
