#' Process LC-MS feature table
#' @description
#' This function processes the feature metadata and intensity matrix from xcms.
#' It renames, filters, and adds summary intensity for downstream annotation workflows.
#'
#' @param result LC-MS result list, typically `result[[2]]` is metadata and `result[[3]]` s intensity.
#' @param rt_range Optional numeric vector of length 2 for RT filtering (e.g. c(60,1200)).
#' @param sort_by Optional column name to sort by.
#' @param sample_cols Either a regex string (e.g. "mzxml") or a character vector of sample column names.
#'
#' @return Processed feature table with intensity_mean column and optional RT filtering/sorting.
#' @export
processFeatureTable <- function(result,
                                rt_range = NULL,
                                sort_by = NULL,
                                sample_cols = "mzxml") {
  stopifnot(length(result) >= 3)

  meta_df <- result[[2]] %>%
    as.data.frame() %>%
    dplyr::select(-peakidx) %>%
    tibble::rownames_to_column("feature_name")

  int_df <- result[[3]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_name")

  sample_expr <- if (is.character(sample_cols) && length(sample_cols) == 1) {
    rlang::expr(dplyr::matches(!!sample_cols))
  } else {
    rlang::expr(dplyr::all_of(!!sample_cols))
  }

  int_df <- int_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(intensity_mean = mean(dplyr::c_across(!!sample_expr), na.rm = TRUE)) %>%
    dplyr::ungroup()

  joined_df <- dplyr::left_join(meta_df, int_df, by = "feature_name")

  if (!is.null(rt_range) && length(rt_range) == 2) {
    joined_df <- dplyr::filter(joined_df, rtmed >= rt_range[1], rtmed <= rt_range[2])
  }

  if (!is.null(sort_by) && sort_by %in% names(joined_df)) {
    joined_df <- dplyr::arrange(joined_df, !!rlang::sym(sort_by))
  }

  return(joined_df)
}
