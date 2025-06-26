#' Join Rcpp annotation results
#'
#' @param dirct_lib_list Required list (library match output)
#' @param lib_anno_list Optional list to join by `Matched_features_to_peaks`
#' @param target_match_list Optional list to join by `Matched_features_to_peaks`
#'
#' @return A data.frame of joined annotation results
#' @export
join_annotation_results <- function(dirct_lib_list, lib_anno_list = NULL, target_match_list = NULL) {
  stopifnot(is.list(dirct_lib_list))

  df_output <- NULL

  dirct_lib_df <- tryCatch({
    dplyr::bind_rows(dirct_lib_list)
  }, error = function(e) {
    stop("Failed to bind dirct_lib_list: ", conditionMessage(e))
  })

  # Validate required columns
  required_input_cols <- c("feature_name", "feature_mz", "feature_rt", "Ion_type")
  missing_cols <- setdiff(required_input_cols, colnames(dirct_lib_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in dirct_lib_df: ", paste(missing_cols, collapse = ", "))
  }

  # Join with lib_anno_list
  if (!is.null(lib_anno_list)) {
    lib_anno_df <- tryCatch({
      dplyr::bind_rows(lib_anno_list)
    }, error = function(e) {
      warning("Failed to bind lib_anno_list: ", conditionMessage(e))
      NULL
    })

    if (!is.null(lib_anno_df)) {
      dirct_lib_df <- tryCatch({
        dplyr::full_join(dirct_lib_df, lib_anno_df, by = c("Matched_features_to_peaks" = "feature_name")) # left_join or full join?
      }, error = function(e) {
        warning("Failed to join dirct_lib_df with lib_anno_df: ", conditionMessage(e))
        dirct_lib_df
      })

      # Rename columns
      dirct_lib_df <- tryCatch({
        dplyr::rename_with(dirct_lib_df, ~ dplyr::case_when(
          .x == "feature_name" ~ "PI_feature_name",
          .x == "feature_mz" ~ "PI_feature_mz",
          .x == "feature_rt" ~ "PI_feature_rt",
          .x == "Ion_type" ~ "PI_Ion_type",
          TRUE ~ .x
        ))
      }, error = function(e) {
        warning("Failed during column renaming in dirct_lib_df: ", conditionMessage(e))
        dirct_lib_df
      })

      expected_cols <- c("PI_feature_name", "PI_feature_mz", "PI_feature_rt", "Compound_name",
                         "Compound_rt", "PI_Ion_type", "Matched_peaks", "Match_peak_score",
                         "Total_library_peaks")
      if (all(expected_cols %in% colnames(dirct_lib_df))) {
        df_output <- tryCatch({
          dirct_lib_df %>%
            dplyr::group_by(PI_feature_name) %>%
            dplyr::mutate(
              PI_feature_mz = dplyr::if_else(dplyr::row_number() == 1, PI_feature_mz, NA_real_),
              PI_feature_rt = dplyr::if_else(dplyr::row_number() == 1, PI_feature_rt, NA_real_),
              Compound_name = dplyr::if_else(dplyr::row_number() == 1, Compound_name, NA_character_),
              Compound_rt = dplyr::if_else(dplyr::row_number() == 1, Compound_rt, NA_real_),
              PI_Ion_type = dplyr::if_else(dplyr::row_number() == 1, PI_Ion_type, NA_character_),
              Matched_peaks = dplyr::if_else(dplyr::row_number() == 1, Matched_peaks, NA_real_),
              Match_peak_score = dplyr::if_else(dplyr::row_number() == 1, Match_peak_score, NA_real_),
              Total_library_peaks = dplyr::if_else(dplyr::row_number() == 1, Total_library_peaks, NA_real_)
            ) %>%
            dplyr::ungroup()
        }, error = function(e) {
          warning("Failed to mutate: ", conditionMessage(e))
          dirct_lib_df
        })
      } else {
        warning("Missing columns: ", paste(setdiff(expected_cols, colnames(dirct_lib_df)), collapse = ", "))
        df_output <- dirct_lib_df
      }
    }
  } else {
    df_output <- dirct_lib_df
  }

  # Join with target_match_list
  if (!is.null(target_match_list)) {
    target_match_df <- tryCatch({
      dplyr::bind_rows(Filter(function(x) {
        !is.null(x) && is.list(x) && "Feature_name" %in% names(x)
      }, target_match_list))
    }, error = function(e) {
      warning("Failed to bind target_match_list: ", conditionMessage(e))
      NULL
    })

    if (!is.null(target_match_df)) {
      if ("Matched_features_to_peaks" %in% colnames(df_output) &&
          "Feature_name" %in% colnames(target_match_df)) {
        df_output <- tryCatch({
          dplyr::left_join(df_output, target_match_df, by = c("Matched_features_to_peaks" = "Feature_name"))
        }, error = function(e) {
          warning("Failed to join df_output with target_match_df: ", conditionMessage(e))
          df_output
        })
      } else {
        warning("Could not join target_match_list: required column(s) missing.")
      }
    }
  }

  return(df_output)
}
