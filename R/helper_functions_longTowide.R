#' Convert long annotations to wide or collapsed format
#'
#' @param df A data frame with annotations.
#' @param col_name1 The feature name column (e.g., Feature_name).
#' @param col_name2 The compound name column (e.g., Comp_name_PAMDB).
#' @param col_name3 The ion name + exact mass column (e.g., PI_name_PAMDB).
#' @param collapse Logical. If TRUE (default), collapse all annotations into one column separated by "; ".
#'
#' @return A data frame with wide or collapsed annotations.
#' @export
long_to_wide.fun <- function(df, col_name1, col_name2, col_name3, collapse = TRUE) {
  # Combine compound and PI into a single annotation
  df <- df %>%
    dplyr::mutate(annotation_combined = paste({{col_name2}}, {{col_name3}})) %>%
    dplyr::distinct({{col_name1}}, annotation_combined) # remove duplicates

  # Wide format
  wide_df <- df %>%
    dplyr::group_by({{col_name1}}) %>%
    dplyr::mutate(
      var_name = paste0("annotation_", dplyr::row_number())
    ) %>%
    tidyr::pivot_wider(names_from = var_name, values_from = annotation_combined) %>%
    dplyr::ungroup() %>%
    dplyr::rename(feature_name = {{col_name1}})

  if (collapse) {
    # Collapse annotations, remove duplicates
    collapsed_df <- wide_df %>%
      tidyr::unite(
        col = "targetList_anno",
        tidyselect::starts_with("annotation_"),
        sep = "; ",
        na.rm = TRUE
      ) %>%
      dplyr::mutate(targetList_anno = sapply(
        strsplit(targetList_anno, "; "),
        function(x) paste(unique(x), collapse = "; ")
      ))
    return(collapsed_df)
  } else {
    return(wide_df)
  }
}

