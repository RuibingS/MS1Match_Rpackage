#' Compute Pairwise Feature Correlation Using Rcpp
#'
#' @param feature_table A data.frame with intensity values and a `feature_name` column
#' @param intensity_cols Column indices for the numeric intensity values (e.g., 3:20)
#' @param feature_pairs A data.frame with columns `feature_1` and `feature_2`
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#'
#' @return A data.frame with feature_1, feature_2, and cor_value columns
#' @export
CorIntensityPair <- function(feature_table, intensity_cols, feature_pairs, method = "pearson") {
  stopifnot("feature_name" %in% colnames(feature_table))

  method <- tolower(method)
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("Invalid method: choose 'pearson', 'spearman', or 'kendall'")
  }


  intensity_matrix <- as.matrix(feature_table[, intensity_cols, drop = FALSE])
  intensity_matrix[intensity_matrix == 0] <- NA
  log_intensity_matrix <- log10(intensity_matrix)


  log_intensity_t <- t(log_intensity_matrix)

  cor_result <- CorPairs_cpp(
    log_intensity_t,
    feature_table$feature_name,
    feature_pairs,
    method = method
  )

  return(cor_result)
}
