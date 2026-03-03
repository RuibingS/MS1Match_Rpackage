#' Compute EIC Correlation for Selected Pairs
#'
#' @description Computes pairwise similarity between selected feature pairs using processed EICs.
#'
#' @param object An `EICData` object containing features and parsed EICs
#' @param feature_pairs A data.frame with columns `feature_1` and `feature_2`
#'                      specifying the pairs of features to compute correlations for
#' @param rt_EIC_thr Retention time tolerance for aligning EICs (default: 1)
#' @param method Correlation method: "pearson", "spearman", or "kendall" (default: "pearson")
#'
#' @return Updated `EICData` object with correlation results stored in `@correlation_df`
#' @export
setMethod("eicCorrelation", "EICData", function(object,
                                                feature_pairs,
                                                rt_EIC_thr = 1,
                                                method = "pearson") {

  # Ensure feature_pairs has required columns
  if (!all(c("feature_1", "feature_2") %in% colnames(feature_pairs))) {
    stop("feature_pairs must contain columns 'feature_1' and 'feature_2'")
  }

  # Call the updated C++ function (debug or final version)
  cor_df <- EIC_cor_mat_pairs(
    object@parsed_EIC,
    feature_pairs,
    rt_EIC_thr,
    method
  )

  # Update slots in EICData object
  object@correlation_method <- method
  object@correlation_df <- cor_df

  return(object)
})
