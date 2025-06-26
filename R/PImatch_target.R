#' Match Features to Target List Using PI Ion
#'
#' @param FT DataFrame of features with columns `feature_name`, `mz`
#' @param Comp_data DataFrame of compounds with `Name`, `MONO_MASS`, and `[M+H]+` or `[M-H]-` columns
#' @param PIon A character vector (e.g., "[M+H]+")
#' @param diff_mz_thr Absolute m/z threshold
#' @param polarity Ion polarity: "pos" or "neg"
#' @param ppm PPM tolerance
#'
#' @return A list of matched features with annotations
#' @export
PImatch_target <- function(FT, Comp_data, PIon = "[M+H]+", diff_mz_thr = 0.005, polarity = "pos", ppm = 10) {
  PImatch_target_rcpp(FT, Comp_data, PIon, diff_mz_thr, polarity, ppm)
}
