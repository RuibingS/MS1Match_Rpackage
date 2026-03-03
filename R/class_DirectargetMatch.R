#' S4 class for direct target list PI matching
#'
#' @slot FT A data.frame of features
#' @slot Comp_data A data frame of parsed compounds target list
#' @slot diff_mz_thr Numeric m/z difference threshold
#' @slot PIon character `[M+H]+` for positive ion mode. `[M-H]-` for negative ion mode
#' @slot diff_rt_thr Numeric RT difference threshold. If rt exist in the target list, it will skip rt threshold
#' @slot ppm Numeric ppm for m/z tolerance
#' @export
setClass("DirectTargetMatch",
         slots = c(
           FT = "data.frame",
           Comp_data = "data.frame",
           PIon = "character",
           diff_mz_thr = "numeric",
           diff_rt_thr = "numeric",
           ppm = "numeric"


                    ))
#' Constructor for DirectTargetMatch
#'
#' @param FT Feature table (data.frame)
#' @param Comp_data Compound target list (data.frame)
#' @param PIon Ion mode, e.g.,`"[M+H]+"` or `"[M-H]-"`.
#' @param diff_mz_thr m/z tolerance (numeric)
#' @param diff_rt_thr retention time tolerance (numeric), default NA
#' @param ppm m/z tolerance in ppm
#'
#' @return A DirectTargetMatch object
#' @export
DirectTargetMatch <- function(FT, Comp_data, PIon,
                              diff_mz_thr, ppm,
                              diff_rt_thr = NA_real_) {
  new("DirectTargetMatch",
      FT = FT,
      Comp_data = Comp_data,
      PIon = PIon,
      diff_mz_thr = diff_mz_thr,
      diff_rt_thr = diff_rt_thr,
      ppm = ppm)
}
