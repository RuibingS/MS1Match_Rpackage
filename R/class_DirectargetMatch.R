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
