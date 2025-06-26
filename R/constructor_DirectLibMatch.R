#' Create a DirectLibMatch object
#'
#' @param FT A data.frame of features
#' @param compounds_list A list of compounds
#' @param diff_mz_thr Numeric m/z difference threshold
#' @param diff_rt_thr Numeric RT difference threshold
#' @param ppm Numeric ppm tolerance
#' @param ms_level Either "MS1" or "MS2"
#' @param polarity Either "pos" or "neg"
#'
#' @return A DirectLibMatch S4 object
#' @export
createDirectLibMatch <- function(FT, compounds_list, diff_mz_thr, diff_rt_thr, ppm,
                                 ms_level = "MS1", polarity = "pos") {
  new("DirectLibMatch",
      FT = FT,
      compounds_list = compounds_list,
      diff_mz_thr = diff_mz_thr,
      diff_rt_thr = diff_rt_thr,
      ppm = ppm,
      ms_level = ms_level,
      polarity = polarity)
}
