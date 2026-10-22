#' Create a DirectLibMatch object
#'
#' @param FT A data.frame of features
#' @param compounds_list A list of compounds
#' @param rt_thr_compound Numeric loose RT threshold for precursor gate vs library RT (e.g. 30 s)
#' @param rt_thr_coelution Numeric strict RT threshold for co-elution pool anchored to
#'   the matched precursor feature (e.g. 5 s). Defaults to \code{5}.
#' @param ppm Numeric ppm tolerance
#' @param ms_level Either \code{"MS1"} or \code{"MS2"}
#' @param polarity Either \code{"pos"} or \code{"neg"}
#'
#' @return A DirectLibMatch S4 object
#' @export
createDirectLibMatch <- function(FT,
                                 compounds_list,
                                 rt_thr_compound,
                                 rt_thr_coelution = 5,
                                 ppm,
                                 ms_level = "MS1",
                                 polarity = "pos") {
  new("DirectLibMatch",
      FT                 = FT,
      compounds_list     = compounds_list,
      rt_thr_compound        = rt_thr_compound,
      rt_thr_coelution = rt_thr_coelution,
      ppm                = ppm,
      ms_level           = ms_level,
      polarity           = polarity,
      results            = list()
  )
}
