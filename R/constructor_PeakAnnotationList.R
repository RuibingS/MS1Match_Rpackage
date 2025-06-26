#' Create a PeakAnnotationList object
#'
#' @param  compounds_list A named \code{list} of compound information
#' @param  FT A \code{data.frame} with feature data (must include mzmed, rtmed)
#' @param  adducts_df A \code{data.frame} parsed from the default adduct table
#' @param  neutral_loss_table A \code{data.frame} parsed from the default neutral loss table
#' @param  polarity A character string: "pos" or "neg"
#' @param  intensity A numeric vector of the mean intensity across samples
#' @export
#'

PeakAnnotationList <- function(compounds_list, FT, adducts_df, neutral_loss_table,
                               polarity, intensity,
                               mz_threshold = 0.005, PI_rt_threshold = 5.0,
                               rt_threshold = 3.0, iso_threshold = 1.003355,
                               max_charge = 3L, max_isotope = 4L) {

  if (!is.list(compounds_list)) stop("`compounds_list` must be a list")
  if (!is.data.frame(FT)) stop("`FT` must be a data frame")
  if (!is.numeric(intensity)) stop("`intensity` must be a numeric vector")
  if (!is.numeric(mz_threshold) || mz_threshold < 0) stop("`mz_threshold` must be >= 0")

  # Validate polarity
  polarity <- validate_polarity(polarity)

  new("PeakAnnotationList",
      compounds_list = compounds_list,
      FT = FT,
      adducts_df = adducts_df,
      neutral_loss_table = neutral_loss_table,
      polarity = polarity,
      intensity = intensity,
      mz_threshold = mz_threshold,
      PI_rt_threshold = PI_rt_threshold,
      rt_threshold = rt_threshold,
      iso_threshold = iso_threshold,
      max_charge = max_charge,
      max_isotope = max_isotope)

}
