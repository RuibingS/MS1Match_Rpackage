#' Create a PeakAnnotationList object
#'
#' @description Constructor for the `PeakAnnotationList` S4 class.
#' Automatically determines if RT-constrained annotation is possible
#' based on the presence of the `rtmed` column in the feature table.
#'
#' @param compounds A named list of parsed compounds
#' @param FT A data.frame of feature data (must include `mzmed`, optionally `rtmed`)
#' @param adducts_df Adduct table
#' @param neutral_loss_table Neutral loss table
#' @param polarity "pos" or "neg"
#' @param intensity Numeric vector of intensities
#' @param ... Additional arguments passed to the `PeakAnnotationList` constructor
#'
#' @return An S4 object of class `PeakAnnotationList`
#' @export
createPeakAnnotation <- function(compounds, FT, adducts_df, neutral_loss_table,
                                 polarity, intensity, ...) {
  if (!is.list(compounds)) {
    stop("`compounds` must be a list.")
  }

  PeakAnnotationList(
    compounds_list = compounds,
    FT = FT,
    adducts_df = adducts_df,
    neutral_loss_table = neutral_loss_table,
    polarity = polarity,
    intensity = intensity,
    ...
  )
}
