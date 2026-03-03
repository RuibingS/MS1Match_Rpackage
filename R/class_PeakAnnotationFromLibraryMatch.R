#' PeakAnnotationFromLibraryMatch class
#'
#' @description
#' Container for feature-level ion annotation derived from library matching results.
#'
#' @slot library_matches list. Output of MS1/MS2 library matching.
#' @slot adducts_df data.frame. Adduct definition table.
#' @slot neutral_loss_table data.frame. Neutral loss table.
#' @slot polarity character. "positive" or "negative".
#' @slot mz_threshold numeric. Absolute m/z tolerance (Da).
#' @slot ppm numeric. ppm tolerance.
#' @slot rt_threshold numeric. RT tolerance (seconds).
#' @slot iso_threshold numeric. Isotope mass delta.
#' @slot max_charge integer. Maximum charge.
#' @slot max_isotope integer. Maximum isotope number.
#'
#' @name PeakAnnotationFromLibraryMatch-class
#' @rdname PeakAnnotationFromLibraryMatch-class
#' @exportClass PeakAnnotationFromLibraryMatch
setClass(
  "PeakAnnotationFromLibraryMatch",
  slots = c(
    library_matches = "list",
    adducts_df = "data.frame",
    neutral_loss_table = "data.frame",
    polarity = "character",
    mz_threshold = "numeric",
    ppm = "numeric",
    rt_threshold = "numeric",
    iso_threshold = "numeric",
    max_charge = "integer",
    max_isotope = "integer"
  )
)
