#' Base class for peak annotation
#'
#' @slot FT Data frame with feature table.
#' @slot compounds_list List of processed compounds.
#' @slot adducts_df Data frame of adduct definitions.
#' @slot neutral_loss_table Data frame of neutral losses.
#' @slot polarity Character. Either "positive" or "negative".
#' @slot mz_threshold Numeric. m/z tolerance for matching.
#' @slot PI_rt_threshold Numeric. RT threshold for PI detection.
#' @slot rt_threshold Numeric. RT threshold for general matching.
#' @slot iso_threshold Numeric. Isotopic mass delta.
#' @slot max_charge Integer. Maximum charge state.
#' @slot max_isotope Integer. Maximum isotope peaks.
#' @export
setClass("PeakAnnotation", slots = c(
  FT = "data.frame",
  compounds_list = "list",
  adducts_df = "data.frame",
  neutral_loss_table = "data.frame",
  polarity = "character",
  mz_threshold = "numeric",
  PI_rt_threshold = "numeric",
  rt_threshold = "numeric",
  iso_threshold = "numeric",
  max_charge = "integer",
  max_isotope = "integer"
))
