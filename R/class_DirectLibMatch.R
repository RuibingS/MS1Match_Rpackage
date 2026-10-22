#' Class for Direct Library Matching
#'
#' @description
#' Represents a configuration for matching detected features against a compound library,
#' using m/z and RT thresholds, ppm tolerance, MS level, and ion polarity.
#'
#' @slot FT A \code{data.frame} of detected features (must include m/z and RT columns)
#' @slot compounds_list A \code{list} of known compounds with metadata
#' @slot rt_thr_compound Numeric loose RT threshold for precursor gate vs library RT (e.g. 30 s)
#' @slot rt_thr_coelution Numeric strict RT threshold for co-elution features to
#'   the matched precursor feature (e.g. 5 s).
#' @slot ppm Numeric parts-per-million tolerance for m/z
#' @slot ms_level Character, either \code{"MS1"} or \code{"MS2"}
#' @slot polarity Character, either \code{"pos"} or \code{"neg"}
#' @slot results List of matched results
#'
#' @export
setClass("DirectLibMatch",
         slots = c(
           FT             = "data.frame",
           compounds_list = "list",
           rt_thr_compound    = "numeric",
           rt_thr_coelution = "numeric",
           ppm            = "numeric",
           ms_level       = "character",
           polarity       = "character",
           results        = "list"
         ))
