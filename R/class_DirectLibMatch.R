#' Class for Direct Library Matching
#'
#' Represents a configuration for matching detected features against a compound library,
#' using m/z and RT thresholds, ppm tolerance, MS level, and ion polarity.
#'
#' @slot FT A \code{data.frame} of detected features (must include m/z and RT columns)
#' @slot compounds_list A \code{list} of known compounds with metadata
#' @slot diff_mz_thr Numeric m/z difference threshold
#' @slot diff_rt_thr Numeric RT difference threshold
#' @slot ppm Numeric parts-per-million tolerance for m/z
#' @slot ms_level Character, either \code{"MS1"} or \code{"MS2"}
#' @slot polarity Character, either \code{"pos"} or \code{"neg"}
#'
#' @export
setClass("DirectLibMatch",
         slots = c(
           FT = "data.frame",
           compounds_list = "list",
           diff_mz_thr = "numeric",
           diff_rt_thr = "numeric",
           ppm = "numeric",
           ms_level = "character",   # "MS1" or "MS2"
           polarity = "character"    # "pos" or "neg"
         ))
