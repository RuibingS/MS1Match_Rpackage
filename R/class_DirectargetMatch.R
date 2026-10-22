#' Class for Direct Target Matching
#'
#' @description
#' Holds all inputs and results for matching a feature table against a
#' compound target list using primary ion m/z (\code{"[M+H]+"} or \code{"[M-H]-"}).
#'
#' The target list (\code{Comp_data}) can be in any supported format —
#' HMDB, plate-based, PA14, etc. Column normalisation is handled
#' automatically inside \code{PImatch_target()}.
#'
#' @slot FT             A \code{data.frame} of detected features (needs
#'   \code{feature_name} and \code{mzmed}/\code{mz} columns; optionally
#'   \code{rtmed}/\code{rt} for RT gating).
#' @slot Comp_data      A \code{data.frame} target list in any supported
#'   format (HMDB, plate-based, PA14, …).
#' @slot PIon           Character. Primary ion selector: \code{"[M+H]+"} for
#'   positive mode, \code{"[M-H]-"} for negative mode.
#' @slot diff_mz_thr    Numeric. Absolute m/z tolerance in Da (floor applied
#'   when ppm tolerance is smaller, e.g. at low m/z).
#' @slot ppm            Numeric. Parts-per-million m/z tolerance.
#' @slot rt_thr_compound Numeric. RT window in seconds for feature vs compound
#'   RT gating. Only applied when the target list contains an RT column and
#'   the feature table has RT values. Set to \code{Inf} to disable.
#' @slot results        List of raw match results (populated by
#'   \code{PImatch_target()}).
#'
#' @export
setClass("DirectTargetMatch",
         slots = c(
           FT              = "data.frame",
           Comp_data       = "data.frame",
           PIon            = "character",
           diff_mz_thr     = "numeric",
           ppm             = "numeric",
           rt_thr_compound = "numeric",
           results         = "list"
         ))
