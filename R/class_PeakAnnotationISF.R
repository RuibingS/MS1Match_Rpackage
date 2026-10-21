# ============================================================
# class_PeakAnnotationISF.R
# S4 class for MS2 fragment-based annotation.
#
# Tier naming convention:
#   1c     — MS2 fragment match, intra-group, Phase 1
#   1c_exp — MS2 fragment match, RT-gated pool, Phase 2
#   1d     — Neutral loss (pairwise), Phase 1
#   1d_exp — Neutral loss (pairwise), Phase 2
#
# Priority: 1a > 1b > 1c > 1d > 2
# ============================================================

#' PeakAnnotationISF S4 class
#'
#' @description
#' Container for MS2 fragment-based annotation of an MS1 feature table.
#'
#' \strong{Phase 1 — Intra-group MS2 fragment (Tier 1c)}
#' For each confirmed parent feature (identified by \code{parent_col} ==
#' \code{parent_val}) with a consensus MS2 spectrum, every other feature
#' in the same library group is checked: if its \code{mzmed} matches any
#' fragment ion within \code{ppm_isf} and \code{rt_isf_threshold} seconds,
#' it receives Tier 1c.
#'
#' \strong{Phase 2 — RT-gated MS2 fragment (Tier 1c_exp)}
#' The full feature pool (features not in any Phase 1 group) is searched
#' within \code{rt_isf_threshold} seconds of the parent.  Matching
#' features receive Tier 1c_exp.
#'
#' @slot phase1_result data.frame. Output of \code{annotate()}.
#' @slot ms2_list list. Named list from \code{\link{read_ms2_consensus}()}.
#' @slot FT data.frame. Full MS1 feature pool.
#' @slot parent_col character. Column identifying parent features.
#'   Default \code{"PI_match"}.
#' @slot parent_val ANY. Value marking a parent.  \code{NULL} = any
#'   non-empty string.  Default \code{NULL}.
#' @slot ppm_isf numeric. ppm tolerance for fragment m/z matching.
#'   Default \code{10}.
#' @slot mz_threshold_isf numeric. Absolute m/z floor (Da).
#'   Default \code{0.005}.
#' @slot rt_isf_threshold numeric. RT window (seconds).  Default \code{5}.
#'
#' @name PeakAnnotationISF-class
#' @rdname PeakAnnotationISF-class
#' @exportClass PeakAnnotationISF
setClass(
  "PeakAnnotationISF",
  slots = c(
    phase1_result    = "data.frame",
    ms2_list         = "list",
    FT          = "data.frame",
    parent_col       = "character",
    parent_val       = "ANY",
    ppm_isf          = "numeric",
    mz_threshold_isf = "numeric",
    rt_isf_threshold = "numeric"
  )
)

#' Constructor for PeakAnnotationISF
#'
#' @inheritParams createPeakAnnotationISF
#' @return A \code{PeakAnnotationISF} object.
#' @seealso \code{\link{createPeakAnnotationISF}()}
#' @export
PeakAnnotationISF <- function(
    phase1_result,
    ms2_list,
    FT,
    parent_col       = "PI_match",
    parent_val       = NULL,
    ppm_isf          = 10,
    mz_threshold_isf = 0.005,
    rt_isf_threshold = 5.0
) {
  createPeakAnnotationISF(
    phase1_result    = phase1_result,
    ms2_list         = ms2_list,
    FT          = FT,
    parent_col       = parent_col,
    parent_val       = parent_val,
    ppm_isf          = ppm_isf,
    mz_threshold_isf = mz_threshold_isf,
    rt_isf_threshold = rt_isf_threshold
  )
}
