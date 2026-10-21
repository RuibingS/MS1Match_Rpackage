#' Annotate features using library matching and cosine similarity
#'
#' @description
#' This function matches features to compounds using direct m/z and RT comparison,
#' followed optionally by cosine similarity scoring. Cosine scoring is computed in parallel
#' using \pkg{BiocParallel}.
#'
#' Two RT thresholds are used:
#' \describe{
#'   \item{\code{rt_thr_compound}}{Loose threshold for the precursor gate: how far a feature's
#'     RT may deviate from the library compound's RT (e.g. 30 s).}
#'   \item{\code{rt_thr_coelution}}{Strict threshold for the co-elution pool: only features
#'     within this window of the matched precursor feature are considered as fragment /
#'     isotope candidates (e.g. 5 s). Applies to both MS1 and MS2 matching.}
#' }
#'
#' @param object A \code{DirectLibMatch} object containing the feature table and compound list
#' @param compute_cosine Logical. Should cosine similarity scores be computed? Default is \code{TRUE}
#' @param ppm Numeric. Parts-per-million (ppm) tolerance used in cosine similarity scoring. Default is \code{10}
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallelization.
#'   Default is \code{BiocParallel::bpparam()}.
#'
#' @return An updated \code{DirectLibMatch} object with added \code{@results} slot
#'   containing cosine similarity and matched features
#'
#' @import BiocParallel
#' @import Spectra
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @export
#' @rdname libAnnotate
setMethod(
  "libAnnotate", "DirectLibMatch",
  function(object, compute_cosine = TRUE, ppm = 10, BPPARAM = BiocParallel::bpparam()) {

    safe_max <- function(x) {
      if (length(x) == 0 || all(is.na(x))) return(1)
      m <- suppressWarnings(max(x, na.rm = TRUE))
      if (!is.finite(m) || m == 0) return(1)
      m
    }

    get_vector <- function(primary, fallback) {
      if (!is.null(primary) && length(primary) > 0) return(primary)
      fallback
    }

    # ---- core annotation (single-threaded, Rcpp) ----
    # Both MS1 and MS2 now receive the strict RT threshold for their
    # co-elution / fragment candidate pool.
    if (object@ms_level == "MS1") {
      annotation_results <- direct_match_MS1Lib_cpp(
        object@FT,
        object@compounds_list,
        object@rt_thr_compound,
        object@rt_thr_coelution,
        object@ppm
      )
    } else if (object@ms_level == "MS2") {
      annotation_results <- direct_match_MS2Lib_cpp(
        object@FT,
        object@compounds_list,
        object@rt_thr_compound,
        object@rt_thr_coelution,
        object@ppm
      )
    } else {
      stop("Unknown MS level in object@ms_level: must be 'MS1' or 'MS2'")
    }

    if (isTRUE(compute_cosine)) {

      if (!methods::is(BPPARAM, "BiocParallelParam")) {
        stop("BPPARAM must be a BiocParallelParam, e.g. SerialParam(), SnowParam(), MulticoreParam().")
      }

      annotation_results <- BiocParallel::bplapply(
        X = seq_along(annotation_results),
        FUN = function(i) {
          match <- annotation_results[[i]]

          if (all(is.na(match$Compound_peak_mz)) || all(is.na(match$Matched_features_mz))) {
            match$inner_cosine    <- NA_real_
            match$outer_cosine    <- NA_real_
            match$intensity_ratio <- NA_real_
            return(match)
          }

          ref_peaks <- data.frame(
            mz            = get_vector(match$Compound_peak_mz_full, match$Compound_peak_mz),
            rel_intensity = get_vector(match$Compound_peak_int_full, match$Compound_peak_int)
          )
          ref_peaks$rel_intensity <- 1000 * ref_peaks$rel_intensity / safe_max(ref_peaks$rel_intensity)

          query_df <- data.frame(
            mzmed         = match$Matched_features_mz,
            rel_intensity = match$Matched_features_int
          )
          query_df$rel_intensity <- 1000 * query_df$rel_intensity / safe_max(query_df$rel_intensity)

          if (nrow(query_df) == 0 || nrow(ref_peaks) == 0) {
            match$inner_cosine    <- NA_real_
            match$outer_cosine    <- NA_real_
            match$intensity_ratio <- NA_real_
            return(match)
          }

          spd <- S4Vectors::DataFrame(
            msLevel  = rep(1L, 2),
            polarity = rep(1L, 2),
            mz       = IRanges::NumericList(query_df$mzmed,         ref_peaks$mz),
            intensity= IRanges::NumericList(query_df$rel_intensity, ref_peaks$rel_intensity)
          )
          sps <- Spectra::Spectra(spd)

          match$inner_cosine <- tryCatch(
            Spectra::compareSpectra(sps, type = "inner", ppm = ppm)[1, 2],
            error = function(e) NA_real_)
          match$outer_cosine <- tryCatch(
            Spectra::compareSpectra(sps, type = "outer", ppm = ppm)[1, 2],
            error = function(e) NA_real_)

          max_query_idx <- which.max(match$Matched_features_int)
          max_ref_idx   <- which.max(match$Compound_peak_int)

          mz_query <- match$Matched_features_mz[max_query_idx]
          mz_ref   <- match$Compound_peak_mz[max_ref_idx]

          if (!is.na(mz_query) && !is.na(mz_ref) && abs(mz_query - mz_ref) <= 0.01) {
            scaled_query          <- 999 * match$Matched_features_int / safe_max(match$Matched_features_int)
            scaled_ref            <- match$Compound_peak_int
            match$intensity_ratio <- scaled_query / scaled_ref
          } else {
            match$intensity_ratio <- rep(NA_real_, length(match$Matched_features_int))
          }

          match
        },
        BPPARAM = BPPARAM
      )
    }

    object@results <- annotation_results
    object
  }
)
