#' Direct match method for DirectLibMatch object
#'
#' Matches features against a compound list using direct m/z and RT thresholding,
#' considering ppm tolerance, MS level, and ion polarity.
#'
#' @param object A \code{DirectLibMatch} object containing features, compound library,
#'        matching thresholds, MS level, and polarity.
#'
#' @return A \code{list} of matched compounds.
#' @param object A \code{DirectLibMatch} object containing features, compound library,
#'        matching thresholds, MS level, and polarity.
#'
#' @return A \code{list} of matched compounds.
#' @export
setMethod("libAnnotate", "DirectLibMatch", function(object) {
  FT <- object@FT
  lib <- object@compounds_list
  ion_mode <- object@polarity
  ms_level <- object@ms_level
  diff_mz <- object@diff_mz_thr
  diff_rt <- object@diff_rt_thr
  ppm <- object@ppm

  if (!is.data.frame(FT)) stop("Feature table must be a data.frame.")
  if (!is.list(lib)) stop("Library must be a list.")

  ms_level <- match.arg(ms_level, choices = c("MS1", "MS2"))

  if (ms_level == "MS2") {
    missing_preion <- sapply(lib, function(comp) {
      is.null(comp[["PreIon"]])
    })

    if (all(missing_preion)) {
      stop("All compounds in the MS2 library are missing 'PreIon'. This field is required for MS2 matching.")
    }

    if (any(missing_preion)) {
      warning(sprintf("%d of %d compounds are missing 'PreIon' and will be skipped.",
                      sum(missing_preion), length(lib)))
    }
  }

  if (ms_level == "MS1") {
    message("Running MS1 library matching...")
    result <- direct_match_MS1Lib_cpp(
      FT = FT,
      compounds_list = lib,
      diff_mz_thr = diff_mz,
      diff_rt_thr = diff_rt,
      ppm = ppm
    )
  } else if (ms_level == "MS2") {
    message("Running MS2 library matching...")
    result <- direct_match_MS2Lib_cpp(
      FT = FT,
      compounds_list = lib,
      diff_mz_thr = diff_mz,
      diff_rt_thr = diff_rt,
      ppm = ppm
    )
  } else {
    stop("Unsupported ms_level: must be 'MS1' or 'MS2'.")
  }

  return(result)
})

