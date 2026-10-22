#' Perform primary-ion matching against a direct target list
#'
#' @description
#' Matches feature table m/z values against a compound target list using
#' \code{[M+H]+} (positive mode) or \code{[M-H]-} (negative mode) ion m/z.
#'
#' The target list is automatically normalised via \code{normalize_target_list()}
#' before matching, so input data frames from different sources (HMDB, plate-based,
#' PA14, etc.) are accepted without manual column renaming.
#'
#' @param object A \code{DirectTargetMatch} object
#' @return The same \code{DirectTargetMatch} object with \code{@results} populated
#' @export
setGeneric("PImatch_target", function(object) standardGeneric("PImatch_target"))

#' @export
setMethod("PImatch_target", "DirectTargetMatch", function(object) {

  # Normalise target list to unified schema regardless of source format
  comp_norm <- normalize_target_list(object@Comp_data)

  polarity <- if (object@PIon == "[M+H]+") "pos" else "neg"

  object@results <- PImatch_target_rcpp(
    FT              = object@FT,
    Comp_data       = comp_norm,
    polarity        = polarity,
    diff_mz_thr     = object@diff_mz_thr,
    ppm             = object@ppm,
    rt_thr_compound = object@rt_thr_compound
  )

  object
})
