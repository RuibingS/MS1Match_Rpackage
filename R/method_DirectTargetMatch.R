#' Perform PI Ion Matching on DirectTargetMatch Object
#'
#' @param object A DirectTargetMatch object
#'
#' @return A list of matched features
#' @export
setGeneric("PImatch_target", function(object) standardGeneric("PImatch_target"))

#' @export
setMethod("PImatch_target", "DirectTargetMatch", function(object) {
  PImatch_target_rcpp(
    FT = object@FT,
    Comp_data = object@Comp_data,
    PIon = object@PIon,
    diff_mz_thr = object@diff_mz_thr,
    polarity = ifelse(object@PIon == "[M+H]+", "pos", "neg"),
    ppm = object@ppm
  )
})
