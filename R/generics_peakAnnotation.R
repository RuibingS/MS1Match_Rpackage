#' Annotate peaks
#'
#' @description
#' S4 generic for running peak annotation.
#' Methods are provided for:
#' \itemize{
#'   \item [PeakAnnotationFromLibraryMatch-class]: annotation of precomputed library matches
#' }
#'
#' @param object Object to annotate.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return Method-dependent annotation results.
#'
#' @export
setGeneric("annotateLibMatch", function(object, ...) standardGeneric("annotateLibMatch"))
