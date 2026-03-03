#' Annotate peaks
#'
#' @description
#' S4 generic for running peak annotation.
#' Methods are provided for:
#' \itemize{
#'   \item [PeakAnnotation-class]: full feature-table × library search
#'   \item [PeakAnnotationFromLibraryMatch-class]: annotation of precomputed library matches
#' }
#'
#' @param object Object to annotate.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return Method-dependent annotation results.
#'
#' @export
setGeneric("annotate", function(object, ...) standardGeneric("annotate"))
