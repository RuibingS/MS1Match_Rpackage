#' Annotate features using library matching and cosine similarity
#'
#' @param object A \code{DirectLibMatch} object
#' @param compute_cosine Logical. Compute cosine similarity scores? Default: TRUE
#' @param ppm Parts-per-million tolerance for cosine scoring
#' @return Updated \code{DirectLibMatch} object
#' @export
setGeneric("libAnnotate", function(object, compute_cosine = TRUE, ppm = 10) standardGeneric("libAnnotate"))
