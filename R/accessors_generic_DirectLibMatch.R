#' Access results from DirectLibMatch object
#'
#' @param object A \code{DirectLibMatch} object
#' @return List of annotated matches
#' @export
setGeneric("results", function(object) standardGeneric("results"))

#' Access features data
#'
#' @param object A \code{DirectLibMatch} object
#' @return Data.frame of detected features
#' @export
setGeneric("features", function(object) standardGeneric("features"))

#' Access library compounds
#'
#' @param object A \code{DirectLibMatch} object
#' @return List of library compounds
#' @export
setGeneric("compounds", function(object) standardGeneric("compounds"))

#' Access strict RT co-elution threshold
#'
#' @param object A \code{DirectLibMatch} object
#' @return Numeric strict RT threshold used for fragment / isotope candidate pool
#' @export
setGeneric("rt_thr_coelution", function(object) standardGeneric("rt_thr_coelution"))

#' Annotate features using library matching and cosine similarity
#'
#' @param object A \code{DirectLibMatch} object
#' @param compute_cosine Logical. Compute cosine similarity scores? Default: TRUE
#' @param ppm Parts-per-million tolerance for cosine scoring
#' @param BPPARAM A \code{BiocParallelParam} object
#' @return Updated \code{DirectLibMatch} object
#' @export
setGeneric(
  "libAnnotate",
  function(object, compute_cosine = TRUE, ppm = 10, BPPARAM = BiocParallel::bpparam())
    standardGeneric("libAnnotate")
)
