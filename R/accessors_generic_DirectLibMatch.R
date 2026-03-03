#' Access results from DirectLibMatch object
#'
#' @param object A \code{DirectLibMatch} object
#' @return List of annotated matches
#' @export
setGeneric("results", function(object) standardGeneric("results"))

#' Access features data
#' @param object A \code{DirectLibMatch} object
#' @return Data.frame of detected features
#' @export
setGeneric("features", function(object) standardGeneric("features"))

#' Access library compounds
#' @param object A \code{DirectLibMatch} object
#' @return List of library compounds
#' @export
setGeneric("compounds", function(object) standardGeneric("compounds"))
