#' @export
setGeneric("rt_similarity", function(object) standardGeneric("rt_similarity"))
setMethod("rt_similarity", "RTSimilarity", function(object) object@sim_matrix)

#' @export
setGeneric("rt_distance", function(object) standardGeneric("rt_distance"))
setMethod("rt_distance", "RTSimilarity", function(object) object@dist_matrix)

#' @export
setGeneric("rt_method", function(object) standardGeneric("rt_method"))
setMethod("rt_method", "RTSimilarity", function(object) object@method)
