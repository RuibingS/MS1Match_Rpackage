#' @rdname results
setMethod("results", "DirectLibMatch", function(object) object@results)

#' @rdname features
setMethod("features", "DirectLibMatch", function(object) object@FT)

#' @rdname compounds
setMethod("compounds", "DirectLibMatch", function(object) object@compounds_list)

#' @rdname rt_thr_coelution
setMethod("rt_thr_coelution", "DirectLibMatch", function(object) object@rt_thr_coelution)
