#' @rdname results
setMethod("results", "DirectLibMatch", function(object) object@results)

#' @rdname features
setMethod("features", "DirectLibMatch", function(object) object@FT)

#' @rdname compounds
setMethod("compounds", "DirectLibMatch", function(object) object@compounds_list)
