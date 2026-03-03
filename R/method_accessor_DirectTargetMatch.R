#' Flatten matched PI result list into data.frame
#'
#' @param object List result from matchPI()
#'
#' @return A data.frame
#' @export
setGeneric("flatten", function(object) standardGeneric("flatten"))

#' @export
setMethod("flatten", "list", function(object) {
  filtered <- Filter(function(x) length(x$Feature_name) > 0, object)
  if (length(filtered) == 0) {
    return(data.frame())
  }
  do.call(rbind, lapply(filtered, as.data.frame))
})
