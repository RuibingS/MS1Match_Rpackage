#' Flatten matched PI result list into a data.frame
#'
#' @description
#' Filters out empty (unmatched) feature slots and row-binds the remainder
#' into a single flat \code{data.frame}. The result includes the new columns
#' added by the updated Rcpp function: \code{exactMass}, \code{Comp_mz},
#' \code{Comp_rt}, \code{mz_diff}, and \code{ppm_error}.
#'
#' @param object A \code{DirectTargetMatch} object, or a raw \code{list}
#'   as returned by \code{PImatch_target_rcpp()}
#' @return A \code{data.frame} of all matched feature–compound pairs,
#'   or an empty \code{data.frame} if nothing matched
#' @export
setGeneric("flatten", function(object) standardGeneric("flatten"))

#' @export
setMethod("flatten", "DirectTargetMatch", function(object) {
  .flatten_results(object@results)
})

#' @export
setMethod("flatten", "list", function(object) {
  .flatten_results(object)
})

# Internal helper shared by both methods
.flatten_results <- function(res) {
  filtered <- Filter(function(x) length(x$Feature_name) > 0, res)
  if (length(filtered) == 0) return(data.frame())
  do.call(rbind, lapply(filtered, as.data.frame, stringsAsFactors = FALSE))
}
