#' Annotate peaks for compound list input
#' @export
setMethod("annotate", "PeakAnnotationList", function(object) {
  has_rt <- any(c("rtmed", "rt", "retention_time") %in% colnames(object@FT)) &&
    is.numeric(object@FT[[intersect(c("rtmed", "rt", "retention_time"), colnames(object@FT))[1]]])

  if (has_rt) {
    message("Using RT-constrained annotation.")
    annotate_has_rt(object)
  } else {
    message("No RT info. Running unconstrained annotation.")
    annotate_no_rt(object)
  }
})


