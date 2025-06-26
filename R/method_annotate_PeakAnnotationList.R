#' @export
setMethod("annotate", "PeakAnnotationList", function(object) {
  rt_cols <- c("rtmed", "rt", "retention_time")
  available_rt_col <- intersect(rt_cols, colnames(object@FT))

  has_rt <- length(available_rt_col) > 0 && is.numeric(object@FT[[available_rt_col[1]]])

  if (has_rt) {
    message("Using RT-constrained annotation.")
    annotate_has_rt(object)
  } else {
    message("No RT info found. Running unconstrained annotation.")
    annotate_no_rt(object)
  }
})
