#' LC peak annotation (R wrapper for rcpp function)
#'
#' @param object A PeakAnnotationList or PeakAnnotationDF object
#' @return A list of annotations
#' @keywords internal
annotate_LC <- function(object) {
  tryCatch({
    if (inherits(object, "PeakAnnotationList")) {
      annotate_LC_cpp(
        compounds_list = object@compounds_list,
        FT = object@FT,
        adducts_df = object@adducts_df,
        neutral_loss_table = object@neutral_loss_table,
        polarity = object@polarity,
        intensity = object@intensity
      )
    } else if (inherits(object, "PeakAnnotationDF")) {
      annotate_LC_cpp(
        compounds_df = object@compounds_df,
        FT = object@FT,
        adducts_df = object@adducts_df,
        neutral_loss_table = object@neutral_loss_table,
        polarity = object@polarity,
        intensity = object@intensity
      )
    } else {
      stop("Unsupported object class for LC annotation.")
    }
  }, error = function(e) {
    stop("LC annotation failed: ", e$message)
  })
}
