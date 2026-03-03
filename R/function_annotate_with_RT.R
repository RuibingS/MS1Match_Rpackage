#'  Internal function RT Constrained Annotation Function
#'
#'@description Matches features to compound list using m/z and RT constraints.
#' @param object A `PeakAnnotation` object
#' @return Annotated data.frame or updates the object
#' @export
annotate_has_rt <- function(object) {
  # Call Rcpp backend
  annotate_with_rt_cpp(
    compounds_list = object@compounds_list,
    FT = object@FT,
    adducts_df = object@adducts_df,
    neutral_loss_table = object@neutral_loss_table,
    polarity = object@polarity,
    mz_threshold = object@mz_threshold,
    PI_rt_threshold = object@PI_rt_threshold,
    rt_threshold = object@rt_threshold,
    iso_threshold = object@iso_threshold,
    max_charge = object@max_charge,
    max_isotope = object@max_isotope
  )
}

