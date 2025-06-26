#' Internal function for non-RT-based annotation
#'
#' @param object A PeakAnnotationList object without RT information
#' @return Annotated peak list (list or data.frame, depending on backend)
annotate_no_rt <- function(object) {
  annotate_without_rt_cpp(
    compounds_list = object@compounds_list,
    FT = object@FT,
    adducts_df = object@adducts_df,
    neutral_loss_table = object@neutral_loss_table,
    polarity = object@polarity,
    intensity = object@intensity,
    mz_threshold = object@mz_threshold,
    iso_threshold = object@iso_threshold,
    max_charge = object@max_charge,
    max_isotope = object@max_isotope
  )
}
