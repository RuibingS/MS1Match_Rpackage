#' @describeIn annotate
#' Annotate features from precomputed spectral library match results.
#' Adduct, neutral-loss, and isotope annotation are performed within each match group.
#' @export
setMethod("annotate", "PeakAnnotationFromLibraryMatch", function(object, ...) {

  annotate_from_librarymatch_cpp(
    library_matches = object@library_matches,
    adducts_df = object@adducts_df,
    neutral_loss_table = object@neutral_loss_table,
    polarity = object@polarity,
    mz_threshold = object@mz_threshold,
    ppm = object@ppm,
    rt_threshold = object@rt_threshold,
    iso_threshold = object@iso_threshold,
    max_charge = object@max_charge,
    max_isotope = object@max_isotope
  )

})
