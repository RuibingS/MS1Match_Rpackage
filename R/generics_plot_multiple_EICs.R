#' Plot multiple EICs from an EICData object
#' @param object An \link{EICData} object.
#' @param feature_indices Integer vector: indices of features to plot.
#' @return A ggplot object showing overlaid EICs.
#' @export
setGeneric("plot_multiple_EICs", function(object, feature_indices) {
  standardGeneric("plot_multiple_EICs")
})

#' Plot the EIC for one feature across selected samples
#'
#' @description
#' For a single feature index `i`, extract the chromatogram from `object@eic_data[i, j]`
#' for each sample index in `j_vec` (or all samples if `j_vec` is `NULL`), and plot
#' the EIC traces with ggplot2.
#'
#' @param object An `EICData` S4 object containing an `@eic_data` matrix-like of
#'   `XChromatogram` objects (features in rows, samples in columns).
#' @param i Integer scalar, row (feature) index to plot.
#' @param j_vec Integer vector of sample indices (columns) to include. If `NULL`,
#'   all columns are used.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' # plot feature 1 for all samples
#' plot_single_feature_EIC(my_eic, i = 1)
#'
#' # plot feature 5 for samples 1:10
#' plot_single_feature_EIC(my_eic, i = 5, j_vec = 1:10)
#' }
#'
#' @export
setGeneric("plot_single_feature_EIC",
           function(object, i, j_vec = NULL) standardGeneric("plot_single_feature_EIC"))
