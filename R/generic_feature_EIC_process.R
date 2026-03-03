#' Process Extracted Ion Chromatograms (EICs) for features
#'
#' @description Parses the EIC data from an \code{EICData} object and prepares it for correlation analysis.
#'
#' @param object An \code{EICData} object containing the raw EIC data and feature table.
#' @param nrows Integer. The number of features (rows) to process from the feature table.
#' @param ncols Integer. The number of samples (columns) in the EIC data matrix.
#'
#' @return An updated \code{EICData} object with the parsed EICs stored in the \code{parsed_EIC} slot.
#' @export
#'
#' @examples
#' obj <- new("EICData", eic_data = result[[4]], features = feature_table)
#' obj <- feature_EIC_process_fun(obj, nrows = 3000, ncols = 9)
setGeneric("feature_EIC_process", function(object, nrows, ncols) {
  standardGeneric("feature_EIC_process")
           }
)
