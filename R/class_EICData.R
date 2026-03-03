#' EICData S4 class
#' @description Holds the XChromatograms object and corresponding feature table.
#' @slot eic_data An XChromatograms object containing EICs.
#' @slot features A data.frame with feature info (e.g. feature_name, rtmed).
#' @slot parsed_EIC A list to store the output from `feature_EIC_process_fun()`
#' @slot correlation_method Character string: correlation method used
#' @slot correlation_matrix A list to store pair wise Eic correlation score
#' @export
setClass("EICData",
         slots = list(
           features = "data.frame",
           eic_data = "XChromatograms",
           parsed_EIC = "list",
           correlation_method = "character",
           correlation_df = "list"
         ))
