#' Parse a Mass Spectral Library
#'
#' @param file_path Path to the library file
#' @param spectrum_type MS1 or MS2 spectrum type (optional)
#' @param ionPolarity Positive or negative ion mode (optional)
#' @param ioniMethod Ionization method (optional)
#' @param colEnergy Collision energy (optional)
#' @param instName Instrument name (optional)
#'
#' @return A parsed library object (list of compounds with metadata and peaks)
#' @import MetaboCoreUtils
#' @export
parse_library <- function(file_path,
                          spectrum_type = NULL,
                          ionPolarity = NULL,
                          ioniMethod = NULL,
                          colEnergy = NULL,
                          instName = NULL) {

  if (!file.exists(file_path)) stop("File not found: ", file_path)

  # Parse the library with Rcpp
  parsed <- parse_library_cpp(
    file_path,
    spectrum_type,
    ionPolarity,
    ioniMethod,
    colEnergy,
    instName
  )

  # cal ExactMass if Formula is present
  parsed <- lapply(parsed, function(compound) {
    if (!is.null(compound$Formula)) {
      exact_mass <- tryCatch(
        MetaboCoreUtils::calculateMass(compound$Formula)[[1]],
        error = function(e) {
          warning("Failed to calculate exact mass for formula: ", compound$Formula)
          NA_real_
        }
      )
      compound$ExactMass <- exact_mass
    } else {
      compound$ExactMass <- NA_real_
    }
    compound
  })

  parsed
}
