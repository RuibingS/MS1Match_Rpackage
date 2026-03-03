#' Import a CSU–PMF Mass Spectral Library (custom MSP format)
#'
#' @description
#' Imports and parses a CSU–PMF `.mspLib` file into R, extracting compound
#' metadata and peak tables. This version supports multiple `SYNON` entries
#' and long one-line peak lists.
#'
#' @param file_path Path to the `.mspLib` file.
#'
#' @return A list of compounds. Each compound is a named list with:
#'   \itemize{
#'     \item Metadata fields (Name, Formula, ExactMass, SYNON, etc.)
#'     \item \code{Peaks}: a data.frame with numeric columns:
#'       \itemize{
#'         \item \code{mz}: m/z values
#'         \item \code{rel_intensity}: relative intensity values
#'       }
#'   }
#'
#' @details
#' - Compatible with CSU–PMF library format (Broeckling et al. 2016).
#' - Parses `SYNON:` lines as a character vector.
#' - Supports peak lists written on a single line or multiple lines.
#' - Calculates \code{ExactMass} if missing but a formula is provided.
#'
#' @examples
#' \dontrun{
#' csu_lib <- import_msp_csu("CSUPMF_C8_lib.mspLib")
#' length(csu_lib)
#' str(csu_lib[[1]])
#' }
#'
#' @import MetaboCoreUtils
#' @export
import_msp_csu <- function(file_path) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)

  # call Rcpp parser
  parsed <- parse_library_CSU_cpp(file_path)

  # post-processing
  parsed <- lapply(parsed, function(compound) {
    # convert ExactMass to numeric or calculate from formula
    if (!is.null(compound$ExactMass)) {
      compound$ExactMass <- suppressWarnings(as.numeric(compound$ExactMass))
    } else if (!is.null(compound$Formula)) {
      compound$ExactMass <- tryCatch(
        MetaboCoreUtils::calculateMass(compound$Formula)[[1]],
        error = function(e) NA_real_
      )
    } else {
      compound$ExactMass <- NA_real_
    }

    # ensure Peaks is a data.frame with numeric columns
    if (!is.null(compound$Peaks)) {
      compound$Peaks$mz <- as.numeric(compound$Peaks$mz)
      compound$Peaks$rel_intensity <- as.numeric(compound$Peaks$rel_intensity)
    }

    compound
  })

  parsed
}
