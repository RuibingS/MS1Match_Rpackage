#' Import a Mass Spectral Library (MSP format)
#'
#' @description
#'  This function imports and parses an MSP file into R, extracting compound
#' metadata and peak tables.
#'
#' @param file_path Path to the `.msp` file.
#' @param spectrum_type Optional filter for spectrum type: `"MS1"` or `"MS2"`.
#' @param ionPolarity Optional filter for ion mode (e.g. `"POSITIVE"`, `"NEGATIVE"`).
#' @param ioniMethod Optional filter for ionization method (e.g. `"ESI-QTOF"`).
#' @param colEnergy Optional filter for collision energy (string match).
#' @param instName Optional filter for instrument name.
#'
#' @return A list of compounds. Each compound is a named list with:
#'   \itemize{
#'     \item Metadata fields (Name, Formula, ExactMass, InChI, SMILES, etc.)
#'     \item \code{Peaks}: a data.frame with two numeric columns:
#'       \itemize{
#'         \item \code{mz}: m/z values
#'         \item \code{rel_intensity}: relative intensity values
#'       }
#'   }
#'
#' @details
#' - Peaks are parsed directly from the lines following `"Num Peaks:"`.
#' - Filters (`spectrum_type`, `ionPolarity`, etc.) are applied during parsing.
#' - If \code{ExactMass} is not present but \code{Formula} is, the exact mass
#'   is calculated using \pkg{MetaboCoreUtils}.
#'
#' @examples
#' \dontrun{
#' msp_data <- import_msp("example.msp", spectrum_type = "MS2", ionPolarity = "POSITIVE")
#' str(msp_data[[1]])
#' }
#'
#' @import MetaboCoreUtils
#' @export
import_msp <- function(file_path,
                       spectrum_type = NULL,
                       ionPolarity = NULL,
                       ioniMethod = NULL,
                       colEnergy = NULL,
                       instName = NULL) {

  if (!file.exists(file_path)) stop("File not found: ", file_path)

  parsed <- parse_msp_library_cpp(
    file_path,
    spectrum_type,
    ionPolarity,
    ioniMethod,
    colEnergy,
    instName
  )

  parsed <- lapply(parsed, function(compound) {
    # Case 1: ExactMass provided in MSP → keep as numeric
    if (!is.null(compound$ExactMass)) {
      compound$ExactMass <- suppressWarnings(as.numeric(compound$ExactMass))

      # Case 2: No ExactMass but Formula available → calculate
    } else if (!is.null(compound$Formula)) {
      compound$ExactMass <- tryCatch(
        MetaboCoreUtils::calculateMass(compound$Formula)[[1]],
        error = function(e) NA_real_
      )

      # Case 3: Neither present
    } else {
      compound$ExactMass <- NA_real_
    }
    compound
  })

  parsed
}
