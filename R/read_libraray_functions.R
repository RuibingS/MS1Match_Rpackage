# =============================================================================
# read_library.R
# Step 1 of the conversion file workflow.
# Generic permissive reader for any MSP-like file format.
# =============================================================================
#' Read any MSP-like spectral library file
#'
#' @description
#' Generic reader for spectral library files regardless of field naming
#' conventions. Accepts Bruker \code{.library}, CSU-PMF \code{.mspLib},
#' MoNA \code{.msp}, NIST \code{.msp}, or any MSP-like plain text format.
#'
#' This function is permissive — it reads whatever fields are present without
#' normalizing or filtering. Use it as \strong{Step 1} of the workflow before
#' \code{\link{convert_to_nist_msp}}.
#'
#' For importing already-clean NIST MSP files into your analysis pipeline,
#' use \code{\link{import_msp}} directly instead.
#'
#' @param file_path Path to the library file
#'   (\code{.msp}, \code{.library}, \code{.mspLib}, \code{.txt}).
#' @param dedup_by Deduplication strategy. Default \code{"none"} — keep all
#'   entries so nothing is lost before conversion. Options:
#'   \code{"name_adduct_mode_energy"}, \code{"name_adduct_mode"},
#'   \code{"name_only"}, \code{"none"}.
#' @param reserve_n Memory hint for expected number of compounds (default 10000).
#' @param max_records Integer. Stop after this many records (-1 = no limit).
#'
#' @return A named list of compounds, each containing raw metadata fields and a
#'   \code{Peaks} data frame (\code{mz}, \code{rel_intensity}).
#'
#' @examples
#' \dontrun{
#' # Step 1a: Read Bruker .library (non-standard fields: IonPolarity, InstName, etc.)
#' bruker_raw <- read_library("inhouse_library.library")
#'
#' # Step 1b: Read CSU-PMF .mspLib (SYNON-encoded metadata)
#' csu_raw <- read_library("CSUPMF_C8_lib.mspLib")
#'
#' # Step 2: Convert both to clean NIST MSP (run once, then use the .msp files)
#' convert_to_nist_msp(bruker_raw, "inhouse_library.msp",
#'                     source_lab = "Your Lab, Your University")
#' convert_to_nist_msp(csu_raw, "CSUPMF_C8.msp",
#'                     source_lab = "Broeckling Lab, Colorado State University")
#' }
#'
#' @seealso \code{\link{convert_to_nist_msp}}, \code{\link{import_msp}}
#' @export
read_library <- function(file_path,
                         dedup_by    = "none",
                         reserve_n   = 10000L,
                         max_records = -1L) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  ext <- tolower(tools::file_ext(file_path))
  if (!ext %in% c("msp", "msplit", "library", "txt", "msplib")) {
    warning("Unexpected file extension '.", ext,
            "' — attempting to parse as MSP format.")
  }
  parsed <- parse_msp_library_cpp(
    file_path,
    spec     = NULL,
    pol      = NULL,
    meth     = NULL,
    cen      = NULL,
    inst     = NULL,
    dedup    = dedup_by,
    reserven = as.integer(reserve_n),
    maxrec   = as.integer(max_records)
  )
  if (length(parsed) == 0) {
    warning("No compounds returned from: ", file_path)
    return(list())
  }
  nms <- vapply(parsed, .safe_extract_name, character(1))
  names(parsed) <- make.unique(nms, sep = ".")
  message("[read_library] Read ", length(parsed),
          " entries from ", basename(file_path))
  parsed
}
