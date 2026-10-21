
#' Export a compound library list to NIST MSP format
#'
#' @description
#' Writes a named list of compounds (as returned by \code{\link{import_msp}})
#' to a NIST-standard \code{.msp} file. Each compound is written as a
#' self-contained block with metadata fields followed by its peak list.
#'
#' Use this function for:
#' \itemize{
#'   \item Saving a filtered or deduplicated library back to disk
#'   \item Exporting processed libraries for use in other tools (e.g. NIST MS
#'     Search, MZmine, MS-DIAL)
#'   \item Archiving subsets produced by \code{\link{import_msp}}
#' }
#'
#' @param library A named list of compounds as returned by
#'   \code{\link{import_msp}}. Each element must contain a \code{Peaks}
#'   data frame with columns \code{mz} and \code{rel_intensity}, plus
#'   any metadata fields to be written.
#' @param file_path Path to the output \code{.msp} file. Parent directory
#'   must exist. Existing files are overwritten.
#' @param fields Character vector of metadata field names to write for each
#'   compound, in the order they will appear in the output. Fields absent
#'   from a compound are silently skipped. Defaults to all standard NIST MSP
#'   fields:
#'   \code{c("Name", "PrecursorMZ", "Adduct", "Formula", "ExactMass",
#'   "InChIKey", "InChI", "SMILES", "IonMode", "CollisionEnergy",
#'   "InstrumentType", "Instrument", "IonizationType", "Spectrum_type",
#'   "Comments")}.
#' @param mz_digits Integer. Number of decimal places for \code{mz} values
#'   in the peak list. Default \code{4}.
#' @param intensity_digits Integer. Number of decimal places for intensity
#'   values in the peak list. Default \code{0} (integer-like output,
#'   consistent with NIST convention).
#'
#' @return Invisibly returns \code{file_path}. Called primarily for its
#'   side effect of writing the \code{.msp} file.
#'
#' @examples
#' \dontrun{
#' # Round-trip: import, filter, export
#' CBIO_library <- import_msp(
#'   system.file("extdata", "CBIO", "CBIO_Standards.msp", package = "MS1Match"),
#'   ionPolarity   = "pos",
#'   spectrum_type = "MS1",
#'   filter        = list(field = "Comments", pattern = "MSMLS")
#' )
#'
#' export_msp(
#'   library   = CBIO_library,
#'   file_path = "CBIO_MSMLS_filtered.msp"
#' )
#'
#' # Custom field order and precision
#' export_msp(
#'   library          = CBIO_library,
#'   file_path        = "CBIO_MSMLS_filtered.msp",
#'   fields           = c("Name", "Formula", "ExactMass", "Adduct",
#'                        "IonMode", "Comments"),
#'   mz_digits        = 6,
#'   intensity_digits = 2
#' )
#' }
#'
#' @seealso \code{\link{import_msp}}, \code{\link{convert_to_nist_msp}}
#' @export
export_msp <- function(library,
                       file_path        = "library_export.msp",
                       fields           = c("Name", "PrecursorMZ", "Adduct",
                                            "Formula", "ExactMass",
                                            "InChIKey", "InChI", "SMILES",
                                            "IonMode", "CollisionEnergy",
                                            "InstrumentType", "Instrument",
                                            "IonizationType", "Spectrum_type",
                                            "Comments"),
                       mz_digits        = 4L,
                       intensity_digits = 0L) {

  # ── Input validation ───────────────────────────────────────────────────────
  if (!is.list(library) || length(library) == 0)
    stop("`library` must be a non-empty named list of compounds.")
  if (!is.character(file_path) || nchar(trimws(file_path)) == 0)
    stop("`file_path` must be a non-empty string.")
  if (!dir.exists(dirname(normalizePath(file_path, mustWork = FALSE))))
    stop("Parent directory does not exist: ", dirname(file_path))
  # ──────────────────────────────────────────────────────────────────────────

  con <- file(file_path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  n_written <- 0L
  n_skipped <- 0L

  for (cmp in library) {

    peaks <- cmp$Peaks

    # Skip entries with no usable peak table
    if (is.null(peaks) || !is.data.frame(peaks) || nrow(peaks) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }

    # ── Metadata fields ──────────────────────────────────────────────────────
    for (fld in fields) {
      val <- cmp[[fld]]
      if (!is.null(val) && !is.na(val) && nchar(trimws(as.character(val))) > 0)
        writeLines(paste0(fld, ": ", trimws(as.character(val))), con)
    }

    # ── Peak list ────────────────────────────────────────────────────────────
    writeLines(paste0("Num Peaks: ", nrow(peaks)), con)

    peak_lines <- paste(
      formatC(peaks$mz,            format = "f", digits = mz_digits),
      formatC(peaks$rel_intensity, format = "f", digits = intensity_digits)
    )
    writeLines(peak_lines, con)

    # Blank line between entries (NIST convention)
    writeLines("", con)

    n_written <- n_written + 1L
  }

  if (n_skipped > 0)
    message("[export_msp] Skipped ", n_skipped,
            " entries with missing or empty peak table.")

  message("[export_msp] Wrote ", n_written, " compounds to: ", file_path)

  invisible(file_path)
}
