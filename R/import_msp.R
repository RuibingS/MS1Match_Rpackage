# =============================================================================
# Import msp file
# Step 3 of the library file conversion workflow.
# Strict importer for clean NIST MSP files with filtering and
# deduplication.
# =============================================================================

#' Import a standardized NIST MSP spectral library
#'
#' @description
#' Full-featured importer for clean NIST-format \code{.msp} files. Supports
#' filtering by spectrum type, ion polarity, ionization method, collision
#' energy and instrument name. Applies smart deduplication — all implemented
#' in C++ for performance.
#'
#' Use this function for:
#' \itemize{
#'   \item MoNA \code{.msp} downloads
#'   \item Files previously converted with \code{\link{convert_to_nist_msp}}
#'   \item Any NIST-standard \code{.msp} file
#' }
#'
#' For raw non-standard files (Bruker \code{.library}, CSU-PMF \code{.mspLib}),
#' use \code{\link{read_library}} → \code{\link{convert_to_nist_msp}} first.
#'
#' @param file_path Path to the NIST \code{.msp} file.
#' @param spectrum_type Optional filter: \code{"MS1"} or \code{"MS2"}.
#'   MoNA entries labelled \code{Spectrum_type: MS} match \code{"MS1"}.
#' @param ionPolarity Optional filter: \code{"POSITIVE"} or \code{"NEGATIVE"}
#'   (also accepts \code{"pos"}, \code{"neg"}, \code{"P"}, \code{"N"}).
#' @param ioniMethod Optional filter for ionization method (string match).
#' @param colEnergy Optional filter for collision energy (string match).
#' @param instName Optional filter for instrument name (string match).
#' @param filter Optional filter for metadata fields. Accepts either:
#'   \itemize{
#'     \item A \strong{single} named list with elements \code{field} and
#'       \code{pattern} — drops entries whose field does not contain the
#'       pattern:
#'       \preformatted{filter = list(field   = "Instrument_type",
#'                    pattern = "LC-ESI-QTOF")}
#'     \item A \strong{list of} such named lists for multiple conditions
#'       (AND logic — all conditions must be satisfied):
#'       \preformatted{filter = list(
#'         list(field = "Instrument_type", pattern = "LC-ESI-QTOF"),
#'         list(field = "Instrument",      pattern = "Agilent 6520 QTOF")
#'       )}
#'   }
#'   Match is case-sensitive and fixed.
#' @param dedup_by Deduplication grouping key:
#'   \describe{
#'     \item{\code{"name_adduct_mode_energy"}}{(default) One entry per
#'       Name + Adduct + IonMode + CollisionEnergy. Preserves fragmentation
#'       diversity across energies — recommended for MoNA MS2.}
#'     \item{\code{"name_adduct_mode"}}{One entry per Name + Adduct + IonMode.}
#'     \item{\code{"name_only"}}{One entry per compound name — recommended
#'       for CSU-PMF and other MS1 libraries where collision energy variants
#'       are not meaningful. Keeps the entry with the most peaks.}
#'     \item{\code{"none"}}{No deduplication — keep all entries.
#'       Recommended for small in-house libraries.}
#'   }
#' @param require_formula Logical. If \code{TRUE}, entries with a missing or
#'   invalid formula (\code{"unk"}, empty) are discarded after parsing.
#'   Default \code{FALSE}. Set to \code{TRUE} when exact mass is required
#'   for downstream matching or annotation — for example, CSU-PMF MS1
#'   libraries where entries without a formula cannot be used.
#' @param reserve_n Memory hint: expected number of compounds (default 60000).
#'   Increase for very large files to reduce internal reallocations.
#' @param max_records Integer. Stop after this many records (-1 = no limit).
#'
#' @return A named list of compounds. Each element contains standardized
#'   metadata fields and a \code{Peaks} data frame (\code{mz},
#'   \code{rel_intensity}).
#'
#' @examples
#' \dontrun{
#' # MoNA MS1 QTOF library — single instrument type filter
#' mona_lib <- import_msp(
#'   system.file("extdata", "Massbank_MS1",
#'               "massbank_result.nist_msp_MS1_TOF.msp", package = "MS1Match"),
#'   spectrum_type = "MS1",
#'   ionPolarity   = "POSITIVE",
#'   dedup_by      = "name_only",
#'   filter        = list(field = "Instrument_type", pattern = "LC-ESI-QTOF")
#' )
#'
#' # Filter by both instrument type AND specific instrument model
#' # Keeps only Agilent 6520 QTOF entries
#' mona_agilent <- import_msp(
#'   "ms1_tof_positive_filtered.msp",
#'   spectrum_type = "MS1",
#'   ionPolarity   = "pos",
#'   filter = list(
#'     list(field = "Instrument_type", pattern = "LC-ESI-QTOF"),
#'     list(field = "Instrument",      pattern = "Agilent 6520 QTOF")
#'   )
#' )
#'
#' # Filter by instrument type only (keeps all LC-ESI-QTOF regardless of model)
#' mona_qtof <- import_msp(
#'   "ms1_tof_positive_filtered.msp",
#'   spectrum_type = "MS1",
#'   ionPolarity   = "pos",
#'   filter        = list(field = "Instrument_type", pattern = "LC-ESI-QTOF")
#' )
#'
#' # No filter — import everything in the pre-filtered MSP
#' mona_all <- import_msp(
#'   "ms1_tof_positive_filtered.msp",
#'   spectrum_type = "MS1",
#'   ionPolarity   = "pos"
#' )
#'
#' # In-house MS1 library (after conversion from Bruker .library)
#' inhouse_lib <- import_msp(
#'   system.file("extdata", "CBIO", "CBIO_Standards.msp", package = "MS1Match"),
#'   spectrum_type = "MS1",
#'   dedup_by      = "name_only"
#' )
#'
#' # CSU-PMF MS1 library — one entry per compound, formula required
#' csu_lib <- import_msp(
#'   system.file("extdata", "CSUPMF_C8_lib", "CSUPMF_library.msp",
#'               package = "MS1Match"),
#'   spectrum_type   = "MS1",
#'   dedup_by        = "name_only",
#'   require_formula = TRUE
#' )
#'
#' # Access individual compounds
#' inhouse_lib[["Caffeine"]]$ExactMass
#' inhouse_lib[["Caffeine"]]$Peaks
#' }
#'
#' @seealso \code{\link{read_library}}, \code{\link{convert_to_nist_msp}},
#'   \code{\link{calculate_exact_mass}}
#' @export
import_msp <- function(file_path,
                       spectrum_type   = NULL,
                       ionPolarity     = NULL,
                       ioniMethod      = NULL,
                       colEnergy       = NULL,
                       instName        = NULL,
                       filter          = NULL,
                       dedup_by        = c("name_adduct_mode_energy",
                                           "name_adduct_mode",
                                           "name_only",
                                           "none"),
                       require_formula = FALSE,
                       reserve_n       = 60000L,
                       max_records     = -1L) {

  if (!file.exists(file_path)) stop("File not found: ", file_path)

  # ── Validate filter argument ─────────────────────────────────────────────
  if (!is.null(filter)) {
    # Normalise: wrap a single filter into a list
    if (is.list(filter) && all(c("field", "pattern") %in% names(filter))) {
      filter <- list(filter)
    }

    for (i in seq_along(filter)) {
      f <- filter[[i]]
      if (!is.list(f) || !all(c("field", "pattern") %in% names(f)))
        stop("Each element of `filter` must be a named list with ",
             "`field` and `pattern`. Problem at index ", i, ".")
      if (!is.character(f$field)   || nchar(trimws(f$field))   == 0)
        stop("`filter[[", i, "]]$field` must be a non-empty string.")
      if (!is.character(f$pattern) || nchar(trimws(f$pattern)) == 0)
        stop("`filter[[", i, "]]$pattern` must be a non-empty string.")
    }
  }


  dedup_by <- match.arg(dedup_by)

  parsed <- parse_msp_library_cpp(
    file_path,
    spec     = spectrum_type,
    pol      = ionPolarity,
    meth     = ioniMethod,
    cen      = colEnergy,
    inst     = instName,
    dedup    = dedup_by,
    reserven = as.integer(reserve_n),
    maxrec   = as.integer(max_records)
  )

  if (length(parsed) == 0) {
    warning("No compounds returned — check filter arguments.")
    return(list())
  }

  # ── pattern filter (multiple conditions) ───────────
  if (!is.null(filter)) {
    for (f in filter) {
      keep <- sapply(parsed, function(x) {
        val <- x[[f$field]]
        !is.null(val) && grepl(f$pattern, val, fixed = TRUE)
      })
      n_dropped <- sum(!keep)
      if (n_dropped > 0)
        message("[import_msp] Dropped ", n_dropped,
                " entries where ", f$field,
                " does not contain \"", f$pattern, "\".")
      parsed <- parsed[keep]
    }
  }

  if (length(parsed) == 0) {
    warning("No compounds returned — check filter arguments.")
    return(list())
  }

  if (require_formula) {
    has_formula <- sapply(parsed, function(x) {
      f <- x$Formula
      !is.null(f) && nchar(trimws(f)) > 0 && tolower(trimws(f)) != "unk"
    })
    n_dropped <- sum(!has_formula)
    if (n_dropped > 0)
      message("[import_msp] Dropped ", n_dropped,
              " entries with missing or invalid formula.")
    parsed <- parsed[has_formula]
  }

  if (length(parsed) == 0) {
    warning("No compounds remaining after formula filter.")
    return(list())
  }

  nms <- vapply(parsed, .safe_extract_name, character(1))
  names(parsed) <- make.unique(nms, sep = ".")

  message("[import_msp] Final library size: ", length(parsed), " compounds.")

  parsed
}
