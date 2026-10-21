# =============================================================================
# convert_to_nist_msp.R
# Step 2 of the conversion workflow.
# One-time normalization of any in-house library to clean NIST MSP format.
#
# Supported source formats:
#   - Bruker .library  (non-standard fields: IonPolarity, InstName, ColEnergy,
#                        InstType, IoniMethod, MSMS, ReserveDWord, etc.)
#   - CSU-PMF .mspLib  (metadata encoded in SYNON $:)
# =============================================================================


#' Convert a parsed library to standard NIST MSP format
#'
#' @description
#' One-time conversion of any in-house library (Bruker \code{.library} or
#' CSU-PMF \code{.mspLib}) to a clean NIST MSP file suitable for Zenodo
#' deposition and downstream use with \code{\link{import_msp}}.
#'
#' What this function does:
#' \itemize{
#'   \item Normalizes field names to NIST standard
#'     (\code{IonPolarity: pos} → \code{Ion_mode: POSITIVE},
#'      \code{InstName} → \code{Instrument}, etc.)
#'   \item Decodes CSU-PMF \code{SYNON} codes to proper fields
#'     (\code{$:11Pos} → \code{Ion_mode: POSITIVE},
#'      \code{$:05} → \code{RetentionTime}, etc.)
#'   \item Infers \code{Spectrum_type: MS1} when the \code{SYNON $:00} value
#'     is \code{in-source} and no \code{MSMS} field is present
#'   \item Calculates \code{ExactMass} from \code{Formula} via
#'     \code{\link{calculate_exact_mass}} if absent (common in Bruker files)
#'   \item Writes peaks one pair per line (NIST standard)
#'   \item Drops Bruker-specific internal fields (\code{ReserveDWord},
#'     \code{AnalName}, \code{BasePeak}, \code{PeakWidth}, etc.)
#'   \item Skips optional identifiers (InChIKey, InChI, SMILES) if absent —
#'     entries do not need to match MoNA completeness
#' }
#'
#' @param library A parsed library list — output of \code{\link{read_library}}.
#' @param output_path Path for the output \code{.msp} file.
#' @param source_lab Optional contributor/lab name appended to
#'   \code{Comments:} (e.g. \code{"Your Lab, Your University"}).
#' @param source_instrument Optional fallback instrument string used when no
#'   instrument field is present in an entry.
#'
#' @return Invisibly returns \code{output_path}.
#'
#' @details
#' \strong{Three-step workflow:}
#' \enumerate{
#'   \item \code{\link{read_library}} — tolerates any MSP-like format
#'   \item \code{convert_to_nist_msp} — one-time normalization; run once and
#'     save the resulting \code{.msp} for all future use
#'   \item \code{\link{import_msp}} — strict NIST MSP importer with filtering,
#'     entropy QC, and deduplication for downstream analysis
#' }
#'
#' @examples
#' \dontrun{
#' # ---- Bruker .library → NIST MSP ----
#' bruker_raw <- read_library("inhouse_library.library")
#' convert_to_nist_msp(bruker_raw, "inhouse_library.msp",
#'                     source_lab = "Your Lab, Your University")
#'
#' # ---- CSU-PMF .mspLib → NIST MSP ----
#' csu_raw <- read_library("CSUPMF_C8_lib.mspLib")
#' convert_to_nist_msp(csu_raw, "CSUPMF_C8.msp",
#'                     source_lab = "Broeckling Lab, Colorado State University")
#'
#' # ---- From now on, always load the clean .msp files via import_msp() ----
#' inhouse_lib <- import_msp("inhouse_library.msp",
#'                            spectrum_type = "MS1",
#'                            dedup_by      = "none")
#' csu_lib     <- import_msp("CSUPMF_C8.msp",
#'                            spectrum_type = "MS1",
#'                            dedup_by      = "none")
#' }
#'
#' @seealso \code{\link{read_library}}, \code{\link{import_msp}},
#'   \code{\link{calculate_exact_mass}}
#' @export
convert_to_nist_msp <- function(library,
                                output_path,
                                source_lab        = NULL,
                                source_instrument = NULL) {

  if (length(library) == 0) stop("Library is empty.")

  # Bruker-specific and internal processing fields
  drop_fields <- c(
    "Peaks", "SYNON", "ReserveDWord", "AnalName",
    "BasePeak", "BasePeakIntensity", "PeakWidth",
    "IonPolarity", "MSMS", "InstType", "InstName", "IoniMethod",
    "ColEnergy", "ColEnergy2", "MassRange", "NumPeaks", "Num Peaks",
    "spectral_entropy", "normalized_entropy", "quality",
    "Ion_mode", "Spectrum_type",
    "ExactMass", "MW", "Formula",
    "Name"
  )

  all_lines <- vector("list", length(library))

  for (i in seq_along(library)) {
    compound <- library[[i]]
    buf <- character(200)
    j   <- 0L
    add <- function(x) { j <<- j + 1L; buf[j] <<- x }

    # Decode CSU-PMF SYNON codes
    synon_meta <- .parse_synon(compound$SYNON)

    # ---- 1. Name (required) ----
    add(paste0("Name: ", compound$Name %||% "Unknown"))

    # ---- 2. Formula ----
    formula <- compound$Formula
    if (!is.null(formula) && !is.na(formula))
      add(paste0("Formula: ", formula))

    # ---- 3. ExactMass----
    em <- NA_real_
    if (!is.null(compound$ExactMass) && !is.na(suppressWarnings(as.numeric(compound$ExactMass)))) {
      em <- as.numeric(compound$ExactMass)
    } else if (!is.null(formula) && !is.na(formula)) {
      em <- tryCatch(calculate_exact_mass(formula)[[1]], error = function(e) NA_real_)
    }
    if (!is.na(em)) add(paste0("ExactMass: ", em))

    # ---- 4. MW ----
    if (!is.null(compound$MW) && !is.na(compound$MW))
      add(paste0("MW: ", compound$MW))

    # ---- 5. Ion mode ----
    ion_mode <- .normalize_ion_mode(
      compound$Ion_mode %||% compound$IonPolarity %||% synon_meta$Ion_mode
    )
    if (!is.null(ion_mode)) add(paste0("Ion_mode: ", ion_mode))

    # ---- 6. Spectrum type: compound field default MS1 ----
    spec_type <- .normalize_spectrum_type(
      compound$Spectrum_type %||%
        (if (!is.null(compound$MSMS)) "MS2" else NULL) %||%
        synon_meta$Spectrum_type %||%
        "MS1"
    )
    add(paste0("Spectrum_type: ", spec_type))

    # ---- 7. Instrument type ----
    inst_type <- compound$Instrument_type %||% compound$InstType %||%
      synon_meta$Instrument_type %||% source_instrument
    if (!is.null(inst_type)) add(paste0("Instrument_type: ", inst_type))

    # ---- 8. Instrument name ----
    inst_name <- compound$Instrument %||% compound$InstName %||%
      synon_meta$Instrument %||% source_instrument
    if (!is.null(inst_name)) add(paste0("Instrument: ", inst_name))

    # ---- 9. Ionization method ----
    ioni <- compound$Ionization_method %||% compound$IoniMethod %||%
      synon_meta$Ionization_method
    if (!is.null(ioni)) add(paste0("Ionization_method: ", ioni))

    # ---- 10. MS2-specific fields ----
    pre_mz <- compound$PrecursorMZ %||% compound$PreIon %||% synon_meta$PrecursorMZ
    if (!is.null(pre_mz)) add(paste0("PrecursorMZ: ", pre_mz))

    pre_type <- compound$Precursor_type
    if (!is.null(pre_type)) add(paste0("Precursor_type: ", pre_type))

    col_en <- compound$Collision_energy %||% compound$ColEnergy %||%
      synon_meta$Collision_energy
    if (!is.null(col_en)) add(paste0("Collision_energy: ", col_en))

    # ---- 11. Retention time ----
    rt <- compound$RetentionTime %||% compound$RetTime %||%
      synon_meta$RetentionTime %||%
      .parse_rt_from_comment(compound$Comment %||% compound$Comments)
    if (!is.null(rt)) add(paste0("RetentionTime: ", rt))

    # ---- 12. Identifiers (only if present) ----
    for (id in c("InChIKey", "InChI", "SMILES", "CAS", "DB#", "Splash")) {
      val <- compound[[id]]
      if (!is.null(val) && !is.na(val) && nchar(trimws(val)) > 0)
        add(paste0(id, ": ", val))
    }

    # ---- 13. Spectral entropy and quality (informational) ----
    if (!is.null(compound$spectral_entropy) && !is.na(compound$spectral_entropy))
      add(paste0("Spectral_entropy: ", round(as.numeric(compound$spectral_entropy), 6)))
    if (!is.null(compound$quality) && !is.na(compound$quality))
      add(paste0("Quality: ", compound$quality))

    # ---- 14. Comments ----
    comment_parts <- c()
    raw_comment <- compound$Comment %||% compound$Comments
    if (!is.null(raw_comment) && nchar(trimws(raw_comment)) > 0)
      comment_parts <- c(comment_parts, trimws(raw_comment))
    if (!is.null(source_lab))
      comment_parts <- c(comment_parts, paste0("Contributor=", source_lab))
    if (length(comment_parts) > 0)
      add(paste0("Comments: ", paste(comment_parts, collapse = "  ")))

    # ---- 15. Peaks — one mz/intensity pair per line (NIST standard) ----
    if (!is.null(compound$Peaks) && nrow(compound$Peaks) > 0) {
      pk <- compound$Peaks
      add(paste0("Num Peaks: ", nrow(pk)))
      peak_lines <- paste(pk$mz, round(pk$rel_intensity, 4))
      all_lines[[i]] <- c(buf[seq_len(j)], peak_lines, "")
    } else {
      all_lines[[i]] <- c(buf[seq_len(j)], "")
    }
  }

  writeLines(unlist(all_lines), output_path, useBytes = FALSE)

  message("[convert_to_nist_msp] Written ", length(library),
          " compounds to ", output_path)
  invisible(output_path)
}
