# =============================================================================
# utils.R
# Internal helpers shared across read_library(), convert_to_nist_msp(),
# and import_msp().
# =============================================================================


# -----------------------------------------------------------------------------
# Null-coalescing operator
# Returns 'a' if it is non-NULL, non-empty, and not all-NA; otherwise 'b'.
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b


# -----------------------------------------------------------------------------
# .safe_extract_name()
# Safely extracts the Name field from a parsed compound list element.
# -----------------------------------------------------------------------------
.safe_extract_name <- function(x) {
  nm <- x$Name
  if (is.null(nm) || length(nm) == 0) return("Unknown")
  if (is.list(nm))  nm <- unlist(nm)        # flatten any accidental nesting
  nm <- as.character(nm)
  nm <- nm[!is.na(nm) & nzchar(nm)]
  if (length(nm) == 0) return("Unknown")
  nm[[1L]]                                  # take first non-empty value
}


# -----------------------------------------------------------------------------
# .parse_synon()
# Decodes CSU-PMF SYNON $:XX codes into named metadata fields.
#
# NIST standard SYNON encoding used by CSU-PMF .mspLib files:
#   $:00 = spectrum type       $:04 = collision energy
#   $:05 = retention time      $:06 = instrument type
#   $:07 = instrument name     $:09 = column / LC method
#   $:10 = ionization method   $:11 = ion polarity
#   $:12 = collision gas       $:14 = mass range
#   $:16 = precursor m/z
# -----------------------------------------------------------------------------
.parse_synon <- function(synon_vec) {
  if (is.null(synon_vec) || length(synon_vec) == 0) return(list())

  out <- list()

  for (s in synon_vec) {
    s <- trimws(s)
    if (!startsWith(s, "$:")) next

    # Extract 2-digit code and value: "$:11Pos" → code="11", value="Pos"
    code  <- substr(s, 3, 4)
    value <- trimws(substr(s, 5, nchar(s)))

    if (nchar(value) == 0) next

    switch(code,
           "00" = { out$Spectrum_type     <- value },
           "04" = { out$Collision_energy  <- value },
           "05" = { out$RetentionTime     <- value },
           "06" = { out$Instrument_type   <- value },
           "07" = { out$Instrument        <- value },
           "09" = { out$Column            <- value },
           "10" = { out$Ionization_method <- value },
           "11" = { out$Ion_mode          <- .normalize_ion_mode(value) },
           "12" = { out$Collision_gas     <- value },
           "14" = { out$Mass_range        <- value },
           "16" = { out$PrecursorMZ       <- value }
    )
  }
  out
}


# -----------------------------------------------------------------------------
# .parse_rt_from_comment()
# Extracts retention time from a Comment string, e.g. "Rt=49.88" or "RT=603.2"
# Returns the value as a character string, or NULL if not found.
# -----------------------------------------------------------------------------
.parse_rt_from_comment <- function(comment_str) {
  if (is.null(comment_str)) return(NULL)
  m <- regmatches(comment_str,
                  regexpr("(?i)\\bRt?\\s*=\\s*([0-9.]+)", comment_str, perl = TRUE))
  if (length(m) == 0) return(NULL)
  val <- regmatches(m, regexpr("[0-9.]+$", m))
  if (length(val) == 0) return(NULL)
  val
}


# -----------------------------------------------------------------------------
# .normalize_ion_mode()
# Normalizes any ion polarity string to "POSITIVE" or "NEGATIVE".
# -----------------------------------------------------------------------------
.normalize_ion_mode <- function(x) {
  if (is.null(x)) return(NULL)
  v <- toupper(trimws(x))
  if (v %in% c("POS", "P", "POSITIVE")) return("POSITIVE")
  if (v %in% c("NEG", "N", "NEGATIVE")) return("NEGATIVE")
  v
}


# -----------------------------------------------------------------------------
# .normalize_spectrum_type()
# Normalizes any spectrum type string to "MS1" or "MS2".
# -----------------------------------------------------------------------------
.normalize_spectrum_type <- function(x) {
  if (is.null(x)) return(NULL)
  v <- toupper(trimws(x))
  if (v %in% c("MS", "MS1", "IN-SOURCE")) return("MS1")
  if (v == "MS2") return("MS2")
  v
}
