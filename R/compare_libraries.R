# =============================================================================
# compare_libraries.R
#
# R wrappers around fast Rcpp backends for spectral library comparison.
# Source this file after compiling the C++ backend:
#
#   Rcpp::sourceCpp("compare_libraries_cpp.cpp")
#   source("compare_libraries.R")
#
# What R handles:              What C++ handles:
#   safe_utf8()   encoding       classify_match_cpp()   row classification loop
#   extract_lib() field extract  compare_libraries_cpp() full pairwise hash-map
#   merge()       formula join
#
# Match priority (highest to lowest confidence):
#   "Formula+InChIKey"              same formula + full 27-char InChIKey  ✓✓✓
#   "Formula+InChIKey_connectivity" same formula + 14-char block-1 match  ✓✓
#   "Formula+Name"                  same formula + exact name match        ✓
#   "Formula+Synonym"               same formula + synonym match            ~
#   "Formula_only"                  same formula only — review manually     ?
#
# Library synonym settings:
#   CSU      → use_synonyms = FALSE  (Synon = instrument metadata)
#   MassBank → use_synonyms = FALSE  (Synon = metadata tokens)
#   CBIO     → use_synonyms = TRUE   (Synon = real compound names)
#
# Library InChIKey availability:
#   CSU      → no InChIKey field
#   MassBank → has InChIKey
#   CBIO     → has InChIKey
# =============================================================================


# -----------------------------------------------------------------------------
# safe_utf8()  —  coerce character vector to valid UTF-8
#
# Prevents "invalid multibyte string" errors from trimws()/tolower() when
# library entries contain Latin-1 / Windows-1252 encoded characters.
# -----------------------------------------------------------------------------
safe_utf8 <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- as.character(x)
  converted <- iconv(x, from = "latin1", to = "UTF-8")
  bad <- is.na(converted)
  if (any(bad))
    converted[bad] <- iconv(x[bad], from = "latin1", to = "UTF-8", sub = "")
  converted
}


# -----------------------------------------------------------------------------
# extract_lib()  —  extract and normalise fields from a parsed MSP library
#
# Returns a data.frame with both display columns (original values) and
# normalised columns (lowercased/trimmed/cleaned) ready for C++ matching.
#
# Parameters:
#   lib          — list of compound entries from import_msp()
#   use_synonyms — if FALSE, synonym field is ignored (set to "")
#
# Columns returned:
#   formula      — original formula string
#   inchikey     — original InChIKey (uppercase, "" if absent/invalid)
#   name_display — original Name (for output)
#   syn_display  — original synonyms joined by ";" (for output)
#   name_norm    — Name lowercased + trimmed (for matching)
#   syn_norm     — synonyms lowercased + trimmed + joined by ";" (for matching)
#   inchikey_norm— InChIKey trimmed uppercase (InChIKey is case-sensitive,
#                  not lowercased)
# -----------------------------------------------------------------------------
extract_lib <- function(lib, use_synonyms = TRUE) {

  get_field <- function(x, field) {
    val <- x[[field]]
    if (is.null(val)) NA_character_ else safe_utf8(val)
  }

  # Validate InChIKey format: 27 chars, XXXXXXXXXXXXXX-XXXXXXXXXX-X
  clean_ik <- function(x) {
    x <- trimws(safe_utf8(x))
    ifelse(!is.na(x) & nchar(x) == 27 &
             grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", x),
           x, NA_character_)
  }

  get_synonyms_str <- function(entry, use_syn) {
    if (!use_syn) return("")
    syn <- unlist(entry[grep("^Synon", names(entry), ignore.case = TRUE)])
    if (length(syn) == 0) return("")
    syn <- safe_utf8(syn)
    syn <- syn[!is.na(syn)]
    # Strip metadata tokens (MassBank/CSU junk)
    syn <- syn[!grepl("^\\$:|^in-source|^\\$:0", syn, ignore.case = TRUE)]
    if (length(syn) == 0) "" else paste(trimws(tolower(syn)), collapse = ";")
  }

  formulas     <- sapply(lib, function(x) {
    f <- get_field(x, "Formula")
    if (is.na(f)) "" else trimws(f)
  })

  ik_raw       <- sapply(lib, function(x) get_field(x, "InChIKey"))
  ik_clean     <- clean_ik(ik_raw)
  ik_display   <- ifelse(is.na(ik_clean), "", ik_clean)

  names_raw    <- sapply(lib, function(x) get_field(x, "Name"))
  names_display <- ifelse(is.na(names_raw), "", names_raw)
  names_norm   <- trimws(tolower(ifelse(is.na(names_raw), "", names_raw)))

  syn_norm     <- sapply(lib, function(x) get_synonyms_str(x, use_synonyms))

  # Display synonyms (original case, for output columns)
  syn_display  <- sapply(lib, function(x) {
    if (!use_synonyms) return("")
    syn <- unlist(x[grep("^Synon", names(x), ignore.case = TRUE)])
    if (length(syn) == 0) return("")
    syn <- safe_utf8(syn)
    syn <- syn[!is.na(syn)]
    syn <- syn[!grepl("^\\$:|^in-source|^\\$:0", syn, ignore.case = TRUE)]
    if (length(syn) == 0) "" else paste(trimws(syn), collapse = "; ")
  })

  data.frame(
    formula       = formulas,
    inchikey      = ik_display,
    name_display  = names_display,
    syn_display   = syn_display,
    name_norm     = names_norm,
    syn_norm      = syn_norm,
    inchikey_norm = ik_display,   # InChIKey not lowercased — case sensitive
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# compare_libraries()
#
# Fast pairwise comparison of two libraries using the C++ hash-map backend.
# Replaces the slow nested R loop — O(n*m) → O(n + m + matches).
#
# Parameters:
#   lib1, lib2         — parsed MSP library lists
#   use_synonyms_lib1  — use Synon field for lib1? (FALSE for CSU/MassBank)
#   use_synonyms_lib2  — use Synon field for lib2? (FALSE for CSU/MassBank)
#
# Returns data.frame with columns:
#   Formula, InChIKey_lib1, InChIKey_lib2,
#   Name_lib1, Synonyms_lib1, Name_lib2, Synonyms_lib2,
#   Match_type
# =============================================================================
compare_libraries <- function(lib1, lib2,
                              use_synonyms_lib1 = TRUE,
                              use_synonyms_lib2 = TRUE) {

  d1 <- extract_lib(lib1, use_synonyms_lib1)
  d2 <- extract_lib(lib2, use_synonyms_lib2)

  result <- compare_libraries_cpp(
    # normalised (for matching)
    formula1 = d1$formula,      formula2 = d2$formula,
    ik1      = d1$inchikey_norm, ik2     = d2$inchikey_norm,
    name1    = d1$name_norm,     name2   = d2$name_norm,
    syn1     = d1$syn_norm,      syn2    = d2$syn_norm,
    # original (for display)
    orig_name1 = d1$name_display, orig_name2 = d2$name_display,
    orig_syn1  = d1$syn_display,  orig_syn2  = d2$syn_display,
    orig_ik1   = d1$inchikey,     orig_ik2   = d2$inchikey
  )

  # Replace empty strings with NA in InChIKey output columns for clarity
  result$InChIKey_lib1[result$InChIKey_lib1 == ""] <- NA
  result$InChIKey_lib2[result$InChIKey_lib2 == ""] <- NA

  result
}


# =============================================================================
# compare_libraries_3()
#
# Compare two or three libraries. Uses merge() for the formula join (fast),
# then calls classify_match_cpp() to classify each merged row — replacing
# the slow apply() loop with a C++ vector operation.
#
# Parameters:
#   lib1, lib2, lib3   — parsed MSP library lists (lib3 optional)
#   lib_names          — display names for output column suffixes
#   use_synonyms       — logical vector of length 3, one per library
#
# Returns data.frame with columns:
#   Formula,
#   InChIKey_<lib>, Name_<lib>, Synonyms_<lib>  (one set per library)
#   In_<lib>                                     (presence flags)
#   Match_type
# =============================================================================
compare_libraries_3 <- function(lib1, lib2, lib3 = NULL,
                                lib_names    = c("Library1", "Library2", "Library3"),
                                use_synonyms = c(TRUE, TRUE, TRUE)) {

  n_libs     <- if (is.null(lib3)) 2L else 3L
  used_names <- lib_names[seq_len(n_libs)]
  libs       <- list(lib1, lib2, lib3)[seq_len(n_libs)]

  # --- Extract and normalise each library ---
  extracted <- mapply(extract_lib, libs, use_synonyms[seq_len(n_libs)],
                      SIMPLIFY = FALSE)

  # --- Build per-library data.frames with suffixed column names ---
  dfs <- lapply(seq_len(n_libs), function(i) {
    d  <- extracted[[i]]
    nm <- used_names[i]
    df <- data.frame(
      Formula                      = d$formula,
      # matching columns (normalised) — used by C++, dropped after merge
      .ik_norm    = d$inchikey_norm,
      .name_norm  = d$name_norm,
      .syn_norm   = d$syn_norm,
      # display columns — kept in output
      stringsAsFactors = FALSE
    )
    # Add display columns with suffixed names
    df[[paste0("InChIKey_",  nm)]] <- d$inchikey
    df[[paste0("Name_",      nm)]] <- d$name_display
    df[[paste0("Synonyms_",  nm)]] <- d$syn_display
    df
  })

  # --- Inner join on Formula (only compounds in ALL libraries) ---
  merged <- Reduce(
    function(x, y) merge(x, y, by = "Formula", all = FALSE,
                         suffixes = paste0(".", seq_len(n_libs)[1:2])),
    dfs
  )

  if (nrow(merged) == 0) {
    message("No formula matches found across all libraries.")
    return(merged)
  }

  # --- Presence flags ---
  for (nm in used_names)
    merged[[paste0("In_", nm)]] <- !is.na(merged[[paste0("Name_", nm)]])

  # --- Identify normalised column names after merge ---
  # After Reduce merge, norm columns get suffixes like .ik_norm.1 / .ik_norm.2
  # We need to retrieve them reliably regardless of suffix pattern
  get_merged_col <- function(merged, prefix, idx) {
    candidates <- grep(paste0("^\\", prefix, "(\\.|$)"), names(merged), value = TRUE)
    if (length(candidates) >= idx) candidates[idx] else rep("", nrow(merged))
  }

  # Build per-library normalised vectors for C++
  ik_list   <- lapply(seq_len(n_libs), function(i)
    merged[[grep(paste0("^\\.ik_norm"),   names(merged))[i]]])
  name_list <- lapply(seq_len(n_libs), function(i)
    merged[[grep(paste0("^\\.name_norm"), names(merged))[i]]])
  syn_list  <- lapply(seq_len(n_libs), function(i)
    merged[[grep(paste0("^\\.syn_norm"),  names(merged))[i]]])

  # Replace any NAs with empty string for C++
  clean_vec <- function(x) ifelse(is.na(x), "", x)
  ik_list   <- lapply(ik_list,   clean_vec)
  name_list <- lapply(name_list, clean_vec)
  syn_list  <- lapply(syn_list,  clean_vec)

  # --- Classify using C++ (handles 2 or 3 libraries) ---
  # For 3-way: require ALL libraries to agree at each level
  if (n_libs == 2L) {
    merged$Match_type <- classify_match_cpp(
      ik1   = ik_list[[1]],   ik2   = ik_list[[2]],
      name1 = name_list[[1]], name2 = name_list[[2]],
      syn1  = syn_list[[1]],  syn2  = syn_list[[2]]
    )
  } else {
    # 3-way: classify 1vs2, 1vs3, 2vs3 and take the MINIMUM confidence level
    match_levels <- c(
      "Formula+InChIKey"              = 1L,
      "Formula+InChIKey_connectivity" = 2L,
      "Formula+Name"                  = 3L,
      "Formula+Synonym"               = 4L,
      "Formula_only"                  = 5L
    )
    level_names <- names(match_levels)

    m12 <- classify_match_cpp(ik_list[[1]], ik_list[[2]],
                               name_list[[1]], name_list[[2]],
                               syn_list[[1]],  syn_list[[2]])
    m13 <- classify_match_cpp(ik_list[[1]], ik_list[[3]],
                               name_list[[1]], name_list[[3]],
                               syn_list[[1]],  syn_list[[3]])
    m23 <- classify_match_cpp(ik_list[[2]], ik_list[[3]],
                               name_list[[2]], name_list[[3]],
                               syn_list[[2]],  syn_list[[3]])

    # Overall match = weakest of the three pairwise matches
    worst <- pmax(match_levels[m12],
                  match_levels[m13],
                  match_levels[m23])
    merged$Match_type <- level_names[worst]
  }

  drop_cols <- grep("^\\.ik_norm|^\\.name_norm|^\\.syn_norm", names(merged))
  if (length(drop_cols) > 0) merged <- merged[, -drop_cols]

  # Replace empty InChIKey strings with NA for clarity
  for (nm in used_names) {
    col <- paste0("InChIKey_", nm)
    if (col %in% names(merged))
      merged[[col]][merged[[col]] == ""] <- NA
  }

  merged
}
