# =============================================================================
# mirror_plot.R
# Mirror plot functions for MS1 + MS2 spectral matching visualization.
# Depends on: utils.R
#
# Changes vs previous version
# ----------------------------
# 1. extract_nl_label_from_adduct() — extracts a clean NL label from pipeline
#    adduct_match strings like "[M+H-H2O]+ 274.131742".
# 2. build_annot_lookup() — stores an nl_label field per feature; prefers NL
#    entries over plain adduct entries when the same feature appears multiple times.
# 3. build_mirror_df() — NL label from pipeline adduct_match shown on fragments;
#    plain [M+H]+ on fragment peaks is suppressed in favour of NL label.
# 4. plot_one_compound() — inject_precursor_feature() removed; precursor
#    injection is now handled in direct_match.cpp.
# 5. annot_res_ms1/annot_res_ms2 are now data frames (e.g.
#    StM14_annotate_res_ms1/ms2), not lists of per-feature records.
#    build_annot_lookup() and the filter_annot_ms1/ms2 helpers inside
#    plot_one_compound() now subset by ROW instead of iterating a list.
#    Field names changed: adduct_match/adduct_ppm/adduct_mda/
#    adduct_confidence -> grp_anno_match/grp_anno_ppm/grp_anno_mda/
#    grp_anno_confidence (Phase1 fields; exp_anno_* are Phase2-only and
#    unused here). NL labels are still parsed out of grp_anno_match only,
#    same as before; the NL_match column is not used for plot labels.
#    All annotation lookups are restricted internally to
#    annotation_phase == "Phase1".
# =============================================================================

#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import ggrepel
NULL
NULL



#' Extract a string for use as a filename
#' @param x Character string
#' @return Character string with special characters replaced by underscores
sanitize_filename <- function(x) {
  x <- gsub("[/\\:*?\"<>|]", "_", x)
  gsub("\\s+", "_", x)
}


#' Extract precursor m/z from a direct_match entry
#'
#' Tries several common field name spellings in order.
#'
#' @param entry A named list (one element of direct_match_res_MS2)
#' @return Numeric precursor m/z, or \code{NA_real_} if not found
extract_precursor_mz <- function(entry) {
  for (nm in c("Precursor_mz", "precursor_mz", "precursor", "precuror")) {
    if (!is.null(entry[[nm]])) {
      val <- suppressWarnings(as.numeric(entry[[nm]][1]))
      if (!is.na(val)) return(val)
    }
  }
  NA_real_
}


#' Apply intensity scaling transformation
#'
#' @param x_rel Numeric vector of relative intensities (0-100 scale)
#' @param method One of \code{"none"}, \code{"sqrt"}, \code{"log10"}
#' @return Transformed numeric vector
scale_intensity <- function(x_rel, method = "none") {
  method <- match.arg(method, c("none", "sqrt", "log10"))
  if (method == "none")  return(x_rel)
  if (method == "sqrt")  return(sqrt(pmax(x_rel, 0)))
  if (method == "log10") return(log10(x_rel + 1))
}


#' Normalize raw query feature intensities to 0-100 scale
#'
#' @param x Numeric vector of raw intensities
#' @return Numeric vector scaled to 0-100
normalize_query <- function(x) {
  if (length(x) == 0 || all(is.na(x)) || max(x, na.rm = TRUE) <= 0)
    return(rep(0, length(x)))
  x / max(x, na.rm = TRUE) * 100
}



#' Extract a clean neutral-loss label from a pipeline adduct_match string
#'
#' Strips trailing mass values and returns the NL portion only when the
#' adduct is a neutral-loss form.
#' \code{"[M+H-H2O]+ 274.131742"} -> \code{"[M+H-H2O]+"}
#' \code{"[M+H]+ 194.057909"}     -> \code{""}  (plain adduct, not NL)
#'
#' @param x Character scalar, possibly with trailing mass
#' @return Character scalar: the NL label (stripped) or \code{""}
extract_nl_label_from_adduct <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return("")
  stripped <- strip_adduct_mass(x)
  if (grepl("\\[M[^-\\]]*-", stripped, perl=TRUE)) stripped else ""
}


# =============================================================================
# Annotation lookup
# =============================================================================

#' Build a per-feature annotation lookup table
#'
#' Indexes an annotate_res data frame by \code{feature_name}.  Only rows
#' with \code{annotation_phase == "Phase1"} are used (Phase2 \code{_exp}
#' rows are dropped here; call this on the full table and it will filter
#' internally). \code{adduct_match} is taken directly from
#' \code{grp_anno_match} with no group-level fallback: in this data model,
#' \code{grp_anno_match} is already populated whenever a feature genuinely
#' has a Tier1a/Tier1b precursor or isotope match, so a feature with an
#' empty \code{grp_anno_match} should get no label — not one inherited
#' from the compound's overall ion type. (A prior version promoted
#' \code{group_ion_type} into \code{adduct_match} for any row lacking its
#' own annotation; that caused every unannotated fragment in a group to
#' share the same "[M+H]+" label text as the true precursor feature, so
#' the plot's per-label intensity-based dedup would silently keep
#' whichever row was tallest — not necessarily the real precursor — and
#' point the leader line at the wrong peak.) A \code{nl_label} field is
#' also stored: it is non-empty only when \code{grp_anno_match} itself
#' already encodes a neutral loss (e.g. \code{"[M+H-H2O]+ 274.131742"} ->
#' \code{"[M+H-H2O]+"}). \code{NL_match} is not used for plot labels.
#'
#' @param annotate_res    Data frame of per-feature annotation rows
#'   (e.g. \code{StM14_annotate_res_ms1}). Expected columns include
#'   \code{feature_name}, \code{grp_anno_match}, \code{grp_anno_ppm},
#'   \code{grp_anno_mda}, \code{grp_anno_confidence}, \code{iso_annotation},
#'   \code{annotation_phase}.
#' @param group_ion_type  Unused (kept for call-site compatibility).
#' @return Named list keyed by \code{feature_name}
build_annot_lookup <- function(annotate_res, group_ion_type = "") {
  if (is.null(annotate_res) || nrow(annotate_res) == 0)
    return(list())

  # Restrict to Phase1 rows — Phase2 (_exp) rows are out of scope here
  if ("annotation_phase" %in% names(annotate_res)) {
    annotate_res <- annotate_res[!is.na(annotate_res$annotation_phase) &
                                  annotate_res$annotation_phase == "Phase1", ,
                                  drop = FALSE]
  }
  if (nrow(annotate_res) == 0) return(list())

  # Is this adduct string a neutral-loss form?
  .is_nl_adduct <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(FALSE)
    grepl("\\[M[^-\\]]*-", strip_adduct_mass(x), perl=TRUE)
  }

  lookup <- list()
  for (i in seq_len(nrow(annotate_res))) {
    item <- annotate_res[i, ]

    fn <- as.character(item$feature_name)
    if (is.null(fn) || length(fn) == 0 || is.na(fn) || fn == "") next

    adduct <- as.character(item$grp_anno_match %||% "")
    if (is.na(adduct)) adduct <- ""

    # Treat bare unbracketed strings (e.g. "M+H") as empty — not valid notation
    if (!grepl("^\\[", trimws(adduct))) adduct <- ""

    # No group-level fallback here — an empty grp_anno_match means this
    # feature genuinely has no primary-ion annotation, and should show no
    # label rather than inheriting the compound's overall ion type.

    iso_val <- as.character(item$iso_annotation %||% "")
    if (is.na(iso_val)) iso_val <- ""

    new_entry <- list(
      adduct_match      = adduct,
      nl_label          = extract_nl_label_from_adduct(adduct),
      adduct_ppm        = as.character(item$grp_anno_ppm        %||% ""),
      adduct_mda        = as.character(item$grp_anno_mda        %||% ""),
      adduct_confidence = as.character(item$grp_anno_confidence %||% ""),
      iso_annotation    = iso_val
    )

    # Only overwrite an existing entry if the new one is more informative:
    # NL adduct beats plain adduct beats empty. This prevents a later plain
    # [M+H]+ entry (same feature_name, different compound context) from
    # silently erasing an NL annotation that was stored first.
    if (is.null(lookup[[fn]])) {
      lookup[[fn]] <- new_entry
    } else if (.is_nl_adduct(adduct) && !.is_nl_adduct(lookup[[fn]]$adduct_match)) {
      lookup[[fn]] <- new_entry   # upgrade plain -> NL
    }
    # else: keep existing (NL already stored, or first plain wins)
  }
  lookup
}


# =============================================================================
# Core data builder
# =============================================================================

#' Build the mirror plot data frame for one spectrum pair
#'
#' Query intensities are normalized to 0-100; reference intensities are
#' used as-is (already 0-100 from the MSP library).
#'
#' NL labels come from the pipeline \code{adduct_match} field only
#' (e.g. \code{"[M+H-H2O]+ 274.131742"} -> \code{"[M+H-H2O]+"}).
#' The precursor peak is never relabelled.
#'
#' @param query_mz            Numeric vector of matched query feature m/z
#' @param query_int           Numeric vector of matched query feature intensities
#'   (raw scale)
#' @param query_feature_names Character vector of feature IDs
#' @param ref_mz_matched      Numeric vector of matched library peak m/z
#' @param ref_mz_full         Numeric vector of all library peak m/z
#' @param ref_int_full        Numeric vector of all library peak intensities
#'   (0-100 scale)
#' @param annot_lookup        Named list from \code{build_annot_lookup()}
#' @param group_ion_type      Character; group-level Ion_type from direct_match
#' @param scale_method        One of \code{"none"}, \code{"sqrt"}, \code{"log10"}
#' @param unmatched_label     Label for unmatched reference peaks
#' @param precursor_mz        Numeric scalar; used to identify the precursor
#'   peak and restrict \code{[M+H]+} labelling to it.  Pass \code{NA_real_}
#'   for MS1.
#' @return Data frame carrying attribute \code{scale_method}
build_mirror_df <- function(query_mz,
                            query_int,
                            query_feature_names,
                            ref_mz_matched,
                            ref_mz_full,
                            ref_int_full,
                            annot_lookup,
                            group_ion_type,
                            scale_method    = "none",
                            unmatched_label = "Unmatched",
                            precursor_mz    = NA_real_) {

  scale_method <- match.arg(scale_method, c("none", "sqrt", "log10"))

  query_int_rel <- normalize_query(query_int)
  n_q <- length(query_mz)


  ion_types   <- character(n_q)
  adduct_mats <- character(n_q)
  adduct_ppms <- character(n_q)
  adduct_mdas <- character(n_q)
  adduct_conf <- character(n_q)
  iso_annots  <- character(n_q)

  # Index of the peak closest to precursor_mz — used to pin [M+H]+ label
  precursor_idx <- NA_integer_
  if (!is.na(precursor_mz) && length(query_mz) > 0)
    precursor_idx <- which.min(abs(query_mz - precursor_mz))

  # TRUE for plain adducts ([M+H]+, [M+Na]+, …) that carry no NL information
  .is_plain_precursor_adduct <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(TRUE)
    !grepl("\\[M[^-\\]]*-", strip_adduct_mass(x), perl=TRUE)
  }

  for (i in seq_len(n_q)) {
    fn          <- if (i <= length(query_feature_names)) query_feature_names[i] else ""
    pipeline_nl <- ""

    if (!is.null(fn) && fn != "" && !is.null(annot_lookup[[fn]])) {
      a              <- annot_lookup[[fn]]
      # Ion_type is no longer a per-feature annotation column (removed from
      # annotate_res); it only ever comes from the group-level direct_match
      # entry, so every feature in this compound shares the same value.
      ion_types[i]   <- group_ion_type
      adduct_mats[i] <- a$adduct_match
      adduct_ppms[i] <- a$adduct_ppm
      adduct_mdas[i] <- a$adduct_mda
      adduct_conf[i] <- a$adduct_confidence
      iso_annots[i]  <- a$iso_annotation
      pipeline_nl    <- a$nl_label %||% ""

      # Clear [M+H]+ promoted from bare "Precursor" on non-precursor peaks
      if (tolower(ion_types[i]) == "precursor" &&
          adduct_mats[i] != "" &&
          !is.na(precursor_idx) &&
          i != precursor_idx) {
        adduct_mats[i] <- ""
      }
    } else {
      ion_types[i] <- group_ion_type
    }

    # Use pipeline NL label on fragment peaks (not the precursor itself)
    is_precursor_peak <- !is.na(precursor_idx) && i == precursor_idx
    if (!is_precursor_peak && pipeline_nl != "") {
      adduct_mats[i] <- pipeline_nl
    }
  }

  query_df <- data.frame(
    mz                = query_mz,
    int_raw           = query_int,
    int_rel           = query_int_rel,
    source            = "Query",
    category          = "Matched",
    match_id          = seq_len(n_q),
    feature_name      = query_feature_names,
    ion_type          = ion_types,
    adduct_match      = adduct_mats,
    adduct_ppm        = adduct_ppms,
    adduct_mda        = adduct_mdas,
    adduct_confidence = adduct_conf,
    iso_annotation    = iso_annots,
    stringsAsFactors  = FALSE
  )

  n_full  <- length(ref_mz_full)
  ref_cat <- rep(unmatched_label, n_full)
  ref_mid <- rep(NA_integer_,     n_full)
  used    <- rep(FALSE,           n_full)

  for (i in seq_along(ref_mz_matched)) {
    j <- which.min(abs(ref_mz_full - ref_mz_matched[i]))
    if (!used[j]) {
      ref_cat[j] <- "Matched"
      ref_mid[j] <- i
      used[j]    <- TRUE
    }
  }

  ref_df <- data.frame(
    mz                = ref_mz_full,
    int_raw           = ref_int_full,
    int_rel           = ref_int_full,
    source            = "Reference",
    category          = ref_cat,
    match_id          = ref_mid,
    feature_name      = NA_character_,
    ion_type          = NA_character_,
    adduct_match      = NA_character_,
    adduct_ppm        = NA_character_,
    adduct_mda        = NA_character_,
    adduct_confidence = NA_character_,
    iso_annotation    = NA_character_,
    stringsAsFactors  = FALSE
  )

  query_df$int_plot  <- scale_intensity(query_df$int_rel,  scale_method)
  ref_df$int_plot    <- scale_intensity(ref_df$int_rel,    scale_method)
  query_df$intensity <-  query_df$int_plot
  ref_df$intensity   <- -ref_df$int_plot

  out <- dplyr::bind_rows(query_df, ref_df)
  attr(out, "scale_method") <- scale_method
  out
}


# =============================================================================
# plot_one_compound
# =============================================================================

#' Plot MS1 and MS2 mirror panels stacked for one compound
#'
#' @param ms1_entry           One element of \code{direct_match_res_MS1}
#' @param ms2_entry           One element of \code{direct_match_res_MS2}
#' @param annot_res_ms1       Full \code{StM14_annotate_res_ms1} data frame.
#'   Filtered internally to \code{annotation_phase == "Phase1"}.
#' @param annot_res_ms2       Full \code{StM14_annotate_res_ms2} data frame.
#'   Filtered internally to \code{annotation_phase == "Phase1"}.
#' @param scale_method        Intensity scale: \code{"none"}, \code{"sqrt"},
#'   \code{"log10"}
#' @param mz_digits           Decimal places for m/z labels
#' @param padding             Fractional x-axis padding
#' @param precursor_extra_frac Extra right-side padding for precursor marker
#' @param shared_xlim         Logical; share x-axis across both panels
#' @param mz_label_mode       One of \code{"all"}, \code{"primary"},
#'   \code{"adduct"}
#' @return A \code{patchwork} object
plot_one_compound <- function(ms1_entry,
                              ms2_entry,
                              annot_res_ms1,
                              annot_res_ms2,
                              scale_method         = "log10",
                              mz_digits            = 1,
                              padding              = 0.05,
                              precursor_extra_frac = 0.04,
                              shared_xlim          = TRUE,
                              mz_label_mode        = "all") {

  cname <- as.character(ms1_entry$Compound_name[1])
  stopifnot(identical(cname, as.character(ms2_entry$Compound_name[1])))

  precursor_mz <- extract_precursor_mz(ms2_entry)

  # ---- Annotation filtering helpers ----------------------------------------
  # UPDATED: annot_res_ms1/ms2 are now data frames (e.g. StM14_annotate_res_ms1),
  # not lists of per-feature records, so these subset by ROW. Phase1 filtering
  # itself happens inside build_annot_lookup(), but restricting here first
  # keeps the Compound_name / feature_name match cheap and unambiguous.
  filter_annot_ms1 <- function(res, cname) {
    if (is.null(res) || nrow(res) == 0) return(res)
    keep <- !is.na(res$Compound_name) &
      trimws(as.character(res$Compound_name)) == trimws(cname)
    res[keep, , drop = FALSE]
  }

  filter_annot_ms2 <- function(res, feature_names) {
    if (is.null(res) || nrow(res) == 0) return(res)
    keep <- !is.na(res$feature_name) & res$feature_name %in% feature_names
    res[keep, , drop = FALSE]
  }

  lookup_ms1 <- build_annot_lookup(
    filter_annot_ms1(annot_res_ms1, cname),
    group_ion_type = ms1_entry$Ion_type %||% ""
  )

  lookup_ms2 <- build_annot_lookup(
    filter_annot_ms2(annot_res_ms2, ms2_entry$Matched_features_to_peaks),
    group_ion_type = ms2_entry$Ion_type %||% ""
  )

  # ---- Build data frames ---------------------------------------------------
  df_ms1 <- build_mirror_df(
    query_mz            = ms1_entry$Matched_features_mz,
    query_int           = ms1_entry$Matched_features_int,
    query_feature_names = ms1_entry$Matched_features_to_peaks,
    ref_mz_matched      = ms1_entry$Compound_peak_mz,
    ref_mz_full         = ms1_entry$Compound_peak_mz_full,
    ref_int_full        = ms1_entry$Compound_peak_int_full,
    annot_lookup        = lookup_ms1,
    group_ion_type      = ms1_entry$Ion_type %||% "",
    scale_method        = scale_method,
    precursor_mz        = NA_real_          # MS1: no precursor filtering
  )

  df_ms2 <- build_mirror_df(
    query_mz            = ms2_entry$Matched_features_mz,
    query_int           = ms2_entry$Matched_features_int,
    query_feature_names = ms2_entry$Matched_features_to_peaks,
    ref_mz_matched      = ms2_entry$Compound_peak_mz,
    ref_mz_full         = ms2_entry$Compound_peak_mz_full,
    ref_int_full        = ms2_entry$Compound_peak_int_full,
    annot_lookup        = lookup_ms2,
    group_ion_type      = ms2_entry$Ion_type %||% "",
    scale_method        = scale_method,
    precursor_mz        = precursor_mz
  )

  xlim <- if (shared_xlim)
    compute_xlim(ms1_entry$Compound_peak_mz_full,
                 ms2_entry$Compound_peak_mz_full,
                 precursor_mz, padding, precursor_extra_frac)
  else NULL

  p1 <- plot_mirror(df_ms1,
                    xlim_range    = xlim,
                    title         = paste0(cname, " \u2014 MS1"),
                    precursor_mz  = precursor_mz,
                    mz_digits     = mz_digits,
                    mz_label_mode = mz_label_mode)

  p2 <- plot_mirror(df_ms2,
                    xlim_range    = xlim,
                    title         = paste0(cname, " \u2014 MS2"),
                    precursor_mz  = precursor_mz,
                    mz_digits     = mz_digits,
                    mz_label_mode = mz_label_mode)

  p1 / p2
}


# =============================================================================
# Remaining helpers — unchanged from original version
# =============================================================================

#' Compute a shared x-axis range for the MS1 + MS2 mirror plot pair
#'
#' @param ref_mz_ms1_full    Numeric vector; all MS1 library peak m/z
#' @param ref_mz_ms2_full    Numeric vector; all MS2 library peak m/z
#' @param precursor_mz       Numeric scalar or \code{NULL}
#' @param padding_frac       Fractional padding added to each side
#' @param precursor_extra_frac Extra right-side padding for precursor marker
#' @return Length-2 numeric vector \code{c(xlo, xhi)}
compute_xlim <- function(ref_mz_ms1_full,
                         ref_mz_ms2_full,
                         precursor_mz         = NULL,
                         padding_frac         = 0.05,
                         precursor_extra_frac = 0.04) {
  all_mz <- c(ref_mz_ms1_full, ref_mz_ms2_full)
  if (!is.null(precursor_mz) && !is.na(precursor_mz))
    all_mz <- c(all_mz, precursor_mz)
  all_mz <- all_mz[!is.na(all_mz)]
  xmin <- min(all_mz); xmax <- max(all_mz)
  xr   <- max(xmax - xmin, 1)
  xlo  <- xmin - xr * padding_frac
  xhi  <- xmax + xr * padding_frac
  if (!is.null(precursor_mz) && !is.na(precursor_mz)) {
    gap <- xr * precursor_extra_frac
    if ((xhi - precursor_mz) < gap) xhi <- precursor_mz + gap
  }
  c(xlo, xhi)
}


#' Build y-axis breaks and labels for the log10(x+1) intensity scale
#'
#' @param ymax Positive numeric; maximum absolute intensity value on the plot
#' @return Named list with elements \code{breaks} and \code{labels}
build_log10_yaxis <- function(ymax) {
  pct_breaks <- c(1, 10, 100)
  brks_pos   <- log10(pct_breaks + 1)
  brks_pos   <- c(0, brks_pos[brks_pos <= ymax * 1.05])
  brks_neg   <- -rev(brks_pos[brks_pos > 0])
  brks       <- unique(c(brks_neg, brks_pos))
  labs <- sapply(brks, function(b) {
    pct <- round(10^abs(b) - 1, 0)
    if (abs(b) < 1e-10) return("0%")
    paste0(format(pct, big.mark = ","), "%")
  })
  list(breaks = brks, labels = labs)
}


#' Find the first entry in a list matching a compound name
#'
#' @param compound_name Character scalar
#' @param res_list      List of entries each containing \code{Compound_name}
#' @return Integer index or \code{NA_integer_}
find_by_name <- function(compound_name, res_list) {
  nms <- vapply(res_list, function(x) {
    if (is.null(x$Compound_name) || length(x$Compound_name) == 0)
      return(NA_character_)
    as.character(x$Compound_name[1])
  }, character(1))
  idx <- which(nms == compound_name)
  if (length(idx) == 0) NA_integer_ else idx[1]
}


#' Strip trailing mass value from a single adduct token
#'
#' \code{"[M+H]+ 194.057909"} -> \code{"[M+H]+"}
#'
#' @param x Character scalar
#' @return Stripped character scalar
.strip_one_adduct <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return("")
  trimws(sub("(\\]\\d*\\+)\\s*[0-9][0-9.]*\\s*$", "\\1", x))
}


#' Strip trailing mass values from an adduct annotation string
#'
#' Handles semicolon-separated multi-adduct strings.
#'
#' @param x Character scalar, possibly semicolon-separated
#' @return Stripped character scalar
strip_adduct_mass <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return("")
  tokens  <- strsplit(x, ";")[[1]]
  tokens  <- trimws(tokens)
  stripped <- vapply(tokens, .strip_one_adduct, character(1))
  paste(stripped, collapse = "; ")
}


#' Test whether a single adduct token is singly charged
#'
#' @param x Character scalar, single stripped adduct token
#' @return Logical
is_singly_charged <- function(x) {
  !grepl("\\]\\d+\\+\\s*$", x)
}


#' Classify whether a feature qualifies as a primary ion
#'
#' @param ion_type     Character scalar
#' @param adduct_match Character scalar, possibly semicolon-separated
#' @return Logical
is_primary_ion <- function(ion_type, adduct_match) {

  raw_am <- if (!is.null(adduct_match) && !is.na(adduct_match) &&
                adduct_match != "") adduct_match else ""
  raw_it <- if (!is.null(ion_type)    && !is.na(ion_type)    &&
                ion_type    != "") ion_type    else ""

  split_strip <- function(s) {
    if (s == "") return(character(0))
    tokens <- trimws(strsplit(s, ";")[[1]])
    tokens <- tokens[tokens != ""]
    vapply(tokens, .strip_one_adduct, character(1))
  }

  am_tokens <- split_strip(raw_am)
  it_tokens <- split_strip(raw_it)
  tokens <- if (length(am_tokens) > 0) am_tokens else it_tokens
  if (length(tokens) == 0) return(FALSE)

  monomer_patterns <- c(
    "^\\[M\\+H\\]\\+$",
    "^\\[M\\+Na\\]\\+$",
    "^\\[M\\+NH4\\]\\+$",
    "^\\[M\\+K\\]\\+$"
  )
  # Matches any singly-charged NL ion [M+X-Y]+ (any loss formula).
  # Matches any singly-charged NL: [M+X-Y]+ (perl=TRUE handles ] in char classes)
  nl_pattern <- "^\\[M[^-\\]]*-.+\\]\\+$"

  any(vapply(tokens, function(tok) {
    if (!is_singly_charged(tok)) return(FALSE)
    if (any(sapply(monomer_patterns, grepl, x = tok,
                   ignore.case = TRUE))) return(TRUE)
    if (grepl(nl_pattern, tok, ignore.case = TRUE, perl = TRUE)) return(TRUE)
    if (grepl("^\\[(2|3|4)M", tok, ignore.case = TRUE)) return(TRUE)
    FALSE
  }, logical(1)))
}


#' Build the label string for a single query peak
#'
#' @param mz            Numeric scalar
#' @param ion_type      Character scalar
#' @param adduct_match  Character scalar, possibly semicolon-separated
#' @param mz_digits     Integer
#' @param mz_label_mode One of \code{"all"}, \code{"primary"}, \code{"adduct"}
#' @return Character scalar label, or \code{""}
build_query_label <- function(mz, ion_type, adduct_match,
                              mz_digits, mz_label_mode) {

  if (mz_label_mode == "all")
    return(format(round(mz, mz_digits), nsmall = mz_digits))

  clean_adduct <- function(s) {
    if (is.null(s) || is.na(s) || s == "") return("")
    tokens   <- trimws(strsplit(s, ";")[[1]])
    tokens   <- tokens[tokens != ""]
    # Require properly bracketed adduct notation — reject bare strings like "M+H"
    tokens   <- tokens[grepl("^\\[", tokens)]
    if (length(tokens) == 0) return("")
    stripped <- vapply(tokens, .strip_one_adduct, character(1))
    if (mz_label_mode == "primary")
      stripped <- stripped[vapply(stripped, is_singly_charged, logical(1))]
    if (length(stripped) == 0) return("")
    paste(stripped, collapse = "; ")
  }

  am <- clean_adduct(adduct_match)
  if (!is.null(am) && !is.na(am) && am != "") return(am)

  it <- clean_adduct(ion_type)
  if (!is.null(it) && !is.na(it) && it != "") return(it)

  ""
}


#' Draw one mirror-plot panel (query up, reference down)
#'
#' @param df           Data frame from \code{build_mirror_df()}
#' @param xlim_range   Length-2 numeric or \code{NULL}
#' @param title        Panel title
#' @param precursor_mz Numeric scalar or \code{NULL}
#' @param mz_digits    Integer
#' @param mz_label_mode One of \code{"all"}, \code{"primary"}, \code{"adduct"}
#' @return A \code{ggplot} object
plot_mirror <- function(df,
                        xlim_range    = NULL,
                        title         = "",
                        precursor_mz  = NULL,
                        mz_digits     = 1,
                        mz_label_mode = "all") {

  scale_method <- attr(df, "scale_method")
  if (is.null(scale_method)) scale_method <- "none"

  y_title <- switch(scale_method,
                    none  = "Relative intensity (%)",
                    sqrt  = "sqrt(Relative intensity (%))",
                    log10 = "Relative intensity (%, log10 scale)"
  )

  if (is.null(xlim_range)) {
    xr         <- diff(range(df$mz, na.rm = TRUE))
    xlim_range <- range(df$mz, na.rm = TRUE) + c(-1, 1) * xr * 0.05
  }

  ymax <- max(abs(df$intensity), na.rm = TRUE) * 1.8

  connector_df <- df %>%
    dplyr::filter(category == "Matched", !is.na(match_id)) %>%
    dplyr::group_by(match_id) %>%
    dplyr::summarise(x = min(mz), xend = max(mz), .groups = "drop") %>%
    dplyr::filter(x != xend)

  p <- ggplot(df, aes(x = mz, y = intensity, colour = category)) +
    geom_segment(aes(xend = mz, yend = 0), linewidth = 0.4, alpha = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_colour_manual(values = c("Matched"   = "#1F78B4",
                                   "Unmatched" = "black")) +
    coord_cartesian(xlim = xlim_range, ylim = c(-ymax, ymax),
                    expand = FALSE, clip = "off") +
    labs(title = title, x = "m/z", y = y_title) +
    theme_bw(base_size = 11) +
    theme(panel.grid      = element_blank(),
          legend.position = "none",
          plot.title      = element_text(face = "bold"),
          axis.title.y    = element_text(margin = margin(r = 8)),
          plot.margin = margin(50, 18, 8, 8))

  if (scale_method == "log10") {
    ax <- build_log10_yaxis(ymax)
    p  <- p + scale_y_continuous(breaks = ax$breaks, labels = ax$labels)
  }

  if (nrow(connector_df) > 0)
    p <- p + geom_segment(
      data        = connector_df,
      aes(x = x, xend = xend, y = 0, yend = 0),
      inherit.aes = FALSE, colour = "grey70", linewidth = 0.3)

  # Build query label data frame without rowwise() to avoid scoping issues
  query_base_df <- df %>%
    dplyr::filter(source == "Query", category == "Matched")

  if (nrow(query_base_df) > 0) {
    # Vectorised label and show_label computation
    labels_vec <- mapply(
      build_query_label,
      mz           = query_base_df$mz,
      ion_type     = query_base_df$ion_type,
      adduct_match = query_base_df$adduct_match,
      MoreArgs     = list(mz_digits = mz_digits, mz_label_mode = mz_label_mode),
      SIMPLIFY     = TRUE, USE.NAMES = FALSE
    )

    show_vec <- if (mz_label_mode == "all") {
      rep(TRUE, nrow(query_base_df))
    } else if (mz_label_mode == "primary") {
      mapply(is_primary_ion,
             ion_type     = query_base_df$ion_type,
             adduct_match = query_base_df$adduct_match,
             SIMPLIFY = TRUE, USE.NAMES = FALSE)
    } else if (mz_label_mode == "adduct") {
      vapply(query_base_df$adduct_match, function(am) {
        !is.null(am) && !is.na(am) && strip_adduct_mass(am) != ""
      }, logical(1), USE.NAMES = FALSE)
    } else {
      rep(TRUE, nrow(query_base_df))
    }

    query_label_df <- query_base_df %>%
      dplyr::mutate(
        label      = labels_vec,
        show_label = show_vec,
        y_lab      = intensity + 0.04 * ymax
      ) %>%
      dplyr::filter(show_label, !is.na(label), label != "") %>%
      # Deduplicate: when multiple peaks share the same label, keep most intense
      dplyr::group_by(label) %>%
      dplyr::slice_max(order_by = intensity, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else {
    query_label_df <- query_base_df %>%
      dplyr::mutate(label = character(0), show_label = logical(0),
                    y_lab = numeric(0))
  }

  if (nrow(query_label_df) > 0)

    # Query labels — horizontal text, long leader lines
    p <- p + ggrepel::geom_text_repel(
      data               = query_label_df,
      aes(x = mz, y = intensity, label = label),
      inherit.aes        = FALSE,
      size               = 2.6,
      colour             = "#1F78B4",
      angle              = 0,                        # ← horizontal (default)
      direction          = "both",                   # ← free to spread x and y
      segment.color      = "#1F78B4",
      segment.size       = 0.3,
      segment.alpha      = 0.6,
      segment.linetype   = 1,
      segment.curvature  = 0,                        # ← straight lines
      point.size         = 0.1,
      point.padding      = unit(0.2, "lines"),
      min.segment.length = 0,                        # ← always draw segment
      box.padding        = 0.5,                      # ← more space around label box
      max.overlaps       = Inf,
      force              = 3.0,                      # ← stronger repulsion between labels
      force_pull         = 0.5,                      # ← weaker pull back to peak (longer segments)
      nudge_y            = ymax * 0.3,               # ← push labels upward initially
      seed               = 42,
      ylim               = c(0, ymax * 1.8),         # ← tall ceiling for label placement
      xlim               = c(xlim_range[1], xlim_range[2])
    )

  ref_label_df <- if (mz_label_mode == "all") {
    df %>%
      dplyr::filter(source == "Reference", category == "Matched") %>%
      dplyr::mutate(
        label = format(round(mz, mz_digits), nsmall = mz_digits),
        y_lab = pmax(intensity - 0.04 * ymax, -ymax * 0.94)
      )
  } else {
    NULL
  }

  if (!is.null(ref_label_df) && nrow(ref_label_df) > 0)
    p <- p + ggrepel::geom_text_repel(
      data               = ref_label_df,
      aes(x = mz, y = intensity, label = label),
      inherit.aes        = FALSE,
      size               = 2.6,
      colour             = "black",
      angle              = 90,
      direction          = "y",
      hjust              = 1,
      segment.color      = "black",
      segment.size       = 0.35,
      segment.alpha      = 0.85,
      segment.linetype   = 1,
      point.size         = 0.1,
      point.padding      = unit(0.4, "lines"),
      min.segment.length = 0,
      box.padding        = 0.15,
      max.overlaps       = Inf,
      force              = 1.2,
      force_pull         = 1.5,
      seed               = 123,
      nudge_y            = -ymax * 0.3,
      ylim               = c(-ymax * 1.8, 0),
      xlim               = c(xlim_range[1], xlim_range[2])
    )
  if (!is.null(precursor_mz) && !is.na(precursor_mz)) {
    x_star <- precursor_mz + diff(xlim_range) * 0.008
    p <- p +
      ggplot2::annotate("segment",
                        x = precursor_mz, xend = precursor_mz, y = 0, yend = -ymax,
                        colour = "red", linetype = "dashed", linewidth = 0.6) +
      ggplot2::annotate("text",
                        x = x_star, y = -ymax,
                        label = "*", colour = "red", size = 7, fontface = "bold")
  }

  p
}


#' Batch-export mirror plots for all compounds in a direct_match list
#'
#' @param input_list_ms1       List; \code{direct_match_StM14_res_MS1}
#' @param input_list_ms2       List; \code{direct_match_StM14_res_MS2}
#' @param annot_res_ms1        Data frame; \code{StM14_annotate_res_ms1}
#'   (only \code{annotation_phase == "Phase1"} rows are used)
#' @param annot_res_ms2        Data frame; \code{StM14_annotate_res_ms2}
#'   (only \code{annotation_phase == "Phase1"} rows are used)
#' @param out_dir              Output directory (created if absent)
#' @param scale_method         Intensity scale
#' @param mz_digits            Decimal places for m/z labels
#' @param padding              Fractional x-axis padding
#' @param precursor_extra_frac Extra right-side padding for precursor marker
#' @param shared_xlim          Logical; share x-axis across MS1 and MS2 panels
#' @param mz_label_mode        One of \code{"all"}, \code{"primary"},
#'   \code{"adduct"}
#' @param width                Plot width in inches
#' @param height               Plot height in inches
#' @param dpi                  Resolution for PNG output
#' @return Invisibly, a list of per-compound output paths
export_mirror_plots <- function(input_list_ms1,
                                input_list_ms2,
                                annot_res_ms1,
                                annot_res_ms2,
                                out_dir              = "mirror_plots",
                                scale_method         = "log10",
                                mz_digits            = 1,
                                padding              = 0.05,
                                precursor_extra_frac = 0.04,
                                shared_xlim          = TRUE,
                                mz_label_mode        = "all",
                                width                = 6.5,
                                height               = 9,
                                dpi                  = 300) {

  stopifnot(mz_label_mode %in% c("all", "primary", "adduct"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  results <- vector("list", length(input_list_ms1))

  for (i in seq_along(input_list_ms1)) {
    ms1 <- input_list_ms1[[i]]
    if (is.null(ms1$Compound_name) || length(ms1$Compound_name) == 0) {
      message("Skipping [[", i, "]]: no Compound_name"); next
    }

    cname   <- as.character(ms1$Compound_name[1])
    ms2_idx <- find_by_name(cname, input_list_ms2)
    if (is.na(ms2_idx)) {
      message("Skipping ", cname, ": no MS2 entry found"); next
    }
    ms2 <- input_list_ms2[[ms2_idx]]

    p <- tryCatch(
      plot_one_compound(
        ms1_entry            = ms1,
        ms2_entry            = ms2,
        annot_res_ms1        = annot_res_ms1,
        annot_res_ms2        = annot_res_ms2,
        scale_method         = scale_method,
        mz_digits            = mz_digits,
        padding              = padding,
        precursor_extra_frac = precursor_extra_frac,
        shared_xlim          = shared_xlim,
        mz_label_mode        = mz_label_mode
      ),
      error = function(e) {
        message("Error plotting ", cname, ": ", e$message); NULL
      }
    )
    if (is.null(p)) next

    base     <- paste0(sanitize_filename(cname), "_", scale_method,
                       "_", mz_label_mode, "_mirror")
    file_png <- file.path(out_dir, paste0(base, ".png"))
    file_pdf <- file.path(out_dir, paste0(base, ".pdf"))

    ggsave(file_png, p, width = width, height = height,
           dpi = dpi, bg = "white")
    ggsave(file_pdf, p, width = width, height = height,
           device = cairo_pdf, bg = "white")

    results[[i]] <- list(compound_name = cname, ms2_index = ms2_idx,
                         png = file_png, pdf = file_pdf)
    message("Exported: ", cname, " [", mz_label_mode, "]")
    rm(p); gc()
  }
  invisible(results)
}
