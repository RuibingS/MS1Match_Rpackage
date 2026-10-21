
setGeneric("annotateLibMatch",
           function(object, cutoff_MH = FALSE, ...) standardGeneric("annotateLibMatch"))

# ── Method ────────────────────────────────────────────────────────────────────

#' Annotate features from spectral library match results
#'
#' @description
#' Two-phase annotation dispatched on a
#' \code{\link{PeakAnnotationFromLibraryMatch-class}} object.
#'
#' \strong{Phase 1} — \code{annotate_librarymatch_phase1_cpp()}:
#' \itemize{
#'   \item Runs \code{Check_Iso_Charge} \emph{once} on the union of all
#'     library-matched features and the full pool.  Produces a global isotope
#'     label table shared by both phases — isotope is never re-computed.
#'   \item For each library-match group: adduct match (Tier 1a), isotope
#'     lookup (Tier 1b), pairwise NL within group (Tier 1c), precursor/MS2
#'     match (Tier 1d), unexplained (Tier 2).
#'   \item RT is used automatically when available; if RT values are
#'     \code{NA} the RT gate simply does not fire.
#'   \item Returns annotated rows plus \code{confirmed_M} (confirmed neutral
#'     masses with reference RT) and \code{iso_labels} (global isotope table)
#'     for Phase 2.
#' }
#'
#' \strong{Phase 2} — \code{annotate_librarymatch_phase2_cpp()}
#' (only when \code{FT} slot is non-empty):
#' \itemize{
#'   \item For each confirmed M from Phase 1, scans pool features not already
#'     annotated within \code{rt_exp_threshold} s of the compound reference RT
#'     and \code{ppm_exp} ppm of any adduct mass (Tier 1a_exp).
#'   \item Looks up iso labels from the Phase 1 global table (Tier 1b_exp).
#'   \item Runs pairwise NL within the same-M expansion group (Tier 1c_exp).
#' }
#'
#' @param object A \code{\link{PeakAnnotationFromLibraryMatch-class}} object.
#' @param cutoff_MH logical.  If \code{TRUE}, adds a \code{cutoff_to_MH}
#'   column to the result that classifies each \code{Tier2_Unexplained} feature
#'   as \code{"above_MH"} or \code{"below_MH"} relative to the theoretical
#'   \eqn{[M+H]^+} / \eqn{[M-H]^-} ion of the matched compound.  All other
#'   tiers and rows with unknown \code{Compound_M} receive \code{NA}.
#'   Default \code{FALSE} — the column is not added unless explicitly requested.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns: \code{feature_name}, \code{mz},
#'   \code{rt}, \code{intensity}, \code{Compound_name}, \code{Lib_group},
#'   \code{Ion_type}, \code{Compound_M}, \code{PI_match},
#'   \code{grp_anno_match}, \code{grp_anno_ppm}, \code{grp_anno_mda},
#'   \code{grp_anno_confidence}, \code{exp_anno_match}, \code{exp_anno_ppm},
#'   \code{exp_anno_mda}, \code{exp_anno_confidence},
#'   \code{NL_match}, \code{iso_annotation}, \code{reference_tier},
#'   \code{annotation_phase}, \code{Match_peak_score}, \code{inner_cosine},
#'   \code{outer_cosine}.  When \code{cutoff_MH = TRUE}, the additional column
#'   \code{cutoff_to_MH} is appended (see \code{cutoff_MH} parameter above).
#'
#' @seealso \code{\link{createPeakAnnotationFromLibraryMatch}()}
#' @export
setMethod("annotateLibMatch", "PeakAnnotationFromLibraryMatch",
          function(object, cutoff_MH = FALSE, ...) {

  # ── Phase 1 ─────────────────────────────────────────────────────────────────
  p1 <- annotate_librarymatch_phase1_cpp(
    library_matches    = object@library_matches,
    FT            = object@FT,
    adducts_df         = object@adducts_df,
    neutral_loss_table = object@neutral_loss_table,
    polarity           = object@polarity,
    mz_threshold       = object@mz_threshold,
    ppm                = object@ppm,
    nl_ppm             = object@nl_ppm,
    rt_threshold       = object@rt_threshold,
    pi_rt_threshold    = object@pi_rt_threshold,
    rt_exp_threshold   = object@rt_exp_threshold,
    iso_threshold      = object@iso_threshold,
    max_charge         = object@max_charge,
    max_isotope        = object@max_isotope
  )

  # ── Phase 2 (skip when no ft provided) ────────────────────────────────────
  if (nrow(object@FT) == 0L) {
    result <- .rows_to_df(p1$rows)
    if (cutoff_MH) result <- .add_cutoff_to_MH(result, polarity = object@polarity)
    return(result)
  }

  #  check if any Phase 2 features were already in Phase 1 groups.

  p2_candidates <- FT_fn <- if (nrow(object@FT) > 0)
    object@FT[["feature_name"]] else character(0L)
  overlap <- intersect(trimws(p1$phase1_features), trimws(p2_candidates))
  if (length(overlap) > 0) {
    message(paste0(
      "[annotateLibMatch] ", length(overlap),
      "features share names with Phase 1 features after trimming — ",
      "they will be correctly skipped in Phase 2. ",
      "If unexpected, check feature name consistency between ",
      "library_matches and FT."
    ))
  }

  p2 <- annotate_librarymatch_phase2_cpp(
    FT            = object@FT,
    confirmed_M_df     = p1$confirmed_M,
    iso_labels_df      = p1$iso_labels,
    phase1_features    = p1$phase1_features,
    adducts_df         = object@adducts_df,
    neutral_loss_table = object@neutral_loss_table,
    polarity           = object@polarity,
    mz_threshold       = object@mz_threshold,
    ppm_exp            = object@ppm_exp,
    nl_ppm             = object@nl_ppm,
    rt_exp_threshold   = object@rt_exp_threshold
  )

  result <- .rows_to_df(c(p1$rows, p2))
  if (cutoff_MH) result <- .add_cutoff_to_MH(result, polarity = object@polarity)
  result
})


# ── Internal helper ───────────────────────────────────────────────────────────
.rows_to_df <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# ── Tier2_Unexplained: [M+H]+ / [M-H]- cutoff column ─────────────────────────
#
# Called only when cutoff_MH = TRUE.
#
# Adds column  cutoff_to_MH  to the result data frame.
# For Tier2_Unexplained rows, classifies observed feature mz relative to the
# theoretical precursor ion of the matched compound:
#
#   Positive mode  ([M+H]+) : precursor_mz = Compound_M + 1.007276
#   Negative mode  ([M-H]-) : precursor_mz = Compound_M - 1.007276
#
# Accepted polarity strings (case-insensitive):
#   positive : "pos", "positive", "+"
#   negative : "neg", "negative", "-"
#   Any other value raises an error.
#
# Two categories:
#   "above_MH"  — feature mz  >  precursor_mz  (higher-mz region, e.g. adducts,
#                                                multimers, in-source dimers)
#   "below_MH"  — feature mz  <= precursor_mz  (lower-mz region, e.g. ISF
#                                                fragments, losses)
#
# NA is assigned for:
#   - all tiers other than Tier2_Unexplained
#   - rows where Compound_M is NA (M unknown)
.add_cutoff_to_MH <- function(result, polarity) {

  if (nrow(result) == 0L || !"reference_tier" %in% colnames(result))
    return(result)

  # ── normalise polarity string ──────────────────────────────────────────────
  pol <- tolower(trimws(polarity))
  if (pol %in% c("pos", "positive", "+")) {
    pol_norm <- "pos"
  } else if (pol %in% c("neg", "negative", "-")) {
    pol_norm <- "neg"
  } else {
    stop(
      "[.add_cutoff_to_MH] unrecognised polarity value: '", polarity,
      "'.  Expected one of: \"pos\", \"positive\", \"+\", ",
      "\"neg\", \"negative\", \"-\"."
    )
  }

  # ── theoretical precursor mz ──────────────────────────────────────────────
  proton <- 1.007276
  precursor_mz <- if (pol_norm == "pos") {
    result$Compound_M + proton   # [M+H]+
  } else {
    result$Compound_M - proton   # [M-H]-
  }

  # ── classify Tier2_Unexplained rows only ──────────────────────────────────
  is_tier2  <- result$reference_tier == "Tier2_Unexplained"
  m_known   <- !is.na(result$Compound_M)
  eligible  <- is_tier2 & m_known

  cutoff <- rep(NA_character_, nrow(result))
  cutoff[eligible] <- ifelse(
    result$mz[eligible] > precursor_mz[eligible],
    "above_MH",
    "below_MH"
  )

  result$cutoff_to_MH <- cutoff
  result
}
