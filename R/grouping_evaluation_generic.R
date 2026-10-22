# =============================================================================
# Tool-agnostic evaluation of feature grouping against the MS1Match
# compound-centric reference set
# =============================================================================
#
# Requirement: MS1Match and the grouping tool must have been run on the SAME
# feature table. Results from tools that built their own feature table cannot
# be evaluated.
# Works for any grouping tool (CAMERA, RAMClustR, ...) as long as the
# tool output can be reduced to one row per feature with:
#   feature_name : the feature name from the feature table MS1Match was run on
#                  (xcms "FT0001"; MZmine id 1 -> "FT0001" via the import)
#   group        : the tool's group / cluster id (NA = feature not grouped)
#   mz           : m/z of the feature (used to verify the feature table)
#   rt           : optional, retention time in seconds (used for rt_window)
#
# Status definitions (see manuscript Table S4). For each compound, the
# Tier 1a [M+H]+ feature is the anchor; the group the tool assigns to the
# anchor is the compound's expected group. Every reference compound has an
# anchor by construction.
#
# Reference features get the FIRST status that applies, in this order:
#   1. Ref-unassigned (no-anchor)  the tool did not group the compound's
#                                  anchor; none of its features can be
#                                  evaluated
#   2. Ref-unassigned              the anchor is grouped, this feature is not
#   3. Ref-confirmed               same group as the anchor
#   4. Ref-conflicting             different group than the anchor
#
# Features outside the reference set:
#   Ref-unmatched, co-grouped      in a group that contains a compound anchor
#   Ref-unmatched, unrelated       in a group without an anchor, or in no group
#
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Collapse MS1Match reference_tier strings to the tier labels used in the benchmark tables
.tier_label <- function(tier) {
  tier   <- as.character(tier)
  is_exp <- grepl("exp", tier, ignore.case = TRUE)
  is_nl  <- grepl("NeutralLoss|^Tier1d", tier, ignore.case = TRUE)
  dplyr::case_when(
    is_nl & is_exp                                  ~ "Tier1d_exp",
    is_nl                                           ~ "Tier1d",
    grepl("^Tier1a", tier, ignore.case = TRUE) & is_exp ~ "Tier1a_exp",
    grepl("^Tier1a", tier, ignore.case = TRUE)          ~ "Tier1a",
    grepl("^Tier1b", tier, ignore.case = TRUE) & is_exp ~ "Tier1b_exp",
    grepl("^Tier1b", tier, ignore.case = TRUE)          ~ "Tier1b",
    grepl("^Tier1c", tier, ignore.case = TRUE) & is_exp ~ "Tier1c_exp",
    grepl("^Tier1c", tier, ignore.case = TRUE)          ~ "Tier1c",
    grepl("^Tier2",  tier, ignore.case = TRUE)          ~ "Tier2",
    TRUE                                                ~ NA_character_
  )
}

.tier_levels <- c("Tier1a", "Tier1a_exp", "Tier1b", "Tier1b_exp",
                  "Tier1c", "Tier1c_exp", "Tier1d", "Tier1d_exp", "Tier2")

.ref_status_levels <- c("Ref-unassigned (no-anchor)", "Ref-unassigned",
                        "Ref-confirmed", "Ref-conflicting",
                        "Not evaluable")

.all_status_levels <- c(.ref_status_levels,
                        "Ref-unmatched, co-grouped", "Ref-unmatched, unrelated")


# -----------------------------------------------------------------------------
# Preparing the grouping table
# -----------------------------------------------------------------------------

#' Standardize a grouping tool's output for evaluation
#'
#' @description
#' Reduces a feature-level grouping result to the columns used by
#' \code{\link{evaluate_grouping}}: \code{feature_name}, \code{group},
#' \code{mz} and \code{rt}. The grouping tool must have been run on the same feature table
#' as MS1Match; this function only renames columns and cannot map features
#' between different feature tables.
#'
#' @param x \code{data.frame} with one row per feature.
#' @param group_col Name of the column holding the group / cluster id.
#' @param feature_col Name of the column holding the feature id (same ids
#'   as in the feature table MS1Match was run on). Plain numbers (1, 2, ...)
#'   are converted to \code{FT0001}, \code{FT0002}, ..., as in
#'   \code{\link{read_feature_table}}.
#'   Default \code{"feature_name"}.
#' @param mz_col Name of the m/z column. Needed to verify that the tool used
#'   the same feature table as MS1Match.
#' @param rt_col Optional name of the retention-time column.
#' @param rt_unit Unit of \code{rt_col}: \code{"sec"} (default) or
#'   \code{"min"}. Minutes are converted to seconds.
#' @param unassigned Values of \code{group_col} that mean "not grouped"
#'   (in addition to \code{NA}). Default \code{c("0", "")}.
#'
#' @return \code{data.frame} with columns \code{feature_name}, \code{group}
#'   (character, \code{NA} = unassigned), \code{mz} and \code{rt}.
#' @export
prepare_grouping_table <- function(x, group_col, feature_col = "feature_name",
                                   mz_col = NULL, rt_col = NULL,
                                   rt_unit = c("sec", "min"),
                                   unassigned = c("0", "")) {
  stopifnot(is.data.frame(x))
  rt_unit <- match.arg(rt_unit)
  if (!feature_col %in% names(x))
    stop("prepare_grouping_table: no '", feature_col, "' column. If your tool ",
         "exported no feature ids, add them with assign_feature_names() first.")
  missing <- setdiff(c(feature_col, group_col, mz_col, rt_col), names(x))
  if (length(missing) > 0L)
    stop("prepare_grouping_table: column(s) not found: ",
         paste(missing, collapse = ", "),
         "\n  Available columns: ", paste(names(x), collapse = ", "))

  grp <- trimws(as.character(x[[group_col]]))
  grp[grp %in% unassigned] <- NA_character_

  out <- data.frame(
    feature_name = .as_feature_name(x[[feature_col]]),
    group        = grp,
    mz           = if (!is.null(mz_col)) as.numeric(x[[mz_col]]) else NA_real_,
    rt           = if (!is.null(rt_col))
                     as.numeric(x[[rt_col]]) * (if (rt_unit == "min") 60 else 1)
                   else NA_real_,
    stringsAsFactors = FALSE
  )

  dup <- duplicated(out$feature_name)
  if (any(dup)) {
    warning("prepare_grouping_table: ", sum(dup), " duplicated feature_name ",
            "value(s); keeping the first occurrence.")
    out <- out[!dup, , drop = FALSE]
  }
  out
}


#' Grouping table from a CAMERA peak list
#'
#' @param peaklist Output of \code{CAMERA::getPeaklist()} (or a saved copy).
#'   Must contain \code{pcgroup}.
#' @param feature_names Optional character vector of xcms feature ids, one per
#'   row of \code{peaklist}. If \code{NULL} and \code{peaklist} has no
#'   \code{feature_name} column, ids are generated in xcms style
#'   (\code{FT0001}, ...), assuming the rows are in xcms feature order.
#'
#' @return See \code{\link{prepare_grouping_table}}.
#' @export
camera_to_grouping <- function(peaklist, feature_names = NULL) {
  stopifnot(is.data.frame(peaklist), "pcgroup" %in% names(peaklist))

  if (!is.null(feature_names)) {
    stopifnot(length(feature_names) == nrow(peaklist))
    peaklist$feature_name <- feature_names
  } else if (!"feature_name" %in% names(peaklist)) {
    n <- nrow(peaklist)
    peaklist$feature_name <- sprintf(
      paste0("FT%0", ceiling(log10(n + 1L)), "d"), seq_len(n))
    message("camera_to_grouping: no feature_name column; generated xcms-style ",
            "ids assuming peaklist rows follow xcms feature order.")
  }

  prepare_grouping_table(
    peaklist,
    group_col = "pcgroup",
    mz_col    = if ("mz" %in% names(peaklist)) "mz" else NULL,
    rt_col    = if ("rt" %in% names(peaklist)) "rt" else NULL
  )
}


#' Grouping table from a RAMClustR result
#'
#' @param RC Output of \code{RAMClustR::ramclustR()}, run on the same feature
#'   table as MS1Match (e.g. the xcms object). \code{RC$frt} then has the same
#'   RT unit as that table, i.e. seconds, and is used as it is.
#'
#' @return See \code{\link{prepare_grouping_table}}. Features with
#'   \code{featclus = 0} (not clustered) get \code{group = NA}.
#' @export
ramclustr_to_grouping <- function(RC) {
  needed  <- c("featnames", "featclus", "cmpd", "fmz", "frt")
  missing <- setdiff(needed, names(RC))
  if (length(missing) > 0L)
    stop("ramclustr_to_grouping: RC object lacks: ", paste(missing, collapse = ", "))

  cl  <- RC$featclus
  grp <- rep(NA_character_, length(cl))
  ok  <- !is.na(cl) & cl > 0
  grp[ok] <- RC$cmpd[cl[ok]]


  data.frame(
    feature_name = .as_feature_name(RC$featnames),
    group        = grp,
    mz           = as.numeric(RC$fmz),
    rt           = as.numeric(RC$frt),
    stringsAsFactors = FALSE
  )
}


#' Grouping table from an MZmine feature-table export
#'
#' @description
#' Converts a grouping stored in an MZmine 3 full feature table (e.g. ion
#' identity networking) to the format used by \code{\link{evaluate_grouping}}.
#' Applies the same conversions as \code{\link{read_mzmine_feature_table}}:
#' MZmine \code{id} becomes \code{feature_name = sprintf("FT\%04d", id)} and RT
#' is converted from minutes to seconds.
#'
#' @param x \code{data.frame} read from the MZmine CSV
#'   (\code{read.csv(file, check.names = FALSE)}). Must contain \code{id}
#'   and \code{group_col}.
#' @param group_col Name of the column holding the group id (check your
#'   export, e.g. the ion identity network id column).
#' @param keep_features Optional character vector of \code{feature_name}
#'   values to keep, typically \code{read_mzmine_feature_table(file)$feature_name},
#'   so that features removed by the import's detection filter are dropped.
#' @param mz_col Name of the m/z column. Default \code{"mz"}.
#' @param rt_col Name of the RT column in minutes. Default \code{"rt"}; set to
#'   \code{NULL} if absent.
#' @param unassigned Values of \code{group_col} meaning "not grouped".
#'   Default \code{""}: MZmine numbers its groups from 0, so 0 is a real group.
#'
#' @return See \code{\link{prepare_grouping_table}}.
#' @export
mzmine_to_grouping <- function(x, group_col, keep_features = NULL,
                               mz_col = "mz", rt_col = "rt",
                               unassigned = "") {
  stopifnot(is.data.frame(x))
  if (!"id" %in% names(x))
    stop("mzmine_to_grouping: no 'id' column; is this an MZmine feature table export?")

  x <- as.data.frame(x)
  x$feature_name <- sprintf("FT%04d", as.integer(x$id))

  if (!is.null(rt_col)) {
    if (!rt_col %in% names(x))
      stop("mzmine_to_grouping: RT column '", rt_col, "' not found.")
    x$.rt_s <- as.numeric(x[[rt_col]]) * 60
  }

  if (!is.null(keep_features)) {
    n_before <- nrow(x)
    x <- x[x$feature_name %in% keep_features, , drop = FALSE]
    message("mzmine_to_grouping: kept ", nrow(x), " of ", n_before,
            " features present in the imported feature table.")
  }

  prepare_grouping_table(
    x,
    group_col  = group_col,
    mz_col     = mz_col,
    rt_col     = if (!is.null(rt_col)) ".rt_s" else NULL,
    unassigned = unassigned
  )
}


# -----------------------------------------------------------------------------
# Feature ids for grouping results that have none
# -----------------------------------------------------------------------------

#' Give a grouping result the feature ids of the MS1Match feature table
#'
#' @description
#' For a grouping result without feature ids. Each row gets the
#' \code{feature_name} of the feature in \code{feature_table} with the same
#' m/z and RT. Because both come from the same feature table, the values are
#' identical; the tolerances only absorb rounding. Rows that match no feature,
#' or more than one, get \code{NA}; two rows matching the same feature both
#' get \code{NA}.
#'
#' @param x Grouping result, one row per feature.
#' @param feature_table The feature table MS1Match was run on
#'   (\code{feature_name}, \code{mzmed}/\code{mz}, \code{rtmed}/\code{rt} in
#'   seconds).
#' @param mz_col,rt_col Names of the m/z and RT columns in \code{x}.
#' @param rt_unit Unit of \code{rt_col}: \code{"sec"} (default) or \code{"min"}.
#' @param mz_tol,rt_tol Rounding tolerances (Da, seconds). Defaults 0.001, 1.
#'
#' @return \code{x} with a \code{feature_name} column and \code{rt_col}
#'   converted to seconds. A message reports how many rows were matched.
#' @export
assign_feature_names <- function(x, feature_table, mz_col = "mz", rt_col = "rt",
                                 rt_unit = c("sec", "min"),
                                 mz_tol = 0.001, rt_tol = 1) {
  stopifnot(is.data.frame(x), is.data.frame(feature_table))
  rt_unit <- match.arg(rt_unit)
  if (!all(c(mz_col, rt_col) %in% names(x)))
    stop("assign_feature_names: x needs columns '", mz_col, "' and '", rt_col, "'.")
  ft_mz_col <- intersect(c("mzmed", "mz"), names(feature_table))[1L]
  ft_rt_col <- intersect(c("rtmed", "rt"), names(feature_table))[1L]
  if (!"feature_name" %in% names(feature_table) || is.na(ft_mz_col) || is.na(ft_rt_col))
    stop("assign_feature_names: feature_table needs feature_name, mzmed and rtmed.")

  x <- as.data.frame(x)
  x[[rt_col]] <- as.numeric(x[[rt_col]]) * (if (rt_unit == "min") 60 else 1)
  .warn_rt_unit(x[[rt_col]], rt_unit, "assign_feature_names")

  ft_id <- as.character(feature_table$feature_name)
  ft_mz <- as.numeric(feature_table[[ft_mz_col]])
  ft_rt <- as.numeric(feature_table[[ft_rt_col]])
  x_mz  <- as.numeric(x[[mz_col]])
  x_rt  <- x[[rt_col]]

  ids <- vapply(seq_len(nrow(x)), function(i) {
    hit <- which(abs(ft_mz - x_mz[i]) <= mz_tol & abs(ft_rt - x_rt[i]) <= rt_tol)
    if (length(hit) == 1L) ft_id[hit] else NA_character_
  }, character(1))

  dup <- !is.na(ids) & (duplicated(ids) | duplicated(ids, fromLast = TRUE))
  ids[dup] <- NA_character_

  n_ok <- sum(!is.na(ids))
  message(sprintf(paste0("assign_feature_names: %d of %d rows matched a feature ",
                         "(%d unmatched or ambiguous)."),
                  n_ok, nrow(x), nrow(x) - n_ok))
  if (n_ok < nrow(x))
    message("  Unmatched rows are not from this feature table; check that the ",
            "tool ran on it and that rt_unit is correct.")

  x$feature_name <- ids
  x
}


# -----------------------------------------------------------------------------
# Feature-name check
# -----------------------------------------------------------------------------

#' Check that MS1Match and a grouping tool used the same feature table
#'
#' @description
#' A yes/no check. The answer is \code{TRUE} only if
#' \enumerate{
#'   \item every MS1Match reference feature exists in \code{feature_table}
#'     with the same m/z, and
#'   \item every feature id shared by the grouping table and
#'     \code{feature_table} has the same m/z (and, where available, the same
#'     RT), and at least one id is shared.
#' }
#' Feature ids alone are not enough: separate xcms runs all name their
#' features \code{FT0001}, \code{FT0002}, ...
#'
#' One table may be a subset of the other, e.g. MS1Match run on an
#' RT-truncated feature table and the grouping tool on the full one. Ids
#' present in only one of them do not fail the check; they are counted and
#' excluded by \code{\link{evaluate_grouping}}.
#'
#' @param anno_result Output of \code{join_annotation_results()}.
#' @param grouping Output of \code{\link{prepare_grouping_table}} or a
#'   converter; must contain \code{mz}.
#' @param feature_table The feature table MS1Match was run on (e.g.
#'   \code{out$ft} for xcms, or \code{read_mzmine_feature_table()} for MZmine).
#'   Must contain \code{feature_name} and \code{mzmed} or \code{mz}.
#' @param mz_tol Allowed m/z difference in Da, only to absorb rounding
#'   (e.g. RAMClustR's \code{mzdec}). Default 0.001.
#' @param rt_tol Allowed RT difference in seconds. Default 1.
#'
#' @return Object of class \code{feature_table_check}: a list with
#'   \code{same_feature_table} (\code{TRUE}/\code{FALSE}), \code{summary}
#'   (counts) and \code{mismatches} (one row per offending feature and reason).
#' @export
check_feature_match <- function(anno_result, grouping, feature_table,
                                mz_tol = 0.001, rt_tol = 1) {

  if (missing(feature_table) || is.null(feature_table))
    stop("check_feature_match: feature_table is required (the table MS1Match was run on).")
  if (!"feature_name" %in% names(feature_table))
    stop("check_feature_match: feature_table needs a feature_name column.")
  ft_mz_col <- intersect(c("mzmed", "mz"), names(feature_table))[1L]
  ft_rt_col <- intersect(c("rtmed", "rt"), names(feature_table))[1L]
  if (is.na(ft_mz_col))
    stop("check_feature_match: feature_table needs an mzmed or mz column.")
  if (!"mz" %in% names(grouping) || all(is.na(grouping$mz)))
    stop("check_feature_match: grouping has no m/z values; use mz_col in ",
         "prepare_grouping_table() so the feature table can be verified.")

  ft <- data.frame(
    feature_name = as.character(feature_table$feature_name),
    ft_mz        = as.numeric(feature_table[[ft_mz_col]]),
    ft_rt        = if (!is.na(ft_rt_col)) as.numeric(feature_table[[ft_rt_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )

  .compare <- function(df, source) {
    m <- merge(df, ft, by = "feature_name", all.x = TRUE)
    reason <- ifelse(is.na(m$ft_mz), "id not in feature table",
              ifelse(!is.na(m$mz) & abs(m$mz - m$ft_mz) > mz_tol, "m/z differs",
              ifelse(!is.na(m$rt) & !is.na(m$ft_rt) & abs(m$rt - m$ft_rt) > rt_tol,
                     "RT differs", NA_character_)))
    data.frame(source = source, feature_name = m$feature_name,
               mz = m$mz, ft_mz = m$ft_mz, rt = m$rt, ft_rt = m$ft_rt,
               reason = reason, stringsAsFactors = FALSE)
  }

  # grouping features
  g <- unique(grouping[, c("feature_name", "mz", "rt")])
  g_res <- .compare(g, "grouping")

  # MS1Match reference features
  summ  <- anno_result$export_feature_compound_summary
  r_mz  <- intersect(c("mz", "mzmed"), names(summ))[1L]
  r_rt  <- intersect(c("rt", "rtmed"), names(summ))[1L]
  r <- unique(data.frame(
    feature_name = as.character(summ$feature_name),
    mz = if (!is.na(r_mz)) as.numeric(summ[[r_mz]]) else NA_real_,
    rt = if (!is.na(r_rt)) as.numeric(summ[[r_rt]]) else NA_real_,
    stringsAsFactors = FALSE))
  r <- r[!is.na(r$feature_name), , drop = FALSE]
  r_res <- .compare(r, "MS1Match reference")

  # "id not in feature table" fails only for MS1Match reference features;
  # for the grouping table it just means the tool's table is larger.
  g_fail <- g_res[g_res$reason %in% c("m/z differs", "RT differs"), , drop = FALSE]
  r_fail <- r_res[!is.na(r_res$reason), , drop = FALSE]
  mism <- rbind(g_fail, r_fail)
  rownames(mism) <- NULL

  n_shared <- sum(g$feature_name %in% ft$feature_name)

  summary <- data.frame(
    n_tool_features           = nrow(g),
    n_shared_features         = n_shared,
    n_tool_not_in_ft          = sum(g_res$reason %in% "id not in feature table"),
    n_tool_mz_mismatch        = sum(g_res$reason %in% "m/z differs"),
    n_tool_rt_mismatch        = sum(g_res$reason %in% "RT differs"),
    n_reference_features      = nrow(r),
    n_reference_not_in_ft     = sum(r_res$reason %in% "id not in feature table"),
    n_reference_mz_mismatch   = sum(r_res$reason %in% "m/z differs"),
    n_reference_not_in_tool   = sum(!r$feature_name %in% g$feature_name),
    n_ft_features             = nrow(ft),
    n_ft_not_in_tool          = sum(!ft$feature_name %in% g$feature_name)
  )

  structure(
    list(same_feature_table = nrow(mism) == 0L && n_shared > 0L,
         summary            = summary,
         mismatches         = mism),
    class = c("feature_table_check", "list")
  )
}


#' @export
print.feature_table_check <- function(x, ...) {
  cat("Same feature table:", if (x$same_feature_table) "YES" else "NO", "\n")
  s <- x$summary
  cat(sprintf("  Shared features:    %d (m/z differ: %d, RT differ: %d)\n",
              s$n_shared_features, s$n_tool_mz_mismatch, s$n_tool_rt_mismatch))
  cat(sprintf("  Only in grouping table: %d | only in feature table: %d (excluded from evaluation)\n",
              s$n_tool_not_in_ft, s$n_ft_not_in_tool))
  cat(sprintf("  MS1Match reference: %d checked, %d not in feature table, %d m/z differ\n",
              s$n_reference_features, s$n_reference_not_in_ft, s$n_reference_mz_mismatch))
  if (s$n_shared_features == 0L)
    cat("  No feature id is shared with the feature table.\n")
  if (!x$same_feature_table) {
    cat("  First mismatches:\n")
    print(utils::head(x$mismatches, 5), row.names = FALSE)
  }
  invisible(x)
}


# -----------------------------------------------------------------------------
# Main evaluation
# -----------------------------------------------------------------------------

#' Evaluate a feature-grouping tool against the MS1Match reference set
#'
#' @description
#' For each library compound, the Tier 1a \code{[M+H]+} feature is the anchor
#' and the group the tool assigns to it defines the compound's expected group.
#' Each reference feature (Phase 1 and Phase 2) gets the first status that
#' applies: Ref-unassigned (no-anchor) if the tool did not group the anchor,
#' Ref-unassigned if the anchor is grouped but the feature is not, otherwise
#' Ref-confirmed or Ref-conflicting. Every other feature of the tool is
#' Ref-unmatched, co-grouped (in a group with a compound anchor) or
#' Ref-unmatched, unrelated (in a group without an anchor, or in no group).
#'
#' @param anno_result Output of \code{join_annotation_results()}. Uses
#'   \code{export_feature_compound_summary} (reference rows) and
#'   \code{export_reduced_anno_table} (anchor search).
#' @param grouping Output of \code{\link{prepare_grouping_table}},
#'   \code{\link{camera_to_grouping}} or \code{\link{ramclustr_to_grouping}}.
#' @param tool_name Label used in printed output and in
#'   \code{\link{compare_grouping}}.
#' @param rt_window Optional \code{c(min, max)} in seconds. Applied only to
#'   features outside the reference set (Ref-unmatched), so their counts are
#'   comparable between tools.
#' @param feature_table The feature table MS1Match was run on. The function
#'   stops unless \code{\link{check_feature_match}} confirms that the
#'   grouping and MS1Match used the same table. If one is a subset of the
#'   other, only features present in both are evaluated.
#'
#' @details
#' The grouping table must list \emph{every} feature the tool processed,
#' including those it did not group (\code{group = NA}). A feature missing
#' from the grouping table is treated as outside the shared feature table and
#' excluded, not counted as Ref-unassigned.
#' @param mz_tol,rt_tol Passed to \code{\link{check_feature_match}}.
#'
#' @return Object of class \code{grouping_evaluation}: a list with
#'   \code{tool}, \code{benchmark} (tier table; \code{pct_grouped} is the
#'   share of reference features the tool grouped at all, i.e. its coverage),
#'   \code{status_counts},
#'   \code{full_table} (one row per reference feature x compound, plus one per
#'   feature outside the reference set), \code{anchors} and
#'   \code{feature_match}.
#' @export
evaluate_grouping <- function(anno_result, grouping, tool_name = "tool",
                              rt_window = NULL, feature_table,
                              mz_tol = 0.001, rt_tol = 1) {

  needed_res <- c("export_feature_compound_summary", "export_reduced_anno_table")
  if (!is.list(anno_result) || !all(needed_res %in% names(anno_result)))
    stop("evaluate_grouping: anno_result must be the output of ",
         "join_annotation_results() (needs ",
         paste(needed_res, collapse = ", "), ").")
  if (!is.data.frame(grouping) ||
      !all(c("feature_name", "group") %in% names(grouping)))
    stop("evaluate_grouping: grouping must have columns feature_name and group; ",
         "use prepare_grouping_table().")
  if (!"rt" %in% names(grouping)) grouping$rt <- NA_real_

  # --- 1. Same-feature-table check -------------------------------
  feature_match <- check_feature_match(anno_result, grouping, feature_table,
                                       mz_tol = mz_tol, rt_tol = rt_tol)
  if (!feature_match$same_feature_table) {
    print(feature_match)
    stop("evaluate_grouping: the grouping table and the MS1Match result do not ",
         "come from the given feature table (see mismatches above). MS1Match and ",
         "the grouping tool must be run on the same feature table.", call. = FALSE)
  }

  # --- 1b. Restrict to features present in both tables -------------------
  common   <- intersect(as.character(feature_table$feature_name),
                        grouping$feature_name)
  grouping <- grouping[grouping$feature_name %in% common, , drop = FALSE]
  n_ref_dropped <- sum(!unique(anno_result$export_feature_compound_summary$feature_name)
                       %in% c(common, NA))
  if (n_ref_dropped > 0L)
    message("evaluate_grouping: ", n_ref_dropped, " reference feature(s) are not ",
            "in the grouping table and are excluded (outside the shared feature table).")

  grp_lookup <- grouping[, c("feature_name", "group")]

  # --- 2. Anchors: one Tier 1a [M+H]+ feature per compound instance ----------
  full_anno  <- anno_result$export_reduced_anno_table
  adduct_col <- intersect(c("grp_anno_match", "adduct_match"), names(full_anno))[1L]
  if (is.na(adduct_col))
    stop("evaluate_grouping: export_reduced_anno_table has no grp_anno_match column.")

  anchors <- full_anno %>%
    dplyr::filter(
      !is.na(Lib_group),
      grepl("^Tier1a_", reference_tier),
      !grepl("exp", reference_tier, ignore.case = TRUE),
      grepl("^\\[M\\+H\\]\\+", .data[[adduct_col]])
    ) %>%
    dplyr::group_by(Lib_group) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(Lib_group, Compound_name, anchor_feature = feature_name) %>%
    dplyr::filter(anchor_feature %in% common) %>%
    dplyr::left_join(
      dplyr::rename(grp_lookup, anchor_feature = feature_name,
                    anchor_group = group),
      by = "anchor_feature")

  if (nrow(anchors) == 0L)
    stop("evaluate_grouping: no Tier 1a [M+H]+ anchors found in anno_result.")

  # --- 3. Reference rows with status ------------------------------------------
  # Phase 1 rows find their anchor via Lib_group; Phase 2 (_exp) rows have
  # Lib_group = NA and find it via Compound_name
  anchors_lib <- dplyr::select(anchors, Lib_group,
                               .af_lib = anchor_feature, .ag_lib = anchor_group)
  anchors_cpd <- anchors %>%
    dplyr::distinct(Compound_name, .keep_all = TRUE) %>%
    dplyr::select(Compound_name, .af_cpd = anchor_feature, .ag_cpd = anchor_group)

  ref <- anno_result$export_feature_compound_summary %>%
    dplyr::filter(!is.na(Compound_name), feature_name %in% common) %>%
    dplyr::mutate(tier = .tier_label(reference_tier)) %>%
    dplyr::left_join(grp_lookup, by = "feature_name") %>%
    dplyr::left_join(anchors_lib, by = "Lib_group") %>%
    dplyr::left_join(anchors_cpd, by = "Compound_name") %>%
    dplyr::mutate(
      .use_lib       = !is.na(.af_lib),
      anchor_feature = dplyr::if_else(.use_lib, .af_lib, .af_cpd),
      anchor_group   = dplyr::if_else(.use_lib, .ag_lib, .ag_cpd),
      # anchor first: without a grouped anchor there is no expected group
      ref_status = dplyr::case_when(
        is.na(anchor_feature) ~ "Not evaluable",               # technical check
        is.na(anchor_group)   ~ "Ref-unassigned (no-anchor)",
        is.na(group)          ~ "Ref-unassigned",
        group == anchor_group ~ "Ref-confirmed",
        TRUE                  ~ "Ref-conflicting"
      )
    ) %>%
    dplyr::select(-.af_lib, -.ag_lib, -.af_cpd, -.ag_cpd, -.use_lib)

  n_not_eval <- sum(ref$ref_status == "Not evaluable")
  if (n_not_eval > 0L)
    warning("evaluate_grouping: ", n_not_eval, " reference row(s) have no Tier 1a ",
            "[M+H]+ anchor in the shared feature table and are 'Not evaluable'. ",
            "This should not happen: check that the anchor is in the grouping ",
            "table (RT window) and that Compound_name matches between phases.")

  # --- 4. Tier benchmark table ------------------------------------------------
  ref_tiered <- ref %>%
    dplyr::filter(!is.na(tier)) %>%
    dplyr::mutate(
      tier       = factor(tier, levels = .tier_levels),
      ref_status = factor(ref_status, levels = .ref_status_levels)
    )

  counts <- ref_tiered %>%
    dplyr::count(tier, ref_status, .drop = FALSE) %>%
    tidyr::pivot_wider(names_from = ref_status, values_from = n,
                       values_fill = 0L) %>%
    dplyr::left_join(
      ref_tiered %>%
        dplyr::group_by(tier, .drop = FALSE) %>%
        dplyr::summarise(n_grouped = sum(!is.na(group)), .groups = "drop"),
      by = "tier") %>%
    dplyr::mutate(tier = as.character(tier))

  total <- counts %>%
    dplyr::summarise(dplyr::across(-tier, sum)) %>%
    dplyr::mutate(tier = "Total")

  pct <- function(x, n) ifelse(n > 0, round(100 * x / n, 1), NA_real_)

  benchmark <- dplyr::bind_rows(counts, total) %>%
    dplyr::mutate(n_ref = rowSums(dplyr::across(dplyr::all_of(.ref_status_levels)))) %>%
    dplyr::filter(n_ref > 0) %>%
    dplyr::transmute(
      Tier            = tier,
      n_ref,
      n_grouped,
      pct_grouped     = pct(n_grouped, n_ref),
      n_confirmed     = `Ref-confirmed`,
      pct_confirmed   = pct(n_confirmed, n_ref),
      n_conflicting   = `Ref-conflicting`,
      pct_conflicting = pct(n_conflicting, n_ref),
      n_unassigned    = `Ref-unassigned`,
      pct_unassigned  = pct(n_unassigned, n_ref),
      n_no_anchor     = `Ref-unassigned (no-anchor)`,
      pct_no_anchor   = pct(n_no_anchor, n_ref),
      n_not_evaluable = `Not evaluable`
    ) %>%
    as.data.frame()

  # --- 5. Features outside the reference set ---------------------------------
  anchor_groups <- unique(stats::na.omit(anchors$anchor_group))

  outside <- grouping %>%
    dplyr::filter(!feature_name %in% ref$feature_name)

  if (!is.null(rt_window)) {
    if (all(is.na(outside$rt))) {
      warning("evaluate_grouping: rt_window ignored, grouping table has no rt.")
    } else {
      outside <- dplyr::filter(outside, !is.na(rt),
                               rt >= rt_window[1], rt <= rt_window[2])
    }
  }

  outside <- outside %>%
    dplyr::mutate(
      .co_grouped  = !is.na(group) & group %in% anchor_groups,
      ref_status   = dplyr::if_else(.co_grouped, "Ref-unmatched, co-grouped",
                                    "Ref-unmatched, unrelated"),
      anchor_group = dplyr::if_else(.co_grouped, group, NA_character_)
    )

  # --- 6. Full table and status counts ---------------------------------------
  full_table <- dplyr::bind_rows(
    ref %>% dplyr::select(feature_name, Compound_name, Lib_group,
                          reference_tier, tier, group, anchor_feature,
                          anchor_group, ref_status),
    outside %>% dplyr::select(feature_name, group, anchor_group, ref_status)
  ) %>%
    as.data.frame()

  status_counts <- full_table %>%
    dplyr::mutate(ref_status = factor(ref_status, levels = .all_status_levels)) %>%
    dplyr::count(ref_status, name = "n", .drop = FALSE) %>%
    dplyr::mutate(pct = round(100 * n / sum(n), 1),
                  ref_status = as.character(ref_status)) %>%
    as.data.frame()

  structure(
    list(
      tool          = tool_name,
      benchmark     = benchmark,
      status_counts = status_counts,
      full_table    = full_table,
      anchors       = as.data.frame(anchors),
      feature_match = feature_match
    ),
    class = c("grouping_evaluation", "list")
  )
}


#' @export
print.grouping_evaluation <- function(x, ...) {
  cat("Grouping evaluation:", x$tool, "\n")
  s <- x$feature_match$summary
  cat("  Same feature table: YES\n")
  cat(sprintf("  Shared features evaluated: %d | reference features excluded: %d\n",
              s$n_shared_features, s$n_reference_not_in_tool))
  cat("  Compound anchors:", nrow(x$anchors), "\n\n")
  print(x$benchmark[, c("Tier", "n_ref", "pct_grouped", "pct_no_anchor",
                        "pct_unassigned", "pct_confirmed", "pct_conflicting")],
        row.names = FALSE)
  invisible(x)
}


#' Compare several grouping evaluations side by side
#'
#' @param ... Objects returned by \code{\link{evaluate_grouping}}.
#' @param metrics Benchmark columns to show per tool.
#'
#' @return \code{data.frame}: one row per tier, the reference size, and one
#'   column per tool x metric.
#' @export
compare_grouping <- function(..., metrics = c("pct_no_anchor",
                                              "pct_unassigned",
                                              "pct_confirmed",
                                              "pct_conflicting")) {
  evals <- list(...)
  stopifnot(length(evals) >= 1L,
            all(vapply(evals, inherits, logical(1), "grouping_evaluation")))

  tabs <- lapply(seq_along(evals), function(i) {
    e <- evals[[i]]
    b <- e$benchmark[, c("Tier", "n_ref", metrics)]
    names(b)[names(b) %in% metrics] <- paste0(e$tool, "_", metrics)
    if (i > 1L) b$n_ref <- NULL
    b
  })

  Reduce(function(a, b) dplyr::full_join(a, b, by = "Tier"), tabs)
}
