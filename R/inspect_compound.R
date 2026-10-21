#' Inspect a compound's CAMERA grouping performance
#'
#' @description
#' For a given compound name, presents a compound-centric evaluation of CAMERA's
#' pcgroup assignments using pre-built evaluation objects from the MS1Match
#' functions. Handles compounds with multiple chromatographic instances (suffix
#' \code{_1}, \code{_2}, etc.) automatically.
#'
#' Unlike the MS1Match functions which operate across all compounds, this
#' function filters the pre-built objects to the queried compound and presents
#' results in a human-readable format. The only computation performed internally
#' is building the feature-level \code{long_table}, which is unique to the
#' per-compound view.
#'
#' @section Tier definitions:
#' \tabular{ll}{
#'   \strong{Tier} \tab \strong{Description} \cr
#'   Tier1a        \tab Confirmed adduct match via direct library matching \cr
#'   Tier1a_exp    \tab Adduct match via Phase 2 full-FT expansion \cr
#'   Tier1b        \tab Isotope peaks (Phase 1) \cr
#'   Tier1b_exp    \tab Isotope peaks found in Phase 2 expansion \cr
#'   Tier1c        \tab In-source MS2 fragments \cr
#'   Tier1d        \tab Neutral loss pairs \cr
#'   Tier2         \tab Unexplained co-eluting features \cr
#' }
#'
#' @param compound_query \code{character(1)}. Compound name to search for.
#'   Matches exact name or name with numeric suffix (e.g. \code{"L-Norvaline"}
#'   matches \code{"L-Norvaline"}, \code{"L-Norvaline_1"}, \code{"L-Norvaline_2"}).
#'   Special regex characters are automatically escaped.
#' @param eval_full \code{data.frame}. Master evaluation table from
#'   \code{\link{build_eval_full}}. Required columns: \code{feature_name},
#'   \code{Compound_name}, \code{Lib_group}, \code{reference_tier},
#'   \code{mz}, \code{rt}, \code{pcgroup}, \code{camera_adduct},
#'   \code{adduct_match}.
#' @param anchors \code{data.frame}. Anchor table from
#'   \code{\link{build_anchors}}. Required columns: \code{Lib_group},
#'   \code{anchor_feature}, \code{anchor_pcgroup}.
#' @param compound_eval_wide \code{data.frame}. Per-compound wide evaluation
#'   table from \code{\link{build_compound_eval}}. Used for tier summary counts.
#' @param anno_accuracy \code{list}. Output of \code{\link{build_anno_accuracy}}.
#'   The \code{per_compound} element is used for annotation accuracy counts.
#' @param purity \code{list}. Output of \code{\link{build_pcgroup_purity}}.
#'   The \code{all_pcgroups} element is used for pcgroup context statistics.
#' @param CAMERA_annotated_pool_FT \code{data.frame}. CAMERA-annotated feature
#'   table. Must contain: \code{feature_name}, \code{pcgroup}.
#'   Optionally: \code{adduct} (or \code{camera_adduct}), \code{isotopes}.
#' @param ppm_tol \code{numeric(1)}. ppm tolerance for annotation accuracy check.
#'   Default \code{10}.
#' @param print_detail \code{logical(1)}. If \code{TRUE} (default), prints a
#'   formatted summary to the console / Rmd output.
#'
#' @return Invisibly returns a named \code{list} with three elements:
#' \describe{
#'   \item{\code{long_table}}{Feature-level detail table for the queried
#'     compound. Columns: \code{Compound_name}, \code{Lib_group},
#'     \code{feature_name}, \code{mz}, \code{rt}, \code{tier_simple},
#'     \code{our_annotation}, \code{camera_pcgroup}, \code{camera_location},
#'     \code{camera_adduct}, \code{camera_isotopes}.}
#'   \item{\code{tier_summary}}{Filtered rows from \code{compound_eval_wide}
#'     with annotation accuracy columns joined from \code{anno_accuracy$per_compound}.}
#'   \item{\code{pcgroup_purity}}{Filtered rows from \code{purity$all_pcgroups}
#'     for the anchor pcgroup(s) of this compound.}
#' }
#'
#' @seealso \code{\link{build_eval_full}}, \code{\link{build_anchors}},
#'   \code{\link{build_compound_eval}}, \code{\link{build_anno_accuracy}},
#'   \code{\link{build_pcgroup_purity}}
#'
#' @importFrom dplyr filter select mutate left_join group_by summarise distinct arrange bind_rows case_when n n_distinct pull
#' @importFrom stringr str_extract_all
#'
#' @examples
#' \dontrun{
#' eval_full          <- build_eval_full(anno_result$export_feature_compound_summary,
#'                                       CAMERA_annotated_pool_FT)
#' anchors            <- build_anchors(eval_full)
#' purity             <- build_pcgroup_purity(CAMERA_annotated_pool_FT,
#'                                            eval_full, anchors)
#' compound_eval_wide <- build_compound_eval(eval_full, anchors, purity)
#' anno_accuracy      <- build_anno_accuracy(eval_full, anchors, ppm_tol = 10)
#'
#' result <- inspect_compound(
#'   compound_query           = "L-Norvaline",
#'   eval_full                = eval_full,
#'   anchors                  = anchors,
#'   compound_eval_wide       = compound_eval_wide,
#'   anno_accuracy            = anno_accuracy,
#'   purity                   = purity,
#'   CAMERA_annotated_pool_FT = CAMERA_annotated_pool_FT
#' )
#'
#' write.csv(result$long_table, "L-Norvaline_detail.csv", row.names = FALSE)
#' }
#'
#' @export
inspect_compound <- function(compound_query,
                             eval_full,
                             anchors,
                             compound_eval_wide,
                             anno_accuracy,
                             purity,
                             CAMERA_annotated_pool_FT,
                             ppm_tol      = 10,
                             print_detail = TRUE) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    is.character(compound_query), length(compound_query) == 1L,
    is.data.frame(eval_full),
    is.data.frame(anchors),
    is.data.frame(compound_eval_wide),
    is.list(anno_accuracy), "per_compound" %in% names(anno_accuracy),
    is.list(purity), "all_pcgroups" %in% names(purity),
    is.data.frame(CAMERA_annotated_pool_FT)
  )

  tier_levels <- c("Tier1a", "Tier1a_exp", "Tier1b", "Tier1b_exp",
                   "Tier1c", "Tier1c_exp", "Tier1d", "Tier1d_exp",
                   "Tier2", "TOTAL")


  .tier_simple <- function(tier) {
    dplyr::case_when(
      grepl("1a.*exp|1a_exp", tier, ignore.case = TRUE) ~ "Tier1a_exp",
      grepl("1a",             tier, ignore.case = TRUE) ~ "Tier1a",
      grepl("1b.*exp|1b_exp", tier, ignore.case = TRUE) ~ "Tier1b_exp",
      grepl("1b",             tier, ignore.case = TRUE) ~ "Tier1b",
      grepl("1c.*exp|1c_exp|MS2Fragment_exp", tier, ignore.case = TRUE) ~ "Tier1c_exp",
      grepl("1c",             tier, ignore.case = TRUE) ~ "Tier1c",
      grepl("1d.*exp|1d_exp", tier, ignore.case = TRUE) ~ "Tier1d_exp",
      grepl("1d",             tier, ignore.case = TRUE) ~ "Tier1d",
      grepl("Tier2|tier2",    tier, ignore.case = TRUE) ~ "Tier2",
      TRUE                                              ~ "Other"
    )
  }

  # Extract numeric neutral mass from adduct string e.g. "[M+H]+ 195.065"
  .extract_mz <- function(x) {
    m <- regmatches(x, regexpr("[0-9]+\\.[0-9]+$", x))
    if (length(m) == 0L || identical(m, character(0))) return(NA_real_)
    as.numeric(m)
  }

  # ── 1. Find matching compound names ────────────────────────────────────────

  query_esc <- gsub("([\\^$.|?*+(){}\\[\\]])", "\\\\\\1", compound_query,
                    perl = TRUE)
  pattern   <- paste0("^", query_esc, "(_[0-9]+)?$")

  matched_names <- unique(eval_full$Compound_name[
    grepl(pattern, eval_full$Compound_name)
  ])

  if (length(matched_names) == 0L) {
    stop("inspect_compound: '", compound_query,
         "' not found in eval_full$Compound_name.\n",
         "  Check spelling or use a partial string.")
  }

  # ── 2. Filter eval_full to this compound ───────────────────────────────────
  ef_cmp <- eval_full %>%
    dplyr::filter(Compound_name %in% matched_names) %>%
    dplyr::mutate(tier_simple = .tier_simple(reference_tier))

  # ── 3. Detect adduct / isotope columns in CAMERA FT ───────────────────────
  cam_adduct_col  <- if ("camera_adduct" %in% names(CAMERA_annotated_pool_FT))
    "camera_adduct"
  else if ("adduct" %in% names(CAMERA_annotated_pool_FT))
    "adduct"
  else
    NULL

  cam_iso_col <- if ("isotopes" %in% names(CAMERA_annotated_pool_FT))
    "isotopes"
  else
    NULL

  # ── 4. Build long_table: one row per feature ───────────────────────────────
  ef_cmp <- ef_cmp %>%
    dplyr::mutate(
      our_annotation = mapply(
        function(grp, exp, nl, iso) {
          parts <- c(grp, exp, nl, iso)
          parts <- parts[!is.na(parts) & nchar(trimws(parts)) > 0L]
          if (length(parts) == 0L) return(NA_character_)
          paste(parts, collapse = "; ")
        },
        if ("grp_anno_match" %in% names(ef_cmp)) ef_cmp[["grp_anno_match"]] else
          if ("adduct_match"    %in% names(ef_cmp)) ef_cmp[["adduct_match"]]     else rep(NA, nrow(ef_cmp)),
        if ("exp_anno_match" %in% names(ef_cmp)) ef_cmp[["exp_anno_match"]] else rep(NA, nrow(ef_cmp)),
        if ("NL_match"         %in% names(ef_cmp)) NL_match         else rep(NA, nrow(ef_cmp)),
        if ("iso_annotation"   %in% names(ef_cmp)) iso_annotation   else rep(NA, nrow(ef_cmp)),
        USE.NAMES = FALSE
      )
    )

  # Get anchor pcgroup for each Lib_group.
  # Phase 2 ISF features have Lib_group = NA but ms2_lib_group set — include
  # those too so their camera_location resolves correctly in long_table.
  all_lib_groups <- unique(c(
    ef_cmp$Lib_group,
    if ("ms2_lib_group" %in% names(ef_cmp)) ef_cmp$ms2_lib_group else NULL
  ))
  all_lib_groups <- all_lib_groups[!is.na(all_lib_groups)]

  anchor_map <- anchors %>%
    dplyr::filter(Lib_group %in% all_lib_groups) %>%
    dplyr::select(Lib_group, anchor_pcgroup)

  # Phase 2 expansion rows (Tier1a_exp_AnnoMatch, Tier1b_exp_Isotope, etc.)

  anchor_map_by_cpd <- anchors %>%
    dplyr::filter(Compound_name %in% ef_cmp$Compound_name) %>%
    dplyr::select(Compound_name, anchor_pcgroup_cpd = anchor_pcgroup) %>%
    dplyr::distinct(Compound_name, .keep_all = TRUE)

  # CAMERA FT columns for this compound's features
  cam_cols <- c("feature_name", "pcgroup",
                cam_adduct_col, cam_iso_col)
  cam_cols <- cam_cols[!is.null(cam_cols)]

  cam_sub <- CAMERA_annotated_pool_FT %>%
    dplyr::filter(feature_name %in% ef_cmp$feature_name) %>%
    dplyr::select(dplyr::any_of(cam_cols)) %>%
    dplyr::rename(camera_pcgroup = pcgroup)

  if (!is.null(cam_adduct_col))
    cam_sub <- dplyr::rename(cam_sub, camera_adduct = dplyr::any_of(cam_adduct_col))
  if (!is.null(cam_iso_col))
    cam_sub <- dplyr::rename(cam_sub, camera_isotopes = dplyr::any_of(cam_iso_col))

  # For Phase 2 ISF features (Lib_group = NA, ms2_lib_group set), join
  # anchor_map via ms2_lib_group
  ef_cmp_join <- ef_cmp %>%
    dplyr::mutate(
      .join_lib = dplyr::coalesce(
        .data$Lib_group,
        if ("ms2_lib_group" %in% names(.)) .data$ms2_lib_group else NA_character_
      )
    )

  long_table <- ef_cmp_join %>%
    dplyr::left_join(anchor_map, by = c(".join_lib" = "Lib_group")) %>%
    dplyr::select(-.join_lib) %>%
    dplyr::left_join(anchor_map_by_cpd, by = "Compound_name") %>%
    dplyr::mutate(
      anchor_pcgroup = dplyr::coalesce(.data$anchor_pcgroup, .data$anchor_pcgroup_cpd)
    ) %>%
    dplyr::select(-anchor_pcgroup_cpd) %>%
    dplyr::left_join(cam_sub,    by = "feature_name") %>%
    dplyr::mutate(
      camera_location = dplyr::case_when(
        is.na(camera_pcgroup)                          ~ "not in CAMERA FT",
        is.na(anchor_pcgroup)                          ~ "no anchor defined",
        camera_pcgroup == anchor_pcgroup               ~ "co-grouped with anchor",
        TRUE                                           ~ "split to other pcgroup"
      )
    ) %>%
    dplyr::select(
      Compound_name, Lib_group, feature_name,
      mz, rt,
      tier_simple, our_annotation,
      camera_pcgroup, camera_location,
      dplyr::any_of(c("camera_adduct", "camera_isotopes",
                      "ms2_lib_group"))
    ) %>%
    dplyr::arrange(Lib_group, tier_simple, mz)

  # ── 5. Tier summary — from pre-built compound_eval_wide ───────────────────
  tier_summary <- compound_eval_wide %>%
    dplyr::filter(Lib_group %in% ef_cmp$Lib_group)


  if (!is.null(anno_accuracy$per_compound) &&
      nrow(anno_accuracy$per_compound) > 0L) {

    tier_col_acc <- if ("tier_simple"     %in% names(anno_accuracy$per_compound))
      "tier_simple"
    else if ("reference_tier" %in% names(anno_accuracy$per_compound))
      "reference_tier"
    else
      NULL

    # Select all accuracy columns that are actually present
    acc_cols_want <- c("Lib_group", tier_col_acc,
                       "n_total", "n_unannotationled", "n_annotationled",
                       "n_correct", "n_wrong", "coverage_pct", "accuracy_pct")
    acc_cols <- intersect(acc_cols_want, names(anno_accuracy$per_compound))

    acc_sub <- anno_accuracy$per_compound %>%
      dplyr::filter(Lib_group %in% ef_cmp$Lib_group) %>%
      dplyr::select(dplyr::all_of(acc_cols))

    tier_summary <- tier_summary %>%
      dplyr::left_join(acc_sub, by = "Lib_group")
  }

  # ── 6. pcgroup purity — from pre-built purity$all_pcgroups ───────────────
  anchor_pcgs <- anchors %>%
    dplyr::filter(Lib_group %in% ef_cmp$Lib_group) %>%
    dplyr::pull(anchor_pcgroup)

  pcgroup_purity <- purity$all_pcgroups %>%
    dplyr::filter(pcgroup %in% anchor_pcgs)

  # ── 7. Optional output ─────────────────────────────────────────
  if (isTRUE(print_detail)) {


    lib_grp_levels <- unique(ef_cmp$Lib_group)
    lib_grp_levels <- lib_grp_levels[!is.na(lib_grp_levels)]

    for (lib_grp in lib_grp_levels) {

      anc   <- anchors %>% dplyr::filter(Lib_group == lib_grp)
      pcg_id <- if (nrow(anc) > 0L) anc$anchor_pcgroup[1L] else NA

      pcg_context <- purity$all_pcgroups %>%
        dplyr::filter(pcgroup == pcg_id)

      has_ms2_lib_lt <- "ms2_lib_group" %in% names(long_table)
      p2_isf_fns <- if (has_ms2_lib_lt)
        long_table$feature_name[
          is.na(long_table$Lib_group) &
          !is.na(long_table$ms2_lib_group) &
          long_table$ms2_lib_group == lib_grp
        ]
      else
        character(0L)


      cmp_lib_groups_here <- unique(ef_cmp$Lib_group[
        !is.na(ef_cmp$Lib_group) &
        ef_cmp$Compound_name %in% ef_cmp$Compound_name[ef_cmp$Lib_group == lib_grp]
      ])
      p2_other_fns <- if (length(cmp_lib_groups_here) == 1L)
        long_table$feature_name[
          is.na(long_table$Lib_group) &
          !(long_table$feature_name %in% p2_isf_fns) &
          long_table$Compound_name %in%
            ef_cmp$Compound_name[ef_cmp$Lib_group == lib_grp]
        ]
      else
        character(0L)

      lt_grp <- long_table %>%
        dplyr::filter(
          (!is.na(Lib_group) & Lib_group == lib_grp) |
          feature_name %in% p2_isf_fns |
          feature_name %in% p2_other_fns
        )

      cat("\n", strrep("=", 72), "\n")
      cat("COMPOUND :", lib_grp, "\n")
      cat(strrep("-", 72), "\n")

      if (nrow(anc) > 0L) {
        cat(sprintf("  Anchor feature : %s\n", anc$anchor_feature[1L]))
        cat(sprintf("  Anchor pcgroup : %s\n", pcg_id))
        cat(sprintf("  Anchor mz      : %.6f\n", anc$anchor_mz[1L]))
      } else {
        cat("  *** No [M+H]+ anchor found for this instance ***\n")
      }

      if (nrow(pcg_context) > 0L) {
        cat(sprintf("  pcgroup size   : %d features, %d compound(s)\n",
                    pcg_context$pcgroup_n_features[1L],
                    pcg_context$pcgroup_n_compounds[1L]))
        if (pcg_context$pcgroup_outside_ref[1L] > 0L)
          cat(sprintf("  Outside ref    : %d features not in reference set\n",
                      pcg_context$pcgroup_outside_ref[1L]))

        # --- Pcgroup purity -----------------

        cat(sprintf("  Pcgroup purity : %d compound(s) anchored in this pcgroup\n",
                    pcg_context$pcgroup_n_compounds[1L]))
        n_owngrp <- sum(lt_grp$camera_location == "co-grouped with anchor",
                        na.rm = TRUE)
        if (pcg_context$pcgroup_n_features[1L] > 0L) {
          refconf_pct <- round(100 * n_owngrp / pcg_context$pcgroup_n_features[1L], 1)
          cat(sprintf("  Ref-confirmed %% : %.1f%% (%d of %d pcgroup features are this compound's own)\n",
                      refconf_pct, n_owngrp, pcg_context$pcgroup_n_features[1L]))
        }
      }

      cat("\n")
      cat(sprintf("  %-12s | %-10s | %8s | %s\n",
                  "Tier", "feature", "Location", "Annotation"))
      cat(sprintf("  %s\n", strrep("-", 65)))

      for (tier in tier_levels[tier_levels != "TOTAL"]) {
        rows <- lt_grp %>% dplyr::filter(tier_simple == tier)
        if (nrow(rows) == 0L) next
        for (i in seq_len(nrow(rows))) {
          r <- rows[i, ]
          loc_short <- dplyr::case_when(
            r$camera_location == "co-grouped with anchor"  ~ "grouped",
            r$camera_location == "split to other pcgroup"  ~ "SPLIT",
            r$camera_location == "not in CAMERA FT"        ~ "not in FT",
            TRUE                                            ~ r$camera_location
          )
          ann <- if (!is.na(r$our_annotation)) r$our_annotation else "—"
          if (nchar(ann) > 45L)
            ann <- paste0(substr(ann, 1L, 42L), "...")
          cat(sprintf("  %-12s | %-10s | %8s | %s\n",
                      if (i == 1L) tier else "",
                      r$feature_name,
                      loc_short,
                      ann))
        }
      }

      cat("\n")
    }
  }

  # ── Return ──────────────────────────────────────────────────────────────────
  invisible(list(
    long_table     = long_table,
    tier_summary   = tier_summary,
    pcgroup_purity = pcgroup_purity
  ))
}
