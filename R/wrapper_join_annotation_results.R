#' Join lists of annotation results to library match output
#'
#' @description
#' Combine direct library match results (MS1/MS2) with downstream annotation
#' results (adduct/isotope/MS2 fragment/NL) onto the feature table.
#'
#' Tier naming convention:
#' \itemize{
#'   \item \strong{1a} — adduct match (Phase 1, within library group)
#'   \item \strong{1b} — isotope (Phase 1, within library group)
#'   \item \strong{1c} — MS2 fragment match (Phase 1: feature already in a
#'     Phase 1 library group)
#'   \item \strong{1d} — neutral loss (Phase 1, pairwise within group)
#'   \item \strong{1a_exp / 1b_exp / 1c_exp / 1d_exp} — Phase 2 expansions
#'     (feature NOT assigned in Phase 1; derived from a Phase 1 confirmed M)
#'   \item \strong{2} — unexplained (in group, no rule fires)
#' }
#'
#'
#' Tiers are mutually exclusive — one tier per feature, highest priority
#' wins: \strong{1a > 1b > 1c > 1d > 2 > 1a_exp > 1b_exp > 1c_exp > 1d_exp}.
#' \code{Tier2} (Unexplained) is a Phase 1 assignment and therefore
#' outranks all Phase 2 expansions.
#'
#' The \code{annotation_phase} column takes exactly two values:
#' \code{"Phase1"} (any Tier 1a/1b/1c/1d or Tier 2 assignment) and
#' \code{"Phase2"} (any \code{_exp} tier).  Sub-labels such as
#' \code{Phase1_MS2frag} or \code{Phase2_IsoExpansion} are not used;
#' \code{reference_tier} already encodes the specific evidence type.
#'
#' @param FT Required data.frame. Feature table.
#' @param direct_lib_list Required list.
#' @param lib_anno_list Optional list/data.frame. Output of annotate().
#' @param target_match_list Optional list/data.frame. Target list annotations.
#' @param isf_result Optional data.frame. Output of \code{annotateISF()} /
#'   \code{run_isf_annotation()}.  Adds Tier 1c / 1c_exp columns and
#'   resolves \code{reference_tier} by priority.
#' @param include_phase2 Logical. If \code{FALSE}, Phase 2 expansion rows
#'   (\code{Tier1a_exp}, \code{Tier1b_exp}, \code{Tier1d_exp}) from
#'   \code{lib_anno_list} are excluded from the output.  Default \code{TRUE}.
#' @param verbose Logical. If \code{FALSE}, suppresses all progress messages.
#'   Default \code{TRUE}.
#'
#' @export
join_annotation_results <- function(FT,
                                    direct_lib_list,
                                    lib_anno_list     = NULL,
                                    target_match_list = NULL,
                                    isf_result        = NULL,
                                    include_phase2    = TRUE,
                                    verbose           = TRUE) {


  stopifnot(is.data.frame(FT))
  stopifnot(is.list(direct_lib_list))

  # ---------------------------------------------------------------------------
  # helpers
  # ---------------------------------------------------------------------------
  safe_max <- function(x) {
    v <- suppressWarnings(max(x, na.rm = TRUE))
    if (is.infinite(v)) NA_real_ else v
  }

  drop_overlap_except_key <- function(y, x_names, key) {
    overlap <- setdiff(intersect(names(y), x_names), key)
    if (length(overlap) == 0L) return(y)
    dplyr::select(y, -dplyr::any_of(overlap))
  }

  target_list_to_df <- function(res_list) {
    if (is.null(res_list))       return(NULL)
    if (is.data.frame(res_list)) return(.rename_target_cols(res_list))
    if (!is.list(res_list))
      stop("target_match_list must be a list or data.frame.")
    res_list <- Filter(function(x) {
      is.list(x) && !is.null(x$Feature_name) && length(x$Feature_name) > 0
    }, res_list)
    if (length(res_list) == 0L) return(NULL)
    out <- do.call(rbind, lapply(res_list, function(x) {
      data.frame(
        Feature_name          = x$Feature_name,
        Feature_mz            = x$Feature_mz,
        Comp_name             = x$Comp_name,
        PI_name               = x$PI_name,
        Comp_mz_targetMatch   = if (!is.null(x$Comp_mz))   x$Comp_mz   else NA_real_,
        mz_diff_targetMatch   = if (!is.null(x$mz_diff))   x$mz_diff   else NA_real_,
        ppm_error_targetMatch = if (!is.null(x$ppm_error)) x$ppm_error else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
    rownames(out) <- NULL
    out
  }

  .rename_target_cols <- function(df) {
    rename_map <- c(
      Comp_mz   = "Comp_mz_targetMatch",
      mz_diff   = "mz_diff_targetMatch",
      ppm_error = "ppm_error_targetMatch"
    )
    for (old in names(rename_map)) {
      new <- rename_map[[old]]
      if (old %in% names(df) && !new %in% names(df))
        names(df)[names(df) == old] <- new
    }
    df
  }

  # ---------------------------------------------------------------------------
  # Tier priority resolve
  # ---------------------------------------------------------------------------
  .tier_priority <- function(tier) {
    tier <- as.character(tier)
    # Actual tier strings emitted by the codebase:
    #   Tier1a_AnnoMatch        (cpp Phase 1 — feature in Phase 1 lib group)
    #   Tier1b_Isotope            (cpp Phase 1 — feature in Phase 1 lib group)
    #   Tier1b_exp_Isotope        (cpp Phase 2 — feature NOT in Phase 1)
    #   Tier1c_MS2Fragment        (annotateISF Phase 1 — feature in Phase 1 lib
    #                              group AND parent is the PI anchor, any adduct)
    #   Tier1c_MS2Fragment_exp    (annotateISF Phase 2 — feature NOT assigned in
    #                              Phase 1, OR parent is not the PI anchor;
    #                              downgraded from Tier1c in either case)
    #   Tier1d_NeutralLoss        (cpp Phase 1 — feature in Phase 1 lib group)
    #   Tier1a_exp_AnnoMatch    (cpp Phase 2 — feature NOT in Phase 1)
    #   Tier1c_exp_NeutralLoss    (cpp Phase 2 — feature NOT in Phase 1)
    #   Tier2_Unexplained         (cpp — in group but no Tier-1 rule fires)
    # Priority: 1a > 1b > 1c > 1d > 2 > 1a_exp > 1b_exp > 1c_exp > 1d_exp
    # Tier2 (Unexplained) is a Phase 1 assignment — the feature IS in a
    # library group but no Tier-1 rule fired.
    dplyr::case_when(
      # ── Phase 1 ─────────────────────────────────────────────────────────
      grepl("^Tier1a_AnnoMatch$|^Tier1a$|^1a$",
                                   tier, ignore.case = TRUE) ~ 1L,
      grepl("^Tier1b",             tier, ignore.case = TRUE) ~ 2L,   # Tier1b_Isotope
      grepl("^Tier1c_MS2Fragment$|^Tier1c$|^1c$",
                                   tier, ignore.case = TRUE) ~ 3L,
      grepl("^Tier1d_NeutralLoss$|^Tier1d$|^1d$",
                                   tier, ignore.case = TRUE) ~ 4L,
      grepl("^Tier2|^2$",          tier, ignore.case = TRUE) ~ 5L,   # Unexplained — Phase 1 group, no rule
      # ── Phase 2 expansions (feature not assigned in Phase 1) ─────────────
      grepl("^Tier1a_exp",         tier, ignore.case = TRUE) ~ 6L,   # Tier1a_exp_AnnoMatch
      grepl("^Tier1b_exp",         tier, ignore.case = TRUE) ~ 7L,   # Tier1b_exp_Isotope
      grepl("Tier1c_MS2Fragment_exp|Tier1c_exp|1c_exp",
                                   tier, ignore.case = TRUE) ~ 8L,   # Tier1c_MS2Fragment_exp
      grepl("^Tier1c_exp_NeutralLoss|1d_exp",
                                   tier, ignore.case = TRUE) ~ 9L,   # Tier1d_exp NeutralLoss
      # ── Catch-all ────────────────────────────────────────────────────────
      grepl("^Tier3|^3$",          tier, ignore.case = TRUE) ~ 10L,
      TRUE                                                    ~ 99L
    )
  }

  # ---------------------------------------------------------------------------
  # Detect FT column names
  # ---------------------------------------------------------------------------
  mz_col  <- if ("mzmed"          %in% names(FT)) "mzmed"          else
             if ("mz"             %in% names(FT)) "mz"             else NULL
  rt_col  <- if ("rtmed"          %in% names(FT)) "rtmed"          else
             if ("rt"             %in% names(FT)) "rt"             else NULL
  int_col <- if ("intensity_mean" %in% names(FT)) "intensity_mean" else
             if ("intensity"      %in% names(FT)) "intensity"      else NULL

  if (is.null(mz_col) || is.null(rt_col))
    stop("FT must contain mz/rt columns (mzmed/mz and rtmed/rt).")

  # ---------------------------------------------------------------------------
  # direct_lib_list -> direct_df
  # ---------------------------------------------------------------------------
  parallel_peak_fields <- c("Compound_peak_mz", "Compound_peak_int")
  full_peak_fields     <- c("Compound_peak_mz_full", "Compound_peak_int_full")

  expand_direct_element <- function(x) {
    matched_feats <- x[["Matched_features_to_peaks"]]
    if (is.null(matched_feats) || length(matched_feats) == 0 ||
        (length(matched_feats) == 1 && is.na(matched_feats)))
      matched_feats <- x[["feature_name"]]
    n   <- length(matched_feats)
    out <- lapply(names(x), function(nm) {
      v <- x[[nm]]
      if (nm %in% parallel_peak_fields) {
        if (length(v) == n) as.character(v) else rep(NA_character_, n)
      } else if (nm %in% full_peak_fields) {
        rep(paste(v[!is.na(v)], collapse = "; "), n)
      } else if (nm == "Matched_features_to_peaks") {
        as.character(matched_feats)
      } else if (length(v) == n) {
        v
      } else {
        rep(if (length(v) == 0) NA else v[[1]], n)
      }
    })
    names(out) <- names(x)
    as.data.frame(out, stringsAsFactors = FALSE)
  }

  direct_df <- tryCatch({
    dplyr::bind_rows(purrr::map(direct_lib_list, expand_direct_element))
  }, error = function(e) {
    stop("Failed to expand/bind direct_lib_list: ", conditionMessage(e))
  })

  req_direct <- c("feature_name", "Matched_features_to_peaks",
                  "Lib_group", "Compound_name")
  miss_direct <- setdiff(req_direct, names(direct_df))
  if (length(miss_direct) > 0)
    stop("Missing required columns in direct match results: ",
         paste(miss_direct, collapse = ", "))

  direct_df <- direct_df %>%
    dplyr::rename(PI_feature_name = feature_name,
                  feature_name    = Matched_features_to_peaks)

  # ── PI feature: one row per (PI_feature_name, Lib_group, Compound_name)
  pi_rows <- direct_df %>%
    dplyr::select(PI_feature_name, Lib_group, Compound_name,
                  dplyr::any_of(c("Compound_M", "Compound_rt",
                                  "Ion_type", "Match_peak_score",
                                  "inner_cosine", "outer_cosine",
                                  "Matched_features_mz", "Matched_features_rt",
                                  "Matched_features_int", "intensity_ratio",
                                  "Compound_peak_mz", "Compound_peak_int",
                                  "Compound_peak_mz_full", "Compound_peak_int_full",
                                  "Total_library_peaks", "Matched_peaks"))) %>%
    dplyr::rename(feature_name = PI_feature_name) %>%
    dplyr::filter(!is.na(feature_name), nchar(feature_name) > 0L) %>%
    # skip if PI feature already appears as a matched pool feature
    dplyr::filter(!feature_name %in% direct_df$feature_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(PI_feature_name = feature_name)

  # ---------------------------------------------------------------------------
  # lib_anno_list to anno_df
  # ---------------------------------------------------------------------------
  anno_df <- NULL
  if (!is.null(lib_anno_list)) {
    anno_df <- if (is.data.frame(lib_anno_list)) {
      lib_anno_list
    } else {
      dplyr::bind_rows(purrr::map(lib_anno_list, function(x)
        as.data.frame(x, stringsAsFactors = FALSE)))
    }
  }
  # ---------------------------------------------------------------------------
  # Merge direct_df + anno_df
  # Phase 1: left_join anno columns onto direct_df rows
  # Phase 2: add expansion rows absent from direct_df
  # ---------------------------------------------------------------------------
  match_df <- direct_df
  if (!is.null(anno_df) && nrow(anno_df) > 0L) {
    anno_clean <- drop_overlap_except_key(anno_df, names(match_df),
                                          "feature_name")
    # join on feature_name,  Lib_group

    tier_only_in_anno <- c("reference_tier", "annotation_phase",
                           "NL_match", "iso_annotation",
                           "grp_anno_match", "grp_anno_ppm",
                           "grp_anno_mda", "grp_anno_confidence",
                           "exp_anno_match", "exp_anno_ppm",
                           "exp_anno_mda", "exp_anno_confidence",
                           "cutoff_to_MH")
    match_df <- dplyr::select(match_df,
                              -dplyr::any_of(tier_only_in_anno))

    anno_clean2 <- drop_overlap_except_key(
      anno_df, names(match_df),
      key = c("feature_name", "Lib_group")
    )
    join_keys <- intersect(c("feature_name", "Lib_group"), names(anno_clean2))
    if (verbose) message("[join_annotation_results] Joining anno on: ",
                         paste(join_keys, collapse = " + "))
    match_df <- dplyr::left_join(match_df, anno_clean2, by = join_keys,
                                 relationship = "many-to-many")

    # ── add PI rows into match_df ────────────────────────────────────────────

    if (nrow(pi_rows) > 0L) {
      pi_rows <- dplyr::left_join(
        pi_rows,
        anno_df %>% dplyr::select(-dplyr::any_of(names(pi_rows)[
          names(pi_rows) %in% names(anno_df) &
          !names(pi_rows) %in% c("feature_name", "Lib_group", "Compound_name")])),
        by = c("feature_name", "Lib_group", "Compound_name"),
        relationship = "many-to-many"
      )

      already_in  <- paste(match_df$feature_name, match_df$Lib_group,
                           match_df$Compound_name, sep = "||")
      pi_key      <- paste(pi_rows$feature_name, pi_rows$Lib_group,
                           pi_rows$Compound_name, sep = "||")
      pi_rows_new <- pi_rows[!pi_key %in% already_in, ]

      if (nrow(pi_rows_new) > 0L) {
        if (verbose) message("[join_annotation_results] Adding ",
                             nrow(pi_rows_new), " PI precursor feature rows.")
        match_df <- dplyr::bind_rows(match_df, pi_rows_new)
      }
    }

    # add Phase 2 expansion features
    phase2_features <- anno_df[
      !anno_df$feature_name %in% direct_df$feature_name &
      grepl("Phase2", anno_df$annotation_phase, ignore.case = TRUE), ,
      drop = FALSE
    ]
    if (include_phase2 && nrow(phase2_features) > 0L) {
      if (verbose) message("[join_annotation_results] Adding ",
                           nrow(phase2_features),
                           " Phase 2 expansion features.")
      match_df <- dplyr::bind_rows(match_df, phase2_features)
    } else if (!include_phase2 && verbose) {
      message("[join_annotation_results] Phase 2 expansion skipped ",
              "(include_phase2 = FALSE).")
    }
  }

  # ---------------------------------------------------------------------------
  # target_match_list
  # ---------------------------------------------------------------------------
  target_df           <- target_list_to_df(target_match_list)
  target_numeric_wide <- NULL
  if (!is.null(target_df) && nrow(target_df) > 0L) {
    target_numeric_wide <- target_df %>%
      dplyr::group_by(Feature_name) %>%
      dplyr::slice_min(order_by = abs(mz_diff_targetMatch), n = 1,
                       with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::rename(feature_name    = Feature_name,
                    targetList_anno = PI_name) %>%
      dplyr::select(feature_name, targetList_anno,
                    dplyr::any_of(c("Comp_mz_targetMatch",
                                    "mz_diff_targetMatch",
                                    "ppm_error_targetMatch")))
  }
  if (!is.null(target_numeric_wide)) {
    match_df <- dplyr::left_join(match_df, target_numeric_wide,
                                 by = "feature_name",
                                 relationship = "many-to-many") %>%
      dplyr::distinct()
  }

  # ---------------------------------------------------------------------------
  # ISF / MS2 fragment result
  # ---------------------------------------------------------------------------
  ms2_cols <- c("parent_feature", "matched_fragment_mz",
                "isf_ppm", "isf_mda", "rt_diff",
                "reference_tier", "annotation_phase", "Lib_group",
                "ms2_lib_group")

  if (!is.null(isf_result) && is.data.frame(isf_result) &&
      nrow(isf_result) > 0L) {


    parent_lib <- unique(anno_df[
      !is.na(anno_df$Lib_group),
      c("feature_name", "Lib_group")])
    names(parent_lib) <- c("parent_feature", "Lib_group_parent")

    # Split ISF result into Phase 1 (Lib_group already set by annotator) and
    # Phase 2 (Lib_group = NA)
    isf_base <- isf_result %>%
      dplyr::select(feature_name, dplyr::any_of(ms2_cols))

    isf_p1 <- isf_base %>%
      dplyr::filter(!is.na(Lib_group))

    isf_p2 <- isf_base %>%
      dplyr::filter(is.na(Lib_group)) %>%
      dplyr::left_join(parent_lib, by = "parent_feature",
                       relationship = "many-to-many") %>%
      dplyr::mutate(

        ms2_lib_group = .data$Lib_group_parent
      ) %>%
      dplyr::select(-dplyr::any_of("Lib_group_parent"))

    isf_slim <- dplyr::bind_rows(isf_p1, isf_p2) %>%
      dplyr::rename(
        ms2_parent_feature      = parent_feature,
        ms2_matched_fragment_mz = matched_fragment_mz,
        ms2_ppm                 = isf_ppm,
        ms2_mda                 = isf_mda,
        ms2_rt_diff             = rt_diff,
        ms2_tier                = reference_tier,
        ms2_annotation_phase    = annotation_phase
      )

    # ── Anchor parent check ───────────────────────────────────────────────────

    if ("reference_tier" %in% names(anno_df)) {
      anchor_parents <- anno_df %>%
        dplyr::filter(!is.na(Lib_group),
                      grepl("^Tier1a_AnnoMatch", reference_tier,
                            ignore.case = FALSE)) %>%
        dplyr::select(ms2_parent_feature = feature_name, Lib_group) %>%
        dplyr::distinct() %>%
        dplyr::mutate(.is_anchor_parent = TRUE)

      isf_slim <- isf_slim %>%
        dplyr::left_join(anchor_parents,
                         by = c("ms2_parent_feature", "Lib_group")) %>%
        dplyr::mutate(
          ms2_tier = dplyr::case_when(
            ms2_tier == "Tier1c_MS2Fragment" & is.na(.is_anchor_parent) ~
              "Tier1c_MS2Fragment_exp",
            TRUE ~ ms2_tier
          ),
          ms2_annotation_phase = dplyr::case_when(
            ms2_tier == "Tier1c_MS2Fragment_exp" &
              ms2_annotation_phase == "Phase1_MS2frag" ~ "Phase2_MS2frag",
            TRUE ~ ms2_annotation_phase
          )
        ) %>%
        dplyr::select(-dplyr::any_of(".is_anchor_parent"))
    }

    # join on feature_name + Lib_group
    match_df <- dplyr::left_join(match_df, isf_slim,
                                 by = c("feature_name", "Lib_group"),
                                 relationship = "many-to-many")

    # ── Add ISF Phase 2 child features not yet in match_df ───────────────────

    if (include_phase2) {
      isf_exp_new <- isf_slim %>%
        dplyr::filter(ms2_tier == "Tier1c_MS2Fragment_exp",
                      !is.na(ms2_lib_group)) %>%
        dplyr::mutate(
          reference_tier   = ms2_tier,
          annotation_phase = ms2_annotation_phase
        )

      if (nrow(isf_exp_new) > 0L) {

        lg_cpd <- unique(anno_df[, c("Lib_group", "Compound_name")])

        # Look up mz / rt / intensity from FT
        ft_cols <- c("feature_name",
                     if ("mzmed" %in% names(FT)) "mzmed" else "mz",
                     if ("rtmed" %in% names(FT)) "rtmed" else "rt",
                     if ("intensity_median" %in% names(FT)) "intensity_median" else NULL)
        ft_sub  <- FT[, intersect(ft_cols, names(FT)), drop = FALSE]

        isf_exp_rows <- isf_exp_new %>%
          dplyr::left_join(lg_cpd,
                           by = c("ms2_lib_group" = "Lib_group")) %>%
          dplyr::left_join(ft_sub,  by = "feature_name") %>%
          dplyr::rename_with(~ sub("^mzmed$", "mz", .),    dplyr::any_of("mzmed")) %>%
          dplyr::rename_with(~ sub("^rtmed$", "rt", .),    dplyr::any_of("rtmed")) %>%
          dplyr::rename_with(~ sub("^intensity_median$", "intensity", .),
                             dplyr::any_of("intensity_median"))

        # Dedup on feature_name only — Lib_group is NA for all exp rows
        already_key  <- match_df$feature_name[is.na(match_df$Lib_group)]
        isf_exp_rows <- isf_exp_rows[
          !isf_exp_rows$feature_name %in% already_key, , drop = FALSE]

        if (nrow(isf_exp_rows) > 0L) {
          if (verbose) message("[join_annotation_results] Adding ",
                               nrow(isf_exp_rows),
                               " Tier1c_MS2Fragment_exp child features.")
          match_df <- dplyr::bind_rows(match_df, isf_exp_rows)
        }
      }
    }

  } else {
    match_df$ms2_parent_feature      <- NA_character_
    match_df$ms2_lib_group           <- NA_character_
    match_df$ms2_matched_fragment_mz <- NA_real_
    match_df$ms2_ppm                 <- NA_real_
    match_df$ms2_mda                 <- NA_real_
    match_df$ms2_rt_diff             <- NA_real_
    match_df$ms2_tier                <- NA_character_
    match_df$ms2_annotation_phase    <- NA_character_
  }



  if ("reference_tier" %in% names(match_df)) {
    match_df <- match_df %>%
      dplyr::mutate(
        .p_lib    = .tier_priority(reference_tier),
        .p_ms2    = .tier_priority(ms2_tier),
        .ms2_wins = !is.na(ms2_tier) & (is.na(reference_tier) | .p_ms2 < .p_lib |
                      (!is.na(reference_tier) & .p_ms2 == .p_lib &
                       reference_tier == ms2_tier)),
        .lib_wins = !is.na(reference_tier) & !.ms2_wins,
        reference_tier = dplyr::case_when(
          is.na(reference_tier) & !is.na(ms2_tier) ~ ms2_tier,
          !is.na(reference_tier) & is.na(ms2_tier) ~ reference_tier,
          .p_ms2 < .p_lib                           ~ ms2_tier,
          TRUE                                      ~ reference_tier
        ),
        annotation_phase = dplyr::case_when(
          is.na(reference_tier)                       ~ NA_character_,
          grepl("_exp", reference_tier, fixed = TRUE) ~ "Phase2",
          TRUE                                         ~ "Phase1"
        )
      )


    lib_evidence_cols <- intersect(
      c("grp_anno_match", "grp_anno_ppm", "grp_anno_mda", "grp_anno_confidence",
        "exp_anno_match", "exp_anno_ppm", "exp_anno_mda", "exp_anno_confidence",
        "NL_match"),
      names(match_df)
    )
    ms2_evidence_cols <- intersect(
      c("ms2_tier", "ms2_parent_feature",
        "ms2_matched_fragment_mz", "ms2_ppm", "ms2_mda", "ms2_rt_diff",
        "ms2_annotation_phase"),
      names(match_df)
    )

    if (length(lib_evidence_cols) > 0L) {
      match_df[match_df$.ms2_wins %in% TRUE, lib_evidence_cols] <- NA
    }
    if (length(ms2_evidence_cols) > 0L) {
      match_df[match_df$.lib_wins %in% TRUE, ms2_evidence_cols] <- NA
    }

    match_df <- match_df %>%
      dplyr::select(-.p_lib, -.p_ms2, -.ms2_wins, -.lib_wins)


    if ("Lib_group" %in% names(match_df)) {
      is_exp <- grepl("_exp", match_df$reference_tier, fixed = TRUE)

      match_df$Lib_group_exp <- dplyr::if_else(
        is_exp & !is.na(match_df$Lib_group),
        paste0(match_df$Lib_group, "_exp"),
        NA_character_
      )
      match_df$Lib_group[is_exp] <- NA_character_
    }
  } else {

    match_df$reference_tier <- ifelse(
      !is.na(match_df$ms2_tier),
      match_df$ms2_tier,
      ifelse(
        !is.na(match_df$Ion_type) & grepl("M\\+H|M\\+Na|M\\+K|M\\+NH4",
                                           match_df$Ion_type, ignore.case = TRUE),
        "Tier1a_AnnoMatch",
        "Tier2_Unexplained"
      )
    )
    match_df$annotation_phase <- ifelse(
      !is.na(match_df$Compound_name),
      ifelse(grepl("_exp", match_df$reference_tier, fixed = TRUE),
             "Phase2", "Phase1"),
      NA_character_
    )
  }



  # ---------------------------------------------------------------------------
  # deduplicate match_df: one row per (feature_name, Compound_name, lib-group)
  # ---------------------------------------------------------------------------

  if (!is.null(anno_df) && nrow(match_df) > 0L && "Compound_name" %in% names(match_df)) {

    .mhas <- function(col) col %in% names(match_df)

    match_df <- match_df %>%
      dplyr::mutate(
        .dedup_key = dplyr::coalesce(
          .data$Lib_group,
          if (.mhas("Lib_group_exp"))
            sub("_exp$", "", .data$Lib_group_exp) else NA_character_,
          if (.mhas("ms2_lib_group"))
            .data$ms2_lib_group else NA_character_,
          .data$Compound_name
        )
      )

    match_anno    <- dplyr::filter(match_df, !is.na(.data$Compound_name))
    match_unanno  <- dplyr::filter(match_df,  is.na(.data$Compound_name))

    grp_counts <- match_anno %>%
      dplyr::count(feature_name, Compound_name, .dedup_key) %>%
      dplyr::filter(n > 1L)

    if (nrow(grp_counts) > 0L) {
      needs_dedup <- match_anno %>%
        dplyr::semi_join(grp_counts,
                         by = c("feature_name", "Compound_name", ".dedup_key"))
      clean_rows  <- match_anno %>%
        dplyr::anti_join(grp_counts,
                         by = c("feature_name", "Compound_name", ".dedup_key"))


      needs_dedup <- needs_dedup %>%
        dplyr::group_by(feature_name, Compound_name, .dedup_key) %>%
        dplyr::mutate(
          .best_row_dd = {
            tiers_dd <- reference_tier
            tiers_dd[is.na(tiers_dd)] <- "Tier_ZZZ"
            idx <- which.min(.tier_priority(tiers_dd))
            idx
          },
          .row_in_grp = dplyr::row_number()
        ) %>%
        dplyr::ungroup()

      deduped <- needs_dedup %>%
        dplyr::group_by(feature_name, Compound_name, .dedup_key) %>%
        dplyr::slice(first(.best_row_dd)) %>%
        dplyr::ungroup() %>%
        { if (.mhas("iso_annotation")) {

            iso_pool <- needs_dedup %>%
              dplyr::group_by(feature_name, Compound_name, .dedup_key) %>%
              dplyr::summarise(
                iso_annotation_pool = {
                  vals <- iso_annotation[!is.na(iso_annotation) &
                                         iso_annotation != ""]
                  if (length(vals) > 0L) vals[1L] else NA_character_
                },
                .groups = "drop"
              )
            dplyr::left_join(., iso_pool,
                             by = c("feature_name", "Compound_name",
                                    ".dedup_key")) %>%
              dplyr::mutate(
                iso_annotation = dplyr::coalesce(
                  dplyr::na_if(iso_annotation, ""),
                  iso_annotation_pool)
              ) %>%
              dplyr::select(-iso_annotation_pool)
          } else . } %>%
        dplyr::select(-.best_row_dd, -.row_in_grp)

      match_anno <- dplyr::bind_rows(clean_rows, deduped)
    }

    match_df <- dplyr::bind_rows(match_anno, match_unanno) %>%
      dplyr::select(-.dedup_key)
  }

  # ---------------------------------------------------------------------------
  # iso lookup
  # ---------------------------------------------------------------------------
  iso_ft <- NULL
  if (!is.null(anno_df) && "iso_annotation" %in% names(anno_df)) {
    iso_ft <- anno_df %>%
      dplyr::select(feature_name, iso_annotation) %>%
      dplyr::filter(!is.na(iso_annotation) & iso_annotation != "") %>%
      dplyr::distinct(feature_name, .keep_all = TRUE) %>%
      dplyr::rename(iso_annotation_global = iso_annotation)
  }

  FT_with_iso <- FT
  if (!is.null(iso_ft)) {
    FT_with_iso <- dplyr::left_join(FT, iso_ft, by = "feature_name",
                                   relationship = "many-to-many")
  } else {
    FT_with_iso$iso_annotation_global <- NA_character_
  }

  FT_full <- FT_with_iso %>%
    dplyr::left_join(match_df, by = "feature_name",
                     relationship = "many-to-many") %>%
    dplyr::mutate(
      iso_annotation = dplyr::coalesce(
        if ("iso_annotation" %in% names(.)) .data$iso_annotation
        else NA_character_,
        .data$iso_annotation_global
      )
    ) %>%
    dplyr::select(-iso_annotation_global) %>%
    dplyr::arrange(.data[[rt_col]]) %>%
    dplyr::distinct()

  # ---------------------------------------------------------------------------
  # Annotation column list
  # ---------------------------------------------------------------------------
  feature_id_cols <- c("feature_name", mz_col, rt_col)

  annotation_cols <- c(
    # Library match
    "PI_feature_name", "Lib_group",
    "Compound_name", "Compound_M", "Compound_rt",
    "Compound_peak_mz", "Compound_peak_int",
    "Compound_peak_mz_full", "Compound_peak_int_full",
    "Total_library_peaks", "Ion_type",
    "Matched_peaks", "Match_peak_score",
    "Matched_features_mz", "Matched_features_rt", "Matched_features_int",
    "inner_cosine", "outer_cosine", "intensity_ratio",
    "PI_match", "cutoff_to_MH",
    "grp_anno_match", "grp_anno_ppm", "grp_anno_mda", "grp_anno_confidence",
    "exp_anno_match", "exp_anno_ppm", "exp_anno_mda", "exp_anno_confidence",
    "NL_match", "iso_annotation",
    # MS2 fragment annotation (Tier 1c / 1c_exp)
    "ms2_tier", "ms2_lib_group", "ms2_parent_feature", "ms2_matched_fragment_mz",
    "ms2_ppm", "ms2_mda", "ms2_rt_diff", "ms2_annotation_phase",
    # final resolved tier
    "reference_tier", "annotation_phase",
    # Target list
    "targetList_anno",
    "Comp_mz_targetMatch", "mz_diff_targetMatch", "ppm_error_targetMatch"
  )

  # ── disambiguate same Compound_name at different RT windows ─────────────────
  # The same compound can appear at multiple RT peaks (e.g. isomers,
  # co-eluting species)
  if ("Compound_rt" %in% names(FT_full)) {
    # Store original unsuffixed name as join key BEFORE any modification.
    FT_full$.orig_compound <- FT_full$Compound_name

    rank_map <- FT_full %>%
      dplyr::filter(!is.na(.orig_compound), !is.na(Lib_group)) %>%
      dplyr::select(.orig_compound, Lib_group) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        .lib_idx = as.integer(sub(".*_(\\d+)$", "\\1", Lib_group))
      ) %>%
      dplyr::group_by(.orig_compound) %>%
      dplyr::mutate(
        .n_grp   = dplyr::n_distinct(Lib_group),
        .rt_rank = dplyr::dense_rank(.lib_idx),
        .suffix  = ifelse(.n_grp > 1L, paste0("_", .rt_rank), "")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(.orig_compound, Lib_group, .suffix)


    lib_suffix_map <- rank_map %>%
      dplyr::select(Lib_group, .suffix_lib = .suffix)

    cpd_suffix_default <- rank_map %>%
      dplyr::group_by(.orig_compound) %>%
      dplyr::slice_min(order_by = .suffix, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(.orig_compound, .suffix_default = .suffix)


    FT_full <- FT_full %>%
      dplyr::left_join(rank_map, by = c(".orig_compound", "Lib_group")) %>%
      dplyr::left_join(cpd_suffix_default, by = ".orig_compound") %>%
      dplyr::mutate(
        .suffix = dplyr::coalesce(.suffix, .suffix_default)
      ) %>%
      dplyr::select(-.suffix_default)


    if ("ms2_lib_group" %in% names(FT_full)) {
      lib_suffix_vec <- setNames(lib_suffix_map$.suffix_lib,
                                 lib_suffix_map$Lib_group)
      isf_p2_mask <- !is.na(FT_full$ms2_lib_group) & is.na(FT_full$Lib_group)
      if (any(isf_p2_mask)) {
        FT_full$.suffix[isf_p2_mask] <-
          lib_suffix_vec[FT_full$ms2_lib_group[isf_p2_mask]]
      }
    }

    FT_full <- FT_full %>%
      dplyr::mutate(
        Compound_name = ifelse(!is.na(.suffix) & .suffix != "",
                               paste0(Compound_name, .suffix),
                               Compound_name)
      ) %>%
      dplyr::select(-.suffix, -.orig_compound)
  }

  # ---------------------------------------------------------------------------
  # export_reduced_table
  # ---------------------------------------------------------------------------
  export_reduced_table <- FT_full %>%
    dplyr::select(dplyr::any_of(c(feature_id_cols, annotation_cols))) %>%
    dplyr::distinct() %>%
    data.frame()

  # ---------------------------------------------------------------------------
  # export_reduced_anno_table
  # ---------------------------------------------------------------------------
  export_reduced_anno_table <- export_reduced_table %>%
    dplyr::filter(!is.na(.data$Compound_name)) %>%
    dplyr::distinct() %>%
    data.frame()

  # ---------------------------------------------------------------------------
  # feature_compound_summary
  # ---------------------------------------------------------------------------
  has_col <- function(df, col) col %in% names(df)


  feature_compound_summary <- FT_full %>%
    dplyr::filter(!is.na(.data$Compound_name)) %>%
    dplyr::mutate(

      .grp_key = dplyr::coalesce(
        .data$Lib_group,
        if ("Lib_group_exp" %in% names(.data))
          sub("_exp$", "", .data$Lib_group_exp) else NA_character_,
        if ("ms2_lib_group" %in% names(.data))
          .data$ms2_lib_group else NA_character_,
        .data$Compound_name
      )
    ) %>%
    dplyr::group_by(feature_name, Compound_name, .grp_key) %>%
    dplyr::summarise(
      n_evidence_rows         = dplyr::n(),
      mz                      = dplyr::first(.data[[mz_col]]),
      rt                      = dplyr::first(.data[[rt_col]]),
      intensity               = if (!is.null(int_col))
                                  dplyr::first(.data[[int_col]])
                                else NA_real_,
      inner_cosine       = if (has_col(dplyr::pick(dplyr::everything()),
                                            "inner_cosine"))
                                  safe_max(inner_cosine) else NA_real_,
      outer_cosine       = if (has_col(dplyr::pick(dplyr::everything()),
                                            "outer_cosine"))
                                  safe_max(outer_cosine) else NA_real_,

      .best_row           = {
        tiers_in_grp <- reference_tier
        tiers_in_grp[is.na(tiers_in_grp)] <- "Tier_ZZZ"
        which.min(.tier_priority(tiers_in_grp))
      },
      reference_tier         = reference_tier[.best_row],
      annotation_phase        = if (has_col(dplyr::pick(dplyr::everything()),
                                            "annotation_phase"))
                                  annotation_phase[.best_row]
                                else NA_character_,
      PI_match                = if (has_col(dplyr::pick(dplyr::everything()),
                                            "PI_match"))
                                  PI_match[.best_row]
                                else NA_character_,
      ms2_tier                = {
        v <- if (has_col(dplyr::pick(dplyr::everything()), "ms2_tier"))
               ms2_tier else rep(NA_character_, dplyr::n())
        best <- v[.best_row]
        if (is.na(best) || best == "") {
          nv <- v[!is.na(v) & v != ""]
          if (length(nv) > 0L) nv[1L] else NA_character_
        } else best
      },
      ms2_parent_feature      = {
        v <- if (has_col(dplyr::pick(dplyr::everything()), "ms2_parent_feature"))
               ms2_parent_feature else rep(NA_character_, dplyr::n())
        best <- v[.best_row]
        if (is.na(best) || best == "") {
          nv <- v[!is.na(v) & v != ""]
          if (length(nv) > 0L) nv[1L] else NA_character_
        } else best
      },
      ms2_lib_group           = {
        v <- if (has_col(dplyr::pick(dplyr::everything()), "ms2_lib_group"))
               ms2_lib_group else rep(NA_character_, dplyr::n())
        best <- v[.best_row]
        if (is.na(best) || best == "") {
          nv <- v[!is.na(v) & v != ""]
          if (length(nv) > 0L) nv[1L] else NA_character_
        } else best
      },
      ms2_matched_fragment_mz = {
        v <- if (has_col(dplyr::pick(dplyr::everything()), "ms2_matched_fragment_mz"))
               ms2_matched_fragment_mz else rep(NA_real_, dplyr::n())
        best <- v[.best_row]
        if (is.na(best)) {
          nv <- v[!is.na(v)]
          if (length(nv) > 0L) nv[1L] else NA_real_
        } else best
      },
      ms2_ppm                 = {
        v <- if (has_col(dplyr::pick(dplyr::everything()), "ms2_ppm"))
               ms2_ppm else rep(NA_real_, dplyr::n())
        best <- v[.best_row]
        if (is.na(best)) {
          nv <- v[!is.na(v)]
          if (length(nv) > 0L) nv[1L] else NA_real_
        } else best
      },
      # all remaining evidence columns indexed from same .best_row
      targetList_anno         = if (has_col(dplyr::pick(dplyr::everything()),
                                            "targetList_anno"))
                                  targetList_anno[.best_row]
                                else NA_character_,
      Comp_mz_targetMatch     = if (has_col(dplyr::pick(dplyr::everything()),
                                            "Comp_mz_targetMatch"))
                                  Comp_mz_targetMatch[.best_row]
                                else NA_real_,
      mz_diff_targetMatch     = if (has_col(dplyr::pick(dplyr::everything()),
                                            "mz_diff_targetMatch"))
                                  mz_diff_targetMatch[.best_row]
                                else NA_real_,
      ppm_error_targetMatch   = if (has_col(dplyr::pick(dplyr::everything()),
                                            "ppm_error_targetMatch"))
                                  ppm_error_targetMatch[.best_row]
                                else NA_real_,
      grp_anno_match        = if (has_col(dplyr::pick(dplyr::everything()),
                                            "grp_anno_match"))
                                  grp_anno_match[.best_row]
                                else NA_character_,
      grp_anno_ppm          = if (has_col(dplyr::pick(dplyr::everything()),
                                            "grp_anno_ppm"))
                                  grp_anno_ppm[.best_row]
                                else NA_character_,
      grp_anno_mda          = if (has_col(dplyr::pick(dplyr::everything()),
                                            "grp_anno_mda"))
                                  grp_anno_mda[.best_row]
                                else NA_character_,
      grp_anno_confidence   = if (has_col(dplyr::pick(dplyr::everything()),
                                            "grp_anno_confidence"))
                                  grp_anno_confidence[.best_row]
                                else NA_character_,
      exp_anno_match        = if (has_col(dplyr::pick(dplyr::everything()),
                                            "exp_anno_match"))
                                  exp_anno_match[.best_row]
                                else NA_character_,
      exp_anno_ppm          = if (has_col(dplyr::pick(dplyr::everything()),
                                            "exp_anno_ppm"))
                                  exp_anno_ppm[.best_row]
                                else NA_character_,
      exp_anno_mda          = if (has_col(dplyr::pick(dplyr::everything()),
                                            "exp_anno_mda"))
                                  exp_anno_mda[.best_row]
                                else NA_character_,
      exp_anno_confidence   = if (has_col(dplyr::pick(dplyr::everything()),
                                            "exp_anno_confidence"))
                                  exp_anno_confidence[.best_row]
                                else NA_character_,

      NL_match                = if (has_col(dplyr::pick(dplyr::everything()),
                                            "NL_match"))
                                  NL_match[.best_row]
                                else NA_character_,

      iso_annotation          = if (has_col(dplyr::pick(dplyr::everything()),
                                            "iso_annotation")) {
                                  best_iso <- iso_annotation[.best_row]
                                  if (is.na(best_iso) || best_iso == "") {
                                    all_iso <- iso_annotation[!is.na(iso_annotation) &
                                                              iso_annotation != ""]
                                    if (length(all_iso) > 0L) all_iso[1L] else NA_character_
                                  } else best_iso
                                } else NA_character_,
      cutoff_to_MH            = if (has_col(dplyr::pick(dplyr::everything()),
                                            "cutoff_to_MH"))
                                  cutoff_to_MH[.best_row]
                                else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::select(-.best_row) %>%

    dplyr::mutate(Lib_group = ifelse(
      grepl("_exp", reference_tier, fixed = TRUE), NA_character_, .grp_key
    )) %>%
    dplyr::select(-.grp_key)


  if ("iso_annotation" %in% names(match_df) &&
      "mz"            %in% names(match_df) &&
      "Lib_group"     %in% names(match_df)) {

    iso_threshold <- 1.003355

    anchors_iso <- match_df %>%
      dplyr::filter(!is.na(iso_annotation),
                    grepl("\\[M\\]", iso_annotation),
                    !grepl("_exp", iso_annotation),
                    !is.na(Lib_group)) %>%
      dplyr::select(feature_name, mz, Lib_group, iso_annotation) %>%
      dplyr::mutate(
        anchor_grp = sub(".*\\[grp[0-9]+_([0-9]+)\\].*", "\\1", iso_annotation),
        lib_num    = sub(".*_([0-9]+)$", "\1", Lib_group)
      )

    relabel <- function(iso_ann, mz_feat, lib_group_feat) {
      if (is.na(iso_ann) || iso_ann == "") return(iso_ann)
      m <- regmatches(iso_ann,
             regexec("^(\\[grp([0-9]+)_([0-9]+)\\])\\[M\\]([0-9]*)\\+_exp$",
                     iso_ann))[[1]]
      if (length(m) == 0) return(iso_ann)
      c_val <- if (nchar(m[5]) == 0) 1L else as.integer(m[5])
      if (c_val < 2L) return(iso_ann)

      lib_num_feat <- sub(".*_([0-9]+)$", "\1", if (is.null(lib_group_feat) || is.na(lib_group_feat)) "" else lib_group_feat)
      if (is.na(lib_num_feat) || lib_num_feat == "") return(iso_ann)

      cands <- anchors_iso[anchors_iso$lib_num == lib_num_feat, ]
      if (nrow(cands) == 0) return(iso_ann)

      for (i in seq_len(nrow(cands))) {
        mz_diff <- mz_feat - cands$mz[i]
        for (mm in 1:4) {
          for (cc in 2:4) {
            expected <- mm * iso_threshold / cc
            if (abs(mz_diff - expected) <= 0.005) {
              c_str <- if (cc > 1) as.character(cc) else ""
              return(paste0("[grp", lib_num_feat, "_", cands$anchor_grp[i],
                            "][M+", mm, "]", c_str, "+_exp"))
            }
          }
        }
      }
      return(iso_ann)
    }

    needs_check <- !is.na(match_df$iso_annotation) &
                   grepl("\\[M\\][0-9]+\\+_exp$", match_df$iso_annotation)

    if (any(needs_check)) {
      match_df$iso_annotation[needs_check] <- mapply(
        relabel,
        match_df$iso_annotation[needs_check],
        match_df$mz[needs_check],
        match_df$Lib_group[needs_check],
        USE.NAMES = FALSE
      )
    }
  }

  match_df <- match_df %>%
    dplyr::arrange(rt) %>%
    dplyr::distinct()

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  list(
    export_full_table               = FT_full,
    export_reduced_table            = export_reduced_table,
    export_reduced_anno_table       = export_reduced_anno_table,
    export_feature_compound_summary = feature_compound_summary,
    export_match_table              = match_df
  )
}
