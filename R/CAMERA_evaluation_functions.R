# =============================================================================
# CAMERA_evaluation_functions.R
#
# Functions for compound-centric CAMERA grouping evaluation.
#
# Functions (in order):
#   1. build_eval_full()        — join tier annotations of MS1Match to CAMERA FT
#   2. build_anchors()          — define [M+H]+ anchor per compound instance
#   3. build_pcgroup_purity()   — assess over-grouping per pcgroup
#   4. build_plot2_data()       — per-pcgroup purity summary for plotting
#   5. build_compound_eval()    — tier-wise grouping recall (wide table format)
#   6. build_anno_accuracy()    — Tier 1a adduct annotation accuracy
#   7. build_benchmark_table()  — aggregate benchmark summary table
#
# See also: inspect_compound.R for per-compound manual inspection.
# See also: CAMERA_evaluation_analysis.Rmd for the full analysis MS1Match.
#
# Dependencies: dplyr, tidyr, stringr
# =============================================================================



# =============================================================================
# classify_ref_status
# =============================================================================

#' Classify features by reference-relative evaluation status
#'
#' @description
#' Evaluation status is defined \emph{relative to the compound-centric
#' reference feature set}, not relative to an assumed ground truth. For each
#' compound, the Tier 1a \code{[M+H]+} feature is the anchor and defines the
#' compound's expected group.
#'
#' \describe{
#'   \item{\strong{Ref-confirmed}}{Feature is assigned to the same group as
#'     its compound's anchor.}
#'   \item{\strong{Ref-conflicting}}{Feature is assigned to a different group
#'     than its compound's anchor.}
#'   \item{\strong{Ref-unmatched, co-grouped}}{Feature has no reference entry
#'     but shares a group with a compound's anchor (possible over-grouping).}
#'   \item{\strong{Ref-unmatched, unrelated}}{Feature has no reference entry,
#'     and its group contains no compound anchor.}
#'   \item{\strong{Ref-unassigned}}{Feature has a reference entry, but the
#'     tool assigned it no group. Zero for CAMERA, which groups every feature.}
#' }
#'
#' Rows that fit none of these statuses (e.g. a reference feature whose
#' compound has no resolvable anchor) get \code{NA} and trigger a warning.
#'
#' @param eval_full \code{data.frame}. Must contain \code{Lib_group},
#'   \code{pcgroup}, \code{CAMERA_grouped}, \code{CAMERA_split}; ideally
#'   \code{anchor_pcgroup}.
#' @param anchors Optional \code{data.frame} (output of
#'   \code{\link{build_anchors}}). Used to define the set of anchor groups for
#'   splitting Ref-unmatched. If \code{NULL}, \code{eval_full$anchor_pcgroup}
#'   is used.
#'
#' @return \code{eval_full} with an ordered factor column \code{ref_status}.
#'
#' @importFrom dplyr mutate case_when
#' @export
classify_ref_status <- function(eval_full, anchors = NULL) {

  stopifnot(is.data.frame(eval_full))

  required <- c("Lib_group", "pcgroup")
  missing  <- setdiff(required, names(eval_full))
  if (length(missing) > 0L)
    stop("classify_ref_status: missing required columns: ",
         paste(missing, collapse = ", "))

  n <- nrow(eval_full)
  grouped <- if ("CAMERA_grouped" %in% names(eval_full))
    eval_full$CAMERA_grouped %in% TRUE else rep(FALSE, n)
  split   <- if ("CAMERA_split" %in% names(eval_full))
    eval_full$CAMERA_split %in% TRUE else rep(FALSE, n)

  # Set of groups that contain a compound anchor
  anchor_set <- if (!is.null(anchors) && "anchor_pcgroup" %in% names(anchors)) {
    unique(stats::na.omit(anchors$anchor_pcgroup))
  } else if ("anchor_pcgroup" %in% names(eval_full)) {
    unique(stats::na.omit(eval_full$anchor_pcgroup))
  } else {
    NULL
  }
  if (is.null(anchor_set))
    warning("classify_ref_status: no anchor_pcgroup available; ",
            "all Ref-unmatched features are labelled 'unrelated'.")
  in_anchor_grp <- !is.na(eval_full$pcgroup) &
                   eval_full$pcgroup %in% anchor_set

  # Reference entry: Phase 1 Lib_group, or Phase 2 feature of a known compound
  is_phase2 <- if (all(c("Compound_name", "annotation_phase") %in% names(eval_full)))
    !is.na(eval_full$Compound_name) &
      grepl("Phase2", eval_full$annotation_phase, ignore.case = TRUE)
  else rep(FALSE, n)
  has_ref <- !is.na(eval_full$Lib_group) | is_phase2

  status_levels <- c("Ref-confirmed", "Ref-conflicting",
                     "Ref-unmatched, co-grouped", "Ref-unmatched, unrelated",
                     "Ref-unassigned")

  ref_status <- dplyr::case_when(
    has_ref & grouped                                  ~ "Ref-confirmed",
    has_ref & split                                    ~ "Ref-conflicting",
    has_ref & is.na(eval_full$pcgroup)                 ~ "Ref-unassigned",
    !has_ref & !is.na(eval_full$pcgroup) & in_anchor_grp  ~ "Ref-unmatched, co-grouped",
    !has_ref & !is.na(eval_full$pcgroup)               ~ "Ref-unmatched, unrelated",
    TRUE                                               ~ NA_character_
  )

  out <- eval_full
  out$ref_status <- factor(ref_status, levels = status_levels, ordered = TRUE)

  n_na <- sum(is.na(out$ref_status))
  if (n_na > 0L)
    warning("classify_ref_status: ", n_na, " rows could not be classified ",
            "(e.g. reference features whose compound has no resolvable anchor).")

  out
}


# =============================================================================
# build_eval_full
# =============================================================================

#' Build the feature × compound evaluation table
#'
#' @description
#' Joins the MS1Match's tier-annotated feature-compound summary onto the CAMERA
#' annotated feature table, so every row carries both the MS1Match annotation
#' tier and the CAMERA pcgroup / adduct assignments.
#'
#' @param feature_compound_summary \code{data.frame}. The
#'   \code{export_feature_compound_summary} element of \code{join_annotation_results()}.
#'   Must contain: \code{feature_name}, \code{Compound_name},
#'   \code{reference_tier}, \code{annotation_phase}, \code{Lib_group},
#'   \code{adduct_match}, \code{rt}.
#' @param CAMERA_ft \code{data.frame}. CAMERA-annotated feature table (e.g.
#'   \code{CAMERA_annotated_pool_FT_sub}). Must contain: \code{feature_name},
#'   \code{pcgroup}, \code{adduct} (CAMERA's adduct/isotope string), \code{mz}.
#'
#' @return A \code{data.frame} with one row per feature × compound combination,
#'   containing all MS1Match annotation columns plus CAMERA's \code{pcgroup},
#'   \code{adduct} (renamed \code{CAMERA_adduct}), and \code{mz}.
#'
#' @importFrom dplyr left_join rename select
#' @export
build_eval_full <- function(feature_compound_summary, CAMERA_ft) {

  stopifnot(
    is.data.frame(feature_compound_summary),
    is.data.frame(CAMERA_ft),
    "feature_name" %in% names(feature_compound_summary),
    "feature_name" %in% names(CAMERA_ft),
    "pcgroup"      %in% names(CAMERA_ft)
  )

  CAMERA_cols <- CAMERA_ft %>%
    dplyr::select(
      feature_name,
      pcgroup,
      CAMERA_mz = mz,
      dplyr::any_of(c("adduct", "isotopes", "rt"))
    ) %>%
    dplyr::rename_with(~"CAMERA_adduct", dplyr::any_of("adduct")) %>%
    dplyr::rename_with(~"CAMERA_rt",     dplyr::any_of("rt"))

  eval_full <- feature_compound_summary %>%
    dplyr::left_join(CAMERA_cols, by = "feature_name")

  eval_full
}


# =============================================================================
# build_anchors
# =============================================================================

#' Define the [M+H]+ anchor feature for each compound instance
#'
#'  @description
#'  For each \code{Lib_group} (compound instance), selects the single Tier 1a
#' \code{[M+H]+} feature as the reference anchor. Where multiple Tier 1a
#' \code{[M+H]+} features exist for one instance, the one with the highest
#' intensity (or first by feature name) is kept.
#'
#' @param eval_full \code{data.frame}. Output of \code{\link{build_eval_full}}.
#'   Must contain: \code{Lib_group}, \code{reference_tier},
#'   \code{adduct_match}, \code{pcgroup}, \code{feature_name}.
#'
#' @return A \code{data.frame} with one row per \code{Lib_group}, containing:
#'   \code{Lib_group}, \code{Compound_name}, \code{anchor_feature},
#'   \code{anchor_pcgroup}, \code{anchor_mz}, \code{anchor_rt}.
#'
#' @importFrom dplyr filter group_by slice_head ungroup select rename
#' @export
build_anchors <- function(eval_full, full_anno = NULL) {


  stopifnot(is.data.frame(eval_full))

  # Anchor source: full_anno (all tier rows per feature) if given, else eval_full
  anchor_src <- if (!is.null(full_anno)) {
    stopifnot(is.data.frame(full_anno))
    full_anno
  } else {
    eval_full
  }

  required <- c("Lib_group", "reference_tier", "feature_name")
  missing  <- setdiff(required, names(anchor_src))
  if (length(missing) > 0L)
    stop("build_anchors: anchor source is missing required columns: ",
         paste(missing, collapse = ", "),
         "\n  Available columns: ", paste(names(anchor_src), collapse = ", "))

  # check mz
  if (!any(c("mz", "mzmed", "mz.x", "CAMERA_mz") %in% names(anchor_src)))
    stop("build_anchors: cannot find any mz column in anchor source.\n",
         "  Available columns: ", paste(names(anchor_src), collapse = ", "))

  # check adduct column
  adduct_col <- intersect(c("adduct_match", "grp_anno_match"), names(anchor_src))[1L]
  if (is.na(adduct_col))
    stop("build_anchors: cannot find adduct column in anchor source. ",
         "Expected 'adduct_match' or 'grp_anno_match'.\n",
         "  Available columns: ", paste(names(anchor_src), collapse = ", "))

  # check pcgroup
  pcgroup_col <- intersect(c("pcgroup", "CAMERA_pcgroup"), names(anchor_src))[1L]
  needs_pcgroup_join <- is.na(pcgroup_col)

  # tier column
  tier_col <- intersect(c("reference_tier", "tier", "Tier"), names(anchor_src))[1L]
  if (is.na(tier_col))
    stop("build_anchors: cannot find tier column in anchor source. ",
         "Expected 'reference_tier'.\n",
         "  Available columns: ", paste(names(anchor_src), collapse = ", "))

  anchors <- anchor_src %>%
    dplyr::filter(
      grepl("^Tier1a_", .data[[tier_col]]) &
        !grepl("exp", .data[[tier_col]], ignore.case = TRUE),
      grepl("^\\[M\\+H\\]\\+", .data[[adduct_col]])
    ) %>%
    dplyr::group_by(Lib_group) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      .mz_use = dplyr::coalesce(
        if ("mz"    %in% names(.)) mz    else NA_real_,
        if ("mzmed" %in% names(.)) mzmed else NA_real_,
        if ("mz.x"  %in% names(.)) mz.x  else NA_real_,
        if ("CAMERA_mz" %in% names(.)) CAMERA_mz else NA_real_
      ),
      .rt_use = dplyr::coalesce(
        if ("rt"    %in% names(.)) rt    else NA_real_,
        if ("rtmed" %in% names(.)) rtmed else NA_real_,
        if ("rt.x"  %in% names(.)) rt.x  else NA_real_
      )
    ) %>%
    dplyr::select(
      Lib_group,
      Compound_name,
      anchor_feature = feature_name,
      dplyr::any_of(c("pcgroup", "CAMERA_pcgroup")),
      anchor_mz      = .mz_use,
      anchor_rt      = .rt_use
    ) %>%
    { df <- .
      if ("CAMERA_pcgroup" %in% names(df))
        dplyr::rename(df, anchor_pcgroup = CAMERA_pcgroup)
      else if ("pcgroup" %in% names(df))
        dplyr::rename(df, anchor_pcgroup = pcgroup)
      else
        df
    } %>%
    # If anchor_src had no pcgroup column,
    # join anchor_pcgroup from eval_full via feature_name.
    { if (needs_pcgroup_join)
        dplyr::left_join(.,
          eval_full %>%
            dplyr::select(feature_name, pcgroup) %>%
            dplyr::distinct() %>%
            dplyr::rename(anchor_feature = feature_name,
                          anchor_pcgroup = pcgroup),
          by = "anchor_feature")
      else . }

  if (nrow(anchors) == 0L)
    warning("build_anchors: no Tier1a [M+H]+ features found. ",
            "Check that reference_tier values match 'Tier1a' exactly and ",
            "that adduct_match contains '[M+H]+'.")

  anchors
}


# =============================================================================
# build_pcgroup_purity
# =============================================================================

#' Assess CAMERA pcgroup purity (over-grouping)
#'
#' @description
#' For every pcgroup that contains at least one anchor feature, counts how many
#' distinct compounds (neutral masses) were merged into the same pcgroup by
#' CAMERA. A pure pcgroup contains exactly one compound.
#'
#' @param CAMERA_ft \code{data.frame}. CAMERA-annotated feature table. Must
#'   contain: \code{feature_name}, \code{pcgroup}.
#' @param eval_full \code{data.frame}. Output of \code{\link{build_eval_full}}.
#' @param anchors \code{data.frame}. Output of \code{\link{build_anchors}}.
#'
#' @return A named \code{list} with two elements:
#' \describe{
#'   \item{\code{all_pcgroups}}{Per-pcgroup summary: \code{pcgroup},
#'     \code{pcgroup_n_features}, \code{pcgroup_n_compounds},
#'     \code{pcgroup_outside_ref}, \code{median_rt}.}
#'   \item{\code{anchor_pcgroups}}{Subset of \code{all_pcgroups} restricted to
#'     pcgroups that contain at least one anchor feature.}
#' }
#'
#' @importFrom dplyr left_join group_by summarise n n_distinct filter ungroup
#' @export
build_pcgroup_purity <- function(CAMERA_ft, eval_full, anchors) {

  stopifnot(
    is.data.frame(CAMERA_ft),
    is.data.frame(eval_full),
    is.data.frame(anchors),
    "pcgroup"      %in% names(CAMERA_ft),
    "feature_name" %in% names(CAMERA_ft),
    "pcgroup"      %in% names(eval_full),
    "anchor_pcgroup" %in% names(anchors)
  )

  # ── MS1Match: feature -> effective Lib_group ─────────────────────────

  feat_compound <- eval_full %>%
    dplyr::mutate(
      .eff_lib = dplyr::coalesce(
        .data$Lib_group,
        if ("ms2_lib_group" %in% names(.)) .data$ms2_lib_group else NA_character_
      )
    ) %>%
    dplyr::select(feature_name, .eff_lib) %>%
    dplyr::distinct()

  ft_annotated <- CAMERA_ft %>%
    dplyr::left_join(feat_compound, by = "feature_name")

  # ── CAMERA: distinct neutral masses per pcgroup ─────────────────────
  # purity = n_distinct_neutralM: the number of distinct neutral masses that
  # CAMERA assigns within the pcgroup.
  # purity = 1 - CAMERA agrees all features share one compound mass.
  # purity > 1 - CAMERA assigned multiple neutral masses - over-grouping.
  # Features with no CAMERA adduct annotation do not contribute a neutral
  # mass and are excluded from this count.
  #
  # Neutral mass extraction mirrors build_plot2_data: parse the numeric value
  # at the end of each CAMERA adduct string, e.g. "[M+H]+ 152.047" - 152.047.
  adduct_col <- intersect(c("adduct", "CAMERA_adduct"), names(CAMERA_ft))[1L]

  .extract_nm <- function(x) {
    entries <- trimws(unlist(strsplit(x, "(?<=\\S)\\s+(?=\\[)", perl = TRUE)))
    masses  <- regmatches(entries, regexpr("[0-9]+\\.[0-9]+$", entries))
    as.numeric(masses)
  }

  if (!is.na(adduct_col)) {
    ft_adduct <- CAMERA_ft %>%
      dplyr::filter(!is.na(.data[[adduct_col]]),
                    nchar(trimws(.data[[adduct_col]])) > 0L) %>%
      dplyr::mutate(
        .nm_list = lapply(.data[[adduct_col]], .extract_nm)
      )

    ft_expanded <- do.call(rbind, lapply(seq_len(nrow(ft_adduct)), function(i) {
      masses <- ft_adduct$.nm_list[[i]]
      if (length(masses) == 0L) return(NULL)
      data.frame(pcgroup      = ft_adduct$pcgroup[i],
                 neutral_mass = round(masses, 4L),
                 stringsAsFactors = FALSE)
    }))

    CAMERA_purity <- if (!is.null(ft_expanded) && nrow(ft_expanded) > 0L) {
      ft_expanded %>%
        dplyr::group_by(pcgroup) %>%
        dplyr::summarise(
          pcgroup_n_compounds = dplyr::n_distinct(neutral_mass),
          .groups = "drop"
        )
    } else {
      data.frame(pcgroup = integer(0), pcgroup_n_compounds = integer(0))
    }
  } else {
    CAMERA_purity <- ft_annotated %>%
      dplyr::group_by(pcgroup) %>%
      dplyr::summarise(
        pcgroup_n_compounds = dplyr::n_distinct(.eff_lib, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # ── per-pcgroup summary ─────────────────────────────────────────
  all_pcgroups <- ft_annotated %>%
    dplyr::group_by(pcgroup) %>%
    dplyr::summarise(
      pcgroup_n_features  = dplyr::n(),
      pcgroup_outside_ref = sum(is.na(.eff_lib)),
      median_rt           = median(rt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(CAMERA_purity, by = "pcgroup") %>%
    dplyr::mutate(
      # Pcgroups with no CAMERA adduct annotation, purity = NA
      pcgroup_n_compounds = dplyr::coalesce(pcgroup_n_compounds, NA_integer_)
    )

  anchor_pcgroups <- all_pcgroups %>%
    dplyr::filter(pcgroup %in% anchors$anchor_pcgroup)

  list(
    all_pcgroups    = all_pcgroups,
    anchor_pcgroups = anchor_pcgroups
  )
}


# =============================================================================
#  build_plot2_data
# =============================================================================

#' Build per-pcgroup purity summary for the over-grouping scatter plot
#'
#' @description
#' Extracts neutral masses from CAMERA's adduct annotation column for every
#' feature, then summarises per pcgroup: how many distinct neutral masses
#' were assigned by CAMERA within that pcgroup, and what is the median RT.
#' The result is the \code{plot2_data} object consumed by the
#' \code{plot-over-grouping} chunk in \code{CAMERA_evaluation_analysis.Rmd}.
#'
#' CAMERA writes adduct annotations as e.g. \code{"[M+H]+ 194.057909"} or
#' \code{"[M+Na]+ 194.057909; [M+K]+ 194.057909"}. This function extracts
#' the trailing numeric mass (the last decimal number on each semicolon-
#' delimited entry), rounds to 4 dp to absorb floating point noise, and
#' counts distinct values per pcgroup.
#'
#' @param CAMERA_ft \code{data.frame}. CAMERA-annotated feature table. Must
#'   contain: \code{pcgroup}, \code{rt}, and a CAMERA adduct column. The adduct
#'   column is detected automatically: the first column whose name matches
#'   \code{adduct} (case-insensitive) is used, falling back to a column named
#'   \code{CAMERA_adduct}. Rows with missing or empty adduct strings are excluded.
#' @param mass_digits \code{integer(1)}. Number of decimal places used when
#'   rounding neutral masses before counting distinct values. Default \code{4}.
#'
#' @return A \code{data.frame} with one row per pcgroup (those with at least
#'   one CAMERA adduct annotation), containing:
#' \describe{
#'   \item{\code{pcgroup}}{CAMERA pcgroup identifier.}
#'   \item{\code{median_rt}}{Median RT (seconds) of all features in the pcgroup.}
#'   \item{\code{n_distinct_neutralM}}{Number of distinct neutral masses found
#'     within the pcgroup.}
#'   \item{\code{grouping_class}}{Ordered factor:
#'     \code{"1 neutral mass"}, \code{"2 neutral masses"},
#'     \code{"3-5 neutral masses"}, \code{">5 neutral masses"}.}
#' }
#'
#' @importFrom dplyr filter mutate group_by summarise n_distinct ungroup
#' @importFrom stringr str_split str_extract
#' @export
build_plot2_data <- function(CAMERA_ft, mass_digits = 4L) {

  stopifnot(
    is.data.frame(CAMERA_ft),
    "pcgroup" %in% names(CAMERA_ft),
    "rt"      %in% names(CAMERA_ft)
  )

  # ── get adduct column ─────────────────────────────────────────────────
  adduct_col <- names(CAMERA_ft)[
    grepl("^adduct$", names(CAMERA_ft), ignore.case = TRUE)
  ]
  if (length(adduct_col) == 0L) {
    adduct_col <- names(CAMERA_ft)[
      grepl("CAMERA_adduct", names(CAMERA_ft), ignore.case = TRUE)
    ]
  }
  if (length(adduct_col) == 0L) {
    stop(
      "build_plot2_data: could not find an adduct column in CAMERA_ft. ",
      "Expected a column named 'adduct' or 'CAMERA_adduct'."
    )
  }
  adduct_col <- adduct_col[1L]

  # ── Helper: extract neutral mass from one adduct entry ───────────────────
  # CAMERA format:e.g. "[M+H]+ 194.057909 [M+Na]+ 194.057909", "[M+H]+ 194.057909"
  .extract_nm <- function(x) {

    entries <- trimws(unlist(strsplit(x, "(?<=\\S)\\s+(?=\\[)",perl = TRUE))) #
    masses  <- stringr::str_extract(entries, "[0-9]+\\.[0-9]+$")
    as.numeric(masses[!is.na(masses)])
  }

  # ── Explode to one row per adduct entry, extract neutral mass ────────────
  ft_adduct <- CAMERA_ft %>%
    dplyr::filter(!is.na(.data[[adduct_col]]),
                  nchar(trimws(.data[[adduct_col]])) > 0L) %>%
    dplyr::mutate(
      neutral_mass_list = lapply(.data[[adduct_col]], .extract_nm)
    )

  ft_expanded <- do.call(rbind, lapply(seq_len(nrow(ft_adduct)), function(i) {
    masses <- ft_adduct$neutral_mass_list[[i]]
    if (length(masses) == 0L) return(NULL)
    data.frame(
      pcgroup      = ft_adduct$pcgroup[i],
      rt           = ft_adduct$rt[i],
      neutral_mass = masses,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(ft_expanded) || nrow(ft_expanded) == 0L) {
    stop("build_plot2_data: no neutral masses could be extracted. ",
         "Check that the adduct column contains entries like '[M+H]+ 194.057909'.")
  }

  ft_expanded$neutral_mass_r <- round(ft_expanded$neutral_mass, mass_digits)

  # ── summary per pcgroup ────────────────────────────────────────────────
  grouping_levels <- c("1 neutral mass", "2 neutral masses",
                       "3-5 neutral masses", ">5 neutral masses")

  plot2_data <- ft_expanded %>%
    dplyr::group_by(pcgroup) %>%
    dplyr::summarise(
      median_rt           = median(rt, na.rm = TRUE),
      n_distinct_neutralM = dplyr::n_distinct(neutral_mass_r),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      grouping_class = dplyr::case_when(
        n_distinct_neutralM == 1L ~ "1 neutral mass",
        n_distinct_neutralM == 2L ~ "2 neutral masses",
        n_distinct_neutralM <= 5L ~ "3-5 neutral masses",
        TRUE                      ~ ">5 neutral masses"
      ),
      grouping_class = factor(grouping_class, levels = grouping_levels)
    )

  plot2_data
}


# =============================================================================
# 5. build_compound_eval
# =============================================================================

#' Build the per-compound tier-wise grouping recall table
#' @description
#' For each compound instance × annotation tier, counts how many features were
#' co-grouped with the \code{[M+H]+} anchor (recall) versus split to a
#' different pcgroup.
#'
#' @param eval_full \code{data.frame}. Output of \code{\link{build_eval_full}}.
#' @param anchors \code{data.frame}. Output of \code{\link{build_anchors}}.
#' @param purity \code{list}. Output of \code{\link{build_pcgroup_purity}}.
#'
#' @return A wide \code{data.frame} with one row per \code{Lib_group}, with
#'   columns \code{n_MS1Match_<Tier>}, \code{n_CAMERA_grouped_<Tier>},
#'   \code{n_CAMERA_split_<Tier>} for each tier.
#'
#' @importFrom dplyr left_join mutate group_by summarise ungroup filter
#' @importFrom tidyr pivot_wider
#' @export
build_compound_eval <- function(eval_full, anchors, purity) {

  stopifnot(
    is.data.frame(eval_full),
    is.data.frame(anchors),
    is.list(purity)
  )

  tier_simple <- function(tier) {
    dplyr::case_when(
      grepl("1a.*exp|1a_exp", tier, ignore.case = TRUE) ~ "Tier1a_exp",
      grepl("1a",             tier, ignore.case = TRUE) ~ "Tier1a",
      grepl("1b.*exp|1b_exp", tier, ignore.case = TRUE) ~ "Tier1b_exp",
      grepl("1b",             tier, ignore.case = TRUE) ~ "Tier1b",
      grepl("1c.*exp|1c_exp|MS2Fragment_exp", tier, ignore.case = TRUE) ~ "Tier1c_exp",
      grepl("1c",             tier, ignore.case = TRUE) ~ "Tier1c",
      grepl("1d",             tier, ignore.case = TRUE) ~ "Tier1d",
      grepl("Tier2|tier2",    tier, ignore.case = TRUE) ~ "Tier2",
      TRUE                                              ~ "Other"
    )
  }


  anchor_by_lib  <- anchors %>%
    dplyr::select(Lib_group, anchor_pcgroup, anchor_rt)

  anchor_by_cpd  <- anchors %>%
    dplyr::select(Compound_name, anchor_pcgroup, anchor_rt) %>%
    dplyr::distinct()

  long <- eval_full %>%
    # Step 1: join by Lib_group for Phase 1 features
    dplyr::left_join(anchor_by_lib, by = "Lib_group") %>%
    # Step 2: for features where Lib_group join gave NA anchor_pcgroup, join by Compound_name
    dplyr::left_join(
      anchor_by_cpd %>%
        dplyr::rename(anchor_pcgroup_cpd = anchor_pcgroup,
                      anchor_rt_cpd      = anchor_rt),
      by = "Compound_name",
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      .use_cpd_anchor = is.na(anchor_pcgroup) & !is.na(anchor_pcgroup_cpd),
      .rt_diff_cpd = abs(rt - anchor_rt_cpd)
    ) %>%
    dplyr::group_by(feature_name, Compound_name, Lib_group) %>%
    dplyr::slice_min(order_by = .rt_diff_cpd, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      anchor_pcgroup = dplyr::if_else(
        .use_cpd_anchor, anchor_pcgroup_cpd, anchor_pcgroup
      )
    ) %>%
    dplyr::select(-.use_cpd_anchor, -.rt_diff_cpd,
                  -anchor_rt, -anchor_pcgroup_cpd, -anchor_rt_cpd) %>%
    dplyr::mutate(
      tier_simple      = tier_simple(reference_tier),
      CAMERA_grouped   = !is.na(pcgroup) & !is.na(anchor_pcgroup) &
                           pcgroup == anchor_pcgroup,
      CAMERA_split     = !is.na(pcgroup) & !is.na(anchor_pcgroup) &
                           pcgroup != anchor_pcgroup,
      CAMERA_unassigned = is.na(pcgroup)
    )

  # Per compound × tier

  long <- long %>%
    dplyr::mutate(
      .eff_lib = dplyr::coalesce(
        .data$Lib_group,
        if ("ms2_lib_group" %in% names(.)) .data$ms2_lib_group else NA_character_
      )
    )

  counts <- long %>%
    dplyr::filter(tier_simple != "Other") %>%
    dplyr::group_by(.eff_lib, Compound_name, tier_simple) %>%
    dplyr::summarise(
      n_MS1Match         = dplyr::n(),
      n_CAMERA_grouped   = sum(CAMERA_grouped,   na.rm = TRUE),
      n_CAMERA_split     = sum(CAMERA_split,     na.rm = TRUE),
      n_CAMERA_unassigned = sum(CAMERA_unassigned, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(Lib_group = .eff_lib)

  totals <- long %>%
    dplyr::group_by(.eff_lib, Compound_name) %>%
    dplyr::summarise(
      tier_simple        = "TOTAL",
      n_MS1Match         = dplyr::n(),
      n_CAMERA_grouped   = sum(CAMERA_grouped,   na.rm = TRUE),
      n_CAMERA_split     = sum(CAMERA_split,     na.rm = TRUE),
      n_CAMERA_unassigned = sum(CAMERA_unassigned, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(Lib_group = .eff_lib)

  all_counts <- dplyr::bind_rows(counts, totals)

  # Pivot to wide
  compound_eval_wide <- all_counts %>%
    tidyr::pivot_wider(
      id_cols     = c(Lib_group, Compound_name),
      names_from  = tier_simple,
      values_from = c(n_MS1Match, n_CAMERA_grouped, n_CAMERA_split,
                      n_CAMERA_unassigned),
      values_fill = 0L
    )

  compound_eval_wide
}


# =============================================================================
#  build_anno_accuracy
# =============================================================================

#' Assess Tier 1a adduct annotation accuracy
#'
#' For each Tier 1a (and optionally Tier 1a_exp) feature, checks whether
#' CAMERA's adduct neutral mass assignment matches the MS1Match's annotated neutral mass within
#' a given ppm tolerance.
#'
#'
#' 1. **Tier-aware adduct column**: Tier1a features are evaluated against
#'    \code{grp_anno_match}; Tier1a_exp features against
#'    \code{exp_anno_match}.
#'
#' 2. **Multi-annotation aware**: CAMERA sometimes assigns multiple adduct
#'    annotations to one feature in a single string, e.g.
#'    \code{"[M+K+HCOOH]+ 63.128 [M+Na+HCOOH]+ 79.102 [M+H+HCOOH]+ 101.084"}.The
#'    updated code splits the string on adduct boundaries and checks whether
#'    \emph{any} of CAMERA's candidate neutral mass values matches the MS1Match neutral mass
#'    within \code{ppm_tol}.
#'
#' @param eval_full \code{data.frame}. Output of \code{\link{build_eval_full}}.
#'   Must contain: \code{reference_tier}, \code{CAMERA_adduct}, and at least
#'   one of \code{grp_anno_match} or \code{exp_anno_match}.
#' @param anchors \code{data.frame}. Output of \code{\link{build_anchors}}.
#' @param ppm_tol \code{numeric(1)}. ppm tolerance for neutral mass matching. Default 10.
#' @param mda_tol \code{numeric(1)}. Absolute mass tolerance floor in millidaltons. A match
#'   is accepted if either the ppm difference is within \code{ppm_tol} OR the absolute
#'   difference is within \code{mda_tol} mDa. Default 5 mDa, which covers the worst-case
#'   truncation error from CAMERA reporting neutral masses to 3 decimal places.
#'
#' @return A named \code{list}:
#' \describe{
#'   \item{\code{per_compound}}{Per \code{Lib_group} annotation accuracy counts
#'     (Tier1a/Tier1a_exp only — the established methodology), with columns
#'     \code{n_total}, \code{n_unannotated}, \code{n_annotated},
#'     \code{n_correct}, \code{n_wrong}, \code{coverage_pct},
#'     \code{accuracy_pct}.}
#'   \item{\code{summary}}{Two-row summary table (one per tier) plus a TOTAL
#'     row, with the same columns. Unaffected by the all-tier extension to
#'     \code{feature_level} below.}
#'   \item{\code{feature_level}}{Row-level table covering EVERY tier (not
#'     only Tier1a/Tier1a_exp), one row per feature x compound instance,
#'     with \code{feature_name}, \code{Lib_group}, \code{Compound_name},
#'     \code{tier_simple}, \code{has_CAMERA_annotation}, and
#'     \code{annotation_match}. Every feature belonging to a compound shares
#'     that compound's one true neutral mass (established from its own
#'     Tier1a/Tier1a_exp evidence), so CAMERA's adduct call on ANY tier's
#'     feature can be checked against it. \code{annotation_match} is
#'     \code{NA} when it cannot be evaluated (no CAMERA annotation present,
#'     or the compound has no resolvable neutral mass at all), and
#'     TRUE/FALSE only when CAMERA's call was actually checked — joinable
#'     into \code{\link{build_compound_instance_table}}.}
#' }
#'
#' @importFrom dplyr filter mutate group_by summarise bind_rows select
#' @export
build_anno_accuracy <- function(eval_full, anchors, ppm_tol = 10, mda_tol = 5) {

  stopifnot(
    is.data.frame(eval_full),
    is.data.frame(anchors),
    is.numeric(ppm_tol), length(ppm_tol) == 1L,
    is.numeric(mda_tol), length(mda_tol) == 1L
  )

  # ---------------------------------------------------------------------------
  # Helper: extract MS1Match neutral mass from a single string "[M+H]+ 195.065"
  # ---------------------------------------------------------------------------
  .extract_mass_single <- function(x) {
    if (is.na(x) || nchar(trimws(x)) == 0L) return(NA_real_)
    m <- regmatches(x, regexpr("[0-9]+\\.[0-9]+$", x))
    if (length(m) == 0L || m == "") return(NA_real_)
    as.numeric(m)
  }

  # ---------------------------------------------------------------------------
  # Helper: check whether any of CAMERA's with multiple adduct neutral mass
  # values matches the MS1Match neutral mass within ppm_tol.

  # ---------------------------------------------------------------------------
  .any_annotation_match <- function(CAMERA_str, MS1Match_mz) {

    if (is.na(MS1Match_mz))                              return(NA)
    if (is.na(CAMERA_str) || nchar(trimws(CAMERA_str)) == 0L) return(NA)

    CAMERA_str_anno <- trimws(unlist(strsplit(CAMERA_str, "(?<=\\S)\\s+(?=\\[)",
                                     perl = TRUE)))
    cam_mass_vals <- vapply(CAMERA_str_anno, function(tok) {
      m <- regmatches(tok, regexpr("[0-9]+\\.[0-9]+$", tok))
      if (length(m) == 0L || m == "") NA_real_ else as.numeric(m)
    }, numeric(1L))

    cam_mass_vals <- cam_mass_vals[!is.na(cam_mass_vals)]
    if (length(cam_mass_vals) == 0L) return(NA)

    # check ppm difference OR absolute difference

    abs_diffs <- abs(cam_mass_vals - MS1Match_mz)
    any(abs_diffs / MS1Match_mz * 1e6 <= ppm_tol | abs_diffs <= mda_tol * 1e-3)
  }

  # ---------------------------------------------------------------------------
  #  required columns
  # ---------------------------------------------------------------------------
  has_grp <- "grp_anno_match" %in% names(eval_full)
  has_exp <- "exp_anno_match" %in% names(eval_full)

  if (!has_grp && !has_exp)
    stop("build_anno_accuracy: eval_full must contain 'grp_anno_match' ",
         "and/or 'exp_anno_match'.")

  if (!"CAMERA_adduct" %in% names(eval_full))
    stop("build_anno_accuracy: eval_full must contain 'CAMERA_adduct'.")

  # ---------------------------------------------------------------------------
  # Tier1a and Tier1a_exp
  # ---------------------------------------------------------------------------
  tier1a_data <- eval_full %>%
    dplyr::filter(grepl("1a", reference_tier, ignore.case = TRUE)) %>%
    dplyr::mutate(
      tier_simple = dplyr::if_else(
        grepl("exp", reference_tier, ignore.case = TRUE),
        "Tier1a_exp", "Tier1a"
      ),
      MS1Match_adduct = dplyr::if_else(
        tier_simple == "Tier1a_exp",
        if (has_exp) exp_anno_match else NA_character_,
        if (has_grp) grp_anno_match else NA_character_
      ),
      has_CAMERA_annotation = !is.na(CAMERA_adduct) &
                           nchar(trimws(CAMERA_adduct)) > 0L,
      MS1Match_mz      = vapply(MS1Match_adduct,
                                 .extract_mass_single, numeric(1L)),
      # Multi-annotation
      annotation_match      = mapply(.any_annotation_match, CAMERA_adduct, MS1Match_mz)
    )

  # ---------------------------------------------------------------------------
  # Per-compound summary
  # ---------------------------------------------------------------------------
  per_compound <- tier1a_data %>%
    dplyr::group_by(Lib_group, Compound_name, tier_simple) %>%
    dplyr::summarise(
      n_total      = dplyr::n(),
      n_unannotated = sum(!has_CAMERA_annotation),
      n_annotated   = sum(has_CAMERA_annotation),
      n_correct    = sum(annotation_match,                    na.rm = TRUE),
      n_wrong      = sum(has_CAMERA_annotation & !annotation_match, na.rm = TRUE),
      coverage_pct = round(100 * n_annotated / n_total,    1),
      accuracy_pct = round(100 * n_correct  / n_annotated, 1),
      .groups = "drop"
    )

  # ---------------------------------------------------------------------------
  # Summary table: one row per tier + TOTAL
  # ---------------------------------------------------------------------------
  .summarise_tier <- function(df, tier_annotation) {
    df %>%
      dplyr::summarise(
        tier         = tier_annotation,
        n_total      = dplyr::n(),
        n_unannotated = sum(!has_CAMERA_annotation),
        n_annotated   = sum(has_CAMERA_annotation),
        n_correct    = sum(annotation_match,                    na.rm = TRUE),
        n_wrong      = sum(has_CAMERA_annotation & !annotation_match, na.rm = TRUE),
        coverage_pct = round(100 * n_annotated / n_total,    1),
        accuracy_pct = round(100 * n_correct  / n_annotated, 1)
      )
  }

  summary_by_tier <- dplyr::bind_rows(
    tier1a_data %>%
      dplyr::filter(tier_simple == "Tier1a") %>%
      .summarise_tier("Tier1a"),
    tier1a_data %>%
      dplyr::filter(tier_simple == "Tier1a_exp") %>%
      .summarise_tier("Tier1a_exp"),
    tier1a_data %>%
      .summarise_tier("TOTAL")
  )

  # ---------------------------------------------------------------------------
  # Row-level correctness (TRUE/FALSE/NA per feature x compound instance)
  # ---------------------------------------------------------------------------
  compound_mz_lookup <- tier1a_data %>%
    dplyr::filter(!is.na(MS1Match_mz)) %>%
    dplyr::mutate(.priority = dplyr::if_else(tier_simple == "Tier1a", 1L, 2L)) %>%
    dplyr::arrange(.priority) %>%
    dplyr::distinct(Compound_name, .keep_all = TRUE) %>%
    dplyr::select(Lib_group, Compound_name, compound_mz = MS1Match_mz)

  lib_mz_lookup <- compound_mz_lookup %>%
    dplyr::filter(!is.na(Lib_group)) %>%
    dplyr::distinct(Lib_group, .keep_all = TRUE) %>%
    dplyr::select(Lib_group, compound_mz)

  cpd_mz_lookup <- compound_mz_lookup %>%
    dplyr::select(Compound_name, compound_mz_cpd = compound_mz)

  feature_level <- eval_full %>%
    dplyr::mutate(.row_id = dplyr::row_number()) %>%
    dplyr::left_join(lib_mz_lookup, by = "Lib_group") %>%
    dplyr::left_join(cpd_mz_lookup, by = "Compound_name") %>%
    dplyr::arrange(.row_id) %>%
    dplyr::mutate(
      compound_mz = dplyr::coalesce(compound_mz, compound_mz_cpd),
      tier_simple = dplyr::case_when(
        grepl("1a.*exp|1a_exp", reference_tier, ignore.case = TRUE) ~ "Tier1a_exp",
        grepl("1a",             reference_tier, ignore.case = TRUE) ~ "Tier1a",
        grepl("1b.*exp|1b_exp", reference_tier, ignore.case = TRUE) ~ "Tier1b_exp",
        grepl("1b",             reference_tier, ignore.case = TRUE) ~ "Tier1b",
        grepl("1c.*exp|1c_exp|MS2Fragment_exp", reference_tier, ignore.case = TRUE) ~ "Tier1c_exp",
        grepl("1c",             reference_tier, ignore.case = TRUE) ~ "Tier1c",
        grepl("1d.*exp|1d_exp", reference_tier, ignore.case = TRUE) ~ "Tier1d_exp",
        grepl("1d",             reference_tier, ignore.case = TRUE) ~ "Tier1d",
        grepl("Tier2|tier2",    reference_tier, ignore.case = TRUE) ~ "Tier2",
        TRUE                                                       ~ "Other"
      ),
      has_CAMERA_annotation = !is.na(CAMERA_adduct) & nchar(trimws(CAMERA_adduct)) > 0L,
      annotation_match       = mapply(.any_annotation_match, CAMERA_adduct, compound_mz)
    ) %>%
    dplyr::select(feature_name, Lib_group, Compound_name, tier_simple,
                  has_CAMERA_annotation, annotation_match)

  list(per_compound = per_compound, summary = summary_by_tier,
      feature_level = feature_level)
}


# =============================================================================
# build_benchmark_table
# =============================================================================

#' Build the tier-stratified benchmark summary table
#'
#' @description
#' Aggregates grouping recall (from \code{build_compound_eval}) and annotation
#' accuracy (from \code{build_anno_accuracy}) into a single table,
#' one row per annotation tier.
#'
#' @param compound_eval_wide \code{data.frame}. Output of
#'   \code{\link{build_compound_eval}}.
#' @param anno_accuracy \code{list}. Output of \code{\link{build_anno_accuracy}}.
#' @param purity \code{list}. Output of \code{\link{build_pcgroup_purity}}.
#'
#' @return A \code{data.frame} with columns: \code{Tier},
#'   \code{n_MS1Match}, \code{n_CAMERA_grouped}, \code{recall_pct},
#'   \code{n_CAMERA_split}, \code{split_pct},
#'   \code{CAMERA_annotation_accuracy_pct}, \code{CAMERA_has_concept}.
#'
#' @importFrom dplyr summarise across starts_with
#' @export
build_benchmark_table <- function(compound_eval_wide, anno_accuracy, purity) {

  stopifnot(
    is.data.frame(compound_eval_wide),
    is.list(anno_accuracy), "summary" %in% names(anno_accuracy),
    is.list(purity)
  )


  tiers <- c("Tier1a", "Tier1a_exp", "Tier1b", "Tier1b_exp",
             "Tier1c", "Tier1c_exp", "Tier1d", "Tier2", "TOTAL")


  CAMERA_concept <- c(
    Tier1a     = "Yes",
    Tier1a_exp = "Yes",
    Tier1b     = "Yes",
    Tier1b_exp = "Yes",
    Tier1c     = "No",
    Tier1c_exp = "No",
    Tier1d     = "No",
    Tier2      = "Yes",
    TOTAL      = "-"
  )

  rows <- lapply(tiers, function(tier) {
    n_col      <- paste0("n_MS1Match_",     tier)
    grp_col    <- paste0("n_CAMERA_grouped_", tier)
    spl_col    <- paste0("n_CAMERA_split_",   tier)



    n_found    <- if (n_col   %in% names(compound_eval_wide))
      sum(compound_eval_wide[[n_col]],   na.rm = TRUE) else NA_integer_
    n_grouped  <- if (grp_col %in% names(compound_eval_wide))
      sum(compound_eval_wide[[grp_col]], na.rm = TRUE) else NA_integer_
    n_split    <- if (spl_col %in% names(compound_eval_wide))
      sum(compound_eval_wide[[spl_col]], na.rm = TRUE) else NA_integer_


    una_col      <- paste0("n_CAMERA_unassigned_", tier)
    n_unassigned <- if (una_col %in% names(compound_eval_wide))
      sum(compound_eval_wide[[una_col]], na.rm = TRUE) else NA_integer_
    unassigned_pct <- if (!is.na(n_found) && n_found > 0)
      round(100 * n_unassigned / n_found, 1) else NA_real_

    recall_pct <- if (!is.na(n_found) && n_found > 0)
      round(100 * n_grouped / n_found, 1) else NA_real_
    split_pct  <- if (!is.na(n_found) && n_found > 0)
      round(100 * n_split   / n_found, 1) else NA_real_

    # Coverage and annotation accuracy — only for Tier1a / Tier1a_exp / TOTAL.

    .get_acc <- function(col) {
      acc <- anno_accuracy$summary
      if (is.null(acc) || nrow(acc) == 0L) return(NA_real_)
      tier_col <- intersect(c("tier", "tier_simple"), names(acc))[1L]
      if (!is.na(tier_col)) {
        idx <- which(acc[[tier_col]] == tier)
        if (length(idx) == 0L) return(NA_real_)
      } else {
        idx <- 1L
      }
      if (!col %in% names(acc)) return(NA_real_)
      acc[[col]][idx[1L]]
    }

    coverage_pct <- if (tier %in% c("Tier1a", "Tier1a_exp", "TOTAL"))
      .get_acc("coverage_pct") else NA_real_

    annotation_acc    <- if (tier %in% c("Tier1a", "Tier1a_exp", "TOTAL"))
      .get_acc("accuracy_pct") else NA_real_

    data.frame(
      Tier                      = tier,
      n_MS1Match                = n_found,
      n_CAMERA_grouped          = n_grouped,
      recall_pct                = recall_pct,
      n_CAMERA_split            = n_split,
      split_pct                 = split_pct,
      CAMERA_coverage_pct       = coverage_pct,
      CAMERA_annotation_accuracy_pct = annotation_acc,
      CAMERA_has_concept        = CAMERA_concept[tier],
      n_CAMERA_unassigned = n_unassigned,
      unassigned_pct      = unassigned_pct,
      stringsAsFactors          = FALSE
    )
  })
  out <- do.call(rbind, rows)

  tier_rows <- out[out$Tier != "TOTAL", ]

  total_n       <- sum(tier_rows$n_MS1Match,       na.rm = TRUE)
  total_grouped <- sum(tier_rows$n_CAMERA_grouped, na.rm = TRUE)
  total_split   <- sum(tier_rows$n_CAMERA_split,   na.rm = TRUE)
  total_unassigned <- sum(tier_rows$n_CAMERA_unassigned, na.rm = TRUE)
  out[out$Tier == "TOTAL", "n_CAMERA_unassigned"] <- total_unassigned
  out[out$Tier == "TOTAL", "unassigned_pct"]      <- round(100 * total_unassigned / total_n, 1)
  out[out$Tier == "TOTAL", "n_MS1Match"]       <- total_n        # 1511
  out[out$Tier == "TOTAL", "n_CAMERA_grouped"] <- total_grouped  # 1230
  out[out$Tier == "TOTAL", "n_CAMERA_split"]   <- total_split    # 261
  out[out$Tier == "TOTAL", "recall_pct"]       <- round(100 * total_grouped / total_n, 1)  # 81.4
  out[out$Tier == "TOTAL", "split_pct"]        <- round(100 * total_split   / total_n, 1)  # 17.3

  out


}


# =============================================================================
# build_compound_instance_table
# =============================================================================

#' Build a flat, all-compounds feature table with CAMERA grouping context
#'
#' @description
#' Long table: one row per feature within a compound instance,
#' across EVERY compound instance present in \code{eval_full} (not one
#' queried compound at a time, unlike \code{\link{inspect_compound}}). Reuses
#' the same per-feature anchor-co-grouping logic as \code{inspect_compound}'s
#' internal \code{long_table}, and adds:
#' \describe{
#'   \item{\code{CAMERA_isotopes}, \code{CAMERA_adduct}}{CAMERA's own native
#'     annotation columns, kept SEPARATE (not combined into one string) so
#'     downstream parsing of CAMERA's adduct hypotheses (e.g. neutral-mass
#'     cross-checks) is not obstructed by a merged column.}
#'   \item{\code{CAMERA_anno_correct}}{TRUE/FALSE: does CAMERA's own adduct
#'     call on THIS feature match the compound's true neutral mass (within
#'     \code{anno_accuracy}'s ppm tolerance), established from that
#'     compound's own Tier1a/Tier1a_exp evidence. Applies across ALL tiers
#'     (a Tier1b/1c/1d/Tier2 feature's CAMERA adduct call is checked against
#'     the same compound-level mass), not only Tier1a/Tier1a_exp. \code{NA}
#'     when it cannot be evaluated: CAMERA assigned no adduct annotation at
#'     all, or the compound has no resolvable neutral mass.}
#'   \item{\code{CAMERA_Anchor_pcgroup_purity}}{Number of distinct reference
#'     compounds whose Tier1a anchor falls in this compound's anchor
#'     pcgroup. How many OTHER compounds is CAMERA merging this one with.
#'     Shown only on the row(s) co-grouped with the anchor.}
#'   \item{\code{CAMERA_Anchor_pcgroup_Refconfirmed_pct}}{This compound's own
#'     anchor-co-grouped feature count (X) / total anchor pcgroup feature
#'     count (Y) x 100. Shown only on the row(s) co-grouped with the anchor.}
#' }
#'
#' @param eval_full \code{data.frame}. Output of \code{\link{build_eval_full}}.
#' @param anchors \code{data.frame}. Output of \code{\link{build_anchors}}.
#' @param purity \code{list}. Output of \code{\link{build_pcgroup_purity}}.
#'   The \code{all_pcgroups} element is used for purity / feature-count context.
#' @param CAMERA_annotated_pool_FT \code{data.frame}. CAMERA-annotated feature
#'   table. Used ONLY as a fallback source for pcgroup/adduct/isotopes if
#'   \code{eval_full} is missing those native columns -- normally
#'   \code{eval_full}'s own \code{pcgroup}, \code{CAMERA_adduct}, and
#'   \code{isotopes} columns are used directly (the same columns
#'   \code{build_anno_accuracy()} reads), so this table and
#'   \code{CAMERA_anno_correct} can never silently disagree for the same row.
#' @param anno_accuracy \code{list}. Output of \code{\link{build_anno_accuracy}}.
#'   Its \code{feature_level} element is joined in to populate
#'   \code{CAMERA_anno_correct}.
#'
#' @return A \code{data.frame}, one row per feature within a compound
#'   instance, with columns: \code{compound_instance}, \code{MS1Match_lib_group},
#'   \code{MS1Match_reference_tier}, \code{feature_name},
#'   \code{MS1Match_anno}, \code{MS1Match_iso}, \code{CAMERA_group_status},
#'   \code{CAMERA_isotopes}, \code{CAMERA_adduct}, \code{CAMERA_anno_correct},
#'   \code{CAMERA_pcgroup}, \code{CAMERA_Anchor_pcgroup},
#'   \code{CAMERA_Anchor_pcgroup_purity}, \code{CAMERA_Anchor_pcgroup_Refconfirmed_pct}.
#'
#' @seealso \code{\link{inspect_compound}}, \code{\link{build_pcgroup_purity}},
#'   \code{\link{build_anno_accuracy}}
#'
#' @importFrom dplyr mutate left_join select case_when group_by ungroup
#'   summarise n filter arrange coalesce any_of distinct rename
#' @export
build_compound_instance_table <- function(eval_full, anchors, purity,
                                           CAMERA_annotated_pool_FT,
                                           anno_accuracy) {

  stopifnot(
    is.data.frame(eval_full),
    is.data.frame(anchors),
    is.list(purity), "all_pcgroups" %in% names(purity),
    is.data.frame(CAMERA_annotated_pool_FT),
    is.list(anno_accuracy), "feature_level" %in% names(anno_accuracy)
  )

  tier_levels <- c("Tier1a", "Tier1a_exp", "Tier1b", "Tier1b_exp",
                   "Tier1c", "Tier1c_exp", "Tier1d", "Tier1d_exp", "Tier2")

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

  .combine_str <- function(...) {
    cols <- list(...)
    n    <- length(cols[[1]])
    out  <- character(n)
    for (i in seq_len(n)) {
      parts <- vapply(cols, function(x) x[i], character(1))
      parts <- parts[!is.na(parts) & nchar(trimws(parts)) > 0L]
      out[i] <- if (length(parts) == 0L) NA_character_ else paste(parts, collapse = "; ")
    }
    out
  }

  ef <- eval_full %>%
    dplyr::mutate(tier_simple = .tier_simple(reference_tier))

  ef <- ef %>%
    dplyr::mutate(
      .ms2_anno = dplyr::if_else(
        grepl("Tier1c", reference_tier, fixed = FALSE) &
          !is.na(ms2_parent_feature) & ms2_parent_feature != "",
        paste0("MS2 fragment (", ms2_parent_feature, ")"),
        NA_character_
      )
    )

  ef <- ef %>%
    dplyr::mutate(
      MS1Match_anno = .combine_str(
        if ("grp_anno_match" %in% names(ef)) ef[["grp_anno_match"]] else
          if ("adduct_match"    %in% names(ef)) ef[["adduct_match"]]     else rep(NA_character_, nrow(ef)),
        if ("exp_anno_match" %in% names(ef)) ef[["exp_anno_match"]] else rep(NA_character_, nrow(ef)),
        if ("NL_match"         %in% names(ef)) ef[["NL_match"]]         else rep(NA_character_, nrow(ef)),
        ef[[".ms2_anno"]]
      ),
      MS1Match_iso = .combine_str(
        if ("iso_annotation"   %in% names(ef)) ef[["iso_annotation"]]   else rep(NA_character_, nrow(ef))
      )
    ) %>%
    dplyr::select(-.ms2_anno)


  cpd_to_libgroup <- anchors %>%
    dplyr::select(Compound_name, Lib_group) %>%
    dplyr::distinct(Compound_name, .keep_all = TRUE) %>%
    dplyr::rename(.lib_group_lookup = Lib_group)

  anchor_by_lib <- anchors %>%
    dplyr::select(Lib_group, anchor_pcgroup)

  anchor_by_cpd <- anchors %>%
    dplyr::select(Compound_name, anchor_pcgroup_cpd = anchor_pcgroup) %>%
    dplyr::distinct(Compound_name, .keep_all = TRUE)


  has_native_pcgroup <- "pcgroup"       %in% names(eval_full)
  has_native_adduct  <- "CAMERA_adduct" %in% names(eval_full)
  has_native_iso     <- "isotopes"      %in% names(eval_full)

  if (!has_native_pcgroup || !has_native_adduct || !has_native_iso) {
    cam_adduct_col <- if ("CAMERA_adduct" %in% names(CAMERA_annotated_pool_FT)) "CAMERA_adduct"
                       else if ("adduct" %in% names(CAMERA_annotated_pool_FT)) "adduct"
                       else NULL
    cam_iso_col <- if ("isotopes" %in% names(CAMERA_annotated_pool_FT)) "isotopes" else NULL

    cam_cols <- c("feature_name", "pcgroup", cam_adduct_col, cam_iso_col)

    cam_sub <- CAMERA_annotated_pool_FT %>%
      dplyr::select(dplyr::any_of(cam_cols)) %>%
      dplyr::rename(.fallback_pcgroup = pcgroup)
    if (!is.null(cam_adduct_col))
      cam_sub <- dplyr::rename(cam_sub, .fallback_adduct = dplyr::any_of(cam_adduct_col))
    if (!is.null(cam_iso_col))
      cam_sub <- dplyr::rename(cam_sub, .fallback_iso = dplyr::any_of(cam_iso_col))
  }

  long <- ef %>%
    dplyr::mutate(
      .join_lib = dplyr::coalesce(
        .data$Lib_group,
        if ("ms2_lib_group" %in% names(.)) .data$ms2_lib_group else NA_character_
      )
    ) %>%
    dplyr::left_join(anchor_by_lib, by = c(".join_lib" = "Lib_group")) %>%
    dplyr::select(-.join_lib) %>%
    dplyr::left_join(anchor_by_cpd, by = "Compound_name") %>%
    dplyr::mutate(
      anchor_pcgroup = dplyr::coalesce(.data$anchor_pcgroup, .data$anchor_pcgroup_cpd)
    ) %>%
    dplyr::select(-anchor_pcgroup_cpd)

  if (!has_native_pcgroup || !has_native_adduct || !has_native_iso) {
    long <- long %>% dplyr::left_join(cam_sub, by = "feature_name")
  }

  long <- long %>%
    dplyr::left_join(cpd_to_libgroup, by = "Compound_name") %>%
    dplyr::mutate(
      compound_instance    = Compound_name,

      MS1Match_lib_group   = dplyr::case_when(
        !is.na(.data$Lib_group) ~ .data$Lib_group,
        !is.na(.data$.lib_group_lookup) ~ paste0(.data$.lib_group_lookup, "_exp"),
        TRUE ~ NA_character_
      ),
      CAMERA_pcgroup_resolved = if (has_native_pcgroup) .data$pcgroup else .data$.fallback_pcgroup,
      CAMERA_group_status  = dplyr::case_when(
        is.na(CAMERA_pcgroup_resolved)            ~ "unassigned",
        is.na(anchor_pcgroup)                     ~ NA_character_,
        CAMERA_pcgroup_resolved == anchor_pcgroup ~ "grouped",
        TRUE                                      ~ "SPLIT"
      ),

      CAMERA_isotopes = if (has_native_iso)    .data$isotopes      else .data$.fallback_iso,
      CAMERA_adduct    = if (has_native_adduct) .data$CAMERA_adduct else .data$.fallback_adduct
    ) %>%
    dplyr::select(-.lib_group_lookup)

  # --- CAMERA annotation correctness (TRUE/FALSE), Tier1a/Tier1a_exp only --
  acc_lookup <- anno_accuracy$feature_level %>%
    dplyr::select(feature_name, Compound_name, CAMERA_anno_correct = annotation_match) %>%
    dplyr::distinct(feature_name, Compound_name, .keep_all = TRUE)

  long <- long %>%
    dplyr::left_join(
      acc_lookup,
      by = c("feature_name" = "feature_name", "compound_instance" = "Compound_name")
    )

  # --- Compound-instance-level metrics: purity & Ref-confirmed pct ---------

  own_grouped <- long %>%
    dplyr::filter(CAMERA_group_status == "grouped") %>%
    dplyr::mutate(.base_lib = sub("_exp$", "", .data$MS1Match_lib_group)) %>%
    dplyr::group_by(.base_lib) %>%
    dplyr::summarise(n_own_grouped = dplyr::n(), .groups = "drop")

  pcgroup_stats <- purity$all_pcgroups %>%
    dplyr::select(pcgroup, pcgroup_n_features, pcgroup_n_compounds)

  metric_lookup <- anchors %>%
    dplyr::select(Lib_group, anchor_pcgroup) %>%
    dplyr::left_join(pcgroup_stats, by = c("anchor_pcgroup" = "pcgroup")) %>%
    dplyr::left_join(own_grouped, by = c("Lib_group" = ".base_lib")) %>%
    dplyr::mutate(
      n_own_grouped = dplyr::coalesce(n_own_grouped, 0L),
      CAMERA_Anchor_pcgroup_purity        = pcgroup_n_compounds,
      CAMERA_Anchor_pcgroup_Refconfirmed_pct = dplyr::if_else(
        !is.na(pcgroup_n_features) & pcgroup_n_features > 0L,
        round(100 * n_own_grouped / pcgroup_n_features, 1),
        NA_real_
      )
    ) %>%
    dplyr::select(Lib_group, CAMERA_Anchor_pcgroup_purity,
                  CAMERA_Anchor_pcgroup_Refconfirmed_pct)

  long <- long %>%
    dplyr::mutate(.join_lib_group = sub("_exp$", "", .data$MS1Match_lib_group)) %>%
    dplyr::left_join(metric_lookup, by = c(".join_lib_group" = "Lib_group")) %>%
    dplyr::select(-.join_lib_group) %>%
    dplyr::mutate(
      CAMERA_Anchor_pcgroup = dplyr::if_else(
        CAMERA_group_status == "grouped", as.numeric(anchor_pcgroup), NA_real_
      ),
      CAMERA_Anchor_pcgroup_purity = dplyr::if_else(
        CAMERA_group_status == "grouped", as.numeric(CAMERA_Anchor_pcgroup_purity), NA_real_
      ),
      CAMERA_Anchor_pcgroup_Refconfirmed_pct = dplyr::if_else(
        CAMERA_group_status == "grouped", as.numeric(CAMERA_Anchor_pcgroup_Refconfirmed_pct), NA_real_
      ),
      tier_simple = factor(tier_simple, levels = tier_levels, ordered = TRUE)
    ) %>%
    dplyr::arrange(compound_instance, tier_simple, mz) %>%
    dplyr::select(
      compound_instance,
      MS1Match_lib_group,
      MS1Match_reference_tier         = tier_simple,
      feature_name,
      mz,
      rt,
      MS1Match_anno,
      MS1Match_iso,
      CAMERA_group_status,
      CAMERA_isotopes,
      CAMERA_adduct,
      CAMERA_anno_correct,
      CAMERA_pcgroup                  = CAMERA_pcgroup_resolved,
      CAMERA_Anchor_pcgroup,
      CAMERA_Anchor_pcgroup_purity,
      CAMERA_Anchor_pcgroup_Refconfirmed_pct
    )

  long
}


# =============================================================================
# build_full_feature_table
# =============================================================================

#' Build the complete CAMERA feature table with reference-relative status
#'
#' @description
#' Extends \code{\link{build_compound_instance_table}} to cover every CAMERA
#' feature in \code{CAMERA_annotated_pool_FT}, not only reference-bearing
#' ones. Reference-bearing rows are produced exactly as
#' \code{build_compound_instance_table} would; every CAMERA feature with no
#' reference entry at all is added as an additional row. Every row — both
#' kinds — carries a single \code{ref_status} column using the our definitions
#' , computed identically to (and consistent with)
#' \code{\link{classify_ref_status}}:
#' \describe{
#'   \item{\strong{Ref-confirmed}}{Co-grouped with the compound's anchor.}
#'   \item{\strong{Ref-conflicting}}{Reference entry exists, split from anchor.}
#'   \item{\strong{Ref-unmatched, co-grouped}}{No reference entry, but sits
#'     in a pcgroup that contains some OTHER compound's anchor — diagnostic
#'     of over-grouping.}
#'   \item{\strong{Ref-unmatched, unrelated}}{No reference entry, and its
#'     pcgroup contains no anchor at all — a reference-set coverage gap, not
#'     a CAMERA grouping outcome.}
#'   \item{\strong{Ref-unassigned}}{Reference entry exists, but the tool
#'     assigned the feature no group (zero for CAMERA).}
#' }
#'
#' Every row has \code{mz} and \code{rt}: reference rows from MS1Match, the
#' other rows from the CAMERA peak list (\code{mz}, \code{rt}).
#'
#' For Ref-unmatched rows, \code{compound_instance}, \code{MS1Match_lib_group},
#' \code{MS1Match_reference_tier}, \code{MS1Match_anno}, and \code{MS1Match_iso} are \code{NA}
#' (no compound owns the row). \code{CAMERA_Anchor_pcgroup_Refconfirmed_pct}
#' and \code{CAMERA_anno_correct} are also left \code{NA} for these rows by
#' design: the former is a specific compound's own share of its anchor
#' pcgroup, the latter requires an MS1Match-expected mass to check CAMERA's
#' call against — neither exists when no compound owns the row.
#'
#' @param eval_full \code{data.frame}. Output of \code{\link{build_eval_full}}.
#' @param anchors \code{data.frame}. Output of \code{\link{build_anchors}}.
#' @param purity \code{list}. Output of \code{\link{build_pcgroup_purity}}.
#' @param CAMERA_annotated_pool_FT \code{data.frame}. The FULL CAMERA-annotated
#'   feature table (every feature, not pre-filtered to reference matches).
#' @param anno_accuracy \code{list}. Output of \code{\link{build_anno_accuracy}}.
#'   Passed through to \code{\link{build_compound_instance_table}} to
#'   populate \code{CAMERA_anno_correct}.
#'
#' @return A \code{data.frame}: all \code{build_compound_instance_table}
#'   columns (including \code{CAMERA_isotopes}, \code{CAMERA_adduct} as
#'   separate columns and \code{CAMERA_anno_correct}), plus \code{ref_status}
#'   (ordered factor, 5 levels), covering every feature in
#'   \code{CAMERA_annotated_pool_FT}.
#'
#' @seealso \code{\link{build_compound_instance_table}},
#'   \code{\link{classify_ref_status}}, \code{\link{build_benchmark_table_from_full}}
#'
#' @importFrom dplyr mutate filter select rename left_join bind_rows case_when
#'   if_else arrange desc any_of n
#' @export
build_full_feature_table <- function(eval_full, anchors, purity,
                                      CAMERA_annotated_pool_FT, anno_accuracy) {

  stopifnot(
    is.data.frame(eval_full), is.data.frame(anchors),
    is.list(purity), "all_pcgroups" %in% names(purity),
    is.data.frame(CAMERA_annotated_pool_FT),
    is.list(anno_accuracy), "feature_level" %in% names(anno_accuracy)
  )

  # --- 1. Reference rows: build_compound_instance_table() ----
  ref_rows <- build_compound_instance_table(
    eval_full, anchors, purity, CAMERA_annotated_pool_FT, anno_accuracy
  ) %>%
    dplyr::mutate(
      MS1Match_reference_tier = as.character(MS1Match_reference_tier),
      ref_status = dplyr::case_when(
        CAMERA_group_status == "grouped"   ~ "Ref-confirmed",
        CAMERA_group_status == "SPLIT"     ~ "Ref-conflicting",
        CAMERA_group_status == "unassigned" ~ "Ref-unassigned",
        TRUE                                ~ NA_character_
      )
    )

  # --- 2. Ref-unmatched rows: CAMERA features with NO reference entry -----
  ref_feature_names <- unique(ref_rows$feature_name)

  cam_adduct_col <- if ("CAMERA_adduct" %in% names(CAMERA_annotated_pool_FT)) "CAMERA_adduct"
                     else if ("adduct" %in% names(CAMERA_annotated_pool_FT)) "adduct"
                     else NULL
  cam_iso_col <- if ("isotopes" %in% names(CAMERA_annotated_pool_FT)) "isotopes" else NULL

  anchor_pcgroup_set <- unique(stats::na.omit(anchors$anchor_pcgroup))

  pcgroup_stats <- purity$all_pcgroups %>%
    dplyr::select(pcgroup, pcgroup_n_compounds)

  unmatched <- CAMERA_annotated_pool_FT %>%
    dplyr::filter(!feature_name %in% ref_feature_names) %>%
    dplyr::select(dplyr::any_of(c("feature_name", "mz", "rt", "pcgroup",
                                  cam_adduct_col, cam_iso_col))) %>%
    dplyr::rename(CAMERA_pcgroup = pcgroup)

  if (!is.null(cam_adduct_col))
    unmatched <- dplyr::rename(unmatched, CAMERA_adduct = dplyr::any_of(cam_adduct_col))
  if (!is.null(cam_iso_col))
    unmatched <- dplyr::rename(unmatched, CAMERA_isotopes = dplyr::any_of(cam_iso_col))
  if (!"CAMERA_adduct" %in% names(unmatched))   unmatched$CAMERA_adduct   <- NA_character_
  if (!"CAMERA_isotopes" %in% names(unmatched)) unmatched$CAMERA_isotopes <- NA_character_
  if (!"mz" %in% names(unmatched)) unmatched$mz <- NA_real_
  if (!"rt" %in% names(unmatched)) unmatched$rt <- NA_real_

  unmatched <- unmatched %>%
    dplyr::mutate(
      .is_co_grouped = CAMERA_pcgroup %in% anchor_pcgroup_set,
      ref_status = dplyr::if_else(
        .is_co_grouped, "Ref-unmatched, co-grouped", "Ref-unmatched, unrelated"
      ),
      compound_instance        = NA_character_,
      MS1Match_lib_group       = NA_character_,
      MS1Match_reference_tier  = NA_character_,
      MS1Match_anno            = NA_character_,
      MS1Match_iso             = NA_character_,
      CAMERA_group_status      = "unmatched",

      CAMERA_anno_correct      = NA,
      CAMERA_Anchor_pcgroup    = dplyr::if_else(.is_co_grouped, as.numeric(CAMERA_pcgroup), NA_real_)
    ) %>%
    dplyr::left_join(pcgroup_stats, by = c("CAMERA_Anchor_pcgroup" = "pcgroup")) %>%
    dplyr::mutate(
      CAMERA_Anchor_pcgroup_purity           = as.numeric(pcgroup_n_compounds),

      CAMERA_Anchor_pcgroup_Refconfirmed_pct = NA_real_
    ) %>%
    dplyr::select(
      compound_instance, MS1Match_lib_group, MS1Match_reference_tier,
      feature_name, mz, rt, MS1Match_anno, MS1Match_iso, CAMERA_group_status,
      CAMERA_isotopes, CAMERA_adduct, CAMERA_anno_correct,
      CAMERA_pcgroup, CAMERA_Anchor_pcgroup, CAMERA_Anchor_pcgroup_purity,
      CAMERA_Anchor_pcgroup_Refconfirmed_pct, ref_status
    )

  # --- 3. Combine -------------------------------------------------------------
  status_levels <- c("Ref-confirmed", "Ref-conflicting",
                     "Ref-unmatched, co-grouped", "Ref-unmatched, unrelated",
                     "Ref-unassigned")

  out <- dplyr::bind_rows(ref_rows, unmatched) %>%
    dplyr::mutate(
      ref_status = factor(ref_status, levels = status_levels, ordered = TRUE)
    ) %>%
    { n_na <- sum(is.na(.$ref_status))
      if (n_na > 0L)
        warning("build_full_feature_table: ", n_na, " reference rows have no ",
                "resolvable anchor and were left with ref_status = NA.")
      . } %>%
    dplyr::arrange(
      dplyr::desc(!is.na(compound_instance)), compound_instance,
      MS1Match_reference_tier
    )

  out
}


# =============================================================================
# build_benchmark_table_from_full
# =============================================================================

#' Recompute the tier-stratified benchmark table from build_full_feature_table()
#'
#' @description
#' Produces the identical column layout as \code{\link{build_benchmark_table}},
#' but every count (\code{n_MS1Match}, \code{n_CAMERA_grouped}/recall,
#' \code{n_CAMERA_split}/split, \code{n_CAMERA_unassigned}) is aggregated
#' directly from \code{ref_status} on \code{full_table} — the same row-level
#' object exported for manual inspection — rather than from the separately
#' computed \code{compound_eval_wide} pathway. Running both
#' \code{build_benchmark_table()} and this function and comparing output is
#' a useful consistency check: any divergence indicates the two pipelines
#' have drifted out of sync somewhere upstream.
#'
#' Annotation coverage/accuracy (\code{CAMERA_coverage_pct},
#' \code{CAMERA_annotation_accuracy_pct}) are NOT re-derived here — they
#' depend on a ppm-tolerance match between CAMERA's own adduct call and the
#' expected mass, a different computation axis from \code{ref_status}
#' entirely — and are still sourced from \code{anno_accuracy$summary}, same
#' as \code{build_benchmark_table()}.
#'
#' @param full_table \code{data.frame}. Output of
#'   \code{\link{build_full_feature_table}}. Must contain \code{ref_status}
#'   and \code{MS1Match_reference_tier}.
#' @param anno_accuracy \code{list}. Output of \code{\link{build_anno_accuracy}}.
#'
#' @return A \code{data.frame}, one row per tier plus TOTAL, with columns:
#'   \code{Tier}, \code{n_MS1Match}, \code{n_CAMERA_grouped}, \code{recall_pct},
#'   \code{n_CAMERA_split}, \code{split_pct}, \code{CAMERA_coverage_pct},
#'   \code{CAMERA_annotation_accuracy_pct}, \code{CAMERA_has_concept},
#'   \code{n_CAMERA_unassigned}, \code{unassigned_pct}.
#'
#' @seealso \code{\link{build_full_feature_table}}, \code{\link{build_benchmark_table}}
#'
#' @importFrom dplyr filter
#' @export
build_benchmark_table_from_full <- function(full_table, anno_accuracy) {

  stopifnot(
    is.data.frame(full_table),
    "ref_status" %in% names(full_table),
    "MS1Match_reference_tier" %in% names(full_table),
    is.list(anno_accuracy), "summary" %in% names(anno_accuracy)
  )

  tiers <- c("Tier1a", "Tier1a_exp", "Tier1b", "Tier1b_exp",
            "Tier1c", "Tier1c_exp", "Tier1d", "Tier2")

  CAMERA_concept <- c(
    Tier1a     = "Yes", Tier1a_exp = "Yes",
    Tier1b     = "Yes", Tier1b_exp = "Yes",
    Tier1c     = "No",  Tier1c_exp = "No",
    Tier1d     = "No",  Tier2      = "Yes",
    TOTAL      = "-"
  )

  .get_acc <- function(tier, col) {
    acc <- anno_accuracy$summary
    if (is.null(acc) || nrow(acc) == 0L) return(NA_real_)
    tier_col <- intersect(c("tier", "tier_simple"), names(acc))[1L]
    if (!is.na(tier_col)) {
      idx <- which(acc[[tier_col]] == tier)
      if (length(idx) == 0L) return(NA_real_)
    } else {
      idx <- 1L
    }
    if (!col %in% names(acc)) return(NA_real_)
    acc[[col]][idx[1L]]
  }

  rows <- lapply(c(tiers, "TOTAL"), function(tier) {

    rows_t <- if (tier == "TOTAL")
      dplyr::filter(full_table, !is.na(MS1Match_reference_tier))
    else
      dplyr::filter(full_table, as.character(MS1Match_reference_tier) == tier)

    n_found     <- nrow(rows_t)
    n_confirmed <- sum(rows_t$ref_status == "Ref-confirmed",   na.rm = TRUE)
    n_conflict  <- sum(rows_t$ref_status == "Ref-conflicting", na.rm = TRUE)
    n_unassigned <- sum(rows_t$ref_status == "Ref-unassigned", na.rm = TRUE)

    confirmed_pct <- if (n_found > 0L) round(100 * n_confirmed / n_found, 1) else NA_real_
    conflict_pct  <- if (n_found > 0L) round(100 * n_conflict  / n_found, 1) else NA_real_
    unassigned_pct <- if (n_found > 0L) round(100 * n_unassigned / n_found, 1) else NA_real_

    coverage_pct <- if (tier %in% c("Tier1a", "Tier1a_exp", "TOTAL"))
      .get_acc(tier, "coverage_pct") else NA_real_
    accuracy_pct <- if (tier %in% c("Tier1a", "Tier1a_exp", "TOTAL"))
      .get_acc(tier, "accuracy_pct") else NA_real_

    data.frame(
      Tier                           = tier,
      n_MS1Match                     = n_found,
      n_CAMERA_grouped               = n_confirmed,
      recall_pct                     = confirmed_pct,
      n_CAMERA_split                 = n_conflict,
      split_pct                      = conflict_pct,
      CAMERA_coverage_pct            = coverage_pct,
      CAMERA_annotation_accuracy_pct = accuracy_pct,
      CAMERA_has_concept             = CAMERA_concept[tier],
      n_CAMERA_unassigned            = n_unassigned,
      unassigned_pct                 = unassigned_pct,
      stringsAsFactors                = FALSE
    )
  })

  do.call(rbind, rows)
}
