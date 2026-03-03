#' Join lists of annotation results to library match output
#'
#' @description
#' Combine direct library match results (MS1/MS2) with downstream annotation results
#' (adduct/NL/isotope), using a match-level key:
#' `feature_name + Lib_group + Compound_name`.
#'
#' This supports cases where the same feature appears in multiple library groups or
#' is assigned to multiple compounds.
#'
#' @param FT Required data.frame. Feature table containing at least `feature_name`
#'   and mz/rt/intensity columns.
#' @param direct_lib_list Required list. Output of direct library matching.
#' @param lib_anno_list Optional list. Output of annotation (e.g., `annotate()`),
#'   either a list of small lists or a data.frame.
#' @param target_match_list Optional list/data.frame. Optional targetlist annotations.
#'
#' Join direct library matches + feature annotations (+ optional target PI matches) onto FT
#'
#' - FT is the base (one row per measured feature).
#' - direct_lib_list is match/peak-level evidence and can contain multiple rows per FT feature.
#' - lib_anno_list is feature-level annotation (PI/adduct/NL/iso) and is joined by FT feature id.
#' - To avoid .x/.y duplicates, we DROP overlapping columns from the right tables before joining.
#'
#' @export
join_annotation_results <- function(FT,
                                    direct_lib_list,
                                    lib_anno_list = NULL,
                                    target_match_list = NULL) {

  stopifnot(is.data.frame(FT))
  stopifnot(is.list(direct_lib_list))

  # -----------------------------
  # helpers
  # -----------------------------
  as_df_safe <- function(x) {
    if (is.data.frame(x)) return(x)
    if (is.list(x) && !is.null(names(x))) return(as.data.frame(x, stringsAsFactors = FALSE))
    stop("Cannot coerce element to data.frame (missing names?).")
  }

  # Flatten target PI list -> df
  target_list_to_df <- function(res_list) {
    if (is.null(res_list)) return(NULL)
    if (is.data.frame(res_list)) return(res_list)
    if (!is.list(res_list)) stop("target_match_list must be a list or data.frame.")

    res_list <- Filter(function(x) {
      is.list(x) && !is.null(x$Feature_name) && length(x$Feature_name) > 0
    }, res_list)

    if (length(res_list) == 0) return(NULL)

    out <- do.call(rbind, lapply(res_list, function(x) {
      data.frame(
        Feature_name = x$Feature_name,
        Feature_mz   = x$Feature_mz,
        Comp_name    = x$Comp_name,
        PI_name      = x$PI_name,
        stringsAsFactors = FALSE
      )
    }))
    rownames(out) <- NULL
    out
  }

  # Drop overlapping columns
  drop_overlap_except_key <- function(y, x_names, key) {
    overlap <- intersect(names(y), x_names)
    overlap <- setdiff(overlap, key)
    if (length(overlap) == 0) return(y)
    dplyr::select(y, -dplyr::any_of(overlap))
  }

  # -----------------------------
  # Check FT column names
  # -----------------------------
  mz_col  <- if ("mzmed" %in% names(FT)) "mzmed" else if ("mz" %in% names(FT)) "mz" else NULL
  rt_col  <- if ("rtmed" %in% names(FT)) "rtmed" else if ("rt" %in% names(FT)) "rt" else NULL
  int_col <- if ("intensity_mean" %in% names(FT)) "intensity_mean" else if ("intensity" %in% names(FT)) "intensity" else NULL
  if (is.null(mz_col) || is.null(rt_col)) stop("FT must contain mz/rt columns (mzmed/mz and rtmed/rt).")

  # -----------------------------
  # direct_lib_list -> direct_df
  # -----------------------------
  direct_df <- tryCatch({
    direct_lib_list %>%
      purrr::map(~ .x[setdiff(names(.x), c("Compound_peak_mz_full", "Compound_peak_int_full"))]) %>%
      purrr::map(as_df_safe) %>%
      dplyr::bind_rows()
  }, error = function(e) {
    stop("Failed to bind direct_lib_list: ", conditionMessage(e))
  })


  req_direct <- c("feature_name", "Matched_features_to_peaks", "Lib_group", "Compound_name")
  miss_direct <- setdiff(req_direct, names(direct_df))
  if (length(miss_direct) > 0) {
    stop("Missing required columns in direct match results: ", paste(miss_direct, collapse = ", "))
  }


  direct_df <- direct_df %>%
    dplyr::rename(PI_feature_name = feature_name,
                  feature_name    = Matched_features_to_peaks)

  # -----------------------------
  # lib_anno_list -> anno_df
  # -----------------------------
  anno_df <- NULL
  if (!is.null(lib_anno_list)) {
    anno_df <- tryCatch({
      if (is.data.frame(lib_anno_list)) lib_anno_list else
        lib_anno_list %>% purrr::map(as_df_safe) %>% dplyr::bind_rows()
    }, error = function(e) {
      warning("Failed to bind lib_anno_list: ", conditionMessage(e))
      NULL
    })
  }

  if (!is.null(anno_df)) {
    if (!"feature_name" %in% names(anno_df)) {
      warning("lib_anno_list missing `feature_name`; skipping lib_anno join.")
      anno_df <- NULL
    } else {

      anno_df <- anno_df %>%
        dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(.x, ""))) %>%
        dplyr::group_by(feature_name) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

      anno_df <- dplyr::select(anno_df, -dplyr::any_of(c("mz", "rt", "mzmed", "rtmed", "intensity", "intensity_mean")))

      anno_df <- drop_overlap_except_key(anno_df, names(direct_df), key = "feature_name")
    }
  }

  # -----------------------------
  # target_match_list -> target_df
  # -----------------------------
  target_wide <- NULL
  if (!is.null(target_match_list)) {
    target_df <- tryCatch(target_list_to_df(target_match_list),
                          error = function(e) {
                            warning("Failed to convert target_match_list: ", conditionMessage(e))
                            NULL
                          })

    if (!is.null(target_df)) {
      # normalize join key
      target_df <- target_df %>% dplyr::rename(feature_name = Feature_name)


      if (exists("long_to_wide.fun", mode = "function") &&
          all(c("feature_name", "Comp_name", "PI_name") %in% names(target_df))) {

        target_wide <- long_to_wide.fun(
          target_df,
          col_name1 = feature_name,
          col_name2 = Comp_name,
          col_name3 = PI_name
        )

        # drop overlaps to avoid .x/.y
        target_wide <- drop_overlap_except_key(target_wide, names(direct_df), key = "feature_name")
      } else {
        warning("target_match_list provided but long_to_wide.fun() not found or required cols missing; skipping target join.")
        target_wide <- NULL
      }
    }
  }

  # -----------------------------
  # match_df
  # -----------------------------
  match_df <- direct_df

  if (!is.null(anno_df)) {
    match_df <- match_df %>% dplyr::left_join(anno_df, by = "feature_name")
  }
  if (!is.null(target_wide)) {
    match_df <- match_df %>% dplyr::left_join(target_wide, by = "feature_name")
  }

  # -----------------------------
  # Join match_df to FT

  # -----------------------------
  FT_full <- FT %>%
    dplyr::select(feature_name, dplyr::all_of(c(mz_col, rt_col)), dplyr::any_of(int_col)) %>%
    dplyr::left_join(match_df, by = "feature_name") %>%
    dplyr::arrange(.data[[rt_col]])

  # -----------------------------
  # Reduced tables
  # -----------------------------
  base_cols <- c("feature_name", mz_col, rt_col, int_col)

  optional_cols <- c(
    "Compound_name", "Lib_group", "Ion_type",
    "Compound_peak_mz", "Compound_peak_int",
    "inner_cosine", "outer_cosine", "intensity_ratio",
    "PI_match", "adduct_match", "NL_match", "iso_annotation"
  )

  reduced <- FT_full %>%
    dplyr::select(dplyr::any_of(c(base_cols, optional_cols))) %>%
    data.frame()

  reduced_anno <- reduced %>%
    { if ("Lib_group" %in% names(.)) dplyr::filter(., !is.na(.data$Lib_group)) else . } %>%
    data.frame()

  feature_compound_summary <- FT_full %>%
    dplyr::filter(!is.na(.data$Compound_name), !is.na(.data$Lib_group)) %>%
    dplyr::group_by(feature_name, Compound_name, Lib_group) %>%
    dplyr::summarise(
      n_evidence_rows = dplyr::n(),
      mz = dplyr::first(.data[[mz_col]]),
      rt = dplyr::first(.data[[rt_col]]),
      intensity = if (!is.null(int_col)) dplyr::first(.data[[int_col]]) else NA_real_,
      best_inner_cosine = if ("inner_cosine" %in% names(dplyr::pick(dplyr::everything())))
        suppressWarnings(max(inner_cosine, na.rm = TRUE)) else NA_real_,
      best_outer_cosine = if ("outer_cosine" %in% names(dplyr::pick(dplyr::everything())))
        suppressWarnings(max(outer_cosine, na.rm = TRUE)) else NA_real_,
      PI_match = if ("PI_match" %in% names(dplyr::pick(dplyr::everything()))) dplyr::first(PI_match) else NA_character_,
      adduct_match = if ("adduct_match" %in% names(dplyr::pick(dplyr::everything()))) dplyr::first(adduct_match) else NA_character_,
      NL_match = if ("NL_match" %in% names(dplyr::pick(dplyr::everything()))) dplyr::first(NL_match) else NA_character_,
      iso_annotation = if ("iso_annotation" %in% names(dplyr::pick(dplyr::everything()))) dplyr::first(iso_annotation) else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(rt)

  list(
    export_match_table             = match_df,
    export_feature_compound_summary= feature_compound_summary,
    export_reduced_table           = reduced,
    export_reduced_anno_table      = reduced_anno,
    export_full_table              = FT_full
  )
}

