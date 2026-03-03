#' Generate primary-ion-anchored feature pairs with labels (and annotations)
#'
#' @param annotate_peaks_df Data frame with columns:
#'   feature_name, mz, rt, PI_match, adduct_match, NL_match, iso_annotation
#' @param rt_threshold Retention time window (seconds) to pair candidate features with the PI
#' @param ppm_tolerance PPM tolerance when comparing adduct-calculated M to PI M
#' @return Tibble with pairs and annotations:
#'   feature_1 (PI), feature_2 (candidate), label (1/0),
#'   feature_1_* and feature_2_* annotation columns, and evidence flags.
#' @importFrom dplyr filter select left_join bind_rows
#' @importFrom tibble tibble
#' @importFrom stringr str_extract str_extract_all str_replace_all
#' @importFrom purrr map_int
generate_PI_pairs <- function(annotate_peaks_df,
                              rt_threshold = 10,
                              ppm_tolerance = 10) {
  stopifnot(is.data.frame(annotate_peaks_df))

  # 1) primary-ion anchors
  primary_ions <- dplyr::filter(annotate_peaks_df, PI_match != "")
  message("Found ", nrow(primary_ions), " primary ion anchors.")
  if (nrow(primary_ions) == 0) {
    return(tibble::tibble())
  }

  # quick lookup by feature_name for attaching annotations later
  annot_cols <- c("feature_name", "mz", "rt", "PI_match", "adduct_match", "NL_match", "iso_annotation")
  annot_lookup <- dplyr::select(annotate_peaks_df, dplyr::all_of(annot_cols))

  out_list <- vector("list", nrow(primary_ions))

  # more permissive regex for trailing numeric (handles integers too)
  trailing_num_re <- "(-?\\d+(?:\\.\\d+)?)$"

  for (i in seq_len(nrow(primary_ions))) {
    pi_row <- primary_ions[i, ]
    pi_feat <- pi_row$feature_name
    pi_rt   <- pi_row$rt

    # candidates in RT window, excluding self
    candidates <- dplyr::filter(
      annotate_peaks_df,
      feature_name != pi_feat,
      abs(rt - pi_rt) <= rt_threshold
    )
    if (nrow(candidates) == 0) next

    # parse PI M (from PI_match's trailing number)
    pi_M <- suppressWarnings(as.numeric(stringr::str_extract(pi_row$PI_match, trailing_num_re)))

    # precompute PI isotope group (fallback: from PI_match if missing)
    pi_group <- stringr::str_extract(pi_row$iso_annotation, "(?<=\\[)\\d+(?=\\])")
    if (is.na(pi_group) || pi_group == "") {
      pi_group <- stringr::str_extract(pi_row$PI_match, "(?<=\\[)\\d+(?=\\])")
    }

    # per-candidate labeling + evidence
    label_vec <- purrr::map_int(seq_len(nrow(candidates)), function(j) {
      cand_row <- candidates[j, ]
      is_positive <- FALSE

      # ----- Adduct evidence -----
      adduct_M <- suppressWarnings(as.numeric(stringr::str_extract(cand_row$adduct_match, trailing_num_re)))
      adduct_ok <- FALSE
      if (!is.na(adduct_M) && !is.na(pi_M) && pi_M != 0) {
        ppm_error <- abs(adduct_M - pi_M) / pi_M * 1e6
        adduct_ok <- ppm_error <= ppm_tolerance
        if (adduct_ok) is_positive <- TRUE
      }

      # ----- Isotope evidence -----
      iso_ok <- FALSE
      if (!isTRUE(is.na(cand_row$iso_annotation)) && cand_row$iso_annotation != "") {
        cand_group <- stringr::str_extract(cand_row$iso_annotation, "(?<=\\[)\\d+(?=\\])")
        if (!is.na(cand_group) && !is.na(pi_group) && cand_group == pi_group) {
          iso_ok <- TRUE
          is_positive <- TRUE
        }
      }

      # ----- Neutral loss evidence -----
      nl_ok <- FALSE
      if (!isTRUE(is.na(cand_row$NL_match)) && cand_row$NL_match != "") {
        referenced <- stringr::str_extract_all(cand_row$NL_match, "\\(FT\\d+\\)") |>
          unlist() |>
          stringr::str_replace_all("[()]", "")
        if (length(referenced) && pi_feat %in% referenced) {
          nl_ok <- TRUE
          is_positive <- TRUE
        }
      }

      # stash evidence into columns by side effect? no — we'll recompute below in one sweep
      ifelse(is_positive, 1L, 0L)
    })

    # Build base pair table
    pair_df <- tibble::tibble(
      feature_1 = pi_feat,
      feature_2 = candidates$feature_name,
      label = label_vec
    )

    # Attach annotations for feature_1 (PI) and feature_2 (candidate)
    pair_df <- pair_df |>
      dplyr::left_join(
        dplyr::rename(annot_lookup, feature_1 = feature_name,
                      feature_1_mz = mz, feature_1_rt = rt,
                      feature_1_PI_match = PI_match,
                      feature_1_adduct_match = adduct_match,
                      feature_1_NL_match = NL_match,
                      feature_1_iso_annotation = iso_annotation),
        by = "feature_1"
      ) |>
      dplyr::left_join(
        dplyr::rename(annot_lookup, feature_2 = feature_name,
                      feature_2_mz = mz, feature_2_rt = rt,
                      feature_2_PI_match = PI_match,
                      feature_2_adduct_match = adduct_match,
                      feature_2_NL_match = NL_match,
                      feature_2_iso_annotation = iso_annotation),
        by = "feature_2"
      )

    # Evidence columns (vectorized over the joined table)
    # adduct ppm agreement between feature_2 adduct M and feature_1 PI M
    pi_M_all   <- suppressWarnings(as.numeric(stringr::str_extract(pair_df$feature_1_PI_match, trailing_num_re)))
    adduct_M_2 <- suppressWarnings(as.numeric(stringr::str_extract(pair_df$feature_2_adduct_match, trailing_num_re)))
    ppm_err <- ifelse(!is.na(pi_M_all) & !is.na(adduct_M_2) & pi_M_all != 0,
                      abs(adduct_M_2 - pi_M_all) / pi_M_all * 1e6, NA_real_)
    adduct_ok <- !is.na(ppm_err) & ppm_err <= ppm_tolerance

    # isotope group match
    pi_grp <- stringr::str_extract(pair_df$feature_1_iso_annotation, "(?<=\\[)\\d+(?=\\])")
    pi_grp[is.na(pi_grp) | pi_grp == ""] <- stringr::str_extract(pair_df$feature_1_PI_match, "(?<=\\[)\\d+(?=\\])")
    cand_grp <- stringr::str_extract(pair_df$feature_2_iso_annotation, "(?<=\\[)\\d+(?=\\])")
    iso_ok <- !is.na(pi_grp) & !is.na(cand_grp) & (pi_grp == cand_grp)

    # NL references feature_1 id
    referenced <- stringr::str_extract_all(pair_df$feature_2_NL_match, "\\(FT\\d+\\)")
    nl_ok <- vapply(seq_along(referenced), function(k) {
      if (length(referenced[[k]]) == 0) return(FALSE)
      any(stringr::str_replace_all(referenced[[k]], "[()]", "") == pair_df$feature_1[k])
    }, logical(1))

    pair_df$evidence_adduct_ppm_ok <- adduct_ok
    pair_df$evidence_adduct_ppm    <- ppm_err
    pair_df$evidence_iso_group_ok  <- iso_ok
    pair_df$evidence_nl_ref_ok     <- nl_ok

    out_list[[i]] <- pair_df
  }

  dplyr::bind_rows(out_list)
}










# generate_PI_pairs <- function(annotate_peaks_df, rt_threshold = 10, ppm_tolerance = 10) {
#
#
#   # 1. Identify primary ions (anchors)
#   primary_ions <- annotate_peaks_df %>%
#     filter(PI_match != "")
#
#   message("Found ", nrow(primary_ions), " primary ion anchors.")
#
#   pair_list <- vector("list", length = nrow(primary_ions))
#
#   # 2. Loop through primary ions
#   for (i in seq_len(nrow(primary_ions))) {
#     pi_feat <- primary_ions$feature_name[i]
#     pi_rt   <- primary_ions$rt[i]
#
#     # Candidate features within RT window (excluding PI itself)
#     candidates <- annotate_peaks_df %>%
#       filter(feature_name != pi_feat, abs(rt - pi_rt) <= rt_threshold)
#
#     if (nrow(candidates) == 0) next
#
#     # Label column: default 0
#     label_vec <- purrr::map_int(seq_len(nrow(candidates)), function(j) {
#       is_positive <- FALSE
#
#       # ----- Adduct match -----
#       adduct_M <- suppressWarnings(as.numeric(str_extract(candidates$adduct_match[j], "[0-9]+\\.[0-9]+$")))
#       pi_M <- suppressWarnings(as.numeric(str_extract(primary_ions$PI_match[i], "[0-9]+\\.[0-9]+$")))
#       if (!is.na(adduct_M) && !is.na(pi_M)) {
#         ppm_error <- abs(adduct_M - pi_M) / pi_M * 1e6
#         if (ppm_error <= ppm_tolerance) is_positive <- TRUE
#       }
#
#       # ----- Isotope match -----
#       if (candidates$iso_annotation[j] != "") {
#         candidate_group <- str_extract(candidates$iso_annotation[j], "(?<=\\[)\\d+(?=\\])")
#         pi_group <- str_extract(primary_ions$iso_annotation[i], "(?<=\\[)\\d+(?=\\])")
#         if (is.na(pi_group) || pi_group == "") {
#           pi_group <- str_extract(primary_ions$PI_match[i], "(?<=\\[)\\d+(?=\\])")
#         }
#         if (!is.na(candidate_group) && !is.na(pi_group) && candidate_group == pi_group) {
#           is_positive <- TRUE
#         }
#       }
#
#       # ----- Neutral Loss match -----
#       if (candidates$NL_match[j] != "") {
#         referenced <- str_extract_all(candidates$NL_match[j], "\\(FT\\d+\\)") %>%
#           unlist() %>% str_replace_all("[()]", "")
#         if (pi_feat %in% referenced) is_positive <- TRUE
#       }
#
#       return(ifelse(is_positive, 1L, 0L))
#     })
#
#     # Build minimal pair table (just the 3 columns)
#     pair_df <- tibble(
#       feature_1 = pi_feat,
#       feature_2 = candidates$feature_name,
#       label = label_vec
#     )
#
#     pair_list[[i]] <- pair_df
#   }
#
#   bind_rows(pair_list)
# }
