# ============================================================
# method_annotateISF.R
#
# S4 generic + method for MS2 fragment-based annotation.
#
# Tier naming:
#   Tier1c_MS2Fragment     — MS2 fragment match, Phase 1 (intra-group)
#   Tier1c_MS2Fragment_exp — MS2 fragment match, Phase 2 (RT-group)
#
# Parent eligibility:
#   The ISF parent has to be the DIRECT LIBRARY-MATCH ANCHOR of its Lib_group —
#   the feature assigned Tier1a_AnnoMatch (any ESI adduct accepted).
#   DDA triggers on the most abundant precursor ion, which is not always [M+H]+.
#
#   Candidate scope determines the tier:
#     Tier1c       — candidates restricted to the same Lib_group
#     Tier1c_exp   — candidates are RT-pool features with no Phase 1 assignment
#
# Priority: 1a > 1b > 1c > 1d > 2 > 1a_exp > 1b_exp > 1c_exp > 1d_exp
# ============================================================

# ── Generic ───────────────────────────────────────────────────────────────────
#' @export
setGeneric("annotateISF",
           function(object, ...) standardGeneric("annotateISF"))

# ── S4 method────────────────────────────────────────────────────

#' Annotate features as MS2 fragment ions (Tier 1c / 1c_exp)
#'
#' @description
#' Dispatched on a \code{\link{PeakAnnotationISF-class}} object.
#' link to \code{\link{run_isf_annotation}()}.
#'
#' \strong{Phase 1 — Intra-group MS2 fragment (Tier 1c)}
#' For each Lib_group, the direct library-match anchor
#' (\code{reference_tier == "Tier1a_AnnoMatch"}) is used as the ISF parent,
#' provided it has an MS2 spectrum.  The anchor may carry any confirmed
#' ESI adduct.  Every other co-group feature is tested: if its \code{mzmed}
#' matches any fragment ion within \code{ppm_isf} and
#' \code{rt_isf_threshold} seconds, it receives \code{Tier1c_MS2Fragment}.
#'
#' \strong{Phase 2 — RT-gated MS2 fragment (Tier 1c_exp)}
#' The same anchor feature searches the full feature pool (features with no
#' Phase 1 assignment) within \code{rt_isf_threshold} seconds.
#' Matching features get \code{Tier1c_MS2Fragment_exp}.
#'
#' @param object A \code{\link{PeakAnnotationISF-class}} object.
#' @param ... Currently unused.
#'
#' @return \code{data.frame} with columns: \code{feature_name},
#'   \code{mz}, \code{rt}, \code{parent_feature}, \code{parent_mz},
#'   \code{parent_rt}, \code{matched_fragment_mz}, \code{isf_ppm},
#'   \code{isf_mda}, \code{rt_diff}, \code{Lib_group},
#'   \code{reference_tier} (\code{"1c"} or \code{"1c_exp"}),
#'   \code{annotation_phase} (\code{"Phase1_MS2frag"} or
#'   \code{"Phase2_MS2frag"}).
#'
#' @seealso \code{\link{createPeakAnnotationISF}()},
#'   \code{\link{run_isf_annotation}()},
#'   \code{\link{read_ms2_consensus}()}
#' @export
setMethod("annotateISF", "PeakAnnotationISF", function(object, ...) {
  run_isf_annotation(
    phase1_result    = object@phase1_result,
    ms2_list         = object@ms2_list,
    FT               = object@FT,
    ppm_isf          = object@ppm_isf,
    mz_threshold_isf = object@mz_threshold_isf,
    rt_isf_threshold = object@rt_isf_threshold
  )
})


#' Run MS2 fragment annotation (Phase 1 + Phase 2)
#'
#' @description
#' R wrapper around \code{isf_annotate_cpp()}.  Identifies features whose
#' \code{mzmed} matches a fragment ion of a confirmed precursor.
#'
#' @param phase1_result \code{data.frame}. Output of \code{annotateLibMatch()}.
#'   Must contain columns \code{feature_name}, \code{Lib_group},
#'   \code{reference_tier}, and \code{Compound_M}.
#' @param ms2_list \code{list}. From \code{\link{read_ms2_consensus}()}.
#' @param FT \code{data.frame}. Full MS1 feature pool.
#' @param ppm_isf Numeric. ppm tolerance for fragment matching.  Default \code{10}.
#' @param mz_threshold_isf Numeric. Absolute m/z floor (Da).  Default \code{0.005}.
#' @param rt_isf_threshold Numeric. RT window (s).  Default \code{30}.
#'
#' @return \code{data.frame} — see \code{\link{annotateISF}()}.
#' @export
run_isf_annotation <- function(
    phase1_result,
    ms2_list,
    FT,
    ppm_isf          = 10,
    mz_threshold_isf = 0.005,
    rt_isf_threshold = 30.0
) {

  mz_col <- if ("mzmed" %in% names(FT)) "mzmed" else "mz"
  rt_col <- if ("rtmed" %in% names(FT)) "rtmed" else
    if ("rt"    %in% names(FT)) "rt"    else NULL

  pool_fn <- FT$feature_name
  pool_mz <- FT[[mz_col]]
  pool_rt <- if (!is.null(rt_col)) FT[[rt_col]] else rep(NA_real_, nrow(FT))

  ft_mz_lookup <- setNames(pool_mz, pool_fn)
  ft_rt_lookup <- setNames(pool_rt, pool_fn)


  if (!"reference_tier" %in% names(phase1_result)) {
    stop("[run_isf_annotation] phase1_result must contain a 'reference_tier' ",
         "column. Check that annotateLibMatch() output is passed as phase1_result.")
  }

  anchor_rows <- phase1_result[
    !is.na(phase1_result$Lib_group) &
    !is.na(phase1_result$reference_tier) &
    grepl("^Tier1a_AnnoMatch", phase1_result$reference_tier,
          ignore.case = FALSE), , drop = FALSE]

  parents <- anchor_rows[
    anchor_rows$feature_name %in% names(ms2_list), , drop = FALSE]
  parents <- parents[!duplicated(paste(parents$feature_name,
                                       parents$Lib_group)), , drop = FALSE]

  if (nrow(parents) == 0L) {
    message("[run_isf_annotation] No anchor features with MS2 found in any ",
            "Lib_group. Returning empty result.")
    return(data.frame())
  }

  message("[run_isf_annotation] ",
          nrow(parents), " anchor features with MS2 (any adduct) | ",
          nrow(FT), " pool features")

  # a Compound_M lookup per Lib_group
  cM_lookup <- NULL
  if ("Compound_M" %in% names(phase1_result)) {
    cM_lookup <- phase1_result %>%
      dplyr::filter(!is.na(Compound_M), !is.na(Lib_group)) %>%
      dplyr::group_by(Lib_group) %>%
      dplyr::summarise(Compound_M_grp = dplyr::first(Compound_M),
                       .groups = "drop")
  }


  phase1_features <- unique(phase1_result$feature_name[
    !(grepl("_exp", phase1_result$reference_tier, fixed = TRUE) %in% TRUE)
  ])

  p1_rows <- list()
  p2_rows <- list()

  for (i in seq_len(nrow(parents))) {

    parent_fn <- parents$feature_name[i]
    lib_group <- parents$Lib_group[i]
    ms2_entry <- ms2_list[[parent_fn]]
    frag_mz   <- ms2_entry$fragment_mz
    if (length(frag_mz) == 0L) next

    #  use FT apex RT for parent
    parent_rt <- if (parent_fn %in% names(ft_rt_lookup) &&
                     !is.na(ft_rt_lookup[[parent_fn]]))
      ft_rt_lookup[[parent_fn]] else NA_real_


    compound_M <- if ("Compound_M" %in% names(parents) &&
                      !is.na(parents$Compound_M[i]))
      parents$Compound_M[i]
    else if (!is.null(cM_lookup)) {
      idx <- match(lib_group, cM_lookup$Lib_group)
      if (!is.na(idx)) cM_lookup$Compound_M_grp[idx] else NA_real_
    } else NA_real_
    guard_mz   <- if (!is.na(compound_M))
      compound_M + 1.007276   # [M+H]+ as ceiling
    else ms2_entry$precursor_mz

    # ── Phase 1 ──────────────────────────

    grp_fns  <- phase1_result$feature_name[
      phase1_result$Lib_group == lib_group &
        phase1_result$feature_name != parent_fn &
        !(grepl("_exp", phase1_result$reference_tier, fixed = TRUE) %in% TRUE)]
    grp_fns  <- unique(grp_fns)
    grp_fns  <- grp_fns[grp_fns %in% pool_fn]   # must be in FT

    grp_mz <- ft_mz_lookup[grp_fns]
    grp_rt <- ft_rt_lookup[grp_fns]

    if (!is.na(parent_rt)) {
      rt_ok  <- !is.na(grp_rt) & abs(grp_rt - parent_rt) <= rt_isf_threshold
      grp_fns <- grp_fns[rt_ok]
      grp_mz  <- grp_mz[rt_ok]
      grp_rt  <- grp_rt[rt_ok]
    }

    if (length(grp_fns) > 0L) {
      hits_p1 <- isf_annotate_cpp(
        pool_mz          = unname(grp_mz),
        pool_ids         = grp_fns,
        fragment_mz      = frag_mz,
        parent_id        = parent_fn,
        parent_mz        = guard_mz,
        parent_rt        = parent_rt,
        pool_rt          = unname(grp_rt),
        ppm_isf          = ppm_isf,
        mz_threshold_isf = mz_threshold_isf
      )
      if (nrow(hits_p1) > 0L) {
        hits_p1$Lib_group        <- lib_group
        hits_p1$reference_tier   <- "Tier1c_MS2Fragment"
        hits_p1$annotation_phase <- "Phase1_MS2frag"
        p1_rows[[length(p1_rows) + 1L]] <- hits_p1
      }
    }

    # ── Phase 2: RT-gated pool (exclude Phase 1 features) ───────────────────
    p2_mask <- !pool_fn %in% phase1_features
    p2_fn   <- pool_fn[p2_mask]
    p2_mz   <- pool_mz[p2_mask]
    p2_rt   <- pool_rt[p2_mask]

    if (!is.na(parent_rt)) {
      rt_ok  <- !is.na(p2_rt) & abs(p2_rt - parent_rt) <= rt_isf_threshold
      p2_fn  <- p2_fn[rt_ok]
      p2_mz  <- p2_mz[rt_ok]
      p2_rt  <- p2_rt[rt_ok]
    }

    if (length(p2_fn) == 0L) next

    # ── Phase 2 ISF: parent's fragments against full p2 pool ────────────────
    hits_p2 <- isf_annotate_cpp(
      pool_mz          = p2_mz,
      pool_ids         = p2_fn,
      fragment_mz      = frag_mz,
      parent_id        = parent_fn,
      parent_mz        = guard_mz,
      parent_rt        = parent_rt,
      pool_rt          = p2_rt,
      ppm_isf          = ppm_isf,
      mz_threshold_isf = mz_threshold_isf
    )

    if (nrow(hits_p2) > 0L) {
      hits_p2$Lib_group        <- NA_character_
      hits_p2$reference_tier   <- "Tier1c_MS2Fragment_exp"
      hits_p2$annotation_phase <- "Phase2_MS2frag"
      p2_rows[[length(p2_rows) + 1L]] <- hits_p2
    }
  }

  all_rows <- c(p1_rows, p2_rows)
  if (length(all_rows) == 0L) {
    message("[run_isf_annotation] No MS2 fragment matches found.")
    return(data.frame())
  }

  result <- do.call(rbind, all_rows)
  rownames(result) <- NULL


  result <- result[order(
    ifelse(result$annotation_phase == "Phase1_MS2frag", 1L, 2L),
    result$isf_ppm
  ), ]

  result <- result[!duplicated(paste(result$feature_name, result$Lib_group)), ]
  result <- result[order(result$annotation_phase, result$feature_name), ]

  message("[run_isf_annotation] Done.",
          "  Tier1c_MS2Fragment: ",
          sum(result$reference_tier == "Tier1c_MS2Fragment"),
          "  Tier1c_MS2Fragment_exp: ",
          sum(result$reference_tier == "Tier1c_MS2Fragment_exp"),
          "  total MS2 fragment features annotated.")

  result
}
