# ============================================================
# constructor_createPeakAnnotationISF.R
# ============================================================

#' Create a PeakAnnotationISF object
#'
#' @description
#' User-facing constructor for \code{\link{PeakAnnotationISF-class}}.
#'
#' Tier naming:
#' \itemize{
#'   \item \strong{Tier 1c} — MS2 fragment match, intra-group (Phase 1)
#'   \item \strong{Tier 1c_ext} — MS2 fragment match, RT-gated features (Phase 2)
#' }
#' Priority: 1a > 1b > 1c > 1d > 2.
#' MS2 tier (1c) only assigned when no higher-priority tier already fired.
#'
#' Parent features are identified by \code{parent_col}/\code{parent_val}:
#' \itemize{
#'   \item \code{parent_val = NULL} (default): any non-empty, non-NA string
#'     in \code{parent_col} is a parent (e.g. \code{PI_match} column
#'     containing compound name strings).
#'   \item Otherwise: exact match \code{parent_col == parent_val}.
#' }
#'
#' @param phase1_result \code{data.frame}. Output of \code{annotateLibMatch()}.
#'   Required columns: \code{feature_name}, \code{Lib_group}, \code{mz},
#'   and \code{parent_col}.  Optional: \code{rt}.
#' @param ms2_list \code{list}. From \code{\link{read_ms2_consensus}()}.
#' @param FT \code{data.frame}. Full MS1 feature tabel.
#'   Required: \code{feature_name}, \code{mzmed}/\code{mz}.
#'   Optional: \code{rtmed}/\code{rt}.
#' @param parent_col Character. Column identifying parent features.
#'   Default \code{"PI_match"}.
#' @param parent_val Value marking a parent.  \code{NULL} = any non-empty
#'   string.  Default \code{NULL}.
#' @param ppm_isf Numeric. ppm tolerance.  Default \code{10}.
#' @param mz_threshold_isf Numeric. Absolute m/z floor (Da).
#'   Default \code{0.005}.
#' @param rt_isf_threshold Numeric. RT window (seconds).  Default \code{5}.
#'
#' @return A \code{\link{PeakAnnotationISF-class}} object.
#'
#' @examples
#' \dontrun{
#' ms2_list <- read_ms2_consensus(ms2_consensus, feature_table = FT)
#'
#' isf_obj <- createPeakAnnotationISF(
#'   phase1_result    = pool_annotate_res_ms1,
#'   ms2_list         = ms2_list,
#'   FT          = FT,
#'   parent_col       = "PI_match",
#'   ppm_isf          = 10,
#'   rt_isf_threshold = 5.0
#' )
#' isf_result <- annotateISF(isf_obj)
#' }
#'
#' @seealso \code{\link{annotateISF}()},
#'   \code{\link{PeakAnnotationISF-class}},
#'   \code{\link{read_ms2_consensus}()}
#' @export
createPeakAnnotationISF <- function(
    phase1_result,
    ms2_list,
    FT,
    parent_col       = "PI_match",
    parent_val       = NULL,
    ppm_isf          = 10,
    mz_threshold_isf = 0.005,
    rt_isf_threshold = 5.0
) {
  # ── Input validation ───────────────────────────────────────────────────────
  if (!is.data.frame(phase1_result))
    stop("`phase1_result` must be a data.frame (output of annotateLibMatch()).")

  required_p1 <- c("feature_name", "Lib_group", "mz")
  missing_p1  <- setdiff(required_p1, names(phase1_result))
  if (length(missing_p1) > 0L)
    stop("`phase1_result` is missing required columns: ",
         paste(missing_p1, collapse = ", "))

  if (!parent_col %in% names(phase1_result))
    stop("`parent_col` '", parent_col,
         "' not found in `phase1_result`. Available columns: ",
         paste(names(phase1_result), collapse = ", "))

  if (!is.list(ms2_list) || length(ms2_list) == 0L)
    stop("`ms2_list` must be a non-empty named list ",
         "from read_ms2_consensus().")

  if (!is.data.frame(FT))
    stop("`FT` must be a data.frame.")

  if (!any(c("mzmed", "mz") %in% names(FT)))
    stop("`FT` must have a column 'mzmed' or 'mz'.")

  if (!"feature_name" %in% names(FT))
    stop("`FT` must have a column 'feature_name'.")

  # ── Count parent features ──────────────────────────────────────────────────
  parent_mask <- .make_parent_mask(phase1_result[[parent_col]], parent_val)
  parents     <- phase1_result[parent_mask, , drop = FALSE]
  n_with_ms2  <- sum(parents$feature_name %in% names(ms2_list))

  message("[createPeakAnnotationISF]",
          " Parent features (", parent_col,
          ifelse(is.null(parent_val), " non-empty",
                 paste0("=", parent_val)), "): ",
          nrow(parents),
          " | with MS2: ", n_with_ms2,
          " | features: ", nrow(FT))

  if (n_with_ms2 == 0L)
    warning("[createPeakAnnotationISF] No parent features have MS2. ",
            "Check that feature_name in phase1_result matches ",
            "names(ms2_list).\n",
            "  phase1_result examples: ",
            paste(head(parents$feature_name, 3), collapse = ", "), "\n",
            "  ms2_list examples: ",
            paste(head(names(ms2_list), 3), collapse = ", "))

  new("PeakAnnotationISF",
      phase1_result    = phase1_result,
      ms2_list         = ms2_list,
      FT          = FT,
      parent_col       = as.character(parent_col),
      parent_val       = parent_val,
      ppm_isf          = as.numeric(ppm_isf),
      mz_threshold_isf = as.numeric(mz_threshold_isf),
      rt_isf_threshold = as.numeric(rt_isf_threshold)
  )
}


# ── Internal: build parent logical mask ──────────────────────────────────────
#' @keywords internal
.make_parent_mask <- function(col_values, parent_val) {
  if (is.null(parent_val)) {
    !is.na(col_values) &
    nchar(trimws(as.character(col_values))) > 0
  } else {
    !is.na(col_values) & col_values == parent_val
  }
}
