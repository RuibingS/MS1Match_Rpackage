#' Flatten annotate result (PeakAnnotation) to a data.frame
#'
#' @param res_list A list returned by annotate(PeakAnnotation_obj)
#' @return A data.frame with one row per list element
#' @export
annotations_df <- function(res_list) {
  if (!is.list(res_list)) stop("`res_list` must be a list.")

  rows <- lapply(res_list, function(x) {
    # each x is a small named list of scalars: feature_name, mz, rt, PI_match, adduct_match, NL_match, iso_annotation
    as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE)
  })

  # fast path: try a direct rbind
  same_cols <- length(unique(vapply(rows, function(df) paste(names(df), collapse = "|"), ""))) == 1
  if (same_cols) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  # fallback: harmonize columns if some elements missed fields
  all_cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(df) {
    miss <- setdiff(all_cols, names(df))
    if (length(miss)) df[miss] <- NA
    df[all_cols]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
