# =============================================================================
# Import a pre-processed feature table (CSV / data.frame) for MS1Match
# =============================================================================
# workflow for software other than xcms / MZmine:
#   CSV -> read_feature_table() -> MS1Match workflow -> anno_result
# Works for a plain feature table (no grouping) as well as one that also
# carries a group / cluster index.
# =============================================================================

#' Read a pre-processed feature table for MS1Match
#'
#' @description
#' Imports a feature table from any software (CSV file or \code{data.frame})
#' into the format MS1Match uses: \code{feature_name}, \code{mzmed},
#' \code{rtmed} (seconds), \code{intensity_mean}, \code{intensity_median},
#' followed by one column per sample. A group / cluster column is not
#' required. Any other columns (e.g. such a group index) are kept unchanged.
#'
#' @param x Path to a CSV file, or a \code{data.frame}.
#' @param mz_col,rt_col Names of the m/z and RT columns.
#' @param rt_unit Unit of \code{rt_col}: \code{"sec"} (default) or
#'   \code{"min"}. Converted to seconds.
#' @param sample_cols Per-sample intensity columns: a character vector of
#'   column names, or a single regular expression matched against the column
#'   names (e.g. \code{"\\\\.mzML$"}).
#' @param id_col Optional feature-id column, used as \code{feature_name}.
#'   Plain numbers (1, 2, ...) are converted to \code{FT0001}, \code{FT0002},
#'   ... other ids are kept
#'   as they are. If \code{NULL}, ids \code{FT0001}, ... are assigned after
#'   sorting by RT then m/z.
#'
#' @return \code{data.frame}, one row per retained feature.
#' @export
read_feature_table <- function(x, mz_col = "mz", rt_col = "rt",
                               rt_unit = c("sec", "min"),
                               sample_cols, id_col = NULL) {

  rt_unit <- match.arg(rt_unit)
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("read_feature_table: file not found: ", x)
    x <- utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE)
  }
  stopifnot(is.data.frame(x))
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  # --- required columns -------------------------------------------------------
  missing <- setdiff(c(mz_col, rt_col, id_col), names(x))
  if (length(missing) > 0L)
    stop("read_feature_table: column(s) not found: ",
         paste(missing, collapse = ", "),
         "\n  Available columns: ", paste(names(x), collapse = ", "))

  if (missing(sample_cols))
    stop("read_feature_table: give the per-sample intensity columns ",
         "(sample_cols = c(...) or a regular expression).")
  if (length(sample_cols) == 1L && !sample_cols %in% names(x))
    sample_cols <- grep(sample_cols, names(x), value = TRUE)
  sample_cols <- setdiff(sample_cols, c(mz_col, rt_col, id_col))
  other_names <- setdiff(sample_cols, names(x))
  if (length(sample_cols) == 0L || length(other_names) > 0L)
    stop("read_feature_table: sample columns not found",
         if (length(other_names)) paste0(": ", paste(other_names, collapse = ", ")) else "", ".")

  # --- m/z, RT ----------------------------------------------------------------
  mzmed <- as.numeric(x[[mz_col]])
  rtmed <- as.numeric(x[[rt_col]]) * (if (rt_unit == "min") 60 else 1)
  .warn_rt_unit(rtmed, rt_unit, "read_feature_table")
  if (anyNA(mzmed) || anyNA(rtmed))
    stop("read_feature_table: m/z or RT contains missing / non-numeric values.")

  # --- intensities -------------------------------------------------------------
  int <- x[, sample_cols, drop = FALSE]
  int[] <- lapply(int, function(v) suppressWarnings(as.numeric(v)))
  int[is.na(int)] <- 0


  # --- feature ids ------------------------------------------------------------
  if (!is.null(id_col)) {
    feature_name <- .as_feature_name(x[[id_col]])
    if (anyDuplicated(feature_name) || anyNA(feature_name))
      stop("read_feature_table: '", id_col, "' has missing or duplicated values.")
  } else {
    ord <- order(rtmed, mzmed)
    x <- x[ord, , drop = FALSE]; int <- int[ord, , drop = FALSE]
    mzmed <- mzmed[ord]; rtmed <- rtmed[ord]
    feature_name <- .as_feature_name(seq_len(nrow(x)))
  }

  # --- assemble ---------------------------------------------------------------
  other <- setdiff(names(x), c(mz_col, rt_col, id_col, sample_cols))
  out <- data.frame(
    feature_name     = feature_name,
    mzmed            = mzmed,
    rtmed            = rtmed,
    intensity_mean   = rowMeans(int),
    intensity_median = apply(int, 1, stats::median),
    x[, other, drop = FALSE],
    int,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  message(sprintf("read_feature_table: %d features, %d samples%s.",
                  nrow(out), length(sample_cols),
                  if (is.null(id_col)) "; ids assigned by RT order" else ""))
  out
}


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Feature ids: plain numbers become "FT" + at least 4 digits (1 -> FT0001,
# 12345 -> FT12345), the same rule as read_mzmine_feature_table(). The rule
# depends only on the id itself, so a feature table and a grouping table give
# the same name even if they contain different numbers of rows. Any other ids
# (e.g. "FT0001", "M123T45") are kept unchanged.
.as_feature_name <- function(id) {
  if (is.numeric(id))          # avoid "1e+05" for large numbers
    id <- ifelse(is.na(id), NA_character_,
                 format(id, scientific = FALSE, trim = TRUE))
  id  <- trimws(as.character(id))
  out <- id
  num <- !is.na(id) & grepl("^[0-9]+$", id)
  if (any(num) && all(num | is.na(id)))
    out[num] <- sprintf("FT%04d", as.integer(id[num]))
  out
}

# Warn if RT values look like the wrong unit
.warn_rt_unit <- function(rt_sec, rt_unit, fn) {
  mx <- suppressWarnings(max(rt_sec, na.rm = TRUE))
  if (is.finite(mx) && rt_unit == "sec" && mx < 60)
    warning(fn, ": all RT values are below 60 s; are they in minutes? ",
            "Set rt_unit = \"min\".", call. = FALSE)
}
