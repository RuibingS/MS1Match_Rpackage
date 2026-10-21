#' Read and Filter an MZmine Full Feature Table
#'
#' @description Imports an MZmine 3 \code{*_feature_table.csv} file, extracts
#' per-sample peak-area columns, optionally renames them to bare sample
#' names, and filters out features that were not detected in a minimum
#' fraction of samples.
#'
#' @param file \code{character(1)}. Path to the MZmine full feature table CSV
#'   (e.g. \code{"test_1_full_feature_table.csv"}).
#'
#' @param min_detect_fraction \code{numeric(1)}. Minimum fraction of samples
#'   (0–1) in which a feature must have a non-zero area to be retained.
#'   Defaults to \code{0.5} (50 \%). Set to \code{0} to keep all features.
#'
#' @param area_col \code{character(1)}. Which per-sample intensity metric to
#'   extract. One of \code{"area"} (default) or \code{"height"}.
#'
#' @param clean_names \code{logical(1)}. If \code{TRUE} (default), strip the
#'   \code{"datafile:"} prefix and \code{".mzXML:area"} (or \code{":height"})
#'   suffix from column names, leaving only the bare sample names.
#'
#' @param keep_meta \code{logical(1)}. If \code{TRUE} (default), prepend the
#'   consensus columns \code{id}, \code{mz}, and \code{rt} to the output.
#'   Set to \code{FALSE} to return a pure sample-by-feature matrix.
#'
#' @param zero_as_na \code{logical(1)}. If \code{TRUE}, replace all zero area
#'   values with \code{NA} before returning. Useful for log-transformation
#'   workflows. Defaults to \code{FALSE}.
#'
#' @param verbose \code{logical(1)}. If \code{TRUE} (default), print a short
#'   summary of features before and after filtering, plus the number of
#'   samples detected.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per retained feature.
#'   If \code{keep_meta = TRUE} the first three columns are \code{id},
#'   \code{mz}, and \code{rt}, followed by one column per sample.
#'   If \code{keep_meta = FALSE} only the sample columns are returned (with
#'   row names set to the feature \code{id} values, if the result is coerced
#'   to a \code{data.frame}).
#'
#' @details
#' MZmine 3 names every per-sample column with the pattern
#' \code{datafile:<sample_filename>.<ext>:<metric>}.  This function uses a
#' regex on the column suffix (\code{:area} or \code{:height}) to identify
#' sample columns, so it is robust to mixed file extensions (\code{.mzML},
#' \code{.mzXML}, \code{.raw}, etc.).
#'
#' Features with \code{feature_state == "ESTIMATED"} (gap-filled) are
#' treated identically to \code{"DETECTED"} features for filtering purposes;
#' only cells containing \code{NA} or \code{0} count as "not detected".
#'
#' @examples
#' \dontrun{
#' # Basic usage — keep features present in >= 50% of samples
#' MZmine_feature_table <- read_mzmine_feature_table( system.file("extdata", "MZmine_feature_table","StM14_MZmine_export","Stm14_mzmine_full_feature_table.csv", package = "MS1Match"))
#'
#' # Stricter threshold — features must appear in >= 80% of samples
#' MZmine_feature_table_strict <- read_mzmine_feature_table(
#'   system.file("extdata", "MZmine_feature_table","StM14_MZmine_export","Stm14_mzmine_full_feature_table.csv", package = "MS1Match")
#'   min_detect_fraction = 0.8
#' )
#'
#' # Keep all features, use height instead of area, return raw matrix
#' MZmine_feature_table_mat <- read_mzmine_feature_table(
#'   system.file("extdata", "MZmine_feature_table","StM14_MZmine_export","Stm14_mzmine_full_feature_table.csv", package = "MS1Match"),
#'   min_detect_fraction = 0,
#'   area_col            = "height",
#'   keep_meta           = FALSE,
#'   verbose             = FALSE
#' )
#'}
#'
#' @importFrom utils read.csv
#' @export
read_mzmine_feature_table <- function(
    file,
    min_detect_fraction = 0.5,
    area_col            = c("area", "height"),
    clean_names         = TRUE,
    keep_meta           = TRUE,
    zero_as_na          = FALSE,
    verbose             = FALSE
) {

  # ---- argument validation --------------------------------------------------

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  area_col <- match.arg(area_col)

  if (!is.numeric(min_detect_fraction) ||
      length(min_detect_fraction) != 1 ||
      min_detect_fraction < 0 ||
      min_detect_fraction > 1) {
    stop("`min_detect_fraction` must be a single number between 0 and 1.")
  }

  # ---- load -----------------------------------------------------------------

  raw <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  # ---- identify sample columns ----------------------------------------------

  # Pattern: ends with ":area" or ":height"
  pattern   <- paste0(":", area_col, "$")
  area_cols <- grep(pattern, colnames(raw), value = TRUE)

  if (length(area_cols) == 0) {
    stop(
      "No columns matching '", pattern, "' found.\n",
      "Check that `area_col` matches your MZmine export settings.\n",
      "Available column suffixes: ",
      paste(unique(sub(".*:", "", colnames(raw))), collapse = ", ")
    )
  }

  n_samples <- length(area_cols)

  if (verbose) {
    message("MZmine feature table loaded: ", nrow(raw), " features, ",
            n_samples, " samples.")
  }

  # ---- extract and coerce area matrix ---------------------------------------

  area_mat <- raw[, area_cols, drop = FALSE]

  area_mat[] <- lapply(area_mat, function(x) {
    if (is.character(x)) suppressWarnings(as.numeric(x)) else x
  })

  area_mat[is.na(area_mat)] <- 0

  # ---- filter by detection frequency ----------------------------------------

  n_features_before <- nrow(area_mat)

  n_detected <- rowSums(area_mat > 0)
  min_n      <- ceiling(min_detect_fraction * n_samples)
  keep_rows  <- n_detected >= min_n

  area_mat <- area_mat[keep_rows, , drop = FALSE]

  n_features_after <- nrow(area_mat)

  if (verbose) {
    message(
      "Detection filter (>= ", round(min_detect_fraction * 100, 1), "% of ",
      n_samples, " samples): ",
      n_features_before, " -> ", n_features_after, " features retained."
    )
  }

  # ---- optional: replace zeros with NA --------------------------------------

  if (zero_as_na) {
    area_mat[area_mat == 0] <- NA
  }

  # ---- optional: clean column names -----------------------------------------

  if (clean_names) {
    # Remove "datafile:" prefix and ".<ext>:<metric>" suffix
    clean_pattern      <- paste0("^datafile:|\\.[^.]+:", area_col, "$")
    colnames(area_mat) <- gsub(clean_pattern, "", colnames(area_mat))
  }
  # ---- assemble output -------------------------------------------------------

  # Compute per-feature mean and median across samples
  intensity_mean   <- rowMeans(area_mat, na.rm = TRUE)
  intensity_median <- apply(area_mat, 1, median, na.rm = TRUE)

  if (keep_meta) {
    meta_cols <- c("id", "mz", "rt")
    meta_cols <- intersect(meta_cols, colnames(raw))
    meta      <- raw[keep_rows, meta_cols, drop = FALSE]

    # Rename to XCMS conventions
    meta[["feature_name"]] <- sprintf("FT%04d", as.integer(meta[["id"]]))
    meta[["id"]] <- NULL
    meta <- meta[, c("feature_name", setdiff(names(meta), "feature_name")), drop = FALSE]
    names(meta)[names(meta) == "mz"] <- "mzmed"
    if ("rt" %in% names(meta)) {
      meta[["rt"]] <- meta[["rt"]] * 60
      names(meta)[names(meta) == "rt"] <- "rtmed"
    }

    meta[["intensity_mean"]]   <- intensity_mean
    meta[["intensity_median"]] <- intensity_median

    out <- cbind(meta, area_mat)

  } else {
    out <- data.frame(
      intensity_mean   = intensity_mean,
      intensity_median = intensity_median,
      area_mat,
      check.names = FALSE
    )
    if ("id" %in% colnames(raw)) {
      rownames(out) <- raw[["id"]][keep_rows]
    }
  }

  if (requireNamespace("tibble", quietly = TRUE)) {
    out <- tibble::as_tibble(out)
  }

  out
}
