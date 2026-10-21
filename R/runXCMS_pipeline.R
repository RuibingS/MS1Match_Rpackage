# =============================================================================
# XCMS-only pipeline: raw files to feature table
# =============================================================================
# Step 0  runXCMSProcessing()   raw mzXML/mzML, XcmsExperiment
# Step 1  getChromPeaks()       chromPeaks as plain data.frame
# Step 2  linkCPtoFeatures()    chromPeaks + feature_name (long format)
# Step 3  buildFeatureTable()   feature defs + per-sample intensity matrix
# Step 4  processFeatureTable() sort by RT, feature_name col, intensity_median
# Step 5  runxcmsFT()         wraps steps 1-4 into one call
# =============================================================================


# -----------------------------------------------------------------------------
# Step 0 — XCMS processing
# -----------------------------------------------------------------------------

#' Run XCMS preprocessing
#'
#' @param files_dir Path to directory containing mzXML or mzML files.
#' @param cwp \code{\link[xcms]{CentWaveParam}}. NULL = defaults for
#'   Bruker maXis HD / UHPLC / positive ESI.
#' @param pdp \code{\link[xcms]{PeakDensityParam}}. NULL = defaults.
#' @param obi \code{\link[xcms]{ObiwarpParam}} for RT alignment.
#'   Pass NULL to skip (single file or pre-aligned data).
#' @param sample_name  Character vector, one name per file. NULL = basenames.
#' @param sample_group Character vector of group labels, parallel to
#'   \code{sample_name}. Example: \code{c(rep("WT",3), rep("KO",3))}.
#' @param BPPARAM \code{\link[BiocParallel]{BiocParallelParam}}.
#'   Default \code{SerialParam()}. Parallelise with
#'   \code{SnowParam(workers = 4)}.
#' @return Named list: \code{x_exp} (XcmsExperiment) and \code{pd}
#'   (sample metadata DataFrame).
#' @export
runXCMSProcessing <- function(
    files_dir,
    cwp          = NULL,
    pdp          = NULL,
    obi          = ObiwarpParam(),
    sample_name  = NULL,
    sample_group = NULL,
    BPPARAM      = BiocParallel::SerialParam()
) {

  stopifnot(is.character(files_dir), length(files_dir) == 1L, dir.exists(files_dir))

  MS_files <- list.files(files_dir, pattern = "\\.(mzXML|mzML)$",
                         full.names = TRUE, ignore.case = TRUE)
  if (!length(MS_files)) stop("No mzXML/mzML files found in: ", files_dir)
  n <- length(MS_files)

  if (is.null(sample_name))  sample_name  <- tools::file_path_sans_ext(basename(MS_files))
  if (is.null(sample_group)) sample_group <- rep("group1", n)

  if (length(sample_name)  != n) stop("`sample_name` must have ", n, " entries")
  if (length(sample_group) != n) stop("`sample_group` must have ", n, " entries")

  if (is.null(cwp))
    cwp <- CentWaveParam(
      peakwidth = c(5, 20), ppm = 5, snthresh = 10,
      mzdiff = 0.001, prefilter = c(3, 1000), noise = 500,
      verboseBetaColumns = TRUE
    )

  if (is.null(pdp))
    pdp <- PeakDensityParam(
      sampleGroups = sample_group, bw = 5, minFraction = 0.5, binSize = 0.01
    )

  pd <- S4Vectors::DataFrame(sample_name = sample_name, sample_group = sample_group)

  message("[1/4] Reading raw data (", n, " files)")
  raw_exp <- MsExperiment::readMsExperiment(spectraFiles = MS_files, sampleData = pd)

  message("[2/4] Detecting peaks (CentWave)")
  x_exp <- findChromPeaks(raw_exp, param = cwp, msLevel = 1L, BPPARAM = BPPARAM)

  if (n > 1 && !is.null(obi)) {
    message("[3/4] Aligning retention times (Obiwarp)")
    x_exp <- adjustRtime(x_exp, obi)
  } else {
    message("[3/4] Skipping RT alignment")
  }

  message("[4/4] Grouping features + filling gaps")
  x_exp <- groupChromPeaks(x_exp, param = pdp)
  x_exp <- fillChromPeaks(x_exp, ChromPeakAreaParam())

  message("Done. ", nrow(featureDefinitions(x_exp)), " features detected.")
  list(x_exp = x_exp, pd = pd)
}


# -----------------------------------------------------------------------------
# Step 1 — xcms raw chromPeaks
# -----------------------------------------------------------------------------

#' Extract chromPeaks as a plain data.frame
#' @param x_exp XcmsExperiment after fillChromPeaks.
#' @return data.frame, one row per detected peak, with \code{sample_name}.
#' @export
getChromPeaks <- function(x_exp) {
  cp <- as.data.frame(chromPeaks(x_exp))
  cp$sample_name <- colnames(featureValues(x_exp, value = "into"))[cp$sample]
  cp
}


# -----------------------------------------------------------------------------
# Step 2 — chromPeaks linked to features (long format)
# -----------------------------------------------------------------------------

#' Link chromPeaks to feature IDs
#' @param x_exp XcmsExperiment.
#' @return data.frame: all chromPeak columns + \code{feature_name} +
#'   \code{sample_name}. Unassigned peaks have \code{NA} in \code{feature_name}.
#' @export
linkCPtoFeatures <- function(x_exp) {
  cp   <- as.data.frame(chromPeaks(x_exp))
  fdef <- featureDefinitions(x_exp)

  feature_name <- rep(NA_character_, nrow(cp))
  for (fi in seq_len(nrow(fdef)))
    feature_name[fdef$peakidx[[fi]]] <- rownames(fdef)[fi]

  cp$feature_name  <- feature_name
  cp$sample_name <- colnames(featureValues(x_exp, value = "into"))[cp$sample]
  cp
}


# -----------------------------------------------------------------------------
# Step 3 — Build feature table
# -----------------------------------------------------------------------------

#' Build simple feature table
#'
#'  @description Feature definitions cbind'd with the per-sample \code{into} intensity matrix.
#' @param x_exp XcmsExperiment after fillChromPeaks.
#' @return data.frame, one row per feature, rownames = feature IDs.
#' @export
buildFeatureTable <- function(x_exp) {
  fdef <- as.data.frame(featureDefinitions(x_exp))
  fint <- featureValues(x_exp, value = "into")
  cbind(fdef, fint)
}


# -----------------------------------------------------------------------------
# Step 4 — Post-process: rename and reorder
# -----------------------------------------------------------------------------

#' Post-process feature table
#'
#' \enumerate{
#'   \item Moves rownames to a leading \code{feature_name} column.
#'   \item Sorts rows by \code{rtmed}.
#'   \item Appends \code{intensity_median} (median of per-sample intensities).
#' }
#' @param ft        data.frame from \code{buildFeatureTable()}.
#' @param x_exp     XcmsExperiment used to produce \code{ft}.
#' @param rt_range  Optional numeric vector of length 2 for RT filtering (e.g. \code{c(60, 1200)}).
#' @param sort_by   Optional column name to sort by (e.g. \code{"rtmed"}, \code{"mzmed"}, \code{"feature_name"}). Defaults to \code{"rtmed"}.
#' @return Processed data.frame.
#' @export
processFeatureTable <- function(ft, x_exp, rt_range = NULL, sort_by = "rtmed") {
  if (!"feature_name" %in% colnames(ft)) {
    ft <- cbind(feature_name = rownames(ft), ft, stringsAsFactors = FALSE)
  }
  rownames(ft) <- NULL

  if (!is.null(rt_range)) {
    if (!is.numeric(rt_range) || length(rt_range) != 2)
      stop("rt_range must be a numeric vector of length 2.")
    rt_col <- if ("rtmed" %in% colnames(ft)) "rtmed" else if ("rt" %in% colnames(ft)) "rt" else
      stop("No RT column ('rtmed' or 'rt') found for rt_range filtering.")
    ft <- ft[ft[[rt_col]] >= rt_range[1] & ft[[rt_col]] <= rt_range[2], ]
  }

  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(ft))
      stop(paste0("sort_by column '", sort_by, "' not found in feature table."))
    ft <- ft[order(ft[[sort_by]]), ]
  }

  fdef_cols <- c("feature_name", colnames(featureDefinitions(x_exp)))
  int_cols  <- setdiff(colnames(ft), fdef_cols)
  if (!length(int_cols)) stop("No intensity columns found.")

  int_mat <- as.matrix(ft[, int_cols, drop = FALSE])
  storage.mode(int_mat) <- "numeric"
  ft$intensity_median <- apply(int_mat, 1, median, na.rm = TRUE)
  ft$intensity_mean   <- apply(int_mat, 1, mean,   na.rm = TRUE)
  ft
}
# -----------------------------------------------------------------------------
# Step 5 — pipeline function
# -----------------------------------------------------------------------------

#' Pipeline function feature table pipeline
#'
#' @param files_dir  Passed to \code{runXCMSProcessing()}.
#' @param cwp,pdp,obi,sample_name,sample_group,BPPARAM
#'   Passed to \code{runXCMSProcessing()}.
#' @return Named list:
#'   \describe{
#'     \item{result}{Output of \code{runXCMSProcessing()} (\code{x_exp}, \code{pd}).}
#'     \item{cp}{Raw chromPeaks data.frame.}
#'     \item{cp_linked}{chromPeaks with \code{feature_name} column.}
#'     \item{ft}{Processed feature table (sorted, \code{intensity_median}).}
#'   }
#' @examples
#' \dontrun{
#' out      <- runxcmsFT("data/mzXML/",
#'                         sample_name  = c("WT_1","WT_2","WT_3","KO_1","KO_2","KO_3"),
#'                         sample_group = c(rep("WT",3), rep("KO",3)))
#' out_ft  <- out$ft
#' }
#' @export
runxcmsFT <- function(
    files_dir,
    cwp          = NULL,
    pdp          = NULL,
    obi          = ObiwarpParam(),
    sample_name  = NULL,
    sample_group = NULL,
    BPPARAM      = BiocParallel::SerialParam()
) {

  result    <- runXCMSProcessing(files_dir, cwp, pdp, obi,
                                 sample_name, sample_group, BPPARAM)
  cp        <- getChromPeaks(result$x_exp)
  cp_linked <- linkCPtoFeatures(result$x_exp)
  ft        <- processFeatureTable(buildFeatureTable(result$x_exp), result$x_exp)

  message(sprintf(
    "Pipeline complete: %d features | %d samples | %d chromPeaks (%d assigned)",
    nrow(ft), length(sample_name %||% basename(list.files(files_dir, pattern = "mzXML|mzML", ignore.case = TRUE))),
    nrow(cp), sum(!is.na(cp_linked$feature_name))
  ))

  list(result = result, cp = cp, cp_linked = cp_linked, ft = ft)
}



# -----------------------------------------------------------------------------
# Export — write feature table to CSV
# -----------------------------------------------------------------------------

#' Export feature table to CSV
#'
#' Writes the processed feature table to a CSV file. The \code{peakidx} list
#' column is collapsed to a single semicolon-separated string per row so the
#' file remains rectangular and readable in Excel / any CSV reader.
#'
#' Any other residual list columns (unexpected) are similarly collapsed as a
#' safety measure.
#'
#' @param ft     data.frame from \code{processFeatureTable()} or \code{runxcmsFT()$ft}.
#' @param file   Output file path, e.g. \code{"results/feature_table.csv"}.
#'   The parent directory is created automatically if it does not exist.
#' @param sep    Field separator. Default \code{","}.
#' @param ...    Additional arguments passed to \code{write.csv()}.
#' @return \code{ft} invisibly (so the call can be piped).
#' @examples
#' \dontrun{
#' out <- runxcmsFT("data/mzXML/", sample_name = c("S1","S2","S3"))
#' exportFeatureTable(out$ft, "results/feature_table.csv")
#' }
#' @export
exportFeatureTable <- function(ft, file, sep = ",", ...) {

  stopifnot(is.data.frame(ft))
  if (!nzchar(file)) stop("`file` must be a non-empty path string.")

  dir <- dirname(file)
  if (!dir.exists(dir)) {
    message("Creating directory: ", dir)
    dir.create(dir, recursive = TRUE)
  }

  # Collapse every list column to a semicolon-separated string
  list_cols <- vapply(ft, is.list, logical(1L))
  if (any(list_cols)) {
    message("Collapsing list column(s): ",
            paste(names(ft)[list_cols], collapse = ", "))
    ft[list_cols] <- lapply(ft[list_cols], function(col) {
      vapply(col,
             function(x) paste(x, collapse = ";"),
             character(1L))
    })
  }

  write.csv(ft, file = file, row.names = FALSE, ...)
  message(sprintf("Feature table written: %s  [%d rows x %d cols]",
                  file, nrow(ft), ncol(ft)))

  invisible(ft)
}

