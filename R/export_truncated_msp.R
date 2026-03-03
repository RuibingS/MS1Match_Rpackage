#' Export a parsed or truncated MSP library to file
#'
#' @param lib_list List of compounds (as returned by parse_any_library_cpp)
#' @param file_path Path to the output .msp file
#' @param overwrite Whether to overwrite existing file (default TRUE)
#'
#' @return Invisibly returns the output file path
#' @export
export_msp <- function(lib_list, file_path, overwrite = TRUE) {
  if (file.exists(file_path) && !overwrite)
    stop("File already exists: ", file_path)

  con <- file(file_path, open = "wt", encoding = "UTF-8")

  on.exit(close(con))

  for (cmp in lib_list) {
    # ---- Write metadata ----
    for (nm in names(cmp)) {
      if (nm == "Peaks") next
      val <- cmp[[nm]]
      if (length(val) == 1 && !is.null(val) && !is.na(val)) {
        writeLines(sprintf("%s: %s", nm, val), con)
      }
    }

    # ---- Write peaks ----
    if (!is.null(cmp$Peaks) && is.data.frame(cmp$Peaks)) {
      peaks <- cmp$Peaks
      writeLines(sprintf("Num Peaks: %d", nrow(peaks)), con)

      # Format m/z and intensity pairs
      peak_lines <- sprintf("%.4f\t%.0f", peaks$mz, peaks$intensity)
      writeLines(peak_lines, con)
      writeLines("", con)
    }
  }

  message("MSP library exported to: ", normalizePath(file_path))
  invisible(file_path)
}
