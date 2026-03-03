#' Plot mirror plots for a list of spectral matches
#'
#' @param object A \code{DirectLibMatch} object containing the compound library
#' @param matches A \code{list} of match results, as returned by \code{libAnnotate()}
#' @param output_dir Optional directory to save PNG files, otherwise plots are shown interactively
#' @param ppm Numeric value for ppm tolerance used in mirror plot matching
#' @import progressr
#' @export
setGeneric("plotMirrorSpectra", function(object, matches, output_dir = NULL, ppm = 10) {
  standardGeneric("plotMirrorSpectra")
})


#' Plot a mirror plot between a query feature and its reference compound
#'
#' @param match A single entry (list) from a DirectLibMatch result.
#' @param object The DirectLibMatch object containing @compounds_list.
#' @param output_dir Optional directory to save the mirror plot as PNG.
#' @param ppm m/z tolerance (default 10 ppm).
#' @return Invisibly returns the ggplot object.
#' @export
#'
plotMirrorSpectraSingle <- function(match, object, output_dir = NULL, ppm = 10) {
  if (is.null(match$Compound_name)) {
    warning("No compound name found in match.")
    return(invisible(NULL))
  }

  comp_name_ref <- match$Compound_name

  # Locate compound (optional)
  compound_idx <- if (!is.null(match$Compound_index) &&
                      length(object@compounds_list) >= match$Compound_index) {
    match$Compound_index
  } else {
    which(vapply(object@compounds_list, function(x)
      isTRUE(x$Name == comp_name_ref), logical(1)))
  }

  if (length(compound_idx) == 0) {
    warning("Compound not found in compound list: ", comp_name_ref)
    return(invisible(NULL))
  }

  # --- Use precomputed scores ---
  inner_cosine <- match$inner_cosine
  outer_cosine <- match$outer_cosine
  match_score  <- match$Match_peak_score
  matched_peaks <- match$Matched_peaks
  total_peaks <- match$Total_library_peaks

  # --- Prepare m/z and intensity ---
  query_mz <- as.numeric(match$Matched_features_mz)
  query_int <- as.numeric(match$Matched_features_int)
  ref_mz <- as.numeric(match$Compound_peak_mz_full)
  ref_int <- as.numeric(match$Compound_peak_int_full)

  if (length(query_mz) == 0 || all(is.na(query_mz)) ||
      length(ref_mz) == 0 || all(is.na(ref_mz))) {
    warning("Missing peak data for: ", comp_name_ref)
    return(invisible(NULL))
  }

  # Normalize intensities
  query_int <- 1000 * query_int / max(query_int, na.rm = TRUE)
  ref_int <- 1000 * ref_int / max(ref_int, na.rm = TRUE)

  # --- Construct proper single-length Spectra objects ---
  query_spectrum <- Spectra::Spectra(
    S4Vectors::DataFrame(
      msLevel = 1L,
      polarity = 1L,
      mz = I(list(query_mz)),        # 'I' ensures it’s one list-column entry
      intensity = I(list(query_int))
    ),
    backend = Spectra::MsBackendMemory()
  )

  ref_spectrum <- Spectra::Spectra(
    S4Vectors::DataFrame(
      msLevel = 1L,
      polarity = 1L,
      mz = I(list(ref_mz)),
      intensity = I(list(ref_int))
    ),
    backend = Spectra::MsBackendMemory()
  )

  # --- Sanity check ---
  if (length(query_spectrum) != 1L || length(ref_spectrum) != 1L) {
    stop("Both query and reference spectra must have length 1.")
  }

  # --- Plot function ---
  plot_fun <- function() {
    Spectra::plotSpectraMirror(
      query_spectrum,
      ref_spectrum,
      tolerance = ppm * 1e-6,
      main = paste0("Mirror Plot: ", comp_name_ref)
    )

    graphics::mtext(
      paste0("Outer Cosine: ", round(outer_cosine, 3),
             " | Inner Cosine: ", round(inner_cosine, 3),
             " | Match score: ", round(match_score, 3),
             " | Matched Peaks: ", matched_peaks, "/", total_peaks),
      line = 0.5, cex = 0.9
    )
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    file_path <- file.path(
      output_dir,
      paste0("MirrorPlot_", gsub("[^A-Za-z0-9_]", "_", comp_name_ref), ".png")
    )
    grDevices::png(file_path, width = 1200, height = 800, res = 300)
    plot_fun()
    grDevices::dev.off()
    message("Saved mirror plot to: ", file_path)
  } else {
    plot_fun()
  }

  invisible(list(
    compound_name = comp_name_ref,
    inner_cosine = inner_cosine,
    outer_cosine = outer_cosine,
    match_score = match_score,
    matched_peaks = matched_peaks,
    total_peaks = total_peaks
  ))
}

#old version
# plotMirrorSpectraSingle <- function(match, object, output_dir = NULL, ppm = 10) {
#   if (is.null(match$Compound_name)) return(NULL)
#
#   comp_name_ref <- match$Compound_name
#
#   # Find compound in reference
#   compound_idx <- which(vapply(object@compounds_list, function(x)
#     isTRUE(x$Name == comp_name_ref), logical(1)))
#
#   if (length(compound_idx) == 0) {
#     warning("Compound not found in compound list: ", comp_name_ref)
#     return(invisible(NULL))
#   }
#
#   compound <- object@compounds_list[[compound_idx]]
#
#   ref_peaks <- data.frame(
#     mz = match$Compound_peak_mz,
#     rel_intensity = match$Compound_peak_int
#   )
#   ref_peaks$rel_intensity <- 1000 * ref_peaks$rel_intensity / max(ref_peaks$rel_intensity, na.rm = TRUE)
#
#   query_df <- data.frame(
#     mzmed = match$Matched_features_mz,
#     rel_intensity = match$Matched_features_int
#   )
#   query_df$rel_intensity <- 1000 * query_df$rel_intensity / max(query_df$rel_intensity, na.rm = TRUE)
#
#   # Create Spectra object
#   spd <- S4Vectors::DataFrame(
#     msLevel = rep(1L, 2),
#     polarity = rep(1L, 2),
#     mz = IRanges::NumericList(query_df$mzmed, ref_peaks$mz),
#     intensity = IRanges::NumericList(query_df$rel_intensity, ref_peaks$rel_intensity)
#   )
#   sps <- Spectra::Spectra(spd)
#
#   inner_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "inner", ppm = ppm)[1, 2], error = function(e) NA)
#   outer_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "outer", ppm = ppm)[1, 2], error = function(e) NA)
#
#   filepath <- NULL
#   if (!is.null(output_dir)) {
#     if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#     filename <- paste0(gsub("[^A-Za-z0-9_]", "_", comp_name_ref), ".png")
#     filepath <- file.path(output_dir, filename)
#
#     grDevices::png(filepath, width = 8, height = 6, units = "in", res = 300)
#     on.exit(grDevices::dev.off(), add = TRUE)
#   }
#
#   Spectra::plotSpectraMirror(sps[[1]], sps[[2]], main = comp_name_ref, ppm = ppm)
#   graphics::mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
#                          " | Inner Cosine: ", round(inner_cosine, 3),
#                          " | Match score: ", round(match$Match_peak_score, 3)),
#                   line = 0.5, cex = 0.9)
#
#   invisible(list(
#     compound_name = comp_name_ref,
#     inner_cosine = inner_cosine,
#     outer_cosine = outer_cosine,
#     plot_path = filepath
#   ))
# }

#' @export
setMethod("plotMirrorSpectra", signature(object = "DirectLibMatch"),
          function(object, matches, output_dir = NULL, ppm = 10) {
            if (!requireNamespace("progressr", quietly = TRUE)) {
              stop("Please install the 'progressr' package to enable progress bars.")
            }

            if (!is.null(output_dir) && !dir.exists(output_dir)) {
              dir.create(output_dir, recursive = TRUE)
            }

            results <- list()
            n <- length(matches)

            progressr::with_progress({
              p <- progressr::progressor(along = matches)

              for (i in seq_along(matches)) {
                match <- matches[[i]]
                p(message = paste("Processing", i, "/", n))

                if (is.null(match$Compound_name)) next

                res <- plotMirrorSpectraSingle(match = match, object = object, output_dir = output_dir, ppm = ppm)
                if (!is.null(res)) {
                  results[[match$Compound_name]] <- res
                }
              }
            })

            invisible(results)
          })


# plotMirrorSpectraSingle <- function(match, object, output_dir = NULL, ppm = 10) {
#   comp_name_ref <- match$Compound_name
#
#   # Find compound in compound list
#   compound_idx <- which(vapply(object@compounds_list, function(x)
#     isTRUE(x$Name == comp_name_ref), logical(1)))
#
#   if (length(compound_idx) == 0) {
#     warning("Compound not found in compound list: ", comp_name_ref)
#     return(invisible(NULL))
#   }
#
#   compound <- object@compounds_list[[compound_idx]]
#
#   # Reference peaks (bottom spectrum)
#   ref_peaks <- data.frame(
#     mz = match$Compound_peak_mz_full,
#     rel_intensity = match$Compound_peak_int_full
#   )
#   if (nrow(ref_peaks) == 0 || all(is.na(ref_peaks$rel_intensity))) {
#     warning("No reference peaks found for compound: ", comp_name_ref)
#     return(invisible(NULL))
#   }
#   ref_peaks$rel_intensity <- 1000 * ref_peaks$rel_intensity / max(ref_peaks$rel_intensity, na.rm = TRUE)
#
#   # Query peaks (top spectrum) — only if found
#   if (all(is.na(match$Matched_features_mz)) || all(is.na(match$Matched_features_int))) {
#     query_df <- data.frame(mzmed = numeric(0), rel_intensity = numeric(0))
#   } else {
#     query_df <- data.frame(
#       mzmed = match$Matched_features_mz,
#       rel_intensity = match$Matched_features_int
#     )
#     query_df$rel_intensity <- 1000 * query_df$rel_intensity / max(query_df$rel_intensity, na.rm = TRUE)
#   }
#
#   # Create Spectra object
#   spd <- S4Vectors::DataFrame(
#     msLevel = rep(1L, 2),
#     polarity = rep(1L, 2),
#     mz = IRanges::NumericList(query_df$mzmed, ref_peaks$mz),
#     intensity = IRanges::NumericList(query_df$rel_intensity, ref_peaks$rel_intensity)
#   )
#   sps <- Spectra::Spectra(spd)
#
#   # Cosine similarity (may be NA if query is empty)
#   inner_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "inner", ppm = ppm)[1, 2], error = function(e) NA)
#   outer_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "outer", ppm = ppm)[1, 2], error = function(e) NA)
#
#   # Output file path
#   filepath <- NULL
#   if (!is.null(output_dir)) {
#     if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#     filename <- paste0(gsub("[^A-Za-z0-9_]", "_", comp_name_ref), ".png")
#     filepath <- file.path(output_dir, filename)
#   }
#
#   # Plot
#   if (!is.null(filepath)) {
#     png(filepath, width = 8, height = 6, units = "in", res = 300)
#     Spectra::plotSpectraMirror(sps[1], sps[2], main = "", ppm = ppm)
#     title(main = comp_name_ref, cex.main = 1.3)
#     if (nrow(query_df) == 0) {
#       mtext("No matched query peaks", line = 0.5, cex = 0.9, col = "red")
#     } else {
#       mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
#                    " | Inner Cosine: ", round(inner_cosine, 3),
#                    " | Match score: ", round(match$Match_peak_score, 3)),
#             line = 0.5, cex = 0.9)
#     }
#     dev.off()
#   } else {
#     Spectra::plotSpectraMirror(sps[1], sps[2], main = comp_name_ref, ppm = ppm)
#     if (nrow(query_df) == 0) {
#       mtext("No matched query peaks", line = 0.5, cex = 0.9, col = "red")
#     } else {
#       mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
#                    " | Inner Cosine: ", round(inner_cosine, 3),
#                    " | Match score: ", round(match$Match_peak_score, 3)),
#             line = 0.5, cex = 0.9)
#     }
#   }
#
#   invisible(list(
#     compound_name = comp_name_ref,
#     inner_cosine = inner_cosine,
#     outer_cosine = outer_cosine,
#     plot_path = filepath
#   ))
# }
