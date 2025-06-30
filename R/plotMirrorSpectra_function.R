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

#' @rdname plotMirrorSpectra
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




# setMethod("plotMirrorSpectra", signature(object = "DirectLibMatch"),
#           function(object, matches, output_dir = NULL, ppm = 10) {
#             if (!is.null(output_dir) && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#
#             results <- list()
#
#             for (i in seq_along(matches)) {
#               match <- matches[[i]]
#               if (is.null(match$Compound_name)) next
#
#               res <- plotMirrorSpectraSingle(match = match, object = object, output_dir = output_dir, ppm = ppm)
#               if (!is.null(res)) {
#                 results[[match$Compound_name]] <- res
#               }
#             }
#
#             invisible(results)
#           })




#' Plot mirror plot for a single match using Spectra
#'
#' @param match A single element of the output list from `libAnnotate()`
#' @param object The DirectLibMatch object used to produce the match
#' @param output_dir Optional directory to save the plot as PNG
#' @param ppm Mass tolerance (in ppm) for matching peaks in cosine and mirror plots
#' @return A named list with cosine scores and plot path
#' @import Spectra
#' @import S4Vectors
#' @import IRanges
#' @export
plotMirrorSpectraSingle <- function(match, object, output_dir = NULL, ppm = 10) {
  comp_name_ref <- match$Compound_name

  compound_idx <- which(vapply(object@compounds_list, function(x)
    isTRUE(x$Name == comp_name_ref), logical(1)))

  if (length(compound_idx) == 0) {
    warning("Compound not found in compound list: ", comp_name_ref)
    return(invisible(NULL))
  }

  compound <- object@compounds_list[[compound_idx]]

  ref_peaks <- data.frame(
    mz = match$Compound_peak_mz_full,
    rel_intensity = match$Compound_peak_int_full
  )
  ref_peaks$rel_intensity <- 1000 * ref_peaks$rel_intensity / max(ref_peaks$rel_intensity, na.rm = TRUE)

  query_df <- data.frame(
    mzmed = match$Matched_features_mz,
    rel_intensity = match$Matched_features_int
  )
  query_df$rel_intensity <- 1000 * query_df$rel_intensity / max(query_df$rel_intensity, na.rm = TRUE)

  # Create Spectra object
  spd <- S4Vectors::DataFrame(
    msLevel = rep(1L, 2),
    polarity = rep(1L, 2),
    mz = IRanges::NumericList(query_df$mzmed, ref_peaks$mz),
    intensity = IRanges::NumericList(query_df$rel_intensity, ref_peaks$rel_intensity)
  )
  sps <- Spectra::Spectra(spd)

  inner_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "inner", ppm = ppm)[1, 2], error = function(e) NA)
  outer_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "outer", ppm = ppm)[1, 2], error = function(e) NA)

  filepath <- NULL
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- paste0(gsub("[^A-Za-z0-9_]", "_", comp_name_ref), ".png")
    filepath <- file.path(output_dir, filename)
  }

  # Plot
  if (!is.null(filepath)) {
    png(filepath, width = 8, height = 6, units = "in", res = 300)
    Spectra::plotSpectraMirror(sps[1], sps[2], main = "", ppm = ppm)
    title(main = comp_name_ref, cex.main = 1.3)
    mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
                 " | Inner Cosine: ", round(inner_cosine, 3),
                 " | Match score: ", round(match$Match_peak_score, 3)),
          line = 0.5, cex = 0.9)
    dev.off()
  } else {
    Spectra::plotSpectraMirror(sps[1], sps[2], main = comp_name_ref, ppm = ppm)
    mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
                 " | Inner Cosine: ", round(inner_cosine, 3),
                 " | Match score: ", round(match$Match_peak_score, 3)),
          line = 0.5, cex = 0.9)
  }

  invisible(list(
    compound_name = comp_name_ref,
    inner_cosine = inner_cosine,
    outer_cosine = outer_cosine,
    plot_path = filepath
  ))
}





