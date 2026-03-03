#' Compare full vs truncated spectra within one library
#'
#' @description
#' This function finds compounds present in both the full and truncated versions
#' of a library (by matching the `$Name` field), computes cosine similarity,
#' and generates mirror plots showing before/after truncation differences.
#'
#' @param full_list List of full (untruncated) spectra.
#' @param trunc_list List of truncated spectra (same structure as full_list).
#' @param output_dir Directory to save PNG plots.
#' @param pdf_name Name of the combined PDF file.
#' @param ppm Numeric, tolerance in ppm for m/z matching.
#'
#' @return
#' A list containing cosine scores, number of peaks before/after truncation,
#' and saved file paths.
#'
#' @import Spectra
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @export
plotMirrorFullVsTrunc <- function(full_list,
                                  trunc_list,
                                  output_dir = "mirror_plots_trunc_compare",
                                  pdf_name = "mirror_plots_trunc_compare.pdf",
                                  ppm = 10) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  pdf_file <- file.path(output_dir, pdf_name)
  grDevices::pdf(pdf_file, width = 8, height = 6)

  full_names  <- vapply(full_list, function(x) x$Name, character(1))
  trunc_names <- vapply(trunc_list, function(x) x$Name, character(1))
  common_names <- intersect(full_names, trunc_names)

  results <- list()
  skipped <- character()

  message("Found ", length(common_names), " overlapping compounds to compare.\n")

  for (nm in common_names) {
    message("Processing: ", nm)

    full_idx  <- which(full_names == nm)[1]
    trunc_idx <- which(trunc_names == nm)[1]
    full_entry  <- full_list[[full_idx]]
    trunc_entry <- trunc_list[[trunc_idx]]

    peaks_full  <- full_entry$Peaks
    peaks_trunc <- trunc_entry$Peaks

    # check structure
    if (is.null(peaks_full) || is.null(peaks_trunc) ||
        !"mz" %in% names(peaks_full) || !"mz" %in% names(peaks_trunc) ||
        !"rel_intensity" %in% names(peaks_full) || !"rel_intensity" %in% names(peaks_trunc)) {
      message("Skipping ", nm, " — missing required columns 'mz' or 'rel_intensity'.")
      skipped <- c(skipped, nm)
      next
    }

    # check empty
    if (nrow(peaks_full) == 0 || nrow(peaks_trunc) == 0) {
      message("Skipping ", nm, " — empty peaks.")
      skipped <- c(skipped, nm)
      next
    }

    # ensure numeric
    peaks_full$rel_intensity  <- as.numeric(peaks_full$rel_intensity)
    peaks_trunc$rel_intensity <- as.numeric(peaks_trunc$rel_intensity)

    # handle invalid intensities
    max_full  <- suppressWarnings(max(peaks_full$rel_intensity, na.rm = TRUE))
    max_trunc <- suppressWarnings(max(peaks_trunc$rel_intensity, na.rm = TRUE))
    if (!is.finite(max_full) || !is.finite(max_trunc) || max_full <= 0 || max_trunc <= 0) {
      message("Skipping ", nm, " — invalid or zero intensities.")
      skipped <- c(skipped, nm)
      next
    }

    # normalize to 1000
    peaks_full$rel_intensity  <- 1000 * peaks_full$rel_intensity  / max_full
    peaks_trunc$rel_intensity <- 1000 * peaks_trunc$rel_intensity / max_trunc

    n_full  <- nrow(peaks_full)
    n_trunc <- nrow(peaks_trunc)
    message(" - Full peaks: ", n_full, " | Trunc peaks: ", n_trunc)

    # create Spectra objects
    spd <- S4Vectors::DataFrame(
      msLevel   = rep(1L, 2),
      polarity  = rep(1L, 2),
      mz        = IRanges::NumericList(peaks_full$mz, peaks_trunc$mz),
      intensity = IRanges::NumericList(peaks_full$rel_intensity, peaks_trunc$rel_intensity)
    )
    sps <- Spectra::Spectra(spd)

    inner_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "inner", ppm = ppm)[1, 2],
                             error = function(e) NA)
    outer_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "outer", ppm = ppm)[1, 2],
                             error = function(e) NA)

    safe_name <- gsub("[^A-Za-z0-9_]", "_", nm)
    png_file <- file.path(output_dir, paste0(safe_name, "_mirror.png"))

    # PNG export
    grDevices::png(png_file, width = 8, height = 6, units = "in", res = 300, type = "cairo")
    Spectra::plotSpectraMirror(sps[1], sps[2],
                               main = paste(nm, "- Full vs Truncated"),
                               ppm = ppm)
    mtext(paste0(
      "Outer Cosine: ", round(outer_cosine, 3),
      " | Inner Cosine: ", round(inner_cosine, 3),
      " | Peaks (full/trunc): ", n_full, "/", n_trunc
    ), line = 0.5, cex = 0.9)
    grDevices::dev.off()

    # PDF export
    Spectra::plotSpectraMirror(sps[1], sps[2],
                               main = paste(nm, "- Full vs Truncated"),
                               ppm = ppm)
    mtext(paste0(
      "Outer Cosine: ", round(outer_cosine, 3),
      " | Inner Cosine: ", round(inner_cosine, 3),
      " | Peaks (full/trunc): ", n_full, "/", n_trunc
    ), line = 0.5, cex = 0.9)

    results[[nm]] <- list(
      compound_name = nm,
      inner_cosine  = inner_cosine,
      outer_cosine  = outer_cosine,
      n_full        = n_full,
      n_trunc       = n_trunc,
      png_path      = png_file
    )
  }

  grDevices::dev.off()
  message("✅ Completed mirror plots. Skipped ", length(skipped), " compounds with missing peaks.\n")

  invisible(results)
}


#' Compare spectra of common compounds across two libraries
#'
#' @description This function generates mirror plots for each compound in a matched data frame by comparing two libraries,
#' saves one PNG per compound, and also saves a single combined PDF with all plots.
#'
#' @param matched_df Data frame with compound name matches between two libraries.
#' @param in_house_list List of in-house compounds.
#' @param public_list List of public compounds.
#' @param col_inhouse Column name in matched_df with in-house compound names.
#' @param col_public Column name in matched_df with public compound names.
#' @param output_dir Directory to save plots.
#' @param pdf_name Name of the pdf.
#' @param ppm Numeric, tolerance in ppm for m/z matching.
#'
#' @return A list of results with cosine scores and saved file paths.
#' @import Spectra
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @export
plotLibraryMirrorCompare <- function(matched_df,
                                     in_house_list,
                                     public_list,
                                     col_inhouse,
                                     col_public,
                                     output_dir = "library_compare_mirror_plots",
                                     pdf_name = "library_compare_mirror_plots.pdf",
                                     ppm = 10) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  pdf_file <- file.path(output_dir, pdf_name)
  grDevices::pdf(pdf_file, width = 8, height = 6)

  results <- list()
  counter_env <- new.env(parent = emptyenv())  # local counter

  for (i in seq_len(nrow(matched_df))) {
    match <- matched_df[i, , drop = FALSE]

    in_name  <- match[[col_inhouse]]
    pub_name <- match[[col_public]]

    in_idx <- which(vapply(in_house_list, function(x) x$Name == in_name, logical(1)))[1]
    pub_idx <- which(vapply(public_list, function(x) x$Name == pub_name, logical(1)))[1]
    if (is.na(in_idx) || is.na(pub_idx)) next

    in_entry  <- in_house_list[[in_idx]]
    pub_entry <- public_list[[pub_idx]]
    peaks1 <- in_entry$Peaks
    peaks2 <- pub_entry$Peaks
    if (nrow(peaks1) == 0 || nrow(peaks2) == 0) next

    # normalize intensities
    peaks1$rel_intensity <- 1000 * peaks1$rel_intensity / max(peaks1$rel_intensity, na.rm = TRUE)
    peaks2$rel_intensity <- 1000 * peaks2$rel_intensity / max(peaks2$rel_intensity, na.rm = TRUE)

    spd <- S4Vectors::DataFrame(
      msLevel   = rep(1L, 2),
      polarity  = rep(1L, 2),
      mz        = IRanges::NumericList(peaks1$mz, peaks2$mz),
      intensity = IRanges::NumericList(peaks1$rel_intensity, peaks2$rel_intensity)
    )
    sps <- Spectra::Spectra(spd)

    inner_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "inner", ppm = ppm)[1, 2], error = function(e) NA)
    outer_cosine <- tryCatch(Spectra::compareSpectra(sps, type = "outer", ppm = ppm)[1, 2], error = function(e) NA)

    # safe filename
    filename <- make_safe_filename(in_name, counter_env = counter_env)
    png_file <- file.path(output_dir, filename)

    # PNG export
    grDevices::png(png_file, width = 8, height = 6, units = "in", res = 300, type = "cairo")
    Spectra::plotSpectraMirror(sps[1], sps[2],
                               main = paste(in_name, " vs ", pub_name),
                               ppm = ppm)
    mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
                 " | Inner Cosine: ", round(inner_cosine, 3)),
          line = 0.5, cex = 0.9)
    grDevices::dev.off()

    # PDF export
    Spectra::plotSpectraMirror(sps[1], sps[2],
                               main = paste(in_name, " vs ", pub_name),
                               ppm = ppm)
    mtext(paste0("Outer Cosine: ", round(outer_cosine, 3),
                 " | Inner Cosine: ", round(inner_cosine, 3)),
          line = 0.5, cex = 0.9)

    results[[i]] <- list(
      compound_inhouse = in_name,
      compound_public  = pub_name,
      inner_cosine     = inner_cosine,
      outer_cosine     = outer_cosine,
      png_path         = png_file
    )
  }

  grDevices::dev.off()
  invisible(results)
}


#' Handles 1:many and many:many mappings in matched_df
#'
#' @param matched_df Data.frame with at least two columns specifying matched compound names
#' @param in_house_list List of in-house library entries, each with $Peaks
#' @param public_list List of public library entries, each with $Peaks
#' @param col_inhouse Column name in matched_df for in-house compound names
#' @param col_public Column name in matched_df for public compound names
#' @param ppm Numeric, tolerance in ppm for m/z matching (default: 10)
#'
#' @return A list of results per compound pair
#' @export
buildVennDataList <- function(matched_df,
                              in_house_list,
                              public_list,
                              col_inhouse = "Name_in_house",
                              col_public  = "Name_public",
                              ppm = 10) {

  make_ids <- function(prefix, n) if (n > 0) paste0(prefix, seq_len(n)) else character(0)

  results <- list()

  for (i in seq_len(nrow(matched_df))) {
    in_name  <- matched_df[[col_inhouse]][i]
    pub_name <- matched_df[[col_public]][i]

    in_idx  <- which(vapply(in_house_list, function(x) x$Name == in_name, logical(1)))
    pub_idx <- which(vapply(public_list,   function(x) x$Name == pub_name, logical(1)))

    if (!length(in_idx) || !length(pub_idx)) {
      warning("Compound not found: ", in_name, " or ", pub_name)
      next
    }

    # loop over *all* matches (1:many or many:many)
    for (ii in in_idx) {
      for (jj in pub_idx) {
        peaks1 <- as.matrix(in_house_list[[ii]]$Peaks[, c("mz","rel_intensity")])
        colnames(peaks1) <- c("mz","intensity")
        peaks2 <- as.matrix(public_list[[jj]]$Peaks[, c("mz","rel_intensity")])
        colnames(peaks2) <- c("mz","intensity")

        inner <- Spectra::joinPeaks(peaks1, peaks2, ppm = ppm, type = "inner")

        both      <- nrow(inner$x)
        in_total  <- nrow(peaks1)
        pub_total <- nrow(peaks2)

        in_only  <- in_total - both
        pub_only <- pub_total - both

        venn_counts <- c(in_only = in_only, pub_only = pub_only, both = both)

        # Build unique IDs. Overlap elements share the SAME labels ("B1..B{both}") in both sets.
        ih_ids <- c(make_ids("IH", in_only), make_ids("B", both))
        pu_ids <- c(make_ids("PU", pub_only), make_ids("B", both))

        # (Optional) sanity checks to catch bugs early
        stopifnot(length(intersect(ih_ids, pu_ids)) == both)
        stopifnot(length(setdiff(ih_ids, pu_ids)) == in_only)
        stopifnot(length(setdiff(pu_ids, ih_ids)) == pub_only)

        venn_data <- list(InHouse = ih_ids, Public = pu_ids)

        results[[length(results) + 1]] <- list(
          in_house_name = in_house_list[[ii]]$Name,
          public_name   = public_list[[jj]]$Name,
          venn_counts   = venn_counts,
          venn_data     = venn_data
        )
      }
    }
  }

  results
}



# helper: sanitize + truncate
#' @description Generate safe filenames from in-house compound names
#' @param in_name character, the in-house compound name
#' @param counter_env environment for keeping counters (so duplicates get suffixes)
#' @param max_len integer, maximum length for base name (default 30)
#' @param max_filename_len integer, maximum total filename length (default 100)
#' @return Safe filename ending in ".png"
#' @noRd
make_safe_filename <- function(in_name,
                               counter_env = .GlobalEnv,
                               max_len = 30,
                               max_filename_len = 100) {


  short_in <- gsub("[^A-Za-z0-9_]", "_", trimws(in_name))
  if (nchar(short_in) > max_len) {
    short_in <- substr(short_in, 1, max_len)
  }

  if (!exists("name_counter", envir = counter_env, inherits = FALSE)) {
    assign("name_counter", list(), envir = counter_env)
  }
  counter <- get("name_counter", envir = counter_env)

  if (!short_in %in% names(counter)) {
    counter[[short_in]] <- 1
  } else {
    counter[[short_in]] <- counter[[short_in]] + 1
  }
  assign("name_counter", counter, envir = counter_env)

  # suffix if duplicated
  suffix <- if (counter[[short_in]] > 1) paste0("_", counter[[short_in]]) else ""

  filename <- paste0(short_in, suffix, ".png")

  # final cutoff
  if (nchar(filename) > max_filename_len) {
    filename <- substr(filename, 1, max_filename_len - 4)
    filename <- paste0(filename, ".png")
  }

  filename
}



#' Export Venn diagrams to PNGs and a combined PDF
#'
#' @param venn_data_list list from buildVennDataList()
#' @param output_dir directory for saving plots
#' @param pdf_file filename for combined PDF (all compounds)
#' @return Invisibly returns NULL
#' @export
exportVennPlots <- function(venn_data_list,
                            output_dir = "venn_plots",
                            pdf_file = "venn_all.pdf") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # local environment for filename counters
  counter_env <- new.env(parent = emptyenv())

  grDevices::pdf(file.path(output_dir, pdf_file), width = 8, height = 6)

  for (res in venn_data_list) {
    if (is.null(res) || is.null(res$venn_data)) next

    venn_data <- res$venn_data
    in_name   <- res$in_house_name
    pub_name  <- res$public_name

    venn_list <- list(
      InHouse = venn_data$InHouse,
      Public  = venn_data$Public
    )

    p <- ggvenn(
      venn_list,
      show_percentage = TRUE,
      fill_color = c("#A6CEE3", "#FFFF99"),
      fill_alpha = 0.6
    ) + ggtitle(paste0(in_name, " vs ", pub_name))

    print(p) # to combined PDF

    filename <- make_safe_filename(in_name, counter_env = counter_env)
    png_file <- file.path(output_dir, filename)

    grDevices::png(png_file, width = 8, height = 6, units = "in", res = 300, type = "cairo")
    print(p)
    grDevices::dev.off()
  }

  grDevices::dev.off()
  message("Saved all Venn plots (PNG + combined PDF) to ", output_dir)
  invisible(NULL)
}

