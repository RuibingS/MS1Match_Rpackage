#' Parse a single compound from a MoNA MSP file
#'
#' @description This function is to parses a NIST format library.
#' It can filter the compound based on the spectrum type, ion mode, instrument_type, Instrument
#' @param compound_lines A character vector of lines for one compound.
#' @param spectrum_type "MS1" or "MS2" to filter by spectrum type.
#' @param ionPolarity "pos" or "neg" to filter polarity.
#' @param ioniMethod Ionization method (optional).
#' @param colEnergy Collision energy (optional).
#' @param instName Instrument name (optional).
#' @import MetaboCoreUtils
#' @return A compound list or NULL if filtered out.
#' @export
parse_single_compound_NIST <- function(compound_lines,
                                       spectrum_type = NULL,
                                       ionPolarity = NULL,
                                       ioniMethod = NULL,
                                       colEnergy = NULL,
                                       instName = NULL) {
  tryCatch({
    compound <- list()

    if (!any(grepl("^Num Peaks: ", compound_lines))) return(NULL)
    if (any(grepl("UVName", compound_lines))) return(NULL)

    is_ms2 <- any(grepl("^MSMS:", compound_lines))
    compound[["spectrum_type"]] <- ifelse(is_ms2, "MS2", "MS1")

    # Early filtering by MS level
    if (!is.null(spectrum_type)) {
      if ((spectrum_type == "MS1" && is_ms2) || (spectrum_type == "MS2" && !is_ms2)) {
        return(NULL)
      }
    }

    for (line in compound_lines) {
      if (grepl(": ", line)) {
        infor_value <- strsplit(line, ": ", fixed = TRUE)[[1]]
        infor <- trimws(infor_value[1])
        value <- trimws(infor_value[2])
        compound[[infor]] <- value
      }

      if (!is.null(ionPolarity) && grepl("^IonPolarity: ", line)) {
        value <- sub("^IonPolarity: ", "", line)
        if (tolower(value) != tolower(ionPolarity)) return(NULL)
      }

      if (!is.null(ioniMethod) && grepl("^IoniMethod: ", line)) {
        value <- sub("^IoniMethod: ", "", line)
        if (value != ioniMethod) return(NULL)
      }

      if (!is.null(colEnergy) && grepl("^ColEnergy: ", line)) {
        value <- sub("^ColEnergy: ", "", line)
        if (value != colEnergy) return(NULL)
      }

      if (!is.null(instName) && grepl("^InstName: ", line)) {
        value <- sub("^InstName: ", "", line)
        if (value != instName) return(NULL)
      }

      if (grepl("^Num Peaks: ", line)) {
        peaks_start <- which(compound_lines == line) + 1
        peaks_line <- paste(compound_lines[peaks_start:length(compound_lines)], collapse = " ")
        peaks_values <- as.numeric(unlist(strsplit(peaks_line, "\\s+")))
        peaks_values <- peaks_values[!is.na(peaks_values)]

        if (length(peaks_values) %% 2 != 0) stop("Uneven number of m/z and intensity values.")

        compound[["Peaks"]] <- data.frame(
          mz = peaks_values[seq(1, length(peaks_values), 2)],
          rel_intensity = peaks_values[seq(2, length(peaks_values), 2)]
        )
      }

      if (grepl("^Formula: ", line)) {
        formula <- trimws(strsplit(line, ": ", fixed = TRUE)[[1]][2])

        # Compute exact mass using MetaboCoreUtils
        exactMass_value <- tryCatch({
          MetaboCoreUtils::calculateMass(formula)[[1]]
        }, error = function(e) {
          message("Error in calculateMass: ", conditionMessage(e))
          NA_real_
        })

        compound[["ExactMass"]] <- exactMass_value
      }
    }

    return(compound)

  }, error = function(e) {
    message("Failed to parse compound: ", conditionMessage(e))
    return(NULL)
  })
}



#' Parses a NIST-format library file using parallel processing.
#'
#' @param file_path Path to the library file.
#' @param spectrum_type Usually "MS1" or "MS2". Set to NULL to disable filtering.
#' @param ionPolarity Either "pos" or "neg". Set to NULL to disable filtering.
#' @param ioniMethod Ionization method (optional).
#' @param colEnergy Collision energy (optional).
#' @param instName Instrument name (optional).
#'
#' @import parallel
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @return A list of parsed and filtered compound spectra.
#' @examples
#' parse_library_parallel("library/Si16.library", ionPolarity = "pos", spectrum_type = "MS1")
#' @export
parse_library_parallel <- function(file_path,
                                   spectrum_type = NULL,
                                   ionPolarity = NULL,
                                   ioniMethod = NULL,
                                   colEnergy = NULL,
                                   instName = NULL) {

  if (!file.exists(file_path)) {
    stop("Error: The file does not exist. Please check the file path.")
  }

  lines <- readLines(file_path, warn = FALSE)
  name_positions <- grep("^Name: ", lines)
  num_compounds <- length(name_positions)

  required_packages <- c("MetaboCoreUtils", "doParallel", "enviPat", "stringr")
  n_cores <- max(1, parallel::detectCores() - 1)

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  chunk_size <- ceiling(num_compounds / n_cores)
  chunk_indices <- split(seq_along(name_positions), ceiling(seq_along(name_positions) / chunk_size))

  # Ensure helper function is available on all nodes
  parallel::clusterExport(cl, varlist = c("parse_single_compound_NIST"), envir = environment())

  # Process chunks in parallel
  compounds <- foreach::foreach(i = chunk_indices, .combine = c, .packages = required_packages) %dopar% {
    res <- list()
    for (idx in i) {
      start <- name_positions[idx]
      end <- if (idx < length(name_positions)) name_positions[idx + 1] - 1 else length(lines)

      compound_lines <- lines[start:end]
      compound <- tryCatch(
        parse_single_compound_NIST(compound_lines,
                                   spectrum_type = spectrum_type,
                                   ionPolarity = ionPolarity,
                                   ioniMethod = ioniMethod,
                                   colEnergy = colEnergy,
                                   instName = instName),
        error = function(e) {
          message("Error parsing compound at index ", idx, ": ", conditionMessage(e))
          return(NULL)
        }
      )

      if (!is.null(compound)) res <- append(res, list(compound))
    }
    res
  }

  parallel::stopCluster(cl)

  return(Filter(Negate(is.null), compounds))
}



















