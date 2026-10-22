# ============================================================
# helper_read_ms2.R
#
# Converts MS2 consensus data into a named list keyed by
# feature_name, ready for run_isf_annotation().
#
# Accepts three input types:
#   1 Spectra object BEFORE empty removal (ms2_consensus)
#   2 Spectra object AFTER empty removal  (ms2_final)
#   3 path to an MGF file (.mgf)
#
# When feature_table is provided, the MS2 scan RT (DDA trigger time)
# is replaced by the MS1 peak apex RT (rtmed) from the feature table.
# This unifies MS2 and MS1 RT at source, avoiding discrepancies caused
# by DDA. (up to ~20s).
# ============================================================

#' Get MS2 peaks with feature_name assigned
#'
#' @description
#' Wrapper around \code{chromPeakSpectra()} that maps
#' \code{chrom_peak_id} to \code{feature_name} and adds it as a
#' spectra variable.
#'
#' @param x_exp A processed \code{XcmsExperiment} object.
#' @param expand_rt Numeric. RT expansion in seconds.  Default \code{10}.
#' @param expand_mz Numeric. m/z expansion in Da.  Default \code{0.02}.
#' @param ... Additional arguments passed to \code{chromPeakSpectra()}.
#'
#' @return A \code{Spectra} object with \code{feature_name} added.
#' @export
get_ms2_peaks <- function(x_exp,
                          expand_rt = 10,
                          expand_mz = 0.02,
                          msms_files = NULL,
                          ...) {
  # ── Step 1: Peak-based assignment via chromPeakSpectra() ───────────────────
  # Assigns MS2 scans to features using xcms-detected chromatographic peaks
  ms2_peaks <- chromPeakSpectra(
    x_exp,
    msLevel     = 2L,
    method      = "all",
    return.type = "Spectra",
    skipFilled  = FALSE,
    expandRt    = expand_rt,
    expandMz    = expand_mz,
    ...
  )

  feat_def <- featureDefinitions(x_exp)
  cp_names <- rownames(chromPeaks(x_exp))

  p2f <- data.frame(
    chrom_peak_id = unlist(lapply(seq_len(nrow(feat_def)),
                                  function(i) cp_names[feat_def$peakidx[[i]]])),
    feature_name  = rep(rownames(feat_def), lengths(feat_def$peakidx)),
    stringsAsFactors = FALSE
  )

  ms2_peaks$feature_name <- p2f$feature_name[
    match(ms2_peaks$chrom_peak_id, p2f$chrom_peak_id)
  ]

  n_assigned_peak <- sum(!is.na(ms2_peaks$feature_name))
  message("[get_ms2_peaks] Peak-based assignments:     ", n_assigned_peak)

  # ── feature lookup ────────────────────────
  fd        <- as.data.frame(feat_def)
  fd$feature_name <- rownames(feat_def)
  fd_mz    <- if ("mzmed" %in% names(fd)) fd$mzmed else (fd$mzmin + fd$mzmax) / 2
  fd_rt    <- if ("rtmed" %in% names(fd)) fd$rtmed else (fd$rtmin + fd$rtmax) / 2
  fd_mzmin <- fd$mzmin - expand_mz
  fd_mzmax <- fd$mzmax + expand_mz
  fd_rtmin <- fd$rtmin - expand_rt
  fd_rtmax <- fd$rtmax + expand_rt
  fd_fn    <- fd$feature_name

  # ── Internal helper: proximity-assign a Spectra object ────────────────────
  .proximity_assign <- function(sp, label = "scans") {
    pmz_vec <- sp$precursorMz
    prt_vec <- sp$rtime
    n       <- length(sp)
    fn_out  <- rep(NA_character_, n)

    for (k in seq_len(n)) {
      pmz <- pmz_vec[k]; prt <- prt_vec[k]
      if (is.na(pmz) || is.na(prt)) next
      cand <- which(pmz >= fd_mzmin & pmz <= fd_mzmax &
                    prt >= fd_rtmin & prt <= fd_rtmax)
      if (length(cand) == 0L) next
      fn_out[k] <- fd_fn[cand[which.min(abs(pmz - fd_mz[cand]))]]
    }

    fn_assigned <- sum(!is.na(fn_out))
    message("[get_ms2_peaks] Proximity matches (", label, "): ", fn_assigned,
            " of ", n, " scans")
    fn_out
  }

  # ── Step 2: fallback for unassigned scans already in ms2_peaks ────────────

  na_mask <- is.na(ms2_peaks$feature_name)
  if (any(na_mask)) {
    sp_na    <- ms2_peaks[na_mask]
    fn_fb    <- .proximity_assign(sp_na, label = "unassigned in ms2_peaks")
    fn_vec   <- ms2_peaks$feature_name
    fn_vec[na_mask] <- fn_fb
    ms2_peaks$feature_name <- fn_vec
    n_fallback <- sum(!is.na(fn_fb))
    message("[get_ms2_peaks] Fallback (unassigned ms2_peaks): ", n_fallback)
  }

  # ── Step 3: Load MSMS files ─────────────────

  if (is.null(msms_files)) {
    # file in x_exp with basename"MSMS" or "msms"
    all_files  <- fileNames(x_exp)
    msms_files <- all_files[grepl("MSMS|msms", basename(all_files))]
  }

  if (length(msms_files) > 0L) {
    message("[get_ms2_peaks] Loading MS2 scans directly from ",
            length(msms_files), " MSMS file(s) ...")

    ms2_direct        <- Spectra::Spectra(msms_files)
    ms2_direct        <- ms2_direct[ms2_direct$msLevel == 2L]

    ms2_direct        <- ms2_direct[lengths(Spectra::mz(ms2_direct)) > 0L]

    fn_direct         <- .proximity_assign(ms2_direct, label = "MSMS direct")
    ms2_direct$feature_name <- fn_direct

    already_covered   <- unique(ms2_peaks$feature_name[
                           !is.na(ms2_peaks$feature_name)])
    keep_mask         <- !is.na(fn_direct) & !fn_direct %in% already_covered

    n_new <- sum(keep_mask)
    message("[get_ms2_peaks] New scans from MSMS: ", n_new,
            " (covering ",
            length(unique(fn_direct[keep_mask])), " additional features)")

    if (n_new > 0L) {
      ms2_new <- ms2_direct[keep_mask]

      sv_peaks <- Spectra::spectraVariables(ms2_peaks)
      sv_new   <- Spectra::spectraVariables(ms2_new)
      for (sv in setdiff(sv_peaks, sv_new)) {
        ex_val <- ms2_peaks[[sv]][1L]
        na_val <- if (is.character(ex_val)) NA_character_ else
                  if (is.numeric(ex_val))   NA_real_      else
                  if (is.integer(ex_val))   NA_integer_   else
                  if (is.logical(ex_val))   NA             else NA
        ms2_new[[sv]] <- rep(na_val, length(ms2_new))
      }
      ms2_peaks <- c(ms2_peaks, ms2_new)
    }
  } else {
    message("[get_ms2_peaks] No MSMS files detected for direct loading.")
  }

  message("[get_ms2_peaks] MS2 scans total:  ", length(ms2_peaks))
  message("[get_ms2_peaks] Unique features:  ",
          length(unique(ms2_peaks$feature_name[
            !is.na(ms2_peaks$feature_name)])))
  ms2_peaks
}


#' Read MS2 consensus data into a named list
#'
#' @description
#' Converts MS2 consensus data into a named list keyed by
#' \code{feature_name}, ready for \code{\link{run_isf_annotation}()}.
#'
#' Three input types accepted:
#' \describe{
#'   \item{Option A}{\code{Spectra} from \code{combineSpectra()} —
#'     empty spectra removed automatically.}
#'   \item{Option B}{\code{Spectra} already filtered
#'     (\code{ms2_final <- ms2_consensus[lengths(mz(ms2_consensus)) > 0]}).}
#'   \item{Option C}{Path to an MGF file with \code{feature_name=} headers.}
#' }
#'
#' When \code{feature_table} is provided, the MS2 scan RT (DDA
#' ) is replaced by the MS1 peak apex RT (\code{rtmed}) from the
#' feature table.  This unifies MS2 and MS1 RT at source, preventing
#' spurious RT differences of up to ~20 s caused by DDA firing on the
#' rising edge of a chromatographic peak.
#'
#' @param x \code{Spectra} object or MGF file path.
#' @param feature_table \code{data.frame}. Optional MS1 feature table with
#'   columns \code{feature_name} and \code{mzmed}/\code{mz} and
#'   \code{rtmed}/\code{rt}.  When provided:
#'   \itemize{
#'     \item MS2 RT is overridden with MS1 apex RT.
#'     \item Precursor m/z is cross-checked within \code{ppm_check}.
#'   }
#' @param ppm_check Numeric. ppm tolerance for precursor m/z check.
#'   Default \code{20}.
#' @param verbose Logical. Print summary.  Default \code{TRUE}.
#'
#' @return Named list keyed by \code{feature_name}.  Each element:
#'   \describe{
#'     \item{\code{feature_name}}{Character.}
#'     \item{\code{precursor_mz}}{Numeric.}
#'     \item{\code{rt}}{Numeric. MS1 apex RT (seconds) if
#'       \code{feature_table} provided; otherwise MS2 scan RT.}
#'     \item{\code{fragment_mz}}{Numeric vector.}
#'     \item{\code{fragment_int}}{Numeric vector.}
#'     \item{\code{n_fragments}}{Integer.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Option A — ms2_consensus directly (empty spectra removed inside)
#' ms2_list <- read_ms2_consensus(ms2_consensus, feature_table = pool_ft)
#'
#' # Option B — pre-filtered ms2_final
#' ms2_final <- ms2_consensus[lengths(mz(ms2_consensus)) > 0]
#' ms2_list  <- read_ms2_consensus(ms2_final, feature_table = pool_ft)
#'
#' # Option C — MGF file
#' ms2_list  <- read_ms2_consensus("xcms_ms2_consensus_final.mgf",
#'                                  feature_table = pool_ft)
#' }
#' @export
read_ms2_consensus <- function(x,
                               feature_table = NULL,
                               ppm_check     = 20,
                               verbose       = TRUE) {

  # ── Build RT lookup from feature_table ──────────────────────────────────
  ft_rt  <- NULL
  ft_mz  <- NULL
  if (!is.null(feature_table)) {
    ft <- feature_table
    if (!"feature_name" %in% names(ft))
      stop("feature_table must have column 'feature_name'.")

    mz_col <- if ("mzmed" %in% names(ft)) "mzmed" else
      if ("mz"    %in% names(ft)) "mz"    else
        stop("feature_table must have column 'mzmed' or 'mz'.")
    rt_col <- if ("rtmed" %in% names(ft)) "rtmed" else
      if ("rt"    %in% names(ft)) "rt"    else NULL

    ft_mz <- setNames(ft[[mz_col]], ft$feature_name)
    ft_rt <- if (!is.null(rt_col))
      setNames(ft[[rt_col]], ft$feature_name)
    else NULL

    if (!is.null(ft_rt) && verbose)
      message("[read_ms2_consensus] MS2 RT will be unified with MS1 ",
              "apex RT from feature_table.")
  }

  # ── Route to correct parser ──────────────────────────────────────────────
  if (inherits(x, "Spectra")) {

    # Remove empty spectra
    n_before <- length(x)
    x        <- x[lengths(Spectra::mz(x)) > 0]
    n_empty  <- n_before - length(x)
    if (verbose && n_empty > 0L)
      message("[read_ms2_consensus] Removed ", n_empty,
              " empty spectra (peaksCount = 0).")

    ms2_list <- .spectra_to_ms2_list(x, ft_rt = ft_rt)

  } else if (is.character(x) && length(x) == 1L) {

    if (!file.exists(x))
      stop("[read_ms2_consensus] File not found: ", x)
    ms2_list <- .mgf_to_ms2_list(x, ft_rt = ft_rt)

  } else {
    stop("[read_ms2_consensus] `x` must be a Spectra object ",
         "(Options A/B) or a path to an MGF file (Option C).")
  }

  # ── Summary ──────────────────────────────────────────────────────────────
  if (verbose) {
    message("[read_ms2_consensus] Loaded ", length(ms2_list),
            " MS2 entries (unique feature_name).")
    n_frags <- vapply(ms2_list, function(e) e$n_fragments, integer(1L))
    message("  Fragment count: median = ", round(median(n_frags), 1),
            ", range = [", min(n_frags), ", ", max(n_frags), "]")
  }

  # ── precursor m/z check ───────────────────────────────────
  if (!is.null(ft_mz)) {
    common <- intersect(names(ms2_list), names(ft_mz))
    if (length(common) == 0L)
      warning("[read_ms2_consensus] No feature_name overlap between MS2 ",
              "list and feature_table — check naming consistency.")

    mismatches <- vapply(common, function(fn) {
      prec <- ms2_list[[fn]]$precursor_mz
      feat <- ft_mz[[fn]]
      if (is.na(prec) || is.na(feat)) return(FALSE)
      abs(prec - feat) / feat * 1e6 > ppm_check
    }, logical(1L))

    n_mm <- sum(mismatches)
    if (n_mm > 0L)
      warning("[read_ms2_consensus] ", n_mm,
              " feature(s) have precursor m/z deviating >", ppm_check,
              " ppm from feature_table: ",
              paste(common[mismatches], collapse = ", "))
    else if (verbose)
      message("  Precursor m/z check passed for ", length(common),
              " features (tolerance: ", ppm_check, " ppm).")
  }

  ms2_list
}


# ── Internal: Spectra to list ──────────────────────────────────────────────────
.spectra_to_ms2_list <- function(sp, ft_rt = NULL) {
  if (!"feature_name" %in% Spectra::spectraVariables(sp))
    stop("Spectra object must have a 'feature_name' spectra variable. ",
         "Use get_ms2_peaks() to add it before calling combineSpectra().")

  fnames  <- sp$feature_name
  prec_mz <- sp$precursorMz
  rts     <- sp$rtime
  peaks   <- Spectra::peaksData(sp)

  entries <- mapply(function(fn, pmz, rt, pk) {
    # Override MS2 scan RT with MS1 apex RT from feature_table
    final_rt <- if (!is.null(ft_rt) &&
                    fn %in% names(ft_rt) &&
                    !is.na(ft_rt[[fn]]))
      ft_rt[[fn]]
    else rt

    list(
      feature_name = fn,
      precursor_mz = pmz,
      rt           = final_rt,
      fragment_mz  = if (nrow(pk) > 0L) pk[, "mz"]       else numeric(0L),
      fragment_int = if (nrow(pk) > 0L) pk[, "intensity"] else numeric(0L),
      n_fragments  = nrow(pk)
    )
  }, fnames, prec_mz, rts, peaks,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)

  setNames(entries, fnames)
}


# ── Internal: MGF file to list ─────────────────────────────────────────────────
.mgf_to_ms2_list <- function(path, ft_rt = NULL) {
  lines    <- readLines(path, warn = FALSE)
  entries  <- list()
  in_block <- FALSE
  current_block <- list(feature_name = NA_character_,
              precursor_mz = NA_real_,
              rt           = NA_real_,
              frag_mz      = numeric(0L),
              frag_int     = numeric(0L))

  for (ln in lines) {
    ln <- trimws(ln)

    if (ln == "BEGIN IONS") {
      in_block <- TRUE
      current_block <- list(feature_name = NA_character_,
                  precursor_mz = NA_real_,
                  rt           = NA_real_,
                  frag_mz      = numeric(0L),
                  frag_int     = numeric(0L))
      next
    }

    if (ln == "END IONS") {
      if (!is.na(current_block$feature_name) && length(current_block$frag_mz) > 0L) {
        fn <- current_block$feature_name

        # Override RT with MS1 apex RT if available
        final_rt <- if (!is.null(ft_rt) &&
                        fn %in% names(ft_rt) &&
                        !is.na(ft_rt[[fn]]))
          ft_rt[[fn]]
        else current_block$rt

        entries[[fn]] <- list(
          feature_name = fn,
          precursor_mz = current_block$precursor_mz,
          rt           = final_rt,
          fragment_mz  = current_block$frag_mz,
          fragment_int = current_block$frag_int,
          n_fragments  = length(current_block$frag_mz)
        )
      }
      in_block <- FALSE
      next
    }

    if (!in_block) next

    if (grepl("^feature_name=", ln, ignore.case = TRUE)) {
      current_block$feature_name <- sub("^feature_name=", "", ln, ignore.case = TRUE)
    } else if (grepl("^feature_id=", ln, ignore.case = TRUE)) {
      current_block$feature_name <- sub("^feature_id=", "", ln, ignore.case = TRUE)
    } else if (grepl("^PEPMASS=", ln)) {
      current_block$precursor_mz <- as.numeric(
        sub("PEPMASS=", "", strsplit(ln, "[ \t]+")[[1L]][1L]))
    } else if (grepl("^RTINSECONDS=", ln)) {
      current_block$rt <- as.numeric(sub("RTINSECONDS=", "", ln))
    } else if (grepl("^[0-9]", ln)) {
      parts <- strsplit(ln, "[ \t]+")[[1L]]
      if (length(parts) >= 2L) {
        current_block$frag_mz  <- c(current_block$frag_mz,  as.numeric(parts[1L]))
        current_block$frag_int <- c(current_block$frag_int, as.numeric(parts[2L]))
      }
    }
  }

  entries
}
