#' summarize CP metrics to features
#'
#' @description
#' Ultra-simple and fast feature ⇄ CP mapping that assumes chromPeaks rownames
#' are of the form "CP00001", "CP00002", ... and correspond to rows 1, 2, ...
#' It summarizes selected CP-level metrics (e.g., beta_cor, beta_snr) over the
#' peak indices listed in \code{featureDefinitions(x)$peakidx}.
#'
#' @param x XCMSnExp or XcmsExperiment with peaks grouped (after groupChromPeaks()).
#' @param metrics character vector of CP metric column names in chromPeaks(x),
#'   default \code{c("beta_cor","beta_snr")}.
#' @param fun summary function: \code{"mean"} (default) or \code{"median"}.
#' @param append_intensity logical, append \code{featureValues(x, "into")} (default TRUE).
#' @param na_rm logical, pass to mean/median (default TRUE).
#'
#' @details
#' This routine is intentionally minimal. It assumes:
#' - \code{rownames(chromPeaks(x))} are "CP" + zero-padded integers that match
#'   row positions (e.g., CP00037 is row 37).
#' - \code{featureDefinitions(x)$peakidx} are integer indices per feature. If
#'   your peakidx are stored as comma-separated character strings, they will be
#'   parsed into integers.
#'
#'
#' @return A DataFrame of featureDefinitions with columns \code{<fun>_<metric>};
#'   if \code{append_intensity=TRUE}, sample intensity columns are cbind-ed.
#'
#' @export
link_feature_CP_id <- function(
    x,
    metrics = c("beta_cor", "beta_snr"),
    fun = c("mean", "median"),
    append_intensity = TRUE,
    na_rm = TRUE
) {
  stopifnot(inherits(x, "XCMSnExp") || inherits(x, "XcmsExperiment"))
  fun <- match.arg(fun)

  cp <- chromPeaks(x)                 # CP table
  fd <- featureDefinitions(x)         # feature table
  pk <- fd$peakidx                    # indices per feature

  metrics <- intersect(metrics, colnames(cp))
  if (!length(metrics)) stop("None of the requested metrics found in chromPeaks(x).")

  # Coerce peakidx to a plain list of integer vectors
  if (!is.list(pk)) pk <- as.list(pk)
  pk <- lapply(pk, function(v) {
    if (is.null(v)) return(integer())
    if (is.character(v) && length(v) == 1L) {
      v <- trimws(unlist(strsplit(v, ",")))
      suppressWarnings(as.integer(v))
    } else {
      suppressWarnings(as.integer(v))
    }
  })


  # per feature
  summarize_vec <- switch(
    fun,
    mean   = function(z) if (length(z)) mean(z, na.rm = na_rm) else NA_real_,
    median = function(z) if (length(z)) median(z, na.rm = na_rm) else NA_real_
  )

  for (m in metrics) {
    vals <- vapply(pk, function(idx) {
      if (!length(idx)) return(NA_real_)

      idx <- idx[idx >= 1 & idx <= nrow(cp)]
      if (!length(idx)) return(NA_real_)
      summarize_vec(cp[idx, m, drop = TRUE])
    }, numeric(1))
    fd[[paste0(fun, "_", m)]] <- vals
  }

  if (append_intensity) {
    fint <- featureValues(x, value = "into")
    fd <- cbind(fd, fint)
  }

  fd
}


#---------------------------------------------------------------------------------
# test
# match_res <- link_feature_CP_id(xdata_MS1)
#
#
# head(match_res)
#
# head(chromPeaks(xdata_MS1))
# head(featureDefinitions(xdata_MS1))$peakidx[1][[1]]

#----------------------------------------------------------------------------------



#' Link beta metrics (beta_cor, beta_snr) from CPs to features
#'
#' Expands CP-level metrics into feature-level summaries (mean/median) and
#' optionally appends intensity values.
#'
#' @param x XcmsExperiment (must have grouped features).
#' @param metrics CP metric names to use (default: beta_cor, beta_snr).
#' @param fun Summary function ("mean" or "median").
#' @param append_intensity Append feature intensity matrix? (default TRUE).
#' @param na_rm Remove NAs in summary functions.
#' @return Data frame with featureDefinitions plus beta metrics.
link_feature_CP_wide <- function(
    x,
    metrics = c("beta_cor", "beta_snr"),
    fun = c("mean", "median"),
    append_intensity = TRUE,
    na_rm = TRUE
) {
  stopifnot(inherits(x, "XcmsExperiment"))
  fun <- match.arg(fun)

  cp <- chromPeaks(x)
  fd <- featureDefinitions(x)
  pk <- fd$peakidx

  metrics <- intersect(metrics, colnames(cp))
  if (!length(metrics)) stop("No requested metrics in chromPeaks(x).")

  # normalize peakidx into integer vectors
  if (!is.list(pk)) pk <- as.list(pk)
  pk <- lapply(pk, function(v) {
    if (is.character(v) && length(v) == 1L) {
      as.integer(trimws(unlist(strsplit(v, ","))))
    } else as.integer(v)
  })

  # summarizer
  summarize_vec <- switch(fun,
                          mean   = function(z) if (length(z)) mean(z, na.rm = na_rm) else NA_real_,
                          median = function(z) if (length(z)) median(z, na.rm = na_rm) else NA_real_
  )

  for (m in metrics) {
    vals <- vapply(pk, function(idx) {
      idx <- idx[idx >= 1 & idx <= nrow(cp)]
      if (!length(idx)) return(NA_real_)
      summarize_vec(cp[idx, m, drop = TRUE])
    }, numeric(1))
    fd[[paste0(fun, "_", m)]] <- vals
  }

  if (append_intensity) {
    fint <- featureValues(x, value = "into")
    fd <- cbind(fd, fint)
  }

  fd
}



# link_feature_CP_keep_vectors_test <- link_feature_CP_wide(xdata_MS1)


# head(link_feature_CP_keep_vectors_test)






#' Link CPC metrics from CPs to features
#'
#' Expands CPC peak-level metrics into wide feature-level columns.
#'
#' @param x XCMSnExp (with peak detection).
#' @param cpc_obj CPC object from filter_xcms_peaklist().
#' @param metrics Which CPC columns to expand (default: all except housekeeping).
#' @param append_intensity Append feature intensity matrix? (default TRUE).
#' @return Data frame with featureDefinitions plus CPC metrics.
link_feature_CPC_wide <- function(
    x,
    cpc_obj,
    metrics = NULL,
    append_intensity = TRUE
) {
  stopifnot(inherits(x, "XCMSnExp"))
  stopifnot(inherits(cpc_obj, "cpc"))

  cp_cpc <- as.data.frame(suppressWarnings(cpc::cpt(cpc_obj)))
  fd <- featureDefinitions(x)
  pk <- fd$peakidx

  if (is.null(metrics)) {
    drop_cols <- c("id","rt","rtmin","rtmax","file","note","exectime")
    metrics <- setdiff(colnames(cp_cpc), drop_cols)
  }

  # normalize peakidx
  if (!is.list(pk)) pk <- as.list(pk)
  pk <- lapply(pk, function(v) {
    if (is.character(v) && length(v) == 1L) {
      as.integer(trimws(unlist(strsplit(v, ","))))
    } else as.integer(v)
  })

  K <- max(vapply(pk, length, integer(1L)), na.rm = TRUE)
  if (!is.finite(K)) K <- 0L

  for (m in metrics) {
    vecs <- lapply(pk, function(idx) if (length(idx)) cp_cpc[idx, m] else numeric(0))

    W <- matrix(NA_real_, nrow = nrow(fd), ncol = K)
    if (K > 0) {
      for (i in seq_len(nrow(fd))) {
        zi <- vecs[[i]]
        if (length(zi)) W[i, seq_len(min(length(zi), K))] <- zi[seq_len(min(length(zi), K))]
      }
      colnames(W) <- paste0(m, "_", seq_len(K))
      fd <- cbind(fd, W, row.names = rownames(fd))
    }
  }

  if (append_intensity) {
    fint <- featureValues(x, value = "into")
    fd <- cbind(fd, fint)
  }

  fd
}

