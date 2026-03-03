#' @rdname plot_single_feature_EIC
#' @export
setMethod("plot_single_feature_EIC",
          signature(object = "EICData"),
          function(object, i, j_vec = NULL) {
            # basic checks
            if (length(i) != 1L || !is.numeric(i) || is.na(i)) {
              stop("'i' must be a single, valid row (feature) index.")
            }
            n_rows <- nrow(object@eic_data)
            n_cols <- ncol(object@eic_data)
            if (i < 1L || i > n_rows) stop("Feature index 'i' is out of bounds.")

            # choose samples
            if (is.null(j_vec)) {
              j_vec <- seq_len(n_cols)
            } else {
              j_vec <- as.integer(j_vec)
              j_vec <- j_vec[!is.na(j_vec)]
              if (!length(j_vec)) stop("No valid sample indices in 'j_vec'.")
              if (any(j_vec < 1L | j_vec > n_cols)) {
                stop("Some 'j_vec' indices are out of bounds [1, ", n_cols, "].")
              }
            }

            # sample labels: use column names if present, else "Sample_<j>"
            coln <- colnames(object@eic_data)
            lab_for <- function(j) if (!is.null(coln)) coln[j] else paste0("Sample_", j)

            # build one data.frame per sample (avoid library/purrr inside package code)
            df_list <- vector("list", length(j_vec))
            for (k in seq_along(j_vec)) {
              j <- j_vec[k]
              chr <- object@eic_data[i, j]
              # chr is an XChromatogram; slots rtime/intensity are numeric vectors
              rt  <- chr@rtime
              inx <- chr@intensity
              keep <- !is.na(rt) & !is.na(inx)
              df_list[[k]] <- data.frame(
                rt = rt[keep],
                intensity = inx[keep],
                sample = lab_for(j),
                stringsAsFactors = FALSE
              )
            }
            df <- do.call(rbind, df_list)
            if (!nrow(df)) stop("No data to plot for the requested feature/samples.")

            # plot
            p <- ggplot2::ggplot(df, ggplot2::aes(x = rt, y = intensity, color = sample)) +
              ggplot2::geom_line() +
              ggplot2::geom_point(size = 0.8) +
              ggplot2::labs(
                title = paste0("EIC for Feature ", i,
                               " across samples: ", paste(j_vec, collapse = ", ")),
                x = "Retention time",
                y = "Intensity",
                color = "Sample"
              ) +
              ggplot2::theme_minimal()

            p
          })
