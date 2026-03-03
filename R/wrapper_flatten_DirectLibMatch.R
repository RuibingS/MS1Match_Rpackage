#' Flatten DirectLibMatch results to a data.frame
#'
#' Convert the `results` slot (a list of per-feature/per-hit lists) from a
#' [DirectLibMatch-class] object into a single `data.frame`. Optionally collapse
#' vector fields `Compound_peak_mz_full` and `Compound_peak_int_full` into
#' comma-separated strings, either overwriting the originals or writing to
#' new `*_str` columns.
#'
#' @param object A [DirectLibMatch-class] object with `results` populated
#'   (e.g., after calling [libAnnotate()]).
#' @param collapse_full Logical; collapse `*_full` vectors to strings? Default is `TRUE`.
#' @param remove_full Logical; if `TRUE`, remove the original `*_full` columns
#'   and write collapsed strings to new columns with `suffix`. If `FALSE`,
#'   overwrite the original `*_full` with the collapsed strings. Default is `FALSE`.
#' @param suffix Character; suffix for new columns when `remove_full = TRUE`.
#'   Default is `"_str"`.
#'
#' @return A `data.frame` with one row per element of `object@results`.
#'
#' @examples
#' \dontrun{
#' obj <- DirectLibMatch(
#'   FT = result_PA14_FT_filtered,
#'   compounds_list = PA14_CBIO_library_MS1,
#'   diff_mz_thr = 0.005, diff_rt_thr = 60, ppm = 10,
#'   ms_level = "MS1", polarity = "pos", results = list()
#' )
#' obj <- libAnnotate(obj)
#' df  <- results_df(obj)
#' df2 <- results_df(obj, remove_full = TRUE)
#' as.data.frame(obj)
#' }
#' @export
setGeneric(
  "results_df",
  function(object, collapse_full = TRUE, remove_full = TRUE, suffix = "_str") {
    standardGeneric("results_df")
  }
)


#' @export
setMethod("results_df", "DirectLibMatch",
          function(object, collapse_full = TRUE, remove_full = TRUE, suffix = "_merge") {

            x <- object@results
            if (!is.list(x) || length(x) == 0L) {
              return(data.frame())
            }

            full_fields <- c("Compound_peak_mz_full", "Compound_peak_int_full")

            rows <- lapply(x, function(el) {
              if (is.null(el)) return(as.data.frame(list(), stringsAsFactors = FALSE))

              # collapse *_full vectors as requested
              if (collapse_full) {
                for (nm in full_fields) {
                  if (!is.null(el[[nm]])) {
                    # collapse as character (works for numeric too)
                    collapsed <- paste(as.character(el[[nm]]), collapse = ", ")
                    if (remove_full) {
                      el[[nm]] <- NULL
                      el[[paste0(nm, suffix)]] <- collapsed
                    } else {
                      el[[nm]] <- collapsed
                    }
                  }
                }
              } else if (remove_full) {

                el[full_fields] <- NULL
              }

              as.data.frame(el, stringsAsFactors = FALSE, optional = TRUE)
            })

            all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
            rows <- lapply(rows, function(df) {
              miss <- setdiff(all_names, names(df))
              if (length(miss)) df[miss] <- NA
              df[all_names]
            })

            out <- do.call(rbind, rows)
            rownames(out) <- NULL
            out
          }
)

#' @export
as.data.frame.DirectLibMatch <- function(x, ..., collapse_full = TRUE,
                                         remove_full = TRUE, suffix = "_merge") {
  results_df(x, collapse_full = collapse_full, remove_full = remove_full, suffix = suffix)
}
