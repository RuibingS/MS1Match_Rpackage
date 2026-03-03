#' Flatten PI target match result list to a data.frame
#'
#' @param res_list A list returned by PImatch_target()
#' @return A data.frame with one row per matched feature-compound pair
#' @export
PImatch_df <- function(res_list) {

  if (!is.list(res_list)) {
    stop("`res_list` must be a list.")
  }

  # keep only non-empty matches
  res_list <- Filter(function(x) {
    !is.null(x$Feature_name) && length(x$Feature_name) > 0
  }, res_list)

  if (length(res_list) == 0) {
    return(data.frame(
      Feature_name = character(),
      Feature_mz   = numeric(),
      Comp_name    = character(),
      PI_name      = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(res_list, function(x) {

    n <- length(x$Feature_name)

    data.frame(
      Feature_name = x$Feature_name,
      Feature_mz   = x$Feature_mz,
      Comp_name    = x$Comp_name,
      PI_name      = x$PI_name,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
