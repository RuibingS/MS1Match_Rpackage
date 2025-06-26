#' Validate and normalize polarity input
#'
#' @param polarity A character string, e.g., "positive", "negative", "pos", or "neg"
#' @return Normalized polarity string: "pos" or "neg"
#' @export
validate_polarity <- function(polarity) {
  polarity <- tolower(polarity)
  if (polarity %in% c("positive", "pos")) return("pos")
  if (polarity %in% c("negative", "neg")) return("neg")
  stop("`polarity` must be 'positive' or 'negative'")
}

