#' PeakAnnotationList class
#'
#' @description S4 subclass of \code{PeakAnnotation} for compound list-based annotation.
#' @slot compounds_list A named list of compounds
#' @export
setClass("PeakAnnotationList",
         slots = c(compounds_list = "list"),
         contains = "PeakAnnotation")

