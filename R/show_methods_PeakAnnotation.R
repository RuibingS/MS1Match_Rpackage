setMethod("show", "PeakAnnotationList", function(object) {
  cat("S4 Object of class 'PeakAnnotationList'\n")
  cat("- Method:", object@method, "\n")
  cat("- Number of features:", nrow(object@FT), "\n")
  cat("- Polarity:", object@polarity, "\n")
})
