setMethod("show", "RTSimilarity", function(object) {
  cat("RTSimilarity object (method:", object@method, ")\n")

  if (is.matrix(object@sim_matrix)) {
    cat(" - Mode: full matrix\n")
    cat(" - Features:", nrow(object@sim_matrix), "\n")
  } else {
    cat(" - Mode: pairwise\n")
    cat(" - Number of pairs:", nrow(object@sim_matrix), "\n")
  }
})
