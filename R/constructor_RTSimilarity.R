#' @export
calRTSimilarity <- function(feature_table,
                            sigma = 5,
                            feature_pairs = NULL,
                            method = c("gaussian", "inverse")) {
  method <- match.arg(method)

  stopifnot(is.data.frame(feature_table))
  stopifnot("feature_name" %in% colnames(feature_table))

  rt_col <- intersect(c("rt", "rtmed"), colnames(feature_table))
  if (length(rt_col) == 0) {
    stop("feature_table must contain either a 'rt' or 'rtmed' column.")
  }

  feature_names <- as.character(feature_table$feature_name)
  rt_values <- as.numeric(feature_table[[rt_col[1]]])

  if (is.null(feature_pairs)) {
    # Full matrix case
    result <- switch(method,
                     gaussian = gaussian_kernel_dist(rt_values, sigma),
                     inverse  = inverse_scaling_dist(rt_values))

    new("RTSimilarity",
        dist_matrix = result$dist_matrix,
        sim_matrix  = result$sim_matrix,
        method = method)

  } else {
    # Pairwise case
    attr(feature_pairs, "feature_names") <- feature_names

    sim_df <- switch(method,
                     gaussian = gaussian_kernel_dist(rt_values, sigma, feature_pairs),
                     inverse  = inverse_scaling_dist(rt_values, feature_pairs))

    feature_index <- setNames(seq_along(feature_names), feature_names)
    sim_df$rt_dist <- abs(rt_values[feature_index[sim_df$feature_1]] -
                            rt_values[feature_index[sim_df$feature_2]])

    new("RTSimilarity",
        dist_matrix = sim_df[c("feature_1", "feature_2", "rt_dist")],
        sim_matrix  = sim_df[c("feature_1", "feature_2", "rt_similarity")],
        method = method)
  }
}
