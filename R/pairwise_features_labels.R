#' Generate Pairwise Labels for All Features
#'
#' @description
#' For all features in the feature table, creates pairs and labels them:
#' 1 if same compound (from compound_list), else 0.
#'
#' @param feature_table A data.frame with at least column `feature_name`
#' @param compound_list A list where each element contains:
#'   - Compound_name: character
#'   - Matched_features_to_peaks: character vector of feature names
#'
#' @return A data.frame with columns: feature_1, feature_2, label (1 or 0)
#' @examples
#' ml_input <- pairwise_FTlabels(result_StM16_FT_filtered, direct_match_Stm16_res_MS1)
#' head(ml_input)
pairwise_FTlabels <- function(feature_table, compound_list) {
  # Flatten compound-feature mapping
  feature_compound_df <- do.call(rbind, lapply(compound_list, function(comp) {
    data.frame(
      feature = comp$Matched_features_to_peaks,
      compound = comp$Compound_name,
      stringsAsFactors = FALSE
    )
  }))

  # Add unmatched features as "Unknown"
  all_features <- feature_table$feature_name
  unmatched_features <- setdiff(all_features, feature_compound_df$feature)
  unmatched_df <- data.frame(
    feature = unmatched_features,
    compound = "Unknown",
    stringsAsFactors = FALSE
  )

  # Combine known and unknown mappings
  full_mapping <- rbind(feature_compound_df, unmatched_df)

  # Generate all feature pairs
  all_pairs <- expand.grid(
    feature_1 = full_mapping$feature,
    feature_2 = full_mapping$feature,
    stringsAsFactors = FALSE
  )

  # Remove self-pairs and duplicate pairs
  all_pairs <- all_pairs[all_pairs$feature_1 < all_pairs$feature_2, ]

  # Merge compound labels
  all_pairs <- merge(all_pairs, full_mapping, by.x = "feature_1", by.y = "feature")
  all_pairs <- merge(all_pairs, full_mapping, by.x = "feature_2", by.y = "feature", suffixes = c("_1", "_2"))

  # Assign labels
  all_pairs$label <- ifelse(
    all_pairs$compound_1 != "Unknown" &
      all_pairs$compound_2 != "Unknown" &
      all_pairs$compound_1 == all_pairs$compound_2,
    1, # Same known compound
    0  # Otherwise
  )

  # Select relevant columns
  result <- all_pairs[, c("feature_1", "feature_2", "label")]

  return(result)
}
