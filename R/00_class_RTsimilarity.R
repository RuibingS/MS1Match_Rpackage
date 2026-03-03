#' Retention Time Similarity Object
#'
#' @description
#' Stores RT distances and similarity information for features.
#'
#' Can handle either:
#'   1. Full `n x n` distance & similarity matrices, OR
#'   2. Pairwise feature RT distance/similarity for specified feature pairs.
#'
#' @slot dist_matrix
#'   - If full calculation: `matrix` of pairwise RT distances (n x n).
#'   - If pairwise calculation: `data.frame` with columns:
#'       * `feature_1`: character, first feature ID
#'       * `feature_2`: character, second feature ID
#'       * `rt_dist`: numeric, absolute RT difference
#' @slot sim_matrix
#'   - If full calculation: `matrix` of RT similarity scores (n x n).
#'   - If pairwise calculation: `data.frame` with columns:
#'       * `feature_1`: character, first feature ID
#'       * `feature_2`: character, second feature ID
#'       * `rt_similarity`: numeric, similarity score
#' @slot method
#'   - Character string indicating the similarity method used (`"gaussian"` or `"inverse"`).
#'
#' @seealso See function `calRTSimilarity()`.
#'
#' @export
setClass("RTSimilarity",
         slots = list(
           dist_matrix = "ANY",
           sim_matrix  = "ANY",
           method      = "character"
         ))
