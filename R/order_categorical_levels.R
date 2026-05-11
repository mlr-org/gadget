#' Order Levels of a Categorical Feature
#'
#' Given x_cat (factor), data, feature, target_feature_name,
#' order_method: builds KxK distance matrix from other features (ECDF diff for numeric, L1 for categorical);
#' embeds in 1D via MDS/PCA/random/raw; reorders levels.
#' Returns factor with reordered levels (or x_cat unchanged if K<=1 or no other features).
#'
#' @param x_cat (`factor()`) \cr
#'   Focal feature values; use \code{droplevels()} first.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Full dataset including target.
#' @param feature (`character(1)`) \cr
#'   Name of the categorical feature in \code{data}.
#' @param target_feature_name (`character(1)`) \cr
#'   Target column name; excluded from distance computation.
#' @param order_method (`character(1)`) \cr
#'   \code{"mds"}, \code{"pca"}, \code{"random"}, or \code{"raw"}.
#'
#' @return (`factor()`) \cr
#'   Same as \code{x_cat} with reordered \code{levels} and \code{ordered = TRUE};
#'   unchanged if \code{nlevels(x_cat) <= 1} or no other features.
#'
#' @details
#' For each pair of levels, a distance is computed from all other features
#' (all columns in \code{data} except \code{feature} and \code{target_feature_name}):
#' \itemize{
#'   \item Numeric: sum over features of max absolute difference of ECDFs (KS-style).
#'   \item Categorical: sum over features of L1 distance of conditional level distributions.
#' }
#' The \eqn{K \times K} distance matrix is embedded in 1D via \code{order_method};
#' that 1D order defines the new level order. Single-level factors or no other
#' features yield \code{x_cat} unchanged.
#'
#' @keywords internal
order_categorical_levels = function(x_cat, data, feature, target_feature_name, order_method = "raw") {
  checkmate::assert_factor(x_cat, .var.name = "x_cat")
  checkmate::assert_data_frame(data, .var.name = "data")
  checkmate::assert_character(feature, len = 1, .var.name = "feature")
  checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
  checkmate::assert_subset(c(feature, target_feature_name), colnames(data),
    .var.name = "feature and target_feature_name")
  checkmate::assert_choice(order_method, c("mds", "pca", "random", "raw"), .var.name = "order_method")
  levels_orig = levels(x_cat)
  K = nlevels(x_cat)
  if (order_method == "raw") {
    return(factor(x_cat, levels = levels_orig, ordered = TRUE))
  }
  # Use all features except the focal feature and the target
  other_feats = setdiff(colnames(data), c(feature, target_feature_name))
  if (K <= 1L || length(other_feats) == 0L) {
    return(x_cat)
  }
  # Build all level pairs
  pairs = expand.grid(from = levels_orig, to = levels_orig, stringsAsFactors = FALSE)
  total_dists = numeric(nrow(pairs))

  for (feat in other_feats) {
    val = data[[feat]]
    if (all(is.na(val))) next
    if (is.numeric(val)) {
      # Numeric other feature: compare ECDFs via KS distance
      grid = stats::quantile(val, probs = seq(0, 1, length.out = 100), na.rm = TRUE, names = FALSE)
      ecdf_list = tapply(val, x_cat, function(v) {
        if (length(stats::na.omit(v)) == 0L) {
          return(numeric(100))
        }
        stats::ecdf(v)(grid)
      })
      ecdf_mat = do.call(cbind, ecdf_list)[, levels_orig, drop = FALSE]
      diffs = abs(ecdf_mat[, pairs$from] - ecdf_mat[, pairs$to])
      total_dists = total_dists + apply(diffs, 2, max)
    } else {
      # Categorical / discrete other feature: compare conditional distributions via L1 distance
      tbl = table(x_cat, val)
      prob_mat = tbl / ifelse(rowSums(tbl) == 0, 1, rowSums(tbl))
      diffs = abs(prob_mat[pairs$from, ] - prob_mat[pairs$to, ])
      total_dists = total_dists + rowSums(diffs)
    }
  }
  # Construct K x K distance matrix
  dist_mat = matrix(total_dists, nrow = K, dimnames = list(levels_orig, levels_orig))

  if (order_method == "random") {
    sorted_lev = sample(levels_orig, length(levels_orig), replace = FALSE)
  } else if (order_method == "mds") {
    mds_ord = stats::cmdscale(dist_mat, k = 1)
    scores = mds_ord[, 1]
    sorted_lev = rownames(dist_mat)[order(scores)]
  } else if (order_method == "pca") {
    pca_res = stats::prcomp(dist_mat, center = TRUE, scale. = FALSE)
    scores = pca_res$x[, 1]
    sorted_lev = rownames(dist_mat)[order(scores)]
  } else {
    # Fallback: keep original order if unknown method is passed
    return(x_cat)
  }

  factor(x_cat, levels = sorted_lev, ordered = TRUE)
}
