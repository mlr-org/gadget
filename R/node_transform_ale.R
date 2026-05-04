#' Node Transform ALE
#'
#' Subsets ALE effect data to the current node's row indices and recomputes
#' per-interval statistics. When \code{split_feature} is non-\code{NULL},
#' forces \code{d_l = 0} for any feature whose values are constant in this node
#' (single unique value).
#'
#' @param Y (`list()`) \cr
#'   ALE effect data per feature.
#' @param idx (`integer()`) \cr
#'   Sample indices in the current node.
#' @param split_feature (`character(1)` or `NULL`) \cr
#'   Feature used for splitting; \code{NULL} = no postprocessing.
#'
#' @return (`list()`) \cr
#'   Transformed ALE effects per feature.
#'
node_transform_ale = function(Y, idx, split_feature = NULL) {
  y_subset = lapply(names(Y), function(feat) {
    y_j = Y[[feat]]
    y_j = y_j[y_j$row_id %in% idx, ]
    y_j[, `:=`(
      int_n  = .N,
      int_s1 = sum(d_l, na.rm = TRUE),
      int_s2 = sum(d_l^2, na.rm = TRUE)
    ), by = interval_index]
    y_j
  })
  names(y_subset) = names(Y)
  if (!is.null(split_feature)) {
    y_processed = lapply(names(y_subset), function(feat) {
      y_j = y_subset[[feat]]
      # Zero out d_l for constant features in this node (ALE undefined when all values equal)
      if (length(unique(y_j$feat_val)) == 1) {
        y_j$d_l = 0
        y_j[, `:=`(
          int_n  = .N,
          int_s1 = sum(d_l, na.rm = TRUE),
          int_s2 = sum(d_l^2, na.rm = TRUE)
        ), by = interval_index]
      }
      y_j
    })
    names(y_processed) = names(y_subset)
    y_processed
  } else {
    y_subset
  }
}
