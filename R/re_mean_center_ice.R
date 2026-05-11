#' Mean-center ICE matrices for a node. Given Y (list of effect matrices), grid (column names per feature),
#' idx: subsets rows by idx; sets non-grid columns to NA; subtracts row means.
#' Returns list of mean-centered matrices.
#'
#' @param Y (`list()`) \cr
#'   Effect matrices per feature.
#' @param grid (`list()`) \cr
#'   Grid column names per feature.
#' @param idx (`integer()`) \cr
#'   Row indices to keep.
#'
#' @return (`list()`) \cr
#'   Mean-centered matrices per feature.
#' @keywords internal
re_mean_center_ice = function(Y, grid, idx) {
  feature_names = names(Y)
  y_centered = lapply(seq_along(Y), function(i) {
    Y[[i]] = Y[[i]][idx, ]
    Y[[i]][, which(!(colnames(Y[[i]]) %in% grid[[i]]))] = NA
    Y[[i]] - rowMeans(Y[[i]], na.rm = TRUE)
  })
  names(y_centered) = feature_names
  y_centered
}
