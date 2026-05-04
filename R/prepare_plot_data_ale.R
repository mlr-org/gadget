utils::globalVariables(c("interval_index", "x_left", "x_right", "d_l", "x_grid", "level"))
interval_index = x_left = x_right = d_l = x_grid = level = NULL

#' Prepare ALE Plot Data for One or More Nodes
#'
#' Given effect (from \code{calculate_ale}), idx (row indices or list of such), features,
#' mean_center: subsets ALE rows by idx; calls \code{mean_center_ale} per feature for cumulative and optional centering.
#' Returns named list of \code{mean_effect} data.tables (or nested list if idx is list).
#'
#' @param effect (`list()`) \cr
#'   List returned by \code{calculate_ale()}.
#' @param idx (`integer()` or `list()` or `NULL`) \cr
#'   Row indices (node subset); list of such vectors; or \code{NULL} for root.
#' @param features (`character()`) \cr
#'   Features to include (default: all in effect).
#' @param mean_center (`logical(1)`) \cr
#'   Whether to mean-center ALE curves.
#'
#' @details
#' Rows are subset with \code{effect[[feat]][row_id \%in\% idx]} when \code{idx} is non-\code{NULL}.
#'
#' \code{mean_center_ale()} builds plot grids from aggregated intervals: means sample-wise \code{d_l}
#' within each \code{(interval_index, x_left, x_right)} group (\code{delta_aggr}), cumulates, then
#' optionally subtracts a weighted scalar \code{f_j0}. Sample-wise \code{d_l == 0} is mapped to
#' \code{NA} before aggregation so exact zeros do not enter group means.
#'
#' For categorical features, each row of \code{mean_effect} corresponds to one row of
#' \code{delta_aggr}; \code{x_grid} uses \code{as.character(delta_aggr$x_left)}, so the number of
#' plotted points follows the number of distinct aggregated intervals after subsetting (not always one
#' row per factor level). Factor \code{x_grid} still carries full \code{levels(feat_val)} for axis
#' ordering.
#'
#' @return (`list()`) \cr
#'   Named list of \code{mean_effect} data.tables per feature; nested if \code{idx} is list.
prepare_plot_data_ale = function(effect, idx = NULL, features = names(effect),
  mean_center = TRUE) {
  if (is.null(effect) || !length(effect)) {
    return(list())
  }
  if (is.list(idx) && length(idx) && !is.atomic(idx)) {
    idx_names = names(idx)
    res = lapply(
      idx,
      function(one_idx) {
        prepare_plot_data_ale(
          effect,
          idx = one_idx,
          features = features,
          mean_center = mean_center
        )
      }
    )
    if (!is.null(idx_names)) names(res) = idx_names
    return(res)
  }
  feats = intersect(features, names(effect))
  if (!length(feats)) {
    return(list())
  }
  out = lapply(feats, function(feat) {
    dt = data.table::as.data.table(effect[[feat]])
    if (!is.null(idx)) {
      dt = dt[dt[["row_id"]] %in% idx]
    }
    mean_dt = mean_center_ale(dt, mean_center = mean_center)
    mean_dt$feature = feat
    list(mean_effect = mean_dt)
  })
  names(out) = feats
  out
}

#' Internal ALE curve computation
#'
#' Given ALE data.table for one feature and mean_center: cumsums d_l by interval;
#' optionally subtracts global mean. Returns data.table with x_grid and .value (cumulative ALE).
#'
#' @param feat (`data.table()`) \cr
#'   Per-interval ALE derivatives and metadata.
#' @param mean_center (`logical(1)`) \cr
#'   Whether to mean-center the ALE curve.
#'
#' @return (`data.table()`) \cr
#'   Cumulative ALE with \code{x_grid} and cumulative values in column \code{d_l}.
#' @keywords internal
mean_center_ale = function(feat, mean_center = TRUE) {
  data.table::setkeyv(feat, c("interval_index"))
  delta_aggr = feat[, list(
    d_l = mean(.SD[[1]], na.rm = TRUE),
    interval_n = sum(!is.na(d_l))
  ), by = c("interval_index", "x_left", "x_right"), .SDcols = "d_l"]

  vals = delta_aggr$d_l
  csum = cumsum_na_as_zero(c(0, vals))
  weights = delta_aggr$interval_n
  mid_vals = (csum[-length(csum)] + csum[-1]) / 2
  denom = sum(weights)
  f_j0 = if (denom > 0) sum(mid_vals * weights) / denom else 0
  if (!isTRUE(mean_center)) f_j0 = 0

  if (is.numeric(feat$feat_val)) {
    x_grid = c(delta_aggr$x_left, tail(delta_aggr$x_right, 1))
    mean_dt = data.table::data.table(
      feature = NA_character_,
      x_grid = x_grid,
      d_l = csum - f_j0
    )
  } else if (is.factor(feat$feat_val)) {
    levs = levels(feat$feat_val)
    if (is.null(levs)) levs = unique(as.character(feat$feat_val))
    level_vals = as.character(delta_aggr$x_left)
    mean_dt = data.table::data.table(
      feature = NA_character_,
      x_grid = level_vals,
      d_l = csum[-1] - f_j0
    )
    mean_dt$x_grid = factor(mean_dt$x_grid, levels = levs)
  } else {
    mean_dt = data.table::data.table(feature = NA_character_, x_grid = numeric(0), d_l = numeric(0))
  }

  mean_dt
}

#' Cumulative sum with NA as zero.
#'
#' @param values (`numeric()`) \cr
#'   Values to cumulatively sum; \code{NA} treated as 0.
#' @return (`numeric()`) \cr
#'   Cumulative sum.
#' @keywords internal
cumsum_na_as_zero = function(values) {
  values[is.na(values)] = 0
  cumsum(values)
}
