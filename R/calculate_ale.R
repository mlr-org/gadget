#' Calculate Accumulated Local Effects (ALE)
#'
#' Given model, data, feature_set, target_feature_name, n_intervals, predict_fun: for each feature,
#' computes finite differences (d_l) and per-interval stats (int_n, int_s1, int_s2).
#' Numeric: quantile intervals; categorical: level-by-level prediction differences.
#' Returns named list of data.tables (row_id, feat_val, d_l, interval_index, int_n, int_s1, int_s2, etc.).
#'
#' @param model (`any`) \cr
#'   Fitted model with predict interface.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data.
#' @param feature_set (`character()`) \cr
#'   Features to compute ALE for.
#' @param target_feature_name (`character(1)`) \cr
#'   Target variable name.
#' @param n_intervals (`integer(1)`) \cr
#'   Equal-frequency intervals for numeric features.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   \code{function(model, data)} returning predictions; \code{NULL} = default.
#'
#' @return (`list()`) \cr
#'   Named list of data.tables per \code{feature_set}. Each has columns:
#'   \item{row_id}{Row index in \code{data}.}
#'   \item{feat_val}{Feature value at that row.}
#'   \item{x_left, x_right}{Interval/category boundaries (numeric) or left/right category (factor).}
#'   \item{d_l}{Local effect (finite difference).}
#'   \item{interval_index}{Interval or category index.}
#'   \item{int_n, int_s1, int_s2}{Per-interval count and sum(d_l), sum(d_l^2) for heterogeneity.}
#'
#' @details
#' Numeric features: builds \code{n_intervals} quantile-based intervals,
#' assigns each row to an interval, and computes finite differences between
#' interval boundaries via \code{predict_fun}.
#'
#' Categorical features: use factor levels as given (typically pre-ordered by
#' \code{order_categorical_levels} in \code{prepare_split_data_ale}). For each
#' row, \code{d_l} is the difference in predictions when the focal feature is
#' set to the next vs. previous level; single-level factors get \code{d_l = 0}.
#'
#' Sample-level columns (\code{row_id}, \code{feat_val}, \code{d_l}, etc.)
#' support subsetting by node and downstream heterogeneity calculation.
#'
#' Downstream plotting (\code{\link{prepare_plot_data_ale}}) aggregates these rows by
#' \code{(interval_index, x_left, x_right)}, cumulates \code{d_l}, and optionally mean-centers the
#' cumulative curve; plot grids for categories derive from aggregated \code{x_left} values, not by
#' re-evaluating \code{calculate_ale}.
#'
#' @keywords internal
calculate_ale = function(model, data, feature_set, target_feature_name, n_intervals = 10, predict_fun = NULL) {
  checkmate::assert_data_frame(data, .var.name = "data")
  checkmate::assert_character(feature_set, min.len = 1, .var.name = "feature_set")
  checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
  checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
  checkmate::assert_subset(feature_set, colnames(data), .var.name = "feature_set")
  checkmate::assert_integerish(n_intervals, len = 1, lower = 1, any.missing = FALSE, .var.name = "n_intervals")
  checkmate::assert_function(predict_fun, null.ok = TRUE, .var.name = "predict_fun")
  if (is.null(predict_fun)) predict_fun = default_predict_fun
  if (data.table::is.data.table(data)) {
    X = data[, setdiff(colnames(data), target_feature_name), with = FALSE]
  } else {
    X = data[, setdiff(colnames(data), target_feature_name), drop = FALSE]
  }
  X_dt = data.table::as.data.table(X)
  n_rows = nrow(X_dt)
  stacked_shared = data.table::rbindlist(list(X_dt, X_dt), use.names = TRUE)
  idx_lower = seq_len(n_rows)
  idx_upper = seq.int(n_rows + 1L, 2L * n_rows)

  eff_list = lapply(feature_set, function(feat) {
    if (is.factor(data[[feat]])) {
      ale_categorical_feature(model = model, data = data, X = X_dt,
        feature = feat, predict_fun = predict_fun,
        stacked = stacked_shared, idx_lower = idx_lower, idx_upper = idx_upper)
    } else {
      ale_numeric_feature(model = model, data = data, X = X_dt,
        feature = feat, n_intervals = n_intervals, predict_fun = predict_fun,
        stacked = stacked_shared, idx_lower = idx_lower, idx_upper = idx_upper)
    }
  })
  names(eff_list) = feature_set
  eff_list
}

#' ALE for a single numeric feature.
#'
#' @param model (`any`) \cr
#'   Fitted model. See \code{\link{calculate_ale}}.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data.
#' @param X (`data.frame()` or `data.table()`) \cr
#'   Features (excl. target).
#' @param feature (`character(1)`) \cr
#'   Feature name.
#' @param target_feature_name (`character(1)`) \cr
#'   Target name.
#' @param n_intervals (`integer(1)`) \cr
#'   Number of intervals.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   Prediction function.
#'
#' @return (`data.table()`) \cr
#'   ALE data with \code{row_id}, \code{feat_val}, \code{d_l}, \code{interval_index}, etc.
#' @param stacked (`NULL` or [data.table::data.table()]) \cr
#'   Shared \code{2n}-row design matrix; omit to allocate internally.
#' @param idx_lower, idx_upper (`integer()` or \code{NULL}) \cr
#'   Lower/upper half row indices inside \code{stacked}.
#' @keywords internal
ale_numeric_feature = function(model, data, X, feature, n_intervals = 10, predict_fun = NULL,
  stacked = NULL, idx_lower = NULL, idx_upper = NULL) {
  x_num = data[[feature]]
  n_rows = nrow(data)
  if (is.null(stacked)) {
    xd = data.table::as.data.table(X)
    stacked = data.table::rbindlist(list(xd, xd), use.names = TRUE)
    idx_lower = seq_len(n_rows)
    idx_upper = seq.int(n_rows + 1L, 2L * n_rows)
  }

  if (length(unique(na.omit(x_num))) <= 1L) {
    return(ale_zero(feat_val = x_num))
  }
  q = stats::quantile(x_num, 0:n_intervals / n_intervals, type = 7, na.rm = TRUE)
  if (length(unique(q)) < 2L) {
    return(ale_zero(feat_val = x_num))
  }
  interval_index = findInterval(x_num, q, left.open = TRUE)
  interval_index[interval_index == 0L] = 1L
  max_id = length(q) - 1L
  interval_index[interval_index > max_id] = max_id
  original = stacked[[feature]][idx_lower]
  data.table::set(stacked, i = idx_lower, j = feature, value = q[interval_index])
  data.table::set(stacked, i = idx_upper, j = feature, value = q[interval_index + 1L])
  pred_raw = predict_fun(model, stacked)
  pred = extract_numeric_prediction(pred_raw, expected_n = 2L * n_rows)
  d_l = pred[idx_upper] - pred[idx_lower]
  data.table::set(stacked, i = idx_lower, j = feature, value = original)
  data.table::set(stacked, i = idx_upper, j = feature, value = original)

  DT = data.table::data.table(
    row_id         = seq_len(n_rows),
    feat_val       = x_num,
    x_left         = q[interval_index],
    x_right        = q[interval_index + 1L],
    d_l             = d_l,
    interval_index = interval_index
  )
  DT[, `:=`(
    int_n  = .N,
    int_s1 = sum(d_l),
    int_s2 = sum(d_l^2)
  ), by = interval_index]
  DT
}

#' ALE for a single categorical feature.
#'
#' @param model (`any`) \cr
#'   Fitted model. See \code{\link{calculate_ale}}.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data.
#' @param X (`data.frame()` or `data.table()`) \cr
#'   Features (excl. target).
#' @param feature (`character(1)`) \cr
#'   Feature name.
#' @param target_feature_name (`character(1)`) \cr
#'   Target name.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   Prediction function.
#'
#' @return (`data.table()`) \cr
#'   ALE data with \code{row_id}, \code{feat_val}, \code{d_l}, \code{interval_index}, etc.
#' @param stacked (`NULL` or [data.table::data.table()]) \cr
#'   Shared \code{2n}-row design matrix for batched categorical ALE (see numeric branch).
#' @param idx_lower, idx_upper (`integer()` or \code{NULL}) \cr
#'   Row halves in \code{stacked}: plus-vector / minus-vector predictions respectively.
#' @keywords internal
ale_categorical_feature = function(model, data, X, feature, predict_fun = NULL,
  stacked = NULL, idx_lower = NULL, idx_upper = NULL) {
  x_cat = droplevels(data[[feature]])
  K = nlevels(x_cat)
  n_rows = nrow(data)
  if (is.null(stacked)) {
    xd = data.table::as.data.table(X)
    stacked = data.table::rbindlist(list(xd, xd), use.names = TRUE)
    idx_lower = seq_len(n_rows)
    idx_upper = seq.int(n_rows + 1L, 2L * n_rows)
  }

  if (K <= 1) {
    return(ale_zero(feat_val = x_cat, interval_index = as.integer(x_cat)))
  }
  levels_orig = levels(x_cat)
  levels_id = as.numeric(x_cat)
  row_ind_plus = seq_len(n_rows)[levels_id < K]
  row_ind_neg  = seq_len(n_rows)[levels_id > 1]

  original = stacked[[feature]][idx_lower]
  cp = X[[feature]]
  cn = X[[feature]]
  cp[row_ind_plus] = levels_orig[levels_id[row_ind_plus] + 1L]
  cn[row_ind_neg] = levels_orig[levels_id[row_ind_neg] - 1L]
  data.table::set(stacked, i = idx_lower, j = feature, value = cp)
  data.table::set(stacked, i = idx_upper, j = feature, value = cn)

  pred_cat = extract_numeric_prediction(
    predict_fun(model, stacked),
    expected_n = 2L * n_rows
  )
  delta = pred_cat[idx_lower] - pred_cat[idx_upper]
  data.table::set(stacked, i = idx_lower, j = feature, value = original)
  data.table::set(stacked, i = idx_upper, j = feature, value = original)

  DT = data.table::data.table(
    row_id         = seq_len(n_rows),
    feat_val       = data[[feature]],
    x_left         = cn,
    x_right        = cp,
    d_l             = delta,
    interval_index = levels_id
  )

  DT[, `:=`(
    int_n  = .N,
    int_s1 = sum(d_l, na.rm = TRUE),
    int_s2 = sum(d_l^2, na.rm = TRUE)
  ), by = interval_index]
  DT
}
