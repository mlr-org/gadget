#' Calculate Partial Dependence Curves
#'
#' Computes ICE (Individual Conditional Expectation) matrices for each feature
#' in \code{feature_set} and returns them in long-format data.tables.
#'
#' @param model (`any`) \cr
#'   Fitted model with a predict interface.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data including the target column.
#' @param target_feature_name (`character(1)`) \cr
#'   Name of the target variable; excluded from feature columns.
#' @param feature_set (`character()` or `NULL`) \cr
#'   Features to compute PD for; \code{NULL} = all non-target columns.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   \code{function(model, data)} returning a numeric vector; \code{NULL} = default.
#' @param n_grid (`integer(1)`) \cr
#'   Number of quantile-based grid points for numeric features.
#' @param pd_engine (`character(1)`) \cr
#'   Backend: \code{"cpp"} (default) or \code{"r"}.
#'
#' @return (`list()`) \cr
#'   Named list with element \code{results}: a named list of data.tables, one per
#'   feature, each with columns \code{.id}, \code{.type}, \code{.feature},
#'   \code{.borders}, \code{.value}.
#'
#' @keywords internal
calculate_pd = function(model, data, target_feature_name, feature_set = NULL,
  predict_fun = NULL, n_grid = 20L, pd_engine = c("cpp", "r")) {
  pd_engine = match.arg(pd_engine)
  features = setdiff(colnames(data), target_feature_name)
  feature_set = resolve_split_features(feature_set, features, "Features")
  if (length(feature_set) == 0L) {
    cli::cli_abort("{.arg feature_set} must contain at least one feature.")
  }
  x_features = if (data.table::is.data.table(data)) {
    data[, features, with = FALSE]
  } else {
    data[, features, drop = FALSE]
  }
  x_features_dt = data.table::as.data.table(x_features)
  x_cols_list = as.list(x_features_dt)
  n_obs = nrow(x_features_dt)

  grids = stats::setNames(
    nm = feature_set,
    lapply(feature_set, function(feat) pd_feature_grid(x_features[[feat]], n_grid = n_grid))
  )

  # Pre-allocate one stacked table for the R path: max_g * n rows, reused per feature.
  max_g_r = max(lengths(grids))
  stacked_pd_cache = NULL
  if (identical(pd_engine, "r")) {
    stacked_pd_cache = list(
      stacked = data.table::as.data.table(lapply(x_features_dt, rep, times = max_g_r)),
      max_g = max_g_r,
      n_obs = n_obs
    )
  }

  results = mlr3misc::map(setNames(nm = feature_set), function(feat) {
    grid = grids[[feat]]
    feat_index = match(feat, names(x_features_dt))
    ice = compute_ice(
      model = model, data = x_features_dt, feature = feat, grid = grid,
      predict_fun = predict_fun, pd_engine = pd_engine,
      base_data_dt = x_features_dt, cols_list = x_cols_list, feature_index = feat_index,
      stacked_pd_cache = stacked_pd_cache
    )
    pd_pack_ice_result(ice, feature = feat, grid = grid)
  })

  list(results = results)
}


#' Compute ICE Matrix (Dispatch)
#'
#' Dispatches ICE computation to the C++ or R backend based on \code{pd_engine}.
#'
#' @param model (`any`) \cr
#'   Fitted model.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Feature data (target column already removed).
#' @param feature (`character(1)`) \cr
#'   Name of the focal feature.
#' @param grid (`atomic vector`) \cr
#'   Grid values for the focal feature.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   Custom predict function; \code{NULL} = default.
#' @param pd_engine (`character(1)`) \cr
#'   \code{"cpp"} or \code{"r"}.
#' @param base_data_dt (`data.table()` or `NULL`) \cr
#'   Pre-converted data.table of \code{data}; avoids repeated conversion.
#' @param cols_list (`list()` or `NULL`) \cr
#'   Pre-extracted column list of \code{base_data_dt}; used by the C++ path.
#' @param feature_index (`integer(1)` or `NULL`) \cr
#'   1-based column index of \code{feature} in \code{base_data_dt}; used by the C++ path.
#' @param stacked_pd_cache (`list()` or `NULL`) \cr
#'   Pre-allocated stacked data.table cache for the R path; \code{NULL} disables caching.
#'
#' @return (`matrix`) \cr
#'   Numeric matrix of shape \code{n_obs x length(grid)} containing ICE predictions.
#'
#' @keywords internal
compute_ice = function(
  model, data, feature, grid, predict_fun = NULL,
  pd_engine = c("cpp", "r"), base_data_dt = NULL,
  cols_list = NULL, feature_index = NULL,
  stacked_pd_cache = NULL
) {
  pd_engine = match.arg(pd_engine)
  if (identical(pd_engine, "cpp")) {
    compute_ice_cpp(
      model = model, data = data, feature = feature, grid = grid,
      predict_fun = predict_fun, base_data_dt = base_data_dt,
      cols_list = cols_list, feature_index = feature_index
    )
  } else {
    compute_ice_r(
      model = model, data = data, feature = feature, grid = grid,
      predict_fun = predict_fun, base_data_dt = base_data_dt,
      stacked_pd_cache = stacked_pd_cache
    )
  }
}


#' Compute ICE Matrix (Pure R)
#'
#' Builds a stacked prediction data.table by repeating each row once per grid
#' value, replaces the focal feature column with each grid value, runs
#' \code{pd_predict}, and reshapes predictions into a matrix.
#'
#' @param model (`any`) \cr
#'   Fitted model.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Feature data (target removed).
#' @param feature (`character(1)`) \cr
#'   Name of the focal feature.
#' @param grid (`atomic vector`) \cr
#'   Grid values for the focal feature.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   Custom predict function; \code{NULL} = default.
#' @param base_data_dt (`data.table()` or `NULL`) \cr
#'   Pre-converted data.table; avoids repeated conversion.
#' @param stacked_pd_cache (`list()` or `NULL`) \cr
#'   Pre-allocated stacked table with elements \code{stacked}, \code{max_g},
#'   \code{n_obs}; \code{NULL} disables caching.
#'
#' @return (`matrix`) \cr
#'   Numeric matrix of shape \code{n_obs x length(grid)}.
#'
#' @keywords internal
compute_ice_r = function(model, data, feature, grid, predict_fun = NULL,
  base_data_dt = NULL, stacked_pd_cache = NULL) {
  checkmate::assert_character(feature, len = 1L, .var.name = "feature")
  checkmate::assert_subset(feature, colnames(data), .var.name = "feature")
  checkmate::assert_atomic_vector(grid, min.len = 1L, .var.name = "grid")

  n_obs = nrow(data)
  grid_len = length(grid)
  dt = if (is.null(base_data_dt)) data.table::as.data.table(data) else base_data_dt
  use_cache = !is.null(stacked_pd_cache)

  if (use_cache) {
    checkmate::assert_list(stacked_pd_cache, names = "named")
    checkmate::assert_names(names(stacked_pd_cache), must.include = c("stacked", "max_g"))
    stacked = stacked_pd_cache$stacked
    if (grid_len > stacked_pd_cache$max_g) {
      cli::cli_abort("Internal PDP cache: grid length ({grid_len}) exceeds max ({stacked_pd_cache$max_g}).")
    }
    if (stacked_pd_cache$n_obs != n_obs) {
      cli::cli_abort("Internal PDP cache row count mismatch.")
    }
    n_take = n_obs * grid_len
    row_take = seq_len(n_take)
  } else {
    stacked = data.table::as.data.table(lapply(dt, rep, times = grid_len))
    n_take = nrow(stacked)
    row_take = seq_len(n_take)
  }

  feature_values = rep(grid, each = n_obs)
  if (is.factor(dt[[feature]])) {
    feature_values = factor(feature_values, levels = levels(dt[[feature]]))
  }
  focal_restore = dt[[feature]]
  if (!use_cache) {
    data.table::set(stacked, j = feature, value = feature_values)
    pred_slice = stacked
  } else {
    data.table::set(stacked, i = row_take, j = feature, value = feature_values)
    pred_slice = stacked[row_take]
  }

  pred = pd_predict(model, pred_slice, predict_fun = predict_fun)
  if (use_cache) {
    data.table::set(stacked, j = feature, value = rep(focal_restore, times = stacked_pd_cache$max_g))
  }

  matrix(pred, nrow = n_obs, ncol = grid_len)
}


#' Compute ICE Matrix (C++ Backend)
#'
#' Uses \code{cpp_pd_stack_newdata} to build the stacked prediction table in C++,
#' then runs \code{pd_predict}. Falls back to \code{compute_ice_r} for
#' character or logical feature columns.
#'
#' @param model (`any`) \cr
#'   Fitted model.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Feature data (target removed).
#' @param feature (`character(1)`) \cr
#'   Name of the focal feature.
#' @param grid (`atomic vector`) \cr
#'   Grid values for the focal feature.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   Custom predict function; \code{NULL} = default.
#' @param base_data_dt (`data.table()` or `NULL`) \cr
#'   Pre-converted data.table; avoids repeated conversion.
#' @param cols_list (`list()` or `NULL`) \cr
#'   Pre-extracted column list of \code{base_data_dt}.
#' @param feature_index (`integer(1)` or `NULL`) \cr
#'   1-based column index of \code{feature} in \code{base_data_dt}.
#'
#' @return (`matrix`) \cr
#'   Numeric matrix of shape \code{n_obs x length(grid)}.
#'
#' @keywords internal
compute_ice_cpp = function(
  model, data, feature, grid, predict_fun = NULL,
  base_data_dt = NULL, cols_list = NULL, feature_index = NULL
) {
  checkmate::assert_character(feature, len = 1L, .var.name = "feature")
  checkmate::assert_subset(feature, colnames(data), .var.name = "feature")
  checkmate::assert_atomic_vector(grid, min.len = 1L, .var.name = "grid")

  dt = if (is.null(base_data_dt)) data.table::as.data.table(data) else base_data_dt
  j = if (is.null(feature_index)) match(feature, names(dt)) else as.integer(feature_index)
  if (is.na(j)) {
    cli::cli_abort("Feature {.val {feature}} not found in {.arg data}.")
  }
  feat_col = dt[[feature]]

  # C kernel handles numeric, integer, factor; fall back to R for character/logical.
  if (is.character(feat_col) || is.logical(feat_col)) {
    return(compute_ice_r(model, data, feature, grid, predict_fun, base_data_dt = dt))
  }

  grid_sexp = if (is.factor(feat_col)) {
    gi = match(as.character(grid), levels(feat_col))
    if (anyNA(gi)) {
      cli::cli_abort("{.arg grid} values must match factor levels of {.field {feature}}.")
    }
    gi
  } else {
    as.numeric(grid)
  }

  cols_shared = if (is.null(cols_list)) as.list(dt) else cols_list
  col_lens = vapply(cols_shared, length, 1L)
  if (any(col_lens != nrow(dt))) {
    cli::cli_abort("All columns in {.arg data} must have length {.val {nrow(dt)}} for stacked PD prediction.")
  }
  stacked_df = cpp_pd_stack_newdata(cols_shared, j - 1L, grid_sexp)
  pred = pd_predict(model, stacked_df, predict_fun = predict_fun)
  n_obs = nrow(data)
  grid_len = length(grid)
  dim(pred) = c(n_obs, grid_len)
  pred
}


#' Build Feature Grid for Partial Dependence
#'
#' Returns grid values for a single feature column:
#' factor levels (after \code{droplevels}), unique sorted values for character,
#' or \code{n_grid} quantile-based numeric values.
#'
#' @param x (`vector`) \cr
#'   Feature column from the training data.
#' @param n_grid (`integer(1)`) \cr
#'   Number of grid points for numeric features; ignored for factor/character.
#'
#' @return (`atomic vector`) \cr
#'   Grid values: \code{character()} for factor/character, \code{numeric()} otherwise.
#'
#' @keywords internal
pd_feature_grid = function(x, n_grid) {
  if (is.factor(x)) return(levels(droplevels(x)))
  if (is.character(x)) {
    u = unique(x[!is.na(x)])
    if (!length(u)) {
      cli::cli_abort("Cannot build PD grid: no non-missing values in {.arg x}.")
    }
    return(sort(u))
  }
  if (!any(is.finite(x))) {
    cli::cli_abort("Cannot build PD grid: no finite values in {.arg x}.")
  }
  probs = seq(0, 1, length.out = as.integer(n_grid))
  g = sort(unique(as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 7))))
  g = g[is.finite(g)]
  if (length(g) < 1L) {
    cli::cli_abort("Cannot build PD grid: no finite quantiles after summarizing {.arg x}.")
  }
  g
}


#' Pack ICE Matrix into Long-Format data.table
#'
#' Converts an \code{n_obs x length(grid)} ICE matrix into a long-format
#' data.table with one row per (observation, grid value) pair.
#'
#' @param ice (`matrix`) \cr
#'   ICE predictions; shape \code{n_obs x length(grid)}.
#' @param feature (`character(1)`) \cr
#'   Name of the focal feature; stored in the \code{.feature} column.
#' @param grid (`atomic vector`) \cr
#'   Grid values used for this feature; stored in the \code{.borders} column.
#'
#' @return (`data.table`) \cr
#'   Columns: \code{.id} (observation index), \code{.type} (\code{"ice"}),
#'   \code{.feature}, \code{.borders}, \code{.value} (prediction).
#'
#' @keywords internal
pd_pack_ice_result = function(ice, feature, grid) {
  n_obs = nrow(ice)
  data.table::data.table(
    .id      = rep(seq_len(n_obs), each = length(grid)),
    .type    = "ice",
    .feature = feature,
    .borders = rep(grid, times = n_obs),
    .value   = as.vector(t(ice))
  )
}


#' Generate Predictions for New Data
#'
#' Calls \code{predict_fun} (or the default predict method) on \code{newdata}
#' and extracts a numeric prediction vector via \code{extract_numeric_prediction}.
#'
#' @param model (`any`) \cr
#'   Fitted model.
#' @param newdata (`data.frame()` or `data.table()`) \cr
#'   New observations to predict.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   \code{function(model, data)} returning predictions; \code{NULL} = default.
#'
#' @return (`numeric()`) \cr
#'   Numeric prediction vector of length \code{nrow(newdata)}.
#'
#' @keywords internal
pd_predict = function(model, newdata, predict_fun = NULL) {
  fun = if (is.null(predict_fun)) pd_select_predict_fun(model) else predict_fun
  pred_raw = fun(model, newdata)
  extract_numeric_prediction(pred_raw, expected_n = nrow(newdata))
}


#' Select Predict Function for a Model
#'
#' Returns a two-argument \code{function(model, data)} suitable for prediction:
#' if \code{model} is itself a function, wraps it directly; otherwise delegates
#' to \code{default_predict_fun}.
#'
#' @param model (`any`) \cr
#'   Fitted model or bare prediction function.
#'
#' @return (`function`) \cr
#'   \code{function(model, data)} returning raw predictions.
#'
#' @keywords internal
pd_select_predict_fun = function(model) {
  if (is.function(model)) {
    return(function(model, data) model(data))
  }
  function(model, data) default_predict_fun(model, data)
}
