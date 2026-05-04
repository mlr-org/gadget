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


pd_predict = function(model, newdata, predict_fun = NULL) {
  fun = if (is.null(predict_fun)) pd_select_predict_fun(model) else predict_fun
  pred_raw = fun(model, newdata)
  extract_numeric_prediction(pred_raw, expected_n = nrow(newdata))
}


pd_select_predict_fun = function(model) {
  if (is.function(model)) {
    return(function(model, data) model(data))
  }
  function(model, data) default_predict_fun(model, data)
}
