run_ale_split_benchmark = function() {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to run this script.", call. = FALSE)
  }
  devtools::load_all(".", quiet = TRUE)

  time_call = function(fun, n_rep = 2L) {
    times = numeric(n_rep)
    for (i in seq_len(n_rep)) {
      gc()
      set.seed(1000L + i)
      times[i] = system.time(fun())[["elapsed"]]
    }
    c(mean = mean(times), min = min(times), max = max(times))
  }

  summarize_fast_root_split = function(result) {
    split_rows = unique(result[result$best_split, c("split_feature", "split_point", "split_objective")])
    split_rows[1L, , drop = FALSE]
  }

  search_best_split_ale_reference = function(Z, effect, min_node_size = 10L, n_quantiles = NULL) {
    st_table = gadget:::build_ale_interval_stats(effect, names(effect))
    split_feature_names = colnames(Z)
    per_feature = lapply(split_feature_names, function(split_feat) {
      res = gadget:::search_best_split_point_ale_with_cpp(
        z = Z[[split_feat]],
        effect = effect,
        st_table = st_table,
        split_feat = split_feat,
        is_categorical = is.factor(Z[[split_feat]]),
        n_quantiles = n_quantiles,
        min_node_size = min_node_size
      )
      data.frame(
        split_feature = split_feat,
        split_point = I(list(res$split_point)),
        split_objective = res$split_objective,
        stringsAsFactors = FALSE
      )
    })
    summary = do.call(rbind, per_feature)
    best_idx = which.min(summary$split_objective)
    list(summary = summary, best_idx = best_idx)
  }

  make_toy_self_effect = function() {
    dt = data.table::data.table(
      row_id = 1:8,
      feat_val = 1:8,
      x_left = 1:8,
      x_right = 1:8,
      d_l = c(0, 0, 0, 0, 10, 10, 10, 10),
      interval_index = rep(1L, 8)
    )
    dt[, int_n := .N, by = interval_index]
    dt[, int_s1 := sum(d_l), by = interval_index]
    dt[, int_s2 := sum(d_l^2), by = interval_index]
    list(x = dt)
  }

  make_interaction_data = function(n = 500L) {
    set.seed(1234)
    x1 = runif(n, -1, 1)
    x2 = runif(n, -1, 1)
    x3 = runif(n, -1, 1)
    data = data.frame(x1, x2, x3)
    predict_fun = function(model, newdata) {
      ifelse(newdata$x3 > 0, 3 * newdata$x1, -3 * newdata$x1) + newdata$x3
    }
    data$y = predict_fun(NULL, data)
    list(data = data, predict_fun = predict_fun)
  }

  make_example_style_data = function(n = 1000L) {
    set.seed(1)
    x1 = round(runif(n, -1, 1), 1)
    x2 = round(runif(n, -1, 1), 3)
    x3 = factor(sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5)))
    x4 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.7, 0.3))
    x5 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
    data = data.frame(x1, x2, x3, x4, x5)
    predict_fun = function(model, newdata) {
      0.2 * newdata$x1 - 8 * newdata$x2 +
        ifelse(newdata$x3 == 0, 16 * newdata$x2, 0) +
        ifelse(newdata$x1 > 0, 8 * newdata$x2, 0)
    }
    data$y = predict_fun(NULL, data)
    list(data = data, predict_fun = predict_fun)
  }

  prepare_interaction_split = function(n) {
    interaction = make_interaction_data(n)
    prepared = gadget:::prepare_split_data_ale(
      model = list(),
      data = interaction$data,
      target_feature_name = "y",
      n_intervals = 10L,
      predict_fun = interaction$predict_fun,
      ale_engine = "cpp"
    )
    st_table = gadget:::build_ale_interval_stats(prepared$Y, names(prepared$Y))
    list(prepared = prepared, st_table = st_table)
  }

  cat("== ALE split regression checks ==\n")

  toy_effect = make_toy_self_effect()
  toy_stats = gadget:::build_ale_interval_stats(toy_effect, "x")
  toy_result = gadget:::search_best_split_point_ale(
    z = 1:8,
    effect = toy_effect,
    st_table = toy_stats,
    split_feat = "x",
    is_categorical = FALSE,
    min_node_size = 2L
  )
  stopifnot(isTRUE(all.equal(toy_result$split_point, 4.5)))
  cat(sprintf("toy self-signal split point: %.1f\n", toy_result$split_point))

  interaction = make_interaction_data()
  prepared_interaction = gadget:::prepare_split_data_ale(
    model = list(),
    data = interaction$data,
    target_feature_name = "y",
    n_intervals = 10L,
    predict_fun = interaction$predict_fun,
    ale_engine = "cpp"
  )
  interaction_result = gadget:::search_best_split_ale(
    Z = prepared_interaction$Z,
    effect = prepared_interaction$Y,
    min_node_size = 10L
  )
  interaction_best = summarize_fast_root_split(interaction_result)
  stopifnot(identical(as.character(interaction_best$split_feature[[1L]]), "x3"))
  cat(sprintf("interaction DGP root split feature: %s\n", interaction_best$split_feature[[1L]]))

  example_style = make_example_style_data()
  prepared_example_style = gadget:::prepare_split_data_ale(
    model = list(),
    data = example_style$data,
    target_feature_name = "y",
    n_intervals = 10L,
    predict_fun = example_style$predict_fun,
    ale_engine = "cpp"
  )
  example_result = gadget:::search_best_split_ale(
    Z = prepared_example_style$Z,
    effect = prepared_example_style$Y,
    min_node_size = 10L
  )
  example_best = summarize_fast_root_split(example_result)
  stopifnot(identical(as.character(example_best$split_feature[[1L]]), "x3"))
  cat(sprintf("example-style DGP root split feature: %s\n", example_best$split_feature[[1L]]))

  cat("\n== ALE split timing benchmark ==\n")

  benchmark_rows = lapply(c(1000L, 2000L, 4000L), function(n) {
    prepared = prepare_interaction_split(n)
    fast_fun = function() {
      gadget:::search_best_split_point_ale(
        z = prepared$prepared$Z$x3,
        effect = prepared$prepared$Y,
        st_table = prepared$st_table,
        split_feat = "x3",
        is_categorical = FALSE,
        min_node_size = 10L
      )
    }
    slow_fun = function() {
      gadget:::search_best_split_point_ale_with_cpp(
        z = prepared$prepared$Z$x3,
        effect = prepared$prepared$Y,
        st_table = prepared$st_table,
        split_feat = "x3",
        is_categorical = FALSE,
        min_node_size = 10L
      )
    }
    fast_time = time_call(fast_fun)
    slow_time = time_call(slow_fun, n_rep = 1L)
    fast_result = fast_fun()
    slow_result = slow_fun()
    data.frame(
      n = n,
      fast_mean_sec = unname(fast_time["mean"]),
      slow_mean_sec = unname(slow_time["mean"]),
      slow_over_fast = unname(slow_time["mean"] / fast_time["mean"]),
      fast_split_point = fast_result$split_point,
      slow_split_point = slow_result$split_point,
      stringsAsFactors = FALSE
    )
  })
  benchmark_table = do.call(rbind, benchmark_rows)
  print(benchmark_table, row.names = FALSE)

  invisible(
    list(
      toy_result = toy_result,
      interaction_best = interaction_best,
      example_best = example_best,
      benchmark_table = benchmark_table,
      slow_reference = search_best_split_ale_reference(
        Z = prepared_interaction$Z,
        effect = prepared_interaction$Y,
        min_node_size = 10L
      )
    )
  )
}

if (sys.nframe() == 0L) {
  run_ale_split_benchmark()
}
