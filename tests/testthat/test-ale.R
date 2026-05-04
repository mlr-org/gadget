test_that("calculate_ale returns named list of data.tables", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  set.seed(1)
  n = 40
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  result = gadget:::calculate_ale(
    model = learner,
    data = data,
    feature_set = c("x1", "x2"),
    target_feature_name = "y",
    n_intervals = 5
  )
  expect_true(is.list(result))
  expect_equal(sort(names(result)), c("x1", "x2"))
  for (nm in names(result)) {
    expect_true(data.table::is.data.table(result[[nm]]))
    expect_true(all(c("row_id", "d_l", "interval_index") %in% names(result[[nm]])))
  }
})

test_that("calculate_ale_heterogeneity_cpp returns numeric", {
  skip_ale_cpp_if_unavailable()
  n = 20
  dt = data.table::data.table(
    row_id = seq_len(n),
    interval_index = rep(1:4, length.out = n),
    d_l = rnorm(n),
    int_n = 5L, int_s1 = 0, int_s2 = 1
  )
  result = gadget:::calculate_ale_heterogeneity_single_cpp(dt$d_l, dt$interval_index)
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(!is.na(result))
  expect_true(result >= 0)
})

test_that("calculate_ale_heterogeneity_list_cpp works with list of ALE data", {
  skip_ale_cpp_if_unavailable()
  n = 15
  dt1 = data.table::data.table(
    row_id = seq_len(n), interval_index = rep(1:3, length.out = n),
    d_l = rnorm(n), int_n = 5L, int_s1 = 0, int_s2 = 1
  )
  dt2 = data.table::data.table(
    row_id = seq_len(n), interval_index = rep(1:5, length.out = n),
    d_l = rnorm(n), int_n = 3L, int_s1 = 0, int_s2 = 1
  )
  Y = list(f1 = dt1, f2 = dt2)
  result = gadget:::calculate_ale_heterogeneity_list_cpp(Y)
  # C++ returns a named list of scalars, not a numeric vector
  expect_true(is.list(result))
  expect_length(result, 2)
  expect_true(all(vapply(result, is.numeric, logical(1))))
  expect_true(all(!is.na(unlist(result))))
  expect_true(all(unlist(result) >= 0))
})

test_that("ALE tree fit stores root and cached effect", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(2)
  n = 50
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1, min_node_size = 15)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  expect_true(!is.null(tree$root))
  strat = tree$strategy
  expect_true(!is.null(strat$effect))
  expect_true(is.list(strat$effect))
})

test_that("AleStrategy defaults to cpp engine when ale_engine is omitted", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(13)
  n = 45
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t_default_engine", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1, min_node_size = 10)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  expect_identical(tree$strategy$ale_engine, "cpp")
})

test_that("prepare_split_data_ale returns Z and Y", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(3)
  n = 30
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  result = gadget:::prepare_split_data_ale(
    model = learner,
    data = data,
    target_feature_name = "y",
    n_intervals = 5,
    split_feature = c("x1", "x2")
  )
  expect_true("Z" %in% names(result))
  expect_true("Y" %in% names(result))
  expect_true(data.table::is.data.table(result$Z))
  expect_true(is.list(result$Y))
})

test_that("prepare_split_data_ale r and cpp engines are numerically aligned", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_ale_cpp_if_unavailable()
  set.seed(11)
  n = 60
  data = data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    y = rnorm(n)
  )
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)

  out_r = gadget:::prepare_split_data_ale(
    model = learner,
    data = data,
    target_feature_name = "y",
    n_intervals = 6,
    feature_set = c("x1", "x2"),
    split_feature = c("x1", "x2"),
    ale_engine = "r"
  )
  out_cpp = gadget:::prepare_split_data_ale(
    model = learner,
    data = data,
    target_feature_name = "y",
    n_intervals = 6,
    feature_set = c("x1", "x2"),
    split_feature = c("x1", "x2"),
    ale_engine = "cpp"
  )

  for (feat in c("x1", "x2")) {
    expect_equal(out_r$Y[[feat]]$d_l, out_cpp$Y[[feat]]$d_l, tolerance = 1e-10)
    expect_equal(out_r$Y[[feat]]$interval_index, out_cpp$Y[[feat]]$interval_index)
  }
})

test_that("prepare_split_data_ale cpp engine works for non-mlr3 model", {
  skip_if_not_installed("ranger")
  skip_ale_cpp_if_unavailable()
  set.seed(12)
  n = 80
  data = data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    y = rnorm(n)
  )
  model = ranger::ranger(y ~ ., data = data, num.trees = 50L, num.threads = 1L, seed = 12L)
  out = gadget:::prepare_split_data_ale(
    model = model,
    data = data,
    target_feature_name = "y",
    n_intervals = 5,
    feature_set = c("x1", "x2"),
    split_feature = c("x1", "x2"),
    ale_engine = "cpp"
  )
  expect_true(data.table::is.data.table(out$Z))
  expect_true(is.list(out$Y))
  expect_true(all(c("x1", "x2") %in% names(out$Y)))
})

test_that("calculate_ale returns zero d_l for constant numeric feature", {
  n = 15L
  data = data.frame(x_const = rep(2.5, n), y = seq_len(n))
  predict_fun = function(model, data) rep(0.0, nrow(data))
  result = gadget:::calculate_ale(
    model = NULL, data = data,
    feature_set = "x_const", target_feature_name = "y",
    predict_fun = predict_fun
  )
  expect_true(all(result$x_const$d_l == 0))
  expect_equal(nrow(result$x_const), n)
})

test_that("calculate_ale returns zero d_l for single-level factor feature", {
  n = 12L
  data = data.frame(x_cat = factor(rep("a", n)), y = seq_len(n))
  predict_fun = function(model, data) rep(0.0, nrow(data))
  result = gadget:::calculate_ale(
    model = NULL, data = data,
    feature_set = "x_cat", target_feature_name = "y",
    predict_fun = predict_fun
  )
  expect_true(all(result$x_cat$d_l == 0))
  expect_equal(nrow(result$x_cat), n)
})

test_that("cpp_ale_numeric_prepare returns zero_effect for constant feature", {
  skip_ale_cpp_if_unavailable()
  result = gadget:::cpp_ale_numeric_prepare(rep(1.0, 20L), n_intervals = 5L)
  expect_true(isTRUE(result$zero_effect))
})

test_that("ale_sweep_cpp returns valid split for simple one-feature case", {
  skip_ale_cpp_if_unavailable()
  n = 6L
  # d_l: all 1 in interval 1 (obs 1-3), all -1 in interval 2 (obs 4-6)
  d_l_mat = matrix(c(1.0, 1.0, 1.0, -1.0, -1.0, -1.0), nrow = 1L)
  interval_idx_mat = matrix(c(1L, 1L, 1L, 2L, 2L, 2L), nrow = 1L)
  result = gadget:::ale_sweep_cpp(
    ord_idx = 1:6,
    d_l_mat = d_l_mat, interval_idx_mat = interval_idx_mat,
    offsets = 0L,
    tot_n = c(3.0, 3.0), tot_s1 = c(3.0, -3.0), tot_s2 = c(3.0, 3.0),
    r_risks = c(0.0, 0.0),
    is_cand = rep(TRUE, 6), min_node_size = 2L,
    split_feat_j = 0L, use_stabilizer = FALSE,
    z_sorted = as.numeric(1:6), n_obs = n
  )
  expect_true(is.list(result))
  expect_true(all(c("best_t", "best_risks_sum") %in% names(result)))
  expect_true(result$best_t >= 2L && result$best_t <= 4L)
  expect_true(is.finite(result$best_risks_sum))
})

test_that("ale_sweep_cpp returns Inf when no candidate splits", {
  skip_ale_cpp_if_unavailable()
  n = 6L
  d_l_mat = matrix(c(1.0, 1.0, 1.0, -1.0, -1.0, -1.0), nrow = 1L)
  interval_idx_mat = matrix(c(1L, 1L, 1L, 2L, 2L, 2L), nrow = 1L)
  result = gadget:::ale_sweep_cpp(
    ord_idx = 1:6,
    d_l_mat = d_l_mat, interval_idx_mat = interval_idx_mat,
    offsets = 0L,
    tot_n = c(3.0, 3.0), tot_s1 = c(3.0, -3.0), tot_s2 = c(3.0, 3.0),
    r_risks = c(0.0, 0.0),
    is_cand = rep(FALSE, 6), min_node_size = 2L,
    split_feat_j = 0L, use_stabilizer = FALSE,
    z_sorted = as.numeric(1:6), n_obs = n
  )
  expect_true(is.na(result$best_t))
  expect_true(is.infinite(result$best_risks_sum) && result$best_risks_sum > 0)
})
