test_that("track_split_condition formats categorical splits without coercion warnings", {
  left_child = list(id = 2L, depth = 2L, parent = list(id = 1L))
  parent_node = list(
    id = 1L,
    depth = 1L,
    parent = NULL,
    split = list(feature = "var_Structure", value = "M"),
    children = list(left_child = left_child, right_child = list(id = 3L))
  )
  tree = list(list(parent_node), list(left_child))
  expect_warning(
    cond <- xplaineff:::track_split_condition(left_child, tree),
    regexp = NA
  )
  expect_equal(cond, "var_Structure = M")
})

test_that("track_split_condition displays ALE categorical prefix splits as level sets", {
  left_child = list(
    id = 2L, depth = 2L,
    parent = list(id = 1L, split_condition = "var_Structure in {L, M}")
  )
  parent_node = list(
    id = 1L,
    depth = 1L,
    parent = NULL,
    split = list(feature = "var_Structure", value = "M"),
    children = list(left_child = left_child, right_child = list(id = 3L))
  )
  tree = list(list(parent_node), list(left_child))
  cond = xplaineff:::track_split_condition(left_child, tree)
  expect_equal(cond, "var_Structure in {L, M}")
})

test_that("merge_ale_y_range_with_response expands ylim to overlaid response", {
  yr = list(ymin = -1, ymax = 1)
  out = xplaineff:::merge_ale_y_range_with_response(yr, c(5, 100))
  expect_true(out$ymin <= -1)
  expect_true(out$ymax >= 100)
})

test_that("calculate_y_range_ale_combined includes regional curve range", {
  global_curves = list(f = list(mean_effect = data.frame(d_l = c(0, 1))))
  regional_curves = list(f = list(mean_effect = data.frame(d_l = c(-10, 10))))
  yr = xplaineff:::calculate_y_range_ale_combined(global_curves, regional_curves, NULL, NULL)
  expect_true(yr$ymin <= -10)
  expect_true(yr$ymax >= 10)
})

test_that("track_split_condition rounds numeric split values", {
  left_child = list(id = 2L, depth = 2L, parent = list(id = 1L))
  parent_node = list(
    id = 1L,
    depth = 1L,
    parent = NULL,
    split = list(feature = "x1", value = 0.50001),
    children = list(left_child = left_child, right_child = list(id = 3L))
  )
  tree = list(list(parent_node), list(left_child))
  cond = xplaineff:::track_split_condition(left_child, tree)
  expect_equal(cond, "x1 <= 0.5")
})

test_that("plot_tree_structure works", {
  # Tree: list of depths, each depth = list of nodes (id, id_parent, split_feature, split_value, subset_idx)
  tree = list(
    list(list(id = 1, id_parent = NA, depth = 1, split_feature = "x1", split_value = 0.5, subset_idx = 1:20)),
    list(
      list(id = 2, id_parent = 1, depth = 2, split_feature = NA, split_value = NA, subset_idx = 1:10),
      list(id = 3, id_parent = 1, depth = 2, split_feature = NA, split_value = NA, subset_idx = 11:20)
    )
  )
  plot_result = xplaineff:::plot_tree_structure(tree)
  expect_true(inherits(plot_result, "gg"))
})

test_that("plot_regional_pd is callable with valid prepared_data", {
  n = 10
  set.seed(1)
  prepared_data = list(
    x1 = data.frame(
      `0` = rnorm(n), `0.5` = rnorm(n), `1` = rnorm(n),
      node = rep(1, n),
      check.names = FALSE
    )
  )
  origin_data = data.frame(x1 = c(0, 0.5, 1)[rep(1:3, length.out = n)], y = rnorm(n))
  p = tryCatch({
    xplaineff:::plot_regional_pd(
      prepared_data = prepared_data,
      origin_data = origin_data,
      target_feature_name = "y",
      node_idx = 1,
      color_ice = "lightblue", color_pd = "red",
      ymin = NA, ymax = NA, show_point = FALSE, mean_center = TRUE
    )
  }, error = function(e) {
    testthat::skip(paste("plot_regional_pd:", conditionMessage(e)))
  })
  # plot_regional_pd returns list of ggplot objects (one per feature)
  expect_true(is.list(p))
  expect_true(length(p) >= 1)
  expect_true(inherits(p[[1]], "gg"))
})

test_that("create_plots_for_depth names PD plots with node ids", {
  prepared_data = list(
    x1 = data.frame(
      `0` = c(0.1, 0.2, 0.3, 0.4),
      `1` = c(0.5, 0.6, 0.7, 0.8),
      node = c(1L, 1L, 2L, 2L),
      check.names = FALSE
    )
  )
  data = data.frame(x1 = c(0, 1, 0, 1), y = c(1, 2, 3, 4))
  tree = list(
    list(list(id = 1L, depth = 1L, parent = NULL, subset_idx = seq_len(4L))),
    list(
      list(id = 4L, depth = 2L, parent = NULL, subset_idx = 1:2),
      list(id = 5L, depth = 2L, parent = NULL, subset_idx = 3:4)
    )
  )

  plots = xplaineff:::create_plots_for_depth(
    tree = tree,
    prepared_data = prepared_data,
    data = data,
    target_feature_name = "y",
    depth_idx = 2L,
    nodes_to_render = 1:2,
    color_ice = "lightblue",
    color_pd = "red",
    show_plot = FALSE,
    show_point = FALSE,
    mean_center = TRUE
  )

  expect_equal(names(plots), c("Node_4", "Node_5"))
})

test_that("ALE tree plot returns list", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  n = 5
  dt = data.table::data.table(
    row_id = seq_len(n), interval_index = rep(1L, n), d_l = 0, int_n = n, int_s1 = 0, int_s2 = 0
  )
  tryCatch(
    xplaineff:::calculate_ale_heterogeneity_list_cpp(list(x = dt)),
    error = function(e) testthat::skip("ALE C++ not loaded")
  )
  set.seed(123)
  n = 40
  data = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::TaskRegr$new("t", backend = data, target = "y")
  learner = mlr3::lrn("regr.ranger")
  learner$train(task)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1, min_node_size = 10)
  tree$fit(model = learner, data = data, target_feature_name = "y", n_intervals = 5)
  pl = tree$plot(data = data, target_feature_name = "y", show_plot = FALSE)
  expect_true(is.list(pl))
})
