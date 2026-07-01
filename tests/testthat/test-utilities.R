test_that("convert_tree_to_list returns depth-based list", {
  skip_if_not(exists("convert_tree_to_list", envir = asNamespace("xplaineff")),
    "convert_tree_to_list not available")
  # Build a minimal root node (no children)
  grid = list(x = 1:5)
  root = xplaineff:::Node$new(
    id = 1, depth = 1, subset_idx = 1:10, grid = grid,
    id_parent = NULL, child_type = NULL
  )
  tree_list = xplaineff:::convert_tree_to_list(root, max_depth = 2)
  expect_true(is.list(tree_list))
  expect_true(length(tree_list) >= 1)
  expect_true(is.list(tree_list[[1]]))
  expect_equal(length(tree_list[[1]]), 1)
  expect_equal(tree_list[[1]][[1]]$id, 1)
})

test_that("extract_split_info works with depth-based tree list", {
  # Tree: list of depths, each depth is list of nodes (new Node structure: split, objective, importance, parent)
  node1 = list(
    id = 1, depth = 1, subset_idx = 1:20,
    split = list(feature = "x", value = 0.5),
    objective = list(value = 1.5, value_j = NULL),
    importance = list(imp = NA, imp_j = NULL),
    parent = NULL,
    improvement_met = FALSE, stop_criterion_met = FALSE,
    children = list(left_child = "dummy", right_child = "dummy")
  )
  node2 = list(
    id = 2, depth = 2, subset_idx = 1:10,
    split = NULL,
    objective = list(value = 0.5, value_j = NULL),
    importance = list(imp = 0.2, imp_j = c(x = 0.2)),
    parent = list(split_feature = "x", split_value = 0.5, objective_value = 1.5, int_imp = NA),
    improvement_met = FALSE, stop_criterion_met = TRUE,
    children = NULL
  )
  node3 = list(
    id = 3, depth = 2, subset_idx = 11:20,
    split = NULL,
    objective = list(value = 0.6, value_j = NULL),
    importance = list(imp = 0.3, imp_j = c(x = 0.3)),
    parent = list(split_feature = "x", split_value = 0.5, objective_value = 1.5, int_imp = NA),
    improvement_met = FALSE, stop_criterion_met = TRUE,
    children = NULL
  )
  tree = list(list(node1), list(node2, node3))
  result = xplaineff:::extract_split_info(tree)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
  expect_true(all(c("depth", "id", "split_feature", "n_obs") %in% names(result)))
  expect_false(any(c("split_levels_left", "split_levels_right", "split_levels_parent") %in% names(result)))
})

test_that("find_node_by_id finds node by id in flat node list", {
  node1 = list(id = 1, depth = 1)
  node2 = list(id = 2, id_parent = 1, depth = 2)
  node3 = list(id = 3, id_parent = 1, depth = 2)
  node_list = list(node1, node2, node3)
  expect_equal(xplaineff:::find_node_by_id(node_list, 2)$id, 2)
  expect_equal(xplaineff:::find_node_by_id(node_list, 1)$id, 1)
  expect_null(xplaineff:::find_node_by_id(node_list, 99))
})

test_that("node_heterogeneity returns non-negative numeric vector", {
  Y = list(
    matrix(rnorm(20), ncol = 2),
    matrix(rnorm(20), ncol = 2)
  )
  result = xplaineff:::node_heterogeneity(Y)
  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(all(!is.na(result)))
  expect_true(all(result >= 0))
})

test_that("order_categorical_levels returns factor with ordered levels", {
  data = data.frame(
    cat = factor(rep(letters[1:3], each = 5)),
    x = rnorm(15),
    y = rnorm(15)
  )
  x_cat = droplevels(data$cat)
  result = xplaineff:::order_categorical_levels(
    x_cat, data, feature = "cat", target_feature_name = "y", order_method = "raw"
  )
  expect_true(is.factor(result))
  expect_equal(length(result), length(x_cat))
  expect_true(all(levels(result) %in% levels(x_cat)))
  result_mds = xplaineff:::order_categorical_levels(
    x_cat, data, feature = "cat", target_feature_name = "y", order_method = "mds"
  )
  expect_true(is.factor(result_mds))
  expect_equal(length(levels(result_mds)), nlevels(x_cat))
})

test_that("order_categorical_levels uses half-L1 for categorical auxiliary features", {
  make_level = function(level, numeric_ones, categorical_v) {
    data.frame(
      cat = rep(level, 10L),
      x_num = c(rep(1, numeric_ones), rep(0, 10L - numeric_ones)),
      x_cat = factor(c(rep("v", categorical_v), rep("u", 10L - categorical_v)), levels = c("u", "v")),
      y = 0,
      stringsAsFactors = FALSE
    )
  }
  data = data.frame(
    rbind(
      make_level("a", numeric_ones = 0L, categorical_v = 0L),
      make_level("b", numeric_ones = 1L, categorical_v = 2L),
      make_level("c", numeric_ones = 2L, categorical_v = 1L)
    )
  )
  data$cat = factor(data$cat, levels = c("a", "b", "c"))
  x_cat = droplevels(data$cat)
  ord = xplaineff:::order_categorical_levels(
    x_cat, data, feature = "cat", target_feature_name = "y", order_method = "mds"
  )
  expect_equal(levels(ord)[2L], "b")
})

test_that("ALE categorical split helpers use ordered-prefix semantics", {
  x = factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"), ordered = TRUE)
  mask = xplaineff:::ordered_categorical_left_mask(x, "b")
  expect_equal(mask, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))

  strategy = list(name = "ale")
  node = xplaineff:::Node$new(
    id = 1L, depth = 1L, subset_idx = seq_along(x), grid = list(cat = x), strategy = strategy
  )
  grids = node$create_child_grids("cat", "b", is_categorical = TRUE)
  expect_equal(as.character(grids$grid_left$cat), c("a", "a", "b", "b"))
  expect_equal(as.character(grids$grid_right$cat), c("c", "c"))
})

test_that("PD categorical child conditions display category sets", {
  x = factor(c("a", "b", "c", "b"), levels = c("a", "b", "c"))
  strategy = list(
    name = "pd",
    get_child_objectives = function(...) {
      list(
        left_objective_value_j = c(cat = 1),
        right_objective_value_j = c(cat = 1),
        left_objective_value = 1,
        right_objective_value = 1
      )
    }
  )
  node = xplaineff:::Node$new(
    id = 1L, depth = 1L, subset_idx = seq_along(x), grid = list(cat = x),
    objective_value = 10, objective_value_j = c(cat = 10), strategy = strategy
  )

  children = node$create_children(
    z_split_feature = x,
    Y = list(),
    split_info = list(split_feature = "cat", split_value = "b", is_categorical = TRUE),
    objective_value_root_j = c(cat = 10),
    objective_value_root = 10,
    impr_par = 0
  )

  expect_equal(children$left_child$parent$split_condition, "cat in {b}")
  expect_equal(children$right_child$parent$split_condition, "cat in {a, c}")
})

test_that("mean_center_ice returns Y and grid from effect with results", {
  set.seed(1)
  effect = list(results = data.frame(
    x = rep(1:3, each = 2),
    .value = rnorm(6),
    .type = "ice",
    .id = rep(1:2, 3)
  ))
  result = xplaineff:::mean_center_ice(effect, feature_set = NULL, mean_center = TRUE)
  expect_true(is.list(result))
  expect_true("Y" %in% names(result))
  expect_true("grid" %in% names(result))
  expect_true(is.list(result$Y))
  expect_true(is.list(result$grid))
})

test_that("prepare_split_data_pd separates effect features from split features", {
  effect = list(results = list(
    x1 = data.frame(
      .id = rep(1:3, each = 2),
      .type = "ice",
      .feature = "x1",
      .borders = rep(c(0, 1), times = 3),
      .value = seq_len(6)
    )
  ))
  data = data.frame(x1 = 1:3, x2 = 4:6, y = 7:9)

  prepared = xplaineff:::prepare_split_data_pd(
    effect = effect,
    data = data,
    target_feature_name = "y",
    feature_set = "x1",
    split_feature = "x2"
  )
  expect_equal(names(prepared$Y), "x1")
  expect_equal(names(prepared$Z), "x2")

  prepared_all_splits = xplaineff:::prepare_split_data_pd(
    effect = effect,
    data = data,
    target_feature_name = "y",
    feature_set = "x1"
  )
  expect_equal(names(prepared_all_splits$Y), "x1")
  expect_equal(names(prepared_all_splits$Z), c("x1", "x2"))
})
