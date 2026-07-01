test_that("PdStrategy can be created", {
  strategy = PdStrategy$new()
  expect_true(inherits(strategy, "PdStrategy"))
  expect_equal(strategy$name, "pd")
})

test_that("PdStrategy find_best_split returns expected structure", {
  tryCatch({
    xplaineff:::search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
  set.seed(1)
  n = 30
  Z = data.frame(x = runif(n), y = runif(n))
  Y = list(
    matrix(rnorm(n * 2), ncol = 2),
    matrix(rnorm(n * 2), ncol = 2)
  )
  strategy = PdStrategy$new()
  res = strategy$find_best_split(Z = Z, Y = Y, min_node_size = 5, n_quantiles = NULL)
  expect_true(is.data.frame(res))
  expect_true(nrow(res) >= 1)
  expect_true(all(c("split_feature", "is_categorical", "split_point",
        "split_objective", "best_split") %in% names(res)))
})

test_that("PdStrategy heterogeneity returns numeric vector", {
  Y = list(matrix(rnorm(20), ncol = 2), matrix(rnorm(20), ncol = 2))
  strategy = PdStrategy$new()
  h = strategy$heterogeneity(Y)
  expect_true(is.numeric(h))
  expect_length(h, 2)
  expect_true(all(!is.na(h)))
  expect_true(all(h >= 0))
})

test_that("AleStrategy can be created", {
  strategy = AleStrategy$new()
  expect_true(inherits(strategy, "AleStrategy"))
  expect_equal(strategy$name, "ale")
})

test_that("AleStrategy heterogeneity returns numeric for ALE-like list", {
  tryCatch({
    dt = data.table::data.table(row_id = 1:5, interval_index = rep(1L, 5), d_l = 0, int_n = 5L, int_s1 = 0, int_s2 = 0)
    xplaineff:::calculate_ale_heterogeneity_list_cpp(list(x = dt))
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("ALE C++ symbols not loaded (install package with compile)")
    }
  })
  n = 20
  dt = data.table::data.table(
    row_id = seq_len(n),
    interval_index = rep(1:4, length.out = n),
    d_l = rnorm(n),
    int_n = 5L, int_s1 = 0, int_s2 = 1
  )
  Y = list(f1 = dt)
  strategy = AleStrategy$new()
  h = strategy$heterogeneity(Y)
  expect_true(is.numeric(h))
  expect_length(h, 1)
  expect_true(!is.na(h))
  expect_true(h >= 0)
})

test_that("AleStrategy child objectives use the selected split when best splits tie", {
  strategy = AleStrategy$new()
  raw = data.frame(
    split_feature = rep(c("x", "z"), each = 2L),
    split_point = rep(c(0.5, 1.5), each = 2L),
    feature = rep(c("f1", "f2"), times = 2L),
    left_objective_value_j = c(1, 2, 100, 200),
    right_objective_value_j = c(3, 4, 300, 400),
    best_split = TRUE
  )
  split_info = list(
    split_feature = "x",
    split_value = 0.5,
    raw_result = raw
  )

  obj = strategy$get_child_objectives(
    Z = NULL, Y = NULL, split_info = split_info,
    idx_left = NULL, idx_right = NULL, grid_left = NULL, grid_right = NULL
  )

  expect_equal(obj$left_objective_value_j, c(f1 = 1, f2 = 2))
  expect_equal(obj$right_objective_value_j, c(f1 = 3, f2 = 4))
  expect_equal(obj$left_objective_value, 3)
  expect_equal(obj$right_objective_value, 7)
})
