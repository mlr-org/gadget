skip_cpp_if_unavailable = function() {
  tryCatch({
    search_best_split_cpp(Z = data.frame(x = 1:5), Y = list(matrix(1:10, ncol = 2)), min_node_size = 2)
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("C++ symbols not loaded (install package with compile)")
    }
  })
}

test_that("search_best_split_cpp works with numeric data", {
  skip_cpp_if_unavailable()
  set.seed(1)
  n = 25
  Z = data.frame(x = runif(n), y = runif(n))
  Y = list(matrix(rnorm(n * 2), ncol = 2))
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 5)
  expect_true(is.data.frame(result))
  expect_true(nrow(result) >= 1)
  expect_true(all(c("split_feature", "is_categorical", "split_point",
        "split_objective", "best_split") %in% names(result)))
  expect_true(all(result$is_categorical %in% c(TRUE, FALSE)))
})

test_that("search_best_split_cpp works with single feature", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:10)
  Y = list(matrix(rnorm(20), ncol = 2))
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 2)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$split_feature[1], "x")
})

test_that("search_best_split_cpp respects min_node_size", {
  skip_cpp_if_unavailable()
  Z = data.frame(x = 1:10)
  Y = list(matrix(rnorm(20), ncol = 2))
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 10)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_false(result$best_split[1])
})

test_that("search_best_split_cpp with multiple Y matrices", {
  skip_cpp_if_unavailable()
  n = 20
  Z = data.frame(x = runif(n))
  Y = list(
    matrix(rnorm(n * 2), ncol = 2),
    matrix(rnorm(n * 2), ncol = 2)
  )
  result = search_best_split_cpp(Z = Z, Y = Y, min_node_size = 5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})
