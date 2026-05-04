test_that("compute_ice cpp matches r (numeric focal feature)", {
  tryCatch(
    gadget:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(2L)
  n = 45L
  d = data.frame(x1 = runif(n), x2 = runif(n), x3 = rnorm(n))
  fit = stats::lm(x3 ~ x1 + x2, data = d)
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = seq(0.05, 0.95, length.out = 11L)
  ice_r = gadget:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = gadget:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("compute_ice cpp matches r (factor focal feature)", {
  tryCatch(
    gadget:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1L),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(3L)
  n = 40L
  d = data.frame(
    x1 = factor(sample(c("low", "high"), n, replace = TRUE)),
    x2 = runif(n),
    y = rnorm(n)
  )
  fit = stats::lm(y ~ x1 + x2, data = d)
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = levels(d$x1)
  ice_r = gadget:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = gadget:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("compute_ice cpp matches r (ranger native model, data= predict)", {
  testthat::skip_if_not_installed("ranger")
  tryCatch(
    gadget:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(4L)
  n = 50L
  d = data.frame(x1 = runif(n), x2 = runif(n), y = rnorm(n))
  fit = ranger::ranger(
    data = d,
    dependent.variable.name = "y",
    num.trees = 50L,
    num.threads = 1L,
    seed = 4L
  )
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = seq(0.1, 0.9, length.out = 9L)
  ice_r = gadget:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = gadget:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("PdStrategy fit aborts when target column is missing from data", {
  d = data.frame(x1 = 1, y = 2)
  fit = stats::lm(y ~ x1, data = d)
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 1L, min_node_size = 1L)
  testthat::expect_error(
    tree$fit(data = d, target_feature_name = "yy", model = fit),
    regexp = "target_feature_name|yy|subset"
  )
})

test_that("AleStrategy fit aborts when target column is missing from data", {
  d = data.frame(x1 = rnorm(20), y = rnorm(20))
  fit = stats::lm(y ~ x1, data = d)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1L, min_node_size = 1L)
  testthat::expect_error(
    tree$fit(data = d, target_feature_name = "yy", model = fit, n_intervals = 3L),
    regexp = "target_feature_name|yy|subset"
  )
})

test_that("pd_feature_grid aborts for all-NA numeric predictor", {
  testthat::expect_error(
    gadget:::pd_feature_grid(c(NA_real_, NA_real_), 5L),
    regexp = "grid|quantiles|finite|NA"
  )
})

test_that("compute_ice cpp matches r for integer focal feature (C++ promotes stacked column to double)", {
  tryCatch(
    gadget:::cpp_pd_stack_newdata(as.list(data.frame(x = 1)), 0L, 1.0),
    error = function(e) {
      if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
        testthat::skip("C++ pd_fast not loaded")
      }
      stop(e)
    }
  )
  set.seed(5L)
  n = 30L
  d = data.frame(x1 = sample(1L:5L, n, replace = TRUE), x2 = runif(n), y = rnorm(n))
  fit = stats::lm(y ~ x1 + x2, data = d)
  x_only = d[, c("x1", "x2"), drop = FALSE]
  grid = sort(unique(d$x1))[1:3]
  ice_r = gadget:::compute_ice_r(fit, x_only, "x1", grid, predict_fun = NULL)
  ice_cpp = gadget:::compute_ice_cpp(fit, x_only, "x1", grid, predict_fun = NULL)
  testthat::expect_equal(ice_cpp, ice_r, tolerance = 1e-10)
})

test_that("extract_numeric_prediction uses response then first prob column for mlr3 Prediction", {
  testthat::skip_if_not_installed("mlr3")
  testthat::skip_if_not_installed("mlr3learners")
  testthat::skip_if_not_installed("ranger")
  task = mlr3::tsk("iris")
  lrn = mlr3::lrn("classif.ranger", predict_type = "prob", num.trees = 30L, num.threads = 1L)
  lrn$train(task)
  pred = lrn$predict_newdata(iris[1:10, ])
  testthat::expect_equal(
    gadget:::extract_numeric_prediction(pred),
    as.numeric(pred$response)
  )
})

test_that("PdStrategy PD path can target one class prob via predict_fun", {
  testthat::skip_if_not_installed("mlr3")
  testthat::skip_if_not_installed("mlr3learners")
  testthat::skip_if_not_installed("ranger")
  set.seed(1L)
  task = mlr3::tsk("iris")
  lrn = mlr3::lrn("classif.ranger", predict_type = "prob", num.trees = 30L, num.threads = 1L)
  lrn$train(task)
  dat = iris[1:20, ]
  pf = function(model, newdata) {
    p = model$predict_newdata(newdata)$prob[, "versicolor"]
    as.numeric(p)
  }
  ice = gadget:::compute_ice_r(
    lrn, dat[, c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")],
    "Petal.Length", grid = c(1.4, 1.5), predict_fun = pf
  )
  testthat::expect_equal(nrow(ice), 20L)
  testthat::expect_equal(ncol(ice), 2L)
})

test_that("pd_feature_grid returns sorted unique levels for character x", {
  g = gadget:::pd_feature_grid(c("b", "a", "a", NA), n_grid = 5L)
  testthat::expect_equal(g, c("a", "b"))
})

test_that("calculate_y_range for PD omits raw target when mean_center is TRUE", {
  prepared_data = list(
    x1 = data.frame(`0` = 0, `1` = 0.1, node = 1L, check.names = FALSE)
  )
  dat = data.frame(y = c(100, 100))
  yr_mc = gadget:::calculate_y_range(prepared_data, dat, "y", mean_center = TRUE)
  yr_raw = gadget:::calculate_y_range(prepared_data, dat, "y", mean_center = FALSE)
  testthat::expect_true(yr_mc$ymax < 50)
  testthat::expect_true(yr_raw$ymax >= 100)
})
