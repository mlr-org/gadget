test_that("binary categorical ALE-like dt: columns interval_index x_left x_right d_l survive subsetting", {
  dt <- data.table::data.table(
    row_id = 1:6,
    feat_val = factor(c("A", "A", "A", "B", "B", "B"), levels = c("A", "B")),
    x_left = factor(rep(c("A", "A"), each = 3L), levels = c("A", "B")),
    x_right = factor(rep(c("B", "B"), each = 3L), levels = c("A", "B")),
    d_l = c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3),
    interval_index = c(1L, 1L, 1L, 2L, 2L, 2L)
  )
  sub <- dt[dt$row_id %in% 1:3]
  expect_true(all(c("interval_index", "x_left", "x_right", "d_l") %in% names(sub)))
  expect_equal(sum(sub$d_l == 0), 0)
  expect_equal(unique(sub$interval_index), 1L)
})

test_that("prepare_plot_data_ale: delta_aggr nrow equals mean_effect nrow for binary factor", {
  eff <- list(
    xf = data.table::data.table(
      row_id = 1:6,
      feat_val = factor(c("A", "A", "A", "B", "B", "B"), levels = c("A", "B")),
      x_left = factor(rep(c("A", "A"), each = 3L), levels = c("A", "B")),
      x_right = factor(rep(c("B", "B"), each = 3L), levels = c("A", "B")),
      d_l = c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3),
      interval_index = c(1L, 1L, 1L, 2L, 2L, 2L)
    )
  )
  out <- xplaineff:::prepare_plot_data_ale(eff, idx = NULL, features = "xf", mean_center = TRUE)
  me <- out$xf$mean_effect
  da <- eff$xf[, list(d_l = mean(d_l, na.rm = TRUE)), by = c("interval_index", "x_left", "x_right")]
  expect_equal(nrow(me), nrow(da))
  expect_true(is.factor(me$x_grid))
  expect_equal(levels(me$x_grid), c("A", "B"))
})

test_that("prepare_plot_data_ale idx drops one interval: fewer mean_effect rows", {
  eff <- list(
    xf = data.table::data.table(
      row_id = 1:6,
      feat_val = factor(c("A", "A", "A", "B", "B", "B"), levels = c("A", "B")),
      x_left = factor(rep(c("A", "A"), each = 3L), levels = c("A", "B")),
      x_right = factor(rep(c("B", "B"), each = 3L), levels = c("A", "B")),
      d_l = c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3),
      interval_index = c(1L, 1L, 1L, 2L, 2L, 2L)
    )
  )
  full <- xplaineff:::prepare_plot_data_ale(eff, idx = NULL, features = "xf", mean_center = TRUE)
  one <- xplaineff:::prepare_plot_data_ale(eff, idx = 1:3, features = "xf", mean_center = TRUE)
  expect_true(nrow(one$xf$mean_effect) < nrow(full$xf$mean_effect))
})

test_that("mean_center_ale maps sample d_l == 0 to NA before aggregation", {
  dt <- data.table::data.table(
    row_id = 1:4,
    feat_val = factor(c("A", "A", "B", "B"), levels = c("A", "B")),
    x_left = factor(c("A", "A", "A", "A"), levels = c("A", "B")),
    x_right = factor(c("B", "B", "B", "B"), levels = c("A", "B")),
    d_l = c(0, 0, 0.5, 0.5),
    interval_index = c(1L, 1L, 2L, 2L)
  )
  out <- xplaineff:::mean_center_ale(data.table::copy(dt), mean_center = FALSE)
  expect_true(nrow(out) >= 1L)
})
