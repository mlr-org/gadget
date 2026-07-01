# Basic behavioral contract tests (no brittle symbol/existence checks)

test_that("strategies satisfy contract: name and GadgetTree creation", {
  library(xplaineff)
  pd = PdStrategy$new()
  ale = AleStrategy$new()
  expect_equal(pd$name, "pd")
  expect_equal(ale$name, "ale")
  expect_s3_class(GadgetTree$new(strategy = pd), "GadgetTree")
  expect_s3_class(GadgetTree$new(strategy = ale), "GadgetTree")
})

test_that("GadgetTree$new rejects non-R6 strategy", {
  library(xplaineff)
  expect_error(
    GadgetTree$new(strategy = "not_r6", n_split = 2),
    "R6"
  )
})

test_that("GadgetTree$fit rejects target_feature_name not in data", {
  library(xplaineff)
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 1)
  data = data.frame(x = 1:5, y = 1:5)
  expect_error(
    tree$fit(data = data, target_feature_name = "z", effect = list()),
    "subset|not.*found|target_feature_name"
  )
})

test_that("PdStrategy fit rejects missing effect", {
  library(xplaineff)
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = 1)
  data = data.frame(x = 1:5, y = 1:5)
  expect_error(
    tree$fit(data = data, target_feature_name = "y"),
    "effect|PdStrategy requires"
  )
})

test_that("AleStrategy fit rejects missing model", {
  library(xplaineff)
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = 1)
  data = data.frame(x = 1:5, y = 1:5)
  expect_error(
    tree$fit(data = data, target_feature_name = "y"),
    "model|AleStrategy requires"
  )
})
