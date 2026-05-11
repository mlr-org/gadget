#!/usr/bin/env Rscript
# GADGET efficiency benchmark: Global PDP, Global ALE, Regional PDP, Regional ALE.
# Two model types: RF (ranger) and toy (analytic DGP).
# Run: Rscript simulation/benchmark_gadget.R [--datadir DIR] [--outdir DIR] [--reps N] [--predict-reps N]
#                                         [--fail-fast true|false|1|0|yes|no|...]

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1) != "Package: gadget") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1) == "Package: gadget") {
    setwd("..")
  } else {
    stop("Run from GADGET package root")
  }
}

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else if (!requireNamespace("gadget", quietly = TRUE)) {
  stop("Install gadget or devtools")
} else {
  library(gadget)
}
library(data.table)
library(ranger)

# Explicit RF configuration used for parity with simulation/benchmark_effector.py.
# Keep this list and Python RF_CONFIG synchronized.
# Notes on non-1:1 mapping:
# - sklearn `min_samples_split` has no exact ranger equivalent (ranger splits are
#   still constrained by `min.node.size` on child nodes).
# - sklearn `min_impurity_decrease` / `ccp_alpha` have no direct ranger analogue.
# - sklearn `min_weight_fraction_leaf` is not used (uniform sample weights here).
rf_config <- list(
  num_trees = 100L,
  mtry_mode = "all_features", # mapped to sklearn max_features = 1.0
  min_node_size = 1L,
  max_depth = NULL, # mapped to sklearn max_depth = None
  replace = TRUE,
  sample_fraction = 1.0,
  splitrule = "variance",
  respect_unordered_factors = "ignore", # numeric-only benchmark data
  # No exact ranger counterpart for sklearn min_samples_split / min_impurity_decrease / ccp_alpha.
  num_threads = 1L,
  seed = 21L
)

args <- commandArgs(trailingOnly = TRUE)
datadir <- "simulation/data/benchmark"
outdir <- "simulation/results/benchmark"
reps <- 5L
predict_reps <- 20L
fixed_N <- 1000L
fixed_D <- 10L
fail_fast <- TRUE

parse_fail_fast_flag <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("true", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "0", "no", "n", "off")) return(FALSE)
  warning("Unrecognized --fail-fast value '", x, "', defaulting to TRUE", immediate. = TRUE)
  TRUE
}

i <- 1
while (i <= length(args)) {
  if (args[i] == "--datadir" && i < length(args)) {
    datadir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--reps" && i < length(args)) {
    reps <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--predict-reps" && i < length(args)) {
    predict_reps <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--fixed-N" && i < length(args)) {
    fixed_N <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--fixed-D" && i < length(args)) {
    fixed_D <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--fail-fast" && i < length(args)) {
    fail_fast <- parse_fail_fast_flag(args[i + 1]); i <- i + 2
  } else {
    i <- i + 1
  }
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Global PDP/ALE use the C fast path only for mlr3 Learners with predict_newdata_fast.
# Plain ranger + custom predict_fun never triggers it; effector compares against sklearn RF.
use_mlr3_rf_for_benchmark <- requireNamespace("mlr3", quietly = TRUE) &&
  requireNamespace("mlr3learners", quietly = TRUE)
if (use_mlr3_rf_for_benchmark) {
  suppressPackageStartupMessages(loadNamespace("mlr3learners"))
  message("RF benchmark: mlr3 regr.ranger (GADGET global PDP/ALE fast path enabled).")
} else {
  message(
    "RF benchmark: ranger + custom predict (GADGET global fast path off; ",
    "install Suggests mlr3 + mlr3learners for comparable timing)."
  )
}

# Defaults match the **small** preset in efficiency_benchmark_plan.md §4.5 / run_benchmark.sh.
# Override via --N-vec / --D-vec / --fixed-N / --fixed-D / --n-grid-vec / --n-int-vec.
N_vec <- c(500L, 1000L, 5000L)
D_vec <- c(5L, 10L, 20L)
n_grid_vec <- c(10L, 20L, 50L)
n_int_vec <- c(10L, 20L, 50L)

parse_int_vec <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
i <- 1
while (i <= length(args)) {
  if (args[i] == "--N-vec" && i < length(args)) {
    N_vec <- parse_int_vec(args[i + 1L]); i <- i + 2L
  } else if (args[i] == "--D-vec" && i < length(args)) {
    D_vec <- parse_int_vec(args[i + 1L]); i <- i + 2L
  } else if (args[i] == "--n-grid-vec" && i < length(args)) {
    n_grid_vec <- parse_int_vec(args[i + 1L]); i <- i + 2L
  } else if (args[i] == "--n-int-vec" && i < length(args)) {
    n_int_vec <- parse_int_vec(args[i + 1L]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

load_data <- function(N, D) {
  f <- file.path(datadir, sprintf("benchmark_N%d_D%d_seed21.csv", N, D))
  if (!file.exists(f)) return(NULL)
  as.data.frame(fread(f))
}

rf_pred_fun <- function(model, newdata) predict(model, newdata)$predictions

fit_rf_gadget_benchmark <- function(dat) {
  p <- max(1L, ncol(dat) - 1L)
  mtry_value <- if (identical(rf_config$mtry_mode, "all_features")) p else rf_config$mtry
  if (use_mlr3_rf_for_benchmark) {
    task <- mlr3::as_task_regr(dat, target = "y", id = "gadget_bench_rf")
    l <- mlr3::lrn(
      "regr.ranger",
      num.trees = rf_config$num_trees,
      mtry = mtry_value,
      min.node.size = rf_config$min_node_size,
      max.depth = rf_config$max_depth,
      replace = rf_config$replace,
      sample.fraction = rf_config$sample_fraction,
      splitrule = rf_config$splitrule,
      respect.unordered.factors = rf_config$respect_unordered_factors,
      num.threads = rf_config$num_threads,
      seed = rf_config$seed
    )
    l$train(task)
    l
  } else {
    ranger(
      y ~ .,
      data = dat,
      num.trees = rf_config$num_trees,
      mtry = mtry_value,
      min.node.size = rf_config$min_node_size,
      max.depth = rf_config$max_depth,
      replace = rf_config$replace,
      sample.fraction = rf_config$sample_fraction,
      splitrule = rf_config$splitrule,
      respect.unordered.factors = rf_config$respect_unordered_factors,
      num.threads = rf_config$num_threads,
      seed = rf_config$seed
    )
  }
}

rf_predict_baseline_batch <- function(model, X) {
  if (inherits(model, "Learner")) {
    model$predict_newdata(as.data.frame(X))
  } else {
    rf_pred_fun(model, X)
  }
}

toy_pred_fun <- function(model, newdata) {
  x1 <- newdata[["x1"]]
  x2 <- newdata[["x2"]]
  x3 <- newdata[["x3"]]
  5 * x1 + 5 * x2 + ifelse(x3 > 0, 10 * x1, 0) - ifelse(x3 > 0, 10 * x2, 0)
}

# ---------------------------------------------------------------------------
# Method runners
# ---------------------------------------------------------------------------

run_global_pdp_engine <- function(dat, model, n_grid, pred_fun, pd_engine = c("cpp", "r")) {
  pd_engine = match.arg(pd_engine)
  tic <- proc.time()
  gadget:::calculate_pd(model, dat, target_feature_name = "y",
    feature_set = NULL, predict_fun = pred_fun, n_grid = n_grid, pd_engine = pd_engine)
  (proc.time() - tic)[["elapsed"]]
}

run_global_ale_engine <- function(dat, model, n_int, pred_fun, ale_engine = c("cpp", "r")) {
  ale_engine = match.arg(ale_engine)
  features = setdiff(colnames(dat), "y")
  tic <- proc.time()
  if (identical(ale_engine, "cpp")) {
    gadget:::calculate_ale_fast(
      model = model,
      data = dat,
      feature_set = features,
      target_feature_name = "y",
      n_intervals = n_int,
      predict_fun = pred_fun
    )
  } else {
    gadget:::calculate_ale(
      model = model,
      data = dat,
      feature_set = features,
      target_feature_name = "y",
      n_intervals = n_int,
      predict_fun = pred_fun
    )
  }
  (proc.time() - tic)[["elapsed"]]
}

# Regional runners return c(split = <tree-only>, total = <global + split>) so that a single
# benchmark run captures both the partitioning cost and the full two-stage cost.
# effector's RegionalPDP/RegionalALE.fit() bundles the global ICE/ALE precompute inside
# itself; the "split" timing isolates our equivalent tree-fitting step for a direct
# algorithmic comparison, while "total" gives a fair pipeline-level comparison.

run_regional_pdp_engine <- function(dat, model, n_grid, pred_fun) {
  old_opts <- options(future.globals.maxSize = 4 * 1024 * 1024^2)
  on.exit(options(old_opts), add = TRUE)

  strat <- PdStrategy$new()
  tree <- GadgetTree$new(strategy = strat, n_split = 2, min_node_size = 50L)
  tree$fit(
    data = dat, target_feature_name = "y",
    model = model, predict_fun = pred_fun, n_grid = n_grid, pd_engine = "cpp"
  )
  t_global <- strat$fit_timing$global
  t_split <- strat$fit_timing$regional
  c(split = t_split, total = t_global + t_split)
}

run_regional_ale_engine <- function(dat, model, n_int, pred_fun) {
  strat <- AleStrategy$new()
  tree <- GadgetTree$new(strategy = strat, n_split = 2, min_node_size = 50L)
  tree$fit(data = dat, target_feature_name = "y", model = model, n_intervals = n_int,
    predict_fun = pred_fun, order_method = "raw", feature_set = NULL, ale_engine = "cpp")
  t_global <- strat$fit_timing$global
  t_split <- strat$fit_timing$regional
  c(split = t_split, total = t_global + t_split)
}

# ---------------------------------------------------------------------------
# Predict baseline (RF only)
# ---------------------------------------------------------------------------

run_predict_baseline <- function() {
  baseline <- list()
  for (N in N_vec) {
    for (D in D_vec) {
      dat <- load_data(N, D)
      if (is.null(dat)) next
      mod <- fit_rf_gadget_benchmark(dat)
      X <- dat[, setdiff(colnames(dat), "y"), drop = FALSE]
      times <- replicate(predict_reps, system.time(rf_predict_baseline_batch(mod, X))[["elapsed"]])
      baseline <- c(baseline, list(data.frame(
        package = "gadget", N = N, D = D,
        predict_time_mean = mean(times), predict_time_sd = sd(times), n_rep = predict_reps
      )))
    }
  }
  if (length(baseline) > 0) {
    df <- do.call(rbind, baseline)
    write.csv(df, file.path(outdir, "predict_baseline_gadget.csv"), row.names = FALSE)
    message("Written: predict_baseline_gadget.csv")
  }
}

# ---------------------------------------------------------------------------
# Record helper
# ---------------------------------------------------------------------------

results_rf <- list()
results_toy <- list()

record <- function(store_name, package_name, impl, method, N, D, n_grid, n_intervals, time_sec, r,
                   sub_experiment = "", status = "ok", error_message = NA_character_) {
  row <- data.frame(
    package = package_name, impl = impl, method = method,
    sub_experiment = sub_experiment,
    N = N, D = D,
    n_grid = if (grepl("pdp", method)) n_grid else NA_integer_,
    n_intervals = if (grepl("ale", method)) n_intervals else NA_integer_,
    repetition = r, time_sec = time_sec,
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  .GlobalEnv[[store_name]] <- c(.GlobalEnv[[store_name]], list(row))
}

# ---------------------------------------------------------------------------
# Run a full method sweep (vs N, vs D, vs resolution)
# ---------------------------------------------------------------------------

run_sweep <- function(store_name, package_name, impl, method, runner, model_type, fixed_N, fixed_D) {
  is_pdp <- grepl("pdp", method)
  res_vec <- if (is_pdp) n_grid_vec else n_int_vec
  default_res <- 20L

  make_model <- function(dat) {
    if (model_type == "rf") fit_rf_gadget_benchmark(dat) else "toy"
  }
  get_pred_fun <- function() {
    if (model_type == "rf") {
      if (use_mlr3_rf_for_benchmark) NULL else rf_pred_fun
    } else {
      toy_pred_fun
    }
  }

  # Handles runners that return either a scalar or a named numeric vector.
  # Named vectors (e.g. c(split=..., total=...)) are expanded to one row per element,
  # with the element name appended to the method name (e.g. "regional_pdp_split").
  warmup <- function(dat, mod, res) {
    tryCatch(
      runner(dat, mod, res, get_pred_fun()),
      error = function(e) {
        msg <- sprintf(
          "[%s] %s warmup failed for package=%s, impl=%s at N=%s, D=%s, res=%s: %s",
          model_type, method, package_name, impl, nrow(dat), ncol(dat) - 1L, res, conditionMessage(e)
        )
        if (isTRUE(fail_fast)) stop(msg, call. = FALSE)
        message("    Warmup error (ignored): ", msg)
        invisible(NULL)
      }
    )
  }

  run_once <- function(dat, mod, res) {
    tryCatch(
      list(ok = TRUE, value = runner(dat, mod, res, get_pred_fun())),
      error = function(e) {
        msg <- sprintf(
          "[%s] %s failed for package=%s, impl=%s at N=%s, D=%s, res=%s: %s",
          model_type, method, package_name, impl, nrow(dat), ncol(dat) - 1L, res, conditionMessage(e)
        )
        if (isTRUE(fail_fast)) stop(msg, call. = FALSE)
        message("    Error: ", msg)
        list(ok = FALSE, value = NA_real_, error = msg)
      }
    )
  }

  record_failure <- function(N_val, D_val, res_val, r_val, sub_exp, err_msg) {
    record(store_name, package_name, impl, method, N_val, D_val,
      if (is_pdp) res_val else NA_integer_,
      if (!is_pdp) res_val else NA_integer_,
      NA_real_, r_val, sub_exp, status = "error", error_message = err_msg)
  }

  do_record <- function(t_raw, N_val, D_val, res_val, r_val, sub_exp) {
    if (is.numeric(t_raw) && length(t_raw) > 1 && !is.null(names(t_raw))) {
      for (nm in names(t_raw)) {
        record(store_name, package_name, impl, paste0(method, "_", nm),
          N_val, D_val, res_val, res_val, t_raw[[nm]], r_val, sub_exp)
      }
    } else {
      record(store_name, package_name, impl, method, N_val, D_val,
        if (is_pdp) res_val else NA_integer_,
        if (!is_pdp) res_val else NA_integer_,
        t_raw, r_val, sub_exp)
    }
  }

  message(sprintf("  [%s] %s — vs N", model_type, method))
  for (N in N_vec) {
    dat <- load_data(N, fixed_D)
    if (is.null(dat)) next
    mod <- make_model(dat)
    warmup(dat, mod, default_res)
    for (r in seq_len(reps)) {
      once <- run_once(dat, mod, default_res)
      if (!once$ok) record_failure(N, fixed_D, default_res, r, "vs_N", once$error)
      else do_record(once$value, N, fixed_D, default_res, r, "vs_N")
    }
  }

  message(sprintf("  [%s] %s — vs D", model_type, method))
  for (D in D_vec) {
    dat <- load_data(fixed_N, D)
    if (is.null(dat)) next
    mod <- make_model(dat)
    warmup(dat, mod, default_res)
    for (r in seq_len(reps)) {
      once <- run_once(dat, mod, default_res)
      if (!once$ok) record_failure(fixed_N, D, default_res, r, "vs_D", once$error)
      else do_record(once$value, fixed_N, D, default_res, r, "vs_D")
    }
  }

  message(sprintf("  [%s] %s — vs resolution", model_type, method))
  for (rv in res_vec) {
    dat <- load_data(fixed_N, fixed_D)
    if (is.null(dat)) next
    mod <- make_model(dat)
    warmup(dat, mod, rv)
    for (r in seq_len(reps)) {
      once <- run_once(dat, mod, rv)
      if (!once$ok) record_failure(fixed_N, fixed_D, rv, r, "vs_res", once$error)
      else do_record(once$value, fixed_N, fixed_D, rv, r, "vs_res")
    }
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

message("=== Predict baseline (RF) ===")
run_predict_baseline()

methods <- list(
  list(package = "gadget", impl = "effect-cpp", name = "global_pdp",
    runner = function(dat, mod, res, pred_fun) run_global_pdp_engine(dat, mod, res, pred_fun, pd_engine = "cpp")),
  list(package = "gadget", impl = "effect-r", name = "global_pdp",
    runner = function(dat, mod, res, pred_fun) run_global_pdp_engine(dat, mod, res, pred_fun, pd_engine = "r")),
  list(package = "gadget", impl = "effect-cpp", name = "global_ale",
    runner = function(dat, mod, res, pred_fun) run_global_ale_engine(dat, mod, res, pred_fun, ale_engine = "cpp")),
  list(package = "gadget", impl = "effect-r", name = "global_ale",
    runner = function(dat, mod, res, pred_fun) run_global_ale_engine(dat, mod, res, pred_fun, ale_engine = "r")),
  # Regional runners return c(split=..., total=...) and are expanded by run_sweep into
  # regional_pdp_split / regional_pdp_total (likewise for ale).
  list(package = "gadget", impl = "cpp", name = "regional_pdp", runner = run_regional_pdp_engine),
  list(package = "gadget", impl = "cpp", name = "regional_ale", runner = run_regional_ale_engine)
)

message("=== RF model benchmarks ===")
for (m in methods) {
  run_sweep("results_rf", m$package, m$impl, m$name, m$runner, "rf", fixed_N, fixed_D)
}

message("=== Toy model benchmarks ===")
for (m in methods) {
  run_sweep("results_toy", m$package, m$impl, m$name, m$runner, "toy", fixed_N, fixed_D)
}

write_results <- function(res_list, filename) {
  if (length(res_list) == 0) {
    message("No results for ", filename)
    return()
  }
  out <- do.call(rbind, res_list)
  write.csv(out, file.path(outdir, filename), row.names = FALSE)
  message("Written: ", filename)
}

write_results(results_rf, "benchmark_gadget_rf.csv")
write_results(results_toy, "benchmark_gadget_toy.csv")
message("Done.")
