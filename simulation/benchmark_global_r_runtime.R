#!/usr/bin/env Rscript
# Global PDP/ALE benchmark against R packages.
# Run from package root:
# Rscript simulation/benchmark_global_r_runtime.R --datadir simulation/data/global_r_runtime --outdir simulation/results/global_r_runtime

Sys.setenv(
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", "1"),
  OMP_THREAD_LIMIT = Sys.getenv("OMP_THREAD_LIMIT", "1"),
  OMP_PROC_BIND = Sys.getenv("OMP_PROC_BIND", "FALSE"),
  KMP_INIT_AT_FORK = Sys.getenv("KMP_INIT_AT_FORK", "FALSE"),
  KMP_AFFINITY = Sys.getenv("KMP_AFFINITY", "disabled"),
  OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", "1"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", "1"),
  VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS", "1"),
  NUMEXPR_NUM_THREADS = Sys.getenv("NUMEXPR_NUM_THREADS", "1"),
  DATATABLE_NUM_THREADS = Sys.getenv("DATATABLE_NUM_THREADS", "1"),
  RCPP_PARALLEL_NUM_THREADS = Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "1")
)

parse_flag = function(x) {
  s = tolower(trimws(as.character(x)))
  if (s %in% c("true", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "0", "no", "n", "off")) return(FALSE)
  FALSE
}

load_gadget_for_benchmark = function() {
  use_load_all = parse_flag(Sys.getenv("GADGET_BENCH_LOAD_ALL", "true"))
  if (use_load_all) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Install pkgload or set GADGET_BENCH_LOAD_ALL=false to use the installed gadget package")
    }
    pkgload::load_all(".", quiet = TRUE)
    message("Loaded gadget from local source with pkgload::load_all().")
  } else if (!requireNamespace("gadget", quietly = TRUE)) {
    stop("Install gadget or set GADGET_BENCH_LOAD_ALL=true")
  } else {
    library(gadget)
    message("Loaded installed gadget package.")
  }
}

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: gadget") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1L) == "Package: gadget") {
    setwd("..")
  } else {
    stop("Run from GADGET package root")
  }
}

load_gadget_for_benchmark()

library(data.table)
setDTthreads(1L)

args = commandArgs(trailingOnly = TRUE)
datadir = "simulation/data/global_r_runtime"
outdir = "simulation/results/global_r_runtime"
reps = 20L
predict_reps = 20L
N_vec = c(1000L, 2500L, 5000L, 10000L)
D_vec = c(5L, 10L, 20L)
fixed_N = 1000L
fixed_D = 10L
default_n_grid = 20L
default_n_intervals = 20L
n_grid_vec = c(10L, 20L, 50L)
n_intervals_vec = c(10L, 20L, 50L)
fail_fast = FALSE
model_types = c("rf", "toy")
sub_experiments = c("vs_N", "vs_D", "vs_res")

parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_chr_vec = function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])

i = 1L
while (i <= length(args)) {
  if (args[i] == "--datadir" && i < length(args)) {
    datadir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--outdir" && i < length(args)) {
    outdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--reps" && i < length(args)) {
    reps = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--predict-reps" && i < length(args)) {
    predict_reps = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--N-vec" && i < length(args)) {
    N_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--D" && i < length(args)) {
    fixed_D = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--D-vec" && i < length(args)) {
    D_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fixed-N" && i < length(args)) {
    fixed_N = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fixed-D" && i < length(args)) {
    fixed_D = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-grid" && i < length(args)) {
    default_n_grid = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-intervals" && i < length(args)) {
    default_n_intervals = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-grid-vec" && i < length(args)) {
    n_grid_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-intervals-vec" && i < length(args)) {
    n_intervals_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-int-vec" && i < length(args)) {
    n_intervals_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fail-fast" && i < length(args)) {
    fail_fast = parse_flag(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--models" && i < length(args)) {
    model_types = strsplit(args[i + 1L], ",", fixed = TRUE)[[1L]]
    model_types = trimws(model_types[nzchar(model_types)])
    i = i + 2L
  } else if (args[i] == "--sub-experiments" && i < length(args)) {
    sub_experiments = parse_chr_vec(args[i + 1L])
    sub_experiments = sub_experiments[nzchar(sub_experiments)]
    i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load_data = function(N, D) {
  path = file.path(datadir, sprintf("benchmark_N%d_D%d_seed21.csv", N, D))
  if (!file.exists(path)) return(NULL)
  as.data.frame(fread(path))
}

toy_pred_fun = function(model, newdata) {
  x1 = newdata[["x1"]]
  x2 = newdata[["x2"]]
  x3 = newdata[["x3"]]
  5 * x1 + 5 * x2 + ifelse(x3 > 0, 10 * x1, 0) - ifelse(x3 > 0, 10 * x2, 0)
}

rf_config = list(
  num_trees = 100L,
  min_node_size = 1L,
  replace = TRUE,
  sample_fraction = 1.0,
  splitrule = "variance",
  respect_unordered_factors = "ignore",
  num_threads = 1L,
  seed = 21L
)

fit_rf = function(dat) {
  p = ncol(dat) - 1L
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install ranger for the RF benchmark")
  }
  ranger::ranger(
    y ~ .,
    data = dat,
    num.trees = rf_config$num_trees,
    mtry = p,
    min.node.size = rf_config$min_node_size,
    replace = rf_config$replace,
    sample.fraction = rf_config$sample_fraction,
    splitrule = rf_config$splitrule,
    respect.unordered.factors = rf_config$respect_unordered_factors,
    num.threads = rf_config$num_threads,
    seed = rf_config$seed
  )
}

fit_mlr3_rf = function(dat) {
  p = ncol(dat) - 1L
  if (!requireNamespace("mlr3", quietly = TRUE) || !requireNamespace("mlr3learners", quietly = TRUE)) {
    stop("Install mlr3 and mlr3learners for the mlr3 RF sensitivity benchmark")
  }
  suppressPackageStartupMessages(loadNamespace("mlr3learners"))
  task = mlr3::as_task_regr(dat, target = "y", id = "global_r_mlr3_rf")
  learner = mlr3::lrn(
    "regr.ranger",
    num.trees = rf_config$num_trees,
    mtry = p,
    min.node.size = rf_config$min_node_size,
    replace = rf_config$replace,
    sample.fraction = rf_config$sample_fraction,
    splitrule = rf_config$splitrule,
    respect.unordered.factors = rf_config$respect_unordered_factors,
    num.threads = rf_config$num_threads,
    seed = rf_config$seed
  )
  learner$train(task)
  learner
}

rf_pred_fun = function(model, newdata) {
  as.numeric(predict(model, data = newdata)$predictions)
}

predict_for_model = function(model_type) {
  if (identical(model_type, "rf")) {
    rf_pred_fun
  } else if (identical(model_type, "mlr3_rf")) {
    NULL
  } else {
    toy_pred_fun
  }
}

fit_model = function(dat, model_type) {
  if (identical(model_type, "rf")) {
    fit_rf(dat)
  } else if (identical(model_type, "mlr3_rf")) {
    fit_mlr3_rf(dat)
  } else {
    "toy"
  }
}

grid_for_feature = function(x, n) {
  sort(unique(as.numeric(stats::quantile(x, probs = seq(0, 1, length.out = n), type = 7, na.rm = TRUE))))
}

make_iml_predictor = function(model, X, y, pred_fun) {
  iml::Predictor$new(
    model = model,
    data = X,
    y = y,
    predict.function = function(object, newdata) pred_fun(object, newdata)
  )
}

make_dalex_explainer = function(model, X, y, pred_fun) {
  DALEX::explain(
    model = model,
    data = X,
    y = y,
    predict_function = function(object, newdata) pred_fun(object, newdata),
    label = "model",
    verbose = FALSE,
    precalculate = FALSE
  )
}

run_gadget_pdp = function(dat, model, pred_fun, engine, n_grid) {
  gadget:::calculate_pd(
    model = model,
    data = dat,
    target_feature_name = "y",
    feature_set = NULL,
    predict_fun = pred_fun,
    n_grid = n_grid,
    pd_engine = engine
  )
}

run_gadget_ale = function(dat, model, pred_fun, engine, n_intervals) {
  features = setdiff(colnames(dat), "y")
  if (identical(engine, "cpp")) {
    gadget:::calculate_ale_fast(
      model = model,
      data = dat,
      feature_set = features,
      target_feature_name = "y",
      n_intervals = n_intervals,
      predict_fun = pred_fun
    )
  } else {
    gadget:::calculate_ale(
      model = model,
      data = dat,
      feature_set = features,
      target_feature_name = "y",
      n_intervals = n_intervals,
      predict_fun = pred_fun
    )
  }
}

run_pdp = function(dat, model, pred_fun, n_grid) {
  X = dat[, setdiff(colnames(dat), "y"), drop = FALSE]
  features = colnames(X)
  for (feat in features) {
    grid = data.frame(value = grid_for_feature(X[[feat]], n_grid))
    names(grid) = feat
    pdp::partial(
      object = model,
      pred.var = feat,
      pred.grid = grid,
      pred.fun = function(object, newdata) pred_fun(object, newdata),
      train = X,
      ice = TRUE,
      plot = FALSE,
      progress = FALSE,
      parallel = FALSE
    )
  }
}

run_iml_pdp = function(dat, model, pred_fun, n_grid) {
  X = dat[, setdiff(colnames(dat), "y"), drop = FALSE]
  predictor = make_iml_predictor(model, X, dat$y, pred_fun)
  for (feat in colnames(X)) {
    iml::FeatureEffect$new(
      predictor = predictor,
      feature = feat,
      method = "ice",
      grid.size = n_grid
    )
  }
}

run_iml_ale = function(dat, model, pred_fun, n_intervals) {
  X = dat[, setdiff(colnames(dat), "y"), drop = FALSE]
  predictor = make_iml_predictor(model, X, dat$y, pred_fun)
  for (feat in colnames(X)) {
    iml::FeatureEffect$new(
      predictor = predictor,
      feature = feat,
      method = "ale",
      grid.size = n_intervals
    )
  }
}

run_ingredients_pdp = function(dat, model, pred_fun, n_grid) {
  X = dat[, setdiff(colnames(dat), "y"), drop = FALSE]
  explainer = make_dalex_explainer(model, X, dat$y, pred_fun)
  ingredients::partial_dependence(
    explainer,
    variables = colnames(X),
    N = nrow(X),
    grid_points = n_grid,
    variable_splits = lapply(X, grid_for_feature, n = n_grid),
    variable_type = "numerical"
  )
}

run_ingredients_ale = function(dat, model, pred_fun, n_intervals) {
  X = dat[, setdiff(colnames(dat), "y"), drop = FALSE]
  explainer = make_dalex_explainer(model, X, dat$y, pred_fun)
  ingredients::accumulated_dependence(
    explainer,
    variables = colnames(X),
    N = nrow(X),
    variable_splits = lapply(X, grid_for_feature, n = n_intervals),
    grid_points = n_intervals,
    variable_type = "numerical"
  )
}

method_specs = list(
  list(package = "gadget", impl = "r", method = "global_pdp", requires = character(),
    model_types = c("rf", "toy", "mlr3_rf"),
    runner = function(dat, model, pred_fun, cell) run_gadget_pdp(dat, model, pred_fun, "r", cell$n_grid)),
  list(package = "gadget", impl = "r", method = "global_ale", requires = character(),
    model_types = c("rf", "toy", "mlr3_rf"),
    runner = function(dat, model, pred_fun, cell) run_gadget_ale(dat, model, pred_fun, "r", cell$n_intervals)),
  list(package = "pdp", impl = "default", method = "global_pdp", requires = "pdp",
    model_types = c("rf", "toy"),
    runner = function(dat, model, pred_fun, cell) run_pdp(dat, model, pred_fun, cell$n_grid)),
  list(package = "iml", impl = "ice", method = "global_pdp", requires = "iml",
    model_types = c("rf", "toy"),
    runner = function(dat, model, pred_fun, cell) run_iml_pdp(dat, model, pred_fun, cell$n_grid)),
  list(package = "iml", impl = "ale", method = "global_ale", requires = "iml",
    model_types = c("rf", "toy"),
    runner = function(dat, model, pred_fun, cell) run_iml_ale(dat, model, pred_fun, cell$n_intervals)),
  list(package = "DALEX/ingredients", impl = "partial_dependence", method = "global_pdp",
    model_types = c("rf", "toy"),
    requires = c("DALEX", "ingredients"),
    runner = function(dat, model, pred_fun, cell) run_ingredients_pdp(dat, model, pred_fun, cell$n_grid)),
  list(package = "DALEX/ingredients", impl = "accumulated_dependence", method = "global_ale",
    model_types = c("rf", "toy"),
    requires = c("DALEX", "ingredients"),
    runner = function(dat, model, pred_fun, cell) run_ingredients_ale(dat, model, pred_fun, cell$n_intervals))
)

make_cells = function(spec) {
  is_pdp = grepl("pdp", spec$method)
  res_vec = if (is_pdp) n_grid_vec else n_intervals_vec
  make_cell = function(sub_experiment, N, D, res) {
    data.frame(
      sub_experiment = sub_experiment,
      N = as.integer(N),
      D = as.integer(D),
      n_grid = if (is_pdp) as.integer(res) else NA_integer_,
      n_intervals = if (!is_pdp) as.integer(res) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  cells = list()
  if ("vs_N" %in% sub_experiments) {
    cells[[length(cells) + 1L]] = do.call(rbind, lapply(N_vec, function(N) make_cell("vs_N", N, fixed_D,
      if (is_pdp) default_n_grid else default_n_intervals)))
  }
  if ("vs_D" %in% sub_experiments) {
    cells[[length(cells) + 1L]] = do.call(rbind, lapply(D_vec, function(D) make_cell("vs_D", fixed_N, D,
      if (is_pdp) default_n_grid else default_n_intervals)))
  }
  if ("vs_res" %in% sub_experiments) {
    cells[[length(cells) + 1L]] = do.call(rbind, lapply(res_vec, function(res) make_cell("vs_res", fixed_N, fixed_D, res)))
  }
  do.call(rbind, cells)
}

record_row = function(
  rows, spec, model_type, cell, time_sec, repetition, status = "ok", error_message = NA_character_
) {
  rows[[length(rows) + 1L]] = data.frame(
    module = "global_r",
    package = spec$package,
    impl = spec$impl,
    method = spec$method,
    model_type = model_type,
    sub_experiment = cell$sub_experiment,
    N = cell$N,
    D = cell$D,
    n_grid = cell$n_grid,
    n_intervals = cell$n_intervals,
    repetition = repetition,
    time_sec = time_sec,
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  rows
}

time_spec = function(spec, dat, model, pred_fun, cell) {
  missing = spec$requires[!vapply(spec$requires, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop(sprintf("missing package(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  gc(FALSE)
  system.time(spec$runner(dat, model, pred_fun, cell))[["elapsed"]]
}

run_model = function(model_type) {
  rows = list()
  model_cache = list()
  pred_fun = predict_for_model(model_type)

  get_model = function(dat, N, D) {
    key = sprintf("%s_N%d_D%d", model_type, N, D)
    if (!key %in% names(model_cache)) {
      model_cache[[key]] <<- fit_model(dat, model_type)
    }
    model_cache[[key]]
  }

  for (spec in method_specs) {
    if (!(model_type %in% spec$model_types)) next
    message(sprintf("=== global R %s: %s / %s / %s ===", model_type, spec$package, spec$impl, spec$method))
    cells = make_cells(spec)
    for (i in seq_len(nrow(cells))) {
      cell = cells[i, ]
      dat = load_data(cell$N, cell$D)
      if (is.null(dat)) next
      model = get_model(dat, cell$N, cell$D)
      res_label = if (grepl("pdp", spec$method)) cell$n_grid else cell$n_intervals
      message(sprintf("  %s: N=%d D=%d res=%d", cell$sub_experiment, cell$N, cell$D, res_label))
      tryCatch(
        time_spec(spec, dat, model, pred_fun, cell),
        error = function(e) {
          if (isTRUE(fail_fast)) stop(e)
          message("    warmup skipped/failed: ", conditionMessage(e))
          invisible(NA_real_)
        }
      )
      for (r in seq_len(reps)) {
        out = tryCatch(
          list(ok = TRUE, time = time_spec(spec, dat, model, pred_fun, cell), error = NA_character_),
          error = function(e) {
            if (isTRUE(fail_fast)) stop(e)
            list(ok = FALSE, time = NA_real_, error = conditionMessage(e))
          }
        )
        rows = record_row(
          rows = rows,
          spec = spec,
          model_type = model_type,
          cell = cell,
          time_sec = out$time,
          repetition = r,
          status = if (out$ok) "ok" else "error",
          error_message = out$error
        )
      }
    }
  }
  rows
}

run_predict_baseline = function() {
  rows = list()
  baseline_cells = unique(rbind(
    data.frame(N = N_vec, D = fixed_D),
    data.frame(N = fixed_N, D = D_vec)
  ))
  for (i in seq_len(nrow(baseline_cells))) {
    cell = baseline_cells[i, ]
    dat = load_data(cell$N, cell$D)
    if (is.null(dat)) next
    model = fit_rf(dat)
    X = dat[, setdiff(colnames(dat), "y"), drop = FALSE]
    times = replicate(predict_reps, system.time(rf_pred_fun(model, X))[["elapsed"]])
    rows[[length(rows) + 1L]] = data.frame(
      package = "r",
      model_type = "rf",
      N = cell$N,
      D = cell$D,
      predict_time_mean = mean(times),
      predict_time_sd = stats::sd(times),
      n_rep = predict_reps
    )
  }
  if (length(rows)) {
    out = do.call(rbind, rows)
    write.csv(out, file.path(outdir, "global_r_predict_baseline.csv"), row.names = FALSE)
    message("Written: global_r_predict_baseline.csv")
  }
}

if ("rf" %in% model_types) {
  message("=== RF prediction baseline ===")
  run_predict_baseline()
  rf_rows = run_model("rf")
  write.csv(do.call(rbind, rf_rows), file.path(outdir, "global_r_runtime_rf.csv"), row.names = FALSE)
  message("Written: global_r_runtime_rf.csv")
}

if ("toy" %in% model_types) {
  toy_rows = run_model("toy")
  write.csv(do.call(rbind, toy_rows), file.path(outdir, "global_r_runtime_toy.csv"), row.names = FALSE)
  message("Written: global_r_runtime_toy.csv")
}

if ("mlr3_rf" %in% model_types) {
  mlr3_rows = run_model("mlr3_rf")
  write.csv(do.call(rbind, mlr3_rows), file.path(outdir, "global_r_runtime_mlr3_rf.csv"), row.names = FALSE)
  message("Written: global_r_runtime_mlr3_rf.csv")
}
