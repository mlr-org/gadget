#!/usr/bin/env Rscript
# Accuracy benchmark (GADGET): split-feature recovery, split-point error, node agreement.
# Plan: simulation/efficiency_benchmark_plan.md §9.
# Requires datasets from: Rscript simulation/generate_accuracy_data.R
# Run from package root: Rscript simulation/benchmark_accuracy_gadget.R

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1) != "Package: gadget") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1) == "Package: gadget") {
    setwd("..")
  } else {
    stop("Run from GADGET package root")
  }
}

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(gadget)
  }
  library(data.table)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
outdir <- "simulation/results/accuracy"
datadir <- "simulation/data/accuracy"
n_seeds <- 30L
N_vec <- c(200L, 500L, 1000L, 5000L)
D_vec <- c(5L, 10L, 20L)
variants <- c("num_0", "num_05", "cat")
n_grid <- 20L
n_split <- 1L
min_node_size <- 50L

i <- 1L
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1L]
    i <- i + 2L
  } else if (args[i] == "--datadir" && i < length(args)) {
    datadir <- args[i + 1L]
    i <- i + 2L
  } else if (args[i] == "--n-seeds" && i < length(args)) {
    n_seeds <- as.integer(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--N-vec" && i < length(args)) {
    N_vec <- as.integer(strsplit(args[i + 1L], ",", fixed = TRUE)[[1L]])
    i <- i + 2L
  } else if (args[i] == "--D-vec" && i < length(args)) {
    D_vec <- as.integer(strsplit(args[i + 1L], ",", fixed = TRUE)[[1L]])
    i <- i + 2L
  } else if (args[i] == "--variants" && i < length(args)) {
    variants <- strsplit(args[i + 1L], ",", fixed = TRUE)[[1L]]
    i <- i + 2L
  } else {
    i <- i + 1L
  }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load shared CSV (same rows as Python benchmark)
# ---------------------------------------------------------------------------

load_accuracy_csv <- function(datadir, N, D, variant, seed) {
  fn <- file.path(datadir, sprintf("acc_N%d_D%d_%s_seed%d.csv", N, D, variant, seed))
  if (!file.exists(fn)) {
    stop("Missing dataset ", fn, ". Run: Rscript simulation/generate_accuracy_data.R", call. = FALSE)
  }
  dat <- as.data.frame(data.table::fread(fn))
  if (variant == "cat" && is.character(dat$x3)) {
    dat$x3 <- factor(dat$x3, levels = c("A", "B"))
  }
  dat
}

toy_pred_fun_factory <- function(variant) {
  function(model, newdata) {
    x1 <- newdata[["x1"]]
    x2 <- newdata[["x2"]]
    if (variant == "cat") {
      m <- newdata[["x3"]] == "B"
    } else {
      th <- if (variant == "num_05") 0.5 else 0
      m <- newdata[["x3"]] > th
    }
    5 * x1 + 5 * x2 + ifelse(m, 10 * x1 - 10 * x2, 0)
  }
}

true_threshold <- function(variant) {
  if (variant == "num_0") {
    0
  } else if (variant == "num_05") {
    0.5
  } else {
    NA_real_
  }
}

true_left_mask <- function(dat, variant) {
  if (variant == "cat") {
    as.character(dat$x3) == "A"
  } else {
    dat$x3 <= true_threshold(variant)
  }
}

root_split_info <- function(si) {
  r <- si[!is.na(si$split_feature) & si$node_type == "root", ]
  if (nrow(r) == 0L) {
    return(list(feature = NA_character_, value = NA))
  }
  list(feature = as.character(r$split_feature[1L]), value = r$split_value[1L])
}

split_feat_ok <- function(sf, variant) {
  !is.na(sf) && sf == "x3"
}

split_point_error <- function(sf, sv, variant) {
  if (variant == "cat") {
    return(NA_real_)
  }
  if (!split_feat_ok(sf, variant) || is.na(sv)) {
    return(NA_real_)
  }
  abs(as.numeric(sv) - true_threshold(variant))
}

node_assignment_acc <- function(dat, variant, sf, sv) {
  if (!split_feat_ok(sf, variant) || is.na(sv)) {
    return(NA_real_)
  }
  tl <- true_left_mask(dat, variant)
  if (inherits(dat[[sf]], "factor")) {
    pl <- dat[[sf]] == sv
    p2 <- dat[[sf]] != sv
  } else {
    sv <- as.numeric(sv)
    pl <- dat[[sf]] <= sv
    p2 <- dat[[sf]] > sv
  }
  max(mean(pl == tl), mean(p2 == tl))
}

run_gadget_pdp <- function(dat, n_grid, pred) {
  eff <- gadget:::calculate_pd(
    model = NULL,
    data = dat,
    target_feature_name = "y",
    predict_fun = pred,
    n_grid = n_grid,
    pd_engine = "cpp"
  )
  strat <- PdStrategy$new()
  tree <- GadgetTree$new(strategy = strat, n_split = n_split, min_node_size = min_node_size)
  tree$fit(data = dat, target_feature_name = "y", effect = eff)
  tree
}

run_gadget_ale <- function(dat, n_grid, pred) {
  strat <- AleStrategy$new()
  tree <- GadgetTree$new(strategy = strat, n_split = n_split, min_node_size = min_node_size)
  tree$fit(
    data = dat,
    target_feature_name = "y",
    model = NULL,
    n_intervals = n_grid,
    predict_fun = pred,
    order_method = "raw",
    feature_set = NULL,
    ale_engine = "cpp"
  )
  tree
}

record_gadget <- function(method_label, variant, N, D, seed, dat, tree) {
  si <- tree$extract_split_info()
  rs <- root_split_info(si)
  sf <- rs$feature
  sv <- rs$value
  hit <- split_feat_ok(sf, variant)
  data.table(
    package = "gadget",
    method = method_label,
    variant = variant,
    N = N,
    D = D,
    seed = seed,
    split_feat_correct = hit,
    split_pt_error = split_point_error(sf, sv, variant),
    node_acc = node_assignment_acc(dat, variant, sf, sv),
    effect_mse_node1 = NA_real_,
    effect_mse_node2 = NA_real_
  )
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

rows <- list()
for (variant in variants) {
  pred <- toy_pred_fun_factory(variant)
  for (N in N_vec) {
    for (D in D_vec) {
      for (s in seq_len(n_seeds)) {
        seed <- 1000L + s
        dat <- load_accuracy_csv(datadir, N, D, variant, seed)
        tr_pdp <- run_gadget_pdp(dat, n_grid, pred)
        rows[[length(rows) + 1L]] <- record_gadget("gadget_pdp", variant, N, D, seed, dat, tr_pdp)
        tr_ale <- run_gadget_ale(dat, n_grid, pred)
        rows[[length(rows) + 1L]] <- record_gadget("gadget_ale", variant, N, D, seed, dat, tr_ale)
      }
    }
  }
}

out <- rbindlist(rows)
fout <- file.path(outdir, "accuracy_gadget.csv")
fwrite(out, fout)
message("Written: ", fout)
