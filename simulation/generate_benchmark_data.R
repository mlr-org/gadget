#!/usr/bin/env Rscript
# Generate shared benchmark data.
# X ~ Uniform(-1,1), y from balanced DGP with x3 as moderator.
# Run: Rscript simulation/generate_benchmark_data.R [--outdir DIR] [--seed N] [--N-vec "..."] [--D-vec "..."]
# Defaults match the **small** preset; `run_benchmark.sh large` passes larger N, D vectors and
# `--outdir simulation/data/benchmark_large` so outputs do not overwrite a small run.

args <- commandArgs(trailingOnly = TRUE)
outdir <- "simulation/data/benchmark"
seed <- 21L
i <- 1
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) { outdir <- args[i + 1]; i <- i + 2 }
  else if (args[i] == "--seed" && i < length(args)) { seed <- as.integer(args[i + 1]); i <- i + 2 }
  else { i <- i + 1 }
}

generate_one <- function(N, D, seed, outdir) {
  set.seed(seed)
  X <- matrix(runif(N * D, -1, 1), nrow = N, ncol = D)
  colnames(X) <- paste0("x", seq_len(D))
  x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
  y_det <- 5 * x1 + 5 * x2 + ifelse(x3 > 0, 10 * x1, 0) - ifelse(x3 > 0, 10 * x2, 0)
  eps <- rnorm(N, 0, 0.1 * sd(y_det))
  dat <- data.frame(X, y = y_det + eps)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fname <- file.path(outdir, sprintf("benchmark_N%d_D%d_seed%d.csv", N, D, seed))
  write.csv(dat, fname, row.names = FALSE)
  message("Written: ", fname)
}

parse_int_vec <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
N_vec <- c(500L, 1000L, 5000L)
D_vec <- c(5L, 10L, 20L)
i <- 1
while (i <= length(args)) {
  if (args[i] == "--N-vec" && i < length(args)) { N_vec <- parse_int_vec(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--D-vec" && i < length(args)) { D_vec <- parse_int_vec(args[i + 1L]); i <- i + 2L }
  else { i <- i + 1L }
}
for (N in N_vec)
  for (D in D_vec)
    generate_one(N, D, seed, outdir)
