#!/usr/bin/env Rscript
# Generate shared CSV datasets for simulation/efficiency_benchmark_plan.md §9 (accuracy).
# Run from package root. Writes: simulation/data/accuracy/acc_N{N}_D{D}_{variant}_seed{seed}.csv

args <- commandArgs(trailingOnly = TRUE)
outdir <- "simulation/data/accuracy"
n_seeds <- 30L
N_vec <- c(200L, 500L, 1000L, 5000L)
D_vec <- c(5L, 10L, 20L)
variants <- c("num_0", "num_05", "cat")

i <- 1L
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1L]
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

make_data <- function(N, D, seed, variant) {
  set.seed(seed)
  X <- matrix(runif(N * D, -1, 1), nrow = N, ncol = D)
  colnames(X) <- paste0("x", seq_len(D))
  dat <- as.data.frame(X)

  if (variant == "cat") {
    dat$x3 <- factor(ifelse(rbinom(N, 1L, 0.5) == 1L, "B", "A"), levels = c("A", "B"))
  }

  x1 <- dat$x1
  x2 <- dat$x2
  y_d <- if (variant == "cat") {
    5 * x1 + 5 * x2 + ifelse(dat$x3 == "B", 10 * x1 - 10 * x2, 0)
  } else {
    th <- if (variant == "num_05") 0.5 else 0
    m <- dat$x3 > th
    5 * x1 + 5 * x2 + ifelse(m, 10 * x1 - 10 * x2, 0)
  }

  sd_eps <- 0.1 * stats::sd(y_d)
  if (!is.finite(sd_eps) || sd_eps == 0) sd_eps <- 0.01
  dat$y <- y_d + stats::rnorm(N, 0, sd_eps)
  dat
}

n_written <- 0L
for (variant in variants) {
  for (N in N_vec) {
    for (D in D_vec) {
      for (s in seq_len(n_seeds)) {
        seed <- 1000L + s
        dat <- make_data(N, D, seed, variant)
        fn <- file.path(outdir, sprintf("acc_N%d_D%d_%s_seed%d.csv", N, D, variant, seed))
        write.csv(dat, fn, row.names = FALSE)
        n_written <- n_written + 1L
      }
    }
  }
}
message("Written ", n_written, " files to ", outdir)
