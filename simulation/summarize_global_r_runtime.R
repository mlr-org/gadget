#!/usr/bin/env Rscript
# Summarize module-level efficiency benchmark results.

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

args = commandArgs(trailingOnly = TRUE)
indir = "simulation/results/global_r_runtime"
figdir = "simulation/results/paper_figures"
fixed_D = 10L
model_types = c("rf", "toy")
include_mlr3 = FALSE

parse_model_vec = function(x) {
  x = strsplit(x, ",", fixed = TRUE)[[1L]]
  trimws(x[nzchar(x)])
}

parse_flag = function(x) {
  s = tolower(trimws(as.character(x)))
  if (s %in% c("true", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "0", "no", "n", "off")) return(FALSE)
  FALSE
}

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--fixed-D" && i < length(args)) {
    fixed_D = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--models" && i < length(args)) {
    model_types = parse_model_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--include-mlr3" && i < length(args)) {
    include_mlr3 = parse_flag(args[i + 1L]); i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)
library(scales)

setDTthreads(1L)

load_csv = function(filename) {
  path = file.path(indir, filename)
  if (!file.exists(path)) return(data.table())
  fread(path)
}

csv_files = list()
if ("rf" %in% model_types) {
  csv_files = c(csv_files, list(load_csv("global_r_runtime_rf.csv")))
}
if ("toy" %in% model_types) {
  csv_files = c(csv_files, list(load_csv("global_r_runtime_toy.csv")))
}
if ("mlr3_rf" %in% model_types || isTRUE(include_mlr3)) {
  csv_files = c(csv_files, list(load_csv("global_r_runtime_mlr3_rf.csv")))
}
dt = rbindlist(csv_files, use.names = TRUE, fill = TRUE)

if (nrow(dt) == 0L) {
  stop("No benchmark data found. Run simulation/run_global_r_runtime.sh first.")
}

if (!("status" %in% names(dt))) dt[, status := "ok"]
dt[is.na(status) | status == "", status := "ok"]
if (!("sub_experiment" %in% names(dt))) dt[, sub_experiment := "vs_N"]
dt[is.na(sub_experiment) | sub_experiment == "", sub_experiment := "vs_N"]
n_err = dt[status == "error", .N]
if (n_err > 0L) {
  message("Warning: ", n_err, " benchmark row(s) have status=error; see error_message in raw CSV.")
}

dt[, n_grid := as.integer(fifelse(is.na(n_grid) | n_grid == "", NA_character_, as.character(n_grid)))]
dt[, n_intervals := as.integer(fifelse(
  is.na(n_intervals) | n_intervals == "",
  NA_character_,
  as.character(n_intervals)
))]
dt[, time_sec := as.numeric(time_sec)]
dt[package == "ingredients", package := "DALEX/ingredients"]
dt = dt[module == "global_r" & !(package == "gadget" & impl == "cpp")]
dt = dt[!(sub_experiment == "vs_N" & N == 500L)]

dt_ok = dt[status == "ok" & is.finite(time_sec)]
if (nrow(dt_ok) == 0L) {
  stop("No successful benchmark timings found.")
}

summary_dt = dt_ok[, .(
  time_median = stats::median(time_sec),
  time_q25 = stats::quantile(time_sec, probs = 0.25, names = FALSE, type = 7),
  time_q75 = stats::quantile(time_sec, probs = 0.75, names = FALSE, type = 7),
  time_min = min(time_sec),
  time_max = max(time_sec),
  time_mean = mean(time_sec),
  time_sd = stats::sd(time_sec),
  n_rep = .N
), by = .(module, package, impl, method, model_type, sub_experiment, N, D, n_grid, n_intervals)]

summary_dt[, label := fifelse(
  package == "gadget" & impl == "r",
  "gadget-r",
  package
)]

write.csv(summary_dt, file.path(indir, "summary.csv"), row.names = FALSE)
message("Written: ", file.path(indir, "summary.csv"))

palette_values = c(
  "gadget-r" = "#1f77b4",
  "pdp" = "#ff7f0e",
  "iml" = "#2ca02c",
  "DALEX/ingredients" = "#9467bd"
)

plot_runtime = function(data, title, filename) {
  if (nrow(data) == 0L) return(invisible(NULL))
  data = copy(data)
  data[, effect := fifelse(grepl("pdp", method), "PDP", "ALE")]
  data[, model_label := fifelse(model_type == "rf", "RF model", "Toy model")]
  data[, panel := ifelse(module == "global_r", paste(effect, model_label, sep = " / "), effect)]
  data[, panel := factor(panel, levels = c(
    "ALE / RF model",
    "ALE / Toy model",
    "PDP / RF model",
    "PDP / Toy model"
  ))]
  data[, sweep_label := fcase(
    sub_experiment == "vs_N", "Sample size N\nD = 10, resolution = 20",
    sub_experiment == "vs_D", "Feature dimension D\nN = 1,000, resolution = 20",
    sub_experiment == "vs_res", "PDP grid / ALE intervals\nN = 1,000, D = 10",
    default = sub_experiment
  )]
  data[, sweep_label := factor(sweep_label, levels = c(
    "Sample size N\nD = 10, resolution = 20",
    "Feature dimension D\nN = 1,000, resolution = 20",
    "PDP grid / ALE intervals\nN = 1,000, D = 10"
  ))]
  data[, x_value := fcase(
    sub_experiment == "vs_N", as.numeric(N),
    sub_experiment == "vs_D", as.numeric(D),
    grepl("pdp", method), as.numeric(n_grid),
    default = as.numeric(n_intervals)
  )]
  data[, time_median_plot := pmax(time_median, .Machine$double.eps)]
  data[, time_q25_plot := pmax(time_q25, .Machine$double.eps)]
  data[, time_q75_plot := pmax(time_q75, .Machine$double.eps)]

  labels = unique(data$label)
  colors = palette_values[labels]
  missing_colors = is.na(colors)
  if (any(missing_colors)) {
    colors[missing_colors] = hue_pal()(sum(missing_colors))
  }

  p = ggplot(
    data,
    aes(x = x_value, y = time_median_plot, color = label, fill = label, group = label)
  ) +
    geom_ribbon(aes(ymin = time_q25_plot, ymax = time_q75_plot), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_grid(panel ~ sweep_label, scales = "free") +
    scale_x_log10(labels = comma, breaks = sort(unique(data$x_value))) +
    scale_y_log10(labels = label_number(accuracy = 0.001, trim = TRUE)) +
    scale_color_manual(values = colors, breaks = labels, drop = FALSE) +
    scale_fill_manual(values = alpha(colors, 0.25), breaks = labels, drop = FALSE, guide = "none") +
    labs(title = title, x = "Value of varied parameter", y = "Median runtime (s, log scale)", color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold")
    )

  ggsave(file.path(figdir, filename), p, width = 11, height = 9, dpi = 220)
  message("Written: ", file.path(figdir, filename))
}

plot_runtime(
  summary_dt[module == "global_r"],
  "Global feature-effect computation in R",
  "global_r_methods.png"
)

message("Done.")
