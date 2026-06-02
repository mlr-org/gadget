#!/usr/bin/env Rscript
# Summarize the archived large split-search benchmark and write the paper-facing figure.

Sys.setenv(
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", "1"),
  OMP_THREAD_LIMIT = Sys.getenv("OMP_THREAD_LIMIT", "1"),
  OMP_PROC_BIND = Sys.getenv("OMP_PROC_BIND", "FALSE"),
  KMP_INIT_AT_FORK = Sys.getenv("KMP_INIT_AT_FORK", "FALSE"),
  KMP_AFFINITY = Sys.getenv("KMP_AFFINITY", "disabled"),
  OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", "1"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", "1"),
  VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS", "1"),
  DATATABLE_NUM_THREADS = Sys.getenv("DATATABLE_NUM_THREADS", "1"),
  RCPP_PARALLEL_NUM_THREADS = Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "1")
)

args = commandArgs(trailingOnly = TRUE)
indir = "simulation/results/split_search_runtime_large"
figdir = "simulation/results/split_search_figures"

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]
    i = i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir = args[i + 1L]
    i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)
library(scales)

setDTthreads(1L)

summary_path = file.path(indir, "summary.csv")
if (!file.exists(summary_path)) {
  stop("No split-search summary found at ", summary_path)
}

dt = fread(summary_path)
dt = dt[method %in% c("regional_pdp_split", "regional_ale_split")]
if (nrow(dt) == 0L) {
  stop("No split-search rows found in ", summary_path)
}

dt[, effect := fifelse(method == "regional_pdp_split", "PDP", "ALE")]
dt[, model_label := fifelse(model_type == "rf", "RF model", "Toy model")]
dt[, method_label := fcase(
  package == "gadget", "gadget-split",
  package == "effector", "effector-split",
  default = package
)]
dt[, resolution := fifelse(effect == "PDP", as.integer(n_grid), as.integer(n_intervals))]
dt[, time := fifelse(is.finite(method_time_mean), method_time_mean, time_mean)]
dt[, time_sd_plot := fifelse(is.finite(method_time_sd), method_time_sd, time_sd)]

make_sweep = function(data, sub_experiment, x_var, sweep_short, fixed_desc) {
  out = copy(data)
  out[, sub_experiment := sub_experiment]
  out[, x_value := as.numeric(get(x_var))]
  out[, panel_title := sprintf("%s / %s - %s\n%s", effect, model_label, sweep_short, fixed_desc)]
  out
}

plot_dt = rbindlist(list(
  make_sweep(
    dt[D == 20L & resolution == 20L],
    "vs_N",
    "N",
    "vs N",
    "D = 20, resolution = 20"
  ),
  make_sweep(
    dt[N == 10000L & resolution == 20L],
    "vs_D",
    "D",
    "vs D",
    "N = 10,000, resolution = 20"
  ),
  make_sweep(
    dt[N == 10000L & D == 20L],
    "vs_res",
    "resolution",
    "vs resolution",
    "N = 10,000, D = 20"
  )
), use.names = TRUE, fill = TRUE)

panel_rows = c("PDP / RF model", "PDP / Toy model", "ALE / RF model", "ALE / Toy model")
panel_cols = c(
  "vs N\nD = 20, resolution = 20",
  "vs D\nN = 10,000, resolution = 20",
  "vs resolution\nN = 10,000, D = 20"
)
panel_levels = unlist(lapply(panel_rows, function(row) paste(row, panel_cols, sep = " - ")))
plot_dt[, panel_title := factor(panel_title, levels = panel_levels)]
plot_dt[, time_lower := pmax(time - time_sd_plot, 0)]
plot_dt[, time_upper := time + time_sd_plot]

palette_values = c(
  "gadget-split" = "#1f77b4",
  "effector-split" = "#2ca02c"
)

format_axis_number = function(x) {
  format(x, trim = TRUE, digits = 3L, scientific = FALSE, big.mark = ",")
}

p = ggplot(
  plot_dt,
  aes(x = x_value, y = time, color = method_label, fill = method_label, group = method_label)
) +
  geom_ribbon(aes(ymin = time_lower, ymax = time_upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ panel_title, ncol = 3L, scales = "free") +
  scale_x_continuous(labels = comma, breaks = sort(unique(plot_dt$x_value))) +
  scale_y_continuous(labels = format_axis_number) +
  scale_color_manual(values = palette_values, breaks = names(palette_values), drop = FALSE) +
  scale_fill_manual(values = alpha(palette_values, 0.25), breaks = names(palette_values), drop = FALSE, guide = "none") +
  labs(
    title = "Regional split-search runtime",
    x = "Value of varied parameter",
    y = "Runtime (s)",
    color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", lineheight = 0.95)
  )

out_path = file.path(figdir, "regional_methods.png")
ggsave(out_path, p, width = 13, height = 10, dpi = 220)
message("Written: ", out_path)
