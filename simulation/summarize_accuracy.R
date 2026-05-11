#!/usr/bin/env Rscript
# Aggregate simulation/results/accuracy/accuracy_{gadget,effector}.csv and write summary + figures.
# Run from package root: Rscript simulation/summarize_accuracy.R

args <- commandArgs(trailingOnly = TRUE)
indir <- "simulation/results/accuracy"
figdir <- "simulation/results/figures"
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir <- args[i + 1L]
    i <- i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir <- args[i + 1L]
    i <- i + 2L
  } else {
    i <- i + 1L
  }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)

g <- fread(file.path(indir, "accuracy_gadget.csv"))
e <- fread(file.path(indir, "accuracy_effector.csv"))
dt <- rbind(g, e, fill = TRUE)

dt[, split_feat_correct := as.logical(split_feat_correct)]
dt[, split_pt_error := as.numeric(split_pt_error)]
dt[, node_acc := as.numeric(node_acc)]

agg <- dt[, .(
  split_feat_hit_mean = mean(split_feat_correct, na.rm = TRUE),
  split_feat_hit_sd = stats::sd(split_feat_correct, na.rm = TRUE),
  split_pt_mae_mean = mean(split_pt_error, na.rm = TRUE),
  split_pt_mae_sd = stats::sd(split_pt_error, na.rm = TRUE),
  node_acc_mean = mean(node_acc, na.rm = TRUE),
  node_acc_sd = stats::sd(node_acc, na.rm = TRUE),
  n_rep = .N
), by = .(method, variant, N, D)]

fwrite(agg, file.path(indir, "accuracy_summary.csv"))
message("Written: ", file.path(indir, "accuracy_summary.csv"))

if (nrow(dt) == 0L) {
  quit(save = "no")
}

p_hit <- ggplot(dt, aes(x = N, y = as.integer(split_feat_correct), color = method)) +
  geom_jitter(height = 0.04, width = 0, alpha = 0.15, size = 0.6) +
  { if (length(unique(dt$N)) > 1L) {
    stat_summary(fun = mean, geom = "line", aes(group = interaction(method, variant)), linewidth = 0.8)
  } else {
    NULL
  } } +
  facet_grid(variant ~ D, labeller = label_both) +
  scale_x_log10() +
  labs(
    title = "Split-feature recovery (first split equals x3)",
    x = "N",
    y = "Correct (0/1)",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figdir, "accuracy_hit_rate.png"), p_hit, width = 9, height = 7, dpi = 150)

dt_num <- dt[variant %in% c("num_0", "num_05") & is.finite(split_pt_error)]
if (nrow(dt_num) > 0L) {
  p_mae <- ggplot(dt_num, aes(x = N, y = split_pt_error, color = method)) +
    geom_point(alpha = 0.12, size = 0.5) +
    { if (length(unique(dt_num$N)) > 1L) {
      stat_summary(fun = median, geom = "line", aes(group = interaction(method, variant)), linewidth = 0.8)
    } else {
      NULL
    } } +
    facet_grid(variant ~ D, labeller = label_both) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Split-point error (numeric moderators)",
      x = "N",
      y = "abs(threshold_hat - threshold*)",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(file.path(figdir, "accuracy_split_point_mae.png"), p_mae, width = 9, height = 6, dpi = 150)
}

p_acc <- ggplot(dt[is.finite(node_acc)], aes(x = N, y = node_acc, color = method)) +
  geom_point(alpha = 0.1, size = 0.5) +
  { if (length(unique(dt$N)) > 1L) {
    stat_summary(fun = mean, geom = "line", aes(group = interaction(method, variant)), linewidth = 0.8)
  } else {
    NULL
  } } +
  facet_grid(variant ~ D, labeller = label_both) +
  scale_x_log10() +
  ylim(0, 1) +
  labs(
    title = "Node assignment agreement (best of <= / > sides)",
    x = "N",
    y = "Accuracy",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(figdir, "accuracy_node_assignment.png"), p_acc, width = 9, height = 7, dpi = 150)

message("Written figures to ", figdir)
