#!/usr/bin/env Rscript
# Aggregate the structural recovery benchmark and write summary figures.

args = commandArgs(trailingOnly = TRUE)
indir = "simulation/results/structural_recovery"
figdir = "simulation/results/paper_figures"
paper_figdir = "paper/figures"

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]
    i = i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir = args[i + 1L]
    i = i + 2L
  } else if (args[i] == "--paper-figdir" && i < length(args)) {
    paper_figdir = args[i + 1L]
    i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(paper_figdir) && nzchar(paper_figdir)) {
  dir.create(paper_figdir, showWarnings = FALSE, recursive = TRUE)
}

library(data.table)
library(ggplot2)

variant_labels = c(
  num_0 = "Numeric, balanced split",
  num_04 = "Numeric, off-center split",
  cat = "Categorical split"
)
method_labels = c(
  gadget_pdp = "gadget PDP",
  gadget_ale = "gadget ALE",
  effector_rpdp = "effector RegionalPDP",
  effector_rale = "effector RegionalALE"
)
method_levels = c("gadget_pdp", "gadget_ale", "effector_rpdp", "effector_rale")

save_plot = function(filename, plot_obj, width, height) {
  ggsave(file.path(figdir, filename), plot_obj, width = width, height = height, dpi = 150)
  if (!is.null(paper_figdir) && nzchar(paper_figdir)) {
    ggsave(file.path(paper_figdir, filename), plot_obj, width = width, height = height, dpi = 150)
  }
}

gadget_res = fread(file.path(indir, "structural_recovery_gadget.csv"))
effector_res = fread(file.path(indir, "structural_recovery_effector.csv"))
dt = rbind(gadget_res, effector_res, fill = TRUE)

dt[, split_feat_correct := as.logical(split_feat_correct)]
dt[, split_pt_error := as.numeric(split_pt_error)]
dt[, node_acc := as.numeric(node_acc)]
dt[, method := factor(method, levels = method_levels)]
dt[, method_label := factor(method_labels[as.character(method)], levels = unname(method_labels))]
dt[, variant_label := factor(variant_labels[variant], levels = unname(variant_labels))]
dt[, d_label := factor(sprintf("D = %s", D), levels = sprintf("D = %s", sort(unique(D))))]

agg = dt[, .(
  split_feat_hit_mean = mean(split_feat_correct, na.rm = TRUE),
  split_feat_hit_sd = stats::sd(split_feat_correct, na.rm = TRUE),
  split_pt_mae_mean = mean(split_pt_error, na.rm = TRUE),
  split_pt_mae_sd = stats::sd(split_pt_error, na.rm = TRUE),
  node_acc_mean = mean(node_acc, na.rm = TRUE),
  node_acc_sd = stats::sd(node_acc, na.rm = TRUE),
  n_rep = .N
), by = .(method, variant, N, D)]

summary_path = file.path(indir, "structural_recovery_summary.csv")
fwrite(agg, summary_path)
message("Written: ", summary_path)

if (nrow(dt) == 0L) {
  quit(save = "no")
}

p_hit = ggplot(dt, aes(x = N, y = as.integer(split_feat_correct), color = method_label)) +
  geom_jitter(height = 0.04, width = 0, alpha = 0.15, size = 0.6) +
  stat_summary(fun = mean, geom = "line", aes(group = interaction(method, variant)), linewidth = 0.8) +
  facet_grid(variant_label ~ d_label) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = expression(paste("Oracle moderator recovery for the regional ", x[2], " effect")),
    subtitle = expression(paste("A hit means that the root split feature is ", x[3])),
    x = "Sample size N",
    y = "Hit rate",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot("accuracy_hit_rate.png", p_hit, width = 9, height = 7)

dt_num = dt[variant %in% c("num_0", "num_04") & is.finite(split_pt_error)]
if (nrow(dt_num) > 0L) {
  p_mae = ggplot(dt_num, aes(x = N, y = split_pt_error, color = method_label)) +
    geom_point(alpha = 0.12, size = 0.5) +
    stat_summary(fun = median, geom = "line", aes(group = interaction(method, variant)), linewidth = 0.8) +
    facet_grid(variant_label ~ d_label) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = expression(paste("Root split-point error on the numeric ", x[3], " variants")),
      subtitle = expression(paste("Absolute error relative to the oracle ", x[3], " threshold")),
      x = "Sample size N",
      y = "Absolute split-point error",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  save_plot("accuracy_split_point_mae.png", p_mae, width = 9, height = 6)
}

p_acc = ggplot(dt[is.finite(node_acc)], aes(x = N, y = node_acc, color = method_label)) +
  geom_point(alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "line", aes(group = interaction(method, variant)), linewidth = 0.8) +
  facet_grid(variant_label ~ d_label) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = expression(paste("Agreement with the oracle ", x[3], " partition")),
    subtitle = "Best label alignment of the selected root split",
    x = "Sample size N",
    y = "Node assignment accuracy",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot("accuracy_node_assignment.png", p_acc, width = 9, height = 7)

message("Written figures to ", figdir)
if (!is.null(paper_figdir) && nzchar(paper_figdir)) {
  message("Synced figures to ", paper_figdir)
}
