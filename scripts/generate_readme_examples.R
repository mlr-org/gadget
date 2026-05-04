# Generate README example outputs
# Run: Rscript scripts/generate_readme_examples.R

options(warn = -1)
dir.create("figures", showWarnings = FALSE)

library(gadget)
library(mlr3)
library(mlr3learners)
library(ISLR2)

# ---- ALE + Bike (README example) ----
cat("\n=== ALE + Bike ===\n")
set.seed(123)
bike = Bikeshare[sample(seq_len(nrow(Bikeshare)), 1000), ]
bike$workingday = as.factor(bike$workingday)
bike_data = bike[, c("hr", "temp", "workingday", "bikers")]
names(bike_data)[names(bike_data) == "bikers"] = "target"

task = TaskRegr$new(id = "bike", backend = bike_data, target = "target")
learner = lrn("regr.ranger")
learner$train(task)

tree_ale_bike = GadgetTree$new(
  strategy = AleStrategy$new(),
  n_split = 2,
  impr_par = 0.01,
  min_node_size = 50
)
tree_ale_bike$fit(
  data = bike_data,
  target_feature_name = "target",
  model = learner,
  n_intervals = 10
)

split_ale_bike = tree_ale_bike$extract_split_info()
cat("Split info (ALE Bike):\n")
print(split_ale_bike)
sink("figures/split_info_ale_bike.txt")
print(split_ale_bike)
sink()

# save tree structure plot
png("figures/ale_bike_tree_structure.png", width = 800, height = 500)
tree_ale_bike$plot_tree_structure()
dev.off()

# collect all regional ALE plots (depth x node)
pl_ale_bike_all = tree_ale_bike$plot(
  data = bike_data,
  target_feature_name = "target",
  features = c("hr", "temp"),
  mean_center = TRUE,
  show_plot = FALSE
)

# save per-node plots
if (length(pl_ale_bike_all) > 0) {
  for (d in seq_along(pl_ale_bike_all)) {
    nodes = pl_ale_bike_all[[d]]
    if (!is.list(nodes) || !length(nodes)) next
    for (k in seq_along(nodes)) {
      fname = sprintf("figures/ale_bike_depth%d_node%d.png", d, k)
      png(fname, width = 800, height = 500)
      print(nodes[[k]])
      dev.off()
    }
  }
}

# ---- PD + Bike (requires iml, README example) ----
cat("\n=== PD + Bike ===\n")
library(iml)
bike_x = bike_data[, c("hr", "temp", "workingday")]
bike_y = bike_data$target
predictor_bike = iml::Predictor$new(model = learner, data = bike_x, y = bike_y)
effect_bike = iml::FeatureEffects$new(predictor_bike, method = "ice", grid.size = 20)

tree_pd_bike = GadgetTree$new(strategy = PdStrategy$new(), n_split = 2, min_node_size = 50)
tree_pd_bike$fit(data = bike_data, target_feature_name = "target", effect = effect_bike)

split_pd_bike = tree_pd_bike$extract_split_info()
cat("Split info (PD Bike):\n")
print(split_pd_bike)
sink("figures/split_info_pd_bike.txt")
print(split_pd_bike)
sink()

# save tree structure plot
png("figures/pd_bike_tree_structure.png", width = 800, height = 500)
tree_pd_bike$plot_tree_structure()
dev.off()

# collect all regional PD/ICE plots (depth x node)
pl_pd_bike_all = tree_pd_bike$plot(
  data = bike_data,
  target_feature_name = "target",
  effect = effect_bike,
  features = c("hr", "temp"),
  show_plot = FALSE
)

# save per-node plots
if (length(pl_pd_bike_all) > 0) {
  for (d in seq_along(pl_pd_bike_all)) {
    nodes = pl_pd_bike_all[[d]]
    if (!is.list(nodes) || !length(nodes)) next
    for (k in seq_along(nodes)) {
      fname = sprintf("figures/pd_bike_depth%d_node%d.png", d, k)
      png(fname, width = 800, height = 500)
      print(nodes[[k]])
      dev.off()
    }
  }
}

cat("\nDone. Outputs in figures/\n")
