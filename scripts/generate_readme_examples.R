# Generate README example outputs
# Run: Rscript scripts/generate_readme_examples.R

options(warn = -1)
dir.create("figures", showWarnings = FALSE)
unlink(Sys.glob("figures/ale_bike_depth*_node*.png"))
unlink(Sys.glob("figures/ale_bike_depth*_id*.png"))
unlink(Sys.glob("figures/pd_bike_depth*_node*.png"))
unlink(Sys.glob("figures/pd_bike_depth*_id*.png"))
unlink(Sys.glob("figures/split_info_*_bike.*"))

devtools::load_all(".", quiet = TRUE)
library(mlr3)
library(mlr3learners)
library(ISLR2)

format_readme_value = function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  if (is.numeric(x)) {
    if (isTRUE(abs(x) >= 1000)) {
      return(format(round(x), big.mark = "", scientific = FALSE, trim = TRUE))
    }
    return(format(round(x, 3L), nsmall = 0L, scientific = FALSE, trim = TRUE))
  }
  as.character(x)
}

write_markdown_table = function(x, path) {
  cols = names(x)
  rows = vapply(seq_len(nrow(x)), function(i) {
    values = vapply(cols, function(col) format_readme_value(x[[col]][i]), character(1L))
    paste0("| ", paste(values, collapse = " | "), " |")
  }, character(1L))
  lines = c(
    paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|"),
    rows
  )
  writeLines(lines, path)
}

plot_file_node_id = function(node_name) {
  as.integer(sub("^Node_", "", node_name))
}

save_readme_plots = function(plot_list, prefix) {
  if (!length(plot_list)) {
    return(invisible(NULL))
  }
  for (depth_name in names(plot_list)) {
    nodes = plot_list[[depth_name]]
    if (!is.list(nodes) || !length(nodes)) next
    depth_id = as.integer(sub("^Depth_", "", depth_name))
    for (node_name in names(nodes)) {
      node_id = plot_file_node_id(node_name)
      fname = sprintf("figures/%s_bike_depth%d_id%d.png", prefix, depth_id, node_id)
      png(fname, width = 800, height = 500)
      print(nodes[[node_name]])
      dev.off()
    }
  }
  invisible(NULL)
}

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
write_markdown_table(split_ale_bike, "figures/split_info_ale_bike.md")

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
save_readme_plots(pl_ale_bike_all, "ale")

# ---- PD + Bike (README example) ----
cat("\n=== PD + Bike ===\n")

tree_pd_bike = GadgetTree$new(strategy = PdStrategy$new(), n_split = 2, min_node_size = 50)
tree_pd_bike$fit(data = bike_data, target_feature_name = "target", model = learner, n_grid = 20L)

split_pd_bike = tree_pd_bike$extract_split_info()
cat("Split info (PD Bike):\n")
print(split_pd_bike)
sink("figures/split_info_pd_bike.txt")
print(split_pd_bike)
sink()
write_markdown_table(split_pd_bike, "figures/split_info_pd_bike.md")

# save tree structure plot
png("figures/pd_bike_tree_structure.png", width = 800, height = 500)
tree_pd_bike$plot_tree_structure()
dev.off()

# collect all regional PD/ICE plots (depth x node)
pl_pd_bike_all = tree_pd_bike$plot(
  data = bike_data,
  target_feature_name = "target",
  features = c("hr", "temp"),
  show_plot = FALSE
)

# save per-node plots
save_readme_plots(pl_pd_bike_all, "pd")

cat("\nDone. Outputs in figures/\n")
