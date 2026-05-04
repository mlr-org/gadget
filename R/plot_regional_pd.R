#' Plot regional PD/ICE for one node.
#'
#' @param prepared_data (`list()`) \cr
#'   Prepared effect matrices per feature.
#' @param origin_data (`data.frame()`) \cr
#'   Original data.
#' @param target_feature_name (`character(1)`) \cr
#'   Target column.
#' @param node_idx (`integer(1)`) \cr
#'   Node index.
#' @param color_ice,color_pd (`character(1)`) \cr
#'   Colors.
#' @param ymin,ymax (`numeric(1)`) \cr
#'   Y-axis limits.
#' @param split_condition (`character(1)` or `NULL`) \cr
#'   Split condition label.
#' @param show_point,mean_center (`logical(1)`) \cr
#'   Plot options.
#'
#' @return (`list()`) \cr
#'   List of ggplot objects per feature.
#' @keywords internal
plot_regional_pd = function(prepared_data, origin_data, target_feature_name, node_idx,
  color_ice, color_pd, ymin, ymax, split_condition = NULL,
  show_point, mean_center) {
  plot = mlr3misc::map(names(prepared_data), function(feat) {
    data = prepared_data[[feat]]
    subset_idx = which(data$node == node_idx)
    data_subset = data[subset_idx, ]
    origin_data_subset = origin_data[subset_idx, ]
    if (feat %in% colnames(origin_data_subset) && is.factor(origin_data_subset[[feat]])) {
      origin_data_subset[[feat]] = factor_to_numeric(origin_data_subset[[feat]])
    }

    # data transformation
    data_subset = data_subset[, -ncol(data_subset), drop = FALSE]
    n_rows = nrow(data_subset)
    n_cols = ncol(data_subset)

    if (n_rows == 0 || n_cols == 0) {
      # return empty plot if no data
      return(ggplot() + theme_bw() + labs(title = "No data"))
    }

    # wide to long conversion: numeric grids use numeric x; categorical (factor levels as colnames) keep labels
    cn = colnames(data_subset)
    grid_as_num = suppressWarnings(as.numeric(cn))
    grid_numeric = length(cn) > 0L && all(!is.na(grid_as_num))
    grid_values = if (grid_numeric) grid_as_num else cn
    value_matrix = as.matrix(data_subset)

    ice_long = data.frame(
      grid = rep(grid_values, each = n_rows),
      value = as.vector(value_matrix),
      id = rep(seq_len(n_rows), times = n_cols),
      type = "ICE",
      stringsAsFactors = FALSE
    )

    valid_values = !is.na(ice_long$value)
    gv = ice_long$grid[valid_values]
    ok_g = !is.na(gv)
    if (sum(valid_values) > 0 && sum(ok_g) > 0) {
      mean_ice = tapply(ice_long$value[valid_values][ok_g], gv[ok_g], mean, na.rm = TRUE)
      gnam = names(mean_ice)
      if (is.null(gnam)) {
        gnam = rep(NA_character_, length(mean_ice))
      }
      pdp_grid = if (grid_numeric) as.numeric(gnam) else gnam
      pdp_centered = data.frame(
        grid = pdp_grid,
        value = as.numeric(mean_ice),
        id = NA,
        type = "PDP",
        stringsAsFactors = FALSE
      )
      plot_data = rbind(
        ice_long[, c("grid", "value", "id", "type")],
        pdp_centered[, c("grid", "value", "id", "type")]
      )
    } else {
      plot_data = ice_long[, c("grid", "value", "id", "type")]
    }

    if (!grid_numeric) {
      plot_data$grid = factor(plot_data$grid, levels = unique(cn))
    }

    # check if we can draw lines
    noline = length(unique(plot_data$grid[!is.na(plot_data$value)])) < 2

    # nolint start: object_usage_linter. (.data, type from ggplot2/rlang NSE)
    p = ggplot(plot_data, aes(x = .data[["grid"]], y = .data[["value"]],
        group = .data[["id"]], color = .data[["type"]]))
    if (!noline) {
      p = p + geom_line(alpha = 0.9, linewidth = 0.5, linetype = "dotted", na.rm = TRUE)
      p = p + geom_line(data = subset(plot_data, type == "PDP"), linewidth = 0.8)
    } else {
      p = p + geom_point(size = 1, shape = 4, na.rm = TRUE)
      p = p + geom_point(data = subset(plot_data, type == "PDP"), size = 3, shape = 4, na.rm = TRUE)
    }
    if (feat %in% colnames(origin_data_subset)) {
      p = p + geom_point(data = origin_data_subset,
        aes(x = .data[[feat]], y = .data[[target_feature_name]]),
        alpha = if (show_point) 0.3 else 0, size = 0.8, inherit.aes = FALSE)
    }
    # nolint end
    ylim_top = ymax
    if (show_point && is.finite(ymin) && is.finite(ymax)) {
      rng = ymax - ymin
      ylim_top = ymax + max(rng, 1e-6) * 0.12
    }
    p = p +
      scale_color_manual(values = c("ICE" = color_ice, "PDP" = color_pd),
        labels = c("ICE" = if (mean_center) "Mean centered ICE" else "ICE",
          "PDP" = if (mean_center) "Mean centered PDP" else "PDP")) +
      coord_cartesian(ylim = c(ymin, ylim_top)) +
      # ylim(ymin, ymax) +
      theme_bw(base_size = 9) +
      labs(
        x = if (!is.null(split_condition)) paste0(feat, " | ", split_condition) else feat,
        y = target_feature_name,
        title = "Partial Dependence Plot",
        color = NULL
      ) +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.6, "lines"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9)
      ) +
      if (show_point) {
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.98),
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill = "white", color = NA),
          legend.box.background = element_rect(fill = "white", color = "grey65", linewidth = 0.2),
          legend.margin = margin(3, 3, 3, 3)
        )
      } else {
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.95, 0.95),
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill = NA, color = NA),
          legend.box.background = element_rect(fill = NA, color = "grey", linewidth = 0.1)
        )
      }
    p
  })
  plot
}
