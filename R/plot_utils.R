#' Build node title for ALE or PD tree plots.
#'
#' @param node (`Node`) \cr
#'   Node object.
#' @param depth_idx (`integer(1)`) \cr
#'   Depth index.
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param style (`character(1)`) \cr
#'   \code{"ale"} or \code{"pd"} for format variant.
#'
#' @return (`character(1)`) \cr
#'   Title string.
#' @keywords internal
build_node_title = function(node, depth_idx, tree, style = c("ale", "pd")) {
  style = match.arg(style)
  n_samples = length(node$subset_idx)
  if (depth_idx == 1) {
    return(paste0("Root node (N = ", n_samples, ")"))
  }
  path_conditions = track_split_condition(node, tree)
  cond_txt = if (length(path_conditions) > 0) paste(path_conditions, collapse = " & ") else NULL
  if (style == "ale") {
    cond_txt = if (is.null(cond_txt) || length(cond_txt) == 0) "Split" else cond_txt
    paste0("Depth ", depth_idx, " - Node ", node$id, " (", cond_txt, ", N = ", n_samples, ")")
  } else {
    cond_txt = if (is.null(cond_txt) || length(cond_txt) == 0) "" else cond_txt
    paste0(depth_idx - 1, ".Split results: ", cond_txt, " (N = ", n_samples, ")")
  }
}

#' Select depth indices to render based on tree, depth, and node_id.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param depth (`integer()` or `NULL`) \cr
#'   Depth levels to consider; \code{NULL} = all.
#' @param node_id (`integer()` or `NULL`) \cr
#'   Node IDs to render; \code{NULL} = all nodes at selected depths.
#'
#' @return (`integer()`) \cr
#'   Depth indices to render.
#' @keywords internal
select_depths_to_render = function(tree, depth, node_id) {
  if (is.null(node_id)) {
    if (is.null(depth)) return(seq_along(tree))
    depths = suppressWarnings(as.integer(depth))
    return(depths[!is.na(depths) & depths > 0 & depths <= length(tree)])
  }
  ids = suppressWarnings(as.integer(node_id))
  depths_hit = integer(0)
  for (depth_idx in seq_along(tree)) {
    nodes = tree[[depth_idx]]
    if (any(mlr3misc::map_lgl(nodes, function(n) !is.null(n) && n$id %in% ids))) {
      depths_hit = c(depths_hit, depth_idx)
    }
  }
  if (!is.null(depth)) {
    depth_filter = suppressWarnings(as.integer(depth))
    depth_filter = depth_filter[!is.na(depth_filter) & depth_filter > 0 & depth_filter <= length(tree)]
    depths_hit = intersect(depths_hit, depth_filter)
  }
  unique(depths_hit)
}

#' Shared implementation for y-axis range from effect values.
#'
#' @param effect_values (`numeric()`) \cr
#'   Effect values (no NA).
#' @param data (`data.frame()` or `NULL`) \cr
#'   Original data; target values are always included in range when available.
#' @param target_feature_name (`character(1)` or `NULL`) \cr
#'   Target column.
#' @return (`list()`) \cr
#'   \code{ymin}, \code{ymax}.
#' @keywords internal
calculate_y_range_impl = function(effect_values, data = NULL, target_feature_name = NULL) {
  effect_values = effect_values[!is.na(effect_values)]
  if (!length(effect_values)) return(list(ymin = 0, ymax = 1))
  ymin = min(effect_values)
  ymax = max(effect_values)
  if (!is.null(data) && !is.null(target_feature_name) && target_feature_name %in% colnames(data)) {
    ymin = min(ymin, min(data[[target_feature_name]], na.rm = TRUE))
    ymax = max(ymax, max(data[[target_feature_name]], na.rm = TRUE))
  }
  rng = ymax - ymin
  list(ymin = ymin, ymax = ymax + 0.2 * rng)
}

#' Compute y-axis range for ALE plots.
#'
#' @param curves (`list()`) \cr
#'   Output of \code{prepare_plot_data_ale}.
#' @param data (`data.frame()` or `NULL`) \cr
#'   Unused; kept for call compatibility. Limits follow \code{d_l} only, not the response column.
#' @param target_feature_name (`character(1)` or `NULL`) \cr
#'   Unused; see \code{data}.
#'
#' @return (`list()`) \cr
#'   \code{ymin}, \code{ymax}.
#' @keywords internal
calculate_y_range_ale = function(curves, data, target_feature_name) {
  effect_values = unlist(mlr3misc::map(curves, function(x) x$mean_effect$d_l), use.names = FALSE)
  calculate_y_range_impl(effect_values, NULL, NULL)
}

# Y limits for regional ALE panels: global + node curves (node-only can lie outside global after subsetting).
calculate_y_range_ale_combined = function(global_curves, regional_curves, data, target_feature_name) {
  gv = unlist(mlr3misc::map(global_curves, function(x) x$mean_effect$d_l), use.names = FALSE)
  rv = unlist(mlr3misc::map(regional_curves, function(x) x$mean_effect$d_l), use.names = FALSE)
  calculate_y_range_impl(c(gv, rv), NULL, NULL)
}

# After effect-only ylim (e.g. mean-centered ICE/PDP or ALE \code{d_l}), expand to include finite overlaid
# response values when \code{show_point = TRUE} (\code{coord_cartesian} clips otherwise).
merge_ale_y_range_with_response = function(y_range, y_obs) {
  checkmate::assert_list(y_range, names = "named", .var.name = "y_range")
  checkmate::assert_numeric(y_obs, .var.name = "y_obs")
  y_obs = y_obs[is.finite(y_obs)]
  if (!length(y_obs)) {
    return(y_range)
  }
  ymin = min(y_range$ymin, min(y_obs))
  ymax = max(y_range$ymax, max(y_obs))
  rng = ymax - ymin
  list(ymin = ymin, ymax = ymax + 0.2 * rng)
}

#' Preprocess PD node data by depth.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param Y,grid_total (`list()`) \cr
#'   Effect matrices and grid.
#' @param mean_center (`logical(1)`) \cr
#'   Whether to mean-center.
#'
#' @return (`list()`) \cr
#'   Processed data per depth.
#' @keywords internal
preprocess_node_data = function(tree, Y, grid_total, mean_center) {
  all_node_data = list()
  for (depth_idx in seq_along(tree)) {
    nodes = tree[[depth_idx]]
    y_processed = mlr3misc::map(names(Y), function(i) {
      y_i = Y[[i]]
      y_i$node = NA_integer_
      grid_i_total = grid_total[[i]]
      for (node_idx in seq_along(nodes)) {
        node = nodes[[node_idx]]
        if (!is.null(node)) {
          subset_idx = node$subset_idx
          grid_i_curr = node$grid[[i]]
          y_i$node[subset_idx] = node_idx

          if (length(grid_i_curr) < length(grid_i_total)) {
            y_i[subset_idx, which(!(grid_i_total %in% grid_i_curr))] = NA
            if (mean_center && length(grid_i_curr) > 1) {
              y_i[subset_idx, grid_i_total] = y_i[subset_idx, grid_i_total] -
                rowMeans(y_i[subset_idx, grid_i_total], na.rm = TRUE)
            }
          }
        }
      }
      y_i
    })
    names(y_processed) = names(Y)
    all_node_data[[depth_idx]] = y_processed
  }
  all_node_data
}

#' Create PD plots for one depth.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param prepared_data (`list()`) \cr
#'   Prepared effect data per feature.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Original data.
#' @param target_feature_name (`character(1)`) \cr
#'   Target column name.
#' @param depth_idx (`integer(1)`) \cr
#'   Depth level.
#' @param nodes_to_render (`integer()`) \cr
#'   Node indices to render.
#' @param color_ice,color_pd (`character(1)`) \cr
#'   Colors for ICE and PDP.
#' @param show_plot,show_point,mean_center (`logical(1)`) \cr
#'   Plot options.
#'
#' @return (`list()`) \cr
#'   Named list of patchwork objects per node.
#' @keywords internal
create_plots_for_depth = function(tree, prepared_data, data, target_feature_name,
  depth_idx, nodes_to_render, color_ice, color_pd,
  show_plot, show_point, mean_center) {
  plots_at_depth = list()
  for (node_idx in nodes_to_render) {
    node = tree[[depth_idx]][[node_idx]]
    if (!is.null(node)) {
      path_conditions = track_split_condition(node, tree)
      split_condition = if (length(path_conditions) > 0) paste(path_conditions, collapse = " & ") else NULL
      title = build_node_title(node, depth_idx, tree, style = "pd")
      y_range = calculate_y_range(prepared_data, data, target_feature_name, mean_center = mean_center)
      if (isTRUE(mean_center) && isTRUE(show_point) && length(node$subset_idx)) {
        y_node = data[[target_feature_name]][node$subset_idx]
        y_range = merge_ale_y_range_with_response(y_range, y_node)
      }

      plots = plot_regional_pd(prepared_data = prepared_data,
        origin_data = data,
        target_feature_name = target_feature_name,
        node_idx = node_idx,
        color_ice = color_ice,
        color_pd = color_pd,
        ymin = y_range$ymin,
        ymax = y_range$ymax,
        split_condition = split_condition,
        show_point = show_point,
        mean_center = mean_center)

      p = patchwork::wrap_plots(plots, ncol = min(length(prepared_data), 3)) +
        patchwork::plot_annotation(title = title) &
        theme(plot.title = element_text(hjust = 0.5))
      if (show_plot) print(p)
      plots_at_depth[[paste0("Node_", node_idx)]] = p
    }
  }
  plots_at_depth
}

#' Compute y-axis range for PD plots.
#'
#' @param prepared_data (`list()`) \cr
#'   Prepared effect data.
#' @param data (`data.frame()`) \cr
#'   Original data.
#' @param target_feature_name (`character(1)`) \cr
#'   Target column.
#' @param mean_center (`logical(1)`) \cr
#'   If \code{TRUE}, ICE/PDP are on a centered scale; the raw target column is not merged into
#'   \code{ylim}. If \code{FALSE}, limits follow ICE/PD values and the target range (for overlays).
#'
#' @return (`list()`) \cr
#'   \code{ymin}, \code{ymax}.
#' @keywords internal
calculate_y_range = function(prepared_data, data, target_feature_name, mean_center = FALSE) {
  checkmate::assert_flag(mean_center, .var.name = "mean_center")
  values_list = mlr3misc::map(prepared_data, function(df) {
    if (is.data.frame(df) && ncol(df) >= 1) {
      cols = setdiff(colnames(df), "node")
      if (length(cols) == 0) return(numeric(0))
      as.numeric(unlist(df[, cols, drop = FALSE]))
    } else {
      numeric(0)
    }
  })
  effect_values = unlist(values_list, use.names = FALSE)
  if (isTRUE(mean_center)) {
    calculate_y_range_impl(effect_values, NULL, NULL)
  } else {
    calculate_y_range_impl(effect_values, data, target_feature_name)
  }
}
