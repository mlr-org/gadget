#' Plot ALE curves for selected tree nodes
#'
#' External helper mirroring \code{plot_tree_pd}, used by \code{AleStrategy$plot}.
#' Produces patchwork objects per node with ALE mean curves and optional observation points.
#' Y-axis limits use global and node-specific cumulative ALE (\code{d_l}). When
#' \code{show_point = TRUE}, the node subset response is merged into \code{ylim} so overlaid
#' points are not clipped by \code{coord_cartesian}. \code{x} limits follow the global grid.
#' Y-axis labels use \code{target_feature_name}, matching PD tree plots.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param effect (`list()`) \cr
#'   Output of \code{calculate_ale()}.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Data with features and target.
#' @param target_feature_name (`character(1)`) \cr
#'   Target variable name.
#' @param depth (`integer()` or `NULL`) \cr
#'   Depth levels to render.
#' @param node_id (`integer()` or `NULL`) \cr
#'   Node IDs to render.
#' @param features (`character()` or `NULL`) \cr
#'   Feature subset.
#' @param color_ale (`character(1)`) \cr
#'   Color for ALE curves.
#' @param show_plot,show_point,mean_center (`logical(1)`) \cr
#'   Plot options.
#' @param ... Additional arguments passed to plotting helpers.
#'
#' @return (`list()`) \cr
#'   Nested list (depth -> node -> patchwork).
#' @keywords internal
plot_tree_ale = function(tree, effect, data, target_feature_name,
  depth = NULL, node_id = NULL, features = NULL,
  color_ale = "lightgreen", show_plot = TRUE,
  show_point = TRUE, mean_center = TRUE, ...) {
  checkmate::assert_list(tree, .var.name = "tree")
  checkmate::assert_true(is.list(effect) || inherits(effect, "R6"), .var.name = "effect")
  checkmate::assert_data_frame(data, .var.name = "data")
  checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
  checkmate::assert_integerish(depth, lower = 1, null.ok = TRUE, .var.name = "depth")
  checkmate::assert_integerish(node_id, lower = 1, null.ok = TRUE, .var.name = "node_id")
  checkmate::assert_character(features, null.ok = TRUE, .var.name = "features")
  checkmate::assert_character(color_ale, len = 1, .var.name = "color_ale")
  checkmate::assert_flag(show_plot, .var.name = "show_plot")
  checkmate::assert_flag(show_point, .var.name = "show_point")
  checkmate::assert_flag(mean_center, .var.name = "mean_center")

  if (is.null(features)) features = names(effect)
  features = intersect(features, names(effect))
  if (!length(features)) cli::cli_abort("No overlapping features between 'features' and effect list.")

  depths_to_render = select_depths_to_render(tree, depth, node_id)
  if (!length(depths_to_render)) {
    cli::cli_warn("No valid depths to render for ALE plots.")
    return(invisible(list()))
  }

  global_curves = prepare_plot_data_ale(effect,
    features = features, mean_center = mean_center)
  x_limits = mlr3misc::map(setNames(nm = features), function(feat) {
    if (is.null(global_curves[[feat]])) return(NULL)
    mean_dt = global_curves[[feat]]$mean_effect
    if (!("x_grid" %in% names(mean_dt))) return(NULL)
    if (is.numeric(mean_dt$x_grid)) {
      range(mean_dt$x_grid, na.rm = TRUE)
    } else if (is.factor(mean_dt$x_grid)) {
      levels(mean_dt$x_grid)
    } else {
      NULL
    }
  })

  plot_list = list()
  for (depth_idx in depths_to_render) {
    nodes = tree[[depth_idx]]
    if (is.null(node_id)) {
      nodes_to_render = seq_along(nodes)
    } else {
      ids = suppressWarnings(as.integer(node_id))
      nodes_to_render = which(mlr3misc::map_lgl(nodes, function(n) n$id %in% ids))
    }
    depth_plots = list()
    for (node_idx in nodes_to_render) {
      node = nodes[[node_idx]]
      curves = prepare_plot_data_ale(effect,
        idx = node$subset_idx, features = features, mean_center = mean_center)
      y_range = calculate_y_range_ale_combined(global_curves, curves, data, target_feature_name)
      point_values = NULL
      if (show_point) {
        y_node = data[[target_feature_name]][node$subset_idx]
        y_range = merge_ale_y_range_with_response(y_range, y_node)
        point_values = mlr3misc::map(setNames(nm = names(curves)), function(feat) {
          data.frame(
            x = data[[feat]][node$subset_idx],
            y = y_node
          )
        })
      }
      feat_panels = plot_regional_ale(curves,
        color_ale = color_ale,
        target_feature_name = target_feature_name,
        mean_center = mean_center,
        ymin = y_range$ymin, ymax = y_range$ymax,
        show_point = show_point,
        point_values = point_values,
        x_limits = x_limits)
      if (!length(feat_panels)) next
      stacked = mlr3misc::map(names(feat_panels), function(feat_name) {
        panel = feat_panels[[feat_name]]
        panel +
          patchwork::plot_annotation(title = feat_name) &
          ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = 0.5))
      })
      if (!length(stacked)) next
      ncol_used = max(1L, min(length(stacked), 3L))
      nrow_used = ceiling(length(stacked) / ncol_used)
      combined = patchwork::wrap_plots(stacked, ncol = ncol_used, nrow = nrow_used) +
        patchwork::plot_annotation(title = build_node_title(node, depth_idx, tree, style = "ale")) &
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      combined$layout$ncol = ncol_used
      combined$layout$nrow = nrow_used
      if (show_plot) print(combined)
      depth_plots[[paste0("Node_", node$id)]] = combined
    }
    plot_list[[paste0("Depth_", depth_idx)]] = depth_plots
  }
  invisible(plot_list)
}
