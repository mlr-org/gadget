#' Visualize the Tree Structure
#'
#' Given tree (depth-list of Node objects): calls
#' \code{prepare_layout_data} to build layout data; creates parent map;
#' builds edge list; creates ggraph plot with nodes labeled by split info
#' and edges representing tree hierarchy. Returns ggplot object.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param label_wrap_width (`integer(1)` or `NULL`) \cr
#'   If not \code{NULL}, wrap each line of node labels to this many characters (see \code{strwrap}).
#' @param node_spread_x,node_spread_y (`numeric(1)`) \cr
#'   Positive multipliers applied to the default \code{"tree"} layout coordinates to separate nodes.
#'
#' @return (ggplot) \cr
#'   Tree structure visualization.
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph create_layout ggraph geom_edge_elbow geom_node_label circle
#' @importFrom ggplot2 aes coord_flip scale_fill_manual theme_void scale_y_reverse theme expansion arrow unit margin
#' @importFrom stats setNames na.omit
#' @importFrom grDevices hcl.colors
#'
#' @keywords internal
plot_tree_structure = function(tree, label_wrap_width = 34L, node_spread_x = 1.55, node_spread_y = 1.12) {
  checkmate::assert_integerish(label_wrap_width, len = 1L, lower = 8L, null.ok = TRUE, .var.name = "label_wrap_width")
  checkmate::assert_numeric(node_spread_x, len = 1L, lower = 0.5, finite = TRUE, .var.name = "node_spread_x")
  checkmate::assert_numeric(node_spread_y, len = 1L, lower = 0.5, finite = TRUE, .var.name = "node_spread_y")
  data = prepare_layout_data(tree)
  if (!is.null(label_wrap_width)) {
    data$label = vapply(
      data$label,
      wrap_tree_label,
      FUN.VALUE = character(1L),
      width = as.integer(label_wrap_width)
    )
  }
  parent_map = setNames(data$id, data$node_id)
  data$parent_id = parent_map[as.character(data$id_parent)]
  edge_list = na.omit(data[, c("parent_id", "id")])
  colnames(edge_list) = c("from", "to")

  g = igraph::graph_from_data_frame(edge_list, vertices = data, directed = TRUE)

  lay = ggraph::create_layout(g, layout = "tree")
  lay$x = lay$x * node_spread_x
  lay$y = lay$y * node_spread_y

  gg = ggraph::ggraph(lay) +
    coord_flip(clip = "off")

  # Only add edges if there are any
  if (nrow(edge_list) > 0) {
    gg = gg + ggraph::geom_edge_elbow(
      arrow = arrow(length = unit(0.05, "cm")),
      end_cap = ggraph::circle(1.5, "mm"),
      edge_colour = "grey40",
      edge_width = 0.4
    )
  }

  gg +
    ggraph::geom_node_label(
      aes(label = label, fill = factor(depth)),
      size = 3.5,
      label.padding = unit(0.25, "lines"),
      # label.size = 0.3,
      label.r = unit(0.1, "lines")
    ) +
    scale_fill_manual(values = hcl.colors(n = length(tree), palette = "Set2")) +
    theme_void() +
    ggplot2::scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_reverse(expand = expansion(mult = c(0.08, 0.08))) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 14, r = 28, b = 14, l = 28, unit = "mm")
    )
}
