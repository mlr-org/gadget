# Format split threshold for path labels: numeric rounded; factors/character as label (no as.numeric).
format_split_condition_value = function(v) {
  if (is.numeric(v)) {
    as.character(round(v, 3L))
  } else {
    as.character(v)
  }
}

format_parent_split_condition = function(parent_node, current_node) {
  if (!is.null(current_node$parent$split_condition)) {
    return(current_node$parent$split_condition)
  }
  op = choose_operator(parent_node, current_node)
  val_txt = format_split_condition_value(parent_node$split$value)
  paste0(parent_node$split$feature, " ", op, " ", val_txt)
}

#' Build path of split conditions from root to node
#'
#' Given node and tree (depth-list): walks parent chain via \code{find_node_by_id};
#' at each step builds condition string (e.g. "x <= 0.5") via \code{choose_operator}.
#' Returns character vector of conditions (root to node).
#'
#' @param node (`Node`) \cr
#'   Node object with \code{parent} and \code{depth}.
#' @param tree (`list()`) \cr
#'   Depth-based list of nodes (from \code{convert_tree_to_list}).
#' @return (`character()`) \cr
#'   Conditions from root to node (e.g. \code{"x <= 0.5"}).
#' @keywords internal
track_split_condition = function(node, tree) {
  path_conditions = character(0)
  current_node = node
  while (!is.null(current_node$parent)) {
    parent_node = find_node_by_id(tree[[current_node$depth - 1]], current_node$parent$id)
    if (is.null(parent_node)) break
    cond = format_parent_split_condition(parent_node, current_node)
    path_conditions = c(cond, path_conditions)
    current_node = parent_node
  }
  path_conditions
}
