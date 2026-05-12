#' Node: Tree Node for Effect-based Decision Trees (R6 class)
#'
#' Represents a single node in an effect-based decision tree, storing split information,
#' effect statistics, and child nodes. Uses grouped fields for clearer structure.
#'
#' @field id (`integer(1)`) \cr
#'   Node identifier within its depth level.
#' @field depth (`integer(1)`) \cr
#'   Depth of the node (root starts at 1).
#' @field subset_idx (`integer()`) \cr
#'   Row indices of data that fall into this node.
#' @field grid (`list()`) \cr
#'   Grid values for each feature in this node.
#' @field parent (`list()` or `NULL`) \cr
#'   Parent info: id, child_type, split_feature, split_value, objective_value, int_imp. NULL for root.
#' @field split (`list()` or `NULL`) \cr
#'   Split info: feature, value. NULL for terminal nodes.
#' @field objective (`list()`) \cr
#'   Objective: value (scalar), value_j (per-feature vector).
#' @field importance (`list()` or `NULL`) \cr
#'   Importance: imp (scalar), imp_j (per-feature). NULL for root and unsplit nodes.
#' @field children (`list()` or `NULL`) \cr
#'   Left and right child nodes (or NULL for terminal nodes).
#' @field stop_criterion_met (`logical(1)`) \cr
#'   Whether the minimal node size or improvement threshold has been reached.
#' @field improvement_met (`logical(1)`) \cr
#'   Whether the improvement-based stop criterion was met (`TRUE` = stop splitting).
#' @field strategy (PdStrategy | AleStrategy) \cr
#'   Strategy for effect-specific operations.
#'
#' @details
#' This class is used internally by GadgetTree and strategy objects to represent
#' and manage nodes in effect-based decision trees. Each node stores split information,
#' effect statistics, and references to its children.
#'
#' @examples
#' # Example: Creating a Node (typically done internally)
#' # node = Node$new(id = 1, depth = 1, subset_idx = 1:100, grid = list(feature1 = 1:10))
#'
#' @importFrom R6 R6Class
#'
#' @keywords internal
Node = R6::R6Class("Node", public = list(
  id = NULL,
  depth = NULL,
  subset_idx = NULL,
  grid = NULL,
  parent = NULL,
  split = NULL,
  objective = NULL,
  importance = NULL,
  children = NULL,
  stop_criterion_met = NULL,
  improvement_met = NULL,
  strategy = NULL,

  #' @description
  #' Create a node from id, depth, subset indices, and grid. Assigns fields and sets
  #' \code{stop_criterion_met = FALSE}. Accepts legacy scalar args for compatibility.
  #' @param id (`integer(1)`) \cr
  #'   Node identifier.
  #' @param depth (`integer(1)` or `NULL`) \cr
  #'   Node depth (root is 1).
  #' @param subset_idx (`integer()`) \cr
  #'   Row indices of data in this node.
  #' @param grid (`list()`) \cr
  #'   Grid values for each feature.
  #' @param id_parent (`integer(1)` or `NULL`) \cr
  #'   Parent node id.
  #' @param child_type (`character(1)` or `NULL`) \cr
  #'   Split direction (\code{"<="}, \code{">"}, \code{"=="}, \code{"!="}, or \code{"in"}).
  #' @param objective_value_parent (`numeric(1)` or `NULL`) \cr
  #'   Parent node's objective value.
  #' @param objective_value_j (`numeric()` or `NULL`) \cr
  #'   Objective values per feature.
  #' @param objective_value (`numeric(1)` or `NULL`) \cr
  #'   Total objective value.
  #' @param improvement_met (`logical(1)`) \cr
  #'   Whether improvement threshold was met.
  #' @param int_imp (`numeric(1)` or `NULL`) \cr
  #'   Interaction importance.
  #' @param int_imp_j (`numeric()` or `NULL`) \cr
  #'   Interaction importance per feature.
  #' @param strategy (PdStrategy | AleStrategy or `NULL`) \cr
  #'   Strategy; \code{NULL} not used in practice.
  initialize = function(id, depth = NULL, subset_idx, grid, id_parent = NULL,
    child_type = NULL, objective_value_parent = NULL, objective_value_j = NULL,
    objective_value = NULL, improvement_met = FALSE, int_imp = NULL, int_imp_j = NULL, strategy = NULL) {

    checkmate::assert_numeric(id, len = 1)
    checkmate::assert_numeric(depth, len = 1, null.ok = TRUE)
    checkmate::assert_numeric(subset_idx, min.len = 1)
    checkmate::assert_numeric(id_parent, len = 1, null.ok = TRUE)
    checkmate::assert_character(child_type, null.ok = TRUE)

    self$id = id
    self$depth = depth
    self$subset_idx = subset_idx
    self$grid = grid
    self$parent = if (is.null(id_parent)) NULL else list(
      id = id_parent,
      child_type = child_type,
      split_feature = NULL,
      split_value = NULL,
      objective_value = objective_value_parent,
      int_imp = NULL
    )
    self$split = NULL
    self$objective = list(value = objective_value, value_j = objective_value_j)
    self$importance = if (is.null(int_imp) && is.null(int_imp_j)) NULL else list(imp = int_imp, imp_j = int_imp_j)
    self$stop_criterion_met = FALSE
    self$improvement_met = improvement_met
    self$strategy = strategy
  },

  #' @description
  #' Given Z, Y, root objectives, and tree params: checks stopping criteria;
  #' finds best split; creates and applies children; recurses into child nodes.
  #' Returns NULL if no valid split.
  #' @param Z (`data.frame()` or `data.table()`) \cr
  #'   Split features.
  #' @param Y (`list()`) \cr
  #'   Effect list.
  #' @param objective_value_root_j (`numeric()`) \cr
  #'   Root objective values per feature.
  #' @param objective_value_root (`numeric(1)`) \cr
  #'   Root total objective value.
  #' @param min_node_size (`integer(1)`) \cr
  #'   Minimum node size.
  #' @param n_quantiles (`integer(1)` or `NULL`) \cr
  #'   Quantiles for candidate split points.
  #' @param impr_par (`numeric(1)`) \cr
  #'   Improvement threshold.
  #' @param depth (`integer(1)`) \cr
  #'   Current node depth.
  #' @param max_depth (`integer(1)`) \cr
  #'   Maximum tree depth.
  #' @return (`NULL`)
  split_node = function(Z, Y, objective_value_root_j, objective_value_root,
    min_node_size, n_quantiles, impr_par, depth, max_depth) {
    t0 = proc.time()
    # 1. Stopping criteria
    if (objective_value_root < 1e-10 || depth >= max_depth ||
        length(self$subset_idx) < min_node_size || isTRUE(self$improvement_met)) {
      self$stop_criterion_met = TRUE
      # Recursion exit: stop splitting at this node
      return(NULL)
    }
    # 2. Find the best split
    split_info = tryCatch({
      y_curr = self$strategy$node_transform(
        Y = Y, idx = self$subset_idx, grid = self$grid,
        is_child = !is.null(self$parent)
      )
      self$find_best_split(Z, y_curr, min_node_size, n_quantiles)
    }, error = function(e) {
      cli::cli_warn("find_best_split error at node {self$id} (depth {self$depth}): {e$message}")
      NULL
    })
    if (is.null(split_info)) {
      self$stop_criterion_met = TRUE
      return(NULL)
    }
    # 3. Create left and right child nodes
    children_info = tryCatch({
      self$create_children(Z[[split_info$split_feature]], Y, split_info,
        objective_value_root_j, objective_value_root, impr_par)
    }, error = function(e) {
      cli::cli_warn("create_children error at node {self$id} (depth {self$depth}): {e$message}")
      NULL
    })
    if (is.null(children_info)) {
      self$stop_criterion_met = TRUE
      return(NULL)
    }
    # 4. Apply the split
    self$apply_split(split_info, children_info)
    t1 = proc.time()
    elapsed = as.numeric((t1 - t0)[3])
    # 5. Record split time
    if (!is.null(self$strategy) && !is.null(self$strategy$tree_ref)) {
      self$strategy$tree_ref$split_benchmark[[length(self$strategy$tree_ref$split_benchmark) + 1]] =
        list(node_id = self$id, depth = self$depth, time = elapsed)
    }
    # 6. Recursively split left and right child nodes
    if (!is.null(self$children$left_child)) {
      self$children$left_child$split_node(
        Z, Y,
        objective_value_root_j,
        objective_value_root,
        min_node_size,
        n_quantiles,
        impr_par,
        depth + 1,
        max_depth)
    }
    if (!is.null(self$children$right_child)) {
      self$children$right_child$split_node(
        Z, Y,
        objective_value_root_j,
        objective_value_root,
        min_node_size,
        n_quantiles,
        impr_par,
        depth + 1,
        max_depth)
    }
  },

  #' @description
  #' Given Z (subset by node indices), y_curr, and params: calls
  #' \code{strategy$find_best_split} and returns list with
  #' \code{split_feature}, \code{split_value}, \code{is_categorical}
  #' (and for AleStrategy: \code{left/right_objective_value_j}).
  #' @param Z (`data.frame()` or `data.table()`) \cr
  #'   Split features.
  #' @param y_curr (`list()`) \cr
  #'   Effect list for current node.
  #' @param min_node_size (`integer(1)`) \cr
  #'   Minimum node size.
  #' @param n_quantiles (`integer(1)` or `NULL`) \cr
  #'   Quantiles for candidate split points.
  #' @return (`list()` or `NULL`) \cr
  #'   Best split info or \code{NULL} if no valid split.
  find_best_split = function(Z, y_curr, min_node_size, n_quantiles) {
    z_subset = Z[self$subset_idx, ]
    split_res = self$strategy$find_best_split(Z = z_subset, Y = y_curr,
      min_node_size = min_node_size, n_quantiles = n_quantiles)
    if (is.null(split_res$best_split) || length(split_res$best_split) == 0 || all(!split_res$best_split)) {
      NULL
    } else {
      list(
        split_feature = split_res$split_feature[split_res$best_split][1],
        split_value = split_res$split_point[split_res$best_split][1],
        is_categorical = split_res$is_categorical[split_res$best_split][1],
        raw_result = split_res
      )
    }
  },

  #' @description
  #' Given Z, Y, split_info, and root objectives: computes idx_left/right,
  #' child grids, objective values; checks improvement threshold;
  #' creates left/right Node instances and sets parent info.
  #' Returns list of \code{left_child}, \code{right_child}, \code{int_imp},
  #' \code{int_imp_j} or NULL if improvement too small.
  #' @param z_split_feature (`numeric()`) \cr
  #'   Numeric vector. Values of the splitting feature of this split.
  #' @param Y (`list()`) \cr
  #'   Effect list.
  #' @param split_info (`list()`) \cr
  #'   Split information.
  #' @param objective_value_root_j (`numeric()`) \cr
  #'   Root objective values per feature.
  #' @param objective_value_root (`numeric(1)`) \cr
  #'   Root total objective value.
  #' @param impr_par (`numeric(1)`) \cr
  #'   Improvement threshold.
  #' @return (`list()`) \cr
  #'   Left/right child nodes and split statistics.
  create_children = function(z_split_feature, Y, split_info, objective_value_root_j, objective_value_root, impr_par) {
    split_feature = split_info$split_feature
    split_value = split_info$split_value
    is_categorical = split_info$is_categorical
    # Get indices for children
    z_sub = z_split_feature[self$subset_idx]
    is_ale_categorical = is_ale_ordered_categorical_split(is_categorical, self$strategy)
    split_groups = NULL
    if (is_ale_categorical) {
      split_groups = ordered_categorical_split_groups(z_split_feature, split_value)
      left_mask = ordered_categorical_left_mask(z_sub, split_value)
      idx_left = self$subset_idx[which(left_mask)]
      idx_right = self$subset_idx[which(!left_mask)]
    } else if (is_categorical) {
      split_groups = one_vs_rest_categorical_split_groups(z_split_feature, split_value)
      idx_left = self$subset_idx[which(z_sub == split_value)]
      idx_right = self$subset_idx[which(z_sub != split_value)]
    } else {
      idx_left = self$subset_idx[which(z_sub <= as.numeric(split_value))]
      idx_right = self$subset_idx[which(z_sub > as.numeric(split_value))]
    }
    if (length(idx_left) == 0 || length(idx_right) == 0) {
      return(NULL)
    }

    grid_info = self$create_child_grids(split_feature, split_value, is_categorical)
    obj = self$strategy$get_child_objectives(
      Z, Y, split_info, idx_left, idx_right,
      grid_info$grid_left, grid_info$grid_right
    )
    left_objective_value_j = obj$left_objective_value_j
    right_objective_value_j = obj$right_objective_value_j
    left_objective_value = obj$left_objective_value
    right_objective_value = obj$right_objective_value
    int_imp_j = (self$objective$value_j - left_objective_value_j - right_objective_value_j) / objective_value_root_j
    int_imp_j[!is.finite(int_imp_j)] = NA_real_
    int_imp = (self$objective$value - left_objective_value - right_objective_value) / objective_value_root

    # Threshold for root node: impr_par; for child node: parent int_imp * impr_par
    threshold = if (is.null(self$parent)) impr_par else self$parent$int_imp * impr_par
    # Check if improvement meets threshold
    if (int_imp < threshold) {
      self$improvement_met = TRUE
      return(NULL) # Improvement not sufficient: stop splitting at this node
    }
    # Create child nodes
    left_child = Node$new(
      id = 2 * self$id, depth = self$depth + 1,
      subset_idx = idx_left, grid = grid_info$grid_left, id_parent = self$id,
      child_type = if (is_categorical) {
        "in"
      } else {
        "<="
      },
      objective_value_parent = self$objective$value,
      objective_value = left_objective_value,
      objective_value_j = left_objective_value_j,
      int_imp = NULL, int_imp_j = NULL,
      improvement_met = self$improvement_met,
      strategy = self$strategy
    )
    right_child = Node$new(
      id = 2 * self$id + 1, depth = self$depth + 1,
      subset_idx = idx_right, grid = grid_info$grid_right, id_parent = self$id,
      child_type = if (is_categorical) {
        "in"
      } else {
        ">"
      },
      objective_value_parent = self$objective$value,
      objective_value = right_objective_value,
      objective_value_j = right_objective_value_j,
      int_imp = NULL, int_imp_j = NULL,
      improvement_met = self$improvement_met,
      strategy = self$strategy
    )
    # Set parent split/int_imp for children
    left_child$parent$split_feature = right_child$parent$split_feature = split_feature
    left_child$parent$split_value = right_child$parent$split_value = split_value
    left_child$parent$int_imp = right_child$parent$int_imp = int_imp
    if (is_categorical) {
      left_child$parent$split_levels = split_groups$left_levels
      right_child$parent$split_levels = split_groups$right_levels
      left_child$parent$split_condition =
        format_categorical_split_condition(split_feature, split_groups$left_levels)
      right_child$parent$split_condition =
        format_categorical_split_condition(split_feature, split_groups$right_levels)
    }

    list(
      left_child = left_child,
      right_child = right_child,
      int_imp = int_imp,
      int_imp_j = int_imp_j
    )
  },

  #' @description
  #' Given split_feature, split_value, and is_categorical: partitions
  #' \code{self$grid[[split_feature]]} into left (<= or ==) and
  #' right (> or !=). Returns list \code{grid_left}, \code{grid_right}.
  #' @param split_feature (`character(1)`) \cr
  #'   Feature used for splitting.
  #' @param split_value (`numeric(1)` or `factor()`) \cr
  #'   Split value.
  #' @param is_categorical (`logical(1)`) \cr
  #'   Whether the split feature is categorical.
  #' @return (`list()`) \cr
  #'   \code{grid_left}, \code{grid_right}.
  create_child_grids = function(split_feature, split_value, is_categorical) {
    grid_left = self$grid
    grid_right = self$grid
    if (split_feature %in% names(self$grid) && length(self$grid[[split_feature]]) > 0L) {
      if (is_ale_ordered_categorical_split(is_categorical, self$strategy)) {
        grid_left_idx = ordered_categorical_left_mask(grid_left[[split_feature]], split_value)
        grid_right_idx = !grid_left_idx
      } else if (is_categorical) {
        grid_left_idx = grid_left[[split_feature]] == split_value
        grid_right_idx = grid_right[[split_feature]] != split_value
      } else {
        grid_left_idx = as.numeric(grid_left[[split_feature]]) <= as.numeric(split_value)
        grid_right_idx = as.numeric(grid_right[[split_feature]]) > as.numeric(split_value)
      }
      grid_left[[split_feature]] = grid_left[[split_feature]][grid_left_idx]
      grid_right[[split_feature]] = grid_right[[split_feature]][grid_right_idx]
    }
    list(grid_left = grid_left, grid_right = grid_right)
  },

  #' @description
  #' Given split_info and children_info: sets \code{split_feature},
  #' \code{split_value}, \code{int_imp}, \code{int_imp_j}, \code{children}.
  #' @param split_info (`list()`) \cr
  #'   Split information.
  #' @param children_info (`list()`) \cr
  #'   Children information.
  #' @return (`NULL`)
  apply_split = function(split_info, children_info) {
    self$split = list(
      feature = split_info$split_feature,
      value = if (split_info$is_categorical) split_info$split_value else as.numeric(split_info$split_value)
    )
    self$importance = list(imp = children_info$int_imp, imp_j = children_info$int_imp_j)
    self$children = list("left_child" = children_info$left_child, "right_child" = children_info$right_child)
  }
))
