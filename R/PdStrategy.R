#' PdStrategy: Generalized Additive Decomposition Based on PD Effects
#'
#' @description
#' PD-based effect strategy (inherits from \link{EffectStrategy}). Given effect or model and data,
#' preprocesses to Z/Y/grid;
#' mean-centers effects per node; computes sum-of-variances heterogeneity;
#' finds best split via C++; fits tree and plots PD/ICE.
#' Character feature columns are coerced to \code{factor} before ICE/PD computation so they match
#' split-matrix treatment and learner conventions (same as \code{prepare_split_data_common}).
#'
#' @usage NULL
#' @format [R6::R6Class] object inheriting from [EffectStrategy].
#'
#' @section Construction:
#' ```
#' s = PdStrategy$new()
#' ```
#'
#' @field effect (`list()` or `R6` or `NULL`) \cr
#'   Cached PD/ICE effect used when \code{$plot()} omits \code{effect}.
#'
#' @details
#' This class is used internally by the GadgetTree framework to implement partial dependence
#' tree growing, splitting, and visualization. It is not intended to be used directly by end users,
#' but can be instantiated for advanced customization.
#'
#' @examples
#' \dontrun{
#' # Example: Fit and plot a PD tree using PdStrategy and GadgetTree
#' # (Assuming effect and data are prepared)
#' pd_strat = PdStrategy$new()
#' tree = GadgetTree$new(strategy = pd_strat, n_split = 2)
#' tree$fit(data = data, target_feature_name = "target", effect = effect)
#' tree$plot(data = data, target_feature_name = "target", effect = effect)
#' }
#'
#' @include EffectStrategy.R
#' @export
PdStrategy = R6::R6Class(
  "PdStrategy",
  inherit = EffectStrategy,
  public = list(
    effect = NULL,

    #' @description
    #' Create a PdStrategy instance (calls \code{super$initialize("pd")}).
    initialize = function() {
      super$initialize("pd")
    },

    #' @description
    #' Preprocess to Z, Y, grid via \code{prepare_split_data_pd}.
    #' @param effect (R6 or `list()`) \cr
    #'   Effect object (e.g. FeatureEffect).
    #' @param data (`data.frame()` or `data.table()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)` or `NULL`) \cr
    #'   Target variable name.
    #' @param feature_set (`character()` or `NULL`) \cr
    #'   Features for effect; \code{NULL} = all.
    #' @param split_feature (`character()` or `NULL`) \cr
    #'   Features for splitting; \code{NULL} = all.
    #' @return (`list()`) \cr
    #'   \code{Z}, \code{Y}, \code{grid}.
    preprocess = function(effect, data, target_feature_name = NULL, feature_set = NULL,
      split_feature = NULL) {
      checkmate::assert_data_frame(data, .var.name = "data")
      checkmate::assert_character(target_feature_name, len = 1, null.ok = TRUE, .var.name = "target_feature_name")
      if (!is.null(target_feature_name)) {
        checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      }
      checkmate::assert_character(feature_set, null.ok = TRUE, .var.name = "feature_set")
      checkmate::assert_character(split_feature, null.ok = TRUE, .var.name = "split_feature")
      checkmate::assert_true(is.null(effect) || is.list(effect) || inherits(effect, "R6"), .var.name = "effect")
      prepare_split_data_pd(effect = effect, data = data, target_feature_name = target_feature_name,
        feature_set = feature_set, split_feature = split_feature)
    },

    #' @description
    #' Subset and mean-center via \code{re_mean_center_ice_cpp}.
    #' @param Y (`list()`) \cr
    #'   Effect matrices per feature.
    #' @param idx (`integer()`) \cr
    #'   Sample indices in the node.
    #' @param grid (`list()`) \cr
    #'   Feature grids; required for PD.
    #' @return (`list()`) \cr
    #'   Mean-centered effect matrices.
    node_transform = function(Y, idx, grid, is_child = FALSE) {
      checkmate::assert_list(Y, .var.name = "Y")
      checkmate::assert_list(grid, .var.name = "grid")
      checkmate::assert_integerish(idx, min.len = 1, .var.name = "idx")
      as_numeric_matrix = function(x) {
        x = as.matrix(x)
        storage.mode(x) = "double"
        x
      }
      # Normalize each effect block before passing to C++.
      y_numeric = lapply(Y, as_numeric_matrix)
      re_mean_center_ice_cpp(Y = y_numeric, grid = grid, idx = idx)
    },

    #' @description
    #' Compute heterogeneity via \code{node_heterogeneity}.
    #' @param Y (`list()`) \cr
    #'   Effect matrices.
    #' @return (`numeric()`) \cr
    #'   Heterogeneity per feature.
    heterogeneity = function(Y) {
      checkmate::assert_list(Y, .var.name = "Y")
      node_heterogeneity(Y)
    },

    #' @description
    #' Compute left/right child objective values via node_transform and heterogeneity.
    #' @param Z,Y,split_info,idx_left,idx_right,grid_left,grid_right \cr
    #'   Node/split context.
    #' @return (`list()`) \cr
    #'   \code{left_objective_value_j}, \code{right_objective_value_j},
    #'   \code{left_objective_value}, \code{right_objective_value}.
    get_child_objectives = function(Z, Y, split_info, idx_left, idx_right, grid_left, grid_right) {
      checkmate::assert_integerish(idx_left, min.len = 1, .var.name = "idx_left")
      checkmate::assert_integerish(idx_right, min.len = 1, .var.name = "idx_right")
      checkmate::assert_list(grid_left, .var.name = "grid_left")
      checkmate::assert_list(grid_right, .var.name = "grid_right")
      y_left = self$node_transform(Y = Y, idx = idx_left, grid = grid_left)
      y_right = self$node_transform(Y = Y, idx = idx_right, grid = grid_right)
      left_objective_value_j = self$heterogeneity(y_left)
      right_objective_value_j = self$heterogeneity(y_right)
      list(
        left_objective_value_j = left_objective_value_j,
        right_objective_value_j = right_objective_value_j,
        left_objective_value = sum(left_objective_value_j, na.rm = TRUE),
        right_objective_value = sum(right_objective_value_j, na.rm = TRUE)
      )
    },

    #' @description
    #' Find best split via \code{search_best_split_cpp}.
    #' @param Z (`data.frame()` or `data.table()`) \cr
    #'   Split features.
    #' @param Y (`list()`) \cr
    #'   Effect matrices.
    #' @param min_node_size (`integer(1)`) \cr
    #'   Minimum node size.
    #' @param n_quantiles (`integer(1)` or `NULL`) \cr
    #'   Quantile candidates.
    #' @return (`data.frame()` or `list()`) \cr
    #'   Best split info.
    find_best_split = function(Z, Y, min_node_size, n_quantiles) {
      checkmate::assert_true(data.table::is.data.table(Z) || is.data.frame(Z), .var.name = "Z")
      checkmate::assert_list(Y)
      checkmate::assert_integerish(min_node_size, len = 1, lower = 1, any.missing = FALSE, .var.name = "min_node_size")
      checkmate::assert_integerish(n_quantiles, len = 1, lower = 1, null.ok = TRUE, .var.name = "n_quantiles")
      search_best_split_cpp(Z = Z, Y = Y, min_node_size = min_node_size, n_quantiles = n_quantiles)
    },

    #' @description
    #' Plot PD/ICE tree via \code{plot_tree_pd}.
    #' @param tree (`list()`) \cr
    #'   Depth-based list of Node objects.
    #' @param effect (R6 or `list()` or `NULL`) \cr
    #'   Effect object.
    #' @param data (`data.frame()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target name.
    #' @param depth (`integer()` or `NULL`) \cr
    #'   Depths to plot.
    #' @param node_id (`integer()` or `NULL`) \cr
    #'   Node IDs to plot.
    #' @param features (`character()` or `NULL`) \cr
    #'   Features to plot.
    #' @param ... Plot arguments.
    #' @return (`list()`) \cr
    #'   Nested list (depth -> node -> patchwork).
    plot = function(
      tree, effect = NULL, data, target_feature_name, depth = NULL,
      node_id = NULL, features = NULL, ...
    ) {
      checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      if (is.null(effect)) {
        effect = self$effect
        if (is.null(effect)) {
          cli::cli_abort("No cached PD effect found. Please pass {.arg effect} or run fit() with {.arg model}.")
        }
      }
      plot_tree_pd(tree = tree, effect = effect, data = data,
        target_feature_name = target_feature_name,
        depth = depth, node_id = node_id, features = features, ...)
    },

    #' @description
    #' Fit tree: preprocess, create root, split recursively.
    #' @param tree (`GadgetTree`) \cr
    #'   Tree instance.
    #' @param effect (R6 or `list()` or `NULL`) \cr
    #'   Optional precomputed effect object.
    #' @param model (`any`) \cr
    #'   Fitted model for internal PD/ICE computation.
    #' @param data (`data.frame()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target name.
    #' @param feature_set,split_feature (`character()` or `NULL`) \cr
    #'   Feature subsets.
    #' @param predict_fun (`function()` or `NULL`) \cr
    #'   Optional prediction function.
    #' @param n_grid (`integer(1)`) \cr
    #'   Number of grid points for numeric features.
    #' @param pd_engine (`character(1)`) \cr
    #'   When computing ICE/PD from \code{model}: \code{"cpp"} (column-wise stacked \code{newdata},
    #'   xplaineff-style) or \code{"r"} (\code{data.table::rbindlist}).
    #' @param ... Ignored.
    #' @return (`GadgetTree`) \cr
    #'   The tree, invisibly.
    fit = function(tree, effect = NULL, model = NULL, data, target_feature_name,
      feature_set = NULL, split_feature = NULL, predict_fun = NULL,
      n_grid = 20L, pd_engine = c("cpp", "r"), ...) {
      checkmate::assert_r6(tree, classes = "GadgetTree", .var.name = "tree")
      checkmate::assert_data_frame(data, .var.name = "data")
      checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
      checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      checkmate::assert_character(feature_set, null.ok = TRUE, .var.name = "feature_set")
      checkmate::assert_character(split_feature, null.ok = TRUE, .var.name = "split_feature")
      checkmate::assert_true(
        is.null(effect) || is.list(effect) || inherits(effect, "R6"),
        .var.name = "effect"
      )
      checkmate::assert_function(predict_fun, null.ok = TRUE, .var.name = "predict_fun")
      checkmate::assert_integerish(n_grid, len = 1L, lower = 2L, .var.name = "n_grid")
      pd_engine = match.arg(pd_engine)

      # After checks: same coercion as prepare_split_data_common (ICE + Z use factor categoricals).
      feat_cols = setdiff(colnames(data), target_feature_name)
      data = ensure_factors(data, feat_cols)

      if (is.null(effect) && is.null(model)) {
        cli::cli_abort("PdStrategy requires either {.arg effect} or {.arg model}.")
      }
      if (!is.null(effect) && !is.null(model)) {
        cli::cli_warn("Both {.arg effect} and {.arg model} were provided; {.arg effect} takes precedence.")
      }
      # --- global part (timed): ICE/PD precompute when needed, then preprocess to Z/Y/grid ---
      t_global = system.time({
        if (is.null(effect)) {
          effect = calculate_pd(
            model = model,
            data = data,
            target_feature_name = target_feature_name,
            feature_set = feature_set,
            predict_fun = predict_fun,
            n_grid = n_grid,
            pd_engine = pd_engine
          )
        }

        self$tree_ref = tree
        self$effect = effect
        prepared_data = self$preprocess(effect = effect, data = data,
          target_feature_name = target_feature_name,
          feature_set = feature_set,
          split_feature = split_feature)
        Z = prepared_data$Z
        Y = prepared_data$Y
        grid = prepared_data$grid
        objective_value_root_j = self$heterogeneity(Y)
        objective_value_root = sum(objective_value_root_j, na.rm = TRUE)
      })[["elapsed"]]

      t_regional = private$fit_tree_internal(tree, Z, Y, grid, objective_value_root_j, objective_value_root)
      self$fit_timing = list(global = t_global, regional = t_regional)
      invisible(tree)
    },

    clean = function() {
      self$tree_ref = NULL
    }
  )
)
