#' AleStrategy: Generalized Additive Decomposition Based on ALE Effects
#'
#' @description
#' ALE-based effect strategy (inherits from \link{EffectStrategy}). Given model and data,
#' preprocesses to Z/Y via \code{prepare_split_data_ale};
#' transforms ALE effects per node; computes ALE-derivative heterogeneity;
#' finds best split via \code{search_best_split_ale}; fits tree and plots ALE curves.
#'
#' @usage NULL
#' @format [R6::R6Class] object inheriting from [EffectStrategy].
#'
#' @section Construction:
#' ```
#' s = AleStrategy$new()
#' ```
#'
#' @field model (`any`) \cr
#'   Fitted model (persistent after \code{$fit()}).
#' @field data (`data.frame()` or `data.table()`) \cr
#'   Data (persistent after \code{$fit()}).
#' @field target_feature_name (`character(1)`) \cr
#'   Target variable name.
#' @field n_intervals (`integer(1)`) \cr
#'   Intervals for numeric ALE.
#' @field predict_fun (`function()`) \cr
#'   \code{function(model, data)} returning predictions.
#' @field order_method (`character(1)`) \cr
#'   Categorical order: \code{"mds"}, \code{"pca"}, \code{"random"}, \code{"raw"}.
#' @field ale_engine (`character(1)`) \cr
#'   ALE backend after \code{$fit()}: \code{"cpp"} or \code{"r"}.
#' @field effect (`list()` or `NULL`) \cr
#'   Cached ALE effect used when \code{$plot()} omits \code{effect}.
#'
#' @details
#' Intended for use through \code{GadgetTree$new(strategy = AleStrategy$new())} and
#' \code{tree$fit(...)}. Can be instantiated directly for custom pipelines.
#'
#' @examples
#' \dontrun{
#' ale_strat = AleStrategy$new()
#' tree = GadgetTree$new(strategy = ale_strat, n_split = 2)
#' tree$fit(model = model, data = data, target_feature_name = "y")
#' }
#'
#' @include EffectStrategy.R
#' @export
AleStrategy = R6::R6Class(
  "AleStrategy",
  inherit = EffectStrategy,
  public = list(
    # Persistent context fields (declared to allow assignment outside initialize)
    model = NULL,
    data = NULL,
    target_feature_name = NULL,
    n_intervals = NULL,
    predict_fun = NULL,
    order_method = "raw",
    ale_engine = "cpp",
    effect = NULL,

    #' @description
    #' Create an AleStrategy instance (calls \code{super$initialize("ale")}).
    initialize = function() {
      super$initialize("ale")
    },

    #' @description
    #' Preprocess to Z and Y via \code{prepare_split_data_ale}.
    #' @param model (`any`) \cr
    #'   Fitted model.
    #' @param effect (`list()` or `NULL`) \cr
    #'   Reserved for future extension. Currently unsupported.
    #' @param data (`data.frame()` or `data.table()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target variable name.
    #' @param n_intervals (`integer(1)`) \cr
    #'   Intervals for numeric ALE.
    #' @param feature_set (`character()` or `NULL`) \cr
    #'   Features for ALE; \code{NULL} = all.
    #' @param split_feature (`character()` or `NULL`) \cr
    #'   Features for splitting; \code{NULL} = all.
    #' @param predict_fun (`function()` or `NULL`) \cr
    #'   Prediction function.
    #' @param order_method (`character(1)`) \cr
    #'   Categorical order: \code{"mds"}, \code{"pca"}, \code{"random"}, or \code{"raw"}.
    #' @param ale_engine (`character(1)`) \cr
    #'   ALE engine: \code{"cpp"} or \code{"r"}; default \code{c("cpp", "r")} resolves to
    #'   \code{"cpp"} via \code{match.arg}.
    #' @return (`list()`) \cr
    #'   \code{Z}: split features; \code{Y}: ALE effect data.tables.
    preprocess = function(model, effect = NULL, data, target_feature_name, n_intervals,
      feature_set = NULL, split_feature = NULL, predict_fun = NULL,
      order_method = "raw", ale_engine = c("cpp", "r")) {
      ale_engine = match.arg(ale_engine)
      checkmate::assert_data_frame(data, .var.name = "data")
      checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
      checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      checkmate::assert_integerish(n_intervals, len = 1, lower = 1, .var.name = "n_intervals")
      checkmate::assert_character(feature_set, null.ok = TRUE, .var.name = "feature_set")
      checkmate::assert_character(split_feature, null.ok = TRUE, .var.name = "split_feature")
      checkmate::assert_function(predict_fun, null.ok = TRUE, .var.name = "predict_fun")
      checkmate::assert_choice(order_method, c("mds", "pca", "random", "raw"), .var.name = "order_method")
      if (!is.null(effect)) {
        cli::cli_abort(
          "Direct ALE {.arg effect} input is not enabled yet; please pass {.arg model} and use gadget ALE."
        )
      }
      if (is.null(model)) {
        cli::cli_abort("{.arg model} is required for AleStrategy$preprocess.")
      }
      prepare_split_data_ale(
        model = model, data = data, target_feature_name = target_feature_name,
        n_intervals = n_intervals, feature_set = feature_set, split_feature = split_feature,
        predict_fun = predict_fun, order_method = order_method, ale_engine = ale_engine
      )
    },

    #' @description
    #' Subset ALE by node indices; handle single-interval and categorical.
    #' @param Y (`list()`) \cr
    #'   ALE effect list from \code{calculate_ale}.
    #' @param idx (`integer()`) \cr
    #'   Row indices in the node.
    #' @param grid (`list()` or `NULL`) \cr
    #'   Ignored for ALE; required by interface.
    #' @param is_child (`logical(1)`) \cr
    #'   Whether the current node is a child node.
    #' @return (`list()`) \cr
    #'   Transformed ALE data.tables.
    node_transform = function(Y, idx, grid = NULL, is_child = FALSE) {
      assert_ale_effect_list(Y)
      checkmate::assert_integerish(idx, min.len = 1, .var.name = "idx")
      checkmate::assert_list(grid, null.ok = TRUE, .var.name = "grid")
      checkmate::assert_flag(is_child, .var.name = "is_child")

      node_transform_ale(
        Y = Y,
        idx = idx,
        is_child = is_child
      )
    },

    #' @description
    #' Compute ALE heterogeneity via \code{calculate_ale_heterogeneity_cpp}.
    #' @param Y (`list()`) \cr
    #'   ALE effect list from \code{calculate_ale}.
    #' @return (`numeric()`) \cr
    #'   Heterogeneity per feature.
    #' @seealso \code{\link{calculate_ale_heterogeneity_cpp}}
    heterogeneity = function(Y) {
      assert_ale_effect_list(Y)
      unlist(calculate_ale_heterogeneity_cpp(Y))
    },

    #' @description
    #' Compute left/right child objective values from split result.
    #' For ALE, extracts from split_info (computed during sweep).
    #' @param Z,Y,split_info,idx_left,idx_right,grid_left,grid_right \cr
    #'   Node/split context.
    #' @return (`list()`) \cr
    #'   \code{left_objective_value_j}, \code{right_objective_value_j},
    #'   \code{left_objective_value}, \code{right_objective_value}.
    get_child_objectives = function(Z, Y, split_info, idx_left, idx_right, grid_left, grid_right) {
      checkmate::assert_list(split_info, .var.name = "split_info")
      raw = split_info$raw_result
      if (is.null(raw) || is.null(raw$left_objective_value_j)) {
        cli::cli_abort("ALE split_info must contain raw_result with left/right objective values.")
      }
      rows = which(raw$best_split & raw$split_feature == split_info$split_feature)
      if (length(rows) == 0L) {
        cli::cli_abort("ALE split_info does not contain objective rows for the selected split.")
      }
      left_objective_value_j = raw$left_objective_value_j[rows]
      right_objective_value_j = raw$right_objective_value_j[rows]
      names(left_objective_value_j) = raw$feature[rows]
      names(right_objective_value_j) = raw$feature[rows]
      list(
        left_objective_value_j = left_objective_value_j,
        right_objective_value_j = right_objective_value_j,
        left_objective_value = sum(left_objective_value_j, na.rm = TRUE),
        right_objective_value = sum(right_objective_value_j, na.rm = TRUE)
      )
    },

    #' @description
    #' Find best split via \code{search_best_split_ale}.
    #' @param Z (`data.frame()` or `data.table()`) \cr
    #'   Split features.
    #' @param Y (`list()`) \cr
    #'   ALE effect from \code{calculate_ale}.
    #' @param min_node_size (`integer(1)`) \cr
    #'   Minimum node size.
    #' @param n_quantiles (`integer(1)` or `NULL`) \cr
    #'   Quantile candidates for numeric.
    #' @return (`list()` or `data.frame()`) \cr
    #'   Best split info: \code{split_feature}, \code{split_point}, etc.
    find_best_split = function(Z, Y, min_node_size, n_quantiles) {
      checkmate::assert_true(data.table::is.data.table(Z) || is.data.frame(Z), .var.name = "Z")
      assert_ale_effect_list(Y)
      checkmate::assert_integerish(min_node_size, len = 1, lower = 1, any.missing = FALSE, .var.name = "min_node_size")
      checkmate::assert_integerish(n_quantiles, len = 1, lower = 1, null.ok = TRUE, .var.name = "n_quantiles")

      search_best_split_ale(
        Z = Z,
        effect = Y,
        min_node_size = min_node_size,
        n_quantiles = n_quantiles
      )
    },

    #' @description
    #' Plot ALE curves via \code{plot_tree_ale}.
    #' @param tree (`list()`) \cr
    #'   Depth-based list of Node objects.
    #' @param effect (`list()` or `NULL`) \cr
    #'   ALE effect; \code{NULL} = use cached \code{effect}.
    #' @param data (`data.frame()` or `data.table()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target variable name.
    #' @param depth (`integer()` or `NULL`) \cr
    #'   Depths to plot.
    #' @param node_id (`integer()` or `NULL`) \cr
    #'   Node IDs to plot.
    #' @param features (`character()` or `NULL`) \cr
    #'   Features to include.
    #' @param show_plot,show_point,mean_center (`logical(1)`) \cr
    #'   Plot options.
    #' @param ... Passed to \code{plot_tree_ale}.
    #' @return (`list()`) \cr
    #'   Nested list (depth -> node -> patchwork).
    plot = function(tree, effect = NULL, data, target_feature_name,
      depth = NULL, node_id = NULL, features = NULL,
      show_plot = TRUE, show_point = TRUE, mean_center = TRUE, ...) {
      checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      if (is.null(effect)) {
        effect = self$effect
        if (is.null(effect)) {
          cli::cli_abort("No cached ALE effect found. Please pass 'effect' or run fit() first.")
        }
      }
      plot_tree_ale(tree = tree, effect = effect, data = data,
        target_feature_name = target_feature_name, depth = depth,
        node_id = node_id, features = features, show_plot = show_plot,
        show_point = show_point, mean_center = mean_center, ...)
    },

    #' @description
    #' Fit tree: preprocess, create root, split recursively.
    #' @param tree (`GadgetTree`) \cr
    #'   Tree instance.
    #' @param model (`any`) \cr
    #'   Fitted model.
    #' @param effect (`list()` or `NULL`) \cr
    #'   Reserved for future extension. Currently unsupported.
    #' @param data (`data.frame()` or `data.table()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target name.
    #' @param n_intervals (`integer(1)`) \cr
    #'   Intervals for numeric ALE.
    #' @param feature_set,split_feature (`character()` or `NULL`) \cr
    #'   Feature subsets.
    #' @param predict_fun (`function()` or `NULL`) \cr
    #'   Prediction function.
    #' @param order_method (`character(1)`) \cr
    #'   Categorical order.
    #' @param ale_engine (`character(1)`) \cr
    #'   ALE engine: \code{"cpp"} or \code{"r"}; default \code{c("cpp", "r")} resolves to
    #'   \code{"cpp"} via \code{match.arg}.
    #' @param ... Ignored.
    #' @return (`GadgetTree`) \cr
    #'   The tree, invisibly.
    fit = function(tree, model, effect = NULL, data, target_feature_name,
      n_intervals = 10, feature_set = NULL, split_feature = NULL,
      predict_fun = NULL, order_method = "raw", ale_engine = c("cpp", "r"), ...) {
      checkmate::assert_r6(tree, classes = "GadgetTree", .var.name = "tree")
      checkmate::assert_data_frame(data, .var.name = "data")
      checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
      checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      checkmate::assert_integerish(n_intervals, len = 1, lower = 1, .var.name = "n_intervals")
      checkmate::assert_list(effect, null.ok = TRUE, .var.name = "effect")
      checkmate::assert_character(feature_set, null.ok = TRUE, .var.name = "feature_set")
      checkmate::assert_character(split_feature, null.ok = TRUE, .var.name = "split_feature")
      checkmate::assert_function(predict_fun, null.ok = TRUE, .var.name = "predict_fun")
      checkmate::assert_choice(order_method, c("mds", "pca", "random", "raw"), .var.name = "order_method")
      ale_engine = match.arg(ale_engine)
      dots = list(...)
      if ("with_stab" %in% names(dots)) {
        cli::cli_abort(
          "{.arg with_stab} has been retired. ALE split search now always uses the fast bias-corrected objective."
        )
      }
      if (is.null(model)) {
        cli::cli_abort("AleStrategy requires {.arg model} to be passed.")
      }
      if (!is.null(effect)) {
        cli::cli_abort(
          "Direct ALE {.arg effect} input is not enabled yet; please pass {.arg model} and use gadget ALE."
        )
      }

      if (is.null(predict_fun)) {
        predict_fun = default_predict_fun
      }
      t_global = system.time({
        self$model = model
        self$data = data
        self$target_feature_name = target_feature_name
        self$n_intervals = n_intervals
        self$predict_fun = predict_fun
        self$order_method = order_method
        self$ale_engine = ale_engine
        self$tree_ref = tree
        prepared_data = self$preprocess(model = model, effect = effect, data = data,
          target_feature_name = target_feature_name, n_intervals = n_intervals,
          feature_set = feature_set, split_feature = split_feature,
          predict_fun = predict_fun, order_method = order_method, ale_engine = ale_engine)
        Z = prepared_data$Z
        Y = prepared_data$Y
        self$effect = Y
        grid = vector("list", length(names(Z)))
        names(grid) = names(Z)
        objective_value_root_j = self$heterogeneity(Y)
        objective_value_root = sum(objective_value_root_j, na.rm = TRUE)
      })[["elapsed"]]

      t_regional = private$fit_tree_internal(tree, Z, Y, grid, objective_value_root_j, objective_value_root)
      self$fit_timing = list(global = t_global, regional = t_regional)
      invisible(tree)
    },
    #' @description
    #' Sets \code{data} and \code{model} to \code{NULL} to free memory after fitting.
    #' \code{effect} is intentionally retained because \code{plot()} requires it post-fit.
    clean = function() {
      self$data = NULL
      self$model = NULL
      self$tree_ref = NULL
    }
  )
)
