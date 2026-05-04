#' Prepare ALE Data for Tree Splitting
#'
#' Given model, data, target_feature_name, n_intervals, and optional feature/split sets:
#' validates features; converts character to factor and orders levels via
#' \code{order_categorical_levels}; builds Z (data.table of split columns);
#' calls \code{calculate_ale} for Y. Returns list \code{Z}, \code{Y}.
#'
#' @param model (`any`) \cr
#'   Fitted model with predict interface.
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data (features and target).
#' @param target_feature_name (`character(1)`) \cr
#'   Name of the target variable.
#' @param n_intervals (`integer(1)`) \cr
#'   Number of intervals for numeric ALE.
#' @param feature_set (`character()` or `NULL`) \cr
#'   Features to compute ALE for; \code{NULL} = all.
#' @param split_feature (`character()` or `NULL`) \cr
#'   Features for splitting; \code{NULL} = all.
#' @param predict_fun (`function()` or `NULL`) \cr
#'   \code{function(model, data)} returning predictions; \code{NULL} = default.
#' @param order_method (`character(1)`) \cr
#'   Categorical level order: \code{"mds"}, \code{"pca"}, \code{"random"}, or \code{"raw"}.
#' @param ale_engine (`character(1)`) \cr
#'   ALE engine: \code{"cpp"} or \code{"r"}; default \code{c("cpp", "r")} resolves to
#'   \code{"cpp"} via \code{match.arg}.
#' @return (`list()`) \cr
#'   \code{Z}: data.table of split features; \code{Y}: list of ALE effect data per feature.
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Resolve \code{feature_set} and \code{split_feature} against feature columns (see
#'     \code{GadgetTree}/\code{AleStrategy} for target-column validation).
#'   \item For \code{union(feature_set, split_feature)}, convert character
#'     columns to factor and order levels via
#'     \code{order_categorical_levels} (using \code{droplevels} internally).
#'   \item Build \code{Z} as \code{data[split_feature]} (data.table).
#'   \item Call \code{calculate_ale(model, data, feature_set, ...)} to get \code{Y}.
#' }
#' Stops with an error if any requested feature is missing from \code{data}.
#'
#' @examples
#' \dontrun{
#' result = prepare_split_data_ale(model, data, "y", n_intervals = 10)
#' result = prepare_split_data_ale(model, data, "y", n_intervals = 10,
#'   feature_set = c("x1", "x2"), split_feature = c("x1", "x2", "x3"))
#' }
#'
#' @keywords internal
prepare_split_data_ale = function(
  model,
  data,
  target_feature_name,
  n_intervals,
  feature_set = NULL,
  split_feature = NULL,
  predict_fun = NULL,
  order_method = "raw",
  ale_engine = c("cpp", "r")
) {
  ale_engine = match.arg(ale_engine)
  common = prepare_split_data_common(data, target_feature_name, feature_set, split_feature)
  data = common$data
  feature_set = common$feature_set
  split_feature = common$split_feature
  for (col in union(feature_set, split_feature)) {
    if (is.factor(data[[col]])) {
      data[[col]] = order_categorical_levels(droplevels(data[[col]]), data, col, target_feature_name, order_method)
    }
  }
  Z = data.table::setDT(take_cols(data, split_feature))
  effect = if (ale_engine == "cpp") {
    calculate_ale_fast(
      model = model,
      data = data,
      target_feature_name = target_feature_name,
      feature_set = feature_set,
      n_intervals = n_intervals,
      predict_fun = predict_fun
    )
  } else {
    calculate_ale(
      model = model,
      data = data,
      target_feature_name = target_feature_name,
      feature_set = feature_set,
      n_intervals = n_intervals,
      predict_fun = predict_fun
    )
  }

  list(Z = Z, Y = effect)
}
