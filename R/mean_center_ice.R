#' Pivot tabular effect data to wide numeric columns
#'
#' Uses \code{data.table::dcast} on \code{data}: rows indexed by \code{id_cols}
#' (all columns except \code{grid_col} and \code{value_col}), columns by \code{grid_col}.
#'
#' @param data (`data.frame`) \cr
#'   Long-ish tabular effect results (ICE-style rows).
#' @param grid_col (`character(1)`) \cr
#'   Column used as column axis after pivot (often feature grid or \code{".borders"}).
#' @param value_col (`character(1)`) \cr
#'   Column holding curve values (typically \code{".value"}).
#' @param drop_cols (`character()`) \cr
#'   Columns to drop from the wide table after pivot (metadata such as \code{".type"}, \code{".id"}).
#' @param mean_center (`logical(1)`) \cr
#'   If \code{TRUE}, subtract the row mean from each row.
#'
#' @return (`data.frame`) \cr
#'   Wide numeric matrix-like frame (columns are grid cells).
#'
#' @keywords internal
pivot_effect_to_wide = function(data, grid_col, value_col = ".value", drop_cols = c(".type", ".id"),
  mean_center = FALSE) {
  if (is.factor(data[[grid_col]])) data[[grid_col]] = factor_to_numeric(data[[grid_col]])
  id_cols = setdiff(colnames(data), c(grid_col, value_col))
  formula = stats::as.formula(paste(paste(id_cols, collapse = " + "), "~", grid_col))
  wide = data.table::dcast(data.table::as.data.table(data), formula,
    value.var = value_col, fill = NA_real_)
  value_cols = setdiff(colnames(wide), drop_cols)
  out = as.data.frame(wide)[, value_cols, drop = FALSE]
  if (mean_center) out = out - rowMeans(out, na.rm = TRUE)
  out
}

#' Build wide effect matrices from iml-style effect containers
#'
#' Expects \code{effect$results} as either a single \code{data.frame} or a named list of
#' per-feature tables. Delegates pivoting to \code{pivot_effect_to_wide}.
#'
#' @param effect (R6 or `list()`) \cr
#'   Effect object with \code{results} field (e.g. from FeatureEffect).
#' @param feature_set (`character()` or `NULL`) \cr
#'   Features to include; \code{NULL} = all.
#' @param mean_center (`logical(1)`) \cr
#'   Whether to mean-center each effect matrix row-wise.
#'
#' @return (`list()`) \cr
#'   \code{Y}: named list of matrices per feature; \code{grid}: column names (grid coordinates) per feature.
#'
#' @details
#' If \code{effect$results} is a \code{data.frame}, one feature is inferred from column names.
#' If it is a named list, each element is pivoted separately.
#'
#' @examples
#' # result = mean_center_ice(effect, feature_set = c("feature1", "feature2"), mean_center = TRUE)
#' # Y = result$Y
#' # grid = result$grid
#'
#' @keywords internal
mean_center_ice = function(effect, feature_set = NULL, mean_center = TRUE) {
  checkmate::assert_true(inherits(effect, "R6") || is.list(effect), .var.name = "effect")
  checkmate::assert_character(feature_set, null.ok = TRUE, .var.name = "feature_set")
  checkmate::assert_flag(mean_center, .var.name = "mean_center")
  effect_results = effect$results

  items = if (inherits(effect_results, "data.frame")) {
    feat = colnames(effect_results)[1]
    setNames(list(list(data = effect_results, grid_col = feat, drop = c(".type", ".id"))), feat)
  } else {
    all_features = names(effect_results)
    feature_set = resolve_split_features(feature_set, all_features, "Features")
    mlr3misc::map(setNames(nm = feature_set), function(f) {
      list(data = effect_results[[f]], grid_col = ".borders",
        drop = c(".type", ".id", ".feature"))
    })
  }

  res = mlr3misc::map(items, function(it) {
    Y = pivot_effect_to_wide(it$data, it$grid_col, drop_cols = it$drop, mean_center = mean_center)
    list(Y = Y, grid = colnames(Y))
  })
  list(
    Y = mlr3misc::map(res, "Y"),
    grid = mlr3misc::map(res, "grid")
  )
}
