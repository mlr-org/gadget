#' Prepare PD Data for Tree Splitting
#'
#' Given effect, data, and optional feature/split sets: resolves features;
#' converts character to factor; builds Z (split columns); calls \code{mean_center_ice} for Y and grid.
#'
#' @param effect (R6 or `list()`) \cr
#'   Effect object (e.g. FeatureEffect).
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Training data.
#' @param target_feature_name (`character(1)` or `NULL`) \cr
#'   Target variable name; \code{NULL} = all columns are features.
#' @param feature_set (`character()` or `NULL`) \cr
#'   Features in effect; \code{NULL} = all non-target columns.
#' @param split_feature (`character()` or `NULL`) \cr
#'   Features for splitting; \code{NULL} = all.
#'
#' @return (`list()`) \cr
#'   \code{Z}: split-feature data.table; \code{Y}: mean-centered effects; \code{grid}: grid list.
#'
#' @keywords internal
prepare_split_data_pd = function(effect, data, target_feature_name = NULL, feature_set = NULL,
  split_feature = NULL) {
  common = prepare_split_data_common(data, target_feature_name, feature_set, split_feature)
  wide_mean_center = mean_center_ice(effect = effect, feature_set = common$split_feature)
  list(Z = common$Z, Y = wide_mean_center$Y, grid = wide_mean_center$grid)
}
