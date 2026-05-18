#' Build order and split candidates for ALE sweep.
#'
#' @param z (`numeric()` or `factor()`) \cr
#'   Split feature values.
#' @param is_categorical (`logical(1)`) \cr
#'   Whether \code{z} is categorical.
#' @param n_quantiles (`integer(1)` or `NULL`) \cr
#'   Quantiles for numeric features.
#'
#' @return (`list()` or `NULL`) \cr
#'   \code{ord_idx}, \code{z_sorted}, \code{n_obs}, \code{is_cand}; for categorical
#'   also \code{boundary_pos}, \code{levels_vec}. \code{NULL} if no valid split.
#' @keywords internal
build_ale_order_and_candidates = function(z, is_categorical, n_quantiles = NULL) {
  if (!is_categorical) {
    ord_idx = order(z, na.last = NA)
    z_sorted = z[ord_idx]
    n_obs = length(z_sorted)
    if (n_obs <= 1L) return(NULL)
    if (!is.null(n_quantiles) && length(unique(z_sorted)) >= n_quantiles) {
      probs = seq(0, 1, length.out = n_quantiles + 2L)[-c(1L, n_quantiles + 2L)]
      splits = unique(as.numeric(quantile(z_sorted, probs, type = 7)))
    } else {
      splits = unique(z_sorted)
    }
    t_idx = findInterval(splits, z_sorted)
    if (length(t_idx) == 0L) return(NULL)
    is_cand = rep(FALSE, n_obs - 1L)
    valid_t = t_idx[t_idx >= 1L & t_idx <= (n_obs - 1L)]
    is_cand[unique(valid_t)] = TRUE
    list(ord_idx = ord_idx, z_sorted = z_sorted, n_obs = n_obs, is_cand = is_cand)
  } else {
    z_fac = droplevels(z)
    z_non_na = which(!is.na(z_fac))
    if (length(z_non_na) <= 1L) return(NULL)
    level_id = as.integer(z_fac[z_non_na])
    ord_idx = z_non_na[order(level_id)]
    n_obs = length(ord_idx)
    counts = tabulate(level_id)
    t_idx = head(cumsum(counts), -1L)
    if (length(t_idx) == 0L) return(NULL)
    is_cand = rep(FALSE, n_obs - 1L)
    is_cand[t_idx] = TRUE
    list(
      ord_idx = ord_idx,
      z_sorted = z_fac[ord_idx],
      n_obs = n_obs,
      is_cand = is_cand,
      boundary_pos = t_idx,
      levels_vec = levels(z_fac)
    )
  }
}

#' Find best ALE split point for one feature.
#'
#' @param z (`numeric()` or `factor()`) \cr
#'   Split feature values.
#' @param effect (`list()`) \cr
#'   ALE effect data (from \code{calculate_ale}).
#' @param st_table (`list()`) \cr
#'   Precomputed interval statistics.
#' @param split_feat (`character(1)`) \cr
#'   Name of split feature.
#' @param is_categorical (`logical(1)`) \cr
#'   Whether \code{z} is categorical.
#' @param n_quantiles (`integer(1)` or `NULL`) \cr
#'   Quantiles for numeric features.
#' @param min_node_size (`integer(1)`) \cr
#'   Minimum observations per child.
#'
#' @return (`list()`) \cr
#'   \code{split_point}, \code{split_objective}, \code{objective_value_j},
#'   \code{left_objective_value_j}, \code{right_objective_value_j}.
#' @keywords internal
search_best_split_point_ale = function(
  z, effect, st_table, split_feat,
  is_categorical,
  n_quantiles = NULL,
  min_node_size = 1L
) {
  feature_names = names(effect)
  p = length(feature_names)
  split_feat_j = match(split_feat, feature_names)
  has_self_ale = !is.na(split_feat_j)

  plan = build_ale_order_and_candidates(z, is_categorical, n_quantiles)
  if (is.null(plan)) {
    return(list(
      split_point = NA_real_,
      split_objective = Inf,
      objective_value_j = rep(NA_real_, p),
      left_objective_value_j = rep(NA_real_, p),
      right_objective_value_j = rep(NA_real_, p)
    ))
  }
  ord_idx = plan$ord_idx
  z_sorted = plan$z_sorted
  n_obs = plan$n_obs
  is_cand_t = plan$is_cand

  split_feat_j_arg = if (has_self_ale) split_feat_j else 0L
  z_sorted_num = if (is.numeric(z_sorted) || is.factor(z_sorted)) as.numeric(z_sorted) else rep(0.0, n_obs)

  cpp_res = ale_sweep_cpp(
    ord_idx = ord_idx,
    d_l_mat = st_table$d_l_mat,
    interval_idx_mat = st_table$interval_idx_mat,
    offsets = st_table$offsets,
    tot_n = st_table$tot_n,
    tot_s1 = st_table$tot_s1,
    tot_s2 = st_table$tot_s2,
    r_risks = st_table$r_risks,
    is_cand = is_cand_t,
    min_node_size = min_node_size,
    split_feat_j = split_feat_j_arg,
    use_stabilizer = FALSE,
    z_sorted = z_sorted_num,
    n_obs = n_obs
  )
  best_t = cpp_res$best_t
  best_risks_sum = cpp_res$best_risks_sum
  best_left_risks = cpp_res$best_left_risks
  best_right_risks = cpp_res$best_right_risks

  if (is.na(best_t) || best_t < 0L) { # in cpp: best_t = -1 means no valid split point found
    return(list(
      split_point = NA_real_,
      split_objective = Inf,
      objective_value_j = rep(NA_real_, p),
      left_objective_value_j = rep(NA_real_, p),
      right_objective_value_j = rep(NA_real_, p)
    ))
  }
  if (!is_categorical) {
    left_value = max(z_sorted[1:best_t])
    right_value = min(z_sorted[-(1:best_t)])
    best_split_point = (left_value + right_value) / 2
  } else {
    k = which(plan$boundary_pos == best_t)[1]
    best_split_point = if (!is.na(k)) plan$levels_vec[k] else NA
  }
  list(
    split_point = best_split_point,
    split_objective = best_risks_sum,
    objective_value_j = best_left_risks + best_right_risks,
    left_objective_value_j = best_left_risks,
    right_objective_value_j = best_right_risks
  )
}
