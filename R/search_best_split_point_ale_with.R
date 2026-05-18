#' Slow reference ALE split point search with boundary stabilizer.
#'
#' Same params and return as \code{search_best_split_point_ale}.
#'
#' @param z (`numeric()` or `factor()`) \cr
#'   Split feature values.
#' @param effect (`list()`) \cr
#'   ALE effect data.
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
#'   Same as \code{search_best_split_point_ale}.
#' @keywords internal
search_best_split_point_ale_with = function(
  z, effect, st_table, split_feat,
  is_categorical,
  n_quantiles = NULL,
  min_node_size = 1L
) {
  feature_names = names(effect)
  p = length(feature_names)
  split_feat_j = match(split_feat, feature_names)
  has_self_ale = !is.na(split_feat_j)

  # Helper: Move one row across all features (updates flattened vectors and per-feature risks)
  move_row_all = function(row_id) {
    d = st_table$d_l_mat[, row_id] # length p, d_l of each feature
    interval_idx = st_table$interval_idx_mat[, row_id] # length p, interval number of each feature
    m = offsets + interval_idx # length p, subscript of each feature in flattened vector

    r_n_old = st_table$r_n[m]
    r_s1_old = st_table$r_s1[m]
    r_s2_old = st_table$r_s2[m]
    r_risk_old = risk_from_stats(r_n_old, r_s1_old, r_s2_old)

    r_n_new = r_n_old - 1.0
    r_s1_new = r_s1_old - d
    r_s2_new = r_s2_old - d * d
    r_risk_new = risk_from_stats(r_n_new, r_s1_new, r_s2_new)

    st_table$r_n[m] <<- r_n_new
    st_table$r_s1[m] <<- r_s1_new
    st_table$r_s2[m] <<- r_s2_new
    r_risks <<- r_risks - r_risk_old + r_risk_new

    l_n_old = st_table$tot_n[m] - r_n_old
    l_s1_old = st_table$tot_s1[m] - r_s1_old
    l_s2_old = st_table$tot_s2[m] - r_s2_old
    l_risk_old = risk_from_stats(l_n_old, l_s1_old, l_s2_old)

    l_n_new = st_table$tot_n[m] - r_n_new
    l_s1_new = st_table$tot_s1[m] - r_s1_new
    l_s2_new = st_table$tot_s2[m] - r_s2_new
    l_risk_new = risk_from_stats(l_n_new, l_s1_new, l_s2_new)

    l_risks <<- l_risks - l_risk_old + l_risk_new
    sum(-l_risk_old - r_risk_old + l_risk_new + r_risk_new)
  }

  #### Old Helper: Boundary stabilizer using cumulative sums (numeric feature only) ####
  # adjust_side_for_feature = function(side, t, w) {
  #   # w = window size (# of points near boundary)
  #   if (w <= 0L) {
  #     return(if (side == "L") l_risks[split_feat_j] else r_risks[split_feat_j])
  #   }
  #   if (side == "L") {
  #     # Left node: ord[1..t]
  #     tot_n = t
  #     if (tot_n <= 1L) {
  #       return(l_risks[split_feat_j])
  #     }
  #     w = min(w, tot_n)
  #     # near: last w samples on the left
  #     s1_near = S1[t] - if (t - w >= 1L) S1[t - w] else 0.0
  #     s2_near = S2[t] - if (t - w >= 1L) S2[t - w] else 0.0
  #
  #     risk_j = l_risks[split_feat_j]
  #     s1_tot = S1[t]
  #     s2_tot = S2[t]
  #   } else {
  #     # Right node: ord[t+1..n_obs]
  #     tot_n = n_obs - t
  #     if (tot_n <= 1L) {
  #       return(r_risks[split_feat_j])
  #     }
  #     w = min(w, tot_n)
  #     # near: first w samples on the right
  #     s1_near = S1[t + w] - S1[t]
  #     s2_near = S2[t + w] - S2[t]
  #
  #     risk_j = r_risks[split_feat_j]
  #     # totals over the right node
  #     s1_tot = S1_tot - S1[t]
  #     s2_tot = S2_tot - S2[t]
  #   }
  #
  #   # far region = node minus near window
  #   n_far = tot_n - w
  #   if (n_far <= 1L) {
  #     return(risk_j)
  #   }
  #   s1_far = s1_tot - s1_near
  #   s2_far = s2_tot - s2_near
  #   var_near = (s2_near - s1_near * s1_near / w) / max(w - 1L, 1L)
  #   var_far = (s2_far - s1_far * s1_far / n_far) / max(n_far - 1L, 1L)
  #   if (is.na(var_far) || is.na(var_near) || var_far <= 0 || var_near < 4.0 * var_far) {
  #     return(risk_j)
  #   }
  #   risk_j - risk_from_stats(w, s1_near, s2_near) + w * var_far
  # }
  #### End old adjust_side_for_feature function ####

  #### New Helper: Boundary stabilizer using cumulative sums (numeric feature only) ####
  adjust_side_for_feature = function(side, t) {
    # 1. Determine range of current side
    is_left = (side == "L")
    if (is_left) {
      idx_start = 1L
      idx_end = t
      original_risk = l_risks[split_feat_j]
    } else {
      idx_start = t + 1L
      idx_end = n_obs
      original_risk = r_risks[split_feat_j]
    }
    n_side = idx_end - idx_start + 1L
    # 2. if sample size <= 20, do not smooth
    if (n_side <= 20L) {
      return(original_risk)
    }
    # 3. window size w: max(10% of side samples, 10)
    w = max(round(0.1 * n_side), 10L)
    # 4. Get interval indices and d_l for current side
    node_int = interval_idx_sorted[idx_start:idx_end]
    node_d_l = d_l_j_sorted[idx_start:idx_end]
    # 5. Get unique set of interval indices within the window
    if (is_left) {
      window_indices = (n_side - w + 1L):n_side # Left window: last w elements
    } else {
      window_indices = 1L:w # Right window: first w elements
    }
    target_intervals = unique(node_int[window_indices])
    # 6. Define Near / Far samples
    # As long as it belongs to target_intervals, treat as Near (regardless of whether inside window w)
    is_near = node_int %in% target_intervals
    # 7. Sample size check
    n_near = sum(is_near)
    n_far = n_side - n_near
    if (n_near < 2 || n_far < 2) {
      return(original_risk)
    }
    # 8. Get d_l and calculate SD
    vals_near = node_d_l[is_near]
    vals_far = node_d_l[!is_near]
    sd_near = sd(vals_near)
    sd_far = sd(vals_far)
    # 9. Condition and replacement
    if (!is.na(sd_near) && !is.na(sd_far) && sd_near > 2.0 * sd_far) {
      node_d_l[is_near] = rnorm(n = n_near, mean = mean(vals_far), sd = sd_far)
      s1_vec = rowsum(node_d_l, group = node_int, reorder = FALSE)
      s2_vec = rowsum(node_d_l^2, group = node_int, reorder = FALSE)
      n_vec = rowsum(rep(1L, length(node_d_l)), group = node_int, reorder = FALSE)
      interval_sse = s2_vec - (s1_vec^2) / n_vec
      return(sum(interval_sse))
    }
    original_risk
  }
  #### End new adjust_side_for_feature function ####

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

  offsets = st_table$offsets
  r_risks = st_table$r_risks
  l_risks = numeric(p)

  #### Code with adjust_side_for_feature ####
  if (!is_categorical && has_self_ale) {
    d_l_j_sorted = st_table$d_l_mat[split_feat_j, ord_idx] # N x 1
    interval_idx_sorted = st_table$interval_idx_mat[split_feat_j, ord_idx]
    #### Code for old adjust_side_for_feature ####
    # S1 = cumsum(d_l_j_sorted)
    # S2 = cumsum(d_l_j_sorted * d_l_j_sorted)
    # S1_tot = S1[n_obs]
    # S2_tot = S2[n_obs]
    #### ####
  }
  #### ####

  risks_sum = sum(r_risks)
  best_risks_sum = Inf
  best_t = NA_integer_
  best_left_risks = NULL
  best_right_risks = NULL

  # Main sweep
  for (t in 1:(n_obs - 1L)) {
    delta = move_row_all(ord_idx[t])
    risks_sum = risks_sum + delta
    if (!is_cand_t[t] || t < min_node_size || (n_obs - t) < min_node_size) next
    curr_l_risks = l_risks
    curr_r_risks = r_risks
    #### Code with adjust_side_for_feature ####
    # Handle single unique value case for the split feature
    l_const = z_sorted[1] == z_sorted[t]
    r_const = z_sorted[t + 1L] == z_sorted[n_obs]
    if (has_self_ale) {
      if (l_const) curr_l_risks[split_feat_j] = 0.0
      if (r_const) curr_r_risks[split_feat_j] = 0.0
    }
    #### ####
    #### Code without adjust_side_for_feature ####
    # When splitting on a feature that has its own ALE, exclude its own risk from the objective
    # This avoids re-fitting the feature's main effect and focuses on interaction effects
    # if (has_self_ale) {
    #   curr_l_risks[split_feat_j] = 0.0
    #   curr_r_risks[split_feat_j] = 0.0
    # }
    #### ####
    total = sum(curr_l_risks) + sum(curr_r_risks)
    #### Code with adjust_side_for_feature ####
    # Boundary Stabilizer
    use_stabilizer = !is_categorical && has_self_ale && !l_const && !r_const && t > 20 && (n_obs - t) > 20
    if (use_stabilizer) {
      risk_t_j = l_risks[split_feat_j] + r_risks[split_feat_j]
      #### Code for old adjust_side_for_feature ####
      # w_l = max(round(0.1 * t), 10L)
      # w_l = min(w_l, t)
      # w_r = max(round(0.1 * (n_obs - t)), 10L)
      # w_r = min(w_r, n_obs - t)
      # adj_l_risk_t_j = adjust_side_for_feature("L", t, w_l)
      # adj_r_risk_t_j = adjust_side_for_feature("R", t, w_r)
      #### ####
      adj_l_risk_t_j = adjust_side_for_feature("L", t)
      adj_r_risk_t_j = adjust_side_for_feature("R", t)
      adj_risk_t_j = adj_l_risk_t_j + adj_r_risk_t_j
      total = risks_sum - risk_t_j + adj_risk_t_j
    }
    #### ####
    if (!is.na(total) && !is.infinite(total) && total < best_risks_sum) {
      best_risks_sum = total
      best_t = t
      best_left_risks = curr_l_risks
      best_right_risks = curr_r_risks
      #### Code with adjust_side_for_feature ####
      if (use_stabilizer) {
        best_left_risks[split_feat_j] = adj_l_risk_t_j
        best_right_risks[split_feat_j] = adj_r_risk_t_j
      }
      #### ####
    }
  }
  if (is.na(best_t)) {
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

#' Slow reference ALE split point search via C++ sweep.
#'
#' Same params and return as \code{search_best_split_point_ale}.
#'
#' @param z (`numeric()` or `factor()`) \cr
#'   Split feature values.
#' @param effect (`list()`) \cr
#'   ALE effect data.
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
#'   Same as \code{search_best_split_point_ale}.
#' @keywords internal
search_best_split_point_ale_with_cpp = function(
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
  use_stab = !is_categorical && has_self_ale
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
    use_stabilizer = use_stab,
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
