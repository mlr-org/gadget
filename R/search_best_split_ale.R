#' Build per-feature interval statistics for ALE effect.
#'
#' @param effect (`list()`) \cr
#'   ALE effect data per feature (from \code{calculate_ale}).
#' @param features (`character()`) \cr
#'   Feature names to include.
#'
#' @return (`list()`) \cr
#'   Statistics: K, offsets, tot_n, tot_s1, tot_s2, r_n, r_s1, r_s2, r_risks, d_l_mat, interval_idx_mat.
#' @keywords internal
build_ale_interval_stats = function(effect, features) {
  p = length(features)
  stats_list = vector("list", p)
  row_pos_list = vector("list", p)
  d_l_list = vector("list", p)
  K = integer(p)
  for (j in seq_len(p)) {
    DT = effect[[features[j]]]
    data.table::setorder(DT, row_id)
    S = unique(DT[, list(interval_index, n = int_n, s1 = int_s1, s2 = int_s2)])
    data.table::setorder(S, interval_index)
    stats_list[[j]] = S
    row_pos_list[[j]] = match(DT$interval_index, S$interval_index)
    d_l_list[[j]] = DT$d_l
    K[j] = nrow(S)
  }
  offsets = c(0L, cumsum(K))
  M = offsets[length(offsets)]
  offsets = offsets[-length(offsets)]

  tot_n = numeric(M)
  tot_s1 = numeric(M)
  tot_s2 = numeric(M)
  r_n = numeric(M)
  r_s1 = numeric(M)
  r_s2 = numeric(M)
  r_risks = numeric(p)
  pos = 1L
  for (j in seq_len(p)) {
    m = K[j]
    S = stats_list[[j]]
    rng = pos:(pos + m - 1L)
    tot_n[rng] = S$n
    tot_s1[rng] = S$s1
    tot_s2[rng] = S$s2
    r_n[rng] = S$n
    r_s1[rng] = S$s1
    r_s2[rng] = S$s2
    r_risks[j] = sum(risk_from_stats(S$n, S$s1, S$s2))
    pos = pos + m
  }
  N = length(effect[[1]]$d_l)
  d_l_mat = matrix(0.0, nrow = p, ncol = N)
  interval_idx_mat = matrix(0L, nrow = p, ncol = N)
  for (j in seq_len(p)) {
    d_l_mat[j, ] = d_l_list[[j]]
    interval_idx_mat[j, ] = row_pos_list[[j]]
  }
  list(K = K, offsets = offsets,
    tot_n = tot_n, tot_s1 = tot_s1, tot_s2 = tot_s2,
    r_n = r_n, r_s1 = r_s1, r_s2 = r_s2, r_risks = r_risks,
    d_l_mat = d_l_mat, interval_idx_mat = interval_idx_mat)
}

#' Find best ALE split across features.
#'
#' @param Z (`data.frame()` or `data.table()`) \cr
#'   Split features.
#' @param effect (`list()`) \cr
#'   ALE effect data per feature (from \code{calculate_ale}).
#' @param min_node_size (`integer(1)`) \cr
#'   Minimum observations per node.
#' @param n_quantiles (`integer(1)` or `NULL`) \cr
#'   Quantiles for numeric split candidates.
#'
#' @return (`data.frame()`) \cr
#'   Best split info with per-feature objective values.
#' @keywords internal
search_best_split_ale = function(
  Z, effect,
  min_node_size = 1L,
  n_quantiles = NULL
) {
  split_feature_names = colnames(Z)
  if (is.null(split_feature_names)) cli::cli_abort("Z (split features) must have column names.")
  st_table = build_ale_interval_stats(effect, names(effect))
  # Per split_feature, compute best split once and capture per-feature vectors
  per_feature_res = lapply(split_feature_names, function(split_feat) {
    res = search_best_split_point_ale(
      z = Z[[split_feat]],
      effect = effect,
      st_table = st_table,
      split_feat = split_feat,
      is_categorical = is.factor(Z[[split_feat]]),
      n_quantiles = n_quantiles,
      min_node_size = min_node_size
    )
    res$split_feature = split_feat
    res$is_categorical = is.factor(Z[[split_feat]])
    res
  })

  # Long format rows
  res = data.table::rbindlist(lapply(per_feature_res, function(res) {
    feature_names = names(effect)
    ovj = res$objective_value_j
    ovjl = res$left_objective_value_j
    ovjr = res$right_objective_value_j
    if (length(ovj) == 1 && is.na(ovj)) ovj = rep(NA_real_, length(feature_names))
    if (length(ovjl) == 1 && is.na(ovjl)) ovjl = rep(NA_real_, length(feature_names))
    if (length(ovjr) == 1 && is.na(ovjr)) ovjr = rep(NA_real_, length(feature_names))
    data.frame(
      split_feature = res$split_feature,
      is_categorical = res$is_categorical,
      split_point = res$split_point,
      split_objective = res$split_objective,
      feature = feature_names,
      objective_value_j = as.numeric(ovj),
      left_objective_value_j = as.numeric(ovjl),
      right_objective_value_j = as.numeric(ovjr),
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)

  min_obj = if (any(is.finite(res$split_objective))) min(res$split_objective, na.rm = TRUE) else Inf
  res$best_split = is.finite(min_obj) & (res$split_objective == min_obj)
  res[, c("split_feature", "is_categorical", "split_point",
      "split_objective", "feature", "objective_value_j",
      "left_objective_value_j", "right_objective_value_j",
      "best_split")]
}
