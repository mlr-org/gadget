/**
 * @file ale_sweep.cpp
 * @brief Fast ALE main sweep for tree splitting - C++ implementation
 *
 * Implements the core loop of ALE-based split search: rows are ordered by the
 * split feature z; we sweep through split positions t=1..n_obs-1, moving one
 * row at a time from the right child to the left. At each position, we compute
 * the total heterogeneity (sum of SSE over intervals and features) and track
 * the best split. Uses sufficient statistics (n, s1, s2 per interval) for O(1)
 * incremental updates. Replaces the R-level for-loop in search_best_split_point_ale
 * and search_best_split_point_ale_with.
 */

#include <Rcpp.h>
#include <set>
#include <map>
#include <cmath>

// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;

/* Minimum observations per child side required to run the boundary stabilizer.
 * With n = 20 the window is ~2 rows (10%), leaving ~18 "far" observations for
 * variance estimation -- the practical lower bound for a meaningful stabilizer. */
static const int kMinSideForStab = 20;

// -----------------------------------------------------------------------------
// risk_from_stats
// Purpose:
//   Compute SSE from sufficient statistics: Var = E[X^2] - E[X]^2, so
//   SSE = sum(x^2) - sum(x)^2/n = s2 - s1^2/n. Used for per-interval risk.
// Notes:
//   Returns 0 when n <= 1 (no variance for single point).
// -----------------------------------------------------------------------------
inline double risk_from_stats(double n, double s1, double s2) {
  if (n <= 1.0) return 0.0;
  return s2 - (s1 * s1) / n;
}

// -----------------------------------------------------------------------------
// adjust_side_cpp
// Purpose:
//   Boundary stabilizer for numeric split feature. Near the split boundary,
//   few samples in edge intervals can inflate heterogeneity. This function
//   identifies a "near" window (last w rows for left child, first w for right),
//   compares d_l variance of near vs "far" intervals. If near is much larger,
//   replaces near d_l with draws from N(mean_far, sd_far) and recomputes SSE.
//   Replicates adjust_side_for_feature logic from R (aleStrategy).
// Inputs:
//   node_d_l: d_l values for the split feature, ordered as in sweep
//   node_int: interval index per row
//   original_risk: raw SSE for this side; returned if stabilizer does not apply
//   is_left, t, n_obs: describe side (left=t rows, right=n_obs-t rows)
// Notes:
//   NA/NaN in node_d_l are skipped. Uses R::rnorm for stochastic replacement;
//   reproducibility depends on R's RNG state.
// -----------------------------------------------------------------------------
double adjust_side_cpp(
    const NumericVector& node_d_l,
    const IntegerVector& node_int,
    double original_risk,
    bool is_left,
    int t,
    int n_obs
) {
  int idx_start = is_left ? 0 : t;
  int idx_end   = is_left ? t - 1 : n_obs - 1;
  int n_side = idx_end - idx_start + 1;
  if (n_side <= kMinSideForStab) return original_risk;

  /* Window size: 10% of side or min 10. */
  int w = std::max(static_cast<int>(0.1 * n_side + 0.5), 10);
  if (w > n_side) w = n_side;

  /* Near = boundary window: last w rows for left child, first w for right. */
  IntegerVector window_indices(w);
  if (is_left) {
    for (int i = 0; i < w; ++i) window_indices[i] = n_side - w + i;
  } else {
    for (int i = 0; i < w; ++i) window_indices[i] = i;
  }

  /* Intervals present in the boundary window = "near". */
  std::set<int> target_intervals;
  for (int i = 0; i < w; ++i) {
    target_intervals.insert(node_int[window_indices[i]]);
  }

  /* Near/Far split; skip NA/NaN in node_d_l. */
  std::vector<double> vals_near, vals_far;
  for (int i = 0; i < n_side; ++i) {
    double v = node_d_l[i];
    if (NumericVector::is_na(v)) continue;
    if (target_intervals.count(node_int[i]))
      vals_near.push_back(v);
    else
      vals_far.push_back(v);
  }
  int n_near = vals_near.size();
  int n_far = vals_far.size();
  if (n_near < 2 || n_far < 2) return original_risk;

  /* Mean and SD of far intervals (used for replacement distribution). */
  double sum_far = 0.0, sum2_far = 0.0;
  for (double v : vals_far) {
    sum_far += v;
    sum2_far += v * v;
  }
  double mean_far = sum_far / n_far;
  double var_far = (sum2_far - sum_far * sum_far / n_far) / std::max(n_far - 1, 1);
  double sd_far = (var_far > 0) ? std::sqrt(var_far) : 0.0;

  /* SD of near intervals (compared to far for stabilizer condition). */
  double sum_near = 0.0, sum2_near = 0.0;
  for (double v : vals_near) {
    sum_near += v;
    sum2_near += v * v;
  }
  double var_near = (sum2_near - sum_near * sum_near / n_near) / std::max(n_near - 1, 1);
  double sd_near = (var_near > 0) ? std::sqrt(var_near) : 0.0;

  /* Only stabilize when near has much higher variance than far. */
  if (R_IsNA(sd_near) || R_IsNA(sd_far) || sd_far <= 0 || sd_near <= 2.0 * sd_far)
    return original_risk;

  /* Replace near d_l with draws from N(mean_far, sd_far); recompute SSE per interval. */
  std::map<int, double> s1_map, s2_map, n_map;
  for (int i = 0; i < n_side; ++i) {
    double v = node_d_l[i];
    if (NumericVector::is_na(v)) continue;
    int k = node_int[i];
    if (target_intervals.count(k))
      v = R::rnorm(mean_far, sd_far);
    s1_map[k] += v;
    s2_map[k] += v * v;
    n_map[k] += 1.0;
  }
  double sse = 0.0;
  for (auto& it : n_map) {
    int k = it.first;
    double nk = n_map[k];
    if (nk <= 1.0) continue;
    double s1 = s1_map[k], s2 = s2_map[k];
    sse += s2 - (s1 * s1) / nk;
  }
  return sse;
}

// -----------------------------------------------------------------------------
// ale_sweep_cpp
// Purpose:
//   Find best split position t by sweeping rows (ordered by split feature) from
//   right to left. At each t, left child = rows 1..t, right = t+1..n_obs. Uses
//   sufficient statistics (n, s1, s2 per interval) for O(1) risk updates when
//   moving a row. For numeric split feature with self-ALE, optionally applies
//   boundary stabilizer to reduce boundary noise.
// Inputs:
//   ord_idx: 1-based row indices in sorted order (by z); length n_obs
//   d_l_mat: p x N, d_l_mat(j,i) = local effect for feature j, sample i
//   interval_idx_mat: p x N, interval index per feature per sample (1-based)
//   offsets: length p, start offset per feature in flattened tot_n/s1/s2 arrays
//   tot_n, tot_s1, tot_s2: full-node totals; r_n = right-node copy, l = tot - r
//   r_risks: initial right risks per feature (SSE per feature)
//   is_cand: length n_obs-1, TRUE where split position is a candidate (e.g. z changes)
//   min_node_size: minimum observations per child
//   split_feat_j: 1-based ALE feature index for split var, or 0 if no self-ALE
//   use_stabilizer: use boundary stabilizer (numeric + has_self_ale only)
//   z_sorted: z values in ord order; for categorical, all 0
//   n_obs: number of observations in node
// Output:
//   List: best_t (1-based), best_risks_sum, best_left_risks, best_right_risks.
//   If no valid split: best_t = NA_INTEGER, best_risks_sum = Inf.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List ale_sweep_cpp(
    IntegerVector ord_idx,
    NumericMatrix d_l_mat,
    IntegerMatrix interval_idx_mat,
    IntegerVector offsets,
    NumericVector tot_n,
    NumericVector tot_s1,
    NumericVector tot_s2,
    NumericVector r_risks,
    LogicalVector is_cand,
    int min_node_size,
    int split_feat_j,
    bool use_stabilizer,
    NumericVector z_sorted,
    int n_obs
) {
  const int p = d_l_mat.nrow();
  const int M = tot_n.size();
  const int j0 = split_feat_j - 1;  /* 0-based index of split feature. */

  /* Working copies: r_n, r_s1, r_s2 = right-node sufficient stats per interval. */
  std::vector<double> r_n(tot_n.begin(), tot_n.end());
  std::vector<double> r_s1(tot_s1.begin(), tot_s1.end());
  std::vector<double> r_s2(tot_s2.begin(), tot_s2.end());
  std::vector<double> left_risks(p, 0.0);
  std::vector<double> right_risks(r_risks.begin(), r_risks.end());

  double risks_sum = 0.0;
  for (int j = 0; j < p; ++j) risks_sum += right_risks[j];

  double best_risks_sum = R_PosInf;
  int best_t = -1;
  std::vector<double> best_left_risks(p), best_right_risks(p);

  bool has_self_ale = (split_feat_j >= 1 && split_feat_j <= p);

  /* Precompute d_l and interval_idx in sweep order for stabilizer. */
  NumericVector d_l_j_sorted(n_obs);
  IntegerVector interval_idx_sorted(n_obs);
  const int N = d_l_mat.ncol();
  if (use_stabilizer && has_self_ale && z_sorted.size() >= (R_xlen_t)n_obs) {
    for (int i = 0; i < n_obs; ++i) {
      int row = ord_idx[i] - 1;
      if (row < 0 || row >= N) stop("ord_idx contains invalid row index");
      d_l_j_sorted[i] = d_l_mat(j0, row);
      interval_idx_sorted[i] = interval_idx_mat(j0, row);
    }
  }
  /* Sweep: at each t, move row ord_idx[t-1] from right to left. */
  for (int t = 1; t <= n_obs - 1; ++t) {
    int row = ord_idx[t - 1] - 1;  /* 1-based to 0-based. */
    if (row < 0 || row >= N) stop("ord_idx contains invalid row index");

    /* Update sufficient stats and risks for each feature. */
    for (int j = 0; j < p; ++j) {
      double d = d_l_mat(j, row);
      int interval_idx_val = interval_idx_mat(j, row);
      int m = offsets[j] + interval_idx_val - 1;  /* Flattened interval index. */
      if (m < 0 || m >= M) continue;

      /* Remove row from right: r_n -= 1, r_s1 -= d, r_s2 -= d^2. */
      double r_n_old = r_n[m];
      double r_s1_old = r_s1[m];
      double r_s2_old = r_s2[m];
      double r_risk_old = risk_from_stats(r_n_old, r_s1_old, r_s2_old);

      double r_n_new = r_n_old - 1.0;
      double r_s1_new = r_s1_old - d;
      double r_s2_new = r_s2_old - d * d;
      double r_risk_new = risk_from_stats(r_n_new, r_s1_new, r_s2_new);

      r_n[m] = r_n_new;
      r_s1[m] = r_s1_new;
      r_s2[m] = r_s2_new;
      right_risks[j] -= r_risk_old;
      right_risks[j] += r_risk_new;

      /* Left stats = tot - right (invariant: tot = left + right). */
      double l_n_old = tot_n[m] - r_n_old;
      double l_s1_old = tot_s1[m] - r_s1_old;
      double l_s2_old = tot_s2[m] - r_s2_old;
      double l_risk_old = risk_from_stats(l_n_old, l_s1_old, l_s2_old);

      double l_n_new = tot_n[m] - r_n_new;
      double l_s1_new = tot_s1[m] - r_s1_new;
      double l_s2_new = tot_s2[m] - r_s2_new;
      double l_risk_new = risk_from_stats(l_n_new, l_s1_new, l_s2_new);

      left_risks[j] -= l_risk_old;
      left_risks[j] += l_risk_new;
      risks_sum += (-l_risk_old - r_risk_old + l_risk_new + r_risk_new);
    }

    /* Skip if not a candidate position or violates min_node_size. */
    if (!is_cand[t - 1] || t < min_node_size || (n_obs - t) < min_node_size)
      continue;

    /* Left/right constant? (all same z value) -> drop self risk. */
    bool l_const = (z_sorted.size() > 0 && std::abs(z_sorted[0] - z_sorted[t - 1]) < 1e-15);
    bool r_const = (z_sorted.size() >= (R_xlen_t)n_obs && std::abs(z_sorted[t] - z_sorted[n_obs - 1]) < 1e-15);
    bool do_stab = use_stabilizer && has_self_ale && !l_const && !r_const && t > kMinSideForStab && (n_obs - t) > kMinSideForStab;

    /* Compute total objective: sum of risks, minus self-ALE risk, plus stabilizer adj. */
    double total = R_PosInf;
    double adj_left = 0.0, adj_right = 0.0;

    if (!has_self_ale) {
      total = risks_sum;
    } else {
      double self_risk = left_risks[j0] + right_risks[j0];
      if (!use_stabilizer) {
        total = risks_sum - self_risk;
      } else if (!do_stab) {
        double drop = 0.0;
        if (l_const) drop += left_risks[j0];
        if (r_const) drop += right_risks[j0];
        total = risks_sum - drop;
      } else {
        NumericVector node_d_l_left(t);
        IntegerVector node_int_left(t);
        for (int i = 0; i < t; ++i) {
          node_d_l_left[i] = d_l_j_sorted[i];
          node_int_left[i] = interval_idx_sorted[i];
        }
        NumericVector node_d_l_right(n_obs - t);
        IntegerVector node_int_right(n_obs - t);
        for (int i = 0; i < n_obs - t; ++i) {
          node_d_l_right[i] = d_l_j_sorted[t + i];
          node_int_right[i] = interval_idx_sorted[t + i];
        }
        adj_left = adjust_side_cpp(node_d_l_left, node_int_left,
            left_risks[j0], true, t, t);
        adj_right = adjust_side_cpp(node_d_l_right, node_int_right,
            right_risks[j0], false, 0, n_obs - t);
        total = risks_sum - self_risk + adj_left + adj_right;
      }
    }

    if (!R_FINITE(total) || total >= best_risks_sum) continue;

    /* New best: store split index and risk vectors. */
    best_risks_sum = total;
    best_t = t;
    best_left_risks = left_risks;
    best_right_risks = right_risks;

    /* For self-ALE feature: store adjusted risks (or 0 if dropped). */
    if (has_self_ale) {
      if (!use_stabilizer) {
        best_left_risks[j0] = 0.0;
        best_right_risks[j0] = 0.0;
      } else if (!do_stab) {
        if (l_const) best_left_risks[j0] = 0.0;
        if (r_const) best_right_risks[j0] = 0.0;
      } else {
        best_left_risks[j0] = adj_left;
        best_right_risks[j0] = adj_right;
      }
    }
  }

  /* No valid split found. */
  if (best_t < 0) {
    return List::create(
      _["best_t"] = NA_INTEGER,
      _["best_risks_sum"] = R_PosInf,
      _["best_left_risks"] = NumericVector(p, NA_REAL),
      _["best_right_risks"] = NumericVector(p, NA_REAL)
    );
  }

  return List::create(
    _["best_t"] = best_t,  /* 1-based: split after row best_t. */
    _["best_risks_sum"] = best_risks_sum,
    _["best_left_risks"] = NumericVector(best_left_risks.begin(), best_left_risks.end()),
    _["best_right_risks"] = NumericVector(best_right_risks.begin(), best_right_risks.end())
  );
}
