/**
 * @file search_best_split.cpp
 * @brief Fast tree splitting for PD strategy - C++/Armadillo implementation.
 *
 * Categorical and numerical splitting with preprocessed effect matrices.
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Helper functions (internal)
// -----------------------------------------------------------------------------

/* Remove consecutive duplicates; assumes x is sorted. O(n). Copies input so caller's x is unchanged. */
inline NumericVector unique_cpp(const NumericVector& x) {
  NumericVector out = Rcpp::clone(x);
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

/* R type-7 quantile for sorted vector (R default). */
inline double quantile_type7(const NumericVector& x, double p) {
  const int n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];
  if (p <= 0.0) return x[0];
  if (p >= 1.0) return x[n - 1];

  const double h = (n - 1) * p;
  const int lo   = static_cast<int>(h);
  const double f = h - lo;

  return x[lo] + f * (x[lo + 1] - x[lo]);
}

/* Convert R matrix to Armadillo (zero-copy when possible). */
inline arma::mat arma_view(SEXP obj) {
  if (!Rf_isMatrix(obj) || TYPEOF(obj) != REALSXP) {
    Rcpp::Function as_matrix("as.matrix");
    obj = as_matrix(obj);
  }
  Rcpp::NumericMatrix M(obj);
  return arma::mat(M.begin(), M.nrow(), M.ncol(), false);
}

// -----------------------------------------------------------------------------
// search_best_split_point_cpp_internal
// Purpose:
//   Find best split for a single feature (categorical or numerical).
//   Accepts preprocessed Ym and S_tot for efficiency.
// Inputs:
//   z: Feature vector (numeric or categorical)
//   Ym: Preprocessed effect matrices (NaN replaced with 0)
//   S_tot: Column sums per matrix
//   n_quantiles, is_categorical, min_node_size
// Output:
//   List: split_point, split_objective
// -----------------------------------------------------------------------------
List search_best_split_point_cpp_internal(
    SEXP              z,
    const std::vector<arma::mat>& Ym,        // & = reference, const = read-only, avoids copying large matrices
    const std::vector<arma::rowvec>& S_tot,  // & = reference, const = read-only, avoids copying large vectors
    Nullable<int>     n_quantiles   = R_NilValue,
    bool              is_categorical = false,
    int               min_node_size  = 1)
{
  const int Ly = Ym.size(); // p
  const int N = Ym[0].n_rows; // n

  double best_obj = R_PosInf, best_split = NA_REAL;
  std::string best_level;

  /* ================= Categorical Variable Splitting ================= */
  if (is_categorical) {
    // Extract factor levels and convert to integer vector
    IntegerVector z_fac(z);
    CharacterVector lev = z_fac.attr("levels");
    const int K = lev.size();

    // If only one level, no valid split possible
    if (K <= 1)
      return List::create(_["split_point"] = R_NaString,
        _["split_objective"] = R_PosInf);

    // Calculate per-level left sums for each effect matrix
    std::vector< std::vector<arma::rowvec> > SumL(Ly);
    std::vector<int> countL(K, 0);
    for (int l = 0; l < Ly; ++l)
      SumL[l].assign(K, arma::rowvec(Ym[l].n_cols, arma::fill::zeros));

    // Accumulate sums for each level (skip NA)
    for (int i = 0; i < N; ++i) {
      if (z_fac[i] == NA_INTEGER) continue;
      int k = z_fac[i] - 1;
      ++countL[k];
      for (int l = 0; l < Ly; ++l) {
        SumL[l][k] += Ym[l].row(i);  // Direct accumulation, NaN already processed in preprocessing
      }
    }

    // Evaluate each level as potential split (one vs rest)
    for (int k = 0; k < K; ++k) {
      int NL = countL[k], NR = N - NL;

      // Skip if either child node would be too small
      if (NL < min_node_size || NR < min_node_size) continue;

      // Calculate objective function for this split
      double obj = 0.0;
      for (int l = 0; l < Ly; ++l) {
        const arma::rowvec SL = SumL[l][k];
        const arma::rowvec SR = S_tot[l] - SL;
        obj += arma::accu( - SL%SL / NL - SR%SR / NR );
      }

      // Update best split if this one is better
      if (obj < best_obj) {
        best_obj = obj;
        best_level = as<std::string>(lev[k]);
      }
    }

    return List::create(_["split_point"]     = best_level,
      _["split_objective"] = best_obj);
  }

  /* ================= Numerical Variable Splitting =================== */

  // Convert to numeric and sort the feature vector once
  NumericVector z_num;
  if (TYPEOF(z) == STRSXP) {
    // Convert character to numeric if possible
    CharacterVector z_char(z);
    z_num = NumericVector(z_char.size());
    for (int i = 0; i < z_char.size(); ++i) {
      if (z_char[i] == NA_STRING) {
        z_num[i] = NA_REAL;
      } else {
        z_num[i] = std::stod(as<std::string>(z_char[i]));
      }
    }
  } else {
    z_num = NumericVector(z);
  }

  // Create sorted index and sorted values
  IntegerVector ord = Rcpp::seq(0, N - 1);
  std::sort(ord.begin(), ord.end(), [&](int i, int j){ return z_num[i] < z_num[j]; });

  NumericVector z_sorted(N);
  for (int i = 0; i < N; ++i) z_sorted[i] = z_num[ord[i]];

  // Generate candidate split points
  NumericVector splits;
  if (n_quantiles.isNotNull()) {
    // Use quantile-based grid
    int nq = Rcpp::as<int>(n_quantiles);
    NumericVector uniq = unique_cpp(z_sorted);
    if (uniq.size() < nq) {
      splits = uniq;
    } else {
      // Calculate quantile positions
      NumericVector q(nq);
      for (int i = 0; i < nq; ++i) q[i] = (i + 1.0)/(nq + 1.0);
      splits = NumericVector(nq);
      for (int i = 0; i < nq; ++i)
        splits[i] = quantile_type7(z_sorted, q[i]);
      splits = unique_cpp(splits);
    }
  } else {
    // Use all unique values as candidates
    splits = unique_cpp(z_sorted);
  }

  // Check if we have any valid split candidates
  if (splits.size() == 0)
    return List::create(_["split_point"]     = NA_REAL,
      _["split_objective"] = R_PosInf);

  // Stream through split candidates and accumulate left sums (incremental; splits are sorted)
  std::vector<arma::rowvec> SL(Ly);
  for (int l = 0; l < Ly; ++l)
    SL[l] = arma::rowvec(Ym[l].n_cols, arma::fill::zeros);

  int idx = 0;
  for (double sp : splits) {
    while (idx < N && z_sorted[idx] <= sp) {
      int r = ord[idx++];
      for (int l = 0; l < Ly; ++l)
        SL[l] += Ym[l].row(r);
    }
    int NL = idx, NR = N - NL;

    // Skip if either child node would be too small
    if (NL < min_node_size || NR < min_node_size) continue;

    // Calculate objective function for this split
    double obj = 0.0;
    for (int l = 0; l < Ly; ++l) {
      const arma::rowvec SR = S_tot[l] - SL[l];
      obj += arma::accu( - SL[l]%SL[l] / NL - SR%SR / NR );
    }
    if (obj < best_obj) { best_obj = obj; best_split = sp; }
  }

  if (best_obj == R_PosInf || R_IsNA(best_split))
    return List::create(_["split_point"]     = NA_REAL,
      _["split_objective"] = R_PosInf);

  // Refine split point to midpoint between adjacent values
  double Lft = -std::numeric_limits<double>::infinity(),
    Rgt =  std::numeric_limits<double>::infinity();
  for (int i = 0; i < N; ++i) {
    double v = z_num[i];
    if (!R_IsNA(v) && v <= best_split && v > Lft) Lft = v;
    if (!R_IsNA(v) && v >  best_split && v < Rgt) Rgt = v;
  }
  double mid = std::isinf(Rgt) ? Lft : (Lft + Rgt) / 2.0;

  return List::create(_["split_point"]     = mid,
    _["split_objective"] = best_obj);
}

// -----------------------------------------------------------------------------
// search_best_split_cpp
// Purpose:
//   Evaluate all features in Z, find best split per feature, return full results.
// Inputs:
//   Z: DataFrame of features
//   Y: List of effect matrices
//   min_node_size, n_quantiles
// Output:
//   DataFrame: split_feature, is_categorical, split_point, split_objective, etc.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame search_best_split_cpp(
    DataFrame       Z,
    List            Y,
    int             min_node_size,
    Nullable<int>   n_quantiles = R_NilValue)
{
  // Initialize output vectors
  const int p = Z.size();
  CharacterVector feat_names = Z.names();

  CharacterVector split_feature(p);
  LogicalVector   is_cat_vec(p);
  CharacterVector split_point_out(p);
  NumericVector   split_obj(p);

  // Preprocess all Y matrices' NaN values to avoid repeated processing
  const int Ly = Y.size();
  std::vector<arma::mat> Ym(Ly);
  std::vector<arma::rowvec> S_tot(Ly);
  for (int l = 0; l < Ly; ++l) {
    Ym[l] = arma_view(Y[l]);
    Ym[l].replace(arma::datum::nan, 0.0);  // Process all NaN values once
    S_tot[l] = arma::sum(Ym[l], 0);        // Precompute column sums
  }
  // Evaluate each feature
  for (int j = 0; j < p; ++j) {
    SEXP z_j  = Z[j];
    bool is_c = Rf_isFactor(z_j);

    // Call internal function directly, passing preprocessed data
    List res = search_best_split_point_cpp_internal(
      z_j, Ym, S_tot, n_quantiles, is_c, min_node_size);

    // Store results
    split_feature[j]   = feat_names[j];
    is_cat_vec[j]      = is_c;
    split_obj[j]       = res["split_objective"];
    split_point_out[j] = as<CharacterVector>(wrap(res["split_point"]))[0];
  }
  // Find best valid split.
  // The internal function always returns R_PosInf when no split is found, so
  // split_obj[j] < best_val (with best_val = R_PosInf) is sufficient: R_PosInf < R_PosInf
  // is false, correctly excluding invalid features without inspecting the split point string.
  double best_val = R_PosInf;
  int best_idx = -1;

  for (int j = 0; j < p; ++j) {
    if (split_obj[j] < best_val) {
      best_val = split_obj[j];
      best_idx = j;
    }
  }
  
  LogicalVector best_split(p, false);
  if (best_idx >= 0) {
    best_split[best_idx] = true;
  }

  // Return results as DataFrame
  return DataFrame::create(
    _["split_feature"]   = split_feature,
    _["is_categorical"]  = is_cat_vec,
    _["split_point"]     = split_point_out,
    _["split_objective"] = split_obj,
    _["best_split"]      = best_split
  );
}
