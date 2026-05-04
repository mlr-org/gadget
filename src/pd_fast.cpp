#include <Rcpp.h>
#include <cstring>

// PDP/ICE stacked newdata construction (xplaineff-style), adapted from xplaineff/src/pdp.c.
// Builds an n*m-row table: for each grid point k, rows1..n are a copy of the original data with
// the focal feature replaced by grid[k]. Non-FOI columns are block-replicated via memcpy where possible.

static R_xlen_t pd_list_nrow(const Rcpp::List& cols) {
  if (cols.size() == 0) {
    return 0;
  }
  return Rf_xlength(cols[0]);
}

static SEXP make_pdp_feature_col_mx(SEXP col, SEXP grid, R_xlen_t n, R_xlen_t m, bool is_foi) {
  const R_xlen_t nm = n * m;
  if (!is_foi) {
    switch (TYPEOF(col)) {
    case REALSXP: {
      SEXP out = PROTECT(Rf_allocVector(REALSXP, nm));
      Rf_copyMostAttrib(col, out);
      const size_t bytes = sizeof(double) * (size_t)n;
      const double* src = REAL(col);
      double* dest = REAL(out);
      for (R_xlen_t k = 0; k < m; k++) {
        std::memcpy(dest + k * n, src, bytes);
      }
      UNPROTECT(1);
      return out;
    }
    case INTSXP: {
      SEXP out = PROTECT(Rf_allocVector(INTSXP, nm));
      Rf_copyMostAttrib(col, out);
      const size_t bytes = sizeof(int) * (size_t)n;
      const int* src = INTEGER(col);
      int* dest = INTEGER(out);
      for (R_xlen_t k = 0; k < m; k++) {
        std::memcpy(dest + k * n, src, bytes);
      }
      UNPROTECT(1);
      return out;
    }
    case LGLSXP: {
      SEXP out = PROTECT(Rf_allocVector(LGLSXP, nm));
      Rf_copyMostAttrib(col, out);
      const size_t bytes = sizeof(int) * (size_t)n;
      const int* src = LOGICAL(col);
      int* dest = LOGICAL(out);
      for (R_xlen_t k = 0; k < m; k++) {
        std::memcpy(dest + k * n, src, bytes);
      }
      UNPROTECT(1);
      return out;
    }
    case STRSXP: {
      SEXP out = PROTECT(Rf_allocVector(STRSXP, nm));
      Rf_copyMostAttrib(col, out);
      for (R_xlen_t k = 0; k < m; k++) {
        for (R_xlen_t i = 0; i < n; i++) {
          SET_STRING_ELT(out, k * n + i, STRING_ELT(col, i));
        }
      }
      UNPROTECT(1);
      return out;
    }
    default:
      Rcpp::stop("Unsupported column type in data.frame (type = %s).", Rf_type2char(TYPEOF(col)));
    }
  }

  SEXP out;
  // Integer storage in `data`: PD grid from R is usually double (quantiles). The old path used
  // INTEGER(grid) only, which does not match REAL grids and would force lossy coercion. When grid
  // is double, promote the stacked focal column to REAL so each grid value is exact for predict().
  if (TYPEOF(col) == INTSXP) {
    if (TYPEOF(grid) == REALSXP) {
      out = PROTECT(Rf_allocVector(REALSXP, nm));
      Rf_copyMostAttrib(col, out);
      double* outp = REAL(out);
      const double* grid_dbl = REAL(grid);
      for (R_xlen_t k = 0; k < m; k++) {
        const double val = grid_dbl[k];
        double* dest = outp + k * n;
        dest[0] = val;
        R_xlen_t filled = 1;
        while (filled * 2 <= n) {
          std::memcpy(dest + filled, dest, sizeof(double) * (size_t)filled);
          filled *= 2;
        }
        if (filled < n) {
          std::memcpy(dest + filled, dest, sizeof(double) * (size_t)(n - filled));
        }
      }
    } else if (TYPEOF(grid) == INTSXP) {
      out = PROTECT(Rf_allocVector(INTSXP, nm));
      Rf_copyMostAttrib(col, out);
      int* outp = INTEGER(out);
      const int* grid_int = INTEGER(grid);
      for (R_xlen_t k = 0; k < m; k++) {
        const int val = grid_int[k];
        int* dest = outp + k * n;
        dest[0] = val;
        R_xlen_t filled = 1;
        while (filled * 2 <= n) {
          std::memcpy(dest + filled, dest, sizeof(int) * (size_t)filled);
          filled *= 2;
        }
        if (filled < n) {
          std::memcpy(dest + filled, dest, sizeof(int) * (size_t)(n - filled));
        }
      }
    } else {
      Rcpp::stop("grid must be numeric or integer for an integer focal feature.");
    }
  } else {
    out = PROTECT(Rf_allocVector(REALSXP, nm));
    Rf_copyMostAttrib(col, out);
    double* outp = REAL(out);
    auto fill_blocks = [&](auto grid_val_at) {
      for (R_xlen_t k = 0; k < m; k++) {
        const double val = static_cast<double>(grid_val_at(k));
        double* dest = outp + k * n;
        dest[0] = val;
        R_xlen_t filled = 1;
        while (filled * 2 <= n) {
          std::memcpy(dest + filled, dest, sizeof(double) * (size_t)filled);
          filled *= 2;
        }
        if (filled < n) {
          std::memcpy(dest + filled, dest, sizeof(double) * (size_t)(n - filled));
        }
      }
    };
    if (TYPEOF(grid) == INTSXP) {
      const int* grid_int = INTEGER(grid);
      fill_blocks([&](R_xlen_t k) { return grid_int[k]; });
    } else if (TYPEOF(grid) == REALSXP) {
      const double* grid_dbl = REAL(grid);
      fill_blocks([&](R_xlen_t k) { return grid_dbl[k]; });
    } else {
      Rcpp::stop("grid must be numeric or integer for a numeric focal feature.");
    }
  }
  UNPROTECT(1);
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_pd_stack_newdata(Rcpp::List data_cols, int feature_index, SEXP grid) {
  const int p = data_cols.size();
  const R_xlen_t n = pd_list_nrow(data_cols);
  const R_xlen_t m = Rf_xlength(grid);
  if (feature_index < 0 || feature_index >= p) {
    Rcpp::stop("feature_index out of bounds.");
  }
  if (m < 1) {
    Rcpp::stop("grid must have length >= 1.");
  }
  SEXP names_attr = Rf_getAttrib(data_cols, R_NamesSymbol);
  if (names_attr == R_NilValue || TYPEOF(names_attr) != STRSXP || (int)XLENGTH(names_attr) != p) {
    Rcpp::stop("data_cols must be a named list (one column per element).");
  }

  Rcpp::List out(p);
  Rf_setAttrib(out, R_NamesSymbol, names_attr);
  for (int col = 0; col < p; col++) {
    SEXP s_data_col = VECTOR_ELT(data_cols, col);
    SEXP s_new = PROTECT(make_pdp_feature_col_mx(s_data_col, grid, n, m, col == feature_index));
    SET_VECTOR_ELT(out, col, s_new);
    UNPROTECT(1);
  }

  const R_xlen_t nm = n * m;
  SEXP s_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(s_rownames)[0] = NA_INTEGER;
  INTEGER(s_rownames)[1] = static_cast<int>(-nm);
  Rf_setAttrib(out, R_RowNamesSymbol, s_rownames);
  UNPROTECT(1);

  SEXP cls = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  UNPROTECT(1);
  return out;
}
