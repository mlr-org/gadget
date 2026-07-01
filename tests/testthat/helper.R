skip_ale_cpp_if_unavailable = function() {
  n = 5
  dt = data.table::data.table(
    row_id = seq_len(n), interval_index = rep(1L, n), d_l = 0, int_n = n, int_s1 = 0, int_s2 = 0
  )
  tryCatch({
    xplaineff:::calculate_ale_heterogeneity_list_cpp(list(x = dt))
  }, error = function(e) {
    if (grepl("not available for .Call", conditionMessage(e), fixed = TRUE)) {
      testthat::skip("ALE C++ symbols not loaded (install package with compile)")
    }
  })
}
