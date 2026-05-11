library(data.table)
dt = fread("simulation/results/benchmark_large/summary.csv")
global = dt[grepl("^global", method)]
global[, pkg := fcase(
  package == "effector", "effector",
  impl == "effect-cpp",  "gadget_cpp",
  impl == "effect-r",    "gadget_r"
)]

ratio_table = function(sub, xvar, res_var, res_val, fixed_other, fixed_other_val) {
  s = sub[get(res_var) == res_val & get(fixed_other) == fixed_other_val,
    .(pkg, model_type, x = get(xvar), time = time_mean)]
  w = dcast(s, x + model_type ~ pkg, value.var = "time", fun.aggregate = mean)
  if (!all(c("gadget_cpp","gadget_r","effector") %in% names(w))) return(w)
  w[, cpp_eff := round(gadget_cpp / effector,   2)]
  w[, r_eff   := round(gadget_r   / effector,   2)]
  w[, r_cpp   := round(gadget_r   / gadget_cpp, 2)]
  setnames(w, "x", xvar)
  w[order(model_type, get(xvar))]
}

ale = global[method == "global_ale"]
pd  = global[method == "global_pdp"]

cat("══════ ALE vs N (D=20, intervals=20) ══════\n")
print(ratio_table(ale, "N", "n_intervals", 20, "D", 20))
cat("\n══════ ALE vs D (N=10000, intervals=20) ══════\n")
print(ratio_table(ale, "D", "n_intervals", 20, "N", 10000))
cat("\n══════ ALE vs Intervals (N=10000, D=20) ══════\n")
print(ratio_table(ale, "n_intervals", "D", 20, "N", 10000))
cat("\n══════ PD vs N (D=20, grid=20) ══════\n")
print(ratio_table(pd, "N", "n_grid", 20, "D", 20))
cat("\n══════ PD vs D (N=10000, grid=20) ══════\n")
print(ratio_table(pd, "D", "n_grid", 20, "N", 10000))
cat("\n══════ PD vs Grid size (N=10000, D=20) ══════\n")
print(ratio_table(pd, "n_grid", "D", 20, "N", 10000))
