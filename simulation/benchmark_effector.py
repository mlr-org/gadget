#!/usr/bin/env python3
"""effector efficiency benchmark: Global PDP, Global ALE, Regional PDP, Regional ALE.
Two model types: RF (sklearn) and toy (analytic DGP).
Run: python simulation/benchmark_effector.py [--datadir DIR] [--outdir DIR] [--reps N] [--predict-reps N]
Requires: pip install effector numpy pandas scikit-learn
"""

import argparse
import csv
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

try:
    import effector
except ImportError:
    sys.stderr.write("effector not installed. Run: pip install effector\n")
    sys.exit(1)

# Defaults match the **small** preset in efficiency_benchmark_plan.md §4.5 / run_benchmark.sh.
# Override via --N-vec / --D-vec / --fixed-N / --fixed-D (large preset uses larger N, D grids).
N_VEC = [500, 1000, 5000]
D_VEC = [5, 10, 20]
N_GRID_VEC = [10, 20, 50]
N_INT_VEC = [10, 20, 50]

# Explicit RF configuration for parity with simulation/benchmark_gadget.R.
# Keep this dict and R rf_config synchronized.
# Notes on non-1:1 mapping:
# - `min_samples_split` has no exact ranger equivalent (closest control is
#   child-node size via ranger `min.node.size`).
# - `min_impurity_decrease` and `ccp_alpha` have no direct ranger analogue.
# - `min_weight_fraction_leaf` is left at 0.0 (uniform sample weights here).
RF_CONFIG = {
    "n_estimators": 100,
    "criterion": "squared_error",
    "max_features": 1.0,
    "max_depth": None,
    "min_samples_split": 2,
    "min_samples_leaf": 1,
    "min_weight_fraction_leaf": 0.0,
    "max_leaf_nodes": None,
    "min_impurity_decrease": 0.0,
    "bootstrap": True,
    "max_samples": None,
    "ccp_alpha": 0.0,
    "random_state": 21,
    "n_jobs": 1,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_data(datadir, N, D, seed=21):
    f = os.path.join(datadir, "benchmark_N{}_D{}_seed{}.csv".format(N, D, seed))
    if not os.path.exists(f):
        return None, None
    df = pd.read_csv(f)
    X = df[["x{}".format(i + 1) for i in range(D)]].values
    y = df["y"].values
    return X, y


def fit_rf(X, y):
    return RandomForestRegressor(**RF_CONFIG).fit(X, y)


def rf_pred_fun(model, X):
    return model.predict(X)


def toy_pred_fun(X):
    x1, x2, x3 = X[:, 0], X[:, 1], X[:, 2]
    return 5.0 * x1 + 5.0 * x2 + np.where(x3 > 0, 10.0 * x1, 0.0) - np.where(x3 > 0, 10.0 * x2, 0.0)


def quantile_points(x, n_points):
    probs = np.linspace(0.0, 1.0, n_points)
    return np.quantile(x, probs)


# ---------------------------------------------------------------------------
# Method measurement functions
# ---------------------------------------------------------------------------

def measure_global_pdp(X, predict, D, n_grid):
    axis_limits = np.array([[-1.0] * D, [1.0] * D])
    m = effector.PDP(data=X, model=predict, axis_limits=axis_limits, nof_instances="all")
    tic = time.time()
    # Align with gadget::compute_pd(): compute effects across all D features.
    m.fit(features="all", centering=True, points_for_centering=n_grid)
    for feat in range(D):
        # Align with gadget grid construction: empirical quantiles per feature.
        xs = quantile_points(X[:, feat], n_grid)
        m.eval(feature=feat, xs=xs, centering=True, heterogeneity=False)
    return time.time() - tic


def measure_global_ale(X, predict, D, n_intervals):
    axis_limits = np.array([[-1.0] * D, [1.0] * D])
    m = effector.ALE(data=X, model=predict, axis_limits=axis_limits, nof_instances="all")
    bm = effector.axis_partitioning.Fixed(nof_bins=n_intervals, min_points_per_bin=0)
    tic = time.time()
    # Align with gadget::calculate_ale_fast(): compute effects across all D features.
    m.fit(features="all", centering=True, points_for_centering=100, binning_method=bm)
    for feat in range(D):
        xs = quantile_points(X[:, feat], 100)
        m.eval(feature=feat, xs=xs, centering=True, heterogeneity=False)
    return time.time() - tic


def measure_regional_pdp(X, predict, D, n_grid):
    """Return {"split": tree-only, "total": Phase1+Phase2}.

    Each phase is timed exactly once and no computation is repeated:

    Phase 1 (t_global) — mirrors gadget's system.time(compute_pd):
      Manually replicate the per-feature PDP precompute that RegionalPDP.fit()
      performs internally: create a PDP object, call pdp.fit() + pdp.eval(),
      and store ICE curves in m.y_ice.

    Phase 2 (t_split) — mirrors gadget's fit_timing$regional:
      Call _create_heterogeneity_function() + _fit_feature() per feature on the
      already-populated m.y_ice.

    total = t_global + t_split  (same definition as gadget-total)
    """
    axis_limits = np.array([[-1.0] * D, [1.0] * D])
    sp = effector.space_partitioning.Best(max_depth=2, min_samples_leaf=50)
    features = list(range(D))

    m = effector.RegionalPDP(data=X, model=predict, axis_limits=axis_limits, nof_instances="all")
    m.y_ice = {}

    # Phase 1: global PDP precompute — runs once.
    tic = time.time()
    for feat in features:
        pdp = effector.PDP(data=X, model=predict, axis_limits=axis_limits, nof_instances="all")
        try:
            pdp.fit(features=feat, centering=True, points_for_centering=n_grid)
        except TypeError:
            pdp.fit(features=feat)
        xx = np.linspace(axis_limits[0, feat], axis_limits[1, feat], n_grid)
        y_ice = pdp.eval(feature=feat, xs=xx, heterogeneity=True, return_all=True)
        m.y_ice["feature_" + str(feat)] = y_ice.T
    t_global = time.time() - tic

    # Phase 2: tree fitting only — runs once.
    tic = time.time()
    for feat in features:
        heter = m._create_heterogeneity_function(foi=feat, min_points=sp.min_points_per_subregion)
        m._fit_feature(feat, heter, sp, "all")
    t_split = time.time() - tic

    return {"split": t_split, "total": t_global + t_split}


def measure_regional_ale(X, predict, D, n_intervals):
    """Return {"split": tree-only, "total": Phase1+Phase2}.

    Each phase is timed exactly once and no computation is repeated:

    Phase 1 (t_global) — mirrors gadget's fit_timing$global:
      Manually replicate the per-feature ALE precompute that RegionalALE.fit()
      performs internally: create an ALE object, call global_ale.fit(), and store
      the effects in m.global_data_effect and m.global_bin_limits.

    Phase 2 (t_split) — mirrors gadget's fit_timing$regional:
      Call _create_heterogeneity_function() + _fit_feature() per feature on the
      already-populated cache.

    total = t_global + t_split  (same definition as gadget-total)
    """
    axis_limits = np.array([[-1.0] * D, [1.0] * D])
    sp = effector.space_partitioning.Best(max_depth=2, min_samples_leaf=50)
    bm = effector.axis_partitioning.Fixed(nof_bins=n_intervals, min_points_per_bin=0)
    features = list(range(D))

    m = effector.RegionalALE(data=X, model=predict, axis_limits=axis_limits, nof_instances="all")
    m.global_data_effect = {}
    m.global_bin_limits = {}

    # Phase 1: global ALE precompute — runs once.
    tic = time.time()
    for feat in features:
        global_ale = effector.ALE(data=X, model=predict, nof_instances="all", axis_limits=axis_limits)
        global_ale.fit(features=feat, binning_method=bm, centering=False)
        m.global_data_effect["feature_" + str(feat)] = global_ale.data_effect_ale["feature_" + str(feat)]
        m.global_bin_limits["feature_" + str(feat)] = global_ale.bin_limits["feature_" + str(feat)]
    t_global = time.time() - tic

    # Phase 2: tree fitting only — runs once.
    tic = time.time()
    for feat in features:
        heter = m._create_heterogeneity_function(feat, sp.min_points_per_subregion, n_intervals)
        m._fit_feature(feat, heter, sp, "all")
    t_split = time.time() - tic

    return {"split": t_split, "total": t_global + t_split}


# ---------------------------------------------------------------------------
# Predict baseline (RF only)
# ---------------------------------------------------------------------------

def run_predict_baseline(datadir, outdir, predict_reps):
    baseline = []
    for N in N_VEC:
        for D in D_VEC:
            X, y = load_data(datadir, N, D)
            if X is None:
                continue
            print("  baseline RF fit N={}, D={}...".format(N, D), flush=True)
            model = fit_rf(X, y)
            times = []
            for _ in range(predict_reps):
                tic = time.time()
                rf_pred_fun(model, X)
                times.append(time.time() - tic)
            baseline.append({
                "package": "effector", "N": N, "D": D,
                "predict_time_mean": np.mean(times),
                "predict_time_sd": np.std(times),
                "n_rep": predict_reps
            })
    if baseline:
        df = pd.DataFrame(baseline)
        df.to_csv(os.path.join(outdir, "predict_baseline_effector.csv"), index=False)
        print("Written: predict_baseline_effector.csv", flush=True)


# ---------------------------------------------------------------------------
# Sweep runner
# ---------------------------------------------------------------------------

METHODS = {
    "global_pdp": {"runner": measure_global_pdp, "is_pdp": True},
    "global_ale": {"runner": measure_global_ale, "is_pdp": False},
    "regional_pdp": {"runner": measure_regional_pdp, "is_pdp": True},
    "regional_ale": {"runner": measure_regional_ale, "is_pdp": False},
}


def run_sweep(method_name, model_type, datadir, reps, fixed_N, fixed_D):
    info = METHODS[method_name]
    runner = info["runner"]
    is_pdp = info["is_pdp"]
    res_vec = N_GRID_VEC if is_pdp else N_INT_VEC
    default_res = 20
    rows = []

    def make_predict(model):
        if model_type == "rf":
            return lambda x: rf_pred_fun(model, x)
        else:
            return toy_pred_fun

    def record(t_raw, N, D, resolution, r):
        """Handle scalar or dict returns from runners.

        Dict runners (regional methods) emit one row per key with the key appended to
        the method name, e.g. "regional_pdp_split" and "regional_pdp_total".
        This mirrors the named-vector expansion in benchmark_gadget.R's run_sweep.
        """
        if isinstance(t_raw, dict):
            for nm, t in t_raw.items():
                rows.append({
                    "package": "effector",
                    "method": "{}_{}".format(method_name, nm),
                    "N": N, "D": D,
                    "n_grid": resolution if is_pdp else "",
                    "n_intervals": resolution if not is_pdp else "",
                    "repetition": r + 1,
                    "time_sec": t,
                })
        else:
            rows.append({
                "package": "effector",
                "method": method_name,
                "N": N, "D": D,
                "n_grid": resolution if is_pdp else "",
                "n_intervals": resolution if not is_pdp else "",
                "repetition": r + 1,
                "time_sec": t_raw,
            })

    print("  [{}] {} -- vs N".format(model_type, method_name), flush=True)
    for N in N_VEC:
        X, y = load_data(datadir, N, fixed_D)
        if X is None:
            continue
        model = fit_rf(X, y) if model_type == "rf" else None
        predict = make_predict(model)
        try:
            runner(X, predict, fixed_D, default_res)  # warmup
        except Exception:
            pass
        for r in range(reps):
            try:
                t = runner(X, predict, fixed_D, default_res)
            except Exception as e:
                print("    Error: {}".format(e), flush=True)
                t = float("nan")
            record(t, N, fixed_D, default_res, r)

    print("  [{}] {} -- vs D".format(model_type, method_name), flush=True)
    for D in D_VEC:
        X, y = load_data(datadir, fixed_N, D)
        if X is None:
            continue
        model = fit_rf(X, y) if model_type == "rf" else None
        predict = make_predict(model)
        try:
            runner(X, predict, D, default_res)  # warmup
        except Exception:
            pass
        for r in range(reps):
            try:
                t = runner(X, predict, D, default_res)
            except Exception as e:
                print("    Error: {}".format(e), flush=True)
                t = float("nan")
            record(t, fixed_N, D, default_res, r)

    print("  [{}] {} -- vs resolution".format(model_type, method_name), flush=True)
    for rv in res_vec:
        X, y = load_data(datadir, fixed_N, fixed_D)
        if X is None:
            continue
        model = fit_rf(X, y) if model_type == "rf" else None
        predict = make_predict(model)
        try:
            runner(X, predict, fixed_D, rv)  # warmup
        except Exception:
            pass
        for r in range(reps):
            try:
                t = runner(X, predict, fixed_D, rv)
            except Exception as e:
                print("    Error: {}".format(e), flush=True)
                t = float("nan")
            record(t, fixed_N, fixed_D, rv, r)

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--datadir", default="simulation/data/benchmark")
    p.add_argument("--outdir", default="simulation/results/benchmark")
    p.add_argument("--reps", type=int, default=5)
    p.add_argument("--predict-reps", type=int, default=20)
    p.add_argument("--fixed-N", type=int, default=1000)
    p.add_argument("--fixed-D", type=int, default=10)
    p.add_argument("--N-vec", default="500,1000,5000")
    p.add_argument("--D-vec", default="5,10,20")
    p.add_argument("--n-grid-vec", default="10,20,50")
    p.add_argument("--n-int-vec", default="10,20,50")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    def parse_int_vec(raw):
        return [int(x) for x in raw.split(",") if x.strip()]

    global N_VEC, D_VEC, N_GRID_VEC, N_INT_VEC
    N_VEC = parse_int_vec(args.N_vec)
    D_VEC = parse_int_vec(args.D_vec)
    N_GRID_VEC = parse_int_vec(args.n_grid_vec)
    N_INT_VEC = parse_int_vec(args.n_int_vec)

    print("=== Predict baseline (RF) ===", flush=True)
    run_predict_baseline(args.datadir, args.outdir, args.predict_reps)

    cols = ["package", "method", "N", "D", "n_grid", "n_intervals", "repetition", "time_sec"]

    def write_results(rows, filename):
        if not rows:
            print("No results for {}".format(filename), flush=True)
            return
        out = os.path.join(args.outdir, filename)
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print("Written: {}".format(out), flush=True)

    print("=== RF model benchmarks ===", flush=True)
    rf_rows = []
    for method_name in METHODS:
        rf_rows.extend(run_sweep(method_name, "rf", args.datadir, args.reps, args.fixed_N, args.fixed_D))
    write_results(rf_rows, "benchmark_effector_rf.csv")

    print("=== Toy model benchmarks ===", flush=True)
    toy_rows = []
    for method_name in METHODS:
        toy_rows.extend(run_sweep(method_name, "toy", args.datadir, args.reps, args.fixed_N, args.fixed_D))
    write_results(toy_rows, "benchmark_effector_toy.csv")

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
